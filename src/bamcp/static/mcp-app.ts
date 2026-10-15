/**
 * BAMCP Alignment Viewer - MCP Apps client
 *
 * Production-quality BAM/CRAM visualization with:
 * - Reference sequence track with semantic zoom
 * - Position ruler with genomic coordinates
 * - MAPQ-based read opacity
 * - Insertion/deletion visualization from CIGAR
 * - Keyboard shortcuts for navigation
 *
 * This module is bundled into viewer.html by Vite.
 */

import { BAMCPClient, DebugInfo, FetchResult } from "./client";
import {
    baseColor,
    displayPos,
    escapeHtml,
    formatNearEndSummary,
    maxOf,
    paddedFetchWindow,
    panViewport,
    svColor,
    svTypeLabel,
    toRegionString,
    zoomViewport,
} from "./constants";
import { variantKey } from "./data-store";
import {
    firstHiddenRow,
    hiddenReadCount,
    hiddenReadsLabel,
    referenceBannerMessage,
} from "./ui-sync";
import { Renderer } from "./renderer";
import { StateManager, DEFAULT_SETTINGS } from "./state";
import {
    ArtifactAssessment,
    ColorBy,
    Read,
    ReadDisplayMode,
    SortBy,
    Variant,
    VariantEvidence,
    ViewerSettings,
} from "./types";

const CONTEXT_UPDATE_DEBOUNCE_MS = 300;
// Fixed height of the in-viewer "N reads not shown" bar (matches viewer.html CSS).
const HIDDEN_READS_BAR_PX = 22;

/** Render per-sample genotypes as `name=allele/allele` (missing alleles as `.`). */
function formatGenotypes(samples?: Record<string, Record<string, unknown>>): string[] {
    if (!samples) return [];
    return Object.entries(samples).map(([name, data]) => {
        const gt = (data as Record<string, unknown> | null)?.GT;
        let gtStr = '.';
        if (Array.isArray(gt)) {
            gtStr = gt.map(a => (a == null ? '.' : String(a))).join('/');
        } else if (gt != null) {
            gtStr = String(gt);
        }
        return `${name}=${gtStr}`;
    });
}

class BAMCPViewer {
    private client: BAMCPClient;
    private state: StateManager;
    private renderer: Renderer;

    // DOM Elements
    private tooltip: HTMLElement;
    private variantTable: HTMLElement;
    private evidencePanel: HTMLElement;
    private readsCanvas: HTMLCanvasElement;
    private sequenceNotice: HTMLElement;
    private debugOverlay: HTMLElement;
    private debugContext: HTMLElement;
    private debugToolCall: HTMLElement;

    // Interaction State
    private _contextUpdateTimer: ReturnType<typeof setTimeout> | null = null;
    private hoverDelayTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly HOVER_DELAY_MS = 200;
    private isFullscreen = false;

    // Tooltip drag state
    private isDraggingTooltip = false;
    private tooltipDragOffset = { x: 0, y: 0 };

    // Sequence loading state
    private pendingSequenceRequest = false;
    private sequenceFallbackTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly BASE_RENDER_SCALE = 10;
    private static readonly SEQUENCE_REGION_THRESHOLD = 500;
    private static readonly SEQUENCE_FETCH_TIMEOUT_MS = 3000;

    // Viewport refetch state — auto-load data when panning/zooming beyond loaded range
    private viewportRefetchTimer: ReturnType<typeof setTimeout> | null = null;
    private variantTableTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly VIEWPORT_REFETCH_OVERLAP = 0.8;     // refetch when <80% overlap
    private static readonly VIEWPORT_REFETCH_DEBOUNCE_MS = 300;
    private fetchErrorTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly FETCH_ERROR_MS = 8000;

    constructor() {
        this.client = new BAMCPClient();
        this.state = new StateManager();
        this.renderer = new Renderer(this.state);
        // Every render re-syncs the DOM chrome, so any path that changes state and
        // renders (host data, pan/zoom, or a direct hydrate that sets
        // settings.displayMode and calls renderer.resize()) leaves the mode select,
        // region box, banner and hidden-reads indicator consistent with the canvas.
        this.renderer.onRendered = () => this.syncUiState();

        // Bind DOM elements
        this.tooltip = document.getElementById('tooltip')!;
        this.variantTable = document.getElementById('variant-table')!;
        this.evidencePanel = document.getElementById('evidence-panel')!;
        this.readsCanvas = document.getElementById('reads-canvas') as HTMLCanvasElement;
        this.sequenceNotice = document.getElementById('sequence-notice')!;
        this.debugOverlay = document.getElementById('debug-overlay')!;
        this.debugContext = document.getElementById('debug-context')!;
        this.debugToolCall = document.getElementById('debug-toolcall')!;

        // Host-initiated data callback (ontoolresult only — NOT app-initiated fetches).
        // App-initiated fetches (viewport refetch, sequence zoom) return data via
        // Promise and call loadTile directly, keeping the paths cleanly separated.
        this.client.setOnDataReceived(data => {
            // Cancel any pending viewport refetch timer (fresh host data arrived)
            if (this.viewportRefetchTimer) {
                clearTimeout(this.viewportRefetchTimer);
                this.viewportRefetchTimer = null;
            }
            this.sequenceNotice.classList.add('hidden');
            this.sequenceNotice.classList.remove('loading');
            this.pendingSequenceRequest = false;
            if (this.sequenceFallbackTimer) {
                clearTimeout(this.sequenceFallbackTimer);
                this.sequenceFallbackTimer = null;
            }
            this.hideFetchError();

            // Host data always resets viewport
            this.state.loadData(data);
            this.syncRegionInput();
            this.renderVariantTable();
            this.renderer.resize();
            this.scheduleContextUpdate();
        });

        // Debug overlay callback — update content but don't auto-show
        this.client.setOnDebugUpdate((info: DebugInfo) => {
            this.debugContext.textContent = info.lastContext || '—';
            this.debugToolCall.textContent = info.lastToolCall || '—';
        });

        document.getElementById('fetch-error-close')?.addEventListener('click', () => {
            this.hideFetchError();
        });

        // Debug close button
        document.getElementById('debug-close')!.addEventListener('click', () => {
            this.debugOverlay.classList.add('hidden');
        });

        // Debug toggle button in toolbar
        document.getElementById('debug-btn')!.addEventListener('click', () => {
            this.debugOverlay.classList.toggle('hidden');
        });

        this.setupEventListeners();
        this.setupKeyboardShortcuts();

        // The reads area scrolls inside #viewer; keep the hidden-reads count current.
        document.getElementById('viewer')?.addEventListener('scroll', () => {
            this.updateHiddenReadsIndicator();
        });

        // Window resize
        window.addEventListener('resize', () => this.renderer.resize());
        // A devicePixelRatio change without a resize (dragging the window to a
        // monitor with a different density, or browser zoom on some platforms)
        // would leave the backing stores at the old density and render blurry.
        this.watchDevicePixelRatio();

        // Init client and update UI based on available modes
        this.client.init().then(() => {
            this.updateDisplayModeButtons();
        });

        // Initial render
        this.renderer.resize();
    }

    /**
     * Re-run resize() whenever devicePixelRatio changes. A resolution media query
     * only matches the DPR it was built for, so re-arm it with the new value on
     * every change.
     */
    private watchDevicePixelRatio(): void {
        if (typeof window.matchMedia !== 'function') return;
        const mql = window.matchMedia(`(resolution: ${window.devicePixelRatio || 1}dppx)`);
        const onChange = () => {
            mql.removeEventListener?.('change', onChange);
            this.renderer.resize();
            this.watchDevicePixelRatio();
        };
        mql.addEventListener?.('change', onChange);
    }

    private updateDisplayModeButtons(): void {
        const availableModes = this.client.getAvailableDisplayModes();
        const pipBtn = document.getElementById('pip-btn') as HTMLButtonElement;
        const fsBtn = document.getElementById('fullscreen-btn') as HTMLButtonElement;

        // Disable PIP button if not supported by host
        if (!availableModes.includes('pip')) {
            pipBtn.disabled = true;
            pipBtn.title = 'PIP mode not supported by host';
            pipBtn.style.opacity = '0.5';
        }

        // Disable fullscreen button if not supported
        if (!availableModes.includes('fullscreen')) {
            fsBtn.disabled = true;
            fsBtn.title = 'Fullscreen not supported by host';
            fsBtn.style.opacity = '0.5';
        }
    }

    private setupEventListeners(): void {
        document.getElementById('go-btn')!.addEventListener('click', () => {
            const region = (document.getElementById('region-input') as HTMLInputElement).value;
            this.client.requestRegion(region);
        });

        // Gene search
        document.getElementById('gene-btn')!.addEventListener('click', () => {
            const gene = (document.getElementById('gene-input') as HTMLInputElement).value.trim();
            if (gene) this.client.searchGene(gene);
        });

        (document.getElementById('gene-input') as HTMLInputElement).addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                const gene = (e.target as HTMLInputElement).value.trim();
                if (gene) this.client.searchGene(gene);
            }
        });

        // Load detail button (for sequences at high zoom)
        document.getElementById('load-detail-btn')!.addEventListener('click', () => {
            this.loadDetailView();
        });

        (document.getElementById('region-input') as HTMLInputElement).addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                const region = (e.target as HTMLInputElement).value.trim();
                if (region) this.client.requestRegion(region);
            }
        });

        // Evidence panel close button
        document.getElementById('close-evidence')!.addEventListener('click', () => {
            this.evidencePanel.classList.remove('visible');
        });

        // Variant filter toggle
        document.getElementById('filter-high')!.addEventListener('click', () => {
            this.state.variantFilter = 'high';
            document.getElementById('filter-high')!.classList.add('active');
            document.getElementById('filter-all')!.classList.remove('active');
            this.renderVariantTable();
        });
        document.getElementById('filter-all')!.addEventListener('click', () => {
            this.state.variantFilter = 'all';
            document.getElementById('filter-all')!.classList.add('active');
            document.getElementById('filter-high')!.classList.remove('active');
            this.renderVariantTable();
        });

        // Variant table header click for sorting
        document.querySelectorAll('#variant-panel th.sortable').forEach(th => {
            th.addEventListener('click', () => {
                const column = (th as HTMLElement).dataset.sort!;
                if (this.state.variantSort.column === column) {
                    this.state.variantSort.direction = this.state.variantSort.direction === 'asc' ? 'desc' : 'asc';
                } else {
                    this.state.variantSort.column = column;
                    this.state.variantSort.direction = 'asc';
                }
                this.renderVariantTable();
            });
        });

        document.getElementById('zoom-in')!.addEventListener('click', () => this.zoom(0.5));
        document.getElementById('zoom-out')!.addEventListener('click', () => this.zoom(2));

        document.getElementById('fullscreen-btn')!.addEventListener('click', () => {
            this.toggleFullscreen();
        });

        document.getElementById('pip-btn')!.addEventListener('click', () => {
            this.togglePip();
        });

        document.getElementById('sync-btn')!.addEventListener('click', () => {
            this.syncToAI();
        });

        // Settings controls
        this.setupSettingsControls();

        // Pan via drag on reads canvas. mouseup clears isDragging before the
        // click event fires, so a drag would end by locking a tooltip; track how
        // far the pointer moved and swallow the click that ends a real drag.
        let isDragging = false;
        let lastX = 0;
        let dragDistance = 0;
        let suppressClick = false;
        const DRAG_CLICK_THRESHOLD_PX = 4;

        this.readsCanvas.addEventListener('mousedown', (e) => {
            isDragging = true;
            lastX = e.clientX;
            dragDistance = 0;
            suppressClick = false;
        });

        window.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const dx = e.clientX - lastX;
            lastX = e.clientX;
            dragDistance += Math.abs(dx);
            if (dx !== 0) this.pan(-dx);
        });

        window.addEventListener('mouseup', () => {
            if (isDragging && dragDistance > DRAG_CLICK_THRESHOLD_PX) suppressClick = true;
            isDragging = false;
        });

        // Tooltips and hover state
        this.readsCanvas.addEventListener('mousemove', (e) => {
            if (isDragging) return;
            const read = this.renderer.getReadAtPosition(e);

            this.showTooltip(e, read);

            // Pair highlighting delay
            if (read !== this.state.pendingHoverRead) {
                this.state.pendingHoverRead = read;

                if (this.hoverDelayTimer) {
                    clearTimeout(this.hoverDelayTimer);
                    this.hoverDelayTimer = null;
                }

                if (read !== this.state.hoveredRead) {
                    this.hoverDelayTimer = setTimeout(() => {
                        this.state.hoveredRead = this.state.pendingHoverRead;
                        this.renderer.render();
                    }, BAMCPViewer.HOVER_DELAY_MS);
                }
            }
        });

        this.readsCanvas.addEventListener('mouseout', () => {
            if (this.hoverDelayTimer) {
                clearTimeout(this.hoverDelayTimer);
                this.hoverDelayTimer = null;
            }
            this.state.pendingHoverRead = null;

            if (!this.state.lockedTooltip) {
                if (this.state.hoveredRead) {
                    this.state.hoveredRead = null;
                    this.renderer.render();
                }
                this.tooltip.style.display = 'none';
            }
        });

        // Click to lock tooltip
        this.readsCanvas.addEventListener('click', (e) => {
            if (isDragging) return;
            if (suppressClick) {
                suppressClick = false;
                return;
            }
            const read = this.renderer.getReadAtPosition(e);
            if (read) {
                this.state.lockedTooltip = { read, x: e.clientX, y: e.clientY };
                this.showTooltipForRead(read, e.clientX, e.clientY);
                this.tooltip.classList.add('locked');

                if (this.hoverDelayTimer) {
                    clearTimeout(this.hoverDelayTimer);
                    this.hoverDelayTimer = null;
                }
                this.state.pendingHoverRead = null;
                this.state.hoveredRead = read;
                this.renderer.render();
            } else {
                this.unlockTooltip();
            }
        });

        // Tooltip drag handling
        this.tooltip.addEventListener('mousedown', (e) => {
            if (!this.state.lockedTooltip) return;
            e.preventDefault();
            e.stopPropagation();

            this.isDraggingTooltip = true;
            this.tooltip.classList.add('dragging');

            const rect = this.tooltip.getBoundingClientRect();
            this.tooltipDragOffset = {
                x: e.clientX - rect.left,
                y: e.clientY - rect.top
            };
        });

        window.addEventListener('mousemove', (e) => {
            if (!this.isDraggingTooltip) return;

            const newX = e.clientX - this.tooltipDragOffset.x;
            const newY = e.clientY - this.tooltipDragOffset.y;

            // Keep tooltip within viewport bounds
            const tooltipRect = this.tooltip.getBoundingClientRect();
            const maxX = window.innerWidth - tooltipRect.width;
            const maxY = window.innerHeight - tooltipRect.height;

            this.tooltip.style.left = Math.max(0, Math.min(newX, maxX)) + 'px';
            this.tooltip.style.top = Math.max(0, Math.min(newY, maxY)) + 'px';
        });

        window.addEventListener('mouseup', () => {
            if (this.isDraggingTooltip) {
                this.isDraggingTooltip = false;
                this.tooltip.classList.remove('dragging');
            }
        });

        // Click outside locked tooltip to dismiss
        document.addEventListener('click', (e) => {
            if (!this.state.lockedTooltip) return;
            if (this.isDraggingTooltip) return;

            // Check if click is inside tooltip or canvas
            const target = e.target as HTMLElement;
            if (this.tooltip.contains(target)) return;
            if (this.readsCanvas.contains(target)) return;

            this.unlockTooltip();
        });

        // Mouse wheel zoom
        this.readsCanvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = e.deltaY > 0 ? 1.2 : 0.8;
            this.zoom(factor);
        });
    }

    private setupKeyboardShortcuts(): void {
        document.addEventListener('keydown', (e) => {
            if (e.target instanceof HTMLInputElement) return;

            switch (e.key) {
                case 'ArrowLeft':
                    this.pan(50);
                    e.preventDefault();
                    break;
                case 'ArrowRight':
                    this.pan(-50);
                    e.preventDefault();
                    break;
                // ... other shortcuts
            }
        });
    }

    private setupSettingsControls(): void {
        // Load saved settings from localStorage
        this.loadSettings();

        // Display mode selector
        const displayModeSelect = document.getElementById('display-mode-select') as HTMLSelectElement;
        displayModeSelect.value = this.state.settings.displayMode;
        displayModeSelect.addEventListener('change', () => {
            this.state.settings.displayMode = displayModeSelect.value as ReadDisplayMode;
            this.state.resortAndRepack();
            this.renderer.resize();
            this.saveSettings();
        });

        // Color-by selector
        const colorBySelect = document.getElementById('color-by-select') as HTMLSelectElement;
        colorBySelect.value = this.state.settings.colorBy;
        colorBySelect.addEventListener('change', () => {
            this.state.settings.colorBy = colorBySelect.value as ColorBy;
            this.renderer.render();
            this.saveSettings();
        });

        // Sort-by selector
        const sortBySelect = document.getElementById('sort-by-select') as HTMLSelectElement;
        sortBySelect.value = this.state.settings.sortBy;
        sortBySelect.addEventListener('change', () => {
            this.state.settings.sortBy = sortBySelect.value as SortBy;
            this.state.resortAndRepack();
            this.renderer.render();
            this.saveSettings();
        });

        // Soft clips toggle
        const showSoftClips = document.getElementById('show-soft-clips') as HTMLInputElement;
        showSoftClips.checked = this.state.settings.showSoftClips;
        showSoftClips.addEventListener('change', () => {
            this.state.settings.showSoftClips = showSoftClips.checked;
            this.renderer.render();
            this.saveSettings();
        });
    }

    private loadSettings(): void {
        try {
            const saved = localStorage.getItem('bamcp-viewer-settings');
            if (saved) {
                const parsed = JSON.parse(saved) as Partial<ViewerSettings>;
                this.state.settings = { ...DEFAULT_SETTINGS, ...parsed };
            }
        } catch {
            // localStorage may be unavailable (e.g. sandboxed iframe, about:blank)
        }
    }

    private saveSettings(): void {
        try {
            localStorage.setItem('bamcp-viewer-settings', JSON.stringify(this.state.settings));
        } catch {
            // localStorage may be unavailable
        }
    }

    private pan(dxPixels: number): void {
        // CSS-pixel width (matches the mouse delta units); readsCanvas.width is the
        // devicePixelRatio-scaled backing store, which would slow panning on HiDPI.
        const width = this.readsCanvas.clientWidth;
        const regionSpan = this.state.viewport.end - this.state.viewport.start;
        const scale = width / regionSpan;
        const bpDelta = dxPixels / scale;
        if (!Number.isFinite(bpDelta)) return;

        this.state.viewport = panViewport(this.state.viewport, bpDelta);
        this.renderer.render();
        this.scheduleContextUpdate();
        this.scheduleVariantTableUpdate();
        this.checkAndFetchViewport();
    }

    private zoom(factor: number): void {
        // Clamped to [MIN_VIEW_SPAN, MAX_VIEW_SPAN] and start >= 0, so a refetch
        // can never ask the server for more than its 1 Mb region limit.
        this.state.viewport = zoomViewport(this.state.viewport, factor);
        this.renderer.render();
        this.scheduleContextUpdate();
        this.scheduleVariantTableUpdate();

        // Check if we need to fetch sequences for base-level view
        this.checkAndRequestSequences();
        this.checkAndFetchViewport();
    }

    /** Debounced in-view highlighting of the variant table during pan/zoom. */
    private scheduleVariantTableUpdate(): void {
        if (this.variantTableTimer) {
            clearTimeout(this.variantTableTimer);
        }
        this.variantTableTimer = setTimeout(() => {
            this.variantTableTimer = null;
            this.updateVariantInView();
        }, 150);
    }

    /**
     * Check if we're zoomed in enough to show bases but missing sequence data.
     * Auto-triggers fetch with file_path for explicit tool invocation.
     * Falls back to button after timeout if LLM doesn't respond.
     */
    private checkAndRequestSequences(): void {
        if (!this.state.data) {
            this.sequenceNotice.classList.add('hidden');
            return;
        }

        // CSS-pixel width so the base-level zoom threshold is compared on the same
        // scale the renderer draws at, not the DPR-scaled backing store.
        const width = this.readsCanvas.clientWidth;
        const viewportSpan = this.state.viewport.end - this.state.viewport.start;
        const scale = width / viewportSpan;

        // Check if we're at base-level zoom
        if (scale < BAMCPViewer.BASE_RENDER_SCALE) {
            this.sequenceNotice.classList.add('hidden');
            return;
        }

        // Check if region is small enough for sequences
        if (viewportSpan > BAMCPViewer.SEQUENCE_REGION_THRESHOLD) {
            this.sequenceNotice.classList.add('hidden');
            return;
        }

        // Check if we're missing sequences (check first read)
        const firstRead = this.state.data.reads[0];
        if (!firstRead || firstRead.sequence) {
            this.sequenceNotice.classList.add('hidden');
            return;
        }

        // Already requesting - don't duplicate
        if (this.pendingSequenceRequest) {
            return;
        }

        // Auto-trigger fetch with loading state
        this.pendingSequenceRequest = true;
        this.sequenceNotice.classList.remove('hidden');
        this.sequenceNotice.classList.add('loading');

        // Build request parameters
        const region = toRegionString(
            this.state.data.contig, this.state.viewport.start, this.state.viewport.end,
        );
        const filePath = this.state.data.file_path;

        // Use direct tool call if file_path available (reliable), else fallback to sendMessage
        if (filePath) {
            // Set fallback timeout - show button if fetch fails
            if (this.sequenceFallbackTimer) {
                clearTimeout(this.sequenceFallbackTimer);
            }
            this.sequenceFallbackTimer = setTimeout(() => {
                if (this.pendingSequenceRequest) {
                    this.sequenceNotice.classList.remove('loading');
                }
            }, BAMCPViewer.SEQUENCE_FETCH_TIMEOUT_MS);

            this.client.fetchRegionDirect(
                filePath, region, this.state.data?.reference,
                this.state.data?.vcf_path, this.state.data?.variant_source,
            ).then(result => {
                this.pendingSequenceRequest = false;
                if (this.sequenceFallbackTimer) {
                    clearTimeout(this.sequenceFallbackTimer);
                    this.sequenceFallbackTimer = null;
                }
                this.sequenceNotice.classList.add('hidden');
                this.sequenceNotice.classList.remove('loading');
                this.applyFetchResult(result);
            });
        } else {
            this.client.requestRegion(region);
            // sendMessage path — response comes via ontoolresult → onDataReceived
            if (this.sequenceFallbackTimer) {
                clearTimeout(this.sequenceFallbackTimer);
            }
            this.sequenceFallbackTimer = setTimeout(() => {
                if (this.pendingSequenceRequest) {
                    this.sequenceNotice.classList.remove('loading');
                }
            }, BAMCPViewer.SEQUENCE_FETCH_TIMEOUT_MS);
        }
    }

    /**
     * Request detailed data for current viewport (with sequences).
     * Used as manual fallback when auto-fetch times out.
     */
    private loadDetailView(): void {
        if (!this.state.data || this.pendingSequenceRequest) return;

        const region = toRegionString(
            this.state.data.contig, this.state.viewport.start, this.state.viewport.end,
        );
        const filePath = this.state.data.file_path;

        this.pendingSequenceRequest = true;
        this.sequenceNotice.classList.add('loading');

        if (filePath) {
            this.client.fetchRegionDirect(
                filePath, region, this.state.data?.reference,
                this.state.data?.vcf_path, this.state.data?.variant_source,
            ).then(result => {
                this.pendingSequenceRequest = false;
                this.sequenceNotice.classList.add('hidden');
                this.sequenceNotice.classList.remove('loading');
                this.applyFetchResult(result);
            });
        } else {
            // sendMessage path — response comes via ontoolresult → onDataReceived
            this.client.requestRegion(region);
        }
    }

    /**
     * Check if viewport has moved beyond the loaded data range and trigger
     * a debounced refetch if needed. Uses trailing-edge debounce so the
     * fetch fires after the user stops panning/zooming.
     */
    private checkAndFetchViewport(): void {
        const data = this.state.data;
        if (!data || !data.file_path) return;

        // Calculate overlap between current viewport and loaded data range
        const vpStart = this.state.viewport.start;
        const vpEnd = this.state.viewport.end;
        const vpSpan = vpEnd - vpStart;

        const overlapStart = Math.max(vpStart, data.start);
        const overlapEnd = Math.min(vpEnd, data.end);
        const overlap = Math.max(0, overlapEnd - overlapStart);
        const overlapRatio = vpSpan > 0 ? overlap / vpSpan : 1;

        // If viewport is mostly within loaded data, no refetch needed
        if (overlapRatio >= BAMCPViewer.VIEWPORT_REFETCH_OVERLAP) return;

        // Check tile cache before fetching — activate instantly if found
        const cached = this.state.store.bestTile(data.contig, vpStart, vpEnd);
        if (cached) {
            this.state.activateTile(cached);
            this.renderer.render();
            this.renderVariantTable();
            return;
        }

        // Trailing-edge debounce: clear previous timer, set new one
        if (this.viewportRefetchTimer) {
            clearTimeout(this.viewportRefetchTimer);
        }

        this.viewportRefetchTimer = setTimeout(async () => {
            this.viewportRefetchTimer = null;

            // Re-validate after delay (user may have panned back into range)
            const d = this.state.data;
            if (!d || !d.file_path) return;
            const vs = this.state.viewport.start;
            const ve = this.state.viewport.end;
            const sp = ve - vs;
            const os = Math.max(vs, d.start);
            const oe = Math.min(ve, d.end);
            const ol = sp > 0 ? Math.max(0, oe - os) / sp : 1;
            if (ol >= BAMCPViewer.VIEWPORT_REFETCH_OVERLAP) return;

            // Re-check tile cache (may have been populated while debouncing)
            const cachedTile = this.state.store.bestTile(d.contig, vs, ve);
            if (cachedTile) {
                this.state.activateTile(cachedTile);
                this.renderer.render();
                this.renderVariantTable();
                return;
            }

            // Pad fetch region by 50% each side to reduce subsequent refetches,
            // capped so the request stays within the server's region limit.
            const win = paddedFetchWindow({ start: vs, end: ve });
            const region = toRegionString(d.contig, win.start, win.end);

            const result = await this.client.fetchRegionDirect(
                d.file_path, region, d.reference, d.vcf_path, d.variant_source,
            );
            this.applyFetchResult(result);
        }, BAMCPViewer.VIEWPORT_REFETCH_DEBOUNCE_MS);
    }

    /**
     * Apply an app-initiated fetch. Stale responses (superseded by a newer request
     * or by host data) are dropped; errors are shown to the user instead of being
     * swallowed; a tile is only activated if it is on the current contig and still
     * overlaps the viewport (the user may have moved on while it was in flight).
     */
    private applyFetchResult(result: FetchResult): void {
        if (result.status === 'stale') return;
        if (result.status === 'error') {
            this.showFetchError(result.message);
            return;
        }
        this.hideFetchError();
        const data = result.data;
        const current = this.state.data;
        const vp = this.state.viewport;
        if (current && data.contig !== current.contig) return;
        if (!(data.end > vp.start && data.start < vp.end)) return;

        this.state.loadTile(data);
        this.renderVariantTable();
        this.renderer.resize();
        this.scheduleContextUpdate();
    }

    private showFetchError(message: string): void {
        const banner = document.getElementById('fetch-error');
        const text = document.getElementById('fetch-error-text');
        if (!banner || !text) return;
        text.textContent = `Could not load this view: ${message}`;
        banner.classList.remove('hidden');
        if (this.fetchErrorTimer) clearTimeout(this.fetchErrorTimer);
        this.fetchErrorTimer = setTimeout(() => this.hideFetchError(), BAMCPViewer.FETCH_ERROR_MS);
    }

    private hideFetchError(): void {
        if (this.fetchErrorTimer) {
            clearTimeout(this.fetchErrorTimer);
            this.fetchErrorTimer = null;
        }
        document.getElementById('fetch-error')?.classList.add('hidden');
    }

    private scheduleContextUpdate(): void {
        if (this._contextUpdateTimer) {
            clearTimeout(this._contextUpdateTimer);
        }
        this._contextUpdateTimer = setTimeout(() => {
            // Update context via client
            if (!this.state.data) return;
            const coverage = this.state.data.coverage;
            const meanCov = coverage.length > 0
                ? (coverage.reduce((a, b) => a + b, 0) / coverage.length).toFixed(1)
                : '0';

            const context = {
                region: toRegionString(
                    this.state.data.contig, this.state.viewport.start, this.state.viewport.end,
                ),
                reads: this.state.data.reads.length,
                meanCoverage: parseFloat(meanCov),
                variantCount: this.state.data.variants?.length ?? 0,
            };
            this.client.updateModelContext(context);
            this.syncRegionInput();
        }, CONTEXT_UPDATE_DEBOUNCE_MS);
    }

    /**
     * Mirror the current viewport into the region box as a 1-based inclusive
     * samtools-style string (what the user would type, and what tools accept).
     * Skipped while the user is editing the box so we never clobber their input.
     */
    private syncRegionInput(): void {
        if (!this.state.data) return;
        const input = document.getElementById('region-input') as HTMLInputElement | null;
        if (!input || document.activeElement === input) return;
        input.value = toRegionString(
            this.state.data.contig, this.state.viewport.start, this.state.viewport.end,
        );
    }

    /** Sync all DOM chrome with the current state (called after every render). */
    private syncUiState(): void {
        this.syncSettingsControls();
        this.syncRegionInput();
        this.updateReferenceBanner();
        this.updateHiddenReadsIndicator();
    }

    /** Make the toolbar selects show the settings actually being rendered. */
    private syncSettingsControls(): void {
        const s = this.state.settings;
        const pairs: Array<[string, string]> = [
            ['display-mode-select', s.displayMode],
            ['color-by-select', s.colorBy],
            ['sort-by-select', s.sortBy],
        ];
        for (const [id, value] of pairs) {
            const el = document.getElementById(id) as HTMLSelectElement | null;
            if (el && el.value !== value && document.activeElement !== el) el.value = value;
        }
        const legend = document.getElementById('legend');
        if (legend && legend.dataset.mode !== s.displayMode) legend.dataset.mode = s.displayMode;
    }

    /** Show/hide the no-reference banner from the payload's applied_filters. */
    private updateReferenceBanner(): void {
        const banner = document.getElementById('reference-banner');
        if (!banner) return;
        const message = referenceBannerMessage(this.state.data);
        const text = document.getElementById('reference-banner-text');
        if (message) {
            if (text && text.textContent !== message) text.textContent = message;
            banner.classList.remove('hidden');
        } else {
            banner.classList.add('hidden');
        }
    }

    /**
     * Count reads that overlap the viewport but sit in packed rows below the
     * visible part of #viewer (or past the canvas cap), and pin an
     * "N reads not shown" bar to the bottom of the visible area.
     */
    private updateHiddenReadsIndicator(): void {
        const el = document.getElementById('hidden-reads-indicator');
        const viewerEl = document.getElementById('viewer');
        if (!el || !viewerEl) return;
        const count = this.hiddenReadsCount(viewerEl);
        if (count <= 0) {
            el.classList.add('hidden');
            el.dataset.count = '0';
            return;
        }
        el.textContent = hiddenReadsLabel(count);
        el.dataset.count = String(count);
        el.style.top = `${viewerEl.scrollTop + viewerEl.clientHeight - HIDDEN_READS_BAR_PX}px`;
        el.classList.remove('hidden');
    }

    private hiddenReadsCount(viewerEl: HTMLElement): number {
        if (!this.state.data || this.state.packedRows.length === 0) return 0;
        const { height, gap } = this.renderer.getReadDimensions();
        // Visible reads area, minus the bar itself (it covers the bottom rows).
        const visibleH = viewerEl.scrollTop + viewerEl.clientHeight
            - this.readsCanvas.offsetTop - HIDDEN_READS_BAR_PX;
        const from = firstHiddenRow(visibleH, height + gap, height, this.renderer.drawnRowCount);
        return hiddenReadCount(this.state.packedRows, this.state.viewport, from);
    }

    private renderVariantTable(): void {
        const variants = this.state.getFilteredAndSortedVariants();
        // variantTable IS the tbody element (id="variant-table")
        this.variantTable.innerHTML = '';

        variants.forEach((v, idx) => {
            const tr = document.createElement('tr');
            const key = variantKey(v);
            tr.dataset.index = idx.toString();
            tr.dataset.key = key;
            tr.dataset.pos = String(v.position);
            tr.className = key === this.state.selectedVariantKey ? 'selected' : '';
            if (v.is_low_confidence) tr.classList.add('low-confidence');

            // VAF color coding: green ≥40%, orange 20-40%, red <20%
            const vafColor = v.vaf >= 0.4 ? '#22c55e' : v.vaf >= 0.2 ? '#f97316' : '#ef4444';

            // Quality color coding: green ≥30, orange 20-30, red <20
            const qual = Number(v.mean_quality) || 0;
            const qualColor = qual >= 30 ? '#22c55e' : qual >= 20 ? '#f97316' : '#ef4444';

            // Strand counts: "5F:3R" format from accumulated evidence
            const evidence = this.state.store.getEvidence(`${v.position}:${v.ref}>${v.alt}`);
            const strandDisplay = escapeHtml(evidence
                ? `${evidence.forward_count}F:${evidence.reverse_count}R`
                : (v.strand_forward != null && v.strand_reverse != null)
                    ? `${v.strand_forward}F:${v.strand_reverse}R`
                    : '-');

            // Structural variants carry a symbolic alt (<DEL>/<DUP>/…); render it as a
            // typed badge and escape all VCF-derived text before injecting it as HTML.
            let altCell: string;
            if (v.variant_kind === 'sv') {
                const svType = svTypeLabel(v.sv_type, v.alt);
                const span = v.sv_end
                    ? `${displayPos(Number(v.position)).toLocaleString()}–${Number(v.sv_end).toLocaleString()}`
                    : '';
                // Badge color matches the canvas SV-span color for that type.
                altCell = `<span class="sv-badge" style="background:${svColor(svType)}" title="${escapeHtml(v.alt)}${span ? ` (${escapeHtml(span)})` : ''}">${escapeHtml(svType)}</span>`;
            } else {
                altCell = `<span style="color:${baseColor(v.alt, '#333')};font-weight:bold">${escapeHtml(v.alt)}</span>`;
            }
            // A multi-sample VCF gets a small sample-count marker on the position cell.
            const sampleCount = Array.isArray(v.sample_names) ? v.sample_names.length : 0;
            const sampleMarker = sampleCount > 1
                ? ` <span class="sample-count" title="${sampleCount} samples genotyped">×${sampleCount}</span>`
                : '';

            tr.innerHTML = `
                <td>${escapeHtml(displayPos(Number(v.position)).toLocaleString())}${sampleMarker}</td>
                <td>${escapeHtml(v.ref)}</td>
                <td>${altCell}</td>
                <td style="color:${vafColor};font-weight:500">${(Number(v.vaf) * 100).toFixed(1)}%</td>
                <td>${escapeHtml(v.depth)}</td>
                <td>${strandDisplay}</td>
                <td style="color:${qualColor}">Q${qual.toFixed(0)}</td>
                <td>
                    <button class="action-btn clinvar">ClinVar</button>
                    <button class="action-btn gnomad">gnomAD</button>
                    <button class="action-btn explain">Explain</button>
                </td>
            `;

            tr.addEventListener('click', (e) => {
                // Handle action buttons
                if ((e.target as HTMLElement).classList.contains('action-btn')) {
                    const btn = e.target as HTMLElement;
                    if (btn.classList.contains('clinvar')) {
                        this.client.lookupClinVar(v);
                    } else if (btn.classList.contains('gnomad')) {
                        this.client.lookupGnomAD(v);
                    } else if (btn.classList.contains('explain')) {
                        this.client.sendVariantMessage(v);
                    }
                    e.stopPropagation();
                    return;
                }

                this.selectVariant(v, key);
            });

            this.variantTable.appendChild(tr);
        });

        this.updateVariantInView();
    }

    /**
     * Refresh the in-view highlighting and badge without rebuilding the rows —
     * this runs on every pan/zoom, where a full rebuild of a large table is costly.
     */
    private updateVariantInView(): void {
        const vpStart = this.state.viewport.start;
        const vpEnd = this.state.viewport.end;
        const rows = this.variantTable.children;
        let inViewCount = 0;
        for (let i = 0; i < rows.length; i++) {
            const tr = rows[i] as HTMLElement;
            const pos = Number(tr.dataset.pos);
            const inView = pos >= vpStart && pos <= vpEnd;
            if (inView) inViewCount++;
            tr.classList.toggle('in-view', inView);
        }
        const totalCount = this.state.store.getAllVariants().length;
        const badge = document.getElementById('variant-count')!;
        badge.textContent = inViewCount < rows.length
            ? `${inViewCount} in view / ${totalCount}`
            : totalCount.toString();
    }

    /** Select a variant row: highlight it, show evidence, and jump the view to it. */
    private selectVariant(v: Variant, key: string): void {
        this.state.selectedVariantKey = key;
        // Mark this variant as the active one for the DeepVariant
        // "supports-variant" channel.
        this.state.settings.activeVariantPosition = v.position;
        this.state.settings.activeVariantAlt = v.alt;
        for (const row of Array.from(this.variantTable.children) as HTMLElement[]) {
            row.classList.toggle('selected', row.dataset.key === key);
        }

        // Show evidence
        this.showVariantEvidence(v);

        // Jump to variant, then load what the new view needs: the data tile (the
        // variant may sit outside the loaded range), sequences at base-level zoom,
        // and the model context — exactly as a pan/zoom would.
        const span = 100;
        const start = Math.max(0, v.position - span / 2);
        this.state.viewport = { start, end: start + span };
        this.renderer.render();
        this.updateVariantInView();
        this.scheduleContextUpdate();
        this.checkAndRequestSequences();
        this.checkAndFetchViewport();
    }

    private showTooltip(e: MouseEvent, read: Read | null): void {
        if (!read || this.state.lockedTooltip) {
            if (!this.state.lockedTooltip) this.tooltip.style.display = 'none';
            return;
        }
        this.showTooltipForRead(read, e.clientX, e.clientY);
    }

    private unlockTooltip(): void {
        this.state.lockedTooltip = null;
        this.tooltip.classList.remove('locked');
        this.tooltip.classList.remove('dragging');
        this.tooltip.style.display = 'none';

        if (this.state.hoveredRead) {
            this.state.hoveredRead = null;
            this.renderer.render();
        }
    }

    private showTooltipForRead(read: Read, x: number, y: number): void {
        this.tooltip.style.display = 'block';

        // MAPQ color coding: green ≥30, orange 10-30, red <10
        const mapqColor = read.mapping_quality >= 30 ? '#22c55e' :
                          read.mapping_quality >= 10 ? '#f97316' : '#ef4444';

        // Format mismatches with ref→alt colors (show first 5)
        let mismatchHtml = '';
        if (read.mismatches.length > 0) {
            const mmList = read.mismatches.slice(0, 5).map(mm =>
                `${escapeHtml(displayPos(Number(mm.pos)).toLocaleString())}: ` +
                `<span style="color:${baseColor(mm.ref, '#666')}">${escapeHtml(mm.ref)}</span>→` +
                `<span style="color:${baseColor(mm.alt, '#666')}">${escapeHtml(mm.alt)}</span>`
            ).join(', ');
            const more = read.mismatches.length > 5 ? ` +${read.mismatches.length - 5} more` : '';
            mismatchHtml = `<br><span style="color:var(--color-text-muted)">Mismatches:</span> ${mmList}${more}`;
        }

        // Paired-end info
        let pairInfo = read.is_paired ? 'Paired' : 'Single';
        if (read.is_paired && read.insert_size) {
            const insertColor = Math.abs(read.insert_size) > 1000 ? 'var(--color-danger)' :
                               Math.abs(read.insert_size) < 100 ? 'var(--color-warning)' : 'var(--color-success)';
            pairInfo += ` | Insert: <span style="color:${insertColor}">${escapeHtml(Math.abs(read.insert_size))}bp</span>`;
        }

        // Every BAM/server-derived value is escaped: a SAM QNAME may legally contain
        // `<`, `>`, `"`, `'` and `&`, so an unescaped name is an XSS vector.
        const pos = displayPos(Number(read.position)).toLocaleString();
        const end = Number(read.end_position).toLocaleString();
        this.tooltip.innerHTML = `
            <strong>${escapeHtml(read.name)}</strong><br>
            Pos: ${escapeHtml(pos)}-${escapeHtml(end)}<br>
            CIGAR: <span style="font-family:monospace">${escapeHtml(read.cigar)}</span><br>
            MAPQ: <span style="color:${mapqColor};font-weight:bold">${escapeHtml(read.mapping_quality)}</span><br>
            ${pairInfo}${mismatchHtml}
        `;

        // Position tooltip using viewport coordinates (since position: fixed)
        // Set initial position to get accurate dimensions
        this.tooltip.style.left = (x + 10) + 'px';
        this.tooltip.style.top = (y + 10) + 'px';

        // Adjust to keep tooltip within viewport
        const tooltipRect = this.tooltip.getBoundingClientRect();
        let left = x + 10;
        let top = y + 10;

        // Keep within right edge
        if (left + tooltipRect.width > window.innerWidth - 10) {
            left = x - tooltipRect.width - 10;
        }
        // Keep within bottom edge
        if (top + tooltipRect.height > window.innerHeight - 10) {
            top = y - tooltipRect.height - 10;
        }
        // Keep within left edge
        if (left < 10) {
            left = 10;
        }
        // Keep within top edge
        if (top < 10) {
            top = 10;
        }

        this.tooltip.style.left = left + 'px';
        this.tooltip.style.top = top + 'px';
    }

    /** Populate the SV-metadata + multi-sample-genotype details block, or hide it. */
    private renderVariantDetails(variant: Variant): void {
        const el = document.getElementById('variant-details');
        if (!el) return;

        const rows: string[] = [];
        if (variant.variant_kind === 'sv') {
            rows.push(
                `<span class="detail-key">Type</span>` +
                `<span class="sv-badge" style="background:${svColor(svTypeLabel(variant.sv_type, variant.alt))}">` +
                `${escapeHtml(svTypeLabel(variant.sv_type, variant.alt))}</span>`
            );
            if (variant.sv_end) {
                const span = `${displayPos(variant.position).toLocaleString()}` +
                    `–${Number(variant.sv_end).toLocaleString()}`;
                rows.push(`<span class="detail-key">Span</span>${escapeHtml(span)}`);
            }
            // sv_len is normalized to this ALT's SVLEN by the backend (load_vcf_variants).
            if (variant.sv_len != null) {
                rows.push(`<span class="detail-key">Length</span>${escapeHtml(String(variant.sv_len))} bp`);
            }
        }

        const genotypes = formatGenotypes(variant.samples);
        if (genotypes.length) {
            const items = genotypes.map(g => `<li>${escapeHtml(g)}</li>`).join('');
            rows.push(`<span class="detail-key">Genotypes</span><ul class="genotype-list">${items}</ul>`);
        }

        if (rows.length) {
            el.innerHTML = rows.map(r => `<div class="detail-row">${r}</div>`).join('');
            el.classList.remove('hidden');
        } else {
            el.innerHTML = '';
            el.classList.add('hidden');
        }
    }

    /** Blank every evidence chart/stat so a no-evidence variant shows nothing stale. */
    private clearEvidenceCharts(): void {
        for (const id of ['strand-chart', 'quality-chart', 'position-chart', 'mapq-chart']) {
            const canvas = document.getElementById(id) as HTMLCanvasElement | null;
            const ctx = canvas?.getContext('2d');
            if (canvas && ctx) ctx.clearRect(0, 0, canvas.width, canvas.height);
        }
        for (const id of ['strand-ratio', 'quality-stats', 'position-stats']) {
            const el = document.getElementById(id);
            if (el) el.textContent = '';
        }
        const warning = document.getElementById('strand-warning');
        if (warning) warning.style.display = 'none';
        const list = document.getElementById('artifact-risk-list');
        if (list) list.innerHTML = '';
    }

    private showVariantEvidence(variant: Variant): void {
        this.evidencePanel.classList.add('visible');
        const evidence = this.state.store.getEvidence(`${variant.position}:${variant.ref}>${variant.alt}`);

        // Update Title
        document.getElementById('evidence-title')!.textContent =
            `Variant Evidence: ${variant.contig}:${displayPos(variant.position)} ${variant.ref}>${variant.alt}`;

        // SV metadata + multi-sample genotypes render for every variant, independent of
        // whether it has BAMCP read-level evidence (SVs/indels don't).
        this.renderVariantDetails(variant);

        if (!evidence) {
            // No BAMCP read-level evidence for this variant (e.g. a pass-through VCF
            // indel/SV). Clear the previous variant's charts so the panel never shows
            // stale evidence for the wrong variant.
            this.clearEvidenceCharts();
            const meter = document.getElementById('artifact-risk-meter');
            if (meter) {
                meter.innerHTML =
                    '<span style="color:var(--color-text-muted)">' +
                    'No BAMCP read-level evidence for this variant.</span>';
            }
            return;
        }

        // Draw Strand Bias Chart
        this.renderStrandChart(evidence);

        // Draw Quality Histogram
        this.renderHistogram(
            'quality-chart',
            evidence.quality_histogram || [],
            ['0-10', '10-20', '20-30', '30-40', '40+'],
            '#3b82f6'
        );
        document.getElementById('quality-stats')!.textContent =
            `Mean: ${evidence.mean_quality.toFixed(1)} | Median: ${evidence.median_quality}`;

        // Draw Position-in-Read Histogram
        this.renderHistogram(
            'position-chart',
            evidence.position_histogram || [],
            ['0-25', '25-50', '50-75', '75-100', '100-150', '150+'],
            '#10b981'
        );
        document.getElementById('position-stats')!.textContent =
            formatNearEndSummary(evidence);

        // Draw MAPQ Histogram
        this.renderHistogram(
            'mapq-chart',
            evidence.mapq_histogram || [],
            ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60+'],
            '#8b5cf6'
        );

        // Render Artifact Risk
        this.renderArtifactRisk(evidence.artifact_risk);
    }

    private renderStrandChart(evidence: VariantEvidence): void {
        const strandCanvas = document.getElementById('strand-chart') as HTMLCanvasElement;
        const sCtx = strandCanvas.getContext('2d')!;
        sCtx.clearRect(0, 0, strandCanvas.width, strandCanvas.height);

        const fwdPct = evidence.forward_count / (evidence.forward_count + evidence.reverse_count || 1);
        const w = strandCanvas.width;
        sCtx.fillStyle = '#60a5fa'; // Forward
        sCtx.fillRect(0, 0, w * fwdPct, 24);
        sCtx.fillStyle = '#a78bfa'; // Reverse
        sCtx.fillRect(w * fwdPct, 0, w * (1 - fwdPct), 24);

        // Add text overlay
        sCtx.fillStyle = '#fff';
        sCtx.font = '10px sans-serif';
        sCtx.textAlign = 'center';
        sCtx.textBaseline = 'middle';
        if (fwdPct > 0.1) sCtx.fillText(evidence.forward_count.toString(), w * fwdPct / 2, 12);
        if (fwdPct < 0.9) sCtx.fillText(evidence.reverse_count.toString(), w * fwdPct + w * (1 - fwdPct) / 2, 12);

        document.getElementById('strand-ratio')!.textContent =
            `Fwd: ${evidence.forward_count} / Rev: ${evidence.reverse_count}`;

        // Show warning if high strand bias
        const warningEl = document.getElementById('strand-warning')!;
        if (evidence.strand_bias > 0.8) {
            warningEl.style.display = 'flex';
        } else {
            warningEl.style.display = 'none';
        }
    }

    private renderHistogram(
        canvasId: string,
        data: number[],
        labels: string[],
        color: string
    ): void {
        const canvas = document.getElementById(canvasId) as HTMLCanvasElement;
        if (!canvas) return;

        const ctx = canvas.getContext('2d')!;
        const width = canvas.width;
        const height = canvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!data || data.length === 0) {
            ctx.fillStyle = '#9ca3af';
            ctx.font = '10px sans-serif';
            ctx.fillText('No data', 10, height / 2);
            return;
        }

        const maxValue = maxOf(data, 1);
        const barWidth = (width - 20) / data.length;
        const barGap = 2;

        // Draw bars
        data.forEach((value, i) => {
            const barHeight = (value / maxValue) * (height - 20);
            const x = 10 + i * barWidth;
            const y = height - 15 - barHeight;

            ctx.fillStyle = color;
            ctx.fillRect(x + barGap / 2, y, barWidth - barGap, barHeight);

            // Value label on bar
            if (value > 0) {
                ctx.fillStyle = '#374151';
                ctx.font = '9px sans-serif';
                ctx.textAlign = 'center';
                ctx.fillText(value.toString(), x + barWidth / 2, y - 2);
            }
        });

        // X-axis labels
        ctx.fillStyle = '#6b7280';
        ctx.font = '8px sans-serif';
        ctx.textAlign = 'center';
        labels.slice(0, data.length).forEach((label, i) => {
            const x = 10 + i * barWidth + barWidth / 2;
            ctx.fillText(label, x, height - 2);
        });
    }

    private renderArtifactRisk(artifactRisk: ArtifactAssessment | undefined): void {
        const meterEl = document.getElementById('artifact-risk-meter')!;
        const listEl = document.getElementById('artifact-risk-list')!;

        if (!artifactRisk) {
            meterEl.innerHTML = '<span style="color:var(--color-text-muted)">No data</span>';
            listEl.innerHTML = '';
            return;
        }

        // Risk meter visualization
        // Clamp to a finite [0, 1] so a malformed score can't inject into the style attr.
        const rawScore = Number(artifactRisk.risk_score);
        const score = Number.isFinite(rawScore) ? Math.min(1, Math.max(0, rawScore)) : 0;
        const likelihood = artifactRisk.artifact_likelihood;
        const color = likelihood === 'high' ? 'var(--color-danger)' :
                      likelihood === 'medium' ? 'var(--color-warning)' : 'var(--color-success)';

        meterEl.innerHTML = `
            <div style="display:flex;align-items:center;gap:8px;">
                <div style="width:100px;height:8px;background:var(--color-border-light);border-radius:4px;overflow:hidden;">
                    <div style="width:${score * 100}%;height:100%;background:${color};"></div>
                </div>
                <span style="color:${color};font-weight:600;text-transform:uppercase;font-size:11px;">
                    ${escapeHtml(likelihood)}
                </span>
            </div>
        `;

        // Risk list
        const risks = Array.isArray(artifactRisk.risks) ? artifactRisk.risks : [];
        if (risks.length > 0) {
            listEl.innerHTML = risks.map(risk => `
                <li style="color:${risk.severity === 'high' ? 'var(--color-danger)' : 'var(--color-warning)'};">
                    ${escapeHtml(risk.description)}
                </li>
            `).join('');
        } else {
            listEl.innerHTML = '<li style="color:var(--color-success);">No concerns identified</li>';
        }
    }

    private async toggleFullscreen(): Promise<void> {
        const btn = document.getElementById('fullscreen-btn')!;
        const container = document.getElementById('container')!;

        // Use MCP Apps display mode API (falls back to browser fullscreen)
        const success = await this.client.toggleFullscreen();

        if (success) {
            this.isFullscreen = this.client.getCurrentDisplayMode() === 'fullscreen';
            btn.classList.toggle('active', this.isFullscreen);
            container.classList.toggle('fullscreen', this.isFullscreen);
        }
    }

    private async togglePip(): Promise<void> {
        const btn = document.getElementById('pip-btn')!;
        const currentMode = this.client.getCurrentDisplayMode();

        // Toggle between pip and inline
        const newMode = currentMode === 'pip' ? 'inline' : 'pip';
        const success = await this.client.requestDisplayMode(newMode);

        if (success) {
            const isPip = this.client.getCurrentDisplayMode() === 'pip';
            btn.classList.toggle('active', isPip);

            // Also update fullscreen button state
            const fsBtn = document.getElementById('fullscreen-btn')!;
            this.isFullscreen = this.client.getCurrentDisplayMode() === 'fullscreen';
            fsBtn.classList.toggle('active', this.isFullscreen);
        }
    }

    private syncToAI(): void {
        if (!this.state.data) return;

        const region = toRegionString(
            this.state.data.contig, this.state.viewport.start, this.state.viewport.end,
        );
        this.client.syncContext(region);
    }
}

// Init — expose the viewer on window for debugging/E2E hooks.
(window as unknown as { viewer: BAMCPViewer }).viewer = new BAMCPViewer();
