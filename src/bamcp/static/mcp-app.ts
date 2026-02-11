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

import { BAMCPClient } from "./client";
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

// Base colors for nucleotides (shared with renderer.ts)
const BASE_COLORS: Record<string, string> = {
    'A': '#22c55e',  // Green
    'T': '#ef4444',  // Red
    'G': '#f97316',  // Orange
    'C': '#3b82f6',  // Blue
    'N': '#9ca3af',  // Gray
};

class BAMCPViewer {
    private client: BAMCPClient;
    private state: StateManager;
    private renderer: Renderer;

    // DOM Elements
    private tooltip: HTMLElement;
    private variantTable: HTMLElement;
    private evidencePanel: HTMLElement;
    private readsCanvas: HTMLCanvasElement;

    // Interaction State
    private _contextUpdateTimer: ReturnType<typeof setTimeout> | null = null;
    private hoverDelayTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly HOVER_DELAY_MS = 200;
    private isFullscreen = false;

    constructor() {
        this.client = new BAMCPClient();
        this.state = new StateManager();
        this.renderer = new Renderer(this.state);

        // Bind DOM elements
        this.tooltip = document.getElementById('tooltip')!;
        this.variantTable = document.getElementById('variant-table')!;
        this.evidencePanel = document.getElementById('evidence-panel')!;
        this.readsCanvas = document.getElementById('reads-canvas') as HTMLCanvasElement;

        // Setup callbacks
        this.client.setOnDataReceived(data => {
            this.state.loadData(data);
            this.renderVariantTable();
            this.renderer.resize();
            this.scheduleContextUpdate();
        });

        this.setupEventListeners();
        this.setupKeyboardShortcuts();

        // Window resize
        window.addEventListener('resize', () => this.renderer.resize());

        // Init client
        this.client.init();

        // Initial render
        this.renderer.resize();
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

        // Settings controls
        this.setupSettingsControls();

        // Pan via drag on reads canvas
        let isDragging = false;
        let lastX = 0;

        this.readsCanvas.addEventListener('mousedown', (e) => {
            isDragging = true;
            lastX = e.clientX;
        });

        window.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const dx = e.clientX - lastX;
            lastX = e.clientX;
            this.pan(-dx);
        });

        window.addEventListener('mouseup', () => isDragging = false);

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
                this.state.lockedTooltip = null;
                this.tooltip.classList.remove('locked');
                this.tooltip.style.display = 'none';

                if (this.state.hoveredRead) {
                    this.state.hoveredRead = null;
                    this.renderer.render();
                }
            }
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
        const saved = localStorage.getItem('bamcp-viewer-settings');
        if (saved) {
            try {
                const parsed = JSON.parse(saved) as Partial<ViewerSettings>;
                this.state.settings = { ...DEFAULT_SETTINGS, ...parsed };
            } catch {
                // Ignore parse errors, use defaults
            }
        }
    }

    private saveSettings(): void {
        localStorage.setItem('bamcp-viewer-settings', JSON.stringify(this.state.settings));
    }

    private pan(dxPixels: number): void {
        const width = this.readsCanvas.width;
        const regionSpan = this.state.viewport.end - this.state.viewport.start;
        const scale = width / regionSpan;
        const bpDelta = dxPixels / scale;

        this.state.viewport.start += bpDelta;
        this.state.viewport.end += bpDelta;
        this.renderer.render();
        this.scheduleContextUpdate();
    }

    private zoom(factor: number): void {
        const span = this.state.viewport.end - this.state.viewport.start;
        const center = (this.state.viewport.start + this.state.viewport.end) / 2;
        const newSpan = span * factor;

        this.state.viewport.start = center - newSpan / 2;
        this.state.viewport.end = center + newSpan / 2;
        this.renderer.render();
        this.scheduleContextUpdate();
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
                region: `${this.state.data.contig}:${this.state.viewport.start.toFixed(0)}-${this.state.viewport.end.toFixed(0)}`,
                reads: this.state.data.reads.length,
                meanCoverage: parseFloat(meanCov),
                variantCount: this.state.data.variants.length,
            };
            this.client.updateModelContext(context);
        }, CONTEXT_UPDATE_DEBOUNCE_MS);
    }

    private renderVariantTable(): void {
        const variants = this.state.getFilteredAndSortedVariants();
        // variantTable IS the tbody element (id="variant-table")
        this.variantTable.innerHTML = '';

        // Update variant count badge
        document.getElementById('variant-count')!.textContent = variants.length.toString();

        variants.forEach((v, idx) => {
            const tr = document.createElement('tr');
            tr.dataset.index = idx.toString();
            tr.className = idx === this.state.selectedVariantIndex ? 'selected' : '';
            if (v.is_low_confidence) tr.classList.add('low-confidence');

            // VAF color coding: green ≥40%, orange 20-40%, red <20%
            const vafColor = v.vaf >= 0.4 ? '#22c55e' : v.vaf >= 0.2 ? '#f97316' : '#ef4444';

            // Quality color coding: green ≥30, orange 20-30, red <20
            const qual = v.mean_quality || 0;
            const qualColor = qual >= 30 ? '#22c55e' : qual >= 20 ? '#f97316' : '#ef4444';

            // Strand counts: "5F:3R" format from evidence if available
            const evidence = this.state.data?.variant_evidence?.[`${v.position}:${v.ref}>${v.alt}`];
            const strandDisplay = evidence
                ? `${evidence.forward_count}F:${evidence.reverse_count}R`
                : (v.alt_count ?? 0) > 0 ? `${v.alt_count}` : '-';

            tr.innerHTML = `
                <td>${v.position.toLocaleString()}</td>
                <td>${v.ref}</td>
                <td style="color:${BASE_COLORS[v.alt] || '#333'};font-weight:bold">${v.alt}</td>
                <td style="color:${vafColor};font-weight:500">${(v.vaf * 100).toFixed(1)}%</td>
                <td>${v.depth}</td>
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
                        // TODO: Implement ClinVar lookup tool call
                    } else if (btn.classList.contains('gnomad')) {
                        // TODO: Implement gnomAD lookup tool call
                    } else if (btn.classList.contains('explain')) {
                        this.client.sendVariantMessage(v);
                    }
                    e.stopPropagation();
                    return;
                }

                this.state.selectedVariantIndex = idx;
                this.renderVariantTable();

                // Show evidence
                this.showVariantEvidence(v);

                // Jump to variant
                const width = this.readsCanvas.width;
                const span = 100;
                this.state.viewport.start = v.position - span / 2;
                this.state.viewport.end = v.position + span / 2;
                this.renderer.render();
            });

            this.variantTable.appendChild(tr);
        });
    }

    private showTooltip(e: MouseEvent, read: Read | null): void {
        if (!read || this.state.lockedTooltip) {
            if (!this.state.lockedTooltip) this.tooltip.style.display = 'none';
            return;
        }
        this.showTooltipForRead(read, e.clientX, e.clientY);
    }

    private showTooltipForRead(read: Read, x: number, y: number): void {
        this.tooltip.style.display = 'block';
        this.tooltip.style.left = (x + 10) + 'px';
        this.tooltip.style.top = (y + 10) + 'px';

        // MAPQ color coding: green ≥30, orange 10-30, red <10
        const mapqColor = read.mapping_quality >= 30 ? '#22c55e' :
                          read.mapping_quality >= 10 ? '#f97316' : '#ef4444';

        // Format mismatches with ref→alt colors (show first 5)
        let mismatchHtml = '';
        if (read.mismatches.length > 0) {
            const mmList = read.mismatches.slice(0, 5).map(mm =>
                `${mm.pos.toLocaleString()}: <span style="color:${BASE_COLORS[mm.ref] || '#666'}">${mm.ref}</span>→` +
                `<span style="color:${BASE_COLORS[mm.alt] || '#666'}">${mm.alt}</span>`
            ).join(', ');
            const more = read.mismatches.length > 5 ? ` +${read.mismatches.length - 5} more` : '';
            mismatchHtml = `<br><span style="color:#9ca3af">Mismatches:</span> ${mmList}${more}`;
        }

        // Paired-end info
        let pairInfo = read.is_paired ? 'Paired' : 'Single';
        if (read.is_paired && read.insert_size) {
            const insertColor = Math.abs(read.insert_size) > 1000 ? '#ef4444' :
                               Math.abs(read.insert_size) < 100 ? '#f97316' : '#22c55e';
            pairInfo += ` | Insert: <span style="color:${insertColor}">${Math.abs(read.insert_size)}bp</span>`;
        }

        this.tooltip.innerHTML = `
            <strong>${read.name}</strong><br>
            Pos: ${read.position.toLocaleString()}-${read.end_position.toLocaleString()}<br>
            CIGAR: <span style="font-family:monospace">${read.cigar}</span><br>
            MAPQ: <span style="color:${mapqColor};font-weight:bold">${read.mapping_quality}</span><br>
            ${pairInfo}${mismatchHtml}
        `;

        // Adjust position to keep tooltip in view
        const containerRect = document.getElementById('container')!.getBoundingClientRect();
        const tooltipRect = this.tooltip.getBoundingClientRect();

        let left = x + 10;
        let top = y + 10;

        if (left + tooltipRect.width > containerRect.right - 10) {
            left = x - tooltipRect.width - 10;
        }
        if (top + tooltipRect.height > containerRect.bottom - 10) {
            top = y - tooltipRect.height - 10;
        }

        this.tooltip.style.left = left + 'px';
        this.tooltip.style.top = top + 'px';
    }

    private showVariantEvidence(variant: Variant): void {
        this.evidencePanel.classList.add('visible');
        const evidence = this.state.data?.variant_evidence?.[`${variant.position}:${variant.ref}>${variant.alt}`];

        // Update Title
        document.getElementById('evidence-title')!.textContent =
            `Variant Evidence: ${variant.contig}:${variant.position} ${variant.ref}>${variant.alt}`;

        if (!evidence) {
            // No detailed evidence available
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
            this.getPositionSummary(evidence.position_histogram);

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

        const maxValue = Math.max(...data, 1);
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

    private getPositionSummary(histogram: number[] | undefined): string {
        if (!histogram || histogram.length === 0) return '';

        const total = histogram.reduce((a, b) => a + b, 0);
        if (total === 0) return '';

        // First two bins (0-25, 25-50) represent near-end positions
        const nearEnd = (histogram[0] || 0) + (histogram[1] || 0);
        const pct = Math.round((nearEnd / total) * 100);

        if (pct > 50) {
            return `Warning: ${pct}% near ends`;
        }
        return `${pct}% near ends`;
    }

    private renderArtifactRisk(artifactRisk: ArtifactAssessment | undefined): void {
        const meterEl = document.getElementById('artifact-risk-meter')!;
        const listEl = document.getElementById('artifact-risk-list')!;

        if (!artifactRisk) {
            meterEl.innerHTML = '<span style="color:#6b7280">No data</span>';
            listEl.innerHTML = '';
            return;
        }

        // Risk meter visualization
        const score = artifactRisk.risk_score;
        const likelihood = artifactRisk.artifact_likelihood;
        const color = likelihood === 'high' ? '#ef4444' :
                      likelihood === 'medium' ? '#f97316' : '#22c55e';

        meterEl.innerHTML = `
            <div style="display:flex;align-items:center;gap:8px;">
                <div style="width:100px;height:8px;background:#e5e7eb;border-radius:4px;overflow:hidden;">
                    <div style="width:${score * 100}%;height:100%;background:${color};"></div>
                </div>
                <span style="color:${color};font-weight:600;text-transform:uppercase;font-size:11px;">
                    ${likelihood}
                </span>
            </div>
        `;

        // Risk list
        if (artifactRisk.risks.length > 0) {
            listEl.innerHTML = artifactRisk.risks.map(risk => `
                <li style="color:${risk.severity === 'high' ? '#ef4444' : '#f97316'};">
                    ${risk.description}
                </li>
            `).join('');
        } else {
            listEl.innerHTML = '<li style="color:#22c55e;">No concerns identified</li>';
        }
    }

    private async toggleFullscreen(): Promise<void> {
        const btn = document.getElementById('fullscreen-btn')!;

        // Use MCP Apps display mode API (falls back to browser fullscreen)
        const success = await this.client.toggleFullscreen();

        if (success) {
            this.isFullscreen = this.client.getCurrentDisplayMode() === 'fullscreen';
            btn.classList.toggle('active', this.isFullscreen);
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
}

// Init
// @ts-ignore
(window as any).viewer = new BAMCPViewer();
