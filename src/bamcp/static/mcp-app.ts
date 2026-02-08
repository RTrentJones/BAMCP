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

import { App } from "@modelcontextprotocol/ext-apps";

// Types for BAM/CRAM data
interface Read {
    name: string;
    sequence: string;
    qualities: number[];
    cigar: string;
    position: number;
    end_position: number;
    mapping_quality: number;
    is_reverse: boolean;
    mismatches: Array<{ pos: number; ref: string; alt: string }>;
    // Paired-end fields
    mate_position?: number | null;
    mate_contig?: string | null;
    insert_size?: number | null;
    is_proper_pair?: boolean;
    is_read1?: boolean;
}

interface Variant {
    contig: string;
    position: number;
    ref: string;
    alt: string;
    vaf: number;
    depth: number;
    // Enhanced fields from tools.py
    strand_forward?: number;
    strand_reverse?: number;
    mean_quality?: number;
    is_low_confidence?: boolean;
}

interface VariantEvidence {
    forward_count: number;
    reverse_count: number;
    strand_bias: number;
    mean_quality: number;
    median_quality: number;
    quality_histogram: Record<string, number>;
    position_histogram: Record<string, number>;
}

interface RegionData {
    contig: string;
    start: number;
    end: number;
    reads: Read[];
    coverage: number[];
    variants: Variant[];
    reference_sequence?: string;
    variant_evidence?: Record<string, VariantEvidence>;
}

interface ToolResultParams {
    structuredContent?: RegionData;
    content?: Array<{ type: string; text?: string }>;
}

interface CigarOp {
    op: string;
    len: number;
}

// Color scheme for nucleotides (IGV-inspired)
const BASE_COLORS: Record<string, string> = {
    'A': '#22c55e',  // Green
    'T': '#ef4444',  // Red
    'G': '#f97316',  // Orange
    'C': '#3b82f6',  // Blue
    'N': '#9ca3af'   // Gray
};

// Track heights
const RULER_HEIGHT = 20;
const REFERENCE_HEIGHT = 24;
const COVERAGE_HEIGHT = 60;
const READ_HEIGHT = 12;
const READ_GAP = 2;
const HEADER_OFFSET = RULER_HEIGHT + REFERENCE_HEIGHT;  // 44px
const CONTEXT_UPDATE_DEBOUNCE_MS = 300;

class BAMCPViewer {
    // Canvas elements
    private rulerCanvas: HTMLCanvasElement;
    private referenceCanvas: HTMLCanvasElement;
    private coverageCanvas: HTMLCanvasElement;
    private readsCanvas: HTMLCanvasElement;

    // Canvas contexts
    private rulerCtx: CanvasRenderingContext2D;
    private referenceCtx: CanvasRenderingContext2D;
    private coverageCtx: CanvasRenderingContext2D;
    private readsCtx: CanvasRenderingContext2D;

    // DOM elements
    private tooltip: HTMLElement;
    private variantTable: HTMLElement;
    private evidencePanel: HTMLElement;
    private strandChart: HTMLCanvasElement;
    private qualityChart: HTMLCanvasElement;
    private positionChart: HTMLCanvasElement;

    // State
    private data: RegionData | null = null;
    private viewport = { start: 0, end: 1000 };
    private packedRows: Read[][] = [];
    private app: App | null = null;
    private isFullscreen = false;
    private _contextUpdateTimer: ReturnType<typeof setTimeout> | null = null;
    private lockedTooltip: { read: Read; x: number; y: number } | null = null;

    // Variant panel state
    private variantFilter: 'high' | 'all' = 'high';
    private variantSort: { column: string; direction: 'asc' | 'desc' } = { column: 'position', direction: 'asc' };
    private selectedVariantIndex: number = -1;
    private expandedVariantIndex: number = -1;

    // Mate pair index: maps read name to its mate (if both are loaded)
    private mateIndex: Map<string, Read> = new Map();
    // Maps read to its row index for Y position lookup
    private readRowIndex: Map<Read, number> = new Map();
    // Currently hovered read (for pair highlighting) - only set after delay or click
    private hoveredRead: Read | null = null;
    // Pending hover (waiting for delay)
    private pendingHoverRead: Read | null = null;
    private hoverDelayTimer: ReturnType<typeof setTimeout> | null = null;
    private static readonly HOVER_DELAY_MS = 600;  // Delay before hover effect kicks in

    constructor() {
        // Initialize canvases
        this.rulerCanvas = document.getElementById('ruler-canvas') as HTMLCanvasElement;
        this.referenceCanvas = document.getElementById('reference-canvas') as HTMLCanvasElement;
        this.coverageCanvas = document.getElementById('coverage-canvas') as HTMLCanvasElement;
        this.readsCanvas = document.getElementById('reads-canvas') as HTMLCanvasElement;

        this.rulerCtx = this.rulerCanvas.getContext('2d')!;
        this.referenceCtx = this.referenceCanvas.getContext('2d')!;
        this.coverageCtx = this.coverageCanvas.getContext('2d')!;
        this.readsCtx = this.readsCanvas.getContext('2d')!;

        this.tooltip = document.getElementById('tooltip')!;
        this.variantTable = document.getElementById('variant-table')!;
        this.evidencePanel = document.getElementById('evidence-panel')!;
        this.strandChart = document.getElementById('strand-chart') as HTMLCanvasElement;
        this.qualityChart = document.getElementById('quality-chart') as HTMLCanvasElement;
        this.positionChart = document.getElementById('position-chart') as HTMLCanvasElement;

        this.setupEventListeners();
        this.setupKeyboardShortcuts();
        this.initApp();
        this.resize();

        window.addEventListener('resize', () => this.resize());
    }

    private async initApp(): Promise<void> {
        try {
            this.app = new App({ name: 'BAMCP Viewer', version: '1.0.0' });

            // Set ontoolresult BEFORE connect() to receive the initial tool result
            this.app.ontoolresult = (params: ToolResultParams) => {
                if (params.structuredContent) {
                    this.loadData(params.structuredContent);
                } else if (params.content) {
                    const text = params.content.find(c => c.type === 'text');
                    if (text?.text) {
                        this.loadData(JSON.parse(text.text));
                    }
                }
            };

            await this.app.connect();
        } catch (e) {
            console.warn('MCP Apps SDK init failed:', e);
        }
    }

    private setupEventListeners(): void {
        document.getElementById('go-btn')!.addEventListener('click', () => {
            const region = (document.getElementById('region-input') as HTMLInputElement).value;
            this.requestRegion(region);
        });

        // Gene search
        document.getElementById('gene-btn')!.addEventListener('click', () => {
            const gene = (document.getElementById('gene-input') as HTMLInputElement).value.trim();
            if (gene) this.searchGene(gene);
        });

        (document.getElementById('gene-input') as HTMLInputElement).addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                const gene = (e.target as HTMLInputElement).value.trim();
                if (gene) this.searchGene(gene);
            }
        });

        // Evidence panel close button
        document.getElementById('close-evidence')!.addEventListener('click', () => {
            this.evidencePanel.classList.remove('visible');
        });

        // Variant filter toggle
        document.getElementById('filter-high')!.addEventListener('click', () => {
            this.setVariantFilter('high');
        });
        document.getElementById('filter-all')!.addEventListener('click', () => {
            this.setVariantFilter('all');
        });

        // Variant table header click for sorting
        document.querySelectorAll('#variant-panel th.sortable').forEach(th => {
            th.addEventListener('click', () => {
                const column = (th as HTMLElement).dataset.sort!;
                this.setVariantSort(column);
            });
        });

        document.getElementById('zoom-in')!.addEventListener('click', () => this.zoom(0.5));
        document.getElementById('zoom-out')!.addEventListener('click', () => this.zoom(2));

        document.getElementById('fullscreen-btn')!.addEventListener('click', () => {
            this.toggleFullscreen();
        });

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

        // Tooltips and hover state for mate highlighting (with delay)
        this.readsCanvas.addEventListener('mousemove', (e) => {
            if (isDragging) return;
            const read = this.getReadAtPosition(e);

            // Tooltip shows immediately (no delay)
            this.showTooltip(e);

            // Pair highlighting has a delay to avoid flickering
            if (read !== this.pendingHoverRead) {
                this.pendingHoverRead = read;

                // Clear existing timer
                if (this.hoverDelayTimer) {
                    clearTimeout(this.hoverDelayTimer);
                    this.hoverDelayTimer = null;
                }

                // Start new timer for hover effect
                if (read !== this.hoveredRead) {
                    this.hoverDelayTimer = setTimeout(() => {
                        this.hoveredRead = this.pendingHoverRead;
                        this.renderReads();
                    }, BAMCPViewer.HOVER_DELAY_MS);
                }
            }
        });

        this.readsCanvas.addEventListener('mouseout', () => {
            // Clear pending hover timer
            if (this.hoverDelayTimer) {
                clearTimeout(this.hoverDelayTimer);
                this.hoverDelayTimer = null;
            }
            this.pendingHoverRead = null;

            // Clear active hover state
            if (this.hoveredRead) {
                this.hoveredRead = null;
                this.renderReads();  // Clear mate highlighting
            }
            if (!this.lockedTooltip) {
                this.tooltip.style.display = 'none';
            }
        });

        // Click to lock/unlock tooltip
        this.readsCanvas.addEventListener('click', (e) => {
            if (isDragging) return;
            const read = this.getReadAtPosition(e);
            if (read) {
                this.lockedTooltip = { read, x: e.clientX, y: e.clientY };
                this.showTooltipForRead(read, e.clientX, e.clientY);
                this.tooltip.classList.add('locked');
            } else {
                this.lockedTooltip = null;
                this.tooltip.classList.remove('locked');
                this.tooltip.style.display = 'none';
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
            // Don't capture if typing in input
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
                case '+':
                case '=':
                    this.zoom(0.5);
                    e.preventDefault();
                    break;
                case '-':
                    this.zoom(2);
                    e.preventDefault();
                    break;
                case 'Home':
                    if (this.data) {
                        this.jumpTo(this.data.start);
                        e.preventDefault();
                    }
                    break;
                case 'End':
                    if (this.data) {
                        this.jumpTo(this.data.end);
                        e.preventDefault();
                    }
                    break;
                // Variant navigation shortcuts
                case 'n':
                case 'N':
                    this.navigateVariant(1);
                    e.preventDefault();
                    break;
                case 'p':
                case 'P':
                    this.navigateVariant(-1);
                    e.preventDefault();
                    break;
                case 'Enter':
                    this.toggleVariantExpansion();
                    e.preventDefault();
                    break;
            }
        });
    }

    private navigateVariant(direction: number): void {
        const variants = this.getFilteredAndSortedVariants();
        if (variants.length === 0) return;

        let newIndex = this.selectedVariantIndex + direction;

        // Wrap around
        if (newIndex < 0) newIndex = variants.length - 1;
        if (newIndex >= variants.length) newIndex = 0;

        this.selectVariant(newIndex);

        const variant = variants[newIndex];
        if (variant) {
            this.jumpTo(variant.position);
            this.renderEvidencePanel(variant);

            // Scroll the variant into view in the table
            const row = this.variantTable.querySelector(`tr[data-index="${newIndex}"]`);
            row?.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
        }
    }

    private toggleVariantExpansion(): void {
        if (this.selectedVariantIndex < 0) return;

        const variants = this.getFilteredAndSortedVariants();
        const variant = variants[this.selectedVariantIndex];

        if (this.expandedVariantIndex === this.selectedVariantIndex) {
            this.expandedVariantIndex = -1;
            this.evidencePanel.classList.remove('visible');
        } else {
            this.expandedVariantIndex = this.selectedVariantIndex;
            if (variant) {
                this.renderEvidencePanel(variant);
            }
        }

        // Update row classes
        this.variantTable.querySelectorAll('tr').forEach(row => {
            const rowIndex = parseInt((row as HTMLElement).dataset.index!);
            row.classList.toggle('expanded', rowIndex === this.expandedVariantIndex);
        });
    }

    private async requestRegion(region: string): Promise<void> {
        if (this.app) {
            try {
                await this.app.sendMessage({
                    role: 'user',
                    content: [{
                        type: 'text',
                        text: `Please browse the genomic region ${region} and show the alignment visualization.`
                    }]
                });
                return;
            } catch (e) {
                console.warn('sendMessage failed:', e);
            }
        }
    }

    private async searchGene(symbol: string): Promise<void> {
        if (this.app) {
            try {
                await this.app.sendMessage({
                    role: 'user',
                    content: [{
                        type: 'text',
                        text: `Search for gene ${symbol} and show the alignment at that location.`
                    }]
                });
            } catch (e) {
                console.warn('Gene search failed:', e);
            }
        }
    }

    loadData(data: RegionData): void {
        this.data = data;
        this.viewport = { start: data.start, end: data.end };

        (document.getElementById('region-input') as HTMLInputElement).value =
            data.contig + ':' + data.start + '-' + data.end;

        this.buildMateIndex();
        this.packReads();
        this.renderVariantTable();
        this.resize();
        this.scheduleContextUpdate();
    }

    private buildMateIndex(): void {
        // Build index of read pairs where both mates are loaded
        this.mateIndex.clear();
        if (!this.data) return;

        // Group reads by name
        const readsByName = new Map<string, Read[]>();
        for (const read of this.data.reads) {
            const existing = readsByName.get(read.name);
            if (existing) {
                existing.push(read);
            } else {
                readsByName.set(read.name, [read]);
            }
        }

        // For pairs where both are loaded, create mate links
        for (const [name, reads] of readsByName) {
            if (reads.length === 2) {
                // Both mates present - link them by position
                const [r1, r2] = reads;
                // Use position as key to find mate
                this.mateIndex.set(`${name}:${r1.position}`, r2);
                this.mateIndex.set(`${name}:${r2.position}`, r1);
            }
        }
    }

    private scheduleContextUpdate(): void {
        if (this._contextUpdateTimer) {
            clearTimeout(this._contextUpdateTimer);
        }
        this._contextUpdateTimer = setTimeout(() => {
            this.updateModelContext();
        }, CONTEXT_UPDATE_DEBOUNCE_MS);
    }

    private async updateModelContext(): Promise<void> {
        if (!this.app || !this.data) return;

        const coverage = this.data.coverage;
        const meanCov = coverage.length > 0
            ? (coverage.reduce((a, b) => a + b, 0) / coverage.length).toFixed(1)
            : '0';

        const context = {
            region: this.data.contig + ':' + this.viewport.start + '-' + this.viewport.end,
            reads: this.data.reads.length,
            meanCoverage: parseFloat(meanCov),
            variantCount: this.data.variants.length,
            variants: this.data.variants.slice(0, 10).map(v => ({
                position: v.contig + ':' + v.position,
                change: v.ref + '>' + v.alt,
                vaf: v.vaf,
                depth: v.depth
            }))
        };

        try {
            await this.app.updateModelContext(context);
        } catch {
            // Ignore errors
        }
    }

    private async sendVariantMessage(variant: Variant): Promise<void> {
        if (!this.app) return;

        const message =
            'Explain the variant at ' + variant.contig + ':' + variant.position +
            ': ' + variant.ref + '>' + variant.alt +
            ', VAF=' + (variant.vaf * 100).toFixed(1) + '%' +
            ', depth=' + variant.depth;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{ type: 'text', text: message }]
            });
        } catch {
            // Ignore errors
        }
    }

    private async toggleFullscreen(): Promise<void> {
        this.isFullscreen = !this.isFullscreen;
        const btn = document.getElementById('fullscreen-btn')!;
        btn.classList.toggle('active', this.isFullscreen);

        const container = document.getElementById('container');
        if (container) {
            try {
                if (this.isFullscreen && !document.fullscreenElement) {
                    await container.requestFullscreen();
                } else if (!this.isFullscreen && document.fullscreenElement) {
                    await document.exitFullscreen();
                }
            } catch {
                // Fullscreen not supported or denied
            }
        }
    }

    private packReads(): void {
        if (!this.data) return;

        this.packedRows = [];
        this.readRowIndex.clear();
        const rowEnds: number[] = [];

        // Group reads by name to identify pairs
        const readsByName = new Map<string, Read[]>();
        for (const read of this.data.reads) {
            const existing = readsByName.get(read.name);
            if (existing) existing.push(read);
            else readsByName.set(read.name, [read]);
        }

        // Sort reads by position, but process pairs together
        const sortedReads: Read[] = [];
        const processed = new Set<string>();
        const allReads = [...this.data.reads].sort((a, b) => a.position - b.position);

        for (const read of allReads) {
            const key = `${read.name}:${read.position}`;
            if (processed.has(key)) continue;

            const pair = readsByName.get(read.name);
            if (pair && pair.length === 2) {
                // Add both reads of the pair consecutively so they pack on same row
                const [r1, r2] = pair.sort((a, b) => a.position - b.position);
                sortedReads.push(r1, r2);
                processed.add(`${r1.name}:${r1.position}`);
                processed.add(`${r2.name}:${r2.position}`);
            } else {
                sortedReads.push(read);
                processed.add(key);
            }
        }

        // Pack into rows (pairs will tend to land on same row due to consecutive processing)
        for (const read of sortedReads) {
            let rowIndex = rowEnds.findIndex(end => end < read.position);

            if (rowIndex === -1) {
                rowIndex = rowEnds.length;
                rowEnds.push(read.end_position);
                this.packedRows.push([]);
            } else {
                rowEnds[rowIndex] = read.end_position;
            }

            this.packedRows[rowIndex].push(read);
            this.readRowIndex.set(read, rowIndex);
        }
    }

    private resize(): void {
        const container = document.getElementById('viewer')!;
        const width = Math.max(container.clientWidth, 300);
        const containerHeight = Math.max(container.clientHeight, HEADER_OFFSET + COVERAGE_HEIGHT + 200);

        // Ruler canvas
        this.rulerCanvas.width = width;
        this.rulerCanvas.height = RULER_HEIGHT;

        // Reference canvas
        this.referenceCanvas.width = width;
        this.referenceCanvas.height = REFERENCE_HEIGHT;

        // Coverage canvas
        this.coverageCanvas.width = width;
        this.coverageCanvas.height = COVERAGE_HEIGHT;

        // Reads canvas - remaining height
        this.readsCanvas.width = width;
        const readsHeight = containerHeight - HEADER_OFFSET - COVERAGE_HEIGHT;
        this.readsCanvas.height = readsHeight;

        // Expand if needed for packed reads
        if (this.packedRows.length > 0) {
            const neededHeight = this.packedRows.length * (READ_HEIGHT + READ_GAP);
            if (neededHeight > this.readsCanvas.height) {
                this.readsCanvas.height = Math.min(neededHeight, 800);
            }
        }

        this.render();
    }

    private getScale(): number {
        return this.readsCanvas.width / (this.viewport.end - this.viewport.start);
    }

    private getZoomLevel(): 'nucleotide' | 'read' | 'overview' {
        const scale = this.getScale();
        if (scale >= 8) return 'nucleotide';
        if (scale >= 0.5) return 'read';
        return 'overview';
    }

    private render(): void {
        this.renderRuler();
        this.renderReference();
        this.renderCoverage();
        this.renderReads();
    }

    // ==================== RULER TRACK ====================

    private renderRuler(): void {
        const ctx = this.rulerCtx;
        const width = this.rulerCanvas.width;
        const height = this.rulerCanvas.height;

        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = '#fafafa';
        ctx.fillRect(0, 0, width, height);

        if (!this.data) return;

        const scale = this.getScale();
        const span = this.viewport.end - this.viewport.start;
        const tickInterval = this.calculateTickInterval(span);

        ctx.fillStyle = '#374151';
        ctx.font = '10px monospace';
        ctx.strokeStyle = '#9ca3af';
        ctx.lineWidth = 1;

        // Draw region label
        ctx.fillStyle = '#6b7280';
        ctx.font = 'bold 10px sans-serif';
        ctx.fillText(this.data.contig, 4, 12);

        ctx.font = '10px monospace';
        ctx.fillStyle = '#374151';

        const startTick = Math.ceil(this.viewport.start / tickInterval) * tickInterval;
        for (let pos = startTick; pos <= this.viewport.end; pos += tickInterval) {
            const x = (pos - this.viewport.start) * scale;

            // Tick mark
            ctx.beginPath();
            ctx.moveTo(x, height - 6);
            ctx.lineTo(x, height);
            ctx.stroke();

            // Label
            const label = this.formatPosition(pos);
            const labelWidth = ctx.measureText(label).width;
            if (x + labelWidth < width - 5) {
                ctx.fillText(label, x + 2, height - 8);
            }
        }

        // Bottom border
        ctx.strokeStyle = '#e5e7eb';
        ctx.beginPath();
        ctx.moveTo(0, height - 0.5);
        ctx.lineTo(width, height - 0.5);
        ctx.stroke();
    }

    private calculateTickInterval(span: number): number {
        const targetTicks = 8;
        const rawInterval = span / targetTicks;
        const magnitude = Math.pow(10, Math.floor(Math.log10(rawInterval)));
        const normalized = rawInterval / magnitude;

        if (normalized <= 1) return magnitude;
        if (normalized <= 2) return 2 * magnitude;
        if (normalized <= 5) return 5 * magnitude;
        return 10 * magnitude;
    }

    private formatPosition(pos: number): string {
        if (pos >= 1000000) return (pos / 1000000).toFixed(1) + 'M';
        if (pos >= 1000) return (pos / 1000).toFixed(1) + 'K';
        return pos.toString();
    }

    // ==================== REFERENCE TRACK ====================

    private renderReference(): void {
        const ctx = this.referenceCtx;
        const width = this.referenceCanvas.width;
        const height = this.referenceCanvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!this.data?.reference_sequence) {
            // No reference available - show placeholder
            ctx.fillStyle = '#f9fafb';
            ctx.fillRect(0, 0, width, height);
            ctx.fillStyle = '#9ca3af';
            ctx.font = '11px sans-serif';
            ctx.fillText('Reference sequence not available', 8, height / 2 + 4);
            return;
        }

        const scale = this.getScale();
        const refSeq = this.data.reference_sequence;

        if (scale >= 8) {
            this.renderReferenceNucleotides(ctx, refSeq, scale, height);
        } else if (scale >= 2) {
            this.renderReferenceBlocks(ctx, refSeq, scale, height);
        } else {
            this.renderReferenceBands(ctx, refSeq, scale, width, height);
        }

        // Bottom border
        ctx.strokeStyle = '#e5e7eb';
        ctx.beginPath();
        ctx.moveTo(0, height - 0.5);
        ctx.lineTo(width, height - 0.5);
        ctx.stroke();
    }

    private renderReferenceNucleotides(
        ctx: CanvasRenderingContext2D,
        refSeq: string,
        scale: number,
        height: number
    ): void {
        ctx.font = 'bold 12px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';

        const startIdx = Math.max(0, Math.floor(this.viewport.start - this.data!.start));
        const endIdx = Math.min(refSeq.length, Math.ceil(this.viewport.end - this.data!.start));

        for (let i = startIdx; i < endIdx; i++) {
            const base = refSeq[i].toUpperCase();
            const x = (this.data!.start + i - this.viewport.start) * scale;

            // Background color
            ctx.fillStyle = BASE_COLORS[base] || '#9ca3af';
            ctx.fillRect(x, 2, Math.max(scale - 1, 1), height - 4);

            // Letter (white on colored background)
            if (scale >= 10) {
                ctx.fillStyle = '#fff';
                ctx.fillText(base, x + scale / 2, height / 2);
            }
        }

        ctx.textAlign = 'start';
        ctx.textBaseline = 'alphabetic';
    }

    private renderReferenceBlocks(
        ctx: CanvasRenderingContext2D,
        refSeq: string,
        scale: number,
        height: number
    ): void {
        const startIdx = Math.max(0, Math.floor(this.viewport.start - this.data!.start));
        const endIdx = Math.min(refSeq.length, Math.ceil(this.viewport.end - this.data!.start));

        for (let i = startIdx; i < endIdx; i++) {
            const base = refSeq[i].toUpperCase();
            const x = (this.data!.start + i - this.viewport.start) * scale;

            ctx.fillStyle = BASE_COLORS[base] || '#9ca3af';
            ctx.fillRect(x, 2, Math.max(scale - 0.5, 1), height - 4);
        }
    }

    private renderReferenceBands(
        ctx: CanvasRenderingContext2D,
        refSeq: string,
        scale: number,
        width: number,
        height: number
    ): void {
        // At very low zoom, aggregate base composition per pixel
        const pixelWidth = 1 / scale;
        const viewStart = Math.max(0, this.viewport.start - this.data!.start);

        for (let px = 0; px < width; px++) {
            const baseStart = Math.floor(viewStart + px * pixelWidth);
            const baseEnd = Math.min(refSeq.length, Math.ceil(viewStart + (px + 1) * pixelWidth));

            if (baseStart >= refSeq.length) break;

            // Count bases in this pixel
            const counts: Record<string, number> = { A: 0, T: 0, G: 0, C: 0 };
            for (let i = baseStart; i < baseEnd; i++) {
                const base = refSeq[i].toUpperCase();
                if (counts[base] !== undefined) counts[base]++;
            }

            // Draw stacked bar for this pixel
            let y = 2;
            const total = baseEnd - baseStart;
            for (const [base, count] of Object.entries(counts)) {
                if (count === 0) continue;
                const h = (count / total) * (height - 4);
                ctx.fillStyle = BASE_COLORS[base];
                ctx.fillRect(px, y, 1, h);
                y += h;
            }
        }
    }

    // ==================== COVERAGE TRACK ====================

    private renderCoverage(): void {
        const ctx = this.coverageCtx;
        const width = this.coverageCanvas.width;
        const height = this.coverageCanvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!this.data) return;

        const coverage = this.data.coverage;
        const maxCov = Math.max(...coverage, 1);
        const scale = this.getScale();
        const lowCoverageThreshold = 10;

        // Draw low coverage warning regions first (background)
        for (let i = 0; i < coverage.length; i++) {
            if (coverage[i] < lowCoverageThreshold) {
                const x = (this.data.start + i - this.viewport.start) * scale;
                if (coverage[i] === 0) {
                    // Zero coverage: bright red stripe
                    ctx.fillStyle = 'rgba(239, 68, 68, 0.4)';
                    ctx.fillRect(x, 0, Math.max(scale, 1), height);
                    // Hatching pattern for emphasis
                    ctx.strokeStyle = 'rgba(239, 68, 68, 0.6)';
                    ctx.lineWidth = 1;
                    for (let ly = 0; ly < height; ly += 4) {
                        ctx.beginPath();
                        ctx.moveTo(x, ly);
                        ctx.lineTo(x + Math.max(scale, 2), ly + 4);
                        ctx.stroke();
                    }
                } else {
                    // Low but non-zero coverage: lighter warning
                    ctx.fillStyle = 'rgba(239, 68, 68, 0.2)';
                    ctx.fillRect(x, 0, Math.max(scale, 1), height);
                }
            }
        }

        // Fill area
        ctx.fillStyle = '#93c5fd';
        ctx.strokeStyle = '#3b82f6';
        ctx.lineWidth = 1;
        ctx.beginPath();

        for (let i = 0; i < coverage.length; i++) {
            const x = (this.data.start + i - this.viewport.start) * scale;
            const h = (coverage[i] / maxCov) * (height - 20);

            if (i === 0) {
                ctx.moveTo(x, height - h);
            } else {
                ctx.lineTo(x, height - h);
            }
        }

        ctx.lineTo((this.data.end - this.viewport.start) * scale, height);
        ctx.lineTo(0, height);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();

        // Draw threshold line (dashed red)
        const thresholdY = height - (lowCoverageThreshold / maxCov) * (height - 20);
        if (thresholdY > 15 && thresholdY < height - 5) {
            ctx.strokeStyle = '#ef4444';
            ctx.setLineDash([4, 2]);
            ctx.beginPath();
            ctx.moveTo(0, thresholdY);
            ctx.lineTo(width, thresholdY);
            ctx.stroke();
            ctx.setLineDash([]);

            // Threshold label
            ctx.fillStyle = '#ef4444';
            ctx.font = '9px sans-serif';
            ctx.fillText(lowCoverageThreshold + 'x', width - 25, thresholdY - 2);
        }

        // Max coverage label
        ctx.fillStyle = '#374151';
        ctx.font = 'bold 10px sans-serif';
        ctx.fillText('Max: ' + maxCov + 'x', 5, 12);

        // Mean coverage
        const meanCov = coverage.length > 0
            ? (coverage.reduce((a, b) => a + b, 0) / coverage.length).toFixed(1)
            : '0';
        ctx.font = '10px sans-serif';
        ctx.fillStyle = '#6b7280';
        ctx.fillText('Mean: ' + meanCov + 'x', 80, 12);
    }

    // ==================== READS TRACK ====================

    private renderReads(): void {
        const ctx = this.readsCtx;
        const width = this.readsCanvas.width;
        const height = this.readsCanvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!this.data) return;

        const zoomLevel = this.getZoomLevel();

        // At very low zoom, show overview message
        if (zoomLevel === 'overview') {
            ctx.fillStyle = '#f9fafb';
            ctx.fillRect(0, 0, width, height);
            ctx.fillStyle = '#6b7280';
            ctx.font = '12px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(
                `${this.data.reads.length} reads - zoom in to view`,
                width / 2,
                height / 2
            );
            ctx.textAlign = 'start';
            return;
        }

        const scale = this.getScale();

        // Determine which reads should show connectors:
        // 1. Hovered pair (if any)
        // 2. SV candidates (abnormal pairs)
        const showConnectorFor = new Set<Read>();

        if (this.hoveredRead) {
            showConnectorFor.add(this.hoveredRead);
            const mateKey = `${this.hoveredRead.name}:${this.hoveredRead.position}`;
            const mate = this.mateIndex.get(mateKey);
            if (mate) showConnectorFor.add(mate);
        }

        // Add all SV candidates
        for (const row of this.packedRows) {
            for (const read of row) {
                if (this.isSVCandidate(read)) {
                    showConnectorFor.add(read);
                }
            }
        }

        // Get the mate of hovered read (for highlighting)
        const hoveredMate = this.hoveredRead
            ? this.mateIndex.get(`${this.hoveredRead.name}:${this.hoveredRead.position}`)
            : null;

        // First pass: render mate connections (so they appear behind reads)
        if (scale >= 0.5) {
            ctx.save();
            for (let rowIdx = 0; rowIdx < this.packedRows.length; rowIdx++) {
                const y = rowIdx * (READ_HEIGHT + READ_GAP);
                if (y > height) break;

                for (const read of this.packedRows[rowIdx]) {
                    // Only render connection from read1 to avoid duplicates
                    if (read.is_read1 && showConnectorFor.has(read)) {
                        const isHovered = read === this.hoveredRead || read === hoveredMate;
                        this.renderMateConnection(ctx, read, y, scale, isHovered);
                    }
                }
            }
            ctx.restore();
        }

        // Second pass: render reads on top (dim non-hovered when hovering)
        for (let rowIdx = 0; rowIdx < this.packedRows.length; rowIdx++) {
            const y = rowIdx * (READ_HEIGHT + READ_GAP);

            if (y > height) break;

            for (const read of this.packedRows[rowIdx]) {
                const isDimmed = this.hoveredRead !== null &&
                    read !== this.hoveredRead &&
                    read !== hoveredMate;

                this.renderRead(ctx, read, y, scale, zoomLevel, isDimmed);
            }
        }
    }

    private renderRead(
        ctx: CanvasRenderingContext2D,
        read: Read,
        y: number,
        scale: number,
        zoomLevel: 'nucleotide' | 'read',
        isDimmed: boolean = false
    ): void {
        const x = (read.position - this.viewport.start) * scale;
        const w = (read.end_position - read.position) * scale;

        if (x + w < 0 || x > this.readsCanvas.width) return;

        // MAPQ-based opacity (MAPQ 60 = full, 0 = 30%), with dimming factor
        // When dimmed, reduce opacity to 15% for clear visual distinction
        const baseDim = isDimmed ? 0.15 : 1.0;
        const opacity = (0.3 + Math.min(read.mapping_quality, 60) / 60 * 0.7) * baseDim;

        // Parse CIGAR for indels
        const cigarOps = this.parseCigar(read.cigar);

        // Draw soft-clipped regions first
        const softClips = this.getSoftClipsFromCigar(cigarOps);
        if (softClips.left > 0) {
            const clipW = softClips.left * scale;
            ctx.fillStyle = read.is_reverse
                ? `rgba(167, 139, 250, ${opacity * 0.4})`
                : `rgba(96, 165, 250, ${opacity * 0.4})`;
            ctx.fillRect(x - clipW, y, clipW, READ_HEIGHT);
            ctx.strokeStyle = '#f59e0b';
            ctx.lineWidth = 1;
            ctx.strokeRect(x - clipW, y, clipW, READ_HEIGHT);
        }
        if (softClips.right > 0) {
            const clipW = softClips.right * scale;
            ctx.fillStyle = read.is_reverse
                ? `rgba(167, 139, 250, ${opacity * 0.4})`
                : `rgba(96, 165, 250, ${opacity * 0.4})`;
            ctx.fillRect(x + w, y, clipW, READ_HEIGHT);
            ctx.strokeStyle = '#f59e0b';
            ctx.lineWidth = 1;
            ctx.strokeRect(x + w, y, clipW, READ_HEIGHT);
        }

        // Draw read body with CIGAR operations
        this.renderReadWithCigar(ctx, read, cigarOps, x, y, scale, opacity, zoomLevel);

        // Direction arrow
        if (w > 15) {
            ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
            ctx.font = '8px sans-serif';
            const arrow = read.is_reverse ? '\u25C0' : '\u25B6';
            ctx.fillText(arrow, read.is_reverse ? x + 2 : x + w - 10, y + 9);
        }
    }

    private renderReadWithCigar(
        ctx: CanvasRenderingContext2D,
        read: Read,
        ops: CigarOp[],
        startX: number,
        y: number,
        scale: number,
        opacity: number,
        zoomLevel: 'nucleotide' | 'read'
    ): void {
        let refPos = read.position;
        let queryPos = 0;

        const baseColor = read.is_reverse ? [167, 139, 250] : [96, 165, 250];

        for (const { op, len } of ops) {
            const x = (refPos - this.viewport.start) * scale;

            switch (op) {
                case 'M':
                case '=':
                case 'X':
                    // Build mismatch lookup for this segment
                    const mismatchMap = new Map<number, { ref: string; alt: string }>();
                    for (const mm of read.mismatches) {
                        if (mm.pos >= refPos && mm.pos < refPos + len) {
                            mismatchMap.set(mm.pos, mm);
                        }
                    }

                    if (zoomLevel === 'nucleotide' && scale >= 8) {
                        // NUCLEOTIDE ZOOM: Show individual bases with quality shading
                        const fontSize = Math.min(scale - 2, 11);
                        ctx.font = `bold ${fontSize}px monospace`;
                        ctx.textAlign = 'center';

                        for (let i = 0; i < len; i++) {
                            const pos = refPos + i;
                            const bx = (pos - this.viewport.start) * scale;

                            // Skip if off-screen
                            if (bx < -scale || bx > this.readsCanvas.width + scale) continue;

                            const base = read.sequence[queryPos + i] || '?';
                            const qual = read.qualities[queryPos + i] || 0;
                            const mismatch = mismatchMap.get(pos);

                            // Quality-based opacity: Q30+ = full, Q10 = 60%, Q0 = 30%
                            const qualAlpha = 0.3 + Math.min(qual, 40) / 40 * 0.7;

                            if (mismatch) {
                                // MISMATCH: Colored background with white text
                                ctx.fillStyle = BASE_COLORS[base] || '#ef4444';
                                ctx.fillRect(bx + 1, y, scale - 2, READ_HEIGHT);

                                // White border for emphasis
                                ctx.strokeStyle = 'rgba(255, 255, 255, 0.8)';
                                ctx.lineWidth = 1.5;
                                ctx.strokeRect(bx + 1, y, scale - 2, READ_HEIGHT);

                                // White letter
                                ctx.fillStyle = '#fff';
                                ctx.fillText(base, bx + scale / 2, y + READ_HEIGHT - 2);
                            } else {
                                // MATCH: Gray background with quality-faded letter
                                // Subtle gray background based on quality
                                ctx.fillStyle = `rgba(${baseColor[0]}, ${baseColor[1]}, ${baseColor[2]}, ${opacity * 0.6})`;
                                ctx.fillRect(bx, y, scale, READ_HEIGHT);

                                // Base letter with quality-based opacity
                                ctx.fillStyle = `rgba(50, 50, 50, ${qualAlpha})`;
                                ctx.fillText(base, bx + scale / 2, y + READ_HEIGHT - 2);
                            }
                        }
                        ctx.textAlign = 'start';
                    } else {
                        // READ/OVERVIEW ZOOM: Solid bar with mismatch markers
                        ctx.fillStyle = `rgba(${baseColor[0]}, ${baseColor[1]}, ${baseColor[2]}, ${opacity})`;
                        ctx.fillRect(x, y, len * scale, READ_HEIGHT);

                        // Draw mismatch markers
                        for (const mm of read.mismatches) {
                            if (mm.pos >= refPos && mm.pos < refPos + len) {
                                const mx = (mm.pos - this.viewport.start) * scale;

                                if (scale > 2) {
                                    // Medium zoom: colored rectangle
                                    ctx.fillStyle = BASE_COLORS[mm.alt] || '#9ca3af';
                                    ctx.fillRect(mx, y, Math.max(scale, 2), READ_HEIGHT);
                                } else {
                                    // Low zoom: prominent vertical line with circle marker
                                    const markerColor = BASE_COLORS[mm.alt] || '#ef4444';

                                    ctx.fillStyle = markerColor;
                                    ctx.fillRect(mx, y, Math.max(2, scale), READ_HEIGHT);

                                    // Vertical line extending above read
                                    ctx.strokeStyle = markerColor;
                                    ctx.lineWidth = 2;
                                    ctx.beginPath();
                                    ctx.moveTo(mx + 1, y - 8);
                                    ctx.lineTo(mx + 1, y + READ_HEIGHT);
                                    ctx.stroke();
                                    ctx.lineWidth = 1;

                                    // Circle marker at top
                                    ctx.fillStyle = markerColor;
                                    ctx.beginPath();
                                    ctx.arc(mx + 1, y - 5, 3, 0, Math.PI * 2);
                                    ctx.fill();

                                    // White outline for visibility
                                    ctx.strokeStyle = '#fff';
                                    ctx.lineWidth = 1;
                                    ctx.stroke();
                                }
                            }
                        }
                    }

                    refPos += len;
                    queryPos += len;
                    break;

                case 'I':
                    // Insertion - draw purple marker
                    ctx.fillStyle = '#9333ea';
                    ctx.beginPath();
                    ctx.moveTo(x - 3, y);
                    ctx.lineTo(x + 3, y);
                    ctx.lineTo(x, y + 5);
                    ctx.closePath();
                    ctx.fill();

                    // Vertical line
                    ctx.strokeStyle = '#9333ea';
                    ctx.lineWidth = 2;
                    ctx.beginPath();
                    ctx.moveTo(x, y);
                    ctx.lineTo(x, y + READ_HEIGHT);
                    ctx.stroke();
                    ctx.lineWidth = 1;

                    queryPos += len;
                    break;

                case 'D':
                case 'N':
                    // Deletion - draw thin dark line
                    ctx.fillStyle = '#1f2937';
                    ctx.fillRect(x, y + READ_HEIGHT / 2 - 1, len * scale, 2);
                    refPos += len;
                    break;

                case 'S':
                    // Soft-clip - already handled above
                    queryPos += len;
                    break;

                case 'H':
                    // Hard-clip - nothing to draw
                    break;
            }
        }
    }

    private parseCigar(cigar: string): CigarOp[] {
        const ops: CigarOp[] = [];
        if (!cigar) return ops;

        const regex = /(\d+)([MIDNSHP=X])/g;
        let match;
        while ((match = regex.exec(cigar)) !== null) {
            ops.push({ len: parseInt(match[1]), op: match[2] });
        }
        return ops;
    }

    private getSoftClipsFromCigar(ops: CigarOp[]): { left: number; right: number } {
        let left = 0, right = 0;
        if (ops.length === 0) return { left, right };
        if (ops[0].op === 'S') left = ops[0].len;
        if (ops.length > 1 && ops[ops.length - 1].op === 'S') right = ops[ops.length - 1].len;
        return { left, right };
    }

    // ==================== NAVIGATION ====================

    private zoom(factor: number): void {
        const center = (this.viewport.start + this.viewport.end) / 2;
        const span = (this.viewport.end - this.viewport.start) * factor;

        // Limit zoom range
        const minSpan = 10;
        const maxSpan = this.data ? (this.data.end - this.data.start) * 10 : 100000;
        const clampedSpan = Math.max(minSpan, Math.min(maxSpan, span));

        this.viewport.start = Math.floor(center - clampedSpan / 2);
        this.viewport.end = Math.ceil(center + clampedSpan / 2);

        this.render();
        this.scheduleContextUpdate();
    }

    private pan(pixels: number): void {
        const scale = this.getScale();
        const basePan = Math.round(pixels / scale);

        this.viewport.start += basePan;
        this.viewport.end += basePan;

        this.render();
        this.scheduleContextUpdate();
    }

    private jumpTo(position: number): void {
        const span = this.viewport.end - this.viewport.start;
        this.viewport.start = position - span / 2;
        this.viewport.end = position + span / 2;
        this.render();
        this.scheduleContextUpdate();
    }

    // ==================== TOOLTIP ====================

    private getReadAtPosition(e: MouseEvent): Read | null {
        if (!this.data) return null;

        const rect = this.readsCanvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;

        const scale = this.getScale();
        const genomePos = Math.floor(this.viewport.start + x / scale);
        const rowIdx = Math.floor(y / (READ_HEIGHT + READ_GAP));

        if (rowIdx >= this.packedRows.length) return null;

        return this.packedRows[rowIdx].find(r =>
            r.position <= genomePos && r.end_position >= genomePos
        ) || null;
    }

    private showTooltipForRead(read: Read, clientX: number, clientY: number): void {
        const mapqColor = read.mapping_quality >= 30 ? '#22c55e' :
                          read.mapping_quality >= 10 ? '#f97316' : '#ef4444';

        let html =
            '<strong>' + read.name + '</strong><br>' +
            'Position: ' + read.position.toLocaleString() + '-' + read.end_position.toLocaleString() + '<br>' +
            'MAPQ: <span style="color:' + mapqColor + '">' + read.mapping_quality + '</span><br>' +
            'CIGAR: <code>' + read.cigar + '</code><br>' +
            'Strand: ' + (read.is_reverse ? '← Reverse' : '→ Forward') + '<br>' +
            'Length: ' + read.sequence.length + ' bp';

        // Add mismatches if present
        if (read.mismatches && read.mismatches.length > 0) {
            const mmList = read.mismatches
                .slice(0, 5)  // Limit to first 5
                .map(mm => `${mm.pos.toLocaleString()} <span style="color:${BASE_COLORS[mm.ref] || '#888'}">${mm.ref}</span>→<span style="color:${BASE_COLORS[mm.alt] || '#ef4444'}">${mm.alt}</span>`)
                .join(', ');
            html += '<br><span style="color:#f97316">Mismatches:</span> ' + mmList;
            if (read.mismatches.length > 5) {
                html += ` <span style="color:#9ca3af">+${read.mismatches.length - 5} more</span>`;
            }
        }

        this.tooltip.innerHTML = html;
        this.tooltip.style.display = 'block';
        this.tooltip.style.left = (clientX + 10) + 'px';
        this.tooltip.style.top = (clientY + 10) + 'px';
    }

    private showTooltip(e: MouseEvent): void {
        // Don't update tooltip if locked
        if (this.lockedTooltip) return;
        if (!this.data) return;

        const read = this.getReadAtPosition(e);

        if (read) {
            this.showTooltipForRead(read, e.clientX, e.clientY);
        } else {
            this.tooltip.style.display = 'none';
        }
    }

    // ==================== VARIANT TABLE ====================

    private setVariantFilter(filter: 'high' | 'all'): void {
        this.variantFilter = filter;

        // Update toggle button styles
        document.getElementById('filter-high')!.classList.toggle('active', filter === 'high');
        document.getElementById('filter-all')!.classList.toggle('active', filter === 'all');

        this.renderVariantTable();
    }

    private setVariantSort(column: string): void {
        if (this.variantSort.column === column) {
            // Toggle direction
            this.variantSort.direction = this.variantSort.direction === 'asc' ? 'desc' : 'asc';
        } else {
            this.variantSort.column = column;
            this.variantSort.direction = 'asc';
        }

        // Update header styles
        document.querySelectorAll('#variant-panel th.sortable').forEach(th => {
            th.classList.remove('sorted-asc', 'sorted-desc');
            if ((th as HTMLElement).dataset.sort === this.variantSort.column) {
                th.classList.add(this.variantSort.direction === 'asc' ? 'sorted-asc' : 'sorted-desc');
            }
        });

        this.renderVariantTable();
    }

    private getFilteredAndSortedVariants(): Variant[] {
        if (!this.data) return [];

        let variants = [...this.data.variants];

        // Filter by confidence
        if (this.variantFilter === 'high') {
            variants = variants.filter(v => !v.is_low_confidence);
        }

        // Sort
        const { column, direction } = this.variantSort;
        const multiplier = direction === 'asc' ? 1 : -1;

        variants.sort((a, b) => {
            let cmp = 0;
            switch (column) {
                case 'position':
                    cmp = a.position - b.position;
                    break;
                case 'vaf':
                    cmp = a.vaf - b.vaf;
                    break;
                case 'depth':
                    cmp = a.depth - b.depth;
                    break;
                case 'quality':
                    cmp = (a.mean_quality || 0) - (b.mean_quality || 0);
                    break;
            }
            return cmp * multiplier;
        });

        return variants;
    }

    private renderVariantTable(): void {
        if (!this.data) return;

        const variants = this.getFilteredAndSortedVariants();

        // Update count badge
        const countBadge = document.getElementById('variant-count');
        if (countBadge) {
            const totalCount = this.data.variants.length;
            const highConfCount = this.data.variants.filter(v => !v.is_low_confidence).length;
            countBadge.textContent = this.variantFilter === 'high'
                ? `${highConfCount}`
                : `${totalCount}`;
        }

        this.variantTable.innerHTML = variants.map((v, index) => {
            const vafColor = v.vaf >= 0.4 ? '#22c55e' : v.vaf >= 0.2 ? '#f97316' : '#ef4444';
            const qualColor = (v.mean_quality || 0) >= 30 ? '#22c55e' :
                              (v.mean_quality || 0) >= 20 ? '#f97316' : '#ef4444';
            const strandF = v.strand_forward ?? 0;
            const strandR = v.strand_reverse ?? 0;
            const classes = [
                v.is_low_confidence ? 'low-confidence' : '',
                index === this.selectedVariantIndex ? 'selected' : '',
                index === this.expandedVariantIndex ? 'expanded' : ''
            ].filter(Boolean).join(' ');

            return `<tr class="${classes}" data-index="${index}" data-pos="${v.position}"
                data-ref="${v.ref}" data-alt="${v.alt}"
                data-vaf="${v.vaf}" data-depth="${v.depth}"
                data-contig="${v.contig}">
                <td>${v.contig}:${v.position.toLocaleString()}</td>
                <td style="font-family:monospace">${v.ref}</td>
                <td style="color:${BASE_COLORS[v.alt] || '#333'};font-family:monospace;font-weight:bold">${v.alt}</td>
                <td style="color:${vafColor}">${(v.vaf * 100).toFixed(1)}%</td>
                <td>${v.depth}</td>
                <td style="font-size:10px;color:#6b7280">${strandF}F:${strandR}R</td>
                <td style="color:${qualColor}">Q${v.mean_quality?.toFixed(0) || '?'}</td>
                <td>
                    <button class="action-btn clinvar" data-action="clinvar" title="Look up in ClinVar">C</button>
                    <button class="action-btn gnomad" data-action="gnomad" title="Look up in gnomAD">g</button>
                </td>
            </tr>`;
        }).join('');

        // Attach row click handlers
        this.variantTable.querySelectorAll('tr').forEach(row => {
            row.addEventListener('click', (e) => {
                const target = e.target as HTMLElement;
                const el = row as HTMLElement;
                const index = parseInt(el.dataset.index!);
                const pos = parseInt(el.dataset.pos!);

                // Handle action button clicks
                if (target.classList.contains('action-btn')) {
                    e.stopPropagation();
                    const action = target.dataset.action;
                    const variant = variants[index];
                    if (action === 'clinvar') {
                        this.lookupClinVar(variant);
                    } else if (action === 'gnomad') {
                        this.lookupGnomAD(variant);
                    }
                    return;
                }

                // Handle row selection
                this.selectVariant(index);
                this.jumpTo(pos);

                const variant = variants[index];
                this.sendVariantMessage(variant);
                this.renderEvidencePanel(variant);
            });
        });
    }

    private selectVariant(index: number): void {
        // Toggle expansion if clicking same row
        if (this.selectedVariantIndex === index) {
            this.expandedVariantIndex = this.expandedVariantIndex === index ? -1 : index;
        } else {
            this.selectedVariantIndex = index;
            this.expandedVariantIndex = -1;
        }

        // Update row classes
        this.variantTable.querySelectorAll('tr').forEach(row => {
            const rowIndex = parseInt((row as HTMLElement).dataset.index!);
            row.classList.toggle('selected', rowIndex === this.selectedVariantIndex);
            row.classList.toggle('expanded', rowIndex === this.expandedVariantIndex);
        });
    }

    private async lookupClinVar(variant: Variant): Promise<void> {
        if (!this.app) return;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{
                    type: 'text',
                    text: `Look up the variant ${variant.contig}:${variant.position} ${variant.ref}>${variant.alt} in ClinVar.`
                }]
            });
        } catch (e) {
            console.warn('ClinVar lookup failed:', e);
        }
    }

    private async lookupGnomAD(variant: Variant): Promise<void> {
        if (!this.app) return;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{
                    type: 'text',
                    text: `Look up the variant ${variant.contig}:${variant.position} ${variant.ref}>${variant.alt} in gnomAD for population frequency.`
                }]
            });
        } catch (e) {
            console.warn('gnomAD lookup failed:', e);
        }
    }

    // ==================== EVIDENCE PANEL ====================

    private renderEvidencePanel(variant: Variant): void {
        const key = `${variant.position}:${variant.ref}>${variant.alt}`;
        const evidence = this.data?.variant_evidence?.[key];

        if (!evidence) {
            this.evidencePanel.classList.remove('visible');
            return;
        }

        // Update title
        document.getElementById('evidence-title')!.textContent =
            `${variant.contig}:${variant.position.toLocaleString()} ${variant.ref}>${variant.alt}`;

        // Render strand bias chart
        this.renderStrandChart(evidence);

        // Render quality histogram
        this.renderQualityChart(evidence);

        // Render position histogram
        this.renderPositionChart(evidence);

        // Update stats and warnings
        const strandRatio = document.getElementById('strand-ratio')!;
        strandRatio.textContent = `${evidence.forward_count}F / ${evidence.reverse_count}R`;

        const strandWarning = document.getElementById('strand-warning')!;
        strandWarning.style.display = evidence.strand_bias > 0.8 ? 'block' : 'none';

        const qualityStats = document.getElementById('quality-stats')!;
        qualityStats.textContent = `Mean: Q${evidence.mean_quality} | Median: Q${evidence.median_quality}`;

        // Show the panel
        this.evidencePanel.classList.add('visible');
    }

    private renderStrandChart(evidence: VariantEvidence): void {
        const ctx = this.strandChart.getContext('2d')!;
        const width = this.strandChart.width;
        const height = this.strandChart.height;

        ctx.clearRect(0, 0, width, height);

        const total = evidence.forward_count + evidence.reverse_count;
        if (total === 0) return;

        const forwardWidth = (evidence.forward_count / total) * width;

        // Forward (blue/green)
        ctx.fillStyle = '#60a5fa';
        ctx.fillRect(0, 0, forwardWidth, height);

        // Reverse (purple)
        ctx.fillStyle = '#a78bfa';
        ctx.fillRect(forwardWidth, 0, width - forwardWidth, height);

        // Border
        ctx.strokeStyle = '#e5e7eb';
        ctx.strokeRect(0, 0, width, height);
    }

    private renderQualityChart(evidence: VariantEvidence): void {
        const ctx = this.qualityChart.getContext('2d')!;
        const width = this.qualityChart.width;
        const height = this.qualityChart.height;

        ctx.clearRect(0, 0, width, height);

        const histogram = evidence.quality_histogram;
        const bins = Object.keys(histogram);
        const values = Object.values(histogram);
        const maxVal = Math.max(...values, 1);

        const barWidth = width / bins.length - 4;
        const barGap = 4;

        ctx.font = '9px sans-serif';
        ctx.textAlign = 'center';

        bins.forEach((bin, i) => {
            const x = i * (barWidth + barGap) + 2;
            const barHeight = (histogram[bin] / maxVal) * (height - 16);

            // Color based on quality range
            if (bin.startsWith('30') || bin.startsWith('40')) {
                ctx.fillStyle = '#22c55e';  // Green for high quality
            } else if (bin.startsWith('20')) {
                ctx.fillStyle = '#f97316';  // Orange for medium
            } else {
                ctx.fillStyle = '#ef4444';  // Red for low
            }

            ctx.fillRect(x, height - 12 - barHeight, barWidth, barHeight);

            // Label
            ctx.fillStyle = '#6b7280';
            ctx.fillText(bin.split('-')[0], x + barWidth / 2, height - 2);
        });
    }

    private renderPositionChart(evidence: VariantEvidence): void {
        const ctx = this.positionChart.getContext('2d')!;
        const width = this.positionChart.width;
        const height = this.positionChart.height;

        ctx.clearRect(0, 0, width, height);

        const histogram = evidence.position_histogram;
        const bins = Object.keys(histogram);
        const values = Object.values(histogram);
        const maxVal = Math.max(...values, 1);

        const barWidth = width / bins.length - 4;
        const barGap = 4;

        ctx.font = '9px sans-serif';
        ctx.textAlign = 'center';

        bins.forEach((bin, i) => {
            const x = i * (barWidth + barGap) + 2;
            const barHeight = (histogram[bin] / maxVal) * (height - 16);

            // Color: edges are more suspicious (potential sequencing artifacts)
            const binStart = parseInt(bin.split('-')[0]);
            if (binStart < 25 || binStart >= 100) {
                ctx.fillStyle = '#f97316';  // Orange for edge positions
            } else {
                ctx.fillStyle = '#3b82f6';  // Blue for middle
            }

            ctx.fillRect(x, height - 12 - barHeight, barWidth, barHeight);

            // Label
            ctx.fillStyle = '#6b7280';
            ctx.fillText(bin.split('-')[0], x + barWidth / 2, height - 2);
        });
    }

    // ==================== PAIRED-END VISUALIZATION ====================

    private getInsertSizeStatus(read: Read): 'normal' | 'short' | 'long' {
        if (!read.insert_size) return 'normal';
        const absSize = Math.abs(read.insert_size);
        // Normal range: 150-500bp for typical WGS libraries
        if (absSize < 100) return 'short';
        if (absSize > 1000) return 'long';
        return 'normal';
    }

    private isSVCandidate(read: Read): boolean {
        // Chimeric (mate on different chromosome)
        if (read.mate_contig && read.mate_contig !== this.data?.contig) return true;

        // Short insert (<100bp) or long insert (>1kb)
        if (read.insert_size) {
            const absSize = Math.abs(read.insert_size);
            if (absSize < 100 || absSize > 1000) return true;
        }

        // Improper pair flag
        if (read.is_proper_pair === false) return true;

        return false;
    }

    private renderMateConnection(
        ctx: CanvasRenderingContext2D,
        read: Read,
        y: number,
        scale: number,
        isHovered: boolean = false
    ): void {
        // Skip if no mate info
        if (read.mate_position == null) return;

        const mateContig = read.mate_contig;
        const midY = y + READ_HEIGHT / 2;
        const width = this.readsCanvas.width;

        // Calculate read's connection point (3' end)
        const readConnectPos = read.is_read1 ? read.end_position : read.position;
        const x1 = (readConnectPos - this.viewport.start) * scale;

        // Skip if this read's connection point is completely off-screen
        if (x1 < -50 || x1 > width + 50) return;

        // Check if mate is actually loaded (not just that position is in viewport)
        const mateKey = `${read.name}:${read.position}`;
        const mate = this.mateIndex.get(mateKey);
        const mateLoaded = mate !== undefined;

        // Determine mate status
        const isChimeric = mateContig !== this.data?.contig;
        const insertStatus = this.getInsertSizeStatus(read);
        const isImproperPair = read.is_proper_pair === false;

        // Color scheme for SV candidates:
        // - Red: insert >1kb (potential deletion)
        // - Orange: insert <100bp or inter-chromosomal
        // - Gray dashed: improper pair flag
        // - Green: normal pair (only shown on hover)
        const getConnectorStyle = (): { color: string; dash: number[] } => {
            if (insertStatus === 'long') return { color: '#ef4444', dash: [] };           // Red
            if (insertStatus === 'short' || isChimeric) return { color: '#f97316', dash: [] };  // Orange
            if (isImproperPair) return { color: '#9ca3af', dash: [3, 3] };                // Gray dashed
            return { color: '#22c55e', dash: [] };                                         // Green (normal)
        };

        const style = getConnectorStyle();
        const color = style.color;
        const isAbnormal = isChimeric || insertStatus !== 'normal' || isImproperPair;

        if (isChimeric) {
            // CHIMERIC: Draw trailing line with chromosome warning
            ctx.strokeStyle = color;
            ctx.setLineDash([2, 2]);
            ctx.lineWidth = 1.5;

            const mateDirection = read.is_read1 ? 1 : -1;
            const trailLength = 25;

            // Trailing line
            ctx.beginPath();
            ctx.moveTo(x1, midY);
            ctx.lineTo(x1 + trailLength * mateDirection, midY);
            ctx.stroke();

            // Warning triangle with chr label
            const triX = x1 + (trailLength + 4) * mateDirection;
            ctx.fillStyle = color;
            ctx.font = 'bold 8px monospace';
            const label = `⚠ ${mateContig || '?'}`;
            ctx.fillText(label, mateDirection > 0 ? triX : triX - ctx.measureText(label).width, midY + 3);

            ctx.setLineDash([]);

        } else if (mateLoaded) {
            // BOTH READS LOADED: Draw curved arc to actual mate position
            // Get mate's connection point (5' end for read2, 3' end for read1)
            const mateConnectPos = read.is_read1 ? mate.position : mate.end_position;
            const x2 = (mateConnectPos - this.viewport.start) * scale;

            // Get mate's actual Y position from row index
            const mateRowIdx = this.readRowIndex.get(mate);
            if (mateRowIdx === undefined) {
                // Mate not in index, show direction indicator
                this.renderOffscreenMateIndicator(ctx, x1, midY, mateConnectPos, color, isAbnormal || isHovered, read);
                return;
            }
            const mateY = mateRowIdx * (READ_HEIGHT + READ_GAP) + READ_HEIGHT / 2;

            // Only draw if mate connection point is in viewport (X axis)
            if (x2 < -50 || x2 > width + 50) {
                // Mate loaded but X position off-screen - show direction indicator
                this.renderOffscreenMateIndicator(ctx, x1, midY, mateConnectPos, color, isAbnormal || isHovered, read);
                return;
            }

            ctx.strokeStyle = color;
            ctx.setLineDash(style.dash);
            ctx.lineWidth = isHovered ? 2 : (isAbnormal ? 1.5 : 1);

            // Arc connects from this read's Y to mate's Y
            const arcMidX = (x1 + x2) / 2;
            const arcMidY = Math.min(midY, mateY);
            const arcHeight = Math.min(18, Math.abs(x2 - x1) / 5 + Math.abs(mateY - midY) / 4);

            ctx.beginPath();
            ctx.moveTo(x1, midY);
            ctx.quadraticCurveTo(arcMidX, arcMidY - arcHeight, x2, mateY);
            ctx.stroke();

            // Insert size label (only if arc is wide enough and reads on same or adjacent rows)
            if (Math.abs(x2 - x1) > 60 && read.insert_size && Math.abs(mateY - midY) < 40) {
                const labelX = arcMidX;
                const labelY = arcMidY - arcHeight - 2;
                const absSize = Math.abs(read.insert_size);
                const sizeLabel = insertStatus === 'normal'
                    ? `${absSize}bp ✓`
                    : `⚠ ${absSize}bp`;

                ctx.font = '8px monospace';
                ctx.fillStyle = color;
                ctx.textAlign = 'center';
                ctx.fillText(sizeLabel, labelX, labelY);
                ctx.textAlign = 'start';
            }

            ctx.setLineDash([]);

        } else {
            // MATE NOT LOADED: Show direction indicator toward mate position
            const matePos = read.mate_position;
            // Use the style color for SV candidates, gray for normal (shouldn't show for normal pairs in default view)
            this.renderOffscreenMateIndicator(ctx, x1, midY, matePos, color, isAbnormal || isHovered, read);
        }
    }

    private renderOffscreenMateIndicator(
        ctx: CanvasRenderingContext2D,
        x1: number,
        midY: number,
        matePos: number,
        color: string,
        showLabel: boolean,
        read: Read
    ): void {
        ctx.strokeStyle = color;
        ctx.setLineDash([4, 3]);
        ctx.lineWidth = 1.5;

        const mateDirection = matePos > this.viewport.end ? 1 : -1;
        const trailLength = 20;

        // Trailing dashed line
        ctx.beginPath();
        ctx.moveTo(x1, midY);
        ctx.lineTo(x1 + trailLength * mateDirection, midY);
        ctx.stroke();

        // Arrow head
        const arrowX = x1 + trailLength * mateDirection;
        ctx.beginPath();
        ctx.moveTo(arrowX, midY);
        ctx.lineTo(arrowX - 4 * mateDirection, midY - 3);
        ctx.moveTo(arrowX, midY);
        ctx.lineTo(arrowX - 4 * mateDirection, midY + 3);
        ctx.stroke();

        // Insert size label (show for abnormal or when hovered)
        if (showLabel && read.insert_size) {
            ctx.font = '9px monospace';
            ctx.fillStyle = color;
            const absSize = Math.abs(read.insert_size);
            const arrow = mateDirection > 0 ? '→' : '←';
            const label = `${absSize}bp ${arrow}`;
            const labelX = mateDirection > 0 ? arrowX + 4 : arrowX - ctx.measureText(label).width - 4;
            ctx.fillText(label, labelX, midY - 4);
        }

        ctx.setLineDash([]);
    }
}


// Initialize the viewer
new BAMCPViewer();
