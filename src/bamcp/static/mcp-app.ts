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
}

interface Variant {
    contig: string;
    position: number;
    ref: string;
    alt: string;
    vaf: number;
    depth: number;
}

interface RegionData {
    contig: string;
    start: number;
    end: number;
    reads: Read[];
    coverage: number[];
    variants: Variant[];
    reference_sequence?: string;
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

    // State
    private data: RegionData | null = null;
    private viewport = { start: 0, end: 1000 };
    private packedRows: Read[][] = [];
    private app: App | null = null;
    private isFullscreen = false;
    private _contextUpdateTimer: ReturnType<typeof setTimeout> | null = null;
    private lockedTooltip: { read: Read; x: number; y: number } | null = null;

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

        // Tooltips
        this.readsCanvas.addEventListener('mousemove', (e) => {
            if (isDragging) return;
            this.showTooltip(e);
        });

        this.readsCanvas.addEventListener('mouseout', () => {
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
            }
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

    loadData(data: RegionData): void {
        this.data = data;
        this.viewport = { start: data.start, end: data.end };

        (document.getElementById('region-input') as HTMLInputElement).value =
            data.contig + ':' + data.start + '-' + data.end;

        this.packReads();
        this.renderVariantTable();
        this.resize();
        this.scheduleContextUpdate();
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
        const rowEnds: number[] = [];

        const sortedReads = [...this.data.reads].sort((a, b) => a.position - b.position);

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
        ctx.fillStyle = 'rgba(239, 68, 68, 0.15)';
        for (let i = 0; i < coverage.length; i++) {
            if (coverage[i] < lowCoverageThreshold) {
                const x = (this.data.start + i - this.viewport.start) * scale;
                ctx.fillRect(x, 0, Math.max(scale, 1), height);
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

        for (let rowIdx = 0; rowIdx < this.packedRows.length; rowIdx++) {
            const y = rowIdx * (READ_HEIGHT + READ_GAP);

            if (y > height) break;

            for (const read of this.packedRows[rowIdx]) {
                this.renderRead(ctx, read, y, scale, zoomLevel);
            }
        }
    }

    private renderRead(
        ctx: CanvasRenderingContext2D,
        read: Read,
        y: number,
        scale: number,
        zoomLevel: 'nucleotide' | 'read'
    ): void {
        const x = (read.position - this.viewport.start) * scale;
        const w = (read.end_position - read.position) * scale;

        if (x + w < 0 || x > this.readsCanvas.width) return;

        // MAPQ-based opacity (MAPQ 60 = full, 0 = 30%)
        const opacity = 0.3 + Math.min(read.mapping_quality, 60) / 60 * 0.7;

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
                    // Match/mismatch segment
                    ctx.fillStyle = `rgba(${baseColor[0]}, ${baseColor[1]}, ${baseColor[2]}, ${opacity})`;
                    ctx.fillRect(x, y, len * scale, READ_HEIGHT);

                    // Draw mismatches - always visible at any zoom level
                    for (const mm of read.mismatches) {
                        if (mm.pos >= refPos && mm.pos < refPos + len) {
                            const mx = (mm.pos - this.viewport.start) * scale;

                            if (zoomLevel === 'nucleotide' || scale > 2) {
                                // Full mismatch rectangle at higher zoom
                                ctx.fillStyle = BASE_COLORS[mm.alt] || '#9ca3af';
                                ctx.fillRect(mx, y, Math.max(scale, 2), READ_HEIGHT);

                                // Draw letter at nucleotide zoom
                                if (zoomLevel === 'nucleotide' && scale >= 10) {
                                    ctx.fillStyle = '#fff';
                                    ctx.font = 'bold 9px monospace';
                                    ctx.textAlign = 'center';
                                    ctx.fillText(mm.alt, mx + scale / 2, y + READ_HEIGHT - 2);
                                    ctx.textAlign = 'start';
                                }
                            } else {
                                // Low zoom: bright tick marker for visibility
                                ctx.fillStyle = BASE_COLORS[mm.alt] || '#ef4444';
                                ctx.fillRect(mx, y, Math.max(2, scale), READ_HEIGHT);

                                // Triangle indicator above read
                                ctx.beginPath();
                                ctx.moveTo(mx, y - 2);
                                ctx.lineTo(mx + 3, y);
                                ctx.lineTo(mx - 1, y);
                                ctx.closePath();
                                ctx.fill();
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

        this.tooltip.innerHTML =
            '<strong>' + read.name + '</strong><br>' +
            'Position: ' + read.position.toLocaleString() + '-' + read.end_position.toLocaleString() + '<br>' +
            'MAPQ: <span style="color:' + mapqColor + '">' + read.mapping_quality + '</span><br>' +
            'CIGAR: <code>' + read.cigar + '</code><br>' +
            'Strand: ' + (read.is_reverse ? '← Reverse' : '→ Forward') + '<br>' +
            'Length: ' + read.sequence.length + ' bp';

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

    private renderVariantTable(): void {
        if (!this.data) return;

        this.variantTable.innerHTML = this.data.variants.map((v) => {
            const vafColor = v.vaf >= 0.4 ? '#22c55e' : v.vaf >= 0.2 ? '#f97316' : '#ef4444';
            return '<tr data-pos="' + v.position + '" ' +
                'data-ref="' + v.ref + '" data-alt="' + v.alt + '" ' +
                'data-vaf="' + v.vaf + '" data-depth="' + v.depth + '" ' +
                'data-contig="' + v.contig + '">' +
                '<td>' + v.contig + ':' + v.position.toLocaleString() + '</td>' +
                '<td style="font-family:monospace">' + v.ref + '</td>' +
                '<td style="color:' + (BASE_COLORS[v.alt] || '#333') + ';font-family:monospace;font-weight:bold">' + v.alt + '</td>' +
                '<td style="color:' + vafColor + '">' + (v.vaf * 100).toFixed(1) + '%</td>' +
                '<td>' + v.depth + '</td>' +
                '</tr>';
        }).join('');

        this.variantTable.querySelectorAll('tr').forEach(row => {
            row.addEventListener('click', () => {
                const el = row as HTMLElement;
                const pos = parseInt(el.dataset.pos!);
                this.jumpTo(pos);

                const variant: Variant = {
                    contig: el.dataset.contig!,
                    position: pos,
                    ref: el.dataset.ref!,
                    alt: el.dataset.alt!,
                    vaf: parseFloat(el.dataset.vaf!),
                    depth: parseInt(el.dataset.depth!)
                };
                this.sendVariantMessage(variant);
            });
        });
    }
}

// Initialize the viewer
new BAMCPViewer();
