/**
 * BAMCP Alignment Viewer - MCP Apps client
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

// Color scheme for nucleotides
const BASE_COLORS: Record<string, string> = {
    'A': '#22c55e',
    'T': '#ef4444',
    'G': '#f97316',
    'C': '#3b82f6',
    'N': '#9ca3af'
};

const READ_HEIGHT = 12;
const READ_GAP = 2;
const COVERAGE_HEIGHT = 60;
const CONTEXT_UPDATE_DEBOUNCE_MS = 300;

class BAMCPViewer {
    private coverageCanvas: HTMLCanvasElement;
    private readsCanvas: HTMLCanvasElement;
    private coverageCtx: CanvasRenderingContext2D;
    private readsCtx: CanvasRenderingContext2D;
    private tooltip: HTMLElement;
    private variantTable: HTMLElement;

    private data: RegionData | null = null;
    private viewport = { start: 0, end: 1000 };
    private packedRows: Read[][] = [];
    private app: App | null = null;
    private isFullscreen = false;
    private _contextUpdateTimer: ReturnType<typeof setTimeout> | null = null;

    constructor() {
        this.coverageCanvas = document.getElementById('coverage-canvas') as HTMLCanvasElement;
        this.readsCanvas = document.getElementById('reads-canvas') as HTMLCanvasElement;
        this.coverageCtx = this.coverageCanvas.getContext('2d')!;
        this.readsCtx = this.readsCanvas.getContext('2d')!;
        this.tooltip = document.getElementById('tooltip')!;
        this.variantTable = document.getElementById('variant-table')!;

        this.setupEventListeners();
        this.initApp();
        this.resize();

        window.addEventListener('resize', () => this.resize());
    }

    private async initApp(): Promise<void> {
        try {
            this.app = new App({ name: 'BAMCP Viewer', version: '1.0.0' });

            // Set ontoolresult BEFORE connect() to receive the initial tool result
            this.app.ontoolresult = (params: ToolResultParams) => {
                // Prefer structuredContent (hidden from LLM, full data for UI)
                if (params.structuredContent) {
                    this.loadData(params.structuredContent);
                } else if (params.content) {
                    // Fall back to parsing from text content
                    const text = params.content.find(c => c.type === 'text');
                    if (text?.text) {
                        this.loadData(JSON.parse(text.text));
                    }
                }
            };

            await this.app.connect();
        } catch (e) {
            // SDK not available or connect failed — viewer still works
            // when data is loaded directly via loadData()
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

        this.readsCanvas.addEventListener('mousemove', (e) => {
            if (isDragging) return;
            this.showTooltip(e);
        });

        this.readsCanvas.addEventListener('mouseout', () => {
            this.tooltip.style.display = 'none';
        });

        this.readsCanvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = e.deltaY > 0 ? 1.2 : 0.8;
            this.zoom(factor);
        });
    }

    private async requestRegion(region: string): Promise<void> {
        if (this.app) {
            try {
                // MCP Apps UI cannot call tools directly - ask the LLM to do it
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
        // Resize adjusts canvas height based on packed rows
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

        // Use browser's native fullscreen API (works in iframe with allow="fullscreen")
        const container = document.getElementById('container');
        if (container) {
            try {
                if (this.isFullscreen && !document.fullscreenElement) {
                    await container.requestFullscreen();
                } else if (!this.isFullscreen && document.fullscreenElement) {
                    await document.exitFullscreen();
                }
            } catch {
                // Fullscreen not supported or denied - just toggle the button state
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
        // Ensure minimum height for reads display (at least 200px for reads)
        const containerHeight = Math.max(container.clientHeight, COVERAGE_HEIGHT + 200);
        const height = containerHeight - COVERAGE_HEIGHT;

        this.coverageCanvas.width = width;
        this.coverageCanvas.height = COVERAGE_HEIGHT;

        this.readsCanvas.width = width;
        this.readsCanvas.height = height;

        // If we have packed reads, ensure canvas is tall enough to show them
        if (this.packedRows.length > 0) {
            const neededHeight = this.packedRows.length * (READ_HEIGHT + READ_GAP);
            if (neededHeight > this.readsCanvas.height) {
                this.readsCanvas.height = Math.min(neededHeight, 800);
            }
        }

        this.render();
    }

    private render(): void {
        this.renderCoverage();
        this.renderReads();
    }

    private renderCoverage(): void {
        const ctx = this.coverageCtx;
        const width = this.coverageCanvas.width;
        const height = this.coverageCanvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!this.data) return;

        const coverage = this.data.coverage;
        const maxCov = Math.max(...coverage, 1);
        const scale = this.getScale();

        ctx.fillStyle = '#93c5fd';
        ctx.strokeStyle = '#3b82f6';
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

        ctx.fillStyle = '#333';
        ctx.font = '10px sans-serif';
        ctx.fillText('Max: ' + maxCov + 'x', 5, 12);
    }

    private renderReads(): void {
        const ctx = this.readsCtx;
        const width = this.readsCanvas.width;
        const height = this.readsCanvas.height;

        ctx.clearRect(0, 0, width, height);

        if (!this.data) return;

        const scale = this.getScale();

        for (let rowIdx = 0; rowIdx < this.packedRows.length; rowIdx++) {
            const y = rowIdx * (READ_HEIGHT + READ_GAP);

            if (y > height) break;

            for (const read of this.packedRows[rowIdx]) {
                this.renderRead(ctx, read, y, scale);
            }
        }
    }

    private renderRead(ctx: CanvasRenderingContext2D, read: Read, y: number, scale: number): void {
        const x = (read.position - this.viewport.start) * scale;
        const w = (read.end_position - read.position) * scale;

        if (x + w < 0 || x > this.readsCanvas.width) return;

        // Parse soft-clips from CIGAR
        const softClips = this.parseSoftClips(read.cigar);

        // Draw soft-clipped regions (semi-transparent with amber border)
        if (softClips.left > 0) {
            const clipW = softClips.left * scale;
            ctx.fillStyle = read.is_reverse ? 'rgba(167,139,250,0.3)' : 'rgba(96,165,250,0.3)';
            ctx.fillRect(x - clipW, y, clipW, READ_HEIGHT);
            ctx.strokeStyle = '#f59e0b';
            ctx.lineWidth = 1;
            ctx.strokeRect(x - clipW, y, clipW, READ_HEIGHT);
        }
        if (softClips.right > 0) {
            const clipW = softClips.right * scale;
            ctx.fillStyle = read.is_reverse ? 'rgba(167,139,250,0.3)' : 'rgba(96,165,250,0.3)';
            ctx.fillRect(x + w, y, clipW, READ_HEIGHT);
            ctx.strokeStyle = '#f59e0b';
            ctx.lineWidth = 1;
            ctx.strokeRect(x + w, y, clipW, READ_HEIGHT);
        }

        // Draw aligned read body
        ctx.fillStyle = read.is_reverse ? '#a78bfa' : '#60a5fa';
        ctx.fillRect(x, y, w, READ_HEIGHT);

        if (scale > 2) {
            for (const mm of read.mismatches) {
                const mx = (mm.pos - this.viewport.start) * scale;
                ctx.fillStyle = BASE_COLORS[mm.alt] || '#9ca3af';
                ctx.fillRect(mx, y, Math.max(scale, 2), READ_HEIGHT);
            }
        }

        ctx.fillStyle = 'white';
        ctx.font = '8px sans-serif';
        const arrow = read.is_reverse ? '\u25C0' : '\u25B6';
        if (w > 15) {
            ctx.fillText(arrow, read.is_reverse ? x + 2 : x + w - 10, y + 9);
        }
    }

    private parseSoftClips(cigar: string): { left: number; right: number } {
        let left = 0, right = 0;
        if (!cigar) return { left, right };
        const ops = cigar.match(/\d+[MIDNSHP=X]/g);
        if (!ops || ops.length === 0) return { left, right };
        if (ops[0].endsWith('S')) left = parseInt(ops[0]);
        if (ops.length > 1 && ops[ops.length - 1].endsWith('S')) right = parseInt(ops[ops.length - 1]);
        return { left, right };
    }

    private getScale(): number {
        const width = this.readsCanvas.width;
        return width / (this.viewport.end - this.viewport.start);
    }

    private zoom(factor: number): void {
        const center = (this.viewport.start + this.viewport.end) / 2;
        const span = (this.viewport.end - this.viewport.start) * factor;

        this.viewport.start = Math.floor(center - span / 2);
        this.viewport.end = Math.ceil(center + span / 2);

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

    private showTooltip(e: MouseEvent): void {
        if (!this.data) return;

        const rect = this.readsCanvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;

        const scale = this.getScale();
        const genomePos = Math.floor(this.viewport.start + x / scale);
        const rowIdx = Math.floor(y / (READ_HEIGHT + READ_GAP));

        if (rowIdx >= this.packedRows.length) {
            this.tooltip.style.display = 'none';
            return;
        }

        const read = this.packedRows[rowIdx].find(r =>
            r.position <= genomePos && r.end_position >= genomePos
        );

        if (read) {
            this.tooltip.innerHTML =
                '<strong>' + read.name + '</strong><br>' +
                'Position: ' + read.position + '-' + read.end_position + '<br>' +
                'MAPQ: ' + read.mapping_quality + '<br>' +
                'CIGAR: ' + read.cigar + '<br>' +
                'Strand: ' + (read.is_reverse ? 'Reverse' : 'Forward');
            this.tooltip.style.display = 'block';
            this.tooltip.style.left = (e.clientX + 10) + 'px';
            this.tooltip.style.top = (e.clientY + 10) + 'px';
        } else {
            this.tooltip.style.display = 'none';
        }
    }

    private renderVariantTable(): void {
        if (!this.data) return;

        this.variantTable.innerHTML = this.data.variants.map((v) => {
            return '<tr data-pos="' + v.position + '" ' +
                'data-ref="' + v.ref + '" data-alt="' + v.alt + '" ' +
                'data-vaf="' + v.vaf + '" data-depth="' + v.depth + '" ' +
                'data-contig="' + v.contig + '">' +
                '<td>' + v.contig + ':' + v.position + '</td>' +
                '<td>' + v.ref + '</td>' +
                '<td style="color: ' + (BASE_COLORS[v.alt] || '#333') + '">' + v.alt + '</td>' +
                '<td>' + (v.vaf * 100).toFixed(1) + '%</td>' +
                '<td>' + v.depth + '</td>' +
                '</tr>';
        }).join('');

        this.variantTable.querySelectorAll('tr').forEach(row => {
            row.addEventListener('click', () => {
                const el = row as HTMLElement;
                const pos = parseInt(el.dataset.pos!);
                this.jumpTo(pos);

                // Send variant explanation request via App SDK
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

    private jumpTo(position: number): void {
        const span = this.viewport.end - this.viewport.start;
        this.viewport.start = position - span / 2;
        this.viewport.end = position + span / 2;
        this.render();
        this.scheduleContextUpdate();
    }
}

// Initialize the viewer
new BAMCPViewer();
