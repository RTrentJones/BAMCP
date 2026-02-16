import { StateManager } from "./state";
import { Read } from "./types";
import {
    BASE_COLORS,
    COLOR_PALETTES,
    DISPLAY_MODE_CONFIGS,
    INSERT_SIZE_THRESHOLDS,
    SOFT_CLIP_STYLE,
} from "./constants";

// Layout constants
const RULER_HEIGHT = 20;
const REFERENCE_HEIGHT = 24;
const COVERAGE_HEIGHT = 60;

export class Renderer {
    private rulerCanvas: HTMLCanvasElement;
    private referenceCanvas: HTMLCanvasElement;
    private coverageCanvas: HTMLCanvasElement;
    private readsCanvas: HTMLCanvasElement;

    private rulerCtx: CanvasRenderingContext2D;
    private referenceCtx: CanvasRenderingContext2D;
    private coverageCtx: CanvasRenderingContext2D;
    private readsCtx: CanvasRenderingContext2D;

    private state: StateManager;

    constructor(state: StateManager) {
        this.state = state;

        this.rulerCanvas = document.getElementById('ruler-canvas') as HTMLCanvasElement;
        this.referenceCanvas = document.getElementById('reference-canvas') as HTMLCanvasElement;
        this.coverageCanvas = document.getElementById('coverage-canvas') as HTMLCanvasElement;
        this.readsCanvas = document.getElementById('reads-canvas') as HTMLCanvasElement;

        this.rulerCtx = this.rulerCanvas.getContext('2d')!;
        this.referenceCtx = this.referenceCanvas.getContext('2d')!;
        this.coverageCtx = this.coverageCanvas.getContext('2d')!;
        this.readsCtx = this.readsCanvas.getContext('2d')!;
    }

    /**
     * Get read dimensions based on current display mode setting.
     */
    private getReadDimensions(): { height: number; gap: number; showLabels: boolean } {
        const mode = this.state.settings.displayMode;
        const config = DISPLAY_MODE_CONFIGS[mode] || DISPLAY_MODE_CONFIGS.compact;
        return {
            height: config.readHeight,
            gap: config.readGap,
            showLabels: config.showLabels,
        };
    }

    /**
     * Get color based on insert size for mate-pair visualization.
     */
    private getInsertSizeColor(read: Read): string {
        const insertSize = read.insert_size ? Math.abs(read.insert_size) : null;
        const isChimeric = read.mate_contig && read.mate_contig !== this.state.data?.contig;
        const palette = COLOR_PALETTES.insertSize;

        if (isChimeric) return palette.chimeric;
        if (!insertSize) return palette.normal;
        if (insertSize < INSERT_SIZE_THRESHOLDS.tooShort) return palette.small;
        if (insertSize > INSERT_SIZE_THRESHOLDS.veryLong) return palette.veryLarge;
        if (insertSize > INSERT_SIZE_THRESHOLDS.tooLong) return palette.large;
        if (insertSize >= INSERT_SIZE_THRESHOLDS.normalMin &&
            insertSize <= INSERT_SIZE_THRESHOLDS.normalMax) return '#22c55e';  // Green - normal
        return palette.normal;
    }

    /**
     * Get read color based on current colorBy setting.
     */
    private getReadColor(read: Read): string {
        const colorBy = this.state.settings.colorBy;

        switch (colorBy) {
            case 'strand':
                return read.is_reverse
                    ? COLOR_PALETTES.strand.reverse
                    : COLOR_PALETTES.strand.forward;

            case 'mapq':
                return this.getColorFromScale(read.mapping_quality, COLOR_PALETTES.mapq);

            case 'insertSize':
                return this.getInsertSizeColor(read);

            case 'baseQuality':
                // Use average base quality if available
                // For now, fall back to MAPQ as a proxy
                return this.getColorFromScale(read.mapping_quality, COLOR_PALETTES.baseQuality);

            default:
                return '#94a3b8';  // Slate gray fallback
        }
    }

    /**
     * Get color from a threshold-based scale.
     */
    private getColorFromScale(
        value: number,
        scale: Array<{ threshold: number; color: string }>
    ): string {
        for (let i = scale.length - 1; i >= 0; i--) {
            if (value >= scale[i].threshold) {
                return scale[i].color;
            }
        }
        return scale[0].color;
    }

    /**
     * Check if a position is within the current viewport.
     */
    private isInViewport(position: number): boolean {
        return position >= this.state.viewport.start && position <= this.state.viewport.end;
    }

    public resize(): void {
        const container = document.getElementById('viewer')!;
        const width = Math.max(container.clientWidth, 300);

        // Ruler canvas
        this.rulerCanvas.width = width;
        this.rulerCanvas.height = RULER_HEIGHT;

        // Reference canvas
        this.referenceCanvas.width = width;
        this.referenceCanvas.height = REFERENCE_HEIGHT;

        // Coverage canvas
        this.coverageCanvas.width = width;
        this.coverageCanvas.height = COVERAGE_HEIGHT;

        // Calculate height needed for actual reads content
        let readsHeight = 200; // Minimum height
        if (this.state.packedRows.length > 0) {
            const { height, gap } = this.getReadDimensions();
            readsHeight = Math.max(readsHeight, this.state.packedRows.length * (height + gap) + 20);
        }

        // Reads canvas - size to content
        this.readsCanvas.width = width;
        this.readsCanvas.height = readsHeight;

        this.render();
    }

    public render(): void {
        this.renderRuler();
        this.renderReference();
        this.renderCoverage();
        this.renderReads();
    }

    private getScale(): number {
        return this.readsCanvas.width / (this.state.viewport.end - this.state.viewport.start);
    }

    // ==================== RULER TRACK ====================

    private renderRuler(): void {
        const ctx = this.rulerCtx;
        const width = this.rulerCanvas.width;
        const height = this.rulerCanvas.height;
        const data = this.state.data;

        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = '#fafafa';
        ctx.fillRect(0, 0, width, height);

        if (!data) return;

        const scale = this.getScale();
        const span = this.state.viewport.end - this.state.viewport.start;
        const tickInterval = this.calculateTickInterval(span);

        ctx.fillStyle = '#374151';
        ctx.font = '10px monospace';
        ctx.strokeStyle = '#9ca3af';
        ctx.lineWidth = 1;

        // Draw region label
        ctx.fillStyle = '#6b7280';
        ctx.font = 'bold 10px sans-serif';
        ctx.fillText(data.contig, 4, 12);

        ctx.font = '10px monospace';
        ctx.fillStyle = '#374151';

        const startTick = Math.ceil(this.state.viewport.start / tickInterval) * tickInterval;
        for (let pos = startTick; pos <= this.state.viewport.end; pos += tickInterval) {
            const x = (pos - this.state.viewport.start) * scale;

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
        const data = this.state.data;

        ctx.clearRect(0, 0, width, height);

        if (!data?.reference_sequence) {
            ctx.fillStyle = '#f9fafb';
            ctx.fillRect(0, 0, width, height);
            ctx.fillStyle = '#9ca3af';
            ctx.font = '11px sans-serif';
            ctx.fillText('Reference sequence not available', 8, height / 2 + 4);
            return;
        }

        const scale = this.getScale();
        const refSeq = data.reference_sequence;

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

        const data = this.state.data!;
        const startIdx = Math.max(0, Math.floor(this.state.viewport.start - data.start));
        const endIdx = Math.min(refSeq.length, Math.ceil(this.state.viewport.end - data.start));

        for (let i = startIdx; i < endIdx; i++) {
            const base = refSeq[i].toUpperCase();
            const x = (data.start + i - this.state.viewport.start) * scale;

            ctx.fillStyle = BASE_COLORS[base] || '#9ca3af';
            ctx.fillRect(x, 2, Math.max(scale - 1, 1), height - 4);

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
        const data = this.state.data!;
        const startIdx = Math.max(0, Math.floor(this.state.viewport.start - data.start));
        const endIdx = Math.min(refSeq.length, Math.ceil(this.state.viewport.end - data.start));

        for (let i = startIdx; i < endIdx; i++) {
            const base = refSeq[i].toUpperCase();
            const x = (data.start + i - this.state.viewport.start) * scale;

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
        const pixelWidth = 1 / scale;
        const data = this.state.data!;
        const viewStart = Math.max(0, this.state.viewport.start - data.start);

        for (let px = 0; px < width; px++) {
            const baseStart = Math.floor(viewStart + px * pixelWidth);
            const baseEnd = Math.min(refSeq.length, Math.ceil(viewStart + (px + 1) * pixelWidth));

            if (baseStart >= refSeq.length) break;

            const counts: Record<string, number> = { A: 0, T: 0, G: 0, C: 0 };
            for (let i = baseStart; i < baseEnd; i++) {
                const base = refSeq[i].toUpperCase();
                if (counts[base] !== undefined) counts[base]++;
            }

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
        const data = this.state.data;

        ctx.clearRect(0, 0, width, height);

        if (!data) return;

        const coverage = data.coverage;
        const maxCov = Math.max(...coverage, 1);
        const scale = this.getScale();
        const lowCoverageThreshold = 10;

        // Draw low coverage warning regions
        for (let i = 0; i < coverage.length; i++) {
            if (coverage[i] < lowCoverageThreshold) {
                const x = (data.start + i - this.state.viewport.start) * scale;
                if (coverage[i] === 0) {
                    ctx.fillStyle = 'rgba(239, 68, 68, 0.4)';
                    ctx.fillRect(x, 0, Math.max(scale, 1), height);
                    ctx.strokeStyle = 'rgba(239, 68, 68, 0.6)';
                    ctx.lineWidth = 1;
                    for (let ly = 0; ly < height; ly += 4) {
                        ctx.beginPath();
                        ctx.moveTo(x, ly);
                        ctx.lineTo(x + Math.max(scale, 2), ly + 4);
                        ctx.stroke();
                    }
                } else {
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
            const x = (data.start + i - this.state.viewport.start) * scale;
            const h = (coverage[i] / maxCov) * (height - 20);
            const y = height - h;

            if (scale >= 1) {
                // Step function: horizontal bar per base, vertical transitions
                if (i === 0) {
                    ctx.moveTo(x, y);
                } else {
                    ctx.lineTo(x, y);
                }
                ctx.lineTo(x + scale, y);
            } else {
                // Smooth polygon for sub-pixel zoom (many bases per pixel)
                if (i === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            }
        }

        const xEnd = (data.end - this.state.viewport.start) * scale;
        const xStart = (data.start - this.state.viewport.start) * scale;
        ctx.lineTo(xEnd, height);
        ctx.lineTo(xStart, height);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();

        // Warning line
        const thresholdY = height - (lowCoverageThreshold / maxCov) * (height - 20);
        if (thresholdY > 15 && thresholdY < height - 5) {
            ctx.strokeStyle = '#ef4444';
            ctx.setLineDash([4, 2]);
            ctx.beginPath();
            ctx.moveTo(0, thresholdY);
            ctx.lineTo(width, thresholdY);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.fillStyle = '#ef4444';
            ctx.font = '9px sans-serif';
            ctx.fillText(lowCoverageThreshold + 'x', width - 25, thresholdY - 2);
        }

        // Labels
        ctx.fillStyle = '#374151';
        ctx.font = 'bold 10px sans-serif';
        ctx.fillText('Max: ' + maxCov + 'x', 5, 12);

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
        const scale = this.getScale();
        const data = this.state.data;
        const { height: READ_HEIGHT, gap: READ_GAP, showLabels } = this.getReadDimensions();

        ctx.clearRect(0, 0, width, height);

        if (!data) return;

        // Clip reads to data range so they don't extend beyond coverage fill
        const clipX = (data.start - this.state.viewport.start) * scale;
        const clipW = (data.end - data.start) * scale;
        ctx.save();
        ctx.beginPath();
        ctx.rect(clipX, 0, clipW, height);
        ctx.clip();

        // Draw reads
        for (let rowIndex = 0; rowIndex < this.state.packedRows.length; rowIndex++) {
            const row = this.state.packedRows[rowIndex];
            const y = rowIndex * (READ_HEIGHT + READ_GAP);

            if (y > height) break;

            for (const read of row) {
                // Calculate display coordinates
                const x1 = (read.position - this.state.viewport.start) * scale;
                const x2 = (read.end_position - this.state.viewport.start) * scale;
                const w = Math.max(x2 - x1, 1);

                // Skip if off-screen
                if (x2 < 0 || x1 > width) continue;

                // Determine styling based on colorBy setting
                let color = this.getReadColor(read);
                let opacity = Math.max(read.mapping_quality / 60, 0.2);
                let stroke = false;

                // Low MAPQ logic - always show MAPQ=0 reads distinctly
                if (read.mapping_quality === 0) {
                    color = '#fff';
                    stroke = true;
                    opacity = 1;
                }

                // Selection / Hover logic
                if (this.state.hoveredRead) {
                    const hovered = this.state.hoveredRead;
                    // Check if this read is the mate of the hovered read
                    // We check both directions in case keying is weird, but state.mateIndex should handle it
                    const isMate = this.state.mateIndex.get(`${hovered.name}:${hovered.position}`) === read ||
                        this.state.mateIndex.get(`${read.name}:${read.position}`) === hovered;

                    if (read === hovered) {
                        ctx.strokeStyle = '#000';
                        ctx.lineWidth = 2;
                        stroke = true;
                        opacity = 1;

                        // Check if mate is off-screen and draw indicator
                        const mate = this.state.mateIndex.get(`${read.name}:${read.position}`);
                        const matePos = read.mate_position;
                        if (matePos && !mate) {
                            // Mate exists but is not in viewport - draw off-screen indicator
                            const arcColor = this.getInsertSizeColor(read);
                            const isChimeric = read.mate_contig && read.mate_contig !== this.state.data?.contig;
                            const direction = matePos > this.state.viewport.end ? 'right' : 'left';
                            const edgeX = direction === 'right' ? width : 0;
                            const myY = y + READ_HEIGHT / 2;
                            const readEndX = direction === 'right' ? x2 : x1;

                            ctx.save();
                            // Dashed line to edge
                            ctx.strokeStyle = arcColor;
                            ctx.setLineDash([5, 3]);
                            ctx.lineWidth = 1.5;
                            ctx.beginPath();
                            ctx.moveTo(readEndX, myY);
                            ctx.lineTo(edgeX, myY);
                            ctx.stroke();

                            // Arrow at edge
                            const arrowDir = direction === 'right' ? -1 : 1;
                            ctx.beginPath();
                            ctx.moveTo(edgeX, myY);
                            ctx.lineTo(edgeX + 6 * arrowDir, myY - 4);
                            ctx.moveTo(edgeX, myY);
                            ctx.lineTo(edgeX + 6 * arrowDir, myY + 4);
                            ctx.stroke();

                            // Position label
                            ctx.fillStyle = arcColor;
                            ctx.font = '9px monospace';
                            const insertLabel = read.insert_size ? `${Math.abs(read.insert_size)}bp` : '';
                            const label = isChimeric
                                ? `⚠ ${read.mate_contig}`
                                : `${insertLabel} → ${matePos.toLocaleString()}`;
                            const labelWidth = ctx.measureText(label).width;
                            const labelX = direction === 'right' ? edgeX - labelWidth - 8 : 8;
                            ctx.fillText(label, labelX, myY - 8);

                            ctx.setLineDash([]);
                            ctx.restore();
                        }
                    } else if (isMate) {
                        // Get insert-size-based color for mate pair
                        const arcColor = this.getInsertSizeColor(hovered);
                        ctx.strokeStyle = arcColor;
                        ctx.lineWidth = 1.5;
                        stroke = true;
                        opacity = 1;

                        // Draw connector line
                        const midX = x1 + w / 2; // Current read (mate) center

                        // Calculate hovered read's center
                        const hX1 = (hovered.position - this.state.viewport.start) * scale;
                        const hX2 = (hovered.end_position - this.state.viewport.start) * scale;
                        const hMidX = hX1 + Math.max(hX2 - hX1, 1) / 2;

                        const hoveredY = (this.state.readRowIndex.get(hovered) || 0) * (READ_HEIGHT + READ_GAP) + READ_HEIGHT / 2;
                        const myY = y + READ_HEIGHT / 2;

                        ctx.save();
                        ctx.beginPath();
                        ctx.moveTo(midX, myY);

                        // Quadratic curve for smooth connection
                        const cpX = (midX + hMidX) / 2;
                        // Offset control point upward for a nicer arc
                        const arcHeight = Math.min(25, Math.abs(hMidX - midX) / 4 + Math.abs(hoveredY - myY) / 3);
                        const cpY = Math.min(myY, hoveredY) - arcHeight;

                        ctx.quadraticCurveTo(cpX, cpY, hMidX, hoveredY);

                        ctx.strokeStyle = arcColor;
                        ctx.lineWidth = 2;
                        ctx.setLineDash([]); // Solid line
                        ctx.stroke();

                        // Add insert size label at arc midpoint
                        if (hovered.insert_size && Math.abs(hMidX - midX) > 50) {
                            const insertSize = Math.abs(hovered.insert_size);
                            const isChimeric = hovered.mate_contig && hovered.mate_contig !== this.state.data?.contig;
                            const label = isChimeric
                                ? `⚠ ${hovered.mate_contig}`
                                : `${insertSize}bp`;

                            ctx.fillStyle = arcColor;
                            ctx.font = '9px monospace';
                            ctx.textAlign = 'center';
                            ctx.fillText(label, cpX, cpY - 2);
                            ctx.textAlign = 'start';
                        }

                        ctx.restore();
                    } else {
                        // Dim unrelated reads
                        opacity = 0.1;
                    }
                }

                // Draw Read Body
                ctx.fillStyle = read.mapping_quality === 0 ? '#fff' : color;
                ctx.globalAlpha = opacity;

                // Arrow shape
                ctx.beginPath();
                if (read.is_reverse) {
                    ctx.moveTo(x1 + 5, y);
                    ctx.lineTo(x2, y);
                    ctx.lineTo(x2, y + READ_HEIGHT);
                    ctx.lineTo(x1 + 5, y + READ_HEIGHT);
                    ctx.lineTo(x1, y + READ_HEIGHT / 2);
                } else {
                    ctx.moveTo(x1, y);
                    ctx.lineTo(x2 - 5, y);
                    ctx.lineTo(x2, y + READ_HEIGHT / 2);
                    ctx.lineTo(x2 - 5, y + READ_HEIGHT);
                    ctx.lineTo(x1, y + READ_HEIGHT);
                }
                ctx.closePath();
                ctx.fill();

                if (stroke) {
                    ctx.globalAlpha = 1;
                    ctx.stroke();
                }

                ctx.globalAlpha = 1; // Reset

                // Render soft clips if enabled
                this.renderSoftClips(ctx, read, x1, x2, y, READ_HEIGHT, scale);

                // High Zoom: Render Bases
                if (scale >= 10 && read.sequence) {
                    ctx.font = '10px monospace';
                    ctx.textAlign = 'center';
                    ctx.textBaseline = 'middle';

                    const seq = read.sequence;
                    // Draw each base with colored background and contrasting text
                    for (let i = 0; i < seq.length; i++) {
                        const baseX = x1 + (i * scale);
                        if (baseX + scale < 0 || baseX > width) continue;

                        const base = seq[i].toUpperCase();
                        const bgColor = BASE_COLORS[base] || '#374151';

                        // Draw colored background for the base
                        ctx.fillStyle = bgColor;
                        ctx.fillRect(baseX, y, Math.max(scale - 0.5, 1), READ_HEIGHT);

                        // Draw white text on colored background for contrast
                        ctx.fillStyle = '#fff';
                        ctx.fillText(base, baseX + scale / 2, y + READ_HEIGHT / 2);
                    }
                }
                // Medium Zoom: Draw Mismatches
                else if (scale > 0.5) {
                    for (const m of read.mismatches) {
                        const mx = (m.pos - this.state.viewport.start) * scale;
                        if (mx >= x1 && mx <= x2) {
                            ctx.fillStyle = BASE_COLORS[m.alt] || '#000';
                            ctx.fillRect(mx, y, Math.max(scale, 1), READ_HEIGHT);
                        }
                    }
                }
            }
        }

        // Draw variant circles on top of reads
        this.renderVariantMarkers(ctx, scale, width);

        ctx.restore(); // Remove data-range clip
    }

    /**
     * Render variant circles directly on reads that carry the variant.
     */
    private renderVariantMarkers(ctx: CanvasRenderingContext2D, scale: number, width: number): void {
        const data = this.state.data;
        if (!data || !data.variants.length) return;

        // Only show markers if zoomed in enough
        if (scale < 0.5) return;

        const { height: READ_HEIGHT, gap: READ_GAP } = this.getReadDimensions();

        // Build a map of variant positions for quick lookup
        const variantMap = new Map<number, { alt: string; ref: string }[]>();
        for (const variant of data.variants) {
            const existing = variantMap.get(variant.position);
            if (existing) {
                existing.push({ alt: variant.alt, ref: variant.ref });
            } else {
                variantMap.set(variant.position, [{ alt: variant.alt, ref: variant.ref }]);
            }
        }

        // Circle size based on zoom level
        const circleRadius = Math.max(3, Math.min(6, scale / 3));

        ctx.save();

        // Draw circles on each read that has a variant mismatch
        for (let rowIdx = 0; rowIdx < this.state.packedRows.length; rowIdx++) {
            const y = rowIdx * (READ_HEIGHT + READ_GAP);
            if (y > ctx.canvas.height) break;

            for (const read of this.state.packedRows[rowIdx]) {
                for (const mm of read.mismatches) {
                    // Check if this mismatch matches a called variant
                    const variantsAtPos = variantMap.get(mm.pos);
                    if (!variantsAtPos) continue;

                    const isVariant = variantsAtPos.some(v => v.alt === mm.alt);
                    if (!isVariant) continue;

                    const mx = (mm.pos - this.state.viewport.start) * scale;

                    // Skip if off-screen
                    if (mx < -circleRadius || mx > width + circleRadius) continue;

                    const centerX = mx + scale / 2;
                    const centerY = y + READ_HEIGHT / 2;

                    // Outer glow/highlight
                    ctx.beginPath();
                    ctx.arc(centerX, centerY, circleRadius + 2, 0, Math.PI * 2);
                    ctx.fillStyle = 'rgba(255, 255, 255, 0.8)';
                    ctx.fill();

                    // Main circle with variant color
                    ctx.beginPath();
                    ctx.arc(centerX, centerY, circleRadius, 0, Math.PI * 2);
                    ctx.fillStyle = BASE_COLORS[mm.alt] || '#ef4444';
                    ctx.fill();

                    // White border for visibility
                    ctx.strokeStyle = '#000';
                    ctx.lineWidth = 1;
                    ctx.stroke();
                }
            }
        }

        ctx.restore();
    }

    /**
     * Render soft clips as dashed rectangles extending beyond read boundaries.
     */
    private renderSoftClips(
        ctx: CanvasRenderingContext2D,
        read: Read,
        x1: number,
        x2: number,
        y: number,
        height: number,
        scale: number
    ): void {
        if (!this.state.settings.showSoftClips || !read.soft_clips?.length) return;

        ctx.save();
        ctx.fillStyle = SOFT_CLIP_STYLE.fillColor;
        ctx.strokeStyle = SOFT_CLIP_STYLE.strokeColor;
        ctx.setLineDash(SOFT_CLIP_STYLE.dashPattern);
        ctx.lineWidth = 1;

        for (const clip of read.soft_clips) {
            const clipWidth = clip.length * scale;

            if (clip.side === 'left') {
                // Draw before the read
                const clipX = x1 - clipWidth;
                ctx.fillRect(clipX, y, clipWidth, height);
                ctx.strokeRect(clipX, y, clipWidth, height);
            } else {
                // Draw after the read
                ctx.fillRect(x2, y, clipWidth, height);
                ctx.strokeRect(x2, y, clipWidth, height);
            }
        }

        ctx.setLineDash([]);
        ctx.restore();
    }

    public getReadAtPosition(event: MouseEvent): Read | null {
        const rect = this.readsCanvas.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;

        const { height: READ_HEIGHT, gap: READ_GAP } = this.getReadDimensions();
        const rowIndex = Math.floor(y / (READ_HEIGHT + READ_GAP));
        if (rowIndex < 0 || rowIndex >= this.state.packedRows.length) return null;

        const row = this.state.packedRows[rowIndex];
        const scale = this.getScale();

        // Binary search or linear scan (linear is fine for row length)
        for (const read of row) {
            const x1 = (read.position - this.state.viewport.start) * scale;
            const x2 = (read.end_position - this.state.viewport.start) * scale;

            if (x >= x1 && x <= x2) {
                return read;
            }
        }
        return null;
    }
}
