/**
 * Shared constants for the BAMCP viewer.
 * Includes color palettes for IGV-style visualization options.
 */

// RegionData payload contract version this viewer build expects. The server stamps
// `schema_version` on every payload; a mismatch means the bundle is out of step with the
// server and may mis-decode, so we warn loudly (see BAMCPClient.handleData).
export const EXPECTED_SCHEMA_VERSION = 1;

// Base colors for nucleotides
export const BASE_COLORS: Record<string, string> = {
    'A': '#22c55e',  // Green
    'T': '#ef4444',  // Red
    'G': '#f97316',  // Orange
    'C': '#3b82f6',  // Blue
    'N': '#9ca3af',  // Gray
};

// Color palettes for different colorBy modes
export const COLOR_PALETTES = {
    strand: {
        forward: '#3b82f6',  // Blue
        reverse: '#ef4444',  // Red (was purple #818cf8)
    },
    mapq: [
        { threshold: 0, color: '#f87171' },   // 0-10: light red
        { threshold: 10, color: '#fbbf24' },  // 10-20: yellow
        { threshold: 20, color: '#a3e635' },  // 20-30: lime
        { threshold: 30, color: '#22c55e' },  // 30-40: green
        { threshold: 40, color: '#14b8a6' },  // 40-50: teal
        { threshold: 50, color: '#0ea5e9' },  // 50-60: sky blue
        { threshold: 60, color: '#6366f1' },  // 60+: indigo
    ],
    insertSize: {
        normal: '#94a3b8',      // Slate - normal (150-500bp)
        small: '#f97316',       // Orange - too short (<100bp)
        large: '#ef4444',       // Red - too long (>1000bp)
        veryLarge: '#7c3aed',   // Purple - very long (>5000bp, likely SV)
        chimeric: '#ec4899',    // Pink - different chromosome
    },
    baseQuality: [
        { threshold: 0, color: '#f87171' },   // 0-10: red
        { threshold: 10, color: '#fbbf24' },  // 10-20: yellow
        { threshold: 20, color: '#a3e635' },  // 20-30: lime
        { threshold: 30, color: '#22c55e' },  // 30+: green
    ],
};

// Display mode configurations
export const DISPLAY_MODE_CONFIGS = {
    squished: { readHeight: 6, readGap: 1, showLabels: false },
    compact: { readHeight: 12, readGap: 2, showLabels: false },
    expanded: { readHeight: 24, readGap: 4, showLabels: true },
    // DeepVariant-inspired pileup modes. Strip mode allocates 6 channels
    // vertically inside each read row; composite mode false-colors a single row.
    'dv-strips': { readHeight: 36, readGap: 4, showLabels: false },
    'dv-composite': { readHeight: 8, readGap: 1, showLabels: false },
};

// Number of channels in DeepVariant-style rendering (base, base_qual, mapq,
// strand, supports_variant, ref_differs).
export const DV_CHANNEL_COUNT = 6;

// Per-channel colors used by the dv-strips rendering for the boolean/categorical
// channels. Quantitative channels (base quality, MAPQ) reuse the threshold
// palettes above.
export const DV_PALETTE = {
    strandForward: '#3b82f6',
    strandReverse: '#ef4444',
    variantSupport: '#22c55e',   // green — read base matches alt allele
    variantReference: '#9ca3af', // gray — read base matches ref allele
    variantNoCall: '#1f2937',    // dark — position outside active variant
    refMatch: '#0f172a',         // near-black where read matches reference
    refDiffers: '#06b6d4',       // cyan where read mismatches reference
};

// Insert size thresholds for color coding
export const INSERT_SIZE_THRESHOLDS = {
    tooShort: 100,
    normalMin: 150,
    normalMax: 500,
    tooLong: 1000,
    veryLong: 5000,
};

// Soft clip styling
export const SOFT_CLIP_STYLE = {
    fillColor: 'rgba(156, 163, 175, 0.4)',  // Semi-transparent gray
    strokeColor: '#6b7280',                  // Gray border
    dashPattern: [2, 2],
};

/**
 * Convert an internal pysam 0-based genomic position to the 1-based value shown
 * to humans and sent to external databases (the VCF / dbSNP / ClinVar / gnomAD /
 * IGV convention). Use ONLY for displayed text (ruler labels, tables, tooltips)
 * and external-lookup coordinates — NEVER for canvas placement, viewport math,
 * or evidence/mate index keys, which all stay in the raw 0-based frame.
 */
export function displayPos(pos: number): number {
    return pos + 1;
}

/**
 * Format a genomic position for a ruler tick, given the spacing between ticks.
 *
 * A fixed `.toFixed(1) + 'K'` rounding collapses every tick finer than ~100 bp to the
 * same string (a 50 bp interval renders "1.1K, 1.1K, 1.2K, 1.2K…"; a 5 bp interval
 * renders "1.1K" five times), so the ruler conveys no information when zoomed in. Pick
 * a K/M unit only when the tick spacing is coarse enough for that unit to stay distinct,
 * with enough decimals that adjacent ticks differ; otherwise use the full comma-grouped
 * position. `tickInterval` is in base pairs.
 */
export function formatRulerLabel(pos: number, tickInterval: number): string {
    const withUnit = (divisor: number, unit: string): string => {
        const step = tickInterval / divisor;
        const decimals = Math.max(0, Math.ceil(-Math.log10(step)));
        return (pos / divisor).toFixed(decimals) + unit;
    };
    if (pos >= 1_000_000 && tickInterval >= 1_000_000) return withUnit(1_000_000, "M");
    if (pos >= 1_000 && tickInterval >= 1_000) return withUnit(1_000, "K");
    return Math.round(pos).toLocaleString("en-US");
}

/** Color per structural-variant type; falls back to SV_COLOR_DEFAULT for unknown types. */
export const SV_COLORS: Record<string, string> = {
    DEL: '#ef4444', // deletion — red
    DUP: '#3b82f6', // duplication — blue
    INV: '#a855f7', // inversion — purple
    INS: '#f97316', // insertion — orange
    CNV: '#14b8a6', // copy-number — teal
    BND: '#6b7280', // breakend / translocation — gray
};
export const SV_COLOR_DEFAULT = '#6b7280';

/** Short SV type label: the VCF SVTYPE, else the symbolic alt (`<DEL>` → `DEL`), else `SV`. */
export function svTypeLabel(svType?: string | null, alt?: string): string {
    if (svType) return svType;
    if (alt?.startsWith('<')) return alt.replace(/[<>]/g, '');
    return 'SV';
}

/** Resolve a rendering color for a structural-variant type (case-insensitive). */
export function svColor(svType: string): string {
    return SV_COLORS[svType.toUpperCase()] ?? SV_COLOR_DEFAULT;
}

/**
 * Escape a string for safe interpolation into innerHTML. VCF-derived fields such as
 * a symbolic SV alt (`<DEL>`, `<DUP>`) contain angle brackets that would otherwise be
 * parsed as markup, so escape before injecting them into any HTML template. Also
 * applies to BAM-derived text: a SAM QNAME may legally contain `<`, `>`, `"`, `'`, `&`.
 * Non-string values are stringified first (null/undefined become empty), so numeric
 * fields from an untrusted payload can be passed through it too.
 */
export function escapeHtml(value: unknown): string {
    return String(value ?? '')
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#39;');
}

/**
 * Build a samtools-style region string (1-based, inclusive on both ends) from the
 * viewer's internal 0-based half-open interval [start0, end0). This is the contract
 * every tool's `region` argument speaks: `chr1:1051-1051` is exactly one base (the
 * internal 0-based position 1050). Fractional viewport edges are widened outward so
 * the whole visible span is covered, a negative start clamps to the first base, and
 * an empty/inverted interval still yields a single base rather than an invalid range.
 */
export function toRegionString(contig: string, start0: number, end0: number): string {
    const s = Math.max(0, Math.floor(start0));
    const e = Math.max(s + 1, Math.ceil(end0));
    return `${contig}:${s + 1}-${e}`;
}

/** Color for a single nucleotide, or `fallback` — never a prototype member like `constructor`. */
export function baseColor(base: string, fallback: string): string {
    return Object.prototype.hasOwnProperty.call(BASE_COLORS, base) ? BASE_COLORS[base] : fallback;
}

/**
 * Largest value in `values`, or `floor` if none is larger (NaN entries are ignored).
 * A loop, not `Math.max(...values)`: spreading passes every element as a call
 * argument, which throws RangeError once an array outgrows the engine's argument
 * limit (~125k in V8) — i.e. any coverage array for a region over ~125 kb.
 */
export function maxOf(values: ArrayLike<number> | null | undefined, floor = -Infinity): number {
    let max = floor;
    if (!values) return max;
    for (let i = 0; i < values.length; i++) {
        const v = values[i];
        if (v > max) max = v;
    }
    return max;
}

// Server-side cap on a single region request (bamcp.core.parsers.MAX_REGION_SIZE).
export const MAX_REGION_SIZE = 1_000_000;
// Navigation limits. Zooming in below ~20 bp shows nothing useful (and a sub-1 bp
// span divides the scale by ~0); zooming out is capped at half the server limit so
// the padded refetch window (viewport + 50% each side) still fits in one request.
export const MIN_VIEW_SPAN = 20;
export const MAX_VIEW_SPAN = MAX_REGION_SIZE / 2;

export interface Span {
    start: number;
    end: number;
}

/**
 * Zoom `vp` about its center by `factor` (<1 zooms in), clamping the span to
 * [minSpan, maxSpan] and the start to >= 0. A view already outside the limits
 * (e.g. a host-loaded 1 Mb region) is never pushed further out of them, but is
 * not forcibly snapped either: zooming in from it still works.
 */
export function zoomViewport(
    vp: Span, factor: number, minSpan = MIN_VIEW_SPAN, maxSpan = MAX_VIEW_SPAN,
): Span {
    const span = Math.max(vp.end - vp.start, 1);
    const center = (vp.start + vp.end) / 2;
    const lo = Math.min(minSpan, span);
    const hi = Math.max(maxSpan, span);
    const newSpan = Math.min(hi, Math.max(lo, span * factor));
    const start = Math.max(0, center - newSpan / 2);
    return { start, end: start + newSpan };
}

/** Shift `vp` by `deltaBp`, never past the start of the contig. */
export function panViewport(vp: Span, deltaBp: number): Span {
    const delta = Math.max(deltaBp, -vp.start);
    return { start: vp.start + delta, end: vp.end + delta };
}

/**
 * Viewport-refetch window: the view padded by `padFrac` of its span on each side,
 * shrunk so the resulting 1-based inclusive request (see toRegionString) stays
 * within `maxSize` bases, and never starting before 0.
 */
export function paddedFetchWindow(
    vp: Span, padFrac = 0.5, maxSize = MAX_REGION_SIZE,
): Span {
    const budget = maxSize - 2; // floor/ceil of fractional edges can add a base each
    const span = vp.end - vp.start;
    if (span >= budget) {
        const start = Math.max(0, vp.start);
        return { start, end: start + budget };
    }
    const pad = Math.max(0, Math.min(span * padFrac, (budget - span) / 2));
    const start = Math.max(0, vp.start - pad);
    return { start, end: Math.min(start + budget, vp.end + pad) };
}

/**
 * Read-position line for the evidence panel. The bias call belongs to the server
 * (it flags `read_position_bias` only when supporting bases sit within N bp of an
 * aligned read end beyond what uniform placement predicts), so prefer its
 * near-end fraction and only say "Warning" when the server raised that risk.
 * Older payloads without the fraction get a neutral histogram description that
 * makes no bias claim of its own.
 */
export function formatNearEndSummary(evidence: {
    near_end_fraction?: number;
    near_end_cutoff_bp?: number;
    near_end_expected_fraction?: number | null;
    position_histogram?: number[];
    artifact_risk?: { risks?: Array<{ type?: string }> };
}): string {
    const flagged = Array.isArray(evidence.artifact_risk?.risks) &&
        evidence.artifact_risk!.risks!.some(r => r?.type === 'read_position_bias');
    const frac = Number(evidence.near_end_fraction);
    if (evidence.near_end_fraction != null && Number.isFinite(frac)) {
        const cutoff = Number(evidence.near_end_cutoff_bp);
        const within = Number.isFinite(cutoff) ? `within ${cutoff}bp of a read end` : 'near a read end';
        let text = `${Math.round(frac * 100)}% ${within}`;
        const expected = Number(evidence.near_end_expected_fraction);
        if (evidence.near_end_expected_fraction != null && Number.isFinite(expected)) {
            text += ` (~${Math.round(expected * 100)}% expected)`;
        }
        return flagged ? `Warning: ${text}` : text;
    }

    const histogram = evidence.position_histogram;
    if (!Array.isArray(histogram) || histogram.length === 0) return '';
    let total = 0;
    for (const v of histogram) total += Number(v) || 0;
    if (total === 0) return '';
    const firstBins = (Number(histogram[0]) || 0) + (Number(histogram[1]) || 0);
    const text = `${Math.round((firstBins / total) * 100)}% in position bins 0-50`;
    return flagged ? `Warning: ${text}` : text;
}
