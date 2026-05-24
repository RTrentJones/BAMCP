/**
 * Shared constants for the BAMCP viewer.
 * Includes color palettes for IGV-style visualization options.
 */

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
