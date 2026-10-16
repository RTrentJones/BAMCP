/**
 * Pure helpers that derive viewer chrome (banner text, hidden-read counts) from
 * state. Kept DOM-free so vitest can cover them directly.
 */

import { Read, RegionData } from "./types";

export const NO_REFERENCE_MESSAGE =
    'No reference genome: variant calls and mismatches are hidden. ' +
    'Pass a reference or use demo:HG002.';

/**
 * Banner text when the payload says SNV detection ran without a reference, or
 * null when no banner is warranted. Every field is optional: older servers send
 * no `applied_filters`, and a missing field must never show (or crash) a banner.
 *
 * - `applied_filters.snv_detection` = "disabled: no reference" → banner.
 * - "disabled: variant_source=vcf …" is a deliberate VCF-authoritative view,
 *   not a missing reference → no banner.
 * - `applied_filters.reference_source` (or a top-level `reference_source`) of
 *   "none" → banner; any other value is appended as context when a banner shows.
 */
export function referenceBannerMessage(data: Partial<RegionData> | null | undefined): string | null {
    if (!data || typeof data !== 'object') return null;
    const filters = (data.applied_filters && typeof data.applied_filters === 'object')
        ? data.applied_filters as Record<string, unknown>
        : {};
    const snv = typeof filters.snv_detection === 'string' ? filters.snv_detection : '';
    const rawSource = filters.reference_source ?? data.reference_source;
    const source = typeof rawSource === 'string' ? rawSource : '';

    // A VCF-authoritative view deliberately skips BAM SNV calling; its variants are shown,
    // so a missing reference there is not "calls hidden" — never banner it.
    if (/^disabled/i.test(snv) && /variant_source\s*=\s*vcf/i.test(snv)) return null;
    const snvNoRef = /^disabled/i.test(snv) && /reference/i.test(snv);
    const sourceNone = /^none$/i.test(source);
    if (!snvNoRef && !sourceNone) return null;
    return NO_REFERENCE_MESSAGE;
}

/**
 * Index of the first packed row that is not fully visible inside a reads area
 * `visibleHeightPx` tall (rows are `rowPitch` apart and `rowHeight` tall), also
 * capped by `drawnRows` (rows past the canvas cap are never drawn at all).
 */
export function firstHiddenRow(
    visibleHeightPx: number,
    rowPitch: number,
    rowHeight: number,
    drawnRows: number,
): number {
    if (!(rowPitch > 0)) return drawnRows;
    const h = Math.max(0, visibleHeightPx);
    // Row r is fully visible when r * pitch + rowHeight <= h.
    const fullyVisible = h < rowHeight ? 0 : Math.floor((h - rowHeight) / rowPitch) + 1;
    return Math.max(0, Math.min(fullyVisible, drawnRows));
}

/**
 * Reads in rows >= `fromRow` that overlap the genomic viewport — i.e. reads the
 * user would see if the reads area were tall enough, but currently can't.
 */
export function hiddenReadCount(
    packedRows: Read[][],
    viewport: { start: number; end: number },
    fromRow: number,
): number {
    let n = 0;
    for (let r = Math.max(0, fromRow); r < packedRows.length; r++) {
        for (const read of packedRows[r]) {
            if (read.end_position > viewport.start && read.position < viewport.end) n++;
        }
    }
    return n;
}

/** Human text for the hidden-reads indicator. */
export function hiddenReadsLabel(count: number): string {
    const noun = count === 1 ? 'read' : 'reads';
    return `▾ ${count.toLocaleString()} ${noun} not shown — scroll down or zoom in`;
}
