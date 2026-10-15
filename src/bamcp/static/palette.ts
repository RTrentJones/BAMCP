/**
 * Read-track palette and the colour maths behind it. Pure (no DOM) so the
 * legibility guarantees are unit-tested (palette.test.ts).
 *
 * Why read bodies are neutral grey: mismatch ticks are coloured by base (A green,
 * C blue, G orange, T red). When read bodies were coloured by strand (forward
 * blue, reverse red) a C on a forward read and a T on a reverse read vanished
 * into the body, so a homozygous-alt C looked like "alt only on reverse reads"
 * (a strand-bias artifact). Bodies are now one light grey for both strands; the
 * strand is carried by a coloured bottom stripe and a coloured 3'-end tip, and
 * every base colour keeps >= 3:1 contrast against the body at every MAPQ fade.
 */

import { BASE_COLORS } from "./constants";

/** Read body fill (both strands) at full mapping quality. */
export const READ_BODY = '#d1d5db';
/** 1px edge around read bodies so faded (low-MAPQ) reads still read as boxes. */
export const READ_EDGE = '#9ca3af';
/** Letter colour on a reference-matching base tile at base zoom. */
export const MATCH_LETTER = '#374151';
/** MAPQ-0 reads: white fill with a dashed outline in this colour. */
export const MAPQ_ZERO_EDGE = '#6b7280';

/** Base qualities strictly below this draw "hollow" (low confidence). */
export const LOW_BASE_QUALITY = 20;
/** How far a low-BQ tick/tile's centre is washed toward white (0..1). */
export const LOW_BQ_WASH = 0.8;

/** Minimum contrast every mismatch colour must keep against any read body. */
export const MIN_TICK_CONTRAST = 3;

/** Opacity-equivalent for a MAPQ: linear 0..60 → 0.2..1, clamped (255 = unavailable). */
export function mapqOpacity(mapq: number): number {
    const v = Number(mapq);
    if (!Number.isFinite(v) || v >= 255) return 1;
    return Math.min(1, Math.max(v / 60, 0.2));
}

function channels(hex: string): [number, number, number] {
    const h = hex.replace('#', '');
    return [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16)) as [number, number, number];
}

/** Mix `a` toward `b` by `t` (0 → a, 1 → b); both `#rrggbb`. */
export function mixColor(a: string, b: string, t: number): string {
    const k = Math.min(1, Math.max(0, t));
    const ca = channels(a);
    const cb = channels(b);
    return '#' + ca.map((v, i) => Math.round(v * (1 - k) + cb[i] * k).toString(16).padStart(2, '0')).join('');
}

/**
 * Paint `color` as if drawn at `opacity` over white, but as an opaque colour, so
 * the result does not depend on the page theme behind the canvas and the
 * contrast numbers below hold in light and dark mode alike.
 */
export function fadeToWhite(color: string, opacity: number): string {
    return mixColor(color, '#ffffff', 1 - opacity);
}

/** WCAG relative luminance of a `#rrggbb` colour. */
export function relativeLuminance(hex: string): number {
    const [r, g, b] = channels(hex).map(v => {
        const c = v / 255;
        return c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
    });
    return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

/** WCAG contrast ratio between two `#rrggbb` colours (1..21). */
export function contrastRatio(a: string, b: string): number {
    const la = relativeLuminance(a);
    const lb = relativeLuminance(b);
    return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
}

/**
 * Opaque read-body fill for a read drawn at `opacity` (MAPQ fade, or hover dim).
 * Identical for both strands by design — see the module comment.
 */
export function readBodyFill(_isReverse: boolean, opacity: number): string {
    return fadeToWhite(READ_BODY, opacity);
}

/** Mismatch / base colour; unknown symbols fall back to the N grey. */
export function tickColor(base: string): string {
    return Object.prototype.hasOwnProperty.call(BASE_COLORS, base) ? BASE_COLORS[base] : BASE_COLORS.N;
}

/** True when a base quality is known and below LOW_BASE_QUALITY. */
export function isLowBaseQuality(qual: number | null | undefined): boolean {
    return qual != null && Number.isFinite(qual) && qual < LOW_BASE_QUALITY;
}
