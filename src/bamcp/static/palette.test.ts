import { describe, expect, it } from "vitest";
import { BASE_COLORS, COLOR_PALETTES } from "./constants";
import {
  LOW_BASE_QUALITY,
  LOW_BQ_WASH,
  MATCH_LETTER,
  MIN_TICK_CONTRAST,
  READ_BODY,
  READ_EDGE,
  contrastRatio,
  fadeToWhite,
  isLowBaseQuality,
  mapqOpacity,
  mixColor,
  readBodyFill,
  relativeLuminance,
  tickColor,
} from "./palette";

const BASES = ["A", "C", "G", "T", "N"];
const STRANDS = [false, true];
// Every MAPQ the fade distinguishes (0 clamps to the 0.2 floor) plus "unavailable".
const MAPQS = [0, 6, 12, 20, 30, 40, 50, 60, 255];

describe("colour maths", () => {
  it("matches WCAG reference values", () => {
    expect(relativeLuminance("#ffffff")).toBeCloseTo(1, 5);
    expect(relativeLuminance("#000000")).toBeCloseTo(0, 5);
    expect(contrastRatio("#000000", "#ffffff")).toBeCloseTo(21, 3);
    expect(contrastRatio("#777777", "#777777")).toBeCloseTo(1, 5);
  });

  it("mixes and fades toward white", () => {
    expect(mixColor("#000000", "#ffffff", 0.5)).toBe("#808080");
    expect(fadeToWhite("#2563eb", 1)).toBe("#2563eb");
    expect(fadeToWhite("#2563eb", 0)).toBe("#ffffff");
  });

  it("flags low base quality only when known and below the threshold", () => {
    expect(isLowBaseQuality(LOW_BASE_QUALITY - 1)).toBe(true);
    expect(isLowBaseQuality(LOW_BASE_QUALITY)).toBe(false);
    expect(isLowBaseQuality(undefined)).toBe(false);
    expect(isLowBaseQuality(null)).toBe(false);
  });
});

describe("read-track legibility (base x strand x MAPQ fade)", () => {
  it("never colours the read body by strand (the old blue-on-blue / red-on-red collision)", () => {
    for (const mapq of MAPQS) {
      const a = mapqOpacity(mapq);
      expect(readBodyFill(false, a)).toBe(readBodyFill(true, a));
    }
  });

  for (const base of BASES) {
    for (const reverse of STRANDS) {
      it(`${base} tick stays >= ${MIN_TICK_CONTRAST}:1 on a ${reverse ? "reverse" : "forward"} read at every MAPQ`, () => {
        for (const mapq of MAPQS) {
          const body = readBodyFill(reverse, mapqOpacity(mapq));
          // Solid (high-BQ) tick, and the border of a hollow (low-BQ) tick: same colour.
          const ratio = contrastRatio(tickColor(base), body);
          expect(ratio, `${base} on ${body} (MAPQ ${mapq})`).toBeGreaterThanOrEqual(MIN_TICK_CONTRAST);
        }
        // MAPQ-0 reads have a white body.
        expect(contrastRatio(tickColor(base), "#ffffff")).toBeGreaterThanOrEqual(MIN_TICK_CONTRAST);
      });
    }
  }

  it("keeps the strand stripe distinct from the body for both strands", () => {
    const body = readBodyFill(false, 1);
    expect(contrastRatio(COLOR_PALETTES.strand.forward, body)).toBeGreaterThanOrEqual(1.5);
    expect(contrastRatio(COLOR_PALETTES.strand.reverse, body)).toBeGreaterThanOrEqual(1.5);
  });

  it("keeps letters readable at base zoom", () => {
    for (const base of BASES) {
      // White letter on a solid mismatch tile.
      expect(contrastRatio("#ffffff", BASE_COLORS[base])).toBeGreaterThanOrEqual(4.5);
      // Coloured letter on a hollow (low-BQ) tile's pale centre.
      const pale = mixColor(BASE_COLORS[base], "#ffffff", LOW_BQ_WASH);
      expect(contrastRatio(BASE_COLORS[base], pale)).toBeGreaterThanOrEqual(MIN_TICK_CONTRAST);
    }
    expect(contrastRatio(MATCH_LETTER, READ_BODY)).toBeGreaterThanOrEqual(4.5);
  });

  it("makes a low-BQ (hollow) tick visibly different from a solid one", () => {
    for (const base of BASES) {
      const pale = mixColor(BASE_COLORS[base], "#ffffff", LOW_BQ_WASH);
      expect(contrastRatio(pale, BASE_COLORS[base])).toBeGreaterThanOrEqual(MIN_TICK_CONTRAST);
    }
  });

  it("gives the read edge a visible outline on white at full MAPQ", () => {
    expect(contrastRatio(READ_EDGE, "#ffffff")).toBeGreaterThanOrEqual(2);
  });
});
