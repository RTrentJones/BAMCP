import { describe, expect, it } from "vitest";
import {
  baseColor,
  BASE_COLORS,
  displayPos,
  escapeHtml,
  EXPECTED_SCHEMA_VERSION,
  formatRulerLabel,
  toRegionString,
} from "./constants";

describe("displayPos — internal 0-based → user-facing 1-based", () => {
  it("adds one to convert the internal frame to the display frame", () => {
    expect(displayPos(0)).toBe(1);
    expect(displayPos(1050)).toBe(1051);
  });

  it("round-trips against the documented convention", () => {
    // Internal 0-based 7676153 is the 1-based 7676154 shown to users / sent to ClinVar.
    expect(displayPos(7676153)).toBe(7676154);
  });
});

describe("formatRulerLabel — ticks must stay distinct at every zoom", () => {
  it("keeps adjacent 50bp ticks distinct (the '1.1K' repeated bug)", () => {
    // 300bp span → 50bp tick interval. The old .toFixed(1)+'K' rendered 1051 and 1101
    // both as '1.1K'; now they must differ.
    const a = formatRulerLabel(1051, 50);
    const b = formatRulerLabel(1101, 50);
    expect(a).not.toBe(b);
    expect(a).toBe("1,051");
    expect(b).toBe("1,101");
  });

  it("keeps adjacent 5bp ticks distinct at base-level zoom", () => {
    expect(formatRulerLabel(1046, 5)).toBe("1,046");
    expect(formatRulerLabel(1051, 5)).toBe("1,051");
    expect(formatRulerLabel(1046, 5)).not.toBe(formatRulerLabel(1051, 5));
  });

  it("uses a K unit only when the tick spacing keeps K labels distinct", () => {
    // 5000bp interval → whole-K ticks.
    expect(formatRulerLabel(5000, 5000)).toBe("5K");
    expect(formatRulerLabel(10000, 5000)).toBe("10K");
    // 1000bp interval → still whole-K.
    expect(formatRulerLabel(2000, 1000)).toBe("2K");
  });

  it("uses an M unit at chromosome scale", () => {
    expect(formatRulerLabel(5_000_000, 1_000_000)).toBe("5M");
    expect(formatRulerLabel(2_500_000, 500_000)).not.toBe(
      formatRulerLabel(3_000_000, 500_000),
    );
  });

  it("never collapses two different tick positions to one label for nice intervals", () => {
    for (const interval of [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000]) {
      const base = 1_000_000;
      const first = formatRulerLabel(base, interval);
      const second = formatRulerLabel(base + interval, interval);
      expect(first).not.toBe(second);
    }
  });
});

describe("payload schema version", () => {
  it("is a positive integer the server can match against", () => {
    expect(Number.isInteger(EXPECTED_SCHEMA_VERSION)).toBe(true);
    expect(EXPECTED_SCHEMA_VERSION).toBeGreaterThan(0);
  });
});

describe("toRegionString — internal 0-based half-open → 1-based inclusive region", () => {
  it("renders a single base as start == end (samtools convention)", () => {
    // Internal [1050, 1051) is the one 1-based base 1051.
    expect(toRegionString("chr1", 1050, 1051)).toBe("chr1:1051-1051");
  });

  it("maps [0, 100) to the first hundred bases 1-100", () => {
    expect(toRegionString("chr1", 0, 100)).toBe("chr1:1-100");
  });

  it("covers the span exactly: inclusive length equals end0 - start0", () => {
    const region = toRegionString("chr2", 999, 2000);
    expect(region).toBe("chr2:1000-2000");
    const [, a, b] = region.match(/:(\d+)-(\d+)$/)!;
    expect(Number(b) - Number(a) + 1).toBe(2000 - 999);
  });

  it("widens fractional viewport edges outward", () => {
    // [10.4, 20.2) touches internal bases 10..20 → 1-based 11..21.
    expect(toRegionString("chr1", 10.4, 20.2)).toBe("chr1:11-21");
  });

  it("clamps a negative start to the first base", () => {
    expect(toRegionString("chr1", -50, 50)).toBe("chr1:1-50");
  });

  it("never emits an empty or inverted range", () => {
    expect(toRegionString("chr1", 100, 100)).toBe("chr1:101-101");
    expect(toRegionString("chr1", 100, 90)).toBe("chr1:101-101");
  });
});

describe("escapeHtml — BAM/VCF/server text is never parsed as markup", () => {
  it("escapes every HTML-significant character a SAM QNAME may carry", () => {
    const qname = `r1<img src=x onerror="alert('x')">&`;
    const out = escapeHtml(qname);
    expect(out).toBe(
      "r1&lt;img src=x onerror=&quot;alert(&#39;x&#39;)&quot;&gt;&amp;",
    );
    expect(out).not.toMatch(/[<>"']/);
  });

  it("stringifies non-string payload values and blanks null/undefined", () => {
    expect(escapeHtml(42)).toBe("42");
    expect(escapeHtml(null)).toBe("");
    expect(escapeHtml(undefined)).toBe("");
  });

  it("leaves a symbolic SV alt readable once rendered", () => {
    expect(escapeHtml("<DEL>")).toBe("&lt;DEL&gt;");
  });
});

describe("baseColor — nucleotide lookup that ignores prototype members", () => {
  it("returns the palette color for a base", () => {
    expect(baseColor("A", "#000")).toBe(BASE_COLORS.A);
  });

  it("falls back for unknown or prototype-named keys", () => {
    expect(baseColor("constructor", "#000")).toBe("#000");
    expect(baseColor("<x>", "#111")).toBe("#111");
  });
});
