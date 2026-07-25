import { describe, expect, it } from "vitest";
import { displayPos, EXPECTED_SCHEMA_VERSION, formatRulerLabel } from "./constants";

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
