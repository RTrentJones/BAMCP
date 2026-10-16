import { describe, expect, it } from "vitest";
import { layoutRead, parseCigar } from "./layout";

describe("parseCigar", () => {
  it("parses a multi-op CIGAR", () => {
    expect(parseCigar("5S10M2I3D")).toEqual([
      { op: "S", len: 5 },
      { op: "M", len: 10 },
      { op: "I", len: 2 },
      { op: "D", len: 3 },
    ]);
  });

  it("returns [] for empty, '*' or malformed input", () => {
    expect(parseCigar("")).toEqual([]);
    expect(parseCigar("*")).toEqual([]);
    expect(parseCigar("10M?3D")).toEqual([]);
    expect(parseCigar(undefined)).toEqual([]);
  });
});

describe("layoutRead — CIGAR-aware base placement", () => {
  it("skips soft-clipped query bases so the first aligned base sits at position", () => {
    const l = layoutRead({ cigar: "3S4M", position: 100, end_position: 104 });
    expect(l.bases[0]).toEqual({ refPos: 100, qpos: 3 });
    expect(l.bases).toHaveLength(4);
  });

  it("records insertions at the junction and advances only the query", () => {
    const l = layoutRead({ cigar: "2M2I2M", position: 10, end_position: 14 });
    expect(l.insertions).toEqual([{ refPos: 12, len: 2 }]);
    expect(l.bases.map(b => [b.refPos, b.qpos])).toEqual([[10, 0], [11, 1], [12, 4], [13, 5]]);
  });

  it("records deletions and N skips as reference spans", () => {
    const l = layoutRead({ cigar: "2M3D1M2N1M", position: 0, end_position: 9 });
    expect(l.deletions).toEqual([
      { start: 2, end: 5, skip: false },
      { start: 6, end: 8, skip: true },
    ]);
    expect(l.bases.map(b => b.refPos)).toEqual([0, 1, 5, 8]);
  });

  it("ignores hard clips", () => {
    const l = layoutRead({ cigar: "5H3M", position: 50, end_position: 53 });
    expect(l.bases[0]).toEqual({ refPos: 50, qpos: 0 });
  });

  it("falls back to a gapless block when the CIGAR is unusable", () => {
    const l = layoutRead({ cigar: "", position: 7, end_position: 10 });
    expect(l.bases.map(b => b.refPos)).toEqual([7, 8, 9]);
    expect(l.insertions).toEqual([]);
    expect(l.deletions).toEqual([]);
  });
});
