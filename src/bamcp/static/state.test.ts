import { describe, expect, it } from "vitest";
import { DataStore } from "./data-store";
import { StateManager, decodeReads, decodeRegionData } from "./state";
import { RegionData } from "./types";

const base = { contig: "chr1", start: 100, end: 200 };

describe("loadData defensiveness", () => {
  it("does not throw when variants is missing", () => {
    const s = new StateManager();
    const data = { ...base, reads: [], coverage: [1, 2, 3] } as unknown as RegionData;
    expect(() => s.loadData(data)).not.toThrow();
    expect(s.data?.variants).toEqual([]);
    expect(s.getFilteredAndSortedVariants()).toEqual([]);
  });

  it("does not throw when reads, coverage and variants are all missing or null", () => {
    const s = new StateManager();
    const data = { ...base, reads: null, coverage: null, variants: null } as unknown as RegionData;
    expect(() => s.loadData(data)).not.toThrow();
    expect(s.data?.reads).toEqual([]);
    expect(s.data?.coverage).toEqual([]);
    expect(s.packedRows).toEqual([]);
    expect(s.viewport).toEqual({ start: 100, end: 200 });
  });

  it("tolerates a non-array variants value and junk entries", () => {
    const s = new StateManager();
    s.loadData({ ...base, reads: [], coverage: [], variants: "oops" } as unknown as RegionData);
    expect(s.data?.variants).toEqual([]);
    s.loadData({ ...base, reads: [], coverage: [], variants: [null, 3] } as unknown as RegionData);
    expect(s.data?.variants).toEqual([]);
  });

  it("loadTile and activateTile are defensive too", () => {
    const s = new StateManager();
    const data = { ...base, reads: [] } as unknown as RegionData;
    expect(() => s.loadTile(data)).not.toThrow();
    expect(() => s.activateTile(data)).not.toThrow();
  });

  it("falls back to a sane span when start/end are missing", () => {
    const d = decodeRegionData({ contig: "chr1", coverage: [0, 0, 0, 0] } as unknown as RegionData);
    expect(d.start).toBe(0);
    expect(d.end).toBe(4);
  });

  it("returns well-formed data unchanged", () => {
    const d = { ...base, reads: [], coverage: [], variants: [] } as unknown as RegionData;
    expect(decodeRegionData(d)).toBe(d);
  });

  it("DataStore.ingest tolerates missing variants directly", () => {
    const store = new DataStore();
    expect(() => store.ingest({ ...base } as unknown as RegionData)).not.toThrow();
    expect(store.getAllVariants()).toEqual([]);
  });
});

describe("decodeReads defensiveness", () => {
  it("fills missing mismatches on object-form reads", () => {
    const reads = decodeReads([
      { name: "a", cigar: "5M", position: 1, end_position: 6, mapping_quality: 60, is_reverse: false },
    ] as never);
    expect(reads[0].mismatches).toEqual([]);
  });

  it("decodes columnar reads with a missing mismatches column and no count", () => {
    const reads = decodeReads({
      name: ["a", "b"], cigar: ["5M", "5M"], position: [1, 2], end_position: [6, 7],
      mapping_quality: [60, 0], is_reverse: [false, true],
    } as never);
    expect(reads).toHaveLength(2);
    expect(reads[1].mismatches).toEqual([]);
    expect(reads[1].is_reverse).toBe(true);
  });
});
