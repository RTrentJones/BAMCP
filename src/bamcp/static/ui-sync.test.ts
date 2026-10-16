import { describe, expect, it } from "vitest";
import {
  NO_REFERENCE_MESSAGE,
  firstHiddenRow,
  hiddenReadCount,
  hiddenReadsLabel,
  referenceBannerMessage,
} from "./ui-sync";
import { Read } from "./types";

describe("referenceBannerMessage", () => {
  it("shows the banner when SNV detection is disabled for lack of a reference", () => {
    expect(referenceBannerMessage({ applied_filters: { snv_detection: "disabled: no reference" } }))
      .toBe(NO_REFERENCE_MESSAGE);
  });

  it("names the fix (pass a reference or use the demo)", () => {
    expect(NO_REFERENCE_MESSAGE).toMatch(/No reference genome/);
    expect(NO_REFERENCE_MESSAGE).toMatch(/demo:HG002/);
  });

  it("shows the banner when reference_source is 'none'", () => {
    expect(referenceBannerMessage({ applied_filters: { reference_source: "none" } })).toBe(NO_REFERENCE_MESSAGE);
    expect(referenceBannerMessage({ reference_source: "none" })).toBe(NO_REFERENCE_MESSAGE);
  });

  it("stays hidden when detection is enabled or the reference was auto-resolved", () => {
    expect(referenceBannerMessage({ applied_filters: { snv_detection: "enabled" } })).toBeNull();
    expect(referenceBannerMessage({
      applied_filters: { snv_detection: "enabled", reference_source: "auto" },
    })).toBeNull();
  });

  it("does not treat a VCF-authoritative view as a missing reference", () => {
    expect(referenceBannerMessage({
      applied_filters: { snv_detection: "disabled: variant_source=vcf (the VCF is authoritative)" },
    })).toBeNull();
  });

  it("does not banner a VCF-authoritative view even when no reference was resolved", () => {
    expect(referenceBannerMessage({
      applied_filters: {
        snv_detection: "disabled: variant_source=vcf (the VCF is authoritative)",
        reference_source: "none",
      },
    })).toBeNull();
  });

  it("copes with absent or malformed fields", () => {
    expect(referenceBannerMessage(null)).toBeNull();
    expect(referenceBannerMessage(undefined)).toBeNull();
    expect(referenceBannerMessage({})).toBeNull();
    expect(referenceBannerMessage({ applied_filters: undefined })).toBeNull();
    expect(referenceBannerMessage({ applied_filters: { snv_detection: 42 } } as never)).toBeNull();
  });
});

function read(position: number, end: number): Read {
  return {
    name: `r${position}`, cigar: `${end - position}M`, position, end_position: end,
    mapping_quality: 60, is_reverse: false, mismatches: [],
  };
}

describe("hidden-reads accounting", () => {
  it("firstHiddenRow counts only fully visible rows", () => {
    // pitch 14 (12 + 2 gap): rows end at 12, 26, 40 …
    expect(firstHiddenRow(40, 14, 12, 100)).toBe(3);
    expect(firstHiddenRow(39, 14, 12, 100)).toBe(2);
    expect(firstHiddenRow(5, 14, 12, 100)).toBe(0);
    expect(firstHiddenRow(-10, 14, 12, 100)).toBe(0);
  });

  it("firstHiddenRow is capped by the rows actually drawn", () => {
    expect(firstHiddenRow(10_000, 14, 12, 7)).toBe(7);
  });

  it("hiddenReadCount counts reads in hidden rows that overlap the viewport", () => {
    const rows = [[read(0, 50)], [read(0, 50), read(60, 90)], [read(200, 250)]];
    const vp = { start: 10, end: 100 };
    expect(hiddenReadCount(rows, vp, 0)).toBe(3);
    expect(hiddenReadCount(rows, vp, 1)).toBe(2);
    expect(hiddenReadCount(rows, vp, 3)).toBe(0);
  });

  it("labels the count", () => {
    expect(hiddenReadsLabel(1)).toMatch(/1 read not shown/);
    expect(hiddenReadsLabel(1234)).toMatch(/1,234 reads not shown/);
  });
});
