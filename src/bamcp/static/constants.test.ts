import { describe, expect, it } from "vitest";
import { displayPos, EXPECTED_SCHEMA_VERSION } from "./constants";

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

describe("payload schema version", () => {
  it("is a positive integer the server can match against", () => {
    expect(Number.isInteger(EXPECTED_SCHEMA_VERSION)).toBe(true);
    expect(EXPECTED_SCHEMA_VERSION).toBeGreaterThan(0);
  });
});
