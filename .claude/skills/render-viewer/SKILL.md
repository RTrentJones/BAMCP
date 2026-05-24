---
name: render-viewer
description: Render the BAMCP viewer to a PNG so you can visually verify a TS/HTML change without a browser. Use this immediately after editing any file under `src/bamcp/static/` to confirm your change does what you intended.
---

# render-viewer

Closes the visual feedback loop on the BAMCP viewer. Edit TS → rebuild → render → `Read` the PNG → iterate. No human-in-the-loop browser needed.

## When to use this skill

Trigger it after **any** of:

- Editing `src/bamcp/static/renderer.ts` — read the PNG to confirm rendering changed as intended.
- Editing `src/bamcp/static/state.ts`, `types.ts`, `constants.ts`, or `mcp-app.ts` and the change could affect output.
- Adding or modifying a display mode (the DV modes especially benefit from visual checks).
- Investigating a viewer bug — render the failing case to see what the user sees.

Skip it for pure-Python changes that don't touch the static directory.

## Workflow

1. **Rebuild the viewer.** `npm run build` regenerates `dist/viewer.html`, which is what the renderer serves. The script warns when source TS files are newer than the dist bundle:

   ```bash
   cd src/bamcp/static && npm run build && cd -
   ```

   If the build emits TS errors, fix those first — they're real type issues the renderer would crash on.

2. **Render the modes you care about.** For renderer.ts changes, default to rendering all three rendering modes side-by-side:

   ```bash
   make render-viewer-all-modes \
       BAM=tests/fixtures/comprehensive.bam \
       REF=tests/fixtures/comprehensive_ref.fa \
       REGION=chr1:1000-1300
   ```

   Output: `.render-dev/expanded.png`, `.render-dev/dv-strips.png`, `.render-dev/dv-composite.png`.

   For a single mode (faster iteration):

   ```bash
   make render-viewer MODE=dv-strips REGION=chr1:1000-1300
   ```

3. **Read each PNG with the `Read` tool.** Claude Code displays PNGs visually. Confirm the change does what you intended. Iterate until correct.

4. **Optional: render a specific variant active.** Add `--active-variant chrom:pos:ref>alt` when invoking the script directly to drive the "supports-variant" channel in DV modes:

   ```bash
   python scripts/render_viewer.py \
       --bam tests/fixtures/comprehensive.bam \
       --reference tests/fixtures/comprehensive_ref.fa \
       --region chr1:1000-1300 \
       --mode dv-strips \
       --active-variant chr1:1050:A>T \
       --out .render-dev/dv-strips-with-variant.png
   ```

## Tips

- The fixture BAMs at `tests/fixtures/comprehensive.bam` (with `comprehensive_ref.fa`) have real variants in `chr1:1000-2000` — good default region.
- Renders take ~1 second per call (per-case fresh Chromium); rendering all three modes runs serially and takes ~3-4 seconds total.
- Output files are gitignored (`.render-dev/` is in `.gitignore`) so you can leave them around between runs.
- If `playwright install chromium` hasn't been run yet, `make eval-vision-setup` does it once.

## Underlying tools

- `scripts/render_viewer.py` — full CLI with `--bam`, `--ui-init-json`, `--mode`, `--active-variant`, `--width`, `--out`, `--print-path-only`.
- `src/bamcp/eval/renderer.py` — the `ViewerRenderer` class that the CLI wraps. Same class powers the visual eval harness, so what you see in `render-viewer` is exactly what the eval LLM sees.
