# BAMCP as an MCP App: Strategic Architecture Plan

## The Reframe

The previous review identified "BAMCP is an MCP server and your vision is a web app" as the core architectural gap. That was wrong — or rather, it was the wrong vision. MCP Apps changes the equation entirely.

The insight: **the LLM IS the genetic curation expert.** BAMCP doesn't need to be a standalone web app where a layman stares at read alignments and somehow figures out what's going on. It needs to be a set of tools and interactive visualizations that plug into an AI assistant, where the AI bridges the gap between raw genomic data and human understanding.

This is a fundamentally better product for a layman audience than IGV could ever be. IGV assumes you know what you're looking at. BAMCP + an LLM assumes you don't.

---

## What MCP Apps Gives You

The MCP Apps extension (spec finalized 2026-01-26, supported in Claude, Claude Desktop, VS Code, ChatGPT, Goose) lets tools return interactive HTML UIs that render inline in the conversation. The key primitives:

**`_meta.ui.resourceUri`** on tool definitions — tells the host "this tool has an interactive UI, fetch it from this `ui://` resource and render it in a sandboxed iframe."

**`app.callServerTool()`** — the UI can call back to the MCP server for more data. This is how the viewer fetches reads as the user pans/zooms, without going through the LLM. Round-trip to server, not model inference.

**`app.updateModelContext()`** — the UI silently tells the LLM what the user is looking at. "User is viewing chr17:7,674,220-7,674,280. 3 reads show C>T at position 7,674,256. Coverage: 45x." The LLM can then proactively look up ClinVar, explain significance, flag concerns.

**`app.sendMessage()`** — the UI can inject a message into the conversation. User clicks a variant → viewer sends "Tell me about this variant at chr17:7,674,256 C>T" → LLM responds with ClinVar lookup and plain-English explanation.

**`visibility: ["app"]`** — tools hidden from the LLM, only callable by the UI. Perfect for high-frequency data fetching (browse a region, get coverage) that the model shouldn't see or decide about.

**Chunked data loading** — documented pattern for streaming large binary data in 500KB chunks via app-only tools. Directly applicable to BAM data.

**Fullscreen mode** — `app.requestDisplayMode("fullscreen")` for when the alignment viewer needs real estate.

---

## Why This Architecture is Right for the Use Case

**The "layman understanding" problem is solved by the LLM, not the UI.** A standalone genomics viewer needs tooltips, legends, help pages, and educational content baked into the interface. An MCP App has an AI sitting right next to the visualization explaining everything in natural language. User: "What's that red thing?" → AI: "That's a deletion — 3 bases are missing from this read compared to the reference genome. This is in the TP53 gene, which is commonly mutated in cancer. Let me check if this is a known variant..." This is not possible in a standalone web app.

**The curation workflow is conversational, not form-based.** Clinical variant curation is a decision process: review evidence, consult databases, weigh pathogenicity criteria. This maps naturally to a conversation with an expert, not a form with dropdowns. The MCP App provides visual evidence (alignment viewer, coverage plots, variant tables) while the LLM provides reasoning.

**Context preservation.** The viewer, the AI's explanation, and the user's questions all live in the same conversation thread. "Go back to that variant we looked at 10 minutes ago" is a natural language request, not a bookmark system.

**Distribution.** Installing an MCP server is `npx` or `uv run`. No hosting infrastructure, no auth system, no cloud costs. A bioinformatician shares a connector URL or a Claude Desktop config JSON, and their collaborator has a genomics curation tool running.

---

## Architecture

### Overview

```
┌─────────────────────────────────────┐
│  MCP Host (Claude / ChatGPT / etc)  │
│                                     │
│  ┌───────────────────────────────┐  │
│  │     LLM (curation expert)     │  │
│  │  - explains variants          │  │
│  │  - calls ClinVar/gnomAD       │  │
│  │  - reasons about pathogenicity│  │
│  └──────────┬────────────────────┘  │
│             │ tool calls             │
│  ┌──────────▼────────────────────┐  │
│  │  Alignment Viewer (MCP App)   │  │
│  │  ┌────────────────────────┐   │  │
│  │  │ Canvas: reads, coverage│   │  │
│  │  │ Variant table          │   │  │
│  │  │ Gene track             │   │  │
│  │  └────────┬───────────────┘   │  │
│  │           │ callServerTool()  │  │
│  │           │ updateModelContext│  │
│  └───────────┼───────────────────┘  │
└──────────────┼──────────────────────┘
               │
    ┌──────────▼──────────────┐
    │   BAMCP MCP Server      │
    │                         │
    │  Python + pysam         │
    │  - parse BAM/CRAM       │
    │  - detect variants      │
    │  - compute coverage     │
    │  - ClinVar API client   │
    │  - gnomAD API client    │
    │                         │
    │  TypeScript (optional)  │
    │  - ext-apps SDK server  │
    │  - serve bundled HTML   │
    └─────────────────────────┘
```

### Language Decision

**Stay Python for the MCP server.** pysam is the critical dependency and only works in CPython. The Python MCP SDK (FastMCP) supports MCP Apps — the `say-server` and `qr-server` examples in the ext-apps repo are both Python with embedded HTML UIs. The pattern is:

- Python server registers tools with `_meta.ui` fields pointing to `ui://` resources
- Python server registers resource handlers that serve bundled HTML
- The HTML UI (which runs in the browser iframe) uses the `@modelcontextprotocol/ext-apps` JS SDK via a `<script>` tag or Vite bundle

The UI is always JS/HTML regardless of server language. Vite + `vite-plugin-singlefile` bundles it into a single HTML string the Python server returns as a resource.

### Tool Design

Two categories: **LLM-facing tools** the model decides when to call, and **App-only tools** (`visibility: ["app"]`) the viewer calls directly for data.

#### LLM-Facing Tools (the model sees and calls these)

| Tool | Purpose | Returns |
|------|---------|---------|
| `load_bam` | Open a BAM/CRAM from URL or path. Validates index, lists contigs, reports basic stats. | Text summary: reference build, contig count, mapped reads. Structured: contig list with lengths. |
| `visualize_region` | Render the alignment viewer for a genomic region. **This is the MCP App tool.** | Text for model ("Showing chr17:7,674,220-7,674,280, 45x coverage, 3 variants detected"). UI renders the full interactive viewer. |
| `lookup_clinvar` | Query ClinVar for a variant (chrom, pos, ref, alt). | Text: clinical significance, review status, conditions, evidence summary. Plain English. |
| `lookup_gnomad` | Query gnomAD for population frequency. | Text: allele frequency across populations, homozygote counts. |
| `lookup_dbsnp` | Get rsID and basic info for a position. | Text: rsID, variant type, MAF. |
| `get_region_summary` | Text summary of a region for the LLM to reason about. Variants, coverage, notable features. | Text-only. No UI. For LLM reasoning. |
| `list_contigs` | List available reference sequences in the loaded BAM. | Text list. |

#### App-Only Tools (viewer calls these, hidden from model)

| Tool | Purpose | Returns |
|------|---------|---------|
| `fetch_reads` | Get alignment data for a viewport region. Chunked. | Structured: read objects with CIGAR, sequence, quality, position, flags. |
| `fetch_coverage` | Get per-base coverage for a region. | Structured: array of coverage values. |
| `fetch_variants` | Get detected variants in a region. | Structured: variant objects with position, ref, alt, depth, VAF. |
| `fetch_reference` | Get reference sequence bases for a region. | Structured: sequence string. |
| `fetch_genes` | Get gene annotations overlapping a region. | Structured: gene/transcript objects with exon boundaries. |

The key insight: app-only tools handle the high-frequency data fetching (every pan/zoom in the viewer) without polluting the LLM's context window or requiring inference. The LLM only sees `visualize_region` (once) and the context updates the viewer pushes.

### The Viewer (MCP App UI)

A single-file HTML bundle (Vite + vite-plugin-singlefile) containing:

- **Canvas alignment renderer** — your existing renderer, adapted. Reads, CIGAR visualization, mismatch highlighting, soft clips, strand coloring. All already built.
- **Coverage track** — bar chart above the alignment view.
- **Reference sequence track** — A/C/G/T colored bases at high zoom.
- **Gene track** — gene models with exons/introns (RefSeq).
- **Variant table** — clickable list of detected variants with basic stats.
- **Navigation controls** — coordinate input, zoom slider, contig selector.
- **Fullscreen toggle** — `app.requestDisplayMode("fullscreen")`.

Communication pattern:
```
User pans viewer right
  → viewer calls fetch_reads({ region: "chr17:7674300-7674400" })  [app-only tool]
  → server returns read data
  → viewer renders new reads on Canvas
  → viewer calls updateModelContext({
      content: [{ type: "text", text: `
        ---
        region: chr17:7,674,300-7,674,400
        coverage: 42x mean
        variants:
          - pos: 7,674,356, C>T, VAF: 0.47, depth: 42
        gene: TP53, exon 7
        ---
        User is viewing alignments in TP53 exon 7.
        One heterozygous variant detected: C>T at 7,674,356.
      `}]
    })

LLM sees context update, proactively responds:
  "I can see you're looking at a C>T change in TP53 exon 7.
   Let me check ClinVar for this..."
  → calls lookup_clinvar({ chrom: "chr17", pos: 7674356, ref: "C", alt: "T" })
  → "This variant (rs28934578) is classified as Pathogenic in ClinVar
     with a 4-star review status..."
```

This loop — user explores visually, AI explains automatically — is the entire product.

### Model Context Updates: The Secret Weapon

`updateModelContext` is what makes this fundamentally different from a standalone viewer. Every significant user interaction in the viewer generates a context update:

- **Navigation**: "User is now viewing [region] in [gene]. Coverage: [Nx]."
- **Variant hover/click**: "User is inspecting variant at [pos]: [ref]>[alt], VAF [x], depth [y]."
- **Zoom level change**: "User zoomed to base-level view" → triggers reference track rendering and more detailed variant descriptions.
- **Selection**: "User selected 5 reads spanning [region], all showing [variant]."

The LLM uses these to drive the conversation. It can:
- Proactively fetch ClinVar/gnomAD annotations
- Warn about low-quality evidence ("I notice this variant has strand bias — 38 of 40 supporting reads are on the forward strand")
- Suggest next regions to examine ("TP53 has several other commonly mutated positions. Want me to check exons 5-8?")
- Explain what the visual patterns mean ("The coverage drop here suggests a possible deletion")

---

## Phased Roadmap

### Phase 1: Core Viewer as MCP App (Weeks 1-3)

**Goal**: User provides BAM URL to their LLM, gets an interactive alignment viewer in the conversation, AI explains what they're seeing.

- [ ] Restructure server to register tools with `_meta.ui` for MCP Apps
- [ ] Create Vite build pipeline for viewer HTML → single-file bundle
- [ ] Adapt existing Canvas renderer to use `App` SDK (`callServerTool` instead of raw `postMessage`)
- [ ] Implement `load_bam` (LLM-facing) and `fetch_reads` (app-only) tools
- [ ] Implement `visualize_region` as the MCP App tool
- [ ] Add `updateModelContext` calls for navigation and viewport changes
- [ ] Basic fullscreen support
- [ ] Test with Claude via cloudflared tunnel or custom connector

**What transfers from current codebase**: parsers.py (pysam wrappers), Canvas renderer, CIGAR parsing, read packing algorithm, mismatch highlighting, fixture generation. All of this works as-is.

**What changes**: tool registration (add `_meta.ui`), viewer communication layer (App SDK instead of raw postMessage), resource serving pattern.

### Phase 2: Reference + Gene Context (Weeks 4-5)

**Goal**: Viewer shows reference sequence at high zoom and gene annotations so variants have biological context.

- [ ] `fetch_reference` tool — serve reference bases from a FASTA (bundled or URL-based for hg19/hg38)
- [ ] Reference track in Canvas renderer — A/C/G/T colored bases
- [ ] `fetch_genes` tool — serve RefSeq gene annotations from a pre-processed BED/SQLite
- [ ] Gene track in Canvas renderer — gene models with exon/intron structure
- [ ] Richer `updateModelContext` — include gene name, exon number, functional context

### Phase 3: ClinVar Integration (Weeks 6-7)

**Goal**: LLM can look up any variant in ClinVar and explain its clinical significance.

- [ ] `lookup_clinvar` tool — hit NCBI E-utilities or a local ClinVar VCF
- [ ] `lookup_gnomad` tool — gnomAD API or local VCF
- [ ] `lookup_dbsnp` tool — rsID lookup
- [ ] Variant annotation overlay in viewer — color-code variants by ClinVar significance
- [ ] `sendMessage` from viewer on variant click → triggers LLM to explain
- [ ] **Miscalling safeguards** (see below)

### Phase 4: Curation Workflow (Weeks 8-10)

**Goal**: A complete variant curation session where a layman can work through a list of variants with AI guidance.

- [ ] VCF overlay — load a called VCF alongside the BAM
- [ ] Variant curation queue — viewer presents variants one-by-one, AI explains each
- [ ] ACMG criteria guidance — LLM walks through evidence categories
- [ ] Export — generate a curation report (markdown or PDF via the LLM)
- [ ] Multi-sample comparison (if applicable)

---

## Miscalling Safeguards (Non-Negotiable)

This is even more important in the MCP App model because the LLM will be interpreting variants for non-experts. The AI's explanation carries authority. Specific requirements:

**In the viewer UI:**
- Every variant displayed shows: "Observed difference from reference — not a validated variant call"
- Color-coded confidence: green (high depth + VAF + balanced strands), yellow (moderate), red (low quality / possible artifact)
- Persistent footer: "BAMCP is a visualization tool. Not validated for clinical use."

**In the LLM's tool responses:**
- `lookup_clinvar` returns ClinVar review stars prominently. A 1-star review is very different from 4-star.
- `get_region_summary` includes quality metrics alongside each variant
- All tool descriptions include: "This tool provides research-grade information. Clinical decisions require validated testing and qualified interpretation."

**In model context updates:**
- Always include quality metrics: "VAF: 0.47, depth: 42, strand ratio: 22F/20R" — so the LLM can reason about data quality
- Flag common artifacts: "Note: this region has high homopolymer content" or "Variant is adjacent to read end, possible alignment artifact"

**Structural safeguard:** The LLM is actually better at communicating uncertainty than a static UI. It can say "This looks like it could be pathogenic based on ClinVar, but I'd note the coverage here is only 12x which is below typical clinical thresholds. A clinical-grade test would provide more confidence."

---

## What Transfers vs. What Changes

### Keeps Working As-Is
- `parsers.py` — all pysam wrappers, read fetching, variant detection, coverage computation
- Canvas renderer logic — read packing, CIGAR rendering, mismatch highlighting, color schemes
- Fixture generation — test BAM files with known variants
- Unit tests for parsers
- Docker infrastructure (may need adjustment for the Vite build step)

### Needs Adaptation
- **Tool registration** — add `_meta.ui.resourceUri` to `visualize_region`, add `visibility: ["app"]` to data-fetching tools
- **Viewer communication** — replace raw `postMessage` handler with `App` SDK. `ontoolresult` replaces the current `message` event listener. `callServerTool` replaces `postMessage` requests.
- **Resource serving** — register a `ui://bamcp/viewer.html` resource that returns the Vite-bundled HTML

### New
- Vite build pipeline for the viewer HTML
- `updateModelContext` calls from the viewer
- ClinVar/gnomAD/dbSNP API clients
- Gene annotation data layer (RefSeq BED/SQLite)
- Reference sequence serving (FASTA)

---

## Technical Decisions

**TypeScript wrapper vs pure Python server?**

Pure Python. The `say-server` example proves Python MCP servers work with MCP Apps. Use FastMCP for tool/resource registration. The HTML bundle is just a string the Python server returns — it doesn't matter what language serves it. Adding a TypeScript wrapper just to use the `registerAppTool` helper isn't worth the complexity when you can set `_meta.ui` directly in Python.

**igv.js vs custom Canvas renderer?**

Custom. Now that you're building an MCP App (not a standalone viewer), igv.js's strengths (standalone operation, file loading UI, session management) don't apply — those are handled by the MCP host and the LLM. Your renderer can be optimized for the MCP App context: simpler chrome, emphasis on clarity over expert features, tight integration with `callServerTool` and `updateModelContext`. igv.js in an iframe-in-an-iframe is also a fragile stack.

**ClinVar: API vs local VCF?**

Start with NCBI E-utilities API (no data to ship, always current). Add local VCF option later for offline use or speed. The API is free, rate-limited to 10 req/sec with an API key (3/sec without). More than enough for interactive curation.

**Reference genomes: bundle vs download?**

Don't bundle. For hg38 reference sequence: use the UCSC DAS server or a pre-built 2bit file the user provides a path to. For gene annotations: ship a compressed RefSeq BED (~5MB) for hg19 and hg38. This is small enough to bundle with the package.

---

## Distribution Model

For Claude Desktop (stdio):
```json
{
  "mcpServers": {
    "bamcp": {
      "command": "uv",
      "args": ["run", "bamcp", "--stdio"]
    }
  }
}
```

For Claude.ai / ChatGPT (remote, via cloudflared or hosted):
```
Custom connector → https://your-tunnel.trycloudflare.com/mcp
```

For npx distribution (if you add a Node wrapper):
```json
{
  "mcpServers": {
    "bamcp": {
      "command": "npx",
      "args": ["-y", "@bamcp/server", "--stdio"]
    }
  }
}
```

The `uv run` path is cleanest for a Python package. Publish to PyPI, users install with `uv` or `pip`, run with `uv run bamcp --stdio`. Zero config.

---

## What Makes This Special

There is no MCP App for genomics today. There are genomics MCP servers (tools that return text about variants) and there are standalone viewers (IGV, JBrowse). Nobody has combined interactive genomic visualization with LLM-mediated interpretation in a single conversational interface.

The closest analog isn't another genomics tool — it's the `map-server` or `cohort-heatmap-server` examples in ext-apps, but for genomics. Interactive visualization that the AI can reason about.

The moat is the combination:
1. Real genomic data parsing (pysam — hard to replicate in JS)
2. Interactive visualization in the conversation (MCP App)
3. AI-mediated interpretation (the LLM explains everything)
4. Database integration (ClinVar, gnomAD — the AI knows how to use them)
5. Your bioinformatics domain knowledge embedded in the tool design

This is a legitimate product. A genetics counselor could use this to review results with patients. A researcher could use it to triage variants. A curious person who got their WES results could actually understand them — because they're not alone staring at aligned reads; they have an AI explaining every detail.