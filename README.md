# BAMCP

**Interactive BAM/CRAM variant visualization for AI assistants via the Model Context Protocol**

[![MCP](https://img.shields.io/badge/MCP-Apps%20Extension-blue)](https://modelcontextprotocol.io)
[![Python](https://img.shields.io/badge/python-3.10+-green)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

BAMCP brings IGV-style alignment visualization directly into your AI conversations. Browse BAM/CRAM files, inspect variant evidence, and navigate genomic regions—all through natural language interaction with Claude or other MCP-enabled assistants.

<p align="center">
  <img src="docs/assets/demo.png" alt="BAMCP Demo" width="800">
</p>

---

## Table of Contents

- [Why BAMCP?](#why-bamcp)
- [Features](#features)
- [Architecture](#architecture)
- [Installation](#installation)
- [Usage](#usage)
- [Tools Reference](#tools-reference)
- [Configuration](#configuration)
- [Implementation Details](#implementation-details)
- [Development](#development)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)

---

## Why BAMCP?

Existing genomics MCP servers (bio-mcp-samtools, AWS HealthOmics MCP) provide command-line operations but return text/JSON. **BAMCP is the first to leverage the [MCP Apps Extension](https://github.com/modelcontextprotocol/ext-apps) to render interactive visualizations inline.**

| Feature | CLI-based MCPs | BAMCP |
|---------|---------------|-------|
| Read BAM/CRAM files | ✅ | ✅ |
| Coverage statistics | ✅ | ✅ |
| **Interactive visualization** | ❌ | ✅ |
| **Variant browsing UI** | ❌ | ✅ |
| **Pan/zoom navigation** | ❌ | ✅ |
| **Inline in chat** | ❌ | ✅ |

---

## Features

- 🔬 **Alignment Viewer** — Visualize reads with color-coded mismatches, insertions, deletions, and soft clips
- 📊 **Coverage Track** — Real-time depth of coverage across the viewing window
- 🎯 **Variant Highlighting** — Automatic detection and highlighting of positions with non-reference alleles
- 🧭 **Region Navigation** — Jump to coordinates, genes, or specific variants
- 📁 **Format Support** — BAM, CRAM (with reference), and indexed remote files (HTTP/S3)
- ⚡ **Canvas Rendering** — High-performance visualization of thousands of reads via HTML5 Canvas/WebGL

---

## Architecture

### High-Level Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        MCP Client (Claude, Cursor, etc.)                    │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  │ MCP Protocol (JSON-RPC over stdio/HTTP)
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          BAMCP MCP Server (Python)                          │
│                                                                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────────┐  │
│  │                 │  │                 │  │                             │  │
│  │  pysam Layer    │  │  Tool Handlers  │  │  UI Resource Provider       │  │
│  │                 │  │                 │  │                             │  │
│  │  • AlignmentFile│  │  • browse_region│  │  • ui://bamcp/viewer        │  │
│  │  • fetch()      │  │  • get_variants │  │  • text/html+mcp            │  │
│  │  • pileup()     │  │  • get_coverage │  │  • Sandboxed iframe content │  │
│  │  • CRAM support │  │  • list_contigs │  │                             │  │
│  │                 │  │  • jump_to      │  │                             │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────────────────┘  │
│                                                                             │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  │ postMessage (JSON-RPC)
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      MCP App UI (Sandboxed iframe)                          │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                     Canvas Rendering Layer                            │  │
│  │                                                                       │  │
│  │  • Read rectangles with CIGAR-based rendering                         │  │
│  │  • Mismatch highlighting (A=green, T=red, G=orange, C=blue)          │  │
│  │  • Insertion/deletion markers                                         │  │
│  │  • Coverage histogram                                                 │  │
│  │  • Soft-clip visualization                                            │  │
│  │                                                                       │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                     Control Layer (Vanilla JS)                        │  │
│  │                                                                       │  │
│  │  • Navigation (pan, zoom, jump-to)                                    │  │
│  │  • Variant table                                                      │  │
│  │  • Read info tooltip                                                  │  │
│  │  • MCP postMessage bridge                                             │  │
│  │                                                                       │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
User: "Show me reads at chr17:7577000-7577500 in tumor.bam"
                              │
                              ▼
                    ┌─────────────────┐
                    │  MCP Client     │
                    │  (Claude)       │
                    └────────┬────────┘
                             │ tools/call: browse_region
                             ▼
                    ┌─────────────────┐
                    │  BAMCP Server   │
                    │                 │
                    │  1. Parse region│
                    │  2. pysam.fetch │
                    │  3. Serialize   │
                    │  4. Return UI   │
                    └────────┬────────┘
                             │ CallToolResult with ui/resourceUri
                             ▼
                    ┌─────────────────┐
                    │  MCP Client     │
                    │  renders iframe │
                    └────────┬────────┘
                             │ postMessage: init with read data
                             ▼
                    ┌─────────────────┐
                    │  Canvas UI      │
                    │                 │
                    │  1. Parse reads │
                    │  2. Pack rows   │
                    │  3. Render      │
                    └─────────────────┘
```

### Component Details

#### 1. MCP Server (Python)

The server implements the MCP protocol and exposes tools for genomic data access.

```python
# src/bamcp/server.py

from mcp.server import Server
from mcp.types import Tool, TextContent, ImageContent

app = Server("bamcp")

@app.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="browse_region",
            description="View aligned reads in a genomic region",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {"type": "string", "description": "Path to BAM/CRAM file"},
                    "region": {"type": "string", "description": "Region (e.g., chr1:1000-2000)"},
                    "reference": {"type": "string", "description": "Reference FASTA (required for CRAM)"}
                },
                "required": ["file_path", "region"]
            }
        ),
        # ... other tools
    ]
```

#### 2. pysam Integration

```python
# src/bamcp/parsers.py

import pysam
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class AlignedRead:
    name: str
    sequence: str
    qualities: List[int]
    cigar: str
    position: int
    end_position: int
    mapping_quality: int
    is_reverse: bool
    mismatches: List[dict]  # [{pos: int, ref: str, alt: str}]

@dataclass
class RegionData:
    contig: str
    start: int
    end: int
    reads: List[AlignedRead]
    coverage: List[int]
    variants: List[dict]
    reference_sequence: Optional[str]

def fetch_region(
    bam_path: str,
    region: str,
    reference_path: Optional[str] = None,
    max_reads: int = 10000,
    min_mapq: int = 0
) -> RegionData:
    """
    Fetch reads from a BAM/CRAM file for a given region.
    
    Args:
        bam_path: Path to BAM/CRAM file (local or remote)
        region: Genomic region (e.g., "chr1:1000-2000")
        reference_path: Path to reference FASTA (required for CRAM)
        max_reads: Maximum reads to return (downsampling if exceeded)
        min_mapq: Minimum mapping quality filter
    
    Returns:
        RegionData with reads, coverage, and detected variants
    """
    # Parse region
    contig, coords = region.replace(",", "").split(":")
    start, end = map(int, coords.split("-"))
    
    # Open alignment file
    mode = "rc" if bam_path.endswith(".cram") else "rb"
    samfile = pysam.AlignmentFile(
        bam_path,
        mode,
        reference_filename=reference_path
    )
    
    reads = []
    coverage = [0] * (end - start)
    base_counts = [{} for _ in range(end - start)]  # For variant detection
    
    # Fetch reference sequence if available
    ref_seq = None
    if reference_path:
        with pysam.FastaFile(reference_path) as fasta:
            ref_seq = fasta.fetch(contig, start, end)
    
    # Iterate over reads
    read_count = 0
    for read in samfile.fetch(contig, start, end):
        if read.mapping_quality < min_mapq:
            continue
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
            
        read_count += 1
        if read_count > max_reads:
            # TODO: Implement smart downsampling
            break
        
        # Extract mismatches
        mismatches = []
        if ref_seq and read.query_sequence:
            aligned_pairs = read.get_aligned_pairs(with_seq=True)
            for qpos, rpos, ref_base in aligned_pairs:
                if rpos is None or qpos is None:
                    continue
                if start <= rpos < end:
                    query_base = read.query_sequence[qpos]
                    if ref_base and query_base != ref_base.upper():
                        mismatches.append({
                            "pos": rpos,
                            "ref": ref_base.upper(),
                            "alt": query_base
                        })
                    # Update base counts for variant detection
                    idx = rpos - start
                    base_counts[idx][query_base] = base_counts[idx].get(query_base, 0) + 1
        
        # Update coverage
        for pos in range(max(read.reference_start, start), min(read.reference_end or read.reference_start, end)):
            coverage[pos - start] += 1
        
        reads.append(AlignedRead(
            name=read.query_name,
            sequence=read.query_sequence or "",
            qualities=list(read.query_qualities or []),
            cigar=read.cigarstring or "",
            position=read.reference_start,
            end_position=read.reference_end or read.reference_start,
            mapping_quality=read.mapping_quality,
            is_reverse=read.is_reverse,
            mismatches=mismatches
        ))
    
    samfile.close()
    
    # Detect variants from pileup
    variants = detect_variants(base_counts, ref_seq, contig, start)
    
    return RegionData(
        contig=contig,
        start=start,
        end=end,
        reads=reads,
        coverage=coverage,
        variants=variants,
        reference_sequence=ref_seq
    )

def detect_variants(
    base_counts: List[dict],
    ref_seq: Optional[str],
    contig: str,
    start: int,
    min_vaf: float = 0.1,
    min_depth: int = 10
) -> List[dict]:
    """Detect variants from base counts."""
    variants = []
    
    if not ref_seq:
        return variants
    
    for i, counts in enumerate(base_counts):
        if not counts:
            continue
            
        total = sum(counts.values())
        if total < min_depth:
            continue
            
        ref_base = ref_seq[i].upper()
        
        for base, count in counts.items():
            if base == ref_base:
                continue
            vaf = count / total
            if vaf >= min_vaf:
                variants.append({
                    "contig": contig,
                    "position": start + i,
                    "ref": ref_base,
                    "alt": base,
                    "vaf": round(vaf, 3),
                    "depth": total,
                    "alt_count": count
                })
    
    return variants
```

#### 3. MCP Apps UI Resource

```python
# src/bamcp/resources.py

from mcp.server import Server
from mcp.types import Resource, ResourceTemplate

def register_ui_resources(app: Server):
    """Register UI resources for the MCP Apps extension."""
    
    @app.list_resources()
    async def list_resources() -> list[Resource]:
        return [
            Resource(
                uri="ui://bamcp/viewer",
                name="BAMCP Alignment Viewer",
                mimeType="text/html+mcp",
                description="Interactive BAM/CRAM alignment visualization"
            )
        ]
    
    @app.read_resource()
    async def read_resource(uri: str) -> str:
        if uri == "ui://bamcp/viewer":
            return get_viewer_html()
        raise ValueError(f"Unknown resource: {uri}")

def get_viewer_html() -> str:
    """Return the HTML content for the alignment viewer."""
    return """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BAMCP Viewer</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
        
        #container { display: flex; flex-direction: column; height: 100vh; }
        
        #toolbar {
            display: flex;
            gap: 8px;
            padding: 8px;
            background: #f5f5f5;
            border-bottom: 1px solid #ddd;
        }
        
        #toolbar input {
            flex: 1;
            padding: 6px 10px;
            border: 1px solid #ccc;
            border-radius: 4px;
        }
        
        #toolbar button {
            padding: 6px 12px;
            background: #4a90d9;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }
        
        #toolbar button:hover { background: #357abd; }
        
        #viewer {
            flex: 1;
            position: relative;
            overflow: hidden;
        }
        
        #coverage-canvas, #reads-canvas {
            position: absolute;
            left: 0;
        }
        
        #coverage-canvas { top: 0; height: 60px; }
        #reads-canvas { top: 60px; }
        
        #tooltip {
            position: absolute;
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 8px;
            border-radius: 4px;
            font-size: 12px;
            pointer-events: none;
            display: none;
            max-width: 300px;
        }
        
        #variant-panel {
            height: 150px;
            border-top: 1px solid #ddd;
            overflow-y: auto;
            font-size: 12px;
        }
        
        #variant-panel table {
            width: 100%;
            border-collapse: collapse;
        }
        
        #variant-panel th, #variant-panel td {
            padding: 4px 8px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }
        
        #variant-panel tr:hover { background: #f0f7ff; cursor: pointer; }
        #variant-panel th { background: #f5f5f5; position: sticky; top: 0; }
    </style>
</head>
<body>
    <div id="container">
        <div id="toolbar">
            <input type="text" id="region-input" placeholder="chr1:1000-2000">
            <button id="go-btn">Go</button>
            <button id="zoom-in">+</button>
            <button id="zoom-out">−</button>
        </div>
        <div id="viewer">
            <canvas id="coverage-canvas"></canvas>
            <canvas id="reads-canvas"></canvas>
            <div id="tooltip"></div>
        </div>
        <div id="variant-panel">
            <table>
                <thead>
                    <tr>
                        <th>Position</th>
                        <th>Ref</th>
                        <th>Alt</th>
                        <th>VAF</th>
                        <th>Depth</th>
                    </tr>
                </thead>
                <tbody id="variant-table"></tbody>
            </table>
        </div>
    </div>
    
    <script>
    // ============================================
    // BAMCP Viewer - Canvas-based alignment renderer
    // ============================================
    
    const BASE_COLORS = {
        'A': '#22c55e', // green
        'T': '#ef4444', // red
        'G': '#f97316', // orange
        'C': '#3b82f6', // blue
        'N': '#9ca3af'  // gray
    };
    
    const READ_HEIGHT = 12;
    const READ_GAP = 2;
    const COVERAGE_HEIGHT = 60;
    
    class BAMCPViewer {
        constructor() {
            this.coverageCanvas = document.getElementById('coverage-canvas');
            this.readsCanvas = document.getElementById('reads-canvas');
            this.coverageCtx = this.coverageCanvas.getContext('2d');
            this.readsCtx = this.readsCanvas.getContext('2d');
            this.tooltip = document.getElementById('tooltip');
            this.variantTable = document.getElementById('variant-table');
            
            this.data = null;
            this.viewport = { start: 0, end: 1000 };
            this.packedRows = [];
            
            this.setupEventListeners();
            this.setupMCPBridge();
            this.resize();
            
            window.addEventListener('resize', () => this.resize());
        }
        
        setupEventListeners() {
            // Navigation
            document.getElementById('go-btn').addEventListener('click', () => {
                const region = document.getElementById('region-input').value;
                this.requestRegion(region);
            });
            
            document.getElementById('zoom-in').addEventListener('click', () => this.zoom(0.5));
            document.getElementById('zoom-out').addEventListener('click', () => this.zoom(2));
            
            // Pan with mouse drag
            let isDragging = false;
            let lastX = 0;
            
            this.readsCanvas.addEventListener('mousedown', (e) => {
                isDragging = true;
                lastX = e.clientX;
            });
            
            window.addEventListener('mousemove', (e) => {
                if (!isDragging) return;
                const dx = e.clientX - lastX;
                lastX = e.clientX;
                this.pan(-dx);
            });
            
            window.addEventListener('mouseup', () => isDragging = false);
            
            // Tooltip on hover
            this.readsCanvas.addEventListener('mousemove', (e) => {
                if (isDragging) return;
                this.showTooltip(e);
            });
            
            this.readsCanvas.addEventListener('mouseout', () => {
                this.tooltip.style.display = 'none';
            });
            
            // Scroll to zoom
            this.readsCanvas.addEventListener('wheel', (e) => {
                e.preventDefault();
                const factor = e.deltaY > 0 ? 1.2 : 0.8;
                this.zoom(factor);
            });
        }
        
        setupMCPBridge() {
            // MCP Apps communicate via postMessage using JSON-RPC
            window.addEventListener('message', (event) => {
                const message = event.data;
                
                if (message.method === 'bamcp/init') {
                    this.loadData(message.params);
                } else if (message.method === 'bamcp/update') {
                    this.loadData(message.params);
                }
            });
            
            // Signal ready to host
            window.parent.postMessage({
                jsonrpc: '2.0',
                method: 'ready'
            }, '*');
        }
        
        requestRegion(region) {
            // Request new region from MCP server via host
            window.parent.postMessage({
                jsonrpc: '2.0',
                method: 'tools/call',
                params: {
                    name: 'browse_region',
                    arguments: { region }
                }
            }, '*');
        }
        
        loadData(data) {
            this.data = data;
            this.viewport = { start: data.start, end: data.end };
            
            document.getElementById('region-input').value = 
                `${data.contig}:${data.start}-${data.end}`;
            
            this.packReads();
            this.renderVariantTable();
            this.render();
        }
        
        packReads() {
            // Pack reads into rows to avoid overlap
            if (!this.data) return;
            
            this.packedRows = [];
            const rowEnds = []; // Track end position of each row
            
            // Sort reads by start position
            const sortedReads = [...this.data.reads].sort((a, b) => a.position - b.position);
            
            for (const read of sortedReads) {
                // Find first row where read fits
                let rowIndex = rowEnds.findIndex(end => end < read.position);
                
                if (rowIndex === -1) {
                    // Need new row
                    rowIndex = rowEnds.length;
                    rowEnds.push(read.end_position);
                    this.packedRows.push([]);
                } else {
                    rowEnds[rowIndex] = read.end_position;
                }
                
                this.packedRows[rowIndex].push(read);
            }
        }
        
        resize() {
            const container = document.getElementById('viewer');
            const width = container.clientWidth;
            const height = container.clientHeight - COVERAGE_HEIGHT;
            
            this.coverageCanvas.width = width;
            this.coverageCanvas.height = COVERAGE_HEIGHT;
            
            this.readsCanvas.width = width;
            this.readsCanvas.height = height;
            
            this.render();
        }
        
        render() {
            this.renderCoverage();
            this.renderReads();
        }
        
        renderCoverage() {
            const ctx = this.coverageCtx;
            const width = this.coverageCanvas.width;
            const height = this.coverageCanvas.height;
            
            ctx.clearRect(0, 0, width, height);
            
            if (!this.data) return;
            
            const coverage = this.data.coverage;
            const maxCov = Math.max(...coverage, 1);
            const scale = this.getScale();
            
            ctx.fillStyle = '#93c5fd';
            ctx.strokeStyle = '#3b82f6';
            ctx.beginPath();
            
            for (let i = 0; i < coverage.length; i++) {
                const x = (this.data.start + i - this.viewport.start) * scale;
                const h = (coverage[i] / maxCov) * (height - 20);
                
                if (i === 0) {
                    ctx.moveTo(x, height - h);
                } else {
                    ctx.lineTo(x, height - h);
                }
            }
            
            // Close path for fill
            ctx.lineTo((this.data.end - this.viewport.start) * scale, height);
            ctx.lineTo(0, height);
            ctx.closePath();
            ctx.fill();
            ctx.stroke();
            
            // Draw max coverage label
            ctx.fillStyle = '#333';
            ctx.font = '10px sans-serif';
            ctx.fillText(`Max: ${maxCov}x`, 5, 12);
        }
        
        renderReads() {
            const ctx = this.readsCtx;
            const width = this.readsCanvas.width;
            const height = this.readsCanvas.height;
            
            ctx.clearRect(0, 0, width, height);
            
            if (!this.data) return;
            
            const scale = this.getScale();
            
            for (let rowIdx = 0; rowIdx < this.packedRows.length; rowIdx++) {
                const y = rowIdx * (READ_HEIGHT + READ_GAP);
                
                if (y > height) break; // Off screen
                
                for (const read of this.packedRows[rowIdx]) {
                    this.renderRead(ctx, read, y, scale);
                }
            }
        }
        
        renderRead(ctx, read, y, scale) {
            const x = (read.position - this.viewport.start) * scale;
            const w = (read.end_position - read.position) * scale;
            
            // Skip if off screen
            if (x + w < 0 || x > this.readsCanvas.width) return;
            
            // Draw read body
            ctx.fillStyle = read.is_reverse ? '#a78bfa' : '#60a5fa';
            ctx.fillRect(x, y, w, READ_HEIGHT);
            
            // Draw mismatches if zoomed in enough
            if (scale > 2) {
                for (const mm of read.mismatches) {
                    const mx = (mm.pos - this.viewport.start) * scale;
                    ctx.fillStyle = BASE_COLORS[mm.alt] || '#9ca3af';
                    ctx.fillRect(mx, y, Math.max(scale, 2), READ_HEIGHT);
                }
            }
            
            // Draw direction indicator
            ctx.fillStyle = 'white';
            ctx.font = '8px sans-serif';
            const arrow = read.is_reverse ? '◀' : '▶';
            if (w > 15) {
                ctx.fillText(arrow, read.is_reverse ? x + 2 : x + w - 10, y + 9);
            }
        }
        
        getScale() {
            const width = this.readsCanvas.width;
            return width / (this.viewport.end - this.viewport.start);
        }
        
        zoom(factor) {
            const center = (this.viewport.start + this.viewport.end) / 2;
            const span = (this.viewport.end - this.viewport.start) * factor;
            
            this.viewport.start = Math.floor(center - span / 2);
            this.viewport.end = Math.ceil(center + span / 2);
            
            this.render();
        }
        
        pan(pixels) {
            const scale = this.getScale();
            const basePan = Math.round(pixels / scale);
            
            this.viewport.start += basePan;
            this.viewport.end += basePan;
            
            this.render();
        }
        
        showTooltip(e) {
            if (!this.data) return;
            
            const rect = this.readsCanvas.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;
            
            const scale = this.getScale();
            const genomePos = Math.floor(this.viewport.start + x / scale);
            const rowIdx = Math.floor(y / (READ_HEIGHT + READ_GAP));
            
            if (rowIdx >= this.packedRows.length) return;
            
            // Find read at position
            const read = this.packedRows[rowIdx].find(r => 
                r.position <= genomePos && r.end_position >= genomePos
            );
            
            if (read) {
                this.tooltip.innerHTML = `
                    <strong>${read.name}</strong><br>
                    Position: ${read.position}-${read.end_position}<br>
                    MAPQ: ${read.mapping_quality}<br>
                    CIGAR: ${read.cigar}<br>
                    Strand: ${read.is_reverse ? 'Reverse' : 'Forward'}
                `;
                this.tooltip.style.display = 'block';
                this.tooltip.style.left = (e.clientX + 10) + 'px';
                this.tooltip.style.top = (e.clientY + 10) + 'px';
            } else {
                this.tooltip.style.display = 'none';
            }
        }
        
        renderVariantTable() {
            if (!this.data) return;
            
            this.variantTable.innerHTML = this.data.variants.map(v => `
                <tr data-pos="${v.position}">
                    <td>${v.contig}:${v.position}</td>
                    <td>${v.ref}</td>
                    <td style="color: ${BASE_COLORS[v.alt]}">${v.alt}</td>
                    <td>${(v.vaf * 100).toFixed(1)}%</td>
                    <td>${v.depth}</td>
                </tr>
            `).join('');
            
            // Click to jump
            this.variantTable.querySelectorAll('tr').forEach(row => {
                row.addEventListener('click', () => {
                    const pos = parseInt(row.dataset.pos);
                    this.jumpTo(pos);
                });
            });
        }
        
        jumpTo(position) {
            const span = this.viewport.end - this.viewport.start;
            this.viewport.start = position - span / 2;
            this.viewport.end = position + span / 2;
            this.render();
        }
    }
    
    // Initialize viewer
    const viewer = new BAMCPViewer();
    </script>
</body>
</html>
"""
```

#### 4. Tool Handlers

```python
# src/bamcp/tools.py

import json
from typing import Any
from mcp.server import Server
from mcp.types import CallToolResult, TextContent

from .parsers import fetch_region, RegionData

def register_tools(app: Server):
    """Register all BAMCP tools."""
    
    @app.call_tool()
    async def call_tool(name: str, arguments: dict[str, Any]) -> CallToolResult:
        if name == "browse_region":
            return await handle_browse_region(arguments)
        elif name == "get_variants":
            return await handle_get_variants(arguments)
        elif name == "get_coverage":
            return await handle_get_coverage(arguments)
        elif name == "list_contigs":
            return await handle_list_contigs(arguments)
        else:
            raise ValueError(f"Unknown tool: {name}")

async def handle_browse_region(args: dict) -> CallToolResult:
    """Handle browse_region tool call."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference")
    
    data = fetch_region(file_path, region, reference)
    
    # Serialize for UI
    payload = {
        "contig": data.contig,
        "start": data.start,
        "end": data.end,
        "reads": [
            {
                "name": r.name,
                "sequence": r.sequence,
                "cigar": r.cigar,
                "position": r.position,
                "end_position": r.end_position,
                "mapping_quality": r.mapping_quality,
                "is_reverse": r.is_reverse,
                "mismatches": r.mismatches
            }
            for r in data.reads
        ],
        "coverage": data.coverage,
        "variants": data.variants,
        "reference_sequence": data.reference_sequence
    }
    
    return CallToolResult(
        content=[TextContent(type="text", text=json.dumps(payload))],
        _meta={
            "ui/resourceUri": "ui://bamcp/viewer",
            "ui/init": payload
        }
    )

async def handle_get_variants(args: dict) -> CallToolResult:
    """Return variants without UI."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference")
    min_vaf = args.get("min_vaf", 0.1)
    min_depth = args.get("min_depth", 10)
    
    data = fetch_region(file_path, region, reference)
    
    # Filter variants
    variants = [
        v for v in data.variants
        if v["vaf"] >= min_vaf and v["depth"] >= min_depth
    ]
    
    return CallToolResult(
        content=[TextContent(
            type="text",
            text=json.dumps({"variants": variants, "count": len(variants)})
        )]
    )

async def handle_get_coverage(args: dict) -> CallToolResult:
    """Return coverage statistics."""
    file_path = args["file_path"]
    region = args["region"]
    reference = args.get("reference")
    
    data = fetch_region(file_path, region, reference)
    
    coverage = data.coverage
    stats = {
        "region": f"{data.contig}:{data.start}-{data.end}",
        "mean": sum(coverage) / len(coverage) if coverage else 0,
        "min": min(coverage) if coverage else 0,
        "max": max(coverage) if coverage else 0,
        "median": sorted(coverage)[len(coverage)//2] if coverage else 0,
        "bases_covered": sum(1 for c in coverage if c > 0),
        "total_bases": len(coverage)
    }
    
    return CallToolResult(
        content=[TextContent(type="text", text=json.dumps(stats))]
    )

async def handle_list_contigs(args: dict) -> CallToolResult:
    """List contigs in a BAM/CRAM file."""
    import pysam
    
    file_path = args["file_path"]
    reference = args.get("reference")
    
    mode = "rc" if file_path.endswith(".cram") else "rb"
    samfile = pysam.AlignmentFile(file_path, mode, reference_filename=reference)
    
    contigs = [
        {"name": name, "length": length}
        for name, length in zip(samfile.references, samfile.lengths)
    ]
    
    samfile.close()
    
    return CallToolResult(
        content=[TextContent(type="text", text=json.dumps({"contigs": contigs}))]
    )
```

---

## Installation

### Prerequisites

- Python 3.10+
- [pysam](https://pysam.readthedocs.io/) (requires htslib)
- Node.js 18+ (for UI development only)

### From PyPI (coming soon)

```bash
pip install bamcp
```

### From Source

```bash
git clone https://github.com/yourusername/bamcp.git
cd bamcp

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows

# Install with dev dependencies
pip install -e ".[dev]"
```

### MCP Client Configuration

#### Claude Desktop

Add to `~/.config/claude/claude_desktop_config.json` (Linux/macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "bamcp": {
      "command": "python",
      "args": ["-m", "bamcp"],
      "env": {
        "BAMCP_REFERENCE": "/path/to/hg38.fa"
      }
    }
  }
}
```

#### Cursor / VS Code

Add to your MCP configuration:

```json
{
  "mcpServers": {
    "bamcp": {
      "command": "uvx",
      "args": ["bamcp"]
    }
  }
}
```

---

## Usage

### Browse a Region

> "Show me the reads at chr17:7577000-7577500 in /data/tumor.bam"

The assistant will call `browse_region` and render an interactive alignment viewer inline.

### Inspect a Variant

> "What's the evidence for the variant at chr7:140453136 in my sample?"

BAMCP will center on the position and highlight reads supporting reference vs. alternate alleles.

### Get Coverage Stats

> "What's the average coverage across TP53 in my exome?"

Returns coverage statistics without visualization.

### List Available Chromosomes

> "What contigs are in this BAM file?"

---

## Tools Reference

| Tool | Description | Required Args | Optional Args |
|------|-------------|---------------|---------------|
| `browse_region` | View aligned reads with interactive UI | `file_path`, `region` | `reference` |
| `get_variants` | Detect and return variants in a region | `file_path`, `region` | `reference`, `min_vaf`, `min_depth` |
| `get_coverage` | Calculate depth statistics | `file_path`, `region` | `reference` |
| `list_contigs` | List chromosomes in BAM/CRAM header | `file_path` | `reference` |

### Region Format

Regions can be specified as:
- `chr1:1000-2000` — Standard format
- `chr1:1,000-2,000` — Commas allowed
- `1:1000-2000` — Without "chr" prefix (depends on BAM header)

---

## Configuration

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `BAMCP_REFERENCE` | Path to reference FASTA (required for CRAM) | None |
| `BAMCP_MAX_READS` | Maximum reads to fetch per region | `10000` |
| `BAMCP_DEFAULT_WINDOW` | Default viewing window size (bp) | `500` |
| `BAMCP_MIN_VAF` | Minimum variant allele frequency to report | `0.1` |
| `BAMCP_MIN_DEPTH` | Minimum depth for variant calls | `10` |
| `BAMCP_MIN_MAPQ` | Minimum mapping quality filter | `0` |

---

## Implementation Details

### Why Canvas over React?

React's virtual DOM reconciliation struggles with thousands of elements. Genomic regions typically contain 1,000-10,000+ reads:

| Approach | Max Reads | DOM Nodes | Frame Rate |
|----------|-----------|-----------|------------|
| React/SVG | ~500 | O(n) | Degrades |
| Canvas 2D | ~10,000 | O(1) | 60fps |
| WebGL | ~100,000+ | O(1) | 60fps |

BAMCP uses Canvas 2D for the rendering layer with vanilla JS for controls.

### Read Packing Algorithm

Reads are packed into rows using a greedy algorithm:

```
For each read (sorted by start position):
    Find first row where read.start > row.end
    If found: add to that row
    Else: create new row
```

This minimizes vertical space while preventing overlaps.

### Variant Detection

BAMCP performs simple pileup-based variant detection:

1. Count bases at each position from aligned reads
2. Calculate variant allele frequency (VAF) = alt_count / total_depth
3. Report positions where VAF ≥ threshold and depth ≥ minimum

This is intentionally basic—for production variant calling, use dedicated tools (GATK, DeepVariant, etc.).

### CRAM Support

CRAM files require a reference FASTA for decoding. Set `BAMCP_REFERENCE` or pass `reference` to each tool call.

---

## Development

### Project Structure

```
bamcp/
├── src/
│   └── bamcp/
│       ├── __init__.py
│       ├── __main__.py       # Entry point
│       ├── server.py         # MCP server setup
│       ├── tools.py          # Tool handlers
│       ├── parsers.py        # pysam wrappers
│       ├── resources.py      # UI resource provider
│       └── config.py         # Configuration
├── tests/
│   ├── test_parsers.py
│   ├── test_tools.py
│   └── fixtures/
│       ├── small.bam
│       ├── small.bam.bai
│       └── ref.fa
├── docs/
│   └── assets/
│       └── demo.png
├── pyproject.toml
├── README.md
├── LICENSE
└── CONTRIBUTING.md
```

### Running Locally

```bash
# Run server directly
python -m bamcp

# Run with MCP inspector
npx @modelcontextprotocol/inspector python -m bamcp
```

### Running Tests

```bash
# All tests
pytest

# With coverage
pytest --cov=bamcp

# Specific test file
pytest tests/test_parsers.py -v
```

### Code Style

```bash
# Format
black src tests
isort src tests

# Lint
ruff check src tests

# Type check
mypy src
```

---

## Roadmap

- [x] **v0.1** — Core viewer with BAM support
- [ ] **v0.2** — CRAM support, coverage track
- [ ] **v0.3** — Variant detection and highlighting
- [ ] **v0.4** — VCF overlay support
- [ ] **v0.5** — Gene annotation track (RefSeq/Ensembl)
- [ ] **v0.6** — Remote file support (HTTP, S3)
- [ ] **v0.7** — Split-screen comparison view
- [ ] **v1.0** — Stable release

---

## Related Projects

- [bio-mcp-samtools](https://github.com/bio-mcp/bio-mcp-samtools) — CLI-based samtools MCP wrapper
- [igv.js](https://github.com/igvteam/igv.js) — JavaScript genome visualization library
- [pysam](https://github.com/pysam-developers/pysam) — Python interface to htslib
- [MCP Apps Extension](https://github.com/modelcontextprotocol/ext-apps) — UI extension for MCP

---

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- The [Model Context Protocol](https://modelcontextprotocol.io) team at Anthropic
- [MCP-UI](https://github.com/idosal/mcp-ui) for pioneering interactive MCP interfaces
- The [pysam](https://github.com/pysam-developers/pysam) and [htslib](https://github.com/samtools/htslib) maintainers
- [IGV](https://igv.org/) for inspiration on genomics visualization

---

<p align="center">
  <b>BAMCP</b> — BAM files + Model Context Protocol<br>
  <sub>Built with 🧬 for the computational biology community</sub>
</p>
