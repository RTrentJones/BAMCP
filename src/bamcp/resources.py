"""UI resource provider for the MCP Apps extension."""

from pathlib import Path


def get_viewer_html() -> str:
    """Return the HTML content for the alignment viewer."""
    return """<!DOCTYPE html>
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
            <button id="zoom-out">\u2212</button>
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
    const BASE_COLORS = {
        'A': '#22c55e',
        'T': '#ef4444',
        'G': '#f97316',
        'C': '#3b82f6',
        'N': '#9ca3af'
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
            document.getElementById('go-btn').addEventListener('click', () => {
                const region = document.getElementById('region-input').value;
                this.requestRegion(region);
            });

            document.getElementById('zoom-in').addEventListener('click', () => this.zoom(0.5));
            document.getElementById('zoom-out').addEventListener('click', () => this.zoom(2));

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

            this.readsCanvas.addEventListener('mousemove', (e) => {
                if (isDragging) return;
                this.showTooltip(e);
            });

            this.readsCanvas.addEventListener('mouseout', () => {
                this.tooltip.style.display = 'none';
            });

            this.readsCanvas.addEventListener('wheel', (e) => {
                e.preventDefault();
                const factor = e.deltaY > 0 ? 1.2 : 0.8;
                this.zoom(factor);
            });
        }

        setupMCPBridge() {
            window.addEventListener('message', (event) => {
                const message = event.data;
                if (message.method === 'bamcp/init') {
                    this.loadData(message.params);
                } else if (message.method === 'bamcp/update') {
                    this.loadData(message.params);
                }
            });

            window.parent.postMessage({
                jsonrpc: '2.0',
                method: 'ready'
            }, '*');
        }

        requestRegion(region) {
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
                data.contig + ':' + data.start + '-' + data.end;

            this.packReads();
            this.renderVariantTable();
            this.render();
        }

        packReads() {
            if (!this.data) return;

            this.packedRows = [];
            const rowEnds = [];

            const sortedReads = [...this.data.reads].sort((a, b) => a.position - b.position);

            for (const read of sortedReads) {
                let rowIndex = rowEnds.findIndex(end => end < read.position);

                if (rowIndex === -1) {
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
            this.readsCanvas.height = Math.max(height, 100);

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

            ctx.lineTo((this.data.end - this.viewport.start) * scale, height);
            ctx.lineTo(0, height);
            ctx.closePath();
            ctx.fill();
            ctx.stroke();

            ctx.fillStyle = '#333';
            ctx.font = '10px sans-serif';
            ctx.fillText('Max: ' + maxCov + 'x', 5, 12);
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

                if (y > height) break;

                for (const read of this.packedRows[rowIdx]) {
                    this.renderRead(ctx, read, y, scale);
                }
            }
        }

        renderRead(ctx, read, y, scale) {
            const x = (read.position - this.viewport.start) * scale;
            const w = (read.end_position - read.position) * scale;

            if (x + w < 0 || x > this.readsCanvas.width) return;

            ctx.fillStyle = read.is_reverse ? '#a78bfa' : '#60a5fa';
            ctx.fillRect(x, y, w, READ_HEIGHT);

            if (scale > 2) {
                for (const mm of read.mismatches) {
                    const mx = (mm.pos - this.viewport.start) * scale;
                    ctx.fillStyle = BASE_COLORS[mm.alt] || '#9ca3af';
                    ctx.fillRect(mx, y, Math.max(scale, 2), READ_HEIGHT);
                }
            }

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

            const read = this.packedRows[rowIdx].find(r =>
                r.position <= genomePos && r.end_position >= genomePos
            );

            if (read) {
                this.tooltip.innerHTML =
                    '<strong>' + read.name + '</strong><br>' +
                    'Position: ' + read.position + '-' + read.end_position + '<br>' +
                    'MAPQ: ' + read.mapping_quality + '<br>' +
                    'CIGAR: ' + read.cigar + '<br>' +
                    'Strand: ' + (read.is_reverse ? 'Reverse' : 'Forward');
                this.tooltip.style.display = 'block';
                this.tooltip.style.left = (e.clientX + 10) + 'px';
                this.tooltip.style.top = (e.clientY + 10) + 'px';
            } else {
                this.tooltip.style.display = 'none';
            }
        }

        renderVariantTable() {
            if (!this.data) return;

            this.variantTable.innerHTML = this.data.variants.map(function(v) {
                return '<tr data-pos="' + v.position + '">' +
                    '<td>' + v.contig + ':' + v.position + '</td>' +
                    '<td>' + v.ref + '</td>' +
                    '<td style="color: ' + (BASE_COLORS[v.alt] || '#333') + '">' + v.alt + '</td>' +
                    '<td>' + (v.vaf * 100).toFixed(1) + '%</td>' +
                    '<td>' + v.depth + '</td>' +
                    '</tr>';
            }).join('');

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

    const viewer = new BAMCPViewer();
    </script>
</body>
</html>"""
