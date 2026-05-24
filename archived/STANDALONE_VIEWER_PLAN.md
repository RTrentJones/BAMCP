# Plan: Standalone Web Viewer at `/standalone`

## Context

BAMCP's viewer currently only works inside MCP hosts (Claude.ai, Claude Desktop) via the `@modelcontextprotocol/ext-apps` SDK. Users want to open `bamcp.rtrentjones.dev/standalone` directly in a browser to visualize BAM/CRAM files — like [igv.org/app](https://igv.org/app/). This requires a thin REST API layer over existing tool handlers and a standalone frontend entry point that uses `fetch()` instead of the MCP Apps SDK.

## Architecture

```
Browser → GET /standalone → standalone.html (self-contained)
       → POST /api/visualize → handle_visualize_region() → RegionData JSON
       → POST /api/gene → GeneClient.search() → coordinates JSON
       → POST /api/clinvar → handle_lookup_clinvar() → JSON
       → POST /api/gnomad → handle_lookup_gnomad() → JSON
       → POST /api/contigs → handle_list_contigs() → JSON
       → POST /api/summary → handle_get_region_summary() → JSON
```

Shareable URLs: `bamcp.rtrentjones.dev/standalone?file=https://...&region=chr1:1000-2000`

## Key Design Decisions

1. **Client abstraction**: Extract `ViewerClient` interface from `BAMCPClient`. Create `StandaloneClient` implementing same interface with `fetch()`. `BAMCPViewer` constructor accepts either client.
2. **REST API**: New `api.py` module — thin Starlette routes wrapping existing `handle_*()` functions from `core/tools.py`. No new business logic.
3. **Auth**: REST API endpoints are public (rate-limited, no OAuth). The standalone page is unauthenticated.
4. **Build**: Two separate Vite builds (singlefile plugin requires `inlineDynamicImports: true`, incompatible with multi-entry). npm script runs both.
5. **Reuse**: renderer.ts, state.ts, types.ts, constants.ts unchanged. All tool handlers unchanged.

## New Files (7)

### `src/bamcp/api.py` — REST API routes

Thin wrappers around existing handlers. Pattern for each endpoint:

```python
async def visualize(request: Request) -> JSONResponse:
    body = await request.json()
    result = await handle_visualize_region(body, config)
    payload = result.get("_meta", {}).get("ui/init", {})
    return JSONResponse(payload)
```

Endpoints: `/api/visualize`, `/api/contigs`, `/api/gene`, `/api/clinvar`, `/api/gnomad`, `/api/summary`

Error handling: `ValueError` → 400, unexpected → 500 with sanitized message.

### `src/bamcp/standalone.py` — Page serving

Serves `dist/standalone.html` via `GET /standalone`. Same fallback chain as `resources.py`.

### `src/bamcp/static/client-interface.ts` — ViewerClient interface

```typescript
export interface ViewerClient {
    setOnDataReceived(cb: (data: RegionData) => void): void;
    init(): Promise<void>;
    fetchRegionDirect(filePath: string, region: string, ref?: string): Promise<void>;
    searchGene(symbol: string): Promise<void>;
    lookupClinVar(variant: Variant): Promise<any>;
    lookupGnomAD(variant: Variant): Promise<any>;
    // ... other methods BAMCPViewer calls
}
```

### `src/bamcp/static/http-client.ts` — StandaloneClient

Implements `ViewerClient` using `fetch()` to `/api/*` endpoints. Reads URL query params on init (`?file=`, `?region=`, `?ref=`). Updates browser URL via `history.replaceState()` on navigation.

LLM-dependent methods (`sendVariantMessage`, `updateModelContext`, `syncContext`) are no-ops.

### `src/bamcp/static/standalone.html` — HTML template

Based on `viewer.html` with:
- **Added**: File URL input bar (BAM URL + optional reference URL + Load button)
- **Added**: BAMCP branding header
- **Removed**: Sync button, Debug button (MCP-specific)
- **Modified**: ClinVar/gnomAD buttons call API directly, show results in evidence panel
- **Entry**: `<script type="module" src="./standalone-app.ts">`

### `src/bamcp/static/standalone-app.ts` — Entry point

Instantiates `BAMCPViewer` with `StandaloneClient` and wires standalone-specific UI (file input bar, Load button).

### `tests/unit/test_api.py` — REST API tests

Starlette `TestClient` tests for each `/api/*` endpoint using existing test fixtures.

## Modified Files (6)

### `src/bamcp/__main__.py`
Mount API routes and standalone page for HTTP transports:
```python
from .api import create_api_routes
from .standalone import standalone_route
app.routes.insert(0, standalone_route)
for route in create_api_routes(config):
    app.routes.insert(0, route)
```

### `src/bamcp/static/mcp-app.ts`
Change constructor to accept optional client parameter:
```typescript
constructor(client?: ViewerClient) {
    this.client = client || new BAMCPClient();
```

### `src/bamcp/static/client.ts`
Add `implements ViewerClient` to `BAMCPClient` class declaration.

### `src/bamcp/static/vite.config.ts`
Keep as viewer-only build. Add separate `vite.standalone.config.ts` for standalone build.

### `src/bamcp/static/package.json`
Update build script: `"build": "vite build && vite build -c vite.standalone.config.ts"`

### `Dockerfile`
Change `COPY --from=frontend .../viewer.html` to `COPY --from=frontend /frontend/dist/ src/bamcp/static/dist/` to capture both built HTML files.

## Unchanged (fully reused)

- `src/bamcp/static/renderer.ts` — Canvas rendering
- `src/bamcp/static/state.ts` — State management
- `src/bamcp/static/types.ts` — Data types
- `src/bamcp/static/constants.ts` — Colors
- `src/bamcp/core/tools.py` — All tool handlers
- `src/bamcp/core/parsers.py` — BAM parsing
- `src/bamcp/clients/*` — ClinVar, gnomAD, gene clients
- `src/bamcp/middleware/*` — Rate limiting already protects all routes

## Implementation Order

1. **Backend API** (`api.py`, `standalone.py`, modify `__main__.py`) — testable independently
2. **Client interface** (`client-interface.ts`, modify `client.ts`) — no behavior change
3. **HTTP client** (`http-client.ts`) — testable against backend
4. **Standalone HTML + entry point** (`standalone.html`, `standalone-app.ts`)
5. **Build config** (`vite.standalone.config.ts`, update `package.json`, `Dockerfile`)
6. **Tests** (`test_api.py`)

## Verification

1. `python -m pytest tests/` — all existing tests pass
2. `npm run build` in `src/bamcp/static/` — produces both `dist/viewer.html` and `dist/standalone.html`
3. Run locally: `BAMCP_TRANSPORT=streamable-http BAMCP_ALLOW_REMOTE_FILES=true python -m bamcp`
4. Open `http://localhost:8000/standalone` — page loads
5. Enter `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam` + region `chr1:1000-2000` → reads render
6. Test shareable URL: `http://localhost:8000/standalone?file=...&region=chr1:1000-2000` — auto-loads
7. `curl -s http://localhost:8000/api/contigs -d '{"file_path": "..."}' -H 'Content-Type: application/json'` — returns JSON
8. Verify MCP viewer still works: `curl -sI https://bamcp.rtrentjones.dev/mcp` — 401 (auth still works)
9. Docker build succeeds with both HTML files included
