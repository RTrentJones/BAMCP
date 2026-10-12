import { App } from "@modelcontextprotocol/ext-apps";
import { RegionData, ToolResultParams, Variant } from "./types";
import { EXPECTED_SCHEMA_VERSION, displayPos } from "./constants";

export type DisplayMode = 'inline' | 'fullscreen' | 'pip';

/**
 * Outcome of an app-initiated fetch. `stale` means a newer request (or host data)
 * superseded this one, so the caller must drop it silently; `error` carries a
 * short, user-presentable message.
 */
export type FetchResult =
    | { status: 'ok'; data: RegionData }
    | { status: 'stale' }
    | { status: 'error'; message: string };

/** First text block of a tool result, trimmed to something safe to show a user. */
function errorText(content: unknown): string {
    if (Array.isArray(content)) {
        for (const c of content) {
            if (c && typeof c === 'object' && 'text' in c && typeof c.text === 'string') {
                return c.text.slice(0, 300);
            }
        }
    }
    return 'the server returned an error';
}

export interface DebugInfo {
    lastContext: string | null;
    lastToolCall: string | null;
}

export class BAMCPClient {
    private app: App | null = null;
    private onDataCallback: ((data: RegionData) => void) | null = null;
    private currentDisplayMode: DisplayMode = 'inline';
    private availableDisplayModes: DisplayMode[] = ['inline'];

    // Debug tracking
    private debugInfo: DebugInfo = { lastContext: null, lastToolCall: null };
    private onDebugCallback: ((info: DebugInfo) => void) | null = null;

    // Monotonic version counter to discard stale fetchRegionDirect responses.
    // Bumped by ontoolresult (host-initiated data) so in-flight app-initiated
    // fetches don't overwrite newer data.
    private dataVersion = 0;
    // Monotonic id of the latest app-initiated fetch. Overlapping fetches (a pan
    // refetch racing a sequence fetch, or two refetches) can resolve out of order;
    // only the newest request's response is applied.
    private requestSeq = 0;

    constructor() {
        // Initialize callback
    }

    public setOnDebugUpdate(callback: (info: DebugInfo) => void): void {
        this.onDebugCallback = callback;
    }

    public getDebugInfo(): DebugInfo {
        return this.debugInfo;
    }

    private updateDebug(field: keyof DebugInfo, value: string): void {
        this.debugInfo[field] = value;
        if (this.onDebugCallback) {
            this.onDebugCallback(this.debugInfo);
        }
    }

    public getAvailableDisplayModes(): DisplayMode[] {
        return this.availableDisplayModes;
    }

    public getCurrentDisplayMode(): DisplayMode {
        return this.currentDisplayMode;
    }

    public setOnDataReceived(callback: (data: RegionData) => void): void {
        this.onDataCallback = callback;
    }

    public async init(): Promise<void> {
        try {
            // Constructor: (appInfo, capabilities, options)
            // Declare support for inline, fullscreen, and PIP modes
            this.app = new App(
                { name: 'BAMCP Viewer', version: '1.0.0' },
                { availableDisplayModes: ['inline', 'fullscreen', 'pip'] },
            );

            // Set ontoolresult BEFORE connect() to receive the initial tool result.
            // Bump dataVersion to invalidate any in-flight fetchRegionDirect calls.
            this.app.ontoolresult = (params: any) => {
                this.dataVersion++;
                const typedParams = params as ToolResultParams;
                if (typedParams.structuredContent) {
                    this.handleData(typedParams.structuredContent);
                } else if (typedParams.content) {
                    const text = typedParams.content.find(c => c.type === 'text');
                    if (text?.text) {
                        try {
                            this.handleData(JSON.parse(text.text));
                        } catch (e) {
                            console.warn('Failed to parse content text:', e);
                        }
                    }
                }
            };

            // Listen for host context updates (available display modes)
            this.app.onhostcontextchanged = (context: any) => {
                if (context?.availableDisplayModes) {
                    this.availableDisplayModes = context.availableDisplayModes;
                }
                if (context?.displayMode) {
                    this.currentDisplayMode = context.displayMode;
                }
            };

            await this.app.connect();

            // Check initial host context for available display modes
            const hostContext = this.app.getHostContext();
            if (hostContext?.availableDisplayModes) {
                this.availableDisplayModes = hostContext.availableDisplayModes as DisplayMode[];
            }
            if (hostContext?.displayMode) {
                this.currentDisplayMode = hostContext.displayMode as DisplayMode;
            }
        } catch (e) {
            console.warn('MCP Apps SDK init failed:', e);
        }
    }

    public async requestDisplayMode(mode: DisplayMode): Promise<boolean> {
        if (!this.app) {
            // Fallback to browser fullscreen API when not connected to MCP host
            return this.fallbackToggleFullscreen(mode);
        }

        try {
            // Use the MCP Apps API to request display mode change
            const result = await this.app.requestDisplayMode({ mode });
            this.currentDisplayMode = result.mode as DisplayMode;
            return result.mode === mode;
        } catch (e) {
            console.warn('Display mode change failed:', e);
            return this.fallbackToggleFullscreen(mode);
        }
    }

    private async fallbackToggleFullscreen(mode: DisplayMode): Promise<boolean> {
        // Fallback to browser fullscreen API
        if (mode === 'fullscreen') {
            try {
                await document.documentElement.requestFullscreen();
                this.currentDisplayMode = 'fullscreen';
                return true;
            } catch {
                return false;
            }
        } else if (mode === 'inline' && document.fullscreenElement) {
            await document.exitFullscreen();
            this.currentDisplayMode = 'inline';
            return true;
        }
        return false;
    }

    public async toggleFullscreen(): Promise<boolean> {
        const newMode = this.currentDisplayMode === 'fullscreen' ? 'inline' : 'fullscreen';
        return this.requestDisplayMode(newMode);
    }

    private handleData(data: RegionData): void {
        if (
            data?.schema_version !== undefined &&
            data.schema_version !== EXPECTED_SCHEMA_VERSION
        ) {
            console.warn(
                `BAMCP payload schema_version ${data.schema_version} != expected ` +
                `${EXPECTED_SCHEMA_VERSION}; viewer bundle may be out of step with the server.`,
            );
        }
        if (this.onDataCallback) {
            this.onDataCallback(data);
        }
    }

    public async requestRegion(region: string, filePath?: string): Promise<void> {
        if (this.app) {
            try {
                // Make the request explicit with file path for more reliable tool invocation
                const message = filePath
                    ? `Call visualize_region for file ${filePath} region ${region}`
                    : `Please browse the genomic region ${region} and show the alignment visualization.`;

                await this.app.sendMessage({
                    role: 'user',
                    content: [{
                        type: 'text',
                        text: message
                    }]
                });
            } catch (e) {
                console.warn('sendMessage failed:', e);
            }
        }
    }

    public async searchGene(symbol: string): Promise<void> {
        if (this.app) {
            try {
                await this.app.sendMessage({
                    role: 'user',
                    content: [{
                        type: 'text',
                        text: `Search for gene ${symbol} and show the alignment at that location.`
                    }]
                });
            } catch (e) {
                console.warn('Gene search failed:', e);
            }
        }
    }

    public async updateModelContext(context: any): Promise<void> {
        if (!this.app) return;
        try {
            // YAML frontmatter format recommended by MCP Apps patterns docs
            const contextText = `---
region: ${context.region}
reads: ${context.reads}
coverage: ${context.meanCoverage}x
variants: ${context.variantCount}
---

User is viewing ${context.region} with ${context.reads} reads at ${context.meanCoverage}x coverage and ${context.variantCount} variants visible.`;
            await this.app.updateModelContext({
                content: [{ type: 'text', text: contextText }]
            });
            this.updateDebug('lastContext', contextText);
        } catch (e) {
            this.updateDebug('lastContext', `ERROR: ${e}`);
            console.warn('updateModelContext failed:', e);
        }
    }

    /**
     * Directly invoke visualize_region tool using callServerTool.
     *
     * Responses are sequenced: if host data arrived (ontoolresult) or a newer
     * fetchRegionDirect was issued while this one was in flight, the result is
     * `stale` and must be dropped, so an older tile can never replace a newer one.
     * Callers decide how to load the data (loadTile vs loadData), keeping this
     * path separate from host-initiated ontoolresult.
     */
    public async fetchRegionDirect(
        filePath: string,
        region: string,
        reference?: string,
        vcfPath?: string,
        variantSource?: string,
    ): Promise<FetchResult> {
        if (!this.app) return { status: 'error', message: 'not connected to an MCP host' };
        const args: Record<string, string> = { file_path: filePath, region: region };
        if (reference) {
            args.reference = reference;
        }
        // Keep VCF-as-primary sourcing on pan/zoom refetches (else tiles revert to auto).
        if (vcfPath) {
            args.vcf_path = vcfPath;
        }
        if (variantSource) {
            args.variant_source = variantSource;
        }
        const toolCallInfo = `visualize_region(${JSON.stringify(args)})`;
        this.updateDebug('lastToolCall', toolCallInfo);

        // Capture version + sequence before await; either moving on means stale.
        const version = this.dataVersion;
        const seq = ++this.requestSeq;
        const isStale = () => this.dataVersion !== version || this.requestSeq !== seq;

        try {
            const result = await this.app.callServerTool({
                name: 'visualize_region',
                arguments: args
            });

            if (isStale()) {
                this.updateDebug('lastToolCall', `${toolCallInfo} → STALE (discarded)`);
                return { status: 'stale' };
            }

            if (result.isError) {
                this.updateDebug('lastToolCall', `${toolCallInfo} → ERROR`);
                console.warn('callServerTool returned error:', result.content);
                return { status: 'error', message: errorText(result.content) };
            }

            // Extract data from structuredContent or content
            let data: RegionData | null = null;
            if (result.structuredContent) {
                data = result.structuredContent as unknown as RegionData;
            } else if (result.content) {
                const text = result.content.find(c => c.type === 'text');
                if (text && 'text' in text) {
                    try {
                        data = JSON.parse(text.text);
                    } catch (e) {
                        console.warn('Failed to parse content text:', e);
                    }
                }
            }
            if (!data) {
                this.updateDebug('lastToolCall', `${toolCallInfo} → ERROR: unparseable`);
                return { status: 'error', message: 'unreadable response from the server' };
            }
            this.updateDebug('lastToolCall', `${toolCallInfo} → OK`);
            return { status: 'ok', data };
        } catch (e) {
            if (isStale()) return { status: 'stale' };
            this.updateDebug('lastToolCall', `${toolCallInfo} → ERROR: ${e}`);
            console.warn('callServerTool failed:', e);
            const message = e instanceof Error ? e.message : String(e);
            return { status: 'error', message: message.slice(0, 300) || 'request failed' };
        }
    }

    public async sendVariantMessage(variant: Variant): Promise<void> {
        if (!this.app) return;

        const message =
            'Explain the variant at ' + variant.contig + ':' + displayPos(variant.position) +
            ': ' + variant.ref + '>' + variant.alt +
            ', VAF=' + (variant.vaf * 100).toFixed(1) + '%' +
            ', depth=' + variant.depth;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{ type: 'text', text: message }]
            });
        } catch {
            // Ignore errors
        }
    }

    public async lookupClinVar(variant: Variant): Promise<void> {
        if (!this.app) return;

        const message =
            `Look up the clinical significance of ${variant.contig}:${displayPos(variant.position)} ` +
            `${variant.ref}>${variant.alt} in ClinVar using the lookup_clinvar tool.`;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{ type: 'text', text: message }]
            });
        } catch {
            // Ignore errors
        }
    }

    public async lookupGnomAD(variant: Variant): Promise<void> {
        if (!this.app) return;

        const message =
            `Look up the population allele frequency of ${variant.contig}:${displayPos(variant.position)} ` +
            `${variant.ref}>${variant.alt} in gnomAD using the lookup_gnomad tool.`;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{ type: 'text', text: message }]
            });
        } catch {
            // Ignore errors
        }
    }

    public async syncContext(region: string): Promise<void> {
        if (!this.app) return;

        const message =
            `Get a summary of the current view at region ${region} ` +
            `using get_region_summary with the same BAM file.`;

        try {
            await this.app.sendMessage({
                role: 'user',
                content: [{ type: 'text', text: message }]
            });
        } catch {
            // Ignore errors
        }
    }
}
