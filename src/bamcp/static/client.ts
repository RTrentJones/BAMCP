import { App } from "@modelcontextprotocol/ext-apps";
import { RegionData, ToolResultParams, Variant } from "./types";

export type DisplayMode = 'inline' | 'fullscreen' | 'pip';

export class BAMCPClient {
    private app: App | null = null;
    private onDataCallback: ((data: RegionData) => void) | null = null;
    private currentDisplayMode: DisplayMode = 'inline';
    private availableDisplayModes: DisplayMode[] = ['inline'];

    constructor() {
        // Initialize callback
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

            // Set ontoolresult BEFORE connect() to receive the initial tool result
            this.app.ontoolresult = (params: any) => {
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
        if (this.onDataCallback) {
            this.onDataCallback(data);
        }
    }

    public async requestRegion(region: string): Promise<void> {
        if (this.app) {
            try {
                await this.app.sendMessage({
                    role: 'user',
                    content: [{
                        type: 'text',
                        text: `Please browse the genomic region ${region} and show the alignment visualization.`
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
            await this.app.updateModelContext(context);
        } catch {
            // Ignore errors
        }
    }

    public async sendVariantMessage(variant: Variant): Promise<void> {
        if (!this.app) return;

        const message =
            'Explain the variant at ' + variant.contig + ':' + variant.position +
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
}
