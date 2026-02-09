import { App } from "@modelcontextprotocol/ext-apps";
import { RegionData, ToolResultParams, Variant } from "./types";

export class BAMCPClient {
    private app: App | null = null;
    private onDataCallback: ((data: RegionData) => void) | null = null;

    constructor() {
        // Initialize callback
    }

    public setOnDataReceived(callback: (data: RegionData) => void): void {
        this.onDataCallback = callback;
    }

    public async init(): Promise<void> {
        try {
            this.app = new App({ name: 'BAMCP Viewer', version: '1.0.0' });

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

            await this.app.connect();
        } catch (e) {
            console.warn('MCP Apps SDK init failed:', e);
        }
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
