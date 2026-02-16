/**
 * Tile-based data cache with variant accumulation.
 *
 * Stores multiple RegionData "tiles" and accumulates variants across
 * fetches so panning/zooming doesn't lose previously discovered variants.
 * Clears everything on contig change.
 */

import { RegionData, Variant, VariantEvidence } from "./types";

interface Tile {
    data: RegionData;
    lastUsed: number;
}

export class DataStore {
    private tiles = new Map<string, Tile>();
    private variants = new Map<string, Variant>();
    private evidence = new Map<string, VariantEvidence>();
    private currentContig: string | null = null;

    static readonly MAX_TILES = 20;
    private static readonly OVERLAP_REPLACE_THRESHOLD = 0.9;
    private static readonly BEST_TILE_MIN_OVERLAP = 0.8;

    private static tileKey(contig: string, start: number, end: number): string {
        return `${contig}:${start}-${end}`;
    }

    private static variantKey(v: Variant): string {
        return `${v.contig}:${v.position}:${v.ref}>${v.alt}`;
    }

    /**
     * Ingest a new RegionData response. Adds to tile cache, merges variants.
     * Clears everything if contig changes.
     */
    public ingest(data: RegionData): void {
        // Contig change — full reset
        if (this.currentContig !== null && data.contig !== this.currentContig) {
            this.clear();
        }
        this.currentContig = data.contig;

        const key = DataStore.tileKey(data.contig, data.start, data.end);

        // Replace overlapping tiles (>90% overlap) to prevent fragmentation
        for (const [existingKey, tile] of this.tiles) {
            if (existingKey === key) continue;
            const overlap = this.overlapRatio(
                tile.data.start, tile.data.end,
                data.start, data.end,
            );
            if (overlap >= DataStore.OVERLAP_REPLACE_THRESHOLD) {
                this.tiles.delete(existingKey);
            }
        }

        // Add tile
        this.tiles.set(key, { data, lastUsed: Date.now() });

        // Evict LRU if over limit
        while (this.tiles.size > DataStore.MAX_TILES) {
            let oldestKey: string | null = null;
            let oldestTime = Infinity;
            for (const [k, t] of this.tiles) {
                if (t.lastUsed < oldestTime) {
                    oldestTime = t.lastUsed;
                    oldestKey = k;
                }
            }
            if (oldestKey) this.tiles.delete(oldestKey);
            else break;
        }

        // Merge variants (newer data wins for same key)
        for (const v of data.variants) {
            this.variants.set(DataStore.variantKey(v), v);
        }

        // Merge variant evidence (keyed by "position:ref>alt")
        if (data.variant_evidence) {
            for (const [key, ev] of Object.entries(data.variant_evidence)) {
                this.evidence.set(key, ev);
            }
        }
    }

    /**
     * Find the best cached tile for a viewport. Returns null if no tile
     * covers at least 80% of the requested range.
     */
    public bestTile(contig: string, start: number, end: number): RegionData | null {
        if (contig !== this.currentContig) return null;

        let best: Tile | null = null;
        let bestOverlap = 0;

        for (const tile of this.tiles.values()) {
            const overlap = this.overlapRatio(start, end, tile.data.start, tile.data.end);
            if (overlap > bestOverlap) {
                bestOverlap = overlap;
                best = tile;
            }
        }

        if (best && bestOverlap >= DataStore.BEST_TILE_MIN_OVERLAP) {
            best.lastUsed = Date.now();
            return best.data;
        }
        return null;
    }

    /** All accumulated variants for the current contig. */
    public getAllVariants(): Variant[] {
        return Array.from(this.variants.values());
    }

    /** Look up accumulated evidence for a variant by key ("position:ref>alt"). */
    public getEvidence(key: string): VariantEvidence | undefined {
        return this.evidence.get(key);
    }

    /** Number of cached tiles. */
    public get tileCount(): number {
        return this.tiles.size;
    }

    /** Full reset. */
    public clear(): void {
        this.tiles.clear();
        this.variants.clear();
        this.evidence.clear();
        this.currentContig = null;
    }

    /**
     * Fraction of [aStart, aEnd) covered by [bStart, bEnd),
     * relative to the first range's span.
     */
    private overlapRatio(
        aStart: number, aEnd: number,
        bStart: number, bEnd: number,
    ): number {
        const span = aEnd - aStart;
        if (span <= 0) return 0;
        const overlapStart = Math.max(aStart, bStart);
        const overlapEnd = Math.min(aEnd, bEnd);
        return Math.max(0, overlapEnd - overlapStart) / span;
    }
}
