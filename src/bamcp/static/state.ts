import { DataStore } from "./data-store";
import { ColumnarReads, Read, RegionData, Variant, ViewerSettings } from "./types";

/**
 * Decode columnar reads (parallel arrays) back into individual Read objects.
 * An already-decoded Read[] is returned unchanged, so direct callers and tests
 * that pass the object form keep working.
 */
export function decodeReads(raw: ColumnarReads | Read[] | null | undefined): Read[] {
    if (!raw) return [];
    if (Array.isArray(raw)) return raw;
    const c = raw;
    const reads: Read[] = new Array(c.count);
    for (let i = 0; i < c.count; i++) {
        const read: Read = {
            name: c.name[i],
            cigar: c.cigar[i],
            position: c.position[i],
            end_position: c.end_position[i],
            mapping_quality: c.mapping_quality[i],
            is_reverse: c.is_reverse[i],
            mismatches: c.mismatches[i] ?? [],
        };
        const seq = c.sequence?.[i];
        if (seq != null) read.sequence = seq;
        const quals = c.qualities?.[i];
        if (quals != null) read.qualities = quals;
        const softClips = c.soft_clips?.[i];
        if (softClips && softClips.length) read.soft_clips = softClips;
        if (c.is_paired?.[i]) {
            read.is_paired = true;
            read.mate_position = c.mate_position?.[i] ?? null;
            read.mate_contig = c.mate_contig?.[i] ?? null;
            read.insert_size = c.insert_size?.[i] ?? null;
            read.is_proper_pair = c.is_proper_pair?.[i] ?? false;
            read.is_read1 = c.is_read1?.[i] ?? false;
        }
        reads[i] = read;
    }
    return reads;
}

/** Return a RegionData whose reads are decoded to Read[] (columnar or already-array). */
export function decodeRegionData(data: RegionData): RegionData {
    const raw = data.reads as unknown as ColumnarReads | Read[];
    const reads = decodeReads(raw);
    return reads === raw ? data : { ...data, reads };
}

// Default viewer settings
export const DEFAULT_SETTINGS: ViewerSettings = {
    displayMode: 'compact',
    colorBy: 'strand',
    sortBy: 'position',
    showSoftClips: false,
    showMismatches: true,
    activeVariantPosition: null,
    activeVariantAlt: null,
};

export class StateManager {
    public data: RegionData | null = null;
    public viewport = { start: 0, end: 1000 };
    public packedRows: Read[][] = [];
    public store = new DataStore();

    // Maps
    public mateIndex: Map<string, Read> = new Map();
    public readRowIndex: Map<Read, number> = new Map();

    // UI State
    public hoveredRead: Read | null = null;
    public pendingHoverRead: Read | null = null;
    public lockedTooltip: { read: Read; x: number; y: number } | null = null;

    // Variant state
    public variantFilter: 'high' | 'all' = 'high';
    public variantSort: { column: string; direction: 'asc' | 'desc' } = { column: 'position', direction: 'asc' };
    public selectedVariantIndex: number = -1;
    public expandedVariantIndex: number = -1;

    // Viewer settings (IGV-style display options)
    public settings: ViewerSettings = { ...DEFAULT_SETTINGS };

    constructor() { }

    /** Load new data from host (ontoolresult) — resets viewport. */
    public loadData(data: RegionData): void {
        data = decodeRegionData(data);
        this.store.ingest(data);
        this.data = data;
        this.viewport = { start: data.start, end: data.end };

        this.buildMateIndex();
        this.packReads();
    }

    /** Load a tile from viewport refetch — preserves current viewport. */
    public loadTile(data: RegionData): void {
        data = decodeRegionData(data);
        this.store.ingest(data);
        this.data = data;

        this.buildMateIndex();
        this.packReads();
    }

    /** Activate an already-cached tile without fetching. */
    public activateTile(data: RegionData): void {
        this.data = decodeRegionData(data);
        this.buildMateIndex();
        this.packReads();
    }

    /**
     * Re-sort and repack reads when sort settings change.
     */
    public resortAndRepack(): void {
        if (!this.data) return;
        this.packReads();
    }

    /**
     * Get sorted reads based on current sort settings.
     */
    private getSortedReads(): Read[] {
        if (!this.data) return [];

        const reads = [...this.data.reads];

        switch (this.settings.sortBy) {
            case 'mapq':
                return reads.sort((a, b) => b.mapping_quality - a.mapping_quality);
            case 'insertSize':
                return reads.sort((a, b) =>
                    Math.abs(a.insert_size || 0) - Math.abs(b.insert_size || 0)
                );
            case 'strand':
                return reads.sort((a, b) =>
                    (a.is_reverse ? 1 : 0) - (b.is_reverse ? 1 : 0)
                );
            default: // position
                return reads.sort((a, b) => a.position - b.position);
        }
    }

    private buildMateIndex(): void {
        this.mateIndex.clear();
        if (!this.data) return;

        // Group reads by name
        const readsByName = new Map<string, Read[]>();
        for (const read of this.data.reads) {
            const existing = readsByName.get(read.name);
            if (existing) {
                existing.push(read);
            } else {
                readsByName.set(read.name, [read]);
            }
        }

        // For pairs where both are loaded, create mate links
        for (const [name, reads] of readsByName) {
            if (reads.length === 2) {
                const [r1, r2] = reads;
                this.mateIndex.set(`${name}:${r1.position}`, r2);
                this.mateIndex.set(`${name}:${r2.position}`, r1);
            }
        }
    }

    private heapPush(heap: Array<{ end: number; rowIndex: number }>, item: { end: number; rowIndex: number }): void {
        heap.push(item);
        let i = heap.length - 1;
        while (i > 0) {
            const parent = Math.floor((i - 1) / 2);
            if (heap[parent].end <= item.end) break;
            heap[i] = heap[parent];
            i = parent;
        }
        heap[i] = item;
    }

    private heapPop(heap: Array<{ end: number; rowIndex: number }>): { end: number; rowIndex: number } | undefined {
        if (heap.length === 0) return undefined;
        const root = heap[0];
        const last = heap.pop()!;
        if (heap.length > 0) {
            let i = 0;
            while (true) {
                const left = 2 * i + 1;
                const right = left + 1;
                if (left >= heap.length) break;
                const child = right < heap.length && heap[right].end < heap[left].end ? right : left;
                if (heap[child].end >= last.end) break;
                heap[i] = heap[child];
                i = child;
            }
            heap[i] = last;
        }
        return root;
    }

    private packReads(): void {
        if (!this.data) return;

        this.packedRows = [];
        this.readRowIndex.clear();
        const rowHeap: Array<{ end: number; rowIndex: number }> = [];

        // Group reads by name to identify pairs (for packing)
        const readsByName = new Map<string, Read[]>();
        for (const read of this.data.reads) {
            const existing = readsByName.get(read.name);
            if (existing) existing.push(read);
            else readsByName.set(read.name, [read]);
        }

        // Get reads sorted according to current settings
        const baseSortedReads = this.getSortedReads();

        // Process reads, keeping pairs together
        const sortedReads: Read[] = [];
        const processed = new Set<string>();

        for (const read of baseSortedReads) {
            const key = `${read.name}:${read.position}`;
            if (processed.has(key)) continue;

            const pair = readsByName.get(read.name);
            if (pair && pair.length === 2) {
                // Add both reads of the pair consecutively
                const [r1, r2] = pair.sort((a, b) => a.position - b.position);
                sortedReads.push(r1, r2);
                processed.add(`${r1.name}:${r1.position}`);
                processed.add(`${r2.name}:${r2.position}`);
            } else {
                sortedReads.push(read);
                processed.add(key);
            }
        }

        // Pack into rows
        for (const read of sortedReads) {
            const reusable = rowHeap.length > 0 && rowHeap[0].end < read.position
                ? this.heapPop(rowHeap)
                : undefined;

            let rowIndex: number;
            if (!reusable) {
                rowIndex = this.packedRows.length;
                this.packedRows.push([]);
            } else {
                rowIndex = reusable.rowIndex;
            }

            this.packedRows[rowIndex].push(read);
            this.readRowIndex.set(read, rowIndex);
            this.heapPush(rowHeap, { end: read.end_position, rowIndex });
        }
    }

    public getFilteredAndSortedVariants(): Variant[] {
        if (!this.data) return [];

        // Use accumulated variants from DataStore (persists across tile fetches)
        let variants = this.store.getAllVariants();

        // Filter: use backend-computed confidence field (high/medium/low)
        if (this.variantFilter === 'high') {
            variants = variants.filter(v =>
                v.confidence === 'high' ||
                v.confidence === 'medium' ||
                v.confidence === undefined  // Include if field not set (backwards compat)
            );
        }

        // Sort
        variants.sort((a, b) => {
            const valA = (a as any)[this.variantSort.column];
            const valB = (b as any)[this.variantSort.column];

            if (typeof valA === 'string') {
                return this.variantSort.direction === 'asc'
                    ? valA.localeCompare(valB)
                    : valB.localeCompare(valA);
            } else {
                return this.variantSort.direction === 'asc'
                    ? valA - valB
                    : valB - valA;
            }
        });

        return variants;
    }
}
