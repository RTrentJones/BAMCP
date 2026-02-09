import { Read, RegionData, Variant } from "./types";

export class StateManager {
    public data: RegionData | null = null;
    public viewport = { start: 0, end: 1000 };
    public packedRows: Read[][] = [];

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

    constructor() { }

    public loadData(data: RegionData): void {
        this.data = data;
        this.viewport = { start: data.start, end: data.end };

        this.buildMateIndex();
        this.packReads();
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

    private packReads(): void {
        if (!this.data) return;

        this.packedRows = [];
        this.readRowIndex.clear();
        const rowEnds: number[] = [];

        // Group reads by name to identify pairs (for packing)
        const readsByName = new Map<string, Read[]>();
        for (const read of this.data.reads) {
            const existing = readsByName.get(read.name);
            if (existing) existing.push(read);
            else readsByName.set(read.name, [read]);
        }

        // Sort reads by position, but process pairs together
        const sortedReads: Read[] = [];
        const processed = new Set<string>();
        const allReads = [...this.data.reads].sort((a, b) => a.position - b.position);

        for (const read of allReads) {
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
            let rowIndex = rowEnds.findIndex(end => end < read.position);

            if (rowIndex === -1) {
                rowIndex = rowEnds.length;
                rowEnds.push(read.end_position);
                this.packedRows.push([]);
            } else {
                rowEnds[rowIndex] = read.end_position;
            }

            this.packedRows[rowIndex].push(read);
            this.readRowIndex.set(read, rowIndex);
        }
    }

    public getFilteredAndSortedVariants(): Variant[] {
        if (!this.data) return [];

        let variants = [...this.data.variants];

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
