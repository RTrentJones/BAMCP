import { Read } from "./types";

type PackRequest = {
    type: "pack-reads";
    reads: Read[];
};

type PackResponse = {
    type: "packed-reads";
    rows: Read[][];
};

type HeapItem = { end: number; rowIndex: number };

function heapPush(heap: HeapItem[], item: HeapItem): void {
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

function heapPop(heap: HeapItem[]): HeapItem | undefined {
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

function packReads(reads: Read[]): Read[][] {
    const rows: Read[][] = [];
    const rowHeap: HeapItem[] = [];

    for (const read of reads) {
        const reusable = rowHeap.length > 0 && rowHeap[0].end < read.position
            ? heapPop(rowHeap)
            : undefined;
        const rowIndex = reusable ? reusable.rowIndex : rows.length;
        if (!reusable) rows.push([]);
        rows[rowIndex].push(read);
        heapPush(rowHeap, { end: read.end_position, rowIndex });
    }

    return rows;
}

self.onmessage = (event: MessageEvent<PackRequest>) => {
    if (event.data.type !== "pack-reads") return;
    const response: PackResponse = {
        type: "packed-reads",
        rows: packReads(event.data.reads),
    };
    self.postMessage(response);
};
