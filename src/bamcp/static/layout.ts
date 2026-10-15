/**
 * CIGAR-aware read layout: where each aligned base, insertion and deletion of a
 * read sits on the reference axis. Pure (no DOM) so it is unit-testable.
 *
 * `read.sequence` is the full query sequence (soft clips included), and
 * `read.position` is the 0-based reference start of the first aligned base.
 * Painting `sequence[i]` at `position + i` (the old behaviour) shifts every base
 * after a soft clip, insertion or deletion; this layout is the correct mapping.
 */

export interface CigarOp {
    op: string;
    len: number;
}

export interface AlignedBase {
    refPos: number;  // 0-based reference position
    qpos: number;    // index into read.sequence / read.qualities
}

export interface ReadLayout {
    bases: AlignedBase[];
    /** Inserted bases sit *before* refPos (between refPos-1 and refPos). */
    insertions: Array<{ refPos: number; len: number }>;
    /** D (deletion) or N (reference skip, e.g. an intron) spans, [start, end). */
    deletions: Array<{ start: number; end: number; skip: boolean }>;
}

const CIGAR_RE = /(\d+)([MIDNSHP=X])/g;

/** Parse a CIGAR string; returns [] for empty/unparseable input ("*", garbage). */
export function parseCigar(cigar: string | null | undefined): CigarOp[] {
    if (!cigar || typeof cigar !== 'string') return [];
    const ops: CigarOp[] = [];
    let consumed = 0;
    for (const m of cigar.matchAll(CIGAR_RE)) {
        ops.push({ op: m[2], len: Number(m[1]) });
        consumed += m[0].length;
    }
    // Any unparsed characters mean the string is malformed — don't half-trust it.
    return consumed === cigar.length ? ops : [];
}

/**
 * Lay a read out on the reference axis. Falls back to a gapless block spanning
 * [position, end_position) when the CIGAR is missing or malformed, so a read
 * always renders.
 */
export function layoutRead(read: {
    cigar?: string | null;
    position: number;
    end_position: number;
}): ReadLayout {
    const out: ReadLayout = { bases: [], insertions: [], deletions: [] };
    const ops = parseCigar(read.cigar);
    if (ops.length === 0) {
        const span = Math.max(1, read.end_position - read.position);
        for (let i = 0; i < span; i++) out.bases.push({ refPos: read.position + i, qpos: i });
        return out;
    }
    let ref = read.position;
    let q = 0;
    for (const { op, len } of ops) {
        switch (op) {
            case 'M':
            case '=':
            case 'X':
                for (let i = 0; i < len; i++) out.bases.push({ refPos: ref + i, qpos: q + i });
                ref += len;
                q += len;
                break;
            case 'I':
                out.insertions.push({ refPos: ref, len });
                q += len;
                break;
            case 'D':
            case 'N':
                out.deletions.push({ start: ref, end: ref + len, skip: op === 'N' });
                ref += len;
                break;
            case 'S':
                q += len;
                break;
            default: // H, P consume neither
                break;
        }
    }
    return out;
}

const layoutCache = new WeakMap<object, ReadLayout>();

/** layoutRead, memoized per read object (render runs on every hover/pan). */
export function cachedLayout(read: {
    cigar?: string | null;
    position: number;
    end_position: number;
}): ReadLayout {
    let layout = layoutCache.get(read);
    if (!layout) {
        layout = layoutRead(read);
        layoutCache.set(read, layout);
    }
    return layout;
}
