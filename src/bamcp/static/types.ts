export interface Read {
    name: string;
    sequence: string;
    qualities: number[];
    cigar: string;
    position: number;
    end_position: number;
    mapping_quality: number;
    is_reverse: boolean;
    mismatches: Array<{ pos: number; ref: string; alt: string }>;
    // Paired-end fields
    mate_position?: number | null;
    mate_contig?: string | null;
    insert_size?: number | null;
    is_proper_pair?: boolean;
    is_read1?: boolean;
    is_paired?: boolean;
}

export interface Variant {
    contig: string;
    position: number;
    ref: string;
    alt: string;
    vaf: number;
    depth: number;
    // Enhanced fields from tools.py
    strand_forward?: number;
    strand_reverse?: number;
    mean_quality?: number;
    is_low_confidence?: boolean;
    alt_count?: number;
}

export interface VariantEvidence {
    forward_count: number;
    reverse_count: number;
    strand_bias: number;
    mean_quality: number;
    median_quality: number;
    quality_histogram: Record<string, number>;
    position_histogram: Record<string, number>;
}

export interface RegionData {
    contig: string;
    start: number;
    end: number;
    reads: Read[];
    coverage: number[];
    variants: Variant[];
    reference_sequence?: string;
    variant_evidence?: Record<string, VariantEvidence>;
}

export interface ToolResultParams {
    structuredContent?: RegionData;
    content?: Array<{ type: string; text?: string }>;
}
