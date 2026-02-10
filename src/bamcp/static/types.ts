export interface Read {
    name: string;
    sequence?: string;  // Only in non-compact mode (high zoom)
    cigar: string;
    position: number;
    end_position: number;
    mapping_quality: number;
    is_reverse: boolean;
    mismatches: Array<{ pos: number; ref: string; alt: string }>;
    // Paired-end fields (only present if is_paired=true)
    mate_position?: number | null;
    mate_contig?: string | null;
    insert_size?: number | null;
    is_proper_pair?: boolean;
    is_read1?: boolean;
    is_paired?: boolean;
}

export interface ArtifactRisk {
    type: 'read_position_bias' | 'strand_bias' | 'low_mapq' | 'homopolymer' | 'low_depth';
    severity: 'low' | 'medium' | 'high';
    description: string;
    value: number;
}

export interface ArtifactAssessment {
    risks: ArtifactRisk[];
    risk_score: number;
    artifact_likelihood: 'low' | 'medium' | 'high';
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
    confidence?: 'high' | 'medium' | 'low';
    is_low_confidence?: boolean;
    alt_count?: number;
    artifact_risk?: ArtifactAssessment;
}

export interface VariantEvidence {
    forward_count: number;
    reverse_count: number;
    strand_bias: number;
    mean_quality: number;
    median_quality: number;
    // Histogram distributions for curation
    quality_histogram?: number[];
    position_histogram?: number[];
    mapq_histogram?: number[];
    // Artifact risk assessment
    artifact_risk?: ArtifactAssessment;
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
