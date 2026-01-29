"""Configuration for BAMCP server, loaded from environment variables."""

import os
from dataclasses import dataclass, field


@dataclass
class BAMCPConfig:
    """Server configuration loaded from environment variables."""

    reference: str | None = None
    max_reads: int = 10000
    default_window: int = 500
    min_vaf: float = 0.1
    min_depth: int = 10
    min_mapq: int = 0

    @classmethod
    def from_env(cls) -> "BAMCPConfig":
        """Create config from environment variables."""
        return cls(
            reference=os.environ.get("BAMCP_REFERENCE"),
            max_reads=int(os.environ.get("BAMCP_MAX_READS", "10000")),
            default_window=int(os.environ.get("BAMCP_DEFAULT_WINDOW", "500")),
            min_vaf=float(os.environ.get("BAMCP_MIN_VAF", "0.1")),
            min_depth=int(os.environ.get("BAMCP_MIN_DEPTH", "10")),
            min_mapq=int(os.environ.get("BAMCP_MIN_MAPQ", "0")),
        )
