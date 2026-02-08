"""Configuration for BAMCP server, loaded from environment variables."""

import os
from dataclasses import dataclass
from pathlib import Path

from .constants import (
    DEFAULT_CACHE_DIR,
    DEFAULT_CACHE_TTL_SECONDS,
    DEFAULT_GENOME_BUILD,
    DEFAULT_GNOMAD_DATASET,
    DEFAULT_HOST,
    DEFAULT_ISSUER_URL,
    DEFAULT_MAX_READS,
    DEFAULT_MIN_DEPTH,
    DEFAULT_MIN_MAPQ,
    DEFAULT_MIN_VAF,
    DEFAULT_PORT,
    DEFAULT_RESOURCE_SERVER_URL,
    DEFAULT_TOKEN_EXPIRY_SECONDS,
    DEFAULT_TRANSPORT,
    DEFAULT_WINDOW_SIZE,
)


@dataclass
class BAMCPConfig:
    """Server configuration loaded from environment variables."""

    # Genomics settings
    reference: str | None = None
    max_reads: int = DEFAULT_MAX_READS
    default_window: int = DEFAULT_WINDOW_SIZE
    min_vaf: float = DEFAULT_MIN_VAF
    min_depth: int = DEFAULT_MIN_DEPTH
    min_mapq: int = DEFAULT_MIN_MAPQ

    # Transport settings
    transport: str = DEFAULT_TRANSPORT
    host: str = DEFAULT_HOST
    port: int = DEFAULT_PORT

    # Auth settings
    auth_enabled: bool = False
    issuer_url: str = DEFAULT_ISSUER_URL
    resource_server_url: str = DEFAULT_RESOURCE_SERVER_URL
    required_scopes: list[str] | None = None
    token_expiry: int = DEFAULT_TOKEN_EXPIRY_SECONDS

    # External database settings
    ncbi_api_key: str | None = None
    clinvar_enabled: bool = True
    gnomad_enabled: bool = True
    gnomad_dataset: str = DEFAULT_GNOMAD_DATASET
    genome_build: str = DEFAULT_GENOME_BUILD

    # Cache settings
    cache_dir: str = ""  # Defaults to ~/.cache/bamcp if empty
    cache_ttl: int = DEFAULT_CACHE_TTL_SECONDS

    def __post_init__(self) -> None:
        """Set default cache_dir if not provided."""
        if not self.cache_dir:
            self.cache_dir = str(DEFAULT_CACHE_DIR)

    @classmethod
    def from_env(cls) -> "BAMCPConfig":
        """Create config from environment variables."""
        env = os.environ

        scopes_raw = env.get("BAMCP_REQUIRED_SCOPES", "")
        scopes = [s.strip() for s in scopes_raw.split(",") if s.strip()] or None

        # Set up cache directory
        cache_dir = env.get("BAMCP_CACHE_DIR") or str(DEFAULT_CACHE_DIR)

        # Ensure cache directory exists
        Path(cache_dir).mkdir(parents=True, exist_ok=True)

        return cls(
            reference=os.environ.get("BAMCP_REFERENCE"),
            max_reads=int(env.get("BAMCP_MAX_READS", str(DEFAULT_MAX_READS))),
            default_window=int(env.get("BAMCP_DEFAULT_WINDOW", str(DEFAULT_WINDOW_SIZE))),
            min_vaf=float(env.get("BAMCP_MIN_VAF", str(DEFAULT_MIN_VAF))),
            min_depth=int(env.get("BAMCP_MIN_DEPTH", str(DEFAULT_MIN_DEPTH))),
            min_mapq=int(env.get("BAMCP_MIN_MAPQ", str(DEFAULT_MIN_MAPQ))),
            transport=env.get("BAMCP_TRANSPORT", DEFAULT_TRANSPORT),
            host=env.get("BAMCP_HOST", DEFAULT_HOST),
            port=int(env.get("BAMCP_PORT", str(DEFAULT_PORT))),
            auth_enabled=env.get("BAMCP_AUTH_ENABLED", "").lower() == "true",
            issuer_url=env.get("BAMCP_ISSUER_URL", DEFAULT_ISSUER_URL),
            resource_server_url=env.get("BAMCP_RESOURCE_SERVER_URL", DEFAULT_RESOURCE_SERVER_URL),
            required_scopes=scopes,
            token_expiry=int(env.get("BAMCP_TOKEN_EXPIRY", str(DEFAULT_TOKEN_EXPIRY_SECONDS))),
            ncbi_api_key=env.get("BAMCP_NCBI_API_KEY"),
            clinvar_enabled=env.get("BAMCP_CLINVAR_ENABLED", "true").lower() == "true",
            gnomad_enabled=env.get("BAMCP_GNOMAD_ENABLED", "true").lower() == "true",
            gnomad_dataset=env.get("BAMCP_GNOMAD_DATASET", DEFAULT_GNOMAD_DATASET),
            genome_build=env.get("BAMCP_GENOME_BUILD", DEFAULT_GENOME_BUILD),
            cache_dir=cache_dir,
            cache_ttl=int(env.get("BAMCP_CACHE_TTL", str(DEFAULT_CACHE_TTL_SECONDS))),
        )
