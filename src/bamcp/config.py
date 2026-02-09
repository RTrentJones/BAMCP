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

    # Security settings
    allowed_directories: list[str] | None = None
    allow_remote_files: bool = False

    def __post_init__(self) -> None:
        """Validate config values and set defaults."""
        # Set default cache_dir if not provided
        if not self.cache_dir:
            self.cache_dir = str(DEFAULT_CACHE_DIR)

        # Validate genomics settings
        if not 0 <= self.min_vaf <= 1:
            raise ValueError(f"min_vaf must be between 0 and 1, got {self.min_vaf}")

        if self.min_depth < 1:
            raise ValueError(f"min_depth must be at least 1, got {self.min_depth}")

        if self.max_reads < 1:
            raise ValueError(f"max_reads must be at least 1, got {self.max_reads}")

        if not 0 <= self.min_mapq <= 255:
            raise ValueError(f"min_mapq must be between 0 and 255, got {self.min_mapq}")

        if self.default_window < 1:
            raise ValueError(f"default_window must be at least 1, got {self.default_window}")

        # Validate transport settings
        valid_transports = ("stdio", "sse", "streamable-http")
        if self.transport not in valid_transports:
            raise ValueError(f"transport must be one of {valid_transports}, got '{self.transport}'")

        if not 1 <= self.port <= 65535:
            raise ValueError(f"port must be between 1 and 65535, got {self.port}")

        # Validate auth settings
        if self.token_expiry < 1:
            raise ValueError(f"token_expiry must be at least 1 second, got {self.token_expiry}")

        # Validate cache settings
        if self.cache_ttl < 0:
            raise ValueError(f"cache_ttl must be non-negative, got {self.cache_ttl}")

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
            allowed_directories=[
                d.strip() for d in env.get("BAMCP_ALLOWED_DIRECTORIES", "").split(",") if d.strip()
            ]
            or None,
            allow_remote_files=env.get("BAMCP_ALLOW_REMOTE_FILES", "false").lower() == "true",
        )
