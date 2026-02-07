"""Configuration for BAMCP server, loaded from environment variables."""

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class BAMCPConfig:
    """Server configuration loaded from environment variables."""

    # Genomics settings
    reference: str | None = None
    max_reads: int = 10000
    default_window: int = 500
    min_vaf: float = 0.1
    min_depth: int = 10
    min_mapq: int = 0

    # Transport settings
    transport: str = "stdio"
    host: str = "0.0.0.0"  # noqa: S104
    port: int = 8000

    # Auth settings
    auth_enabled: bool = False
    issuer_url: str = "http://localhost:8000"
    resource_server_url: str = "http://localhost:8000"
    required_scopes: list[str] | None = None
    token_expiry: int = 3600

    # External database settings
    ncbi_api_key: str | None = None
    clinvar_enabled: bool = True
    gnomad_enabled: bool = True
    gnomad_dataset: str = "gnomad_r4"
    genome_build: str = "GRCh38"

    # Cache settings
    cache_dir: str = ""  # Set in from_env()
    cache_ttl: int = 86400  # 24 hours

    @classmethod
    def from_env(cls) -> "BAMCPConfig":
        """Create config from environment variables."""
        scopes_raw = os.environ.get("BAMCP_REQUIRED_SCOPES", "")
        scopes = [s.strip() for s in scopes_raw.split(",") if s.strip()] or None

        # Set up cache directory
        cache_dir = os.environ.get("BAMCP_CACHE_DIR") or str(Path.home() / ".cache" / "bamcp")

        # Ensure cache directory exists
        Path(cache_dir).mkdir(parents=True, exist_ok=True)

        return cls(
            reference=os.environ.get("BAMCP_REFERENCE"),
            max_reads=int(os.environ.get("BAMCP_MAX_READS", "10000")),
            default_window=int(os.environ.get("BAMCP_DEFAULT_WINDOW", "500")),
            min_vaf=float(os.environ.get("BAMCP_MIN_VAF", "0.1")),
            min_depth=int(os.environ.get("BAMCP_MIN_DEPTH", "10")),
            min_mapq=int(os.environ.get("BAMCP_MIN_MAPQ", "0")),
            transport=os.environ.get("BAMCP_TRANSPORT", "stdio"),
            host=os.environ.get("BAMCP_HOST", "0.0.0.0"),  # noqa: S104
            port=int(os.environ.get("BAMCP_PORT", "8000")),
            auth_enabled=os.environ.get("BAMCP_AUTH_ENABLED", "").lower() == "true",
            issuer_url=os.environ.get("BAMCP_ISSUER_URL", "http://localhost:8000"),
            resource_server_url=os.environ.get(
                "BAMCP_RESOURCE_SERVER_URL", "http://localhost:8000"
            ),
            required_scopes=scopes,
            token_expiry=int(os.environ.get("BAMCP_TOKEN_EXPIRY", "3600")),
            ncbi_api_key=os.environ.get("BAMCP_NCBI_API_KEY"),
            clinvar_enabled=os.environ.get("BAMCP_CLINVAR_ENABLED", "true").lower() == "true",
            gnomad_enabled=os.environ.get("BAMCP_GNOMAD_ENABLED", "true").lower() == "true",
            gnomad_dataset=os.environ.get("BAMCP_GNOMAD_DATASET", "gnomad_r4"),
            genome_build=os.environ.get("BAMCP_GENOME_BUILD", "GRCh38"),
            cache_dir=cache_dir,
            cache_ttl=int(os.environ.get("BAMCP_CACHE_TTL", "86400")),
        )
