"""Configuration settings for the AWS converter service."""

from __future__ import annotations

from functools import lru_cache

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # AWS Configuration
    aws_region: str = "us-east-1"
    s3_bucket_name: str = "expression-atlas-converter"

    # Security
    converter_api_key: str | None = None  # Optional API key auth
    allowed_domains: list[str] = [
        "ftp.ebi.ac.uk",
        "www.ebi.ac.uk",
        "ebi.ac.uk",
    ]

    # Limits
    max_rdata_size_mb: int = 500
    signed_url_expiry_minutes: int = 60
    download_timeout_seconds: int = 300
    r_timeout_seconds: int = 600

    # Paths
    temp_dir: str = "/tmp"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


@lru_cache()
def get_settings() -> Settings:
    """Get cached settings instance."""
    return Settings()
