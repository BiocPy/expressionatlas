"""Pydantic schemas for request/response validation."""

from __future__ import annotations

from datetime import datetime
from enum import Enum

from pydantic import BaseModel, Field, HttpUrl, field_validator


class OutputFormat(str, Enum):
    """Supported output formats."""

    MTX_BUNDLE = "mtx_bundle"
    TSV_BUNDLE = "tsv_bundle"


class ConvertRequest(BaseModel):
    """Request schema for /convert endpoint."""

    rdata_url: str = Field(
        ...,
        description="URL to the .RData file (ftp:// or https://)",
        max_length=2000,
    )
    accession: str = Field(
        ...,
        description="Experiment accession (e.g., E-MTAB-7841)",
        pattern=r"^E-\w{4}-\d+$",
    )
    output_format: OutputFormat = Field(
        default=OutputFormat.MTX_BUNDLE,
        description="Output format for the bundle",
    )
    assay_name: str | None = Field(
        default=None,
        description="Specific assay to extract (uses first if not specified)",
    )
    force: bool = Field(
        default=False,
        description="Force re-conversion even if cached",
    )

    @field_validator("rdata_url")
    @classmethod
    def validate_url_scheme(cls, v: str) -> str:
        """Ensure URL uses allowed schemes."""
        if not v.startswith(("https://", "ftp://")):
            raise ValueError("URL must use https:// or ftp:// scheme")
        if not v.endswith(".Rdata") and not v.endswith(".RData"):
            raise ValueError("URL must point to an .Rdata file")
        return v


class DatasetInfo(BaseModel):
    """Information about a single dataset in the bundle."""

    name: str
    class_type: str
    dimensions: tuple[int, int]
    assay_names: list[str]
    row_data_columns: list[str]
    col_data_columns: list[str]


class ConversionMeta(BaseModel):
    """Metadata about the conversion result."""

    accession: str
    source_url: str
    cache_key: str
    datasets: list[DatasetInfo]
    r_version: str | None = None
    bioconductor_version: str | None = None
    converted_at: datetime


class ConvertResponse(BaseModel):
    """Response schema for successful conversion."""

    status: str = "success"
    signed_url: str
    cache_hit: bool
    meta: ConversionMeta
    expires_at: datetime


class ErrorResponse(BaseModel):
    """Response schema for errors."""

    status: str = "error"
    error: str
    detail: str | None = None
