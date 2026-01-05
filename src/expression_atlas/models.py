"""Data models for Expression Atlas."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import pandas as pd


class ExperimentType(str, Enum):
    """Valid Expression Atlas experiment types."""

    TRANSCRIPTION_PROFILING_ARRAY = "transcription profiling by array"
    MICRORNA_PROFILING_ARRAY = "microRNA profiling by array"
    ANTIGEN_PROFILING = "antigen profiling"
    PROTEOMIC_PROFILING = "proteomic profiling by mass spectrometer"
    RNASEQ_CODING = "RNA-seq of coding RNA"
    RNASEQ_NONCODING = "RNA-seq of non coding RNA"
    RNASEQ_TOTAL = "RNA-seq of total RNA"
    RNASEQ_SINGLE_CELL_CODING = "RNA-seq of coding RNA from single cells"
    RNASEQ_SINGLE_CELL_NONCODING = "RNA-seq of non coding RNA from single cells"

    @classmethod
    def is_rnaseq(cls, exp_type: str) -> bool:
        """Check if experiment type is RNA-seq."""
        return "rna-seq" in exp_type.lower()

    @classmethod
    def is_microarray(cls, exp_type: str) -> bool:
        """Check if experiment type is microarray."""
        return "array" in exp_type.lower()

    @classmethod
    def get_eligible_types(cls) -> list[str]:
        """Return list of all eligible experiment type values."""
        return [e.value for e in cls]


@dataclass
class SearchResult:
    """Container for search results from BioStudies API."""

    accession: str
    species: str | None
    experiment_type: str | None
    title: str | None
    connection_error: bool = False

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "Accession": self.accession,
            "Species": self.species,
            "Type": self.experiment_type,
            "Title": self.title,
        }


@dataclass
class ExperimentSummary:
    """
    Container for downloaded Expression Atlas experiment data.

    For RNA-seq experiments, data is stored in `rnaseq` attribute.
    For microarray experiments, data is stored by array design accession.
    """

    accession: str
    experiment_type: str | None = None
    species: str | None = None
    title: str | None = None

    # RNA-seq data (AnnData object with counts)
    rnaseq: Any | None = None

    # Microarray data: dict mapping array design accession -> AnnData
    microarray: dict[str, Any] = field(default_factory=dict)

    # Raw metadata from download
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def is_rnaseq(self) -> bool:
        """Check if this is an RNA-seq experiment."""
        return self.rnaseq is not None

    @property
    def is_microarray(self) -> bool:
        """Check if this is a microarray experiment."""
        return len(self.microarray) > 0

    @property
    def array_designs(self) -> list[str]:
        """Get list of array design accessions for microarray experiments."""
        return list(self.microarray.keys())

    def __getitem__(self, key: str) -> Any:
        """
        Access experiment data by key.

        For RNA-seq: use 'rnaseq'
        For microarray: use array design accession (e.g., 'A-AFFY-126')
        """
        if key == "rnaseq":
            return self.rnaseq
        return self.microarray.get(key)

    def __repr__(self) -> str:
        data_type = "RNA-seq" if self.is_rnaseq else "Microarray"
        if self.is_microarray:
            data_type += f" ({len(self.microarray)} array design(s))"
        return f"ExperimentSummary(accession='{self.accession}', type='{data_type}')"


def search_results_to_dataframe(results: list[SearchResult]) -> pd.DataFrame:
    """Convert list of SearchResult objects to a pandas DataFrame."""
    if not results:
        return pd.DataFrame(columns=["Accession", "Species", "Type", "Title"])

    data = [r.to_dict() for r in results if not r.connection_error]
    df = pd.DataFrame(data)

    # Sort by Species, Type, then Accession (matching R package behavior)
    if not df.empty:
        df = df.sort_values(["Species", "Type", "Accession"]).reset_index(drop=True)

    return df
