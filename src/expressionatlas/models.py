"""Data models for Expression Atlas."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

from biocframe import BiocFrame


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


def search_results_to_biocframe(results: list[SearchResult]) -> BiocFrame:
    """Convert list of SearchResult objects to a BiocFrame."""
    columns = ["Accession", "Species", "Type", "Title"]
    if not results:
        return BiocFrame({col: [] for col in columns}, column_names=columns)

    valid_results = [r for r in results if not r.connection_error]
    
    # Sort by Species, Type, then Accession (matching R package behavior)
    valid_results.sort(
        key=lambda r: (
            r.species if r.species is not None else "",
            r.experiment_type if r.experiment_type is not None else "",
            r.accession,
        )
    )

    data = {col: [] for col in columns}
    for r in valid_results:
        d = r.to_dict()
        for col in columns:
            data[col].append(d[col])

    return BiocFrame(data, column_names=columns)
