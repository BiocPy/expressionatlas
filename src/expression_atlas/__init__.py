"""
Expression Atlas Python Client

A Python client for searching and downloading gene expression datasets
from EMBL-EBI Expression Atlas.

Full R compatibility: Data structures match the R package exactly:
- SummarizedExperiment for RNA-seq (genes × samples matrix)
- ExpressionSet for microarray (probes × samples matrix)
- SimpleList for experiment containers
"""

from expression_atlas.client import ExpressionAtlasClient
from expression_atlas.download import (
    get_atlas_data,
    get_atlas_experiment,
    has_converter_available,
    has_r_available,
    has_tsv_files,
)
from expression_atlas.exceptions import (
    APIError,
    DownloadError,
    ExpressionAtlasError,
    InvalidAccessionError,
)
from expression_atlas.models import SearchResult
from expression_atlas.rcompat import (
    ExpressionSet,
    SimpleList,
    SummarizedExperiment,
)

__version__ = "0.1.0"
__all__ = [
    # Main client
    "ExpressionAtlasClient",
    # R-compatible functions (same names as R package)
    "get_atlas_experiment",  # getAtlasExperiment()
    "get_atlas_data",  # getAtlasData()
    # Utility functions
    "has_tsv_files",  # Check if experiment has TSV files
    "has_r_available",  # Check if R is available
    "has_converter_available",  # Check if cloud converter is configured
    # R-compatible data structures
    "SimpleList",
    "SummarizedExperiment",
    "ExpressionSet",
    # Other
    "SearchResult",
    "ExpressionAtlasError",
    "InvalidAccessionError",
    "DownloadError",
    "APIError",
]
