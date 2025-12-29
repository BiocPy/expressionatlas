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
from expression_atlas.download import get_atlas_experiment, get_atlas_data
from expression_atlas.rcompat import (
    SimpleList,
    SummarizedExperiment,
    ExpressionSet,
)
from expression_atlas.models import SearchResult
from expression_atlas.exceptions import (
    ExpressionAtlasError,
    InvalidAccessionError,
    DownloadError,
    APIError,
)

__version__ = "0.1.0"
__all__ = [
    # Main client
    "ExpressionAtlasClient",
    # R-compatible functions (same names as R package)
    "get_atlas_experiment",  # getAtlasExperiment()
    "get_atlas_data",  # getAtlasData()
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
