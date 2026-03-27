"""
Expression Atlas Python Client

A Python client for searching and downloading gene expression datasets
from EMBL-EBI Expression Atlas.

Full BiocPy compatibility: Data structures use BiocPy ecosystem:
- SummarizedExperiment for RNA-seq and microarray (genes × samples matrix)
- NamedList for experiment containers
"""

import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from expressionatlas.client import ExpressionAtlasClient
from expressionatlas.download import (
    get_atlas_data,
    get_atlas_experiment,
    has_converter_available,
    has_tsv_files,
)
from expressionatlas.exceptions import (
    APIError,
    DownloadError,
    ExpressionAtlasError,
    InvalidAccessionError,
)
from expressionatlas.models import SearchResult
