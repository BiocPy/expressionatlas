"""
R-compatible data structures for Expression Atlas.

These classes mirror the Bioconductor R classes:
- SummarizedExperiment: For RNA-seq data (genes × samples matrix)
- ExpressionSet: For microarray data (probes × samples matrix)

The Python implementations maintain the same data orientation and access patterns
as the R versions, ensuring full compatibility with the R package's output.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class SummarizedExperiment:
    """
    Python equivalent of Bioconductor's RangedSummarizedExperiment.

    Stores RNA-seq data in the same format as the R package:
    - assays: Dict of matrices (genes × samples), typically {"counts": matrix}
    - rowData: DataFrame of gene annotations (index = gene IDs)
    - colData: DataFrame of sample annotations (index = sample IDs)
    - metadata: Dict of experiment-level metadata

    Access patterns match R:
    - assays(sumexp)$counts  ->  sumexp.assays["counts"]
    - colData(sumexp)        ->  sumexp.colData
    - rowData(sumexp)        ->  sumexp.rowData
    - metadata(sumexp)       ->  sumexp.metadata

    Example
    -------
    >>> # R: assays(sumexp)$counts[1:5, 1:3]
    >>> sumexp.assays["counts"][:5, :3]
    >>> # R: colData(sumexp)
    >>> sumexp.colData
    >>> # R: rowData(sumexp)$gene_name
    >>> sumexp.rowData["gene_name"]
    """

    # Expression matrix: genes (rows) × samples (columns)
    # Dict allows multiple assays: {"counts": raw_counts, "tpm": tpm_values}
    assays: dict[str, np.ndarray] = field(default_factory=dict)

    # Gene annotations - DataFrame with gene IDs as index
    rowData: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Sample annotations - DataFrame with sample IDs as index
    colData: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Experiment-level metadata
    metadata: dict[str, Any] = field(default_factory=dict)

    # Row names (gene IDs) - kept separate for fast access
    rownames: list[str] = field(default_factory=list)

    # Column names (sample IDs) - kept separate for fast access
    colnames: list[str] = field(default_factory=list)

    @property
    def counts(self) -> np.ndarray | None:
        """Shortcut to access the counts matrix (most common assay)."""
        return self.assays.get("counts")

    @property
    def shape(self) -> tuple[int, int]:
        """Return (n_genes, n_samples) shape."""
        if "counts" in self.assays:
            shape = self.assays["counts"].shape
            return (int(shape[0]), int(shape[1]))
        return (len(self.rownames), len(self.colnames))

    @property
    def n_genes(self) -> int:
        """Number of genes (rows)."""
        return self.shape[0]

    @property
    def n_samples(self) -> int:
        """Number of samples (columns)."""
        return self.shape[1]

    def __repr__(self) -> str:
        assay_names = list(self.assays.keys())
        return (
            f"SummarizedExperiment with {self.n_genes} genes × {self.n_samples} samples\n"
            f"  assays: {assay_names}\n"
            f"  rowData: {list(self.rowData.columns)}\n"
            f"  colData: {list(self.colData.columns)}\n"
            f"  metadata: {list(self.metadata.keys())}"
        )


@dataclass
class ExpressionSet:
    """
    Python equivalent of Bioconductor's ExpressionSet.

    Stores microarray data in the same format as the R package:
    - exprs: Matrix of normalized intensities (probes × samples)
    - phenoData: DataFrame of sample annotations
    - featureData: DataFrame of probe/gene annotations
    - experimentData: Dict with experiment metadata (preprocessing info, etc.)

    Access patterns match R:
    - exprs(eset)           ->  eset.exprs
    - pData(eset)           ->  eset.phenoData
    - fData(eset)           ->  eset.featureData
    - preproc(experimentData(eset))  ->  eset.experimentData["preprocessing"]

    Example
    -------
    >>> # R: exprs(eset)[1:5, 1:3]
    >>> eset.exprs[:5, :3]
    >>> # R: pData(eset)$condition
    >>> eset.phenoData["condition"]
    """

    # Expression matrix: probes/genes (rows) × samples (columns)
    exprs: np.ndarray = field(default_factory=lambda: np.array([]))

    # Sample annotations (phenotype data)
    phenoData: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Feature (probe/gene) annotations
    featureData: pd.DataFrame = field(default_factory=pd.DataFrame)

    # Experiment metadata including preprocessing info
    experimentData: dict[str, Any] = field(default_factory=dict)

    # Row names (probe/gene IDs)
    featureNames: list[str] = field(default_factory=list)

    # Column names (sample IDs)
    sampleNames: list[str] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int]:
        """Return (n_features, n_samples) shape."""
        return self.exprs.shape if self.exprs.size > 0 else (0, 0)

    @property
    def n_features(self) -> int:
        """Number of features/probes (rows)."""
        return self.shape[0]

    @property
    def n_samples(self) -> int:
        """Number of samples (columns)."""
        return self.shape[1]

    # Aliases to match R accessor function names
    @property
    def pData(self) -> pd.DataFrame:
        """Alias for phenoData (matches R's pData() function)."""
        return self.phenoData

    @property
    def fData(self) -> pd.DataFrame:
        """Alias for featureData (matches R's fData() function)."""
        return self.featureData

    def __repr__(self) -> str:
        return (
            f"ExpressionSet with {self.n_features} features × {self.n_samples} samples\n"
            f"  phenoData: {list(self.phenoData.columns)}\n"
            f"  featureData: {list(self.featureData.columns)}\n"
            f"  experimentData: {list(self.experimentData.keys())}"
        )


@dataclass
class SimpleList(dict):
    """
    Python equivalent of S4Vectors::SimpleList.

    A dict-like container that holds experiment data, matching R's SimpleList behavior.
    For RNA-seq experiments, access via ["rnaseq"].
    For microarray experiments, access via array design accession (e.g., ["A-AFFY-126"]).

    Example
    -------
    >>> # R: experiment$rnaseq
    >>> experiment["rnaseq"]
    >>> # R: experiment[["A-AFFY-126"]]
    >>> experiment["A-AFFY-126"]
    """

    def __repr__(self) -> str:
        items = list(self.keys())
        return f"SimpleList with {len(items)} element(s): {items}"


def r_dataframe_to_pandas(r_df: Any) -> pd.DataFrame:
    """
    Convert R DataFrame (from pyreadr) to pandas DataFrame.

    Handles the conversion of R data structures to pandas, preserving
    row names, column names, and data types as closely as possible.
    """
    if isinstance(r_df, pd.DataFrame):
        return r_df

    # pyreadr returns pandas DataFrames directly in most cases
    if hasattr(r_df, "to_pandas"):
        result: pd.DataFrame = r_df.to_pandas()
        return result

    return pd.DataFrame(r_df)


def r_matrix_to_numpy(r_matrix: Any) -> np.ndarray:
    """
    Convert R matrix to numpy array.

    Preserves the genes × samples orientation from R.
    """
    if isinstance(r_matrix, np.ndarray):
        return r_matrix

    if isinstance(r_matrix, pd.DataFrame):
        result: np.ndarray = r_matrix.values
        return result

    arr: np.ndarray = np.asarray(r_matrix)
    return arr
