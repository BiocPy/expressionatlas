"""
FTP download functionality for Expression Atlas experiments.

Provides full compatibility with the R package:
- If R is installed: Uses rpy2 to load .Rdata files (exact R behavior)
- Fallback: Downloads TSV files from FTP server

The data structures (SummarizedExperiment, ExpressionSet) match R exactly.
"""

import io
import logging
import tempfile
from pathlib import Path
from typing import Any, TYPE_CHECKING
from urllib.request import urlopen
from urllib.error import URLError

import numpy as np
import pandas as pd

from expression_atlas.exceptions import DownloadError
from expression_atlas.rcompat import ExpressionSet, SimpleList, SummarizedExperiment
from expression_atlas.validation import validate_accession

logger = logging.getLogger(__name__)

# FTP base URL for Expression Atlas experiment data
FTP_BASE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

# Lazy check for rpy2/R availability
_rpy2_checked = False
_has_rpy2 = False


def _check_rpy2() -> bool:
    """Lazily check if rpy2 and R are available."""
    global _rpy2_checked, _has_rpy2

    if _rpy2_checked:
        return _has_rpy2

    _rpy2_checked = True
    try:
        import rpy2.robjects as ro
        # Try to actually use R
        ro.r("1+1")
        _has_rpy2 = True
        logger.debug("rpy2 and R available - using native .Rdata loading")
    except Exception as e:
        _has_rpy2 = False
        logger.debug(f"rpy2/R not available ({e}) - using TSV fallback")

    return _has_rpy2


def get_atlas_experiment(experiment_accession: str) -> SimpleList | None:
    """
    Download and return the SimpleList object representing a single Expression Atlas experiment.

    This is the Python equivalent of R's getAtlasExperiment() function.

    Parameters
    ----------
    experiment_accession : str
        Valid ArrayExpress/BioStudies accession (e.g., "E-MTAB-1624").

    Returns
    -------
    SimpleList or None
        For RNA-seq: SimpleList with key "rnaseq" containing SummarizedExperiment
        For microarray: SimpleList with array design accessions as keys, each containing ExpressionSet
        Returns None if download fails.

    Data Structure
    --------------
    The returned data matches R exactly:

    RNA-seq (SummarizedExperiment):
        exp["rnaseq"].assays["counts"]  # genes × samples matrix (same as R!)
        exp["rnaseq"].colData           # sample annotations
        exp["rnaseq"].rowData           # gene annotations
        exp["rnaseq"].rownames          # gene IDs
        exp["rnaseq"].colnames          # sample IDs

    Microarray (ExpressionSet):
        exp["A-AFFY-126"].exprs         # probes × samples matrix (same as R!)
        exp["A-AFFY-126"].phenoData     # sample annotations (pData in R)
        exp["A-AFFY-126"].featureData   # probe annotations (fData in R)

    Examples
    --------
    >>> exp = get_atlas_experiment("E-MTAB-1625")
    >>> sumexp = exp["rnaseq"]
    >>> counts = sumexp.assays["counts"]  # numpy array, genes × samples
    >>> print(f"Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
    """
    validate_accession(experiment_accession)

    # Build URL (same as R package)
    atlas_file = f"{experiment_accession}-atlasExperimentSummary.Rdata"
    full_url = f"{FTP_BASE_URL}/{experiment_accession}/{atlas_file}"

    logger.info(f"Downloading Expression Atlas experiment summary from:\n {full_url}")

    try:
        if _check_rpy2():
            # Use R to load .Rdata files (full compatibility)
            experiment_summary = _download_and_load_rdata_rpy2(full_url, experiment_accession)
        else:
            # Fallback: download TSV files
            experiment_summary = _download_tsv_fallback(experiment_accession)

        if experiment_summary:
            logger.info(f"Successfully downloaded experiment summary object for {experiment_accession}")
        return experiment_summary

    except Exception as e:
        logger.warning(
            f"Error encountered while trying to download experiment summary for {experiment_accession}:\n"
            f"{e}\n"
            f"There may not currently be an Expression Atlas experiment summary available for {experiment_accession}.\n"
            f"Please try again later, check the website at http://www.ebi.ac.uk/gxa/experiments/{experiment_accession},\n"
            f"or contact us at https://www.ebi.ac.uk/about/contact/support/gxa"
        )
        return None


def get_atlas_data(experiment_accessions: list[str]) -> SimpleList:
    """
    Download SimpleList objects for one or more Expression Atlas experiments.

    This is the Python equivalent of R's getAtlasData() function.

    Parameters
    ----------
    experiment_accessions : list[str]
        List of experiment accessions to download.

    Returns
    -------
    SimpleList
        Dictionary-like object mapping accession to experiment data.
        Failed downloads are excluded from the result.

    Examples
    --------
    >>> all_exps = get_atlas_data(["E-MTAB-1624", "E-MTAB-1625"])
    >>> # Access RNA-seq data
    >>> sumexp = all_exps["E-MTAB-1625"]["rnaseq"]
    >>> # Access microarray data
    >>> eset = all_exps["E-MTAB-1624"]["A-AFFY-126"]
    """
    from expression_atlas.validation import filter_valid_accessions

    if not experiment_accessions:
        raise ValueError("Please provide a vector of experiment accessions to download.")

    valid_accessions = filter_valid_accessions(experiment_accessions)

    if not valid_accessions:
        raise ValueError(
            "None of the accessions passed are valid ArrayExpress/BioStudies accessions. Cannot continue."
        )

    results = SimpleList()
    for accession in valid_accessions:
        experiment = get_atlas_experiment(accession)
        if experiment is not None:
            results[accession] = experiment

    return results


# =============================================================================
# RPY2 IMPLEMENTATION (when R is available)
# =============================================================================

def _download_and_load_rdata_rpy2(url: str, accession: str) -> SimpleList:
    """Load .Rdata file using R via rpy2 (full R compatibility)."""
    if not _check_rpy2():
        raise RuntimeError("rpy2 not available")

    import rpy2.robjects as ro

    # Download to temp file
    with tempfile.NamedTemporaryFile(suffix=".Rdata", delete=False) as tmp:
        try:
            with urlopen(url, timeout=120) as response:
                tmp.write(response.read())
        except Exception as e:
            raise DownloadError(accession, str(e)) from e
        tmp_path = Path(tmp.name)

    try:
        # Load using R via rpy2
        r_path = str(tmp_path).replace("\\", "/")
        ro.r(f'load("{r_path}")')

        # Convert R SimpleList to Python
        return _convert_r_simplelist_rpy2(accession)

    finally:
        tmp_path.unlink()


def _convert_r_simplelist_rpy2(accession: str) -> SimpleList:
    """Convert R SimpleList to Python using rpy2."""
    import rpy2.robjects as ro

    result = SimpleList()

    # Get the names of the SimpleList
    names = list(ro.r("names(experiment_summary)"))

    for name in names:
        try:
            if name == "rnaseq":
                result[name] = _convert_summarized_experiment_rpy2(name, accession)
            elif name.startswith("A-"):
                result[name] = _convert_expression_set_rpy2(name, accession)
        except Exception as e:
            logger.warning(f"Failed to convert element '{name}': {e}")

    return result


def _convert_summarized_experiment_rpy2(name: str, accession: str) -> SummarizedExperiment:
    """Convert R SummarizedExperiment to Python using rpy2."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    se = SummarizedExperiment()

    try:
        # Get counts matrix
        counts = ro.r(f'as.matrix(assays(experiment_summary${name})$counts)')
        se.assays["counts"] = np.array(counts)

        # Get row/column names
        se.rownames = list(ro.r(f'rownames(assays(experiment_summary${name})$counts)') or [])
        se.colnames = list(ro.r(f'colnames(assays(experiment_summary${name})$counts)') or [])

        # Get colData
        try:
            coldata_r = ro.r(f'as.data.frame(colData(experiment_summary${name}))')
            with localconverter(ro.default_converter + pandas2ri.converter):
                se.colData = ro.conversion.rpy2py(coldata_r)
        except Exception:
            se.colData = pd.DataFrame(index=se.colnames)

        # Get rowData
        try:
            rowdata_r = ro.r(f'as.data.frame(rowData(experiment_summary${name}))')
            with localconverter(ro.default_converter + pandas2ri.converter):
                se.rowData = ro.conversion.rpy2py(rowdata_r)
        except Exception:
            se.rowData = pd.DataFrame(index=se.rownames)

        se.metadata["accession"] = accession

    except Exception as e:
        logger.warning(f"Error converting SummarizedExperiment: {e}")

    return se


def _convert_expression_set_rpy2(name: str, accession: str) -> ExpressionSet:
    """Convert R ExpressionSet to Python using rpy2."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    eset = ExpressionSet()

    try:
        # Get expression matrix
        exprs = ro.r(f'as.matrix(exprs(experiment_summary[["{name}"]]))')
        eset.exprs = np.array(exprs)

        # Get names
        eset.featureNames = list(ro.r(f'featureNames(experiment_summary[["{name}"]])') or [])
        eset.sampleNames = list(ro.r(f'sampleNames(experiment_summary[["{name}"]])') or [])

        # Get phenoData
        try:
            pdata_r = ro.r(f'pData(experiment_summary[["{name}"]])')
            with localconverter(ro.default_converter + pandas2ri.converter):
                eset.phenoData = ro.conversion.rpy2py(pdata_r)
        except Exception:
            eset.phenoData = pd.DataFrame(index=eset.sampleNames)

        # Get featureData
        try:
            fdata_r = ro.r(f'fData(experiment_summary[["{name}"]])')
            with localconverter(ro.default_converter + pandas2ri.converter):
                eset.featureData = ro.conversion.rpy2py(fdata_r)
        except Exception:
            eset.featureData = pd.DataFrame(index=eset.featureNames)

        # Get preprocessing info
        try:
            preproc = ro.r(f'preproc(experimentData(experiment_summary[["{name}"]]))')
            eset.experimentData["preprocessing"] = list(preproc)
        except Exception:
            pass

        eset.experimentData["accession"] = accession

    except Exception as e:
        logger.warning(f"Error converting ExpressionSet: {e}")

    return eset


# =============================================================================
# TSV FALLBACK IMPLEMENTATION (when R is not available)
# =============================================================================

def _download_tsv_fallback(accession: str) -> SimpleList:
    """
    Download experiment data from TSV files when R is not available.

    TSV files provide the same data in a portable format.
    """
    result = SimpleList()
    base_url = f"{FTP_BASE_URL}/{accession}"

    # Try RNA-seq first (raw counts)
    counts_url = f"{base_url}/{accession}-raw-counts.tsv"
    design_url = f"{base_url}/{accession}-experiment-design.tsv"

    try:
        counts_df = _download_tsv(counts_url)
        design_df = _try_download_tsv(design_url)
        result["rnaseq"] = _create_summarized_experiment_from_tsv(counts_df, design_df, accession)
        logger.info(f"Downloaded RNA-seq data from TSV for {accession}")
        return result
    except URLError:
        logger.debug(f"No raw counts TSV for {accession}")

    # Try normalized expressions (microarray or RNA-seq TPM)
    norm_url = f"{base_url}/{accession}-normalized-expressions.tsv"
    try:
        norm_df = _download_tsv(norm_url)
        design_df = _try_download_tsv(design_url)
        # Could be microarray or RNA-seq normalized
        result["normalized"] = _create_expression_set_from_tsv(norm_df, design_df, accession)
        logger.info(f"Downloaded normalized data from TSV for {accession}")
        return result
    except URLError:
        logger.debug(f"No normalized TSV for {accession}")

    raise DownloadError(
        accession,
        "No TSV data files found. Install R and rpy2 for full .Rdata support, "
        "or check if the experiment has downloadable data."
    )


def _download_tsv(url: str) -> pd.DataFrame:
    """Download and parse a TSV file from URL."""
    logger.debug(f"Downloading: {url}")
    with urlopen(url, timeout=60) as response:
        content = response.read().decode("utf-8")
    return pd.read_csv(io.StringIO(content), sep="\t")


def _try_download_tsv(url: str) -> pd.DataFrame | None:
    """Try to download TSV, return None if fails."""
    try:
        return _download_tsv(url)
    except Exception:
        return None


def _create_summarized_experiment_from_tsv(
    counts_df: pd.DataFrame,
    design_df: pd.DataFrame | None,
    accession: str,
) -> SummarizedExperiment:
    """Create SummarizedExperiment from TSV data."""
    se = SummarizedExperiment()

    if counts_df.empty:
        return se

    # First column is gene IDs
    gene_col = counts_df.columns[0]
    sample_cols = list(counts_df.columns[1:])

    se.rownames = counts_df[gene_col].tolist()
    se.colnames = sample_cols

    # Expression matrix: genes × samples (same orientation as R!)
    se.assays["counts"] = counts_df[sample_cols].values.astype(np.float64)

    # Create rowData
    se.rowData = pd.DataFrame(index=se.rownames)
    se.rowData.index.name = "gene_id"

    # Create colData from design file if available
    if design_df is not None and not design_df.empty:
        # Design file usually has sample ID as first column
        se.colData = design_df.set_index(design_df.columns[0])
    else:
        se.colData = pd.DataFrame(index=se.colnames)
    se.colData.index.name = "sample_id"

    se.metadata["accession"] = accession
    se.metadata["source"] = "tsv"

    return se


def _create_expression_set_from_tsv(
    exprs_df: pd.DataFrame,
    design_df: pd.DataFrame | None,
    accession: str,
) -> ExpressionSet:
    """Create ExpressionSet from TSV data."""
    eset = ExpressionSet()

    if exprs_df.empty:
        return eset

    # First column is probe/gene IDs
    probe_col = exprs_df.columns[0]
    sample_cols = list(exprs_df.columns[1:])

    eset.featureNames = exprs_df[probe_col].tolist()
    eset.sampleNames = sample_cols

    # Expression matrix: probes × samples (same orientation as R!)
    eset.exprs = exprs_df[sample_cols].values.astype(np.float64)

    # Create featureData
    eset.featureData = pd.DataFrame(index=eset.featureNames)
    eset.featureData.index.name = "probe_id"

    # Create phenoData from design file if available
    if design_df is not None and not design_df.empty:
        eset.phenoData = design_df.set_index(design_df.columns[0])
    else:
        eset.phenoData = pd.DataFrame(index=eset.sampleNames)
    eset.phenoData.index.name = "sample_id"

    eset.experimentData["accession"] = accession
    eset.experimentData["source"] = "tsv"

    return eset


# Backwards-compatible aliases
download_experiment = get_atlas_experiment
download_experiments = get_atlas_data
