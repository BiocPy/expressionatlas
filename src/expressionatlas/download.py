"""
FTP download functionality for Expression Atlas experiments.

Provides compatibility with the BiocPy ecosystem:
- Uses rds2py to load .rds files if available (replaces rpy2 and .Rdata)
- Fallback: Downloads TSV files from FTP server

The data structures use biocutils, biocframe, and summarizedexperiment.
"""

from __future__ import annotations

import csv
import io
import logging
import os
import tempfile
from pathlib import Path
from typing import Any
from urllib.request import urlopen

import gzip
import numpy as np
import scipy.io
from pybiocfilecache import BiocFileCache

from biocframe import BiocFrame
from biocutils import NamedList
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

from .exceptions import DownloadError
from .validation import validate_accession

logger = logging.getLogger(__name__)

# FTP base URL for Expression Atlas experiment data
FTP_BASE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
FTP_SC_BASE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments"


def has_tsv_files(accession: str) -> bool:
    """Check if an experiment has TSV files available for download.

    Args:
        accession:
            Valid ArrayExpress/BioStudies accession (e.g., "E-MTAB-1624").

    Returns:
        True if TSV files are available, False otherwise.
    """
    validate_accession(accession)
    try:
        with urlopen(f"{FTP_BASE_URL}/{accession}/", timeout=10) as response:
            content = response.read().decode("utf-8")
            return any(
                x in content for x in ["-raw-counts.tsv", "-tpms.tsv", "-fpkms.tsv", "-normalized-expressions.tsv"]
            )
    except Exception:
        return False


def has_converter_available() -> bool:
    """Check if the cloud converter service is configured."""
    return bool(os.environ.get("CONVERTER_URL", ""))


_BFC_INSTANCE: BiocFileCache | None = None

def _get_cache() -> BiocFileCache:
    """Get or create the BiocFileCache instance for Expression Atlas downloads."""
    global _BFC_INSTANCE
    if _BFC_INSTANCE is None:
        cache_dir = Path.home() / ".cache" / "expressionatlas_bfc"
        cache_dir.mkdir(parents=True, exist_ok=True)
        _BFC_INSTANCE = BiocFileCache(cache_dir)
    return _BFC_INSTANCE

def set_cache_dir(cache_dir: str | Path) -> None:
    """Set the BiocFileCache directory globally.
    
    Args:
        cache_dir: Path to the new cache directory.
    """
    global _BFC_INSTANCE
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    _BFC_INSTANCE = BiocFileCache(cache_path)

def _get_filepath(bfc: BiocFileCache, resource: Any) -> str:
    """Extract file path from BiocFileCache resource record."""
    if hasattr(resource, "rpath"):
        rel_path = str(resource.rpath)
    elif hasattr(resource, "get"):
        rel_path = str(resource.get("rpath"))
    else:
        raise RuntimeError("Failed to resolve cache path.")
    return str(Path(bfc.config.cache_dir) / rel_path)

def _cached_download(url: str, key: str) -> str:
    """Download a URL and store it in BiocFileCache, or return cached path."""
    if url.startswith("file://"):
        from urllib.request import url2pathname
        return url2pathname(url[7:])

    bfc = _get_cache()
    try:
        existing = bfc.get(key)
        if existing:
            path = _get_filepath(bfc, existing)
            if os.path.exists(path) and os.path.getsize(path) > 0:
                logger.debug(f"Using cached file for {key}: {path}")
                return path
    except Exception:
        pass

    logger.info(f"Downloading {url} to cache...")
    resource = bfc.add(key, url, rtype="web", download=True)
    path = _get_filepath(bfc, resource)
    
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        try:
            bfc.remove(key)
        except Exception:
            pass
        raise RuntimeError(f"Download failed for {url}")
        
    return path



def get_atlas_experiment(experiment_accession: str) -> NamedList | None:
    """Download and return the data representing a single Expression Atlas experiment.

    Args:
        experiment_accession:
            Valid ArrayExpress/BioStudies accession (e.g., "E-MTAB-1624" or "E-MTAB-6945").

    Returns:
        For RNA-seq (bulk): NamedList with key "rnaseq" containing SummarizedExperiment
        For microarray (bulk): NamedList with array design accessions as keys, each containing SummarizedExperiment
        For Single-cell: SingleCellExperiment object
        Returns None if download fails.
    """
    validate_accession(experiment_accession)

    atlas_file = f"{experiment_accession}-atlasExperimentSummary.rds"
    full_url = f"{FTP_BASE_URL}/{experiment_accession}/{atlas_file}"

    logger.info(f"Downloading Expression Atlas experiment summary from:\n {full_url}")

    try:
        try:
            experiment_summary = _download_and_load_rds(full_url, experiment_accession)
        except DownloadError:
            try:
                rdata_file = f"{experiment_accession}-atlasExperimentSummary.Rdata"
                rdata_url = f"{FTP_BASE_URL}/{experiment_accession}/{rdata_file}"
                logger.info(f"RDS not available, trying direct RData download from:\n {rdata_url}")
                experiment_summary = _download_and_load_rds(rdata_url, experiment_accession)
            except DownloadError:
                logger.info("RData not available, trying TSV fallback...")

                try:
                    experiment_summary = _download_tsv_fallback(experiment_accession)
                except DownloadError:
                    logger.info("Bulk RData/TSV not available, trying single cell endpoint...")
                    try:
                        experiment_summary = _download_sc_experiment(experiment_accession)
                    except DownloadError:
                        if has_converter_available():
                            logger.info("TSV/SC not available, trying cloud converter service on RData...")
                            experiment_summary = _download_via_converter(
                                f"{FTP_BASE_URL}/{experiment_accession}/{experiment_accession}-atlasExperimentSummary.Rdata",
                                experiment_accession,
                            )
                        else:
                            raise

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


def get_atlas_data(experiment_accessions: list[str]) -> NamedList:
    """Download NamedList objects for one or more Expression Atlas experiments.

    Args:
        experiment_accessions:
            List of experiment accessions to download.

    Returns:
        Dictionary-like object mapping accession to experiment data.
    """
    from expressionatlas.validation import filter_valid_accessions

    if not experiment_accessions:
        raise ValueError("Please provide a vector of experiment accessions to download.")

    valid_accessions = filter_valid_accessions(experiment_accessions)

    if not valid_accessions:
        raise ValueError("None of the accessions passed are valid ArrayExpress/BioStudies accessions. Cannot continue.")

    results = NamedList()
    for accession in valid_accessions:
        experiment = get_atlas_experiment(accession)
        if experiment is not None:
            results[accession] = experiment

    return results


def _download_and_load_rds(url: str, accession: str) -> NamedList:
    """Download and load RDS or RData/rda file using rds2py."""
    import rds2py

    parsed_url = url.split("?")[0]
    suffix = ".rds"

    if parsed_url.lower().endswith(".rdata"):
        suffix = ".rdata"
    elif parsed_url.lower().endswith(".rda"):
        suffix = ".rda"

    try:
        path = _cached_download(url, url)
    except Exception as e:
        raise DownloadError(accession, str(e)) from e

    try:
        if suffix in [".rdata", ".rda"]:
            data = rds2py.read_rda(path)
        else:
            data = rds2py.read_rds(path)
    except Exception as e:
        raise DownloadError(accession, f"Failed to parse {suffix} file: {e}") from e

    result = NamedList()
    if isinstance(data, dict):
        for k, v in data.items():
            result[k] = v
    else:
        result["data"] = data

    return result


def _download_tsv_fallback(accession: str) -> NamedList:
    """Download experiment data from TSV files when RDS is not available."""
    result = NamedList()
    base_url = f"{FTP_BASE_URL}/{accession}"

    try:
        with urlopen(f"{base_url}/", timeout=20) as response:
            ftp_listing = response.read().decode("utf-8")
    except Exception as e:
        raise DownloadError(accession, f"FTP directory not accessible: {e}") from e

    files = [line.split()[-1] for line in ftp_listing.strip().split("\n") if line]

    sdrf_file = next((f for f in files if f.endswith(".condensed-sdrf.tsv")), f"{accession}.condensed-sdrf.tsv")
    design_df = _try_download_sdrf(f"{base_url}/{sdrf_file}")

    rnaseq_files = []
    for suffix, assay_name in [
        ("-raw-counts.tsv", "counts"),
        ("-raw-counts.tsv.undecorated", "counts"),
        ("-tpms.tsv", "tpms"),
        ("-fpkms.tsv", "fpkms"),
    ]:
        for f in files:
            if f == f"{accession}{suffix}":
                rnaseq_files.append((f, assay_name))

    if rnaseq_files:
        f, assay_name = rnaseq_files[0]
        try:
            df = _download_tsv(f"{base_url}/{f}")
            result["rnaseq"] = _create_summarized_experiment_from_tsv(df, design_df, accession, assay_name)
            logger.info(f"Downloaded RNA-seq data ({assay_name}) from TSV for {accession}")
        except Exception as e:
            logger.debug(f"Failed to download or parse {f}: {e}")

    microarray_files = []
    for f in files:
        if f.startswith(f"{accession}_") and f.endswith("-normalized-expressions.tsv"):
            design = f[len(accession) + 1 :].split("-normalized-expressions")[0]
            microarray_files.append((f, design))

    for f, design in microarray_files:
        try:
            df = _download_tsv(f"{base_url}/{f}")
            result[design] = _create_summarized_experiment_from_tsv(df, design_df, accession, "exprs")
            logger.info(f"Downloaded microarray data ({design}) from TSV for {accession}")
        except Exception as e:
            logger.debug(f"Failed to download or parse {f}: {e}")

    if len(result) == 0:
        raise DownloadError(accession, "No TSV data files found in FTP directory.")

    return result


def _download_tsv(url: str) -> dict[str, list[str]]:
    """Download and parse a TSV file from URL into a column-oriented dictionary."""
    path = _cached_download(url, url)
    logger.debug(f"Reading: {path}")
    
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        data = {h: [] for h in header}

        for row in reader:
            for i, h in enumerate(header):
                val = row[i] if i < len(row) else None
                data[h].append(val)

    return data


def _try_download_sdrf(url: str) -> dict[str, dict[str, str]] | None:
    try:
        path = _cached_download(url, url)
        logger.debug(f"Reading sample annotations: {path}")
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()

        lines = content.strip().split("\n")
        if not lines:
            return None

        first_cols = lines[0].split("\t")
        if len(first_cols) >= 6 and first_cols[1] == "":
            sample_idx, attr_idx, value_idx = 2, 4, 5
        else:
            sample_idx, attr_idx, value_idx = 1, 3, 4

        records = {}
        for line in lines:
            parts = line.split("\t")
            if len(parts) > max(sample_idx, attr_idx, value_idx):
                sample_id = parts[sample_idx]
                attr_name = parts[attr_idx]
                attr_value = parts[value_idx]

                if sample_id not in records:
                    records[sample_id] = {}

                if attr_name not in records[sample_id]:
                    records[sample_id][attr_name] = attr_value

        return records
    except Exception as e:
        logger.debug(f"Could not download sample annotations: {e}")
        return None


def _create_summarized_experiment_from_tsv(
    df_data: dict[str, list[str]],
    design_data: dict[str, dict[str, str]] | None,
    accession: str,
    assay_name: str = "counts",
) -> SummarizedExperiment:
    """Create SummarizedExperiment from TSV data."""
    if not df_data:
        return SummarizedExperiment()

    all_cols = list(df_data.keys())
    if not all_cols:
        return SummarizedExperiment()

    numeric_cols = []
    annotation_cols = []

    for col in all_cols:
        vals = df_data[col]
        is_num = True
        for v in vals:
            if v is not None and v.strip() != "" and v.strip().lower() != "na":
                try:
                    float(v)
                except ValueError:
                    is_num = False
                    break
        if is_num:
            numeric_cols.append(col)
        else:
            annotation_cols.append(col)

    if not numeric_cols:
        logger.warning("No numeric columns found in TSV")
        return SummarizedExperiment()

    gene_col = annotation_cols[0] if annotation_cols else all_cols[0]
    sample_cols = numeric_cols

    rownames = df_data[gene_col]
    colnames = sample_cols

    matrix_data = []
    for c in sample_cols:
        col_float = []
        for v in df_data[c]:
            if v is None or v.strip() == "" or v.strip().lower() == "na":
                col_float.append(np.nan)
            else:
                col_float.append(float(v))
        matrix_data.append(col_float)

    matrix = np.array(matrix_data, dtype=np.float64).T
    assays = {assay_name: matrix}

    row_data = {}
    for col in annotation_cols:
        if col != gene_col:
            row_data[col] = df_data[col]

    row_bioc = BiocFrame(row_data, row_names=rownames)

    col_data = {}
    if design_data is not None and len(design_data) > 0:
        all_attrs = set()
        for s in colnames:
            if s in design_data:
                all_attrs.update(design_data[s].keys())

        all_attrs = sorted(list(all_attrs))

        for attr in all_attrs:
            col_data[attr] = []
            for s in colnames:
                val = design_data.get(s, {}).get(attr, None)
                col_data[attr].append(val)

    col_bioc = BiocFrame(col_data, row_names=colnames)

    metadata = {"accession": accession, "source": "tsv"}

    return SummarizedExperiment(assays=assays, row_data=row_bioc, column_data=col_bioc, metadata=metadata)


def _download_via_converter(rdata_url: str, accession: str) -> NamedList:
    """Download experiment data via cloud converter."""
    from expressionatlas.converter import ConverterClient, ConverterError

    client = ConverterClient()

    try:
        bundles = client.convert_and_load(rdata_url, accession)
    except ConverterError as e:
        raise DownloadError(accession, f"Cloud converter failed: {e}") from e

    result = NamedList()

    for name, bundle in bundles.items():
        key = name.replace("dataset_", "") if name.startswith("dataset_") else name

        # We need to make sure bundle.genes and bundle.samples return dictionaries of column names mapping to lists of values
        row_bioc = BiocFrame(bundle.genes, row_names=bundle.rownames)
        col_bioc = BiocFrame(bundle.samples, row_names=bundle.colnames)

        assays = {}
        if bundle.matrix is not None:
            assays["counts" if key == "rnaseq" else "exprs"] = bundle.matrix

        meta = bundle.meta.copy()
        meta["source"] = "converter"

        se = SummarizedExperiment(assays=assays, row_data=row_bioc, column_data=col_bioc, metadata=meta)
        result[key] = se

    return result


def _download_sc_experiment(accession: str) -> SingleCellExperiment:
    """Download and construct a SingleCellExperiment from SC Expression Atlas."""
    base_url = f"{FTP_SC_BASE_URL}/{accession}"
    logger.info(f"Trying single cell FTP for {accession}: {base_url}/")

    try:
        mtx_url = f"{base_url}/{accession}.aggregated_filtered_counts.mtx.gz"
        logger.debug(f"Downloading mtx.gz: {mtx_url}")
        mtx_path = _cached_download(mtx_url, mtx_url)
            
        logger.debug("Parsing mtx...")
        with open(mtx_path, "rb") as f:
            matrix = scipy.io.mmread(io.BytesIO(gzip.decompress(f.read())))
        
        logger.debug("Downloading mtx rows and cols...")
        rows_url = f"{base_url}/{accession}.aggregated_filtered_counts.mtx_rows"
        rows_path = _cached_download(rows_url, rows_url)
        with open(rows_path, "r", encoding="utf-8") as f:
            rows = [line.split()[-1] for line in f.read().strip().split('\n')]
            
        cols_url = f"{base_url}/{accession}.aggregated_filtered_counts.mtx_cols"
        cols_path = _cached_download(cols_url, cols_url)
        with open(cols_path, "r", encoding="utf-8") as f:
            cols = [line.strip() for line in f.read().strip().split('\n')]
            
    except Exception as e:
        raise DownloadError(accession, f"Failed to download single cell MTX components: {e}") from e

    design_data = _try_download_sdrf(f"{base_url}/{accession}.condensed-sdrf.tsv")
    col_data = {}
    if design_data is not None and len(design_data) > 0:
        all_attrs = set()
        for s in cols:
            if s in design_data:
                all_attrs.update(design_data[s].keys())
        all_attrs = sorted(list(all_attrs))
        for attr in all_attrs:
            col_data[attr] = []
            for s in cols:
                val = design_data.get(s, {}).get(attr, None)
                col_data[attr].append(val)

    row_bioc = BiocFrame({}, row_names=rows)
    col_bioc = BiocFrame(col_data, row_names=cols)
    metadata = {"accession": accession, "source": "sc_mtx"}

    return SingleCellExperiment(
        assays={"counts": matrix},
        row_data=row_bioc,
        column_data=col_bioc,
        metadata=metadata,
    )


download_experiment = get_atlas_experiment
download_experiments = get_atlas_data
