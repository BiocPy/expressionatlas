"""
FTP download functionality for Expression Atlas experiments.

Provides compatibility with the BiocPy ecosystem:
- Uses rds2py to load .rds files if available (replaces rpy2 and .Rdata)
- Fallback: Downloads TSV files from FTP server

The data structures use biocutils, biocframe, and summarizedexperiment.
"""

from __future__ import annotations

import io
import logging
import tempfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd

from biocutils import NamedList
from biocframe import BiocFrame
from summarizedexperiment import SummarizedExperiment

from expressionatlas.exceptions import DownloadError
from expressionatlas.validation import validate_accession

logger = logging.getLogger(__name__)

# FTP base URL for Expression Atlas experiment data
FTP_BASE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"


def has_tsv_files(accession: str) -> bool:
    """
    Check if an experiment has TSV files available for download.

    Parameters
    ----------
    accession : str
        Valid ArrayExpress/BioStudies accession (e.g., "E-MTAB-1624").

    Returns
    -------
    bool
        True if TSV files are available, False otherwise.
    """
    validate_accession(accession)
    base_url = f"{FTP_BASE_URL}/{accession}"

    counts_url = f"{base_url}/{accession}-raw-counts.tsv"
    norm_url = f"{base_url}/{accession}-normalized-expressions.tsv"

    for url in [counts_url, norm_url]:
        try:
            with urlopen(url, timeout=10) as response:
                response.read(100)
                return True
        except Exception:
            continue

    return False


def has_converter_available() -> bool:
    """Check if the cloud converter service is configured."""
    import os

    return bool(os.environ.get("CONVERTER_URL", ""))


def get_atlas_experiment(experiment_accession: str) -> NamedList | None:
    """
    Download and return the data representing a single Expression Atlas experiment.

    Parameters
    ----------
    experiment_accession : str
        Valid ArrayExpress/BioStudies accession (e.g., "E-MTAB-1624").

    Returns
    -------
    NamedList or None
        For RNA-seq: NamedList with key "rnaseq" containing SummarizedExperiment
        For microarray: NamedList with array design accessions as keys, each containing SummarizedExperiment
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
            logger.info("RDS not available, trying TSV fallback...")
            try:
                experiment_summary = _download_tsv_fallback(experiment_accession)
            except DownloadError:
                if has_converter_available():
                    logger.info("TSV not available, trying cloud converter service...")
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
    """
    Download NamedList objects for one or more Expression Atlas experiments.

    Parameters
    ----------
    experiment_accessions : list[str]
        List of experiment accessions to download.

    Returns
    -------
    NamedList
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
    """Download and load RDS file using rds2py."""
    import rds2py

    with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as tmp:
        try:
            with urlopen(url, timeout=120) as response:
                tmp.write(response.read())
        except Exception as e:
            raise DownloadError(accession, str(e)) from e
        tmp_path = Path(tmp.name)

    try:
        data = rds2py.read_rds(str(tmp_path))

        result = NamedList()
        if isinstance(data, dict):
            for k, v in data.items():
                result[k] = v
        else:
            # Fallback if the top level object is not a dict
            result["data"] = data

        return result
    finally:
        tmp_path.unlink()


def _download_tsv_fallback(accession: str) -> NamedList:
    """Download experiment data from TSV files when RDS is not available."""
    result = NamedList()
    base_url = f"{FTP_BASE_URL}/{accession}"

    sdrf_url = f"{base_url}/{accession}.condensed-sdrf.tsv"
    design_df = _try_download_sdrf(sdrf_url)

    counts_url = f"{base_url}/{accession}-raw-counts.tsv"
    try:
        counts_df = _download_tsv(counts_url)
        result["rnaseq"] = _create_summarized_experiment_from_tsv(counts_df, design_df, accession, "counts")
        logger.info(f"Downloaded RNA-seq data from TSV for {accession}")
        return result
    except URLError:
        logger.debug(f"No raw counts TSV for {accession}")

    norm_url = f"{base_url}/{accession}-normalized-expressions.tsv"
    try:
        norm_df = _download_tsv(norm_url)
        # mapped to SummarizedExperiment as per instructions
        result["normalized"] = _create_summarized_experiment_from_tsv(norm_df, design_df, accession, "exprs")
        logger.info(f"Downloaded normalized data from TSV for {accession}")
        return result
    except URLError:
        logger.debug(f"No normalized TSV for {accession}")

    raise DownloadError(accession, "No TSV or RDS data files found.")


def _download_tsv(url: str) -> pd.DataFrame:
    """Download and parse a TSV file from URL."""
    logger.debug(f"Downloading: {url}")
    with urlopen(url, timeout=60) as response:
        content = response.read().decode("utf-8")
    return pd.read_csv(io.StringIO(content), sep="\t")


def _try_download_sdrf(url: str) -> pd.DataFrame | None:
    try:
        logger.debug(f"Downloading sample annotations: {url}")
        with urlopen(url, timeout=60) as response:
            content = response.read().decode("utf-8")

        lines = content.strip().split("\n")
        if not lines:
            return None

        first_cols = lines[0].split("\t")
        if len(first_cols) >= 6 and first_cols[1] == "":
            sample_idx, attr_idx, value_idx = 2, 4, 5
        else:
            sample_idx, attr_idx, value_idx = 1, 3, 4

        records: list[tuple[str, str, str]] = []
        for line in lines:
            parts = line.split("\t")
            if len(parts) > max(sample_idx, attr_idx, value_idx):
                sample_id = parts[sample_idx]
                attr_name = parts[attr_idx]
                attr_value = parts[value_idx]
                records.append((sample_id, attr_name, attr_value))

        if not records:
            return None

        df = pd.DataFrame(records, columns=["sample_id", "attribute", "value"])

        result = df.pivot_table(
            index="sample_id",
            columns="attribute",
            values="value",
            aggfunc="first",
        )

        result.columns.name = None
        result.index.name = "sample_id"

        return result
    except Exception as e:
        logger.debug(f"Could not download sample annotations: {e}")
        return None


def _create_summarized_experiment_from_tsv(
    df_data: pd.DataFrame, design_df: pd.DataFrame | None, accession: str, assay_name: str = "counts"
) -> SummarizedExperiment:
    """Create SummarizedExperiment from TSV data."""
    if df_data.empty:
        return SummarizedExperiment()

    numeric_cols = df_data.select_dtypes(include=[np.number]).columns.tolist()
    annotation_cols = [c for c in df_data.columns if c not in numeric_cols]

    if not numeric_cols:
        logger.warning("No numeric columns found in TSV")
        return SummarizedExperiment()

    gene_col = annotation_cols[0] if annotation_cols else df_data.columns[0]
    sample_cols = numeric_cols

    rownames = df_data[gene_col].tolist()
    colnames = sample_cols

    assays = {assay_name: df_data[sample_cols].values.astype(np.float64)}

    row_data = {}
    for col in annotation_cols:
        if col != gene_col:
            row_data[col] = df_data[col].values.tolist()

    row_bioc = BiocFrame(row_data, row_names=rownames)

    col_data = {}
    if design_df is not None and not design_df.empty:
        reindexed_df = design_df.reindex(colnames)
        for col in reindexed_df.columns:
            col_data[col] = reindexed_df[col].values.tolist()

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

        row_bioc = BiocFrame(bundle.genes.to_dict("list"), row_names=bundle.rownames)
        col_bioc = BiocFrame(bundle.samples.to_dict("list"), row_names=bundle.colnames)

        assays = {}
        if bundle.matrix is not None:
            assays["counts" if key == "rnaseq" else "exprs"] = bundle.matrix

        meta = bundle.meta.copy()
        meta["source"] = "converter"

        se = SummarizedExperiment(assays=assays, row_data=row_bioc, column_data=col_bioc, metadata=meta)
        result[key] = se

    return result


download_experiment = get_atlas_experiment
download_experiments = get_atlas_data
