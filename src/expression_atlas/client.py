"""Main client interface for Expression Atlas."""

from __future__ import annotations

import logging
from collections.abc import Sequence

import pandas as pd

from expression_atlas.api import BioStudiesAPI
from expression_atlas.download import get_atlas_data, get_atlas_experiment
from expression_atlas.models import search_results_to_dataframe
from expression_atlas.rcompat import SimpleList
from expression_atlas.validation import validate_accession

logger = logging.getLogger(__name__)


class ExpressionAtlasClient:
    """
    Client for searching and downloading Expression Atlas data.

    This is the main entry point for interacting with Expression Atlas.
    It provides methods equivalent to the R package's exported functions:
    - search_experiments() -> searchAtlasExperiments()
    - get_experiment() -> getAtlasExperiment()
    - get_experiments() -> getAtlasData()

    Data is returned in R-compatible formats:
    - RNA-seq: SummarizedExperiment (genes × samples matrix, same orientation as R)
    - Microarray: ExpressionSet (probes × samples matrix, same orientation as R)

    Examples
    --------
    >>> client = ExpressionAtlasClient()
    >>> # Search for experiments
    >>> results = client.search_experiments(["cancer"], species="homo sapiens")
    >>> # Download a single experiment
    >>> exp = client.get_experiment("E-MTAB-1624")
    >>> # Download multiple experiments
    >>> exps = client.get_experiments(["E-MTAB-1624", "E-MTAB-1625"])
    """

    def __init__(self, timeout: int = 30) -> None:
        """
        Initialize Expression Atlas client.

        Parameters
        ----------
        timeout : int
            Request timeout in seconds (default: 30).
        """
        self.timeout = timeout
        self._api: BioStudiesAPI | None = None

    @property
    def api(self) -> BioStudiesAPI:
        """Lazy-loaded BioStudies API client."""
        if self._api is None:
            self._api = BioStudiesAPI(timeout=self.timeout)
        return self._api

    def search_experiments(
        self,
        properties: str | Sequence[str],
        species: str | None = None,
    ) -> pd.DataFrame:
        """
        Search for Expression Atlas experiments matching given criteria.

        Equivalent to R function: searchAtlasExperiments()

        Parameters
        ----------
        properties : str or list of str
            Search terms (e.g., "cancer" or ["cancer", "breast"]).
        species : str, optional
            Species to filter by (e.g., "homo sapiens", "mus musculus").
            If not provided, searches across all species.

        Returns
        -------
        pandas.DataFrame
            DataFrame with columns: Accession, Species, Type, Title.
            Sorted by Species, Type, then Accession.

        Raises
        ------
        ValueError
            If no search properties provided.
        APIError
            If the BioStudies API request fails.

        Examples
        --------
        >>> client = ExpressionAtlasClient()
        >>> # Search for salt stress experiments in rice
        >>> results = client.search_experiments("salt", species="oryza sativa")
        >>> # Search with multiple terms
        >>> results = client.search_experiments(["cancer", "breast"], species="homo sapiens")
        """
        if not properties:
            raise ValueError("Please provide at least one search term.")

        # Normalize to list
        if isinstance(properties, str):
            properties = [properties]

        if species is None:
            logger.info("No species provided. Searching across all available species.")

        results = self.api.search(properties=list(properties), species=species)

        # Filter out connection errors and convert to DataFrame
        df = search_results_to_dataframe(results)

        # Log warning if any connection errors occurred
        error_count = sum(1 for r in results if r.connection_error)
        if error_count > 0:
            logger.warning(
                f"{error_count} experiment(s) excluded due to connection errors."
            )

        return df

    def get_experiment(self, accession: str) -> SimpleList | None:
        """
        Download a single Expression Atlas experiment.

        Equivalent to R function: getAtlasExperiment()

        Parameters
        ----------
        accession : str
            ArrayExpress/BioStudies experiment accession (e.g., "E-MTAB-1624").

        Returns
        -------
        SimpleList or None
            The downloaded experiment data, or None if download fails.
            For RNA-seq: access via ["rnaseq"] to get SummarizedExperiment
            For microarray: access via array design (e.g., ["A-AFFY-126"]) to get ExpressionSet

        Raises
        ------
        InvalidAccessionError
            If the accession format is invalid.

        Examples
        --------
        >>> client = ExpressionAtlasClient()
        >>> # RNA-seq experiment
        >>> exp = client.get_experiment("E-MTAB-1625")
        >>> sumexp = exp["rnaseq"]  # SummarizedExperiment
        >>> sumexp.assays["counts"]  # counts matrix (genes × samples)
        >>> sumexp.colData  # sample annotations
        >>>
        >>> # Microarray experiment
        >>> exp = client.get_experiment("E-MTAB-1624")
        >>> eset = exp["A-AFFY-126"]  # ExpressionSet
        >>> eset.exprs  # expression matrix (probes × samples)
        >>> eset.phenoData  # sample annotations
        """
        validate_accession(accession)
        return get_atlas_experiment(accession)

    def get_experiments(
        self,
        accessions: Sequence[str],
        skip_invalid: bool = True,
    ) -> SimpleList:
        """
        Download multiple Expression Atlas experiments.

        Equivalent to R function: getAtlasData()

        Parameters
        ----------
        accessions : list of str
            List of experiment accessions to download.
        skip_invalid : bool
            If True (default), skip invalid accessions with a warning.
            If False, raise an error on invalid accessions.

        Returns
        -------
        SimpleList
            Dictionary-like object mapping accession to experiment data (SimpleList).
            Failed downloads are excluded from the result.

        Raises
        ------
        ValueError
            If no valid accessions provided.
        InvalidAccessionError
            If skip_invalid is False and an invalid accession is found.

        Examples
        --------
        >>> client = ExpressionAtlasClient()
        >>> results = client.search_experiments("cancer", species="homo sapiens")
        >>> # Download all RNA-seq experiments from search results
        >>> rnaseq_accessions = results[results["Type"].str.contains("RNA-seq", na=False)]["Accession"]
        >>> experiments = client.get_experiments(rnaseq_accessions.tolist())
        >>> # Access: experiments["E-MTAB-XXXX"]["rnaseq"].assays["counts"]
        """
        return get_atlas_data(list(accessions))

    def close(self) -> None:
        """Close the client and release resources."""
        if self._api is not None:
            self._api.close()
            self._api = None

    def __enter__(self) -> ExpressionAtlasClient:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def __repr__(self) -> str:
        return f"ExpressionAtlasClient(timeout={self.timeout})"
