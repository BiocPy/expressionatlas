"""BioStudies API client for Expression Atlas."""

from __future__ import annotations

import logging
from typing import Any
from urllib.parse import quote

import requests

from expression_atlas.exceptions import APIError
from expression_atlas.models import ExperimentType, SearchResult

logger = logging.getLogger(__name__)

# API endpoints
BIOSTUDIES_API_BASE = "http://www.ebi.ac.uk/biostudies/api/v1"
BIOSTUDIES_SEARCH_URL = f"{BIOSTUDIES_API_BASE}/search"
BIOSTUDIES_STUDY_URL = f"{BIOSTUDIES_API_BASE}/studies"

# Default page size for pagination
DEFAULT_PAGE_SIZE = 100


class BioStudiesAPI:
    """Client for BioStudies API to search Expression Atlas experiments."""

    def __init__(self, timeout: int = 30) -> None:
        """
        Initialize BioStudies API client.

        Parameters
        ----------
        timeout : int
            Request timeout in seconds (default: 30).
        """
        self.timeout = timeout
        self.session = requests.Session()

    def search(
        self,
        properties: list[str],
        species: str | None = None,
        page_size: int = DEFAULT_PAGE_SIZE,
    ) -> list[SearchResult]:
        """
        Search for Expression Atlas experiments.

        Parameters
        ----------
        properties : list[str]
            Search terms (e.g., ["cancer", "breast"]).
        species : str, optional
            Species to filter by (e.g., "homo sapiens").
        page_size : int
            Number of results per page (default: 100).

        Returns
        -------
        list[SearchResult]
            List of search results with experiment metadata.

        Raises
        ------
        APIError
            If the API request fails.
        """
        # Build query URL
        query_terms = "".join(quote(p) for p in properties)
        url = f"{BIOSTUDIES_SEARCH_URL}?query={query_terms}&gxa=TRUE&pageSize={page_size}"

        if species:
            url += f"&organism={quote(species)}"

        logger.info("Searching for Expression Atlas experiments...")

        # Initial request to get total count
        response = self._request(url)
        data = response.json()

        total_hits = int(data.get("totalHits", 0))
        if total_hits == 0:
            logger.info("No results found.")
            return []

        logger.info(f"Found {total_hits} experiments matching query.")

        # Verify page size
        if int(data.get("pageSize", 0)) != page_size:
            logger.warning("Page size mismatch in API response.")
            return []

        if not data.get("isTotalHitsExact", True):
            logger.warning("Total hits count from BioStudies is not exact.")

        # Paginate through all results
        all_accessions = self._paginate_results(url, total_hits, page_size)

        if len(all_accessions) != total_hits:
            logger.warning(
                f"Expected {total_hits} accessions, got {len(all_accessions)}."
            )

        # Fetch metadata for each experiment
        logger.info(f"Retrieving metadata for {len(all_accessions)} experiments...")
        results = self._fetch_experiment_metadata(all_accessions)
        logger.info("Metadata retrieval completed.")

        return results

    def _paginate_results(
        self, base_url: str, total_hits: int, page_size: int
    ) -> list[str]:
        """Paginate through search results to collect all accessions."""
        all_accessions: list[str] = []

        # Calculate number of pages
        num_pages = (total_hits + page_size - 1) // page_size

        for page_num in range(1, num_pages + 1):
            page_url = f"{base_url}&page={page_num}"
            response = self._request(page_url)
            data = response.json()

            accessions = data.get("hits", [])
            if isinstance(accessions, list) and accessions:
                # Extract accession from each hit
                for hit in accessions:
                    if isinstance(hit, dict):
                        acc = hit.get("accession")
                    else:
                        acc = hit
                    if acc:
                        all_accessions.append(acc)

        return all_accessions

    def _fetch_experiment_metadata(self, accessions: list[str]) -> list[SearchResult]:
        """Fetch detailed metadata for each experiment."""
        results: list[SearchResult] = []

        for accession in accessions:
            url = f"{BIOSTUDIES_STUDY_URL}/{accession}"

            try:
                response = self._request(url)
                data = response.json()
                result = self._parse_study_metadata(accession, data)
            except (APIError, requests.RequestException) as e:
                logger.warning(f"Failed to fetch metadata for {accession}: {e}")
                result = SearchResult(
                    accession=accession,
                    species=None,
                    experiment_type=None,
                    title=None,
                    connection_error=True,
                )

            results.append(result)

        return results

    def _parse_study_metadata(self, accession: str, data: dict[str, Any]) -> SearchResult:
        """Parse study metadata from BioStudies API response."""
        attributes = data.get("section", {}).get("attributes", [])

        title = self._extract_attribute(attributes, "Title")
        species = self._extract_attribute(attributes, "Organism")
        exp_type = self._extract_attribute(attributes, "Study type")

        # If multiple study types, find the first eligible one
        if exp_type is None:
            exp_types = self._extract_all_attributes(attributes, "Study type")
            exp_type = self._get_eligible_experiment_type(exp_types)

        return SearchResult(
            accession=accession,
            species=species,
            experiment_type=exp_type,
            title=title,
            connection_error=False,
        )

    def _extract_attribute(
        self, attributes: list[dict[str, Any]], name: str
    ) -> str | None:
        """Extract first matching attribute value."""
        for attr in attributes:
            if attr.get("name") == name:
                return attr.get("value")
        return None

    def _extract_all_attributes(
        self, attributes: list[dict[str, Any]], name: str
    ) -> list[str]:
        """Extract all matching attribute values."""
        values = []
        for attr in attributes:
            if attr.get("name") == name:
                val = attr.get("value")
                if val:
                    values.append(val)
        return values

    def _get_eligible_experiment_type(self, exp_types: list[str]) -> str | None:
        """Find the first eligible Expression Atlas experiment type."""
        eligible = ExperimentType.get_eligible_types()
        for exp_type in exp_types:
            if exp_type in eligible:
                return exp_type
        # Return first if no eligible match found
        return exp_types[0] if exp_types else None

    def _request(self, url: str) -> requests.Response:
        """Make HTTP request with error handling."""
        try:
            response = self.session.get(url, timeout=self.timeout)
        except requests.RequestException as e:
            raise APIError(0, str(e)) from e

        if response.status_code != 200:
            raise APIError(response.status_code)

        return response

    def close(self) -> None:
        """Close the session."""
        self.session.close()

    def __enter__(self) -> BioStudiesAPI:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()
