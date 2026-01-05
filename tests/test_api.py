"""Tests for BioStudies API client."""

import pytest
import responses

from expression_atlas.api import BIOSTUDIES_SEARCH_URL, BIOSTUDIES_STUDY_URL, BioStudiesAPI
from expression_atlas.exceptions import APIError


class TestBioStudiesAPI:
    """Tests for BioStudiesAPI class."""

    @responses.activate
    def test_search_no_results(self) -> None:
        """Should return empty list when no results found."""
        responses.add(
            responses.GET,
            BIOSTUDIES_SEARCH_URL,
            json={"totalHits": 0, "pageSize": 100, "hits": []},
            status=200,
        )

        api = BioStudiesAPI()
        results = api.search(properties=["nonexistent"])
        assert results == []

    @responses.activate
    def test_search_api_error(self) -> None:
        """Should raise APIError on HTTP error."""
        responses.add(
            responses.GET,
            BIOSTUDIES_SEARCH_URL,
            status=500,
        )

        api = BioStudiesAPI()
        with pytest.raises(APIError) as exc_info:
            api.search(properties=["cancer"])
        assert exc_info.value.status_code == 500

    @responses.activate
    def test_search_single_result(self) -> None:
        """Should return single result correctly."""
        # Mock search endpoint
        responses.add(
            responses.GET,
            BIOSTUDIES_SEARCH_URL,
            json={
                "totalHits": 1,
                "pageSize": 100,
                "isTotalHitsExact": True,
                "hits": [{"accession": "E-MTAB-1624"}],
            },
            status=200,
        )

        # Mock study details endpoint
        responses.add(
            responses.GET,
            f"{BIOSTUDIES_STUDY_URL}/E-MTAB-1624",
            json={
                "accession": "E-MTAB-1624",
                "section": {
                    "attributes": [
                        {"name": "Title", "value": "Test Experiment"},
                        {"name": "Organism", "value": "Homo sapiens"},
                        {"name": "Study type", "value": "RNA-seq of coding RNA"},
                    ]
                },
            },
            status=200,
        )

        api = BioStudiesAPI()
        results = api.search(properties=["test"])

        assert len(results) == 1
        assert results[0].accession == "E-MTAB-1624"
        assert results[0].title == "Test Experiment"
        assert results[0].species == "Homo sapiens"
        assert results[0].experiment_type == "RNA-seq of coding RNA"

    @responses.activate
    def test_search_with_species_filter(self) -> None:
        """Should include species in query URL."""
        responses.add(
            responses.GET,
            BIOSTUDIES_SEARCH_URL,
            json={"totalHits": 0, "pageSize": 100, "hits": []},
            status=200,
        )

        api = BioStudiesAPI()
        api.search(properties=["cancer"], species="homo sapiens")

        # Check that the request URL included the organism parameter
        assert "organism=homo%20sapiens" in responses.calls[0].request.url

    def test_context_manager(self) -> None:
        """Should work as context manager."""
        with BioStudiesAPI() as api:
            assert api.session is not None
        # Session should be closed after exiting context
