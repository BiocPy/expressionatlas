"""Integration tests for Expression Atlas client.

These tests require network access and hit real APIs.
Run with: pytest -m integration
"""

import pytest

from expression_atlas import ExpressionAtlasClient
from expression_atlas.validation import is_valid_accession


@pytest.mark.integration
class TestExpressionAtlasClientIntegration:
    """Integration tests for ExpressionAtlasClient."""

    def test_search_cancer_human(self) -> None:
        """Search for cancer datasets in human should return results."""
        client = ExpressionAtlasClient()
        results = client.search_experiments(properties=["cancer"], species="homo sapiens")

        assert len(results) > 0
        assert "Accession" in results.columns
        assert "Species" in results.columns
        assert "Type" in results.columns
        assert "Title" in results.columns

        # All accessions should be valid
        for acc in results["Accession"]:
            assert is_valid_accession(acc)

    def test_search_salt_oryza(self) -> None:
        """Search for salt stress in rice should return results."""
        client = ExpressionAtlasClient()
        results = client.search_experiments(properties=["salt"], species="oryza sativa")

        assert len(results) > 0

    def test_download_single_experiment(self) -> None:
        """Download a single experiment should succeed."""
        client = ExpressionAtlasClient()
        # E-MTAB-1624 is used in the R package tests
        exp = client.get_experiment("E-MTAB-1624")

        # May return None if no download method available
        # If returns SimpleList, check it has expected keys
        if exp is not None:
            # SimpleList is dict-like, check it has data
            assert len(exp) > 0, "Expected at least one dataset in result"
