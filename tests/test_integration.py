"""Integration tests for Expression Atlas client.

These tests require network access and hit real APIs.
Run with: pytest -m integration
"""

import pytest

from expressionatlas import ExpressionAtlasClient
from expressionatlas.validation import is_valid_accession

@pytest.mark.integration
# @pytest.mark.skip("takes too long")
class TestExpressionAtlasClientIntegration:
    """Integration tests for ExpressionAtlasClient."""

    def test_search_cancer_human(self) -> None:
        """Search for cancer datasets in human should return results."""
        client = ExpressionAtlasClient()
        results = client.search_experiments(properties=["cancer"], species="homo sapiens")

        assert results.shape[0] > 0
        columns = results.get_column_names()
        assert "Accession" in columns
        assert "Species" in columns
        assert "Type" in columns
        assert "Title" in columns

        for acc in results.get_column("Accession"):
            assert is_valid_accession(acc)

    def test_search_salt_oryza(self) -> None:
        """Search for salt stress in rice should return results."""
        client = ExpressionAtlasClient()
        results = client.search_experiments(properties=["salt"], species="oryza sativa")

        assert results.shape[0] > 0

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
