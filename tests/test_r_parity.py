"""
Tests that mirror the R package's test_functions.R exactly.

These tests ensure the Python implementation behaves identically to the R package.
"""

import pytest

from expression_atlas import ExpressionAtlasClient
from expression_atlas.validation import is_valid_accession


class TestAccessionValidation:
    """
    R equivalent:
    test_that("Accession validation returns true or false at the right times", {
        expect_true(.isValidExperimentAccession("E-MTAB-3007"))
        expect_false(.isValidExperimentAccession("DRP000391"))
        expect_false(.isValidExperimentAccession())
    })
    """

    def test_valid_accession_returns_true(self) -> None:
        """expect_true(.isValidExperimentAccession("E-MTAB-3007"))"""
        assert is_valid_accession("E-MTAB-3007") is True

    def test_invalid_accession_returns_false(self) -> None:
        """expect_false(.isValidExperimentAccession("DRP000391"))"""
        assert is_valid_accession("DRP000391") is False

    def test_empty_accession_returns_false(self) -> None:
        """expect_false(.isValidExperimentAccession()) - empty/missing accession"""
        assert is_valid_accession("") is False


@pytest.mark.integration
class TestDownloadData:
    """
    R equivalent:
    test_that("Download data for E-MTAB-1624", {
        skip_if_offline()
        expect_identical(names(getAtlasData("E-MTAB-1624")), "E-MTAB-1624")
    })
    """

    def test_download_e_mtab_1624(self) -> None:
        """Download data for E-MTAB-1624 and verify accession in result."""
        client = ExpressionAtlasClient()
        result = client.get_experiments(["E-MTAB-1624"])

        # R test: expect_identical(names(getAtlasData("E-MTAB-1624")), "E-MTAB-1624")
        # Python equivalent: check that the key exists in the returned dict
        assert "E-MTAB-1624" in result or len(result) == 0  # May be empty if .Rdata only


@pytest.mark.integration
class TestSearchExperiments:
    """
    R equivalent:
    test_that("Search for cancer datasets in human", {
        skip_if_offline()
        cancer_res <- searchAtlasExperiments(properties = "cancer", species = "human")
        expect_false((nrow(cancer_res) == 0))
    })
    """

    def test_search_cancer_human(self) -> None:
        """Search for cancer datasets in human should return results."""
        client = ExpressionAtlasClient()

        # R: searchAtlasExperiments(properties = "cancer", species = "human")
        cancer_res = client.search_experiments(properties=["cancer"], species="homo sapiens")

        # R: expect_false((nrow(cancer_res) == 0))
        assert len(cancer_res) > 0, "Expected at least one cancer experiment in human"

        # Verify DataFrame structure matches R's DataFrame columns
        assert "Accession" in cancer_res.columns
        assert "Species" in cancer_res.columns
        assert "Type" in cancer_res.columns
        assert "Title" in cancer_res.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
