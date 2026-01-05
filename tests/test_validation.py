"""Tests for validation utilities."""

import pytest

from expression_atlas.exceptions import InvalidAccessionError
from expression_atlas.validation import (
    filter_valid_accessions,
    is_valid_accession,
    validate_accession,
)


class TestIsValidAccession:
    """Tests for is_valid_accession function."""

    def test_valid_mtab_accession(self) -> None:
        """Valid E-MTAB accession should return True."""
        assert is_valid_accession("E-MTAB-3007") is True

    def test_valid_geod_accession(self) -> None:
        """Valid E-GEOD accession should return True."""
        assert is_valid_accession("E-GEOD-11175") is True

    def test_valid_mexp_accession(self) -> None:
        """Valid E-MEXP accession should return True."""
        assert is_valid_accession("E-MEXP-1234") is True

    def test_invalid_drp_accession(self) -> None:
        """Non-ArrayExpress accession should return False."""
        assert is_valid_accession("DRP000391") is False

    def test_invalid_format(self) -> None:
        """Malformed accession should return False."""
        assert is_valid_accession("invalid") is False

    def test_empty_string(self) -> None:
        """Empty string should return False."""
        assert is_valid_accession("") is False

    def test_missing_prefix(self) -> None:
        """Accession without E- prefix should return False."""
        assert is_valid_accession("MTAB-1234") is False

    def test_wrong_separator(self) -> None:
        """Accession with wrong separator should return False."""
        assert is_valid_accession("E_MTAB_1234") is False


class TestValidateAccession:
    """Tests for validate_accession function."""

    def test_valid_accession_returns_same(self) -> None:
        """Valid accession should be returned unchanged."""
        assert validate_accession("E-MTAB-1624") == "E-MTAB-1624"

    def test_invalid_accession_raises(self) -> None:
        """Invalid accession should raise InvalidAccessionError."""
        with pytest.raises(InvalidAccessionError) as exc_info:
            validate_accession("invalid")
        assert "invalid" in str(exc_info.value)


class TestFilterValidAccessions:
    """Tests for filter_valid_accessions function."""

    def test_filters_invalid_accessions(self) -> None:
        """Should filter out invalid accessions."""
        accessions = ["E-MTAB-1624", "invalid", "E-GEOD-11175", "DRP000391"]
        result = filter_valid_accessions(accessions)
        assert result == ["E-MTAB-1624", "E-GEOD-11175"]

    def test_empty_list(self) -> None:
        """Empty list should return empty list."""
        assert filter_valid_accessions([]) == []

    def test_all_invalid(self) -> None:
        """All invalid accessions should return empty list."""
        result = filter_valid_accessions(["invalid", "DRP000391"])
        assert result == []

    def test_raise_on_invalid(self) -> None:
        """Should raise on invalid when raise_on_invalid=True."""
        with pytest.raises(InvalidAccessionError):
            filter_valid_accessions(["E-MTAB-1624", "invalid"], raise_on_invalid=True)
