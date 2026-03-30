"""Tests for data models."""


from expressionatlas.models import (
    ExperimentType,
    SearchResult,
    search_results_to_biocframe,
)


class TestExperimentType:
    """Tests for ExperimentType enum."""

    def test_is_rnaseq_true(self) -> None:
        """RNA-seq types should return True."""
        assert ExperimentType.is_rnaseq("RNA-seq of coding RNA") is True
        assert ExperimentType.is_rnaseq("RNA-seq of total RNA") is True

    def test_is_rnaseq_false(self) -> None:
        """Non-RNA-seq types should return False."""
        assert ExperimentType.is_rnaseq("transcription profiling by array") is False

    def test_is_microarray_true(self) -> None:
        """Microarray types should return True."""
        assert ExperimentType.is_microarray("transcription profiling by array") is True
        assert ExperimentType.is_microarray("microRNA profiling by array") is True

    def test_is_microarray_false(self) -> None:
        """Non-microarray types should return False."""
        assert ExperimentType.is_microarray("RNA-seq of coding RNA") is False

    def test_get_eligible_types(self) -> None:
        """Should return all eligible experiment types."""
        types = ExperimentType.get_eligible_types()
        assert "RNA-seq of coding RNA" in types
        assert "transcription profiling by array" in types
        assert len(types) == 9


class TestSearchResult:
    """Tests for SearchResult dataclass."""

    def test_to_dict(self) -> None:
        """Should convert to dictionary with expected keys."""
        result = SearchResult(
            accession="E-MTAB-1624",
            species="Homo sapiens",
            experiment_type="RNA-seq of coding RNA",
            title="Test experiment",
        )
        d = result.to_dict()
        assert d["Accession"] == "E-MTAB-1624"
        assert d["Species"] == "Homo sapiens"
        assert d["Type"] == "RNA-seq of coding RNA"
        assert d["Title"] == "Test experiment"

    def test_default_connection_error(self) -> None:
        """Default connection_error should be False."""
        result = SearchResult(
            accession="E-MTAB-1624",
            species=None,
            experiment_type=None,
            title=None,
        )
        assert result.connection_error is False


class TestSearchResultsToBiocframe:
    """Tests for search_results_to_biocframe function."""

    def test_empty_list(self) -> None:
        """Empty list should return empty BiocFrame with correct columns."""
        bf = search_results_to_biocframe([])
        assert list(bf.get_column_names()) == ["Accession", "Species", "Type", "Title"]
        assert bf.shape[0] == 0

    def test_filters_connection_errors(self) -> None:
        """Should exclude results with connection errors."""
        results = [
            SearchResult("E-MTAB-1624", "Human", "RNA-seq", "Test 1"),
            SearchResult("E-MTAB-1625", None, None, None, connection_error=True),
        ]
        bf = search_results_to_biocframe(results)
        assert bf.shape[0] == 1
        assert bf.get_column("Accession")[0] == "E-MTAB-1624"

    def test_sorts_by_species_type_accession(self) -> None:
        """Should sort by Species, Type, then Accession."""
        results = [
            SearchResult("E-MTAB-2", "Zebra", "RNA-seq", "Test 2"),
            SearchResult("E-MTAB-1", "Human", "Array", "Test 1"),
            SearchResult("E-MTAB-3", "Human", "RNA-seq", "Test 3"),
        ]
        bf = search_results_to_biocframe(results)
        # Human Array, Human RNA-seq, Zebra RNA-seq
        ids = bf.get_column("Accession")
        assert ids[0] == "E-MTAB-1"
        assert ids[1] == "E-MTAB-3"
        assert ids[2] == "E-MTAB-2"
