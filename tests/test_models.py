"""Tests for data models."""


from expression_atlas.models import (
    ExperimentSummary,
    ExperimentType,
    SearchResult,
    search_results_to_dataframe,
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


class TestExperimentSummary:
    """Tests for ExperimentSummary dataclass."""

    def test_is_rnaseq_with_data(self) -> None:
        """Should return True when rnaseq data present."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        exp.rnaseq = "mock_adata"
        assert exp.is_rnaseq is True

    def test_is_rnaseq_without_data(self) -> None:
        """Should return False when no rnaseq data."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        assert exp.is_rnaseq is False

    def test_is_microarray_with_data(self) -> None:
        """Should return True when microarray data present."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        exp.microarray["A-AFFY-126"] = "mock_adata"
        assert exp.is_microarray is True

    def test_array_designs(self) -> None:
        """Should return list of array design accessions."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        exp.microarray["A-AFFY-126"] = "mock1"
        exp.microarray["A-AFFY-127"] = "mock2"
        assert set(exp.array_designs) == {"A-AFFY-126", "A-AFFY-127"}

    def test_getitem_rnaseq(self) -> None:
        """Should access rnaseq data via ['rnaseq']."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        exp.rnaseq = "mock_adata"
        assert exp["rnaseq"] == "mock_adata"

    def test_getitem_microarray(self) -> None:
        """Should access microarray data via array design accession."""
        exp = ExperimentSummary(accession="E-MTAB-1624")
        exp.microarray["A-AFFY-126"] = "mock_adata"
        assert exp["A-AFFY-126"] == "mock_adata"


class TestSearchResultsToDataframe:
    """Tests for search_results_to_dataframe function."""

    def test_empty_list(self) -> None:
        """Empty list should return empty DataFrame with correct columns."""
        df = search_results_to_dataframe([])
        assert list(df.columns) == ["Accession", "Species", "Type", "Title"]
        assert len(df) == 0

    def test_filters_connection_errors(self) -> None:
        """Should exclude results with connection errors."""
        results = [
            SearchResult("E-MTAB-1624", "Human", "RNA-seq", "Test 1"),
            SearchResult("E-MTAB-1625", None, None, None, connection_error=True),
        ]
        df = search_results_to_dataframe(results)
        assert len(df) == 1
        assert df.iloc[0]["Accession"] == "E-MTAB-1624"

    def test_sorts_by_species_type_accession(self) -> None:
        """Should sort by Species, Type, then Accession."""
        results = [
            SearchResult("E-MTAB-2", "Zebra", "RNA-seq", "Test 2"),
            SearchResult("E-MTAB-1", "Human", "Array", "Test 1"),
            SearchResult("E-MTAB-3", "Human", "RNA-seq", "Test 3"),
        ]
        df = search_results_to_dataframe(results)
        # Human Array, Human RNA-seq, Zebra RNA-seq
        assert df.iloc[0]["Accession"] == "E-MTAB-1"
        assert df.iloc[1]["Accession"] == "E-MTAB-3"
        assert df.iloc[2]["Accession"] == "E-MTAB-2"
