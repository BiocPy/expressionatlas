# Expression Atlas Python Client - AI Agent Guide

## Project Overview

Python client for EMBL-EBI Expression Atlas providing **full R package compatibility**. The defining characteristic is maintaining identical data structures and behavior to the Bioconductor R package. Not a typical REST API wrapper—downloads experiment data via FTP with dual implementation strategy.

## Critical Architecture Decisions

### Dual Implementation Strategy (R Compatibility First)
The package has two parallel download paths in [download.py](../src/expression_atlas/download.py):
1. **Primary**: `rpy2` to load native `.Rdata` files when R is available
2. **Fallback**: TSV parsing when R is unavailable

Check `_check_rpy2()` pattern—lazy evaluation with global state to avoid repeated R environment checks.

### Data Structures Mirror R Bioconductor Classes
Python classes in [rcompat.py](../src/expression_atlas/rcompat.py) replicate R Bioconductor objects exactly:
- `SummarizedExperiment`: RNA-seq data (genes × samples matrix)
- `ExpressionSet`: Microarray data (probes × samples matrix)
- `SimpleList`: Container matching S4Vectors::SimpleList

**Matrix orientation matters**: genes/probes in rows, samples in columns (matching R). Property aliases provided (`pData`/`phenoData`, `fData`/`featureData`) to match R accessor patterns.

### Three-Layer API Design
1. **User-facing**: `ExpressionAtlasClient` class ([client.py](../src/expression_atlas/client.py))
2. **Low-level API**: `BioStudiesAPI` ([api.py](../src/expression_atlas/api.py))—handles pagination, metadata parsing
3. **Download layer**: FTP operations with dual R/TSV strategy ([download.py](../src/expression_atlas/download.py))

## Key Patterns

### Accession Validation
Strict regex pattern enforced: `E-XXXX-####` (e.g., `E-MTAB-1624`). See [validation.py](../src/expression_atlas/validation.py) line 9. Always validate before FTP requests.

### Error Handling
Custom exceptions in [exceptions.py](../src/expression_atlas/exceptions.py) include helpful URLs:
- `InvalidAccessionError`: Malformed accession strings
- `APIError`: BioStudies API failures (includes status code)
- `DownloadError`: FTP/data parsing failures (includes accession and reason)

All include contact links: `https://www.ebi.ac.uk/about/contact/support/gxa`

### Logging Strategy
Verbose logging throughout—users track long downloads. Pattern: `logger.info()` for progress, `logger.warning()` for skipped items. See API search pagination in [api.py](../src/expression_atlas/api.py) lines 75-85.

### R Parity Tests
[test_r_parity.py](../tests/test_r_parity.py) mirrors R package test suite exactly. Each test includes R equivalent as docstring. Example:
```python
def test_valid_accession_returns_true(self) -> None:
    """expect_true(.isValidExperimentAccession("E-MTAB-3007"))"""
    assert is_valid_accession("E-MTAB-3007") is True
```

## Development Workflow

### Running Tests
```bash
pytest                    # Unit tests (mocked, fast)
pytest -m integration     # Integration tests (require network)
pytest --cov             # With coverage (configured in pyproject.toml)
```

Integration tests marked with `@pytest.mark.integration` and skipped by default (see [pytest.ini](../pytest.ini)).

### Type Checking & Linting
Project enforces strict typing (mypy configuration in [pyproject.toml](../pyproject.toml)):
- `disallow_untyped_defs = true`
- All functions must have return type hints
- Use `SomeType | None` not `Optional[SomeType]` (Python 3.9+ union syntax)

Ruff linting enabled with line length 100. Imports sorted with `isort` profile.

### Testing with responses Library
API tests use `responses` library to mock HTTP. Pattern in [test_api.py](../tests/test_api.py):
```python
@responses.activate
def test_search_single_result(self) -> None:
    responses.add(responses.GET, BIOSTUDIES_SEARCH_URL, json={...})
    # test code
```

## External Dependencies

- **BioStudies API**: `http://www.ebi.ac.uk/biostudies/api/v1` - search and metadata
- **Expression Atlas FTP**: `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/{accession}/` - data files
- **Optional R integration**: Requires R installation + rpy2 for `.Rdata` loading

## When Modifying Code

1. **Adding new data types**: Extend `ExperimentType` enum in [models.py](../src/expression_atlas/models.py) and update `is_rnaseq()`/`is_microarray()` class methods
2. **Changing data structures**: Verify R compatibility in [rcompat.py](../src/expression_atlas/rcompat.py)—check matrix orientation and property names
3. **API changes**: Update both search and pagination logic in [api.py](../src/expression_atlas/api.py)—they're coupled
4. **New exceptions**: Add to [exceptions.py](../src/expression_atlas/exceptions.py) and include contact URL in message

## Common Gotchas

- Don't use `Optional[T]`—use `T | None` (enforced by ruff UP rules)
- Matrix indexing: `[rows, columns]` means `[genes/probes, samples]`—never transpose
- FTP URLs must include full path: `{FTP_BASE_URL}/{accession}/{filename}`
- BioStudies pagination starts at page 1, not 0
- `SimpleList` is a dict subclass—use dict operations, not list operations
