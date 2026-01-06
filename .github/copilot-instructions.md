# Expression Atlas Python Client - AI Agent Guide

## Project Overview

Python client for EMBL-EBI Expression Atlas with **full R Bioconductor package compatibility**. Not a typical REST API wrapper—downloads experiment data via FTP with a tri-fallback strategy. The defining constraint: data structures and behavior must match the R package exactly.

## Architecture

### Tri-Fallback Download Strategy (`download.py`)
```
1. rpy2 + R → Load native .Rdata files (full fidelity)
2. TSV parsing → Download TSV files from FTP (no R needed)
3. Cloud converter → CONVERTER_URL env var service (no R, no TSV)
```
The `_check_rpy2()` function uses lazy evaluation with global state (`_rpy2_checked`, `_has_rpy2`) to avoid repeated R environment checks.

### Three-Layer Design
| Layer | Module | Purpose |
|-------|--------|---------|
| User API | `client.py` | `ExpressionAtlasClient` class |
| Low-level | `api.py` | BioStudies API, pagination, metadata |
| Download | `download.py` | FTP ops, R/TSV/converter fallbacks |

### R-Compatible Data Structures (`rcompat.py`)
- `SummarizedExperiment`: RNA-seq (genes × samples matrix)
- `ExpressionSet`: Microarray (probes × samples matrix)  
- `SimpleList`: Dict subclass matching S4Vectors::SimpleList

**Critical**: Matrix orientation is `[genes/probes, samples]`—never transpose. Property aliases (`pData`/`phenoData`) match R accessor patterns.

## Code Conventions

### Type Hints (Strict)
```python
# CORRECT - Python 3.9+ union syntax
def func(value: str | None) -> dict[str, Any]: ...

# WRONG - don't use Optional
def func(value: Optional[str]) -> Dict[str, Any]: ...
```

### Accession Validation
Pattern: `E-XXXX-####` (e.g., `E-MTAB-1624`). Always call `validate_accession()` before FTP requests.

### Exception Pattern
All exceptions in `exceptions.py` include EBI support URL:
```python
raise DownloadError(accession, f"Cloud converter failed: {e}")
# Message includes: "contact https://www.ebi.ac.uk/about/contact/support/gxa"
```

### Logging
Use `logger.info()` for progress, `logger.warning()` for skipped items. Users rely on logs to track long downloads.

## Development Commands

```bash
pytest                     # Unit tests (mocked, fast)
pytest -m integration      # Network tests (BioStudies API, FTP)
pytest --cov               # Coverage report
python -m ruff check src/  # Linting (line-length=100)
python -m mypy src/        # Type checking (disallow_untyped_defs=true)
```

## Testing Patterns

### Mock HTTP with `responses`
```python
@responses.activate
def test_api_call(self) -> None:
    responses.add(responses.GET, BIOSTUDIES_SEARCH_URL, json={...})
```

### R Parity Tests (`test_r_parity.py`)
Tests mirror R package exactly—include R equivalent in docstring:
```python
def test_valid_accession_returns_true(self) -> None:
    """expect_true(.isValidExperimentAccession("E-MTAB-3007"))"""
    assert is_valid_accession("E-MTAB-3007") is True
```

## External Services

| Service | URL Pattern | Purpose |
|---------|-------------|---------|
| BioStudies API | `http://www.ebi.ac.uk/biostudies/api/v1` | Search, metadata |
| Atlas FTP | `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/{accession}/` | Data files |
| Converter | `CONVERTER_URL` env var | .Rdata → portable format |

## Modification Checklist

- **New experiment type**: Add to `ExperimentType` enum in `models.py`, update `is_rnaseq()`/`is_microarray()`
- **Data structure changes**: Verify R compatibility in `rcompat.py`—check matrix orientation
- **New exception**: Add to `exceptions.py` with EBI support URL
- **API changes**: Update both search and pagination in `api.py` (they're coupled)

## Gotchas

- BioStudies pagination starts at page 1, not 0
- `SimpleList` is a dict subclass—use dict operations
- FTP URLs need full path: `{FTP_BASE_URL}/{accession}/{filename}`
- ruff ignores `N802`/`N815` for R-compatible naming (`pData`, `fData`)
