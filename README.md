# Expression Atlas Python Client

[![PyPI version](https://badge.fury.io/py/expression-atlas.svg)](https://badge.fury.io/py/expression-atlas)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Python client for searching and downloading gene expression datasets from [EMBL-EBI Expression Atlas](https://www.ebi.ac.uk/gxa), providing full compatibility with the [R Bioconductor package](https://bioconductor.org/packages/ExpressionAtlas/).

## Overview

Expression Atlas is a comprehensive resource of gene and protein expression data across species and biological conditions. This Python package provides programmatic access to:

- **Search**: Query thousands of curated RNA-seq and microarray experiments
- **Download**: Retrieve experiment data with automatic format handling
- **Analyze**: Work with R-compatible data structures in Python

### Key Features

- **Full R Parity**: Data structures mirror Bioconductor's `SummarizedExperiment` and `ExpressionSet`
- **Download Strategy**: R loading, TSV parsing, or optional cloud converter
- **Type Safety**: Comprehensive type hints for IDE autocomplete
- **Extensive Tests**: Includes R parity tests
- **BioStudies Integration**: Search across all Expression Atlas experiments

---

## Installation

### Basic Installation

```bash
pip install expression-atlas
```

### With R Support (Recommended)

For full compatibility and optimal performance:

```bash
pip install expression-atlas[r]
```

**Requirements**: R must be installed separately. See [rpy2 installation guide](https://rpy2.github.io/doc/latest/html/overview.html#installation).

### Development Installation

```bash
git clone https://github.com/gdeol4/expression-atlas.git
cd expression-atlas
pip install -e ".[dev]"
```

---

## Quick Start

### Basic Usage

```python
from expression_atlas import ExpressionAtlasClient

# Initialize client
client = ExpressionAtlasClient()

# Search for experiments
results = client.search_experiments(
    properties=["cancer", "breast"],
    species="homo sapiens"
)
print(results.head())
#   Accession        Species                Type  ...
# 0 E-MTAB-1624  homo sapiens  microarray data  ...
```

### Download RNA-seq Data

```python
# Download a single experiment
exp = client.get_experiment("E-MTAB-7041")

# Access RNA-seq data (SummarizedExperiment)
rnaseq = exp["rnaseq"]
counts = rnaseq.assays["counts"]  # numpy array: genes × samples

print(f"Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
# Shape: 58735 genes × 48 samples

# Sample metadata
sample_info = rnaseq.colData
print(sample_info.columns)
# Index(['Sample', 'Assay', 'treatment', 'time', ...])

# Gene annotations
gene_info = rnaseq.rowData
print(gene_info.head())
```

### Download Microarray Data

```python
exp = client.get_experiment("E-MTAB-1624")

# Microarray data is keyed by array design
array_design = "A-AFFY-126"
eset = exp[array_design]

# Expression matrix (probes × samples)
intensities = eset.exprs
print(intensities.shape)
# (54675, 96)

# Phenotype data (R: pData)
sample_annotations = eset.phenoData  # or eset.pData
print(sample_annotations.head())

# Feature annotations (R: fData)
probe_annotations = eset.featureData  # or eset.fData
```

### Batch Downloads

```python
# Download multiple experiments
accessions = results["Accession"].head(10).tolist()
experiments = client.get_experiments(accessions)

# Access individual experiments
for acc, exp in experiments.items():
    if exp is not None:
        print(f"{acc}: {exp['rnaseq'].shape if 'rnaseq' in exp else 'microarray'}")
```

---

## Architecture

### Download Strategy

The package employs multiple paths for maximum compatibility:

```
┌─────────────────────────────────────────────────────┐
│ 1. R + rpy2 (if available)                         │
│    → Load native .Rdata files                       │
│    → Full fidelity, exact R behavior                │
└─────────────────────────────────────────────────────┘
                    ↓ (fallback)
┌─────────────────────────────────────────────────────┐
│ 2. TSV Parser (no R required)                       │
│    → Download from EBI FTP server                   │
│    → Parse TSV files directly                       │
└─────────────────────────────────────────────────────┘
                    ↓ (fallback)
┌─────────────────────────────────────────────────────┐
│ 3. Cloud Converter Service (optional)               │
│    → Convert .Rdata via AWS Lambda/App Runner       │
│    → Returns portable MTX/CSV bundles               │
└─────────────────────────────────────────────────────┘
```

**Benefits**:
- Works without R via TSV or converter
- R path preserves exact Bioconductor structures
- Extensible to add new backends without breaking code

### Data Structures

#### `SummarizedExperiment` (RNA-seq)

Mirrors [Bioconductor's SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/):

```python
sumexp.assays["counts"]  # Primary matrix (genes × samples)
sumexp.colData           # Sample annotations (DataFrame)
sumexp.rowData           # Gene annotations (DataFrame)
sumexp.metadata          # Experiment-level metadata (dict)
sumexp.rownames          # Gene IDs (list)
sumexp.colnames          # Sample IDs (list)
```

**Matrix Orientation**: `[genes, samples]` — matches R exactly (not transposed).

#### `ExpressionSet` (Microarray)

Mirrors [Biobase's ExpressionSet](https://bioconductor.org/packages/Biobase/):

```python
eset.exprs             # Expression matrix (probes × samples)
eset.phenoData         # Sample annotations (DataFrame)
eset.featureData       # Probe/gene annotations (DataFrame)
eset.experimentData    # Preprocessing metadata (dict)

# R-style aliases
eset.pData             # Alias for phenoData
eset.fData             # Alias for featureData
```

#### `SimpleList`

Dict-like container matching [S4Vectors::SimpleList](https://bioconductor.org/packages/S4Vectors/):

```python
# RNA-seq experiments
exp["rnaseq"]  # → SummarizedExperiment

# Microarray experiments (by array design)
exp["A-AFFY-126"]  # → ExpressionSet
```

---

## Advanced Usage

### Cloud Converter Service

For environments without R, you can deploy an AWS-based converter service:

#### Service Architecture

```
┌──────────────┐      ┌───────────────┐      ┌──────────────┐
│  Python      │─────→│  AWS Lambda   │─────→│  S3 Bucket   │
│  Client      │      │  / App Runner │      │  (Cache)     │
└──────────────┘      └───────────────┘      └──────────────┘
                            │
                            └──→ R + rpy2 in container
                                 Converts .Rdata → MTX/CSV
```

#### Setup

1. **Deploy the converter** (AWS Lambda or App Runner):
   - Docker image with R + rpy2 + conversion scripts
   - REST API: `POST /convert` with `{"rdata_url": "...", "accession": "..."}`
   - Returns signed S3 URL to converted bundle

2. **Configure client**:
   ```python
   import os
   os.environ["CONVERTER_URL"] = "https://your-converter-service.region.amazonaws.com"
   
   from expression_atlas import ExpressionAtlasClient
   client = ExpressionAtlasClient()
   ```

3. **Authentication** (choose one):
   - **API Key**: `os.environ["CONVERTER_API_KEY"] = "your-key"`
   - **IAM (SigV4)**: Use AWS credentials for secure access

#### Converter Client API

```python
from expression_atlas.converter import ConverterClient

converter = ConverterClient(
    service_url="https://your-converter.example.com",
    use_iam_auth=False,  # Set True for AWS IAM
    cache_dir=Path.home() / ".atlas_cache"
)

# Convert and load
bundles = converter.convert_and_load(
    rdata_url="ftp://ftp.ebi.ac.uk/.../E-MTAB-7041-atlasExperimentSummary.Rdata",
    accession="E-MTAB-7041"
)
```

**Deployment Guide**: See `docs/CONVERTER_SETUP.md` (if you create this) or contact maintainers for AWS CloudFormation templates.

### Checking Availability

```python
from expression_atlas.download import (
    has_r_available,
    has_tsv_files,
    has_converter_available
)

# Check which backends are available
print(f"R available: {has_r_available()}")
print(f"TSV files for E-MTAB-1624: {has_tsv_files('E-MTAB-1624')}")
print(f"Converter configured: {has_converter_available()}")
```

### Error Handling

```python
from expression_atlas import (
    ExpressionAtlasClient,
    InvalidAccessionError,
    DownloadError,
    APIError
)

try:
    exp = client.get_experiment("INVALID-ACCESSION")
except InvalidAccessionError as e:
    print(f"Bad accession format: {e}")
except DownloadError as e:
    print(f"Download failed: {e}")
except APIError as e:
    print(f"API error (HTTP {e.status_code}): {e}")
```

---

## API Reference

### `ExpressionAtlasClient`

Main interface for searching and downloading experiments.

#### `search_experiments(properties, species=None)`

Search BioStudies for Expression Atlas experiments.

**Parameters**:
- `properties` (str | list[str]): Search terms (e.g., `"cancer"` or `["breast", "cancer"]`)
- `species` (str, optional): Species filter (e.g., `"homo sapiens"`, `"mus musculus"`)

**Returns**: `pandas.DataFrame` with columns `Accession`, `Species`, `Type`, `Title`

**Example**:
```python
results = client.search_experiments(
    properties=["rnaseq", "heart"],
    species="mus musculus"
)
```

#### `get_experiment(accession)`

Download a single experiment.

**Parameters**:
- `accession` (str): Valid ArrayExpress/BioStudies accession (format: `E-XXXX-####`)

**Returns**: `SimpleList` or `None` if download fails

**Raises**:
- `InvalidAccessionError`: If accession format is invalid
- `DownloadError`: If download fails after all fallbacks

#### `get_experiments(accessions, skip_invalid=True)`

Batch download multiple experiments.

**Parameters**:
- `accessions` (list[str]): List of experiment accessions
- `skip_invalid` (bool): If `True`, skip invalid accessions with warning

**Returns**: `SimpleList` mapping accession → experiment data

---

## Development

### Project Structure

```
expression-atlas/
├── src/
│   └── expression_atlas/
│       ├── client.py          # Main API (ExpressionAtlasClient)
│       ├── api.py             # BioStudies API wrapper
│       ├── download.py        # FTP + fallback logic
│       ├── rcompat.py         # R-compatible data structures
│       ├── converter.py       # Cloud converter client
│       ├── models.py          # Data models & enums
│       ├── validation.py      # Accession validation
│       └── exceptions.py      # Custom exceptions
├── tests/
│   ├── test_api.py            # BioStudies API tests
│   ├── test_r_parity.py       # R equivalence tests
│   ├── test_validation.py     # Accession validation tests
│   └── test_integration.py    # Network-dependent tests
├── .github/
│   └── workflows/
│       └── publish.yml        # PyPI publishing via OIDC
└── pyproject.toml             # Package metadata
```

### Running Tests

```bash
# Unit tests (fast, mocked)
pytest -m "not integration"

# Integration tests (requires network)
pytest -m integration

# All tests with coverage
pytest --cov=expression_atlas --cov-report=html
```

### Code Quality

```bash
# Linting
python -m ruff check src/

# Type checking
python -m mypy src/

# Format check
python -m ruff format --check src/
```

 

## Compatibility

Supports Python 3.9+ and works with or without R. When R is available (via rpy2), native .Rdata files are loaded; otherwise TSV parsing or the optional cloud converter path is used.

---

## Comparison with R Package

| Feature | R Package | Python Package | Notes |
|---------|-----------|----------------|-------|
| Search experiments | `searchAtlasExperiments()` | `search_experiments()` | Identical behavior |
| Download single | `getAtlasExperiment()` | `get_experiment()` | Same return structure |
| Download multiple | `getAtlasData()` | `get_experiments()` | Same return structure |
| Data structures | `SummarizedExperiment`, `ExpressionSet` | Same | Full parity |
| Matrix orientation | genes × samples | genes × samples | Preserved |
| Accession validation | `.isValidExperimentAccession()` | `is_valid_accession()` | Same regex |


---

## Citation

If you use this package in published research, please cite:

**Expression Atlas**:
> Moreno, P., et al. (2022). Expression Atlas update: gene and protein expression in multiple species. *Nucleic Acids Research*, 50(D1), D129-D140. [doi:10.1093/nar/gkab1030](https://doi.org/10.1093/nar/gkab1030)

**This Python Package**:
```
Expression Atlas Python Client (v0.1.0)
https://github.com/gdeol4/expression-atlas
```

---



## License

GNU General Public License v3.0 or later (GPL-3.0-or-later)

See [LICENSE](LICENSE) for full text.

---



## Links

- **PyPI**: [pypi.org/project/expression-atlas](https://pypi.org/project/expression-atlas)
- **GitHub**: [github.com/gdeol4/expression-atlas](https://github.com/gdeol4/expression-atlas)
- **Expression Atlas**: [ebi.ac.uk/gxa](https://www.ebi.ac.uk/gxa)
- **BioStudies**: [ebi.ac.uk/biostudies](https://www.ebi.ac.uk/biostudies)
- **R Package**: [bioconductor.org/packages/ExpressionAtlas](https://bioconductor.org/packages/ExpressionAtlas/)
