# Expression Atlas Python Client

Python client for searching and downloading gene expression datasets from [EMBL-EBI Expression Atlas](https://www.ebi.ac.uk/gxa), mirroring the R Bioconductor package.

## Features

- Search for Expression Atlas experiments by properties and species
- Download RNA-seq and microarray experiment data
- R-compatible data structures: `SummarizedExperiment` (RNA-seq) and `ExpressionSet` (microarray) wrapped in `SimpleList`
- Sync API with full type hints

## Installation

```bash
pip install expression-atlas
```

For development:

```bash
pip install -e ".[dev]"
```

## Quick Start

```python
from expression_atlas import ExpressionAtlasClient

client = ExpressionAtlasClient()

# Search experiments (DataFrame with Accession/Species/Type/Title)
results = client.search_experiments(properties=["cancer"], species="homo sapiens")

# Download a single experiment (SimpleList)
exp = client.get_experiment("E-MTAB-1624")

# RNA-seq example
rnaseq = exp["rnaseq"]  # SummarizedExperiment
counts = rnaseq.assays["counts"]  # numpy array genes × samples
sample_annotations = rnaseq.colData

# Microarray example
eset = exp["A-AFFY-126"]  # ExpressionSet
exprs = eset.exprs  # probes × samples
pheno = eset.phenoData

# Multiple experiments
exps = client.get_experiments(["E-MTAB-1624", "E-MTAB-1625"])
```

## Data Structures

### RNA-seq Data
RNA-seq experiments are returned as `SummarizedExperiment` objects containing:
- `assays["counts"]`: genes × samples matrix (orientation matches R package)
- `colData`: sample annotations
- `rowData`: gene annotations

### Microarray Data
Microarray experiments are returned as `ExpressionSet` objects containing:
- `exprs`: probes × samples matrix (orientation matches R package)
- `phenoData`: sample annotations
- `featureData`: probe annotations

## API Reference

### `ExpressionAtlasClient`

#### `search_experiments(properties, species=None)`
Search for experiments matching given properties.

**Parameters:**
- `properties`: List of search terms (e.g., `["cancer", "breast"]`)
- `species`: Optional species filter (e.g., `"homo sapiens"`)

**Returns:** `pandas.DataFrame` with columns: Accession, Species, Type, Title

#### `get_experiment(accession)`
Download a single experiment.

**Parameters:**
- `accession`: ArrayExpress/BioStudies accession (e.g., `"E-MTAB-1624"`)

**Returns:** `ExperimentSummary` object

#### `get_experiments(accessions)`
Download multiple experiments.

**Parameters:**
- `accessions`: List of accessions

**Returns:** Dictionary mapping accessions to `ExperimentSummary` objects

## License

GPL-3.0-or-later

## Links

- [Expression Atlas](https://www.ebi.ac.uk/gxa)
- [BioStudies](https://www.ebi.ac.uk/biostudies)
- [Contact Support](https://www.ebi.ac.uk/about/contact/support/gxa)
