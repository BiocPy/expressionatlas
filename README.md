# Expression Atlas Python Client

A Python client for searching and downloading gene expression datasets from [EMBL-EBI Expression Atlas](https://www.ebi.ac.uk/gxa).

## Features

- Search for Expression Atlas experiments by properties and species
- Download RNA-seq and microarray experiment data
- Returns data as pandas DataFrames and AnnData objects
- Full type hints and async support

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

# Initialize client
client = ExpressionAtlasClient()

# Search for experiments
results = client.search_experiments(properties=["cancer"], species="homo sapiens")
print(results)

# Download a single experiment
experiment = client.get_experiment("E-MTAB-1624")

# Download multiple experiments
experiments = client.get_experiments(["E-MTAB-1624", "E-MTAB-1625"])
```

## Data Structures

### RNA-seq Data
RNA-seq experiments are returned as `AnnData` objects containing:
- Raw counts matrix
- Sample annotations in `.obs`
- Gene annotations in `.var`
- Experiment metadata in `.uns`

### Microarray Data
Microarray experiments are returned as dictionaries keyed by array design accession,
with each value being an `AnnData` object containing normalized intensities.

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
