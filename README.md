[![PyPI-Server](https://img.shields.io/pypi/v/expressionatlas.svg)](https://pypi.org/project/expressionatlas/)
![Unit tests](https://github.com/biocpy/expressionatlas/actions/workflows/run-tests.yml/badge.svg)

# expressionatlas

A Python client for searching and downloading gene expression datasets from [EMBL-EBI Expression Atlas](https://www.ebi.ac.uk/gxa), providing full compatibility with the [R Bioconductor package](https://bioconductor.org/packages/ExpressionAtlas/).

> [!NOTE]
> This package is a fork of [expression-atlas](https://github.com/gdeol4/expression-atlas) to use BiocPy data structures. 

## Install

To get started, install the package from [PyPI](https://pypi.org/project/expressionatlas/)

```bash
pip install expressionatlas
```

## Get Started

Expression Atlas is a comprehensive resource of gene and protein expression data across species and biological conditions. This Python package provides programmatic access to:

- **Search**: Query thousands of curated RNA-seq and microarray experiments
- **Download**: Retrieve experiment data with automatic format handling
- **Analyze**: Work with R-compatible data structures in Python

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
counts = rnaseq.assay("counts")  # numpy array: genes × samples

print(f"Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
# Shape: 58735 genes × 48 samples

# Sample metadata (BiocFrame)
sample_info = rnaseq.get_column_data()
print(sample_info.get_column_names())

# Gene annotations (BiocFrame)
gene_info = rnaseq.get_row_data()
print(gene_info.shape)
```

### Download Microarray Data

```python
exp = client.get_experiment("E-MTAB-1624")

# Microarray data is keyed by array design
array_design = "A-AFFY-126"
eset = exp[array_design] # This is also a SummarizedExperiment now

# Expression matrix (probes × samples)
intensities = eset.assay("exprs")
print(intensities.shape)
# (54675, 96)

# Sample metadata (BiocFrame)
sample_annotations = eset.get_column_data()
print(sample_annotations.shape)

# Feature annotations (BiocFrame)
probe_annotations = eset.get_row_data()
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

<!-- biocsetup-notes -->

## Note

This project has been set up using [BiocSetup](https://github.com/biocpy/biocsetup)
and [PyScaffold](https://pyscaffold.org/).
