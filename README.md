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
from expressionatlas import ExpressionAtlasClient

# Initialize client (optionally specify a custom cache directory)
client = ExpressionAtlasClient(cache_dir="~/.cache/my_custom_cache")

# Search for experiments
results = client.search_experiments(
    properties=["cancer", "breast"],
    species="homo sapiens"
)
print(results)
```
    BiocFrame with 208 rows and 4 columns
            Accession      Species                    Type                   Title
                <list>       <list>                  <list>                  <list>
    [0]  E-MTAB-8198          None                    None Functional effect of...
    [1]  E-MTAB-8532          None                    None DNA microarray studi...
    [2] E-GEOD-43306          None                    None Translating transcri...
                ...          ...                     ...                     ...
    [205]   E-MTAB-779        None                    None OncomiRs like let-7 ...
    [206]  E-TABM-1118        None                    None Transcrption profili...
    [207]   E-TABM-601        None                    None Transcription profil...

### Fetch Full Metadata

The initial search is optimized for speed and does not fetch full metadata. To retrieve complete details (including `Species` and `Type`), use `fetch_experiment_metadata`:

```python
# Fetch full metadata for specific experiments
metadata = client.fetch_experiment_metadata(["E-MTAB-8198", "E-MTAB-8532"])
print(metadata)
```
    BiocFrame with 2 rows and 4 columns
            Accession      Species                    Type                   Title
               <list>       <list>                  <list>                  <list>
    [0]   E-MTAB-8198 Homo sapiens Cell line - High-thr... Functional effect of...
    [1]   E-MTAB-8532 Homo sapiens Human - One-color mi... DNA microarray studi...

### Download RNA-seq Data

```python
# Download a single experiment
exp = client.get_experiment("E-MTAB-1625")

# Access RNA-seq data (SummarizedExperiment)
rnaseq = exp["rnaseq"]
counts = rnaseq.assay("counts")  # numpy array: genes × samples

print(f"Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
# Shape: 58735 genes × 24 samples

# Sample metadata (BiocFrame)
sample_info = rnaseq.get_column_data()
print(sample_info.get_column_names())
# ['cell line', 'compound', 'developmental stage', 'disease', 'dose', 'genotype', 'organism', 'organism part']

# Gene annotations (BiocFrame)
gene_info = rnaseq.get_row_data()
print(gene_info.shape)
# (58735, 1)

print(rnaseq)
```
    class: SummarizedExperiment
    dimensions: (58735, 24)
    assays(1): ['counts']
    row_data columns(1): ['Gene Name']
    row_names(58735): ['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419', ..., 'ENSG00000285992', 'ENSG00000285993', 'ENSG00000285994']
    column_data columns(8): ['cell line', 'compound', 'developmental stage', 'disease', 'dose', 'genotype', 'organism', 'organism part']
    column_names(24): ['ERR3456453', 'ERR3456442', 'ERR3456443', ..., 'ERR3456450', 'ERR3456459', 'ERR3456444']
    metadata(2): accession source


### Batch Downloads

```python
# Download multiple experiments
accessions = results.get_column("Accession")[:10]
experiments = client.get_experiments(accessions)

# Access individual experiments
for acc, exp in experiments.items():
    if exp is not None:
        print(f"{acc}: {exp['rnaseq'].shape if 'rnaseq' in exp else 'microarray'}")
```

### Caching Mechanism

To optimize performance and reduce load on the FTP servers, all data downloads are automatically cached locally using `pyBiocFileCache`. 

- By default, the cache is stored at `~/.cache/expressionatlas_bfc`.
- You can customize this location when initializing the client by passing the `cache_dir` argument: `client = ExpressionAtlasClient(cache_dir="/path/to/custom/cache")`.

### Direct RData / rda Support

The client automatically downloads and parses both `.rds` and `.Rdata` / `.rda` files directly without relying on a cloud converter service:
- Tries downloading and parsing `.rds` file using `rds2py.read_rds`.
- If the `.rds` file is not available, falls back to direct download and loading of the `.Rdata` file using `rds2py.read_rda`.

<!-- biocsetup-notes -->

## Note

This project has been set up using [BiocSetup](https://github.com/biocpy/biocsetup)
and [PyScaffold](https://pyscaffold.org/).
