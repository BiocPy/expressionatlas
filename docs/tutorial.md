# Tutorial

The `pyexpressionatlas` package bridges the gap between Python and the [EMBL-EBI Expression Atlas](https://www.ebi.ac.uk/gxa), making it simple to search, download, and analyze curated gene expression datasets in Python. This client mirrors the [Bioconductor package](https://www.bioconductor.org/packages/release/data/experiment/html/ExpressionAtlas.html), and plugs right into the BiocPy ecosystem (`BiocFrame`, `SummarizedExperiment`, and `SingleCellExperiment`).

---

## Getting Started

First, make sure you have the package installed:

```bash
pip install pyexpressionatlas
```

Then, import the client:

```python
from pyexpressionatlas import ExpressionAtlasClient

# Initialize the client.
# By default, files are cached in ~/.cache/pyexpressionatlas_bfc.
# You can customize this directory by passing the cache_dir parameter:
client = ExpressionAtlasClient(timeout=30, cache_dir="/my/custom/cache/path")
```

## Searching for Datasets

Expression Atlas hosts thousands of curated experiments. You can search for terms related to diseases, cell lines, developmental stages, and more.

Let's say we're looking for breast cancer datasets in humans:

```python
results = client.search_experiments(
    properties=["breast", "cancer"],
    species="homo sapiens"
)

print(results)
```

This returns a `BiocFrame` containing all matches.

If you want the complete metadata for your hits, you can fetch it on-demand:

```python
# Grab the first 5 accessions from our search
accessions = results.get_column("Accession")[:5]

# Fetch complete metadata for just these experiments
metadata = client.fetch_experiment_metadata(accessions)
print(metadata)
```

## Downloading Bulk Experiments

Let's download `E-MTAB-1625` as an example. When you download a bulk experiment, it usually comes back as a `NamedList` containing `SummarizedExperiment` objects (since an accession might contain multiple array designs or assays).

```python
exp = client.get_experiment("E-MTAB-1625")

# For RNA-seq datasets, the data lives under the "rnaseq" key
rnaseq = exp["rnaseq"]
print(rnaseq)
```

This is a `SummarizedExperiment` object. You can access the count matrix, sample annotations, and gene annotations using standard accessors:

```python
# The expression matrix (numpy array: genes x samples)
counts = rnaseq.assay("counts")

# The sample annotations (BiocFrame)
sample_metadata = rnaseq.get_column_data()

# The gene annotations (BiocFrame)
gene_metadata = rnaseq.get_row_data()
```

It works exactly the same for Microarray data, except the keys in the `NamedList` will correspond to the array design (e.g., `A-AFFY-126`), and the assay is usually called `"exprs"`.

## Downloading Single-Cell Experiments

If you pass an accession from the **Single Cell Expression Atlas**, the client automatically detects it, falls back to the single-cell FTP, and grabs the Matrix Market components. It returns a `SingleCellExperiment` object.

```python
# Let's download a single-cell dataset
sc_exp = client.get_experiment("E-MTAB-6945")

print(type(sc_exp))
# <class 'singlecellexperiment.SingleCellExperiment.SingleCellExperiment'>

print(sc_exp)
```

You can interact with `sc_exp` exactly like a bulk `SummarizedExperiment`.

## Batch Downloads

If you're doing meta-analysis across dozens of experiments, grabbing them one by one is tedious. Use `get_experiments` to download a whole batch at once:

```python
# Download a batch of accessions
batch_results = client.get_experiments(["E-MTAB-1624", "E-MTAB-1625", "E-MTAB-6945"])

# Iterate through the downloaded objects
for accession, data in batch_results.items():
    print(f"Loaded {accession} successfully!")
```

And that's it! You're ready to start pulling down Expression Atlas data and throwing it straight into your Python pipelines. Happy analyzing!
