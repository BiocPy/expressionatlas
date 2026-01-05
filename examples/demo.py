#!/usr/bin/env python
"""
Expression Atlas Python Client - Demo Script

This script demonstrates the main features of the expression-atlas package:
1. Searching for experiments
2. Checking TSV availability (for systems without R)
3. Downloading experiment data
4. Accessing expression matrices and annotations
"""

from expression_atlas import (
    ExpressionAtlasClient,
    has_r_available,
    has_tsv_files,
)

# ANSI colors for nice terminal output
BOLD = "\033[1m"
CYAN = "\033[96m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
MAGENTA = "\033[95m"
RESET = "\033[0m"


def header(text: str) -> None:
    """Print a formatted header."""
    print(f"\n{BOLD}{CYAN}{'=' * 60}{RESET}")
    print(f"{BOLD}{CYAN}{text}{RESET}")
    print(f"{BOLD}{CYAN}{'=' * 60}{RESET}\n")


def subheader(text: str) -> None:
    """Print a formatted subheader."""
    print(f"\n{YELLOW}>>> {text}{RESET}\n")


def main() -> None:
    """Run the demo."""
    print(f"\n{BOLD}{MAGENTA}🧬 Expression Atlas Python Client Demo 🧬{RESET}")

    # =========================================================================
    # 0. CHECK ENVIRONMENT
    # =========================================================================
    header("0. Environment Check")

    r_available = has_r_available()
    print(f"R + rpy2 available: {GREEN if r_available else RED}{r_available}{RESET}")

    if r_available:
        print("  → Can download ANY experiment (.Rdata files)")
    else:
        print("  → Can only download experiments with TSV files")
        print("  → Install R for full functionality")

    # Initialize the client
    client = ExpressionAtlasClient()

    # =========================================================================
    # 1. SEARCH FOR EXPERIMENTS
    # =========================================================================
    header("1. Searching for Experiments")

    subheader("Search for 'cancer' experiments in Homo sapiens")
    results = client.search_experiments(properties=["cancer"], species="homo sapiens")

    print(f"{GREEN}Found {len(results)} experiments!{RESET}\n")
    print(f"DataFrame columns: {list(results.columns)}")
    print(f"\nFirst 5 results:")
    print(results[["Accession", "Type"]].head().to_string(index=False))

    # =========================================================================
    # 2. CHECK TSV AVAILABILITY
    # =========================================================================
    header("2. Checking TSV File Availability")

    # Check a few experiments for TSV files
    test_accessions = results["Accession"].head(3).tolist()

    print(f"Checking {len(test_accessions)} experiments for TSV files...\n")

    experiments_with_tsv = []
    for acc in test_accessions:
        has_tsv = has_tsv_files(acc)
        status = f"{GREEN}✓ TSV{RESET}" if has_tsv else f"{RED}✗ No TSV{RESET}"
        print(f"  {acc}: {status}")
        if has_tsv:
            experiments_with_tsv.append(acc)

    print(f"\n{GREEN}{len(experiments_with_tsv)}{RESET} of {len(test_accessions)} have TSV files available")

    # =========================================================================
    # 3. DOWNLOAD AN EXPERIMENT
    # =========================================================================
    header("3. Downloading Experiment Data")

    # Choose which experiment to download
    if r_available:
        # With R, we can download any experiment
        accession = test_accessions[0]
        print(f"Using R to download: {accession}")
    elif experiments_with_tsv:
        # Without R, use one with TSV files
        accession = experiments_with_tsv[0]
        print(f"Using TSV fallback to download: {accession}")
    else:
        # Try some known experiments that often have TSV files
        known_tsv_experiments = ["E-MTAB-513", "E-GEOD-26284", "E-MTAB-62"]
        accession = None
        print("No TSV files found in search results. Trying known experiments...")
        for acc in known_tsv_experiments:
            if has_tsv_files(acc):
                accession = acc
                print(f"  Found: {accession}")
                break

    if accession is None:
        print(f"\n{RED}Could not find any experiment with TSV files.{RESET}")
        print("Install R + rpy2 for full functionality.")
        return

    subheader(f"Downloading: {accession}")
    experiment = client.get_experiment(accession)

    if experiment is None:
        print(f"{RED}Download failed.{RESET}")
        return

    print(f"{GREEN}Successfully downloaded!{RESET}\n")
    print(f"Experiment container: {experiment}")
    print(f"Available keys: {list(experiment.keys())}")

    # =========================================================================
    # 4. ACCESS THE DATA (RNA-seq: SummarizedExperiment)
    # =========================================================================
    header("4. Accessing RNA-seq Data (SummarizedExperiment)")

    # Check what type of data we have
    if "rnaseq" in experiment:
        sumexp = experiment["rnaseq"]

        print(f"{BOLD}Object structure (mirrors R's SummarizedExperiment):{RESET}")
        print(f"{sumexp}\n")

        # ----- RAW COUNTS MATRIX -----
        subheader("4a. Raw Counts Matrix - assays['counts']")
        print("This is the UN-normalized count data for DESeq2/edgeR analysis\n")

        counts = sumexp.assays.get("counts")
        if counts is not None and counts.size > 0:
            print(f"  Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
            print(f"  Data type: {counts.dtype}")
            print(f"  Total counts: {counts.sum():,.0f}")
            print(f"\n  First 5x5 of counts matrix:")
            print(f"  {counts[:5, :5]}")

        # ----- GENE ANNOTATIONS (rowData) -----
        subheader("4b. Gene Annotations - rowData")
        print("Gene-level metadata (IDs, names, etc.)\n")

        print(f"  Gene IDs (rownames): {len(sumexp.rownames)} genes")
        print(f"  First 5: {sumexp.rownames[:5]}")

        if not sumexp.rowData.empty:
            print(f"\n  rowData columns: {list(sumexp.rowData.columns)}")
            print(f"\n  First 5 rows:")
            print(sumexp.rowData.head().to_string())
        else:
            print("  (No additional gene annotations in TSV mode)")

        # ----- SAMPLE ANNOTATIONS (colData) -----
        subheader("4c. Sample Annotations - colData")
        print("Sample metadata for DEG design (condition, batch, sex, etc.)")
        print("This is what you use as design matrix in DESeq2/edgeR!\n")

        print(f"  Sample IDs (colnames): {len(sumexp.colnames)} samples")
        print(f"  First 5: {sumexp.colnames[:5]}")

        if not sumexp.colData.empty:
            print(f"\n  colData columns: {list(sumexp.colData.columns)}")
            print(f"\n  First 5 rows:")
            print(sumexp.colData.head().to_string())
        else:
            print("\n  (Sample annotations not available in TSV-only mode)")
            print("  Install R for full colData from .Rdata files")

        # ----- EXPERIMENT METADATA -----
        subheader("4d. Experiment Metadata - metadata")
        print("Experiment-level information\n")

        if sumexp.metadata:
            for key, value in sumexp.metadata.items():
                print(f"  {key}: {value}")
        else:
            print("  (No experiment metadata)")

    # =========================================================================
    # 5. MICROARRAY DATA (ExpressionSet) - if available
    # =========================================================================
    array_keys = [k for k in experiment.keys() if k.startswith("A-")]
    if array_keys or "normalized" in experiment:
        header("5. Accessing Microarray Data (ExpressionSet)")

        if array_keys:
            eset = experiment[array_keys[0]]
            print(f"Array design: {array_keys[0]}\n")
        else:
            eset = experiment["normalized"]

        print(f"{BOLD}Object structure (mirrors R's ExpressionSet):{RESET}")
        print(f"{eset}\n")

        subheader("5a. Normalized Expression Matrix - exprs")
        print("Pre-normalized intensities (NOT raw counts)\n")
        if eset.exprs.size > 0:
            print(f"  Shape: {eset.exprs.shape[0]} probes × {eset.exprs.shape[1]} samples")

        subheader("5b. Sample Annotations - phenoData / pData")
        print("R-style accessors: eset.pData is alias for eset.phenoData\n")
        if not eset.phenoData.empty:
            print(f"  Columns: {list(eset.phenoData.columns)}")
        print(f"  eset.pData is eset.phenoData: {eset.pData is eset.phenoData}")

        subheader("5c. Probe Annotations - featureData / fData")
        print("R-style accessors: eset.fData is alias for eset.featureData\n")
        if not eset.featureData.empty:
            print(f"  Columns: {list(eset.featureData.columns)}")
        print(f"  eset.fData is eset.featureData: {eset.fData is eset.featureData}")

    # =========================================================================
    # 6. SUMMARY - What You Get
    # =========================================================================
    header("6. Summary: What You Get from Expression Atlas")

    print(f"""
{BOLD}For RNA-seq experiments (SummarizedExperiment):{RESET}

  {GREEN}sumexp.assays['counts']{RESET}  → Raw counts matrix (genes × samples)
                             Use this for DESeq2/edgeR DEG analysis!

  {GREEN}sumexp.colData{RESET}           → Sample annotations (condition, batch, etc.)
                             Use as design matrix for DEG analysis!

  {GREEN}sumexp.rowData{RESET}           → Gene annotations (gene names, biotypes)

  {GREEN}sumexp.rownames{RESET}          → Gene IDs (Ensembl IDs)
  {GREEN}sumexp.colnames{RESET}          → Sample IDs

  {GREEN}sumexp.metadata{RESET}          → Experiment-level info

{BOLD}For Microarray experiments (ExpressionSet):{RESET}

  {GREEN}eset.exprs{RESET}               → Normalized intensity matrix
  {GREEN}eset.phenoData / pData{RESET}   → Sample annotations
  {GREEN}eset.featureData / fData{RESET} → Probe/gene annotations
  {GREEN}eset.experimentData{RESET}      → Experiment metadata

{BOLD}Search Results DataFrame:{RESET}

  {GREEN}Accession{RESET}  → E-MTAB-xxxx, E-GEOD-xxxx
  {GREEN}Species{RESET}    → Organism
  {GREEN}Type{RESET}       → RNA-seq, microarray, etc.
  {GREEN}Title{RESET}      → Experiment title

{BOLD}This demo:{RESET}
  • Found {len(results)} experiments matching 'cancer' in human
  • Downloaded {accession} with {sumexp.n_genes if 'rnaseq' in experiment else 0} genes × {sumexp.n_samples if 'rnaseq' in experiment else 0} samples
""")

    print(f"{BOLD}{MAGENTA}✨ Demo complete! ✨{RESET}\n")


if __name__ == "__main__":
    main()
