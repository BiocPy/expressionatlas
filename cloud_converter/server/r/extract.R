#!/usr/bin/env Rscript
#
# Expression Atlas .RData Extraction Script
#
# Extracts expression matrices and metadata from Expression Atlas .RData files
# and converts them to portable formats (Matrix Market, CSV).
#
# Usage:
#   Rscript extract.R --input /path/to/file.RData --outdir /path/to/output --format mtx_bundle
#
# Output:
#   outdir/
#     dataset_<name>/
#       matrix.mtx (or counts.tsv.gz for small matrices)
#       genes.csv
#       samples.csv
#       barcodes.tsv
#       features.tsv
#     meta.json

suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(jsonlite)
})

# Try to load Bioconductor packages (may not be installed)
tryCatch({
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(S4Vectors)
  })
}, error = function(e) {
  message("Note: Bioconductor packages not fully available")
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  result <- list(
    input = NULL,
    outdir = NULL,
    format = "mtx_bundle",
    assay_name = NULL
  )
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" && i < length(args)) {
      result$input <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--outdir" && i < length(args)) {
      result$outdir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--format" && i < length(args)) {
      result$format <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--assay_name" && i < length(args)) {
      result$assay_name <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  return(result)
}

params <- parse_args(args)

if (is.null(params$input) || is.null(params$outdir)) {
  stop("Usage: Rscript extract.R --input <file.RData> --outdir <output_dir> [--format mtx_bundle] [--assay_name <name>]")
}

message("=== Expression Atlas RData Extraction ===")
message(sprintf("Input: %s", params$input))
message(sprintf("Output: %s", params$outdir))
message(sprintf("Format: %s", params$format))

# Create output directory
dir.create(params$outdir, recursive = TRUE, showWarnings = FALSE)

# Load the RData file
message("Loading RData file...")
env <- new.env()
load(params$input, envir = env)
objects <- ls(env)
message(sprintf("Found objects: %s", paste(objects, collapse = ", ")))

# Find the main experiment object
find_main_object <- function(env) {
  objects <- ls(env)
  
  # Prefer 'experimentSummary' or 'experiment_summary'
  preferred <- c("experimentSummary", "experiment_summary", "atlasExperimentSummary")
  for (name in preferred) {
    if (name %in% objects) {
      message(sprintf("Using preferred object: %s", name))
      return(list(name = name, obj = get(name, envir = env)))
    }
  }
  
  # Otherwise find the largest object or known container type
  for (name in objects) {
    obj <- get(name, envir = env)
    cls <- class(obj)[1]
    if (cls %in% c("SimpleList", "SummarizedExperiment", "RangedSummarizedExperiment", "ExpressionSet")) {
      message(sprintf("Using object '%s' of class %s", name, cls))
      return(list(name = name, obj = obj))
    }
  }
  
  # Fallback: use first object
  if (length(objects) > 0) {
    name <- objects[1]
    message(sprintf("Fallback: using first object '%s'", name))
    return(list(name = name, obj = get(name, envir = env)))
  }
  
  stop("No suitable object found in RData file")
}

main_obj <- find_main_object(env)

# Normalize to list of datasets
normalize_to_datasets <- function(obj, name) {
  cls <- class(obj)[1]
  
  if (cls == "SimpleList" || cls == "list") {
    # Iterate each element
    datasets <- list()
    elem_names <- names(obj)
    if (is.null(elem_names)) elem_names <- paste0("dataset_", seq_along(obj))
    
    for (i in seq_along(obj)) {
      elem <- obj[[i]]
      elem_name <- elem_names[i]
      datasets[[elem_name]] <- elem
    }
    return(datasets)
  } else {
    # Single object
    return(list(dataset_1 = obj))
  }
}

datasets <- normalize_to_datasets(main_obj$obj, main_obj$name)
message(sprintf("Found %d dataset(s): %s", length(datasets), paste(names(datasets), collapse = ", ")))

# Helper: write matrix to file
write_matrix <- function(mat, outdir, format, name) {
  if (format == "mtx_bundle") {
    # Write as Matrix Market format
    mtx_file <- file.path(outdir, "matrix.mtx")
    
    # Convert to sparse if not already
    if (!inherits(mat, "sparseMatrix")) {
      mat <- as(mat, "dgCMatrix")
    }
    
    writeMM(mat, mtx_file)
    message(sprintf("  Wrote matrix.mtx (%d x %d)", nrow(mat), ncol(mat)))
    
    # Write barcodes (column names) and features (row names)
    if (!is.null(colnames(mat))) {
      write.table(colnames(mat), file.path(outdir, "barcodes.tsv"), 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    if (!is.null(rownames(mat))) {
      write.table(rownames(mat), file.path(outdir, "features.tsv"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  } else {
    # Write as gzipped TSV
    tsv_file <- file.path(outdir, "counts.tsv.gz")
    gz <- gzfile(tsv_file, "w")
    write.table(as.matrix(mat), gz, sep = "\t", quote = FALSE)
    close(gz)
    message(sprintf("  Wrote counts.tsv.gz (%d x %d)", nrow(mat), ncol(mat)))
  }
}

# Helper: write data.frame to CSV
write_df <- function(df, filepath) {
  if (is.null(df) || nrow(df) == 0) {
    # Write empty file with just header
    write.csv(data.frame(row.names = character(0)), filepath, row.names = TRUE)
  } else {
    write.csv(as.data.frame(df), filepath, row.names = TRUE)
  }
}

# Process each dataset
meta_info <- list(
  datasets = list(),
  r_version = R.version.string,
  bioconductor_version = tryCatch(as.character(BiocManager::version()), error = function(e) NULL)
)

for (ds_name in names(datasets)) {
  message(sprintf("\nProcessing dataset: %s", ds_name))
  ds <- datasets[[ds_name]]
  cls <- class(ds)[1]
  
  # Create dataset output directory
  ds_outdir <- file.path(params$outdir, paste0("dataset_", gsub("[^a-zA-Z0-9_-]", "_", ds_name)))
  dir.create(ds_outdir, recursive = TRUE, showWarnings = FALSE)
  
  ds_meta <- list(
    name = ds_name,
    class_type = cls,
    dimensions = c(0, 0),
    assay_names = character(0),
    row_data_columns = character(0),
    col_data_columns = character(0)
  )
  
  tryCatch({
    if (cls %in% c("SummarizedExperiment", "RangedSummarizedExperiment")) {
      message("  Class: SummarizedExperiment")
      
      # Get assay names
      assay_names <- assayNames(ds)
      ds_meta$assay_names <- as.list(assay_names)
      message(sprintf("  Assays: %s", paste(assay_names, collapse = ", ")))
      
      # Pick assay
      if (!is.null(params$assay_name) && params$assay_name %in% assay_names) {
        selected_assay <- params$assay_name
      } else {
        selected_assay <- assay_names[1]
      }
      message(sprintf("  Using assay: %s", selected_assay))
      
      # Extract matrix
      mat <- assay(ds, selected_assay)
      ds_meta$dimensions <- dim(mat)
      write_matrix(mat, ds_outdir, params$format, ds_name)
      
      # Extract rowData (gene annotations)
      rd <- rowData(ds)
      if (!is.null(rd)) {
        rd_df <- as.data.frame(rd)
        ds_meta$row_data_columns <- as.list(colnames(rd_df))
        write_df(rd_df, file.path(ds_outdir, "genes.csv"))
        message(sprintf("  Wrote genes.csv (%d rows, %d cols)", nrow(rd_df), ncol(rd_df)))
      }
      
      # Extract colData (sample annotations)
      cd <- colData(ds)
      if (!is.null(cd)) {
        cd_df <- as.data.frame(cd)
        ds_meta$col_data_columns <- as.list(colnames(cd_df))
        write_df(cd_df, file.path(ds_outdir, "samples.csv"))
        message(sprintf("  Wrote samples.csv (%d rows, %d cols)", nrow(cd_df), ncol(cd_df)))
      }
      
    } else if (cls == "ExpressionSet") {
      message("  Class: ExpressionSet")
      
      # Get expression matrix
      mat <- exprs(ds)
      ds_meta$assay_names <- list("exprs")
      ds_meta$dimensions <- dim(mat)
      write_matrix(mat, ds_outdir, params$format, ds_name)
      
      # Extract fData (feature annotations)
      fd <- fData(ds)
      if (!is.null(fd)) {
        ds_meta$row_data_columns <- as.list(colnames(fd))
        write_df(fd, file.path(ds_outdir, "genes.csv"))
        message(sprintf("  Wrote genes.csv (%d rows, %d cols)", nrow(fd), ncol(fd)))
      }
      
      # Extract pData (sample annotations)
      pd <- pData(ds)
      if (!is.null(pd)) {
        ds_meta$col_data_columns <- as.list(colnames(pd))
        write_df(pd, file.path(ds_outdir, "samples.csv"))
        message(sprintf("  Wrote samples.csv (%d rows, %d cols)", nrow(pd), ncol(pd)))
      }
      
    } else if (is.matrix(ds) || inherits(ds, "Matrix")) {
      message("  Class: matrix/Matrix")
      ds_meta$assay_names <- list("matrix")
      ds_meta$dimensions <- dim(ds)
      write_matrix(ds, ds_outdir, params$format, ds_name)
      
    } else {
      message(sprintf("  Warning: Unsupported class '%s', skipping", cls))
    }
    
  }, error = function(e) {
    message(sprintf("  Error processing dataset %s: %s", ds_name, e$message))
    ds_meta$error <- e$message
  })
  
  meta_info$datasets[[ds_name]] <- ds_meta
}

# Write meta.json
meta_file <- file.path(params$outdir, "meta.json")
write(toJSON(meta_info, auto_unbox = TRUE, pretty = TRUE), meta_file)
message(sprintf("\nWrote meta.json"))

message("\n=== Extraction complete ===")
