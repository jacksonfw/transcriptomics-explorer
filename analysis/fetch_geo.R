# =============================================================================
# fetch_geo.R
# Get a small RNA-seq count matrix ready for DESeq2/edgeR analysis
#
# This script has TWO options:
#
# OPTION A (recommended for testing) — uses the built-in 'airway' Bioconductor
#   dataset. 8 samples, 63k genes, no internet download needed. This is the
#   dataset used to produce the included de_results.csv.
#   Experiment: dexamethasone drug treatment vs untreated, human airway cells.
#
# OPTION B — download directly from GEO using GEOquery.
#   Works for any GSE accession that provides a supplementary count matrix file.
#   NOTE: GEO datasets vary in format — you may need to adjust the count file
#   parsing depending on the dataset.
#
# Output files:
#   counts.csv    — raw count matrix (genes × samples)
#   metadata.csv  — sample information (condition labels)
# =============================================================================

# ── Set which option to use ──────────────────────────────────────────────────
USE_AIRWAY <- TRUE    # TRUE = Option A (airway), FALSE = Option B (GEO)
GEO_ID     <- "GSE96870"   # Only used when USE_AIRWAY = FALSE

# ── Install packages if missing ──────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# ════════════════════════════════════════════════════════════════════════════
# OPTION A — airway built-in dataset (fast, no download required)
# ════════════════════════════════════════════════════════════════════════════
if (USE_AIRWAY) {

  for (pkg in c("airway", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE)
  }

  library(airway)
  library(SummarizedExperiment)

  data(airway)
  cat("Loaded airway dataset:", nrow(airway), "genes x", ncol(airway), "samples\n")
  cat("Experiment: dexamethasone treatment vs untreated, human airway smooth muscle cells\n")

  # Raw count matrix
  counts <- assay(airway, "counts")

  # Sample metadata
  meta <- as.data.frame(colData(airway))
  meta$sample_id <- rownames(meta)

  cat("\nSample groups:\n")
  print(table(meta$dex))

  write.csv(as.data.frame(counts), "counts.csv",   row.names = TRUE)
  write.csv(meta,                  "metadata.csv",  row.names = FALSE)

  cat("\nSaved counts.csv and metadata.csv\n")
  cat("\nNext steps:\n")
  cat("  - Open deseq2_analysis.R\n")
  cat("  - Set CONDITION_COL <- 'dex'\n")
  cat("  - Set GROUP_A <- 'untrt'  (reference)\n")
  cat("  - Set GROUP_B <- 'trt'    (treated)\n")
  cat("  - Run deseq2_analysis.R → produces de_results.csv\n")

# ════════════════════════════════════════════════════════════════════════════
# OPTION B — download from GEO
# ════════════════════════════════════════════════════════════════════════════
} else {

  if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery", ask = FALSE)
  library(GEOquery)

  cat("Fetching", GEO_ID, "from GEO...\n")

  # Download series matrix (metadata)
  gse <- getGEO(GEO_ID, GSEMatrix = TRUE, destdir = ".", AnnotGPL = FALSE)
  if (is.list(gse)) gse <- gse[[1]]

  # Save metadata
  meta <- pData(phenoData(gse))
  keep_cols <- grep(
    "title|source|characteristics|condition|group|treatment|tissue|sex|geo_accession",
    colnames(meta), ignore.case = TRUE, value = TRUE
  )
  meta_clean <- meta[, c("geo_accession", keep_cols), drop = FALSE]
  write.csv(meta_clean, "metadata.csv", row.names = FALSE)
  cat("Saved metadata.csv —", nrow(meta_clean), "samples\n")

  # Print column names and group values so user can set CONDITION_COL correctly
  cat("\nMetadata columns:\n")
  print(colnames(meta_clean))
  cat("\nSample titles:\n")
  print(meta_clean$title)

  # Download supplementary count file
  cat("\nDownloading supplementary count matrix...\n")
  supp_files <- getGEOSuppFiles(GEO_ID, makeDirectory = FALSE)
  supp_path  <- rownames(supp_files)
  count_file <- supp_path[grepl("count|raw|matrix", basename(supp_path),
                                 ignore.case = TRUE)][1]

  if (is.na(count_file)) {
    stop(
      "Could not auto-detect count file. Downloaded: ",
      paste(basename(supp_path), collapse = ", "),
      "\nSet count_file manually."
    )
  }

  cat("Reading:", basename(count_file), "\n")
  counts_raw <- read.table(gzfile(count_file), header = TRUE, sep = "\t",
                            row.names = 1, check.names = FALSE)
  cat("Dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")

  # Filter low-count genes
  keep <- rowSums(counts_raw >= 10) >= 3
  counts_filtered <- counts_raw[keep, ]
  cat("After filtering:", nrow(counts_filtered), "genes retained\n")

  write.csv(counts_filtered, "counts.csv", row.names = TRUE)
  cat("Saved counts.csv\n")

  cat("\nNext steps:\n")
  cat("  - Open metadata.csv and check the column names and group labels\n")
  cat("  - Set CONDITION_COL, GROUP_A, GROUP_B in deseq2_analysis.R accordingly\n")
  cat("  - Run deseq2_analysis.R → produces de_results.csv\n")
}
