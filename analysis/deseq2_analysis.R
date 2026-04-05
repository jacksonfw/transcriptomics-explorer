# =============================================================================
# deseq2_analysis.R
# Differential expression analysis using DESeq2, with gene symbol annotation
#
# Input:
#   counts.csv    — raw integer count matrix (genes × samples)
#                   produced by fetch_geo.R, or your own count matrix
#   metadata.csv  — sample info table with at least:
#                   - a column matching sample names
#                   - a column for the condition/group label
#
# Output:
#   de_results.csv   — DESeq2 results with gene symbols → upload to dashboard
#   ma_plot.pdf      — MA plot QC
#   volcano_plot.pdf — Volcano plot preview
#
# Usage:
#   Rscript deseq2_analysis.R
#   or open in RStudio and run interactively
# =============================================================================

# ── 0. Install packages if missing ──────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("DESeq2", "ashr", "org.Hs.eg.db", "AnnotationDbi")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}
for (pkg in c("ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ── 1. USER SETTINGS — edit these ───────────────────────────────────────────

# Path to count matrix (rows = genes, columns = samples)
COUNTS_FILE <- "counts.csv"

# Path to sample metadata
METADATA_FILE <- "metadata.csv"

# Which column in metadata contains the sample IDs
# (must match the column names in counts.csv exactly)
# For airway dataset: use "sample_id"
# For GEO datasets:   use "geo_accession"
SAMPLE_ID_COL <- "sample_id"

# Which column in metadata contains the condition/group label
# For airway dataset: use "dex"
# For GEO datasets:   check your metadata.csv column names
CONDITION_COL <- "dex"

# The two conditions to compare (exactly as they appear in CONDITION_COL)
# GROUP_B is compared against GROUP_A (positive LFC = higher in GROUP_B)
# For airway dataset:
GROUP_A <- "untrt"   # reference / control
GROUP_B <- "trt"     # treatment / case

# What type of gene IDs are in your count matrix row names?
# "ENSEMBL"    → ENSG00000... (most common, used by airway and STAR/HISAT2)
# "ENTREZID"   → numeric IDs like 7157
# "SYMBOL"     → already gene symbols like TP53 (skip annotation step)
GENE_ID_TYPE <- "ENSEMBL"

# Which organism annotation database to use for gene symbol mapping
# Human:  org.Hs.eg.db  (loaded above)
# Mouse:  org.Mm.eg.db  (install separately: BiocManager::install("org.Mm.eg.db"))
# Other organisms: see https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
ANNOTATION_DB <- org.Hs.eg.db

# Significance thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0   # minimum absolute log2 fold change

# Output file name
OUTPUT_FILE <- "de_results.csv"

# ── 2. Load data ─────────────────────────────────────────────────────────────
cat("Loading count matrix from:", COUNTS_FILE, "\n")
counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))
counts[counts < 0] <- 0
cat("Count matrix:", nrow(counts), "genes ×", ncol(counts), "samples\n")

cat("Loading metadata from:", METADATA_FILE, "\n")
meta <- read.csv(METADATA_FILE)

# Keep only the two groups we want to compare
meta <- meta[meta[[CONDITION_COL]] %in% c(GROUP_A, GROUP_B), ]
cat("Samples:", nrow(meta),
    "(", sum(meta[[CONDITION_COL]] == GROUP_A), "×", GROUP_A,
    "/",  sum(meta[[CONDITION_COL]] == GROUP_B), "×", GROUP_B, ")\n")

# Align count matrix columns to metadata rows
sample_ids <- meta[[SAMPLE_ID_COL]]
counts <- counts[, colnames(counts) %in% sample_ids, drop = FALSE]
counts <- counts[, meta[[SAMPLE_ID_COL]], drop = FALSE]

# ── 3. Build DESeq2 object ───────────────────────────────────────────────────
coldata <- data.frame(
  row.names = meta[[SAMPLE_ID_COL]],
  condition = factor(meta[[CONDITION_COL]], levels = c(GROUP_A, GROUP_B))
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

# Pre-filter: remove genes with very low counts
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
cat("Genes after pre-filter:", nrow(dds), "\n")

# ── 4. Run DESeq2 ────────────────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeq(dds)

res <- results(dds,
               contrast = c("condition", GROUP_B, GROUP_A),
               alpha    = PADJ_CUTOFF)

# LFC shrinkage — reduces noise for lowly expressed genes (better for visualisation)
res_shrunk <- lfcShrink(dds,
                         contrast = c("condition", GROUP_B, GROUP_A),
                         type     = "ashr",
                         res      = res)

cat("\nDESeq2 summary:\n")
summary(res_shrunk, alpha = PADJ_CUTOFF)

# ── 5. Annotate Ensembl IDs → gene symbols ───────────────────────────────────
res_df <- as.data.frame(res_shrunk)
res_df$ensembl_id <- rownames(res_df)
rownames(res_df) <- NULL

if (GENE_ID_TYPE != "SYMBOL") {
  cat("\nMapping", GENE_ID_TYPE, "IDs to gene symbols...\n")

  symbols <- mapIds(
    ANNOTATION_DB,
    keys      = res_df$ensembl_id,
    column    = "SYMBOL",
    keytype   = GENE_ID_TYPE,
    multiVals = "first"   # if one Ensembl ID maps to multiple symbols, take first
  )

  res_df$gene <- ifelse(is.na(symbols), res_df$ensembl_id, symbols)
  n_mapped <- sum(!is.na(symbols))
  cat("Mapped", n_mapped, "of", nrow(res_df), "genes to symbols\n")
  cat("Unmapped genes will keep their", GENE_ID_TYPE, "ID\n")
} else {
  # Gene IDs are already symbols
  res_df$gene <- res_df$ensembl_id
}

# ── 6. Format and save results ───────────────────────────────────────────────
res_df <- res_df %>%
  dplyr::select(gene, ensembl_id, baseMean, log2FoldChange, pvalue, padj) %>%
  dplyr::mutate(
    regulation = dplyr::case_when(
      !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up",
      !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

n_up   <- sum(res_df$regulation == "Up",   na.rm = TRUE)
n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
cat("\nSignificant DE genes:", n_up + n_down,
    "(", n_up, "up,", n_down, "down )\n")

write.csv(res_df, OUTPUT_FILE, row.names = FALSE)
cat("Saved", OUTPUT_FILE, "\n")

# ── 7. QC plots ──────────────────────────────────────────────────────────────
pdf("ma_plot.pdf", width = 7, height = 5)
plotMA(res_shrunk, alpha = PADJ_CUTOFF,
       main = paste("MA plot:", GROUP_B, "vs", GROUP_A))
dev.off()
cat("Saved ma_plot.pdf\n")

res_plot <- res_df %>% dplyr::filter(!is.na(padj))
pdf("volcano_plot.pdf", width = 7, height = 5)
ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), colour = regulation)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_colour_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#aaaaaa")) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", colour = "black") +
  labs(title = paste("Volcano:", GROUP_B, "vs", GROUP_A),
       x = "log2 Fold Change", y = "-log10(padj)") +
  theme_bw()
dev.off()
cat("Saved volcano_plot.pdf\n")

cat("\nAll done!\n")
cat("Upload", OUTPUT_FILE, "to the Transcriptomics Explorer dashboard.\n")
