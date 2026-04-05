# =============================================================================
# edger_analysis.R
# Differential expression analysis using edgeR, with gene symbol annotation
#
# Alternative to deseq2_analysis.R — use either one, not both.
# edgeR tends to be more sensitive for experiments with fewer replicates (< 3).
# DESeq2 is generally preferred when you have 3+ replicates per group.
#
# Input:
#   counts.csv    — raw integer count matrix (genes × samples)
#   metadata.csv  — sample info with condition column
#
# Output:
#   de_results.csv   — edgeR results with gene symbols → upload to dashboard
#   ma_plot.pdf      — MA plot QC
#   volcano_plot.pdf — Volcano plot preview
#
# Usage:
#   Rscript edger_analysis.R
#   or open in RStudio and run interactively
# =============================================================================

# ── 0. Install packages if missing ──────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("edgeR", "limma", "org.Hs.eg.db", "AnnotationDbi")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}
for (pkg in c("ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(edgeR)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ── 1. USER SETTINGS — edit these ───────────────────────────────────────────

COUNTS_FILE   <- "counts.csv"
METADATA_FILE <- "metadata.csv"

# Which column in metadata contains the sample IDs
# For airway dataset: use "sample_id"
# For GEO datasets:   use "geo_accession"
SAMPLE_ID_COL <- "sample_id"

# Which column in metadata contains the condition/group label
# For airway dataset: use "dex"
# For GEO datasets:   check your metadata.csv column names
CONDITION_COL <- "dex"

# For airway dataset:
GROUP_A <- "untrt"   # reference / control
GROUP_B <- "trt"     # treatment / case

# Gene ID type in your count matrix row names
# "ENSEMBL"  → ENSG00000... (most common)
# "ENTREZID" → numeric IDs
# "SYMBOL"   → already gene symbols (skip annotation)
GENE_ID_TYPE <- "ENSEMBL"

# Annotation database (human by default)
# For mouse: install org.Mm.eg.db and set ANNOTATION_DB <- org.Mm.eg.db
ANNOTATION_DB <- org.Hs.eg.db

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

OUTPUT_FILE <- "de_results.csv"

# ── 2. Load data ─────────────────────────────────────────────────────────────
cat("Loading count matrix from:", COUNTS_FILE, "\n")
counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))
counts[counts < 0] <- 0
cat("Count matrix:", nrow(counts), "genes ×", ncol(counts), "samples\n")

cat("Loading metadata from:", METADATA_FILE, "\n")
meta <- read.csv(METADATA_FILE)

meta <- meta[meta[[CONDITION_COL]] %in% c(GROUP_A, GROUP_B), ]
cat("Samples:", nrow(meta),
    "(", sum(meta[[CONDITION_COL]] == GROUP_A), "×", GROUP_A,
    "/",  sum(meta[[CONDITION_COL]] == GROUP_B), "×", GROUP_B, ")\n")

sample_ids <- meta[[SAMPLE_ID_COL]]
counts <- counts[, colnames(counts) %in% sample_ids, drop = FALSE]
counts <- counts[, meta[[SAMPLE_ID_COL]], drop = FALSE]

# ── 3. Build edgeR DGEList ───────────────────────────────────────────────────
group <- factor(meta[[CONDITION_COL]], levels = c(GROUP_A, GROUP_B))
dge   <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(dge, group = group, min.count = 10, min.total.count = 15)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes after filtering:", nrow(dge), "\n")

# ── 4. Normalise and estimate dispersion ─────────────────────────────────────
dge    <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ group)
dge    <- estimateDisp(dge, design, robust = TRUE)
cat("Common BCV:", round(sqrt(dge$common.dispersion), 3), "\n")

# ── 5. Fit model and test ─────────────────────────────────────────────────────
fit <- glmQLFit(dge, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)

res    <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "PValue")
res_df <- as.data.frame(res$table)
res_df$ensembl_id <- rownames(res_df)
rownames(res_df) <- NULL

# ── 6. Annotate gene IDs → symbols ───────────────────────────────────────────
if (GENE_ID_TYPE != "SYMBOL") {
  cat("\nMapping", GENE_ID_TYPE, "IDs to gene symbols...\n")

  symbols <- mapIds(
    ANNOTATION_DB,
    keys      = res_df$ensembl_id,
    column    = "SYMBOL",
    keytype   = GENE_ID_TYPE,
    multiVals = "first"
  )

  res_df$gene <- ifelse(is.na(symbols), res_df$ensembl_id, symbols)
  n_mapped <- sum(!is.na(symbols))
  cat("Mapped", n_mapped, "of", nrow(res_df), "genes to symbols\n")
} else {
  res_df$gene <- res_df$ensembl_id
}

# ── 7. Format and save results ────────────────────────────────────────────────
# edgeR doesn't output baseMean — use average CPM as proxy
res_df$baseMean <- rowMeans(cpm(dge, log = FALSE))[res_df$ensembl_id]

res_df <- res_df %>%
  dplyr::rename(log2FoldChange = logFC, pvalue = PValue, padj = FDR) %>%
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

# ── 8. QC plots ───────────────────────────────────────────────────────────────
pdf("ma_plot.pdf", width = 7, height = 5)
plotMD(qlf, main = paste("MA plot:", GROUP_B, "vs", GROUP_A))
abline(h = c(-LFC_CUTOFF, LFC_CUTOFF), col = "blue", lty = 2)
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
