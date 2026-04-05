---
title: Transcriptomics Explorer
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_file: app.py
pinned: false
---

# Transcriptomics Explorer

An interactive, browser-based dashboard for exploring **RNA-seq differential expression results**. Upload your DESeq2 or edgeR output and instantly visualise it across six linked plots. No database required — just upload a CSV and explore.

> **Demo mode:** The app loads 800 simulated genes on startup so you can explore all features immediately without uploading anything.

---

## What It Does

| Tab | What you see | Data needed |
|-----|-------------|-------------|
| **Volcano Plot** | Which genes are significantly up/down regulated | DESeq2/edgeR CSV |
| **MA Plot** | Expression magnitude vs. fold change (QC view) | DESeq2/edgeR CSV |
| **PCA** | How your samples cluster by gene expression | DESeq2/edgeR CSV |
| **Sample Correlation** | Pearson correlation between samples | DESeq2/edgeR CSV |
| **Heatmap** | Top variable genes across samples (Z-scored) | DESeq2/edgeR CSV |
| **Pathway Enrichment** | Which biological pathways are activated/suppressed | DESeq2/edgeR CSV |

All six tabs update automatically when you upload your data — no page refresh needed.

### Key features
- Upload DESeq2 or edgeR output (CSV or TSV) — all tabs update instantly
- Volcano plot thresholds (log₂FC, adjusted p-value) are interactive sliders
- Box-select points on the Volcano plot → linked data table updates instantly
- Export significant genes as CSV with one click
- Data table is searchable, sortable, and filterable
- Pathway enrichment bar chart shows Up vs. Down gene separation per pathway

---

## Files

```
transcriptomics-explorer/
├── app.py                   # The entire dashboard (single file)
├── requirements.txt         # Python dependencies
├── Dockerfile               # For Hugging Face Spaces deployment
├── test_data_1000genes.csv  # Sample file to try the upload feature
├── GUIDE.md                 # Full user guide with interpretation help
├── README.md
└── analysis/                # R scripts for generating de_results.csv
    ├── fetch_geo.R          # Get count matrix (airway dataset or GEO)
    ├── deseq2_analysis.R    # Run DESeq2 → de_results.csv
    ├── edger_analysis.R     # Alternative: run edgeR → de_results.csv
    ├── counts.csv           # Count matrix (generated)
    ├── metadata.csv         # Sample info (generated)
    ├── de_results.csv       # Final output → upload to dashboard
    └── README.md            # Step-by-step analysis instructions
```

---

## Setup & Run

### Option A — pip (any Python environment)

```bash
cd transcriptomics-explorer
pip install -r requirements.txt
python app.py
```

Then open **http://localhost:7860** in your browser.

### Option B — conda (recommended if you use Anaconda/Miniconda)

```bash
conda create -n transcriptomics python=3.11
conda activate transcriptomics
pip install -r requirements.txt
python app.py
```

> **Why pip inside conda?** These packages aren't all on `conda-forge` with matching versions, so `pip` is the safest install method regardless of environment manager.

Press `Ctrl + C` in the terminal to stop the app.

---

## Uploading Your Own Data

The app accepts **DESeq2 or edgeR results** as CSV or TSV. Column names are flexible — common variants are auto-detected.

| Column | Accepted names | Example |
|--------|---------------|---------|
| Gene symbol | `gene`, `gene_id`, `geneId`, `symbol` | `TP53` |
| Log₂ fold change | `log2FoldChange`, `log2fc`, `lfc` | `2.41` |
| Raw p-value | `pvalue`, `pval`, `p.value` | `0.000032` |
| Adjusted p-value | `padj`, `p.adj`, `fdr` | `0.0012` |
| Base mean (optional) | `baseMean`, `mean_expr`, `avgExpr` | `340.5` |

A sample file (`test_data_1000genes.csv`) is included — upload it to try all features.

---

## Dependency Compatibility

`requirements.txt` uses minimum version ranges (`>=`) rather than exact pins, so it works alongside whatever you already have installed.

- Python 3.9 – 3.12
- Dash 2.14+
- Pandas 1.5+
- NumPy 1.23+
- SciPy 1.9+
