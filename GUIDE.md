# Transcriptomics Explorer — User Guide

This dashboard is a visualisation tool for exploring **RNA-seq differential expression results**. It is designed for researchers, students, and bioinformaticians. No coding required — just upload your data and interact with the plots.

**Default behaviour:** The app loads 800 simulated genes on startup so you can explore all features immediately. Upload your own DESeq2/edgeR CSV to replace the demo data across all six tabs at once.

---

## How to Run

1. Open Terminal
2. Navigate to this folder: `cd path/to/transcriptomics-explorer`
3. Install dependencies (first time only): `pip install -r requirements.txt`
4. Start the app: `python app.py`
5. Open your browser: **http://localhost:7860**
6. Press `Ctrl + C` in the terminal to stop.

---

## Tab 1 — Volcano Plot

### What is it?
The core feature of the dashboard. A Volcano plot answers the question:
**"Which genes are significantly turned ON or OFF between two conditions?"**
(e.g. cancer vs. normal tissue, drug-treated vs. untreated cells)

Each dot = one gene.

### How to read it

**X-axis — log₂ Fold Change** (how much expression changed)
- Positive (right) = gene is MORE active in condition B
- Negative (left) = gene is LESS active in condition B
- Value of 1 = 2× more expressed / Value of 2 = 4× more expressed

**Y-axis — -log₁₀(adjusted p-value)** (how confident we are)
- Higher = more statistically significant
- ~1.3 = p-value of 0.05 (standard significance cutoff)

**Colour coding**
- Red = Significantly UPREGULATED (top-right)
- Blue = Significantly DOWNREGULATED (top-left)
- Grey = Not significant

**What to look for:** The most interesting genes are in the top corners.

### Interactive controls
- **log₂FC slider** — raise the fold change threshold
- **p-value dropdown** — tighten or loosen significance
- **Gene search box** — type gene names (e.g. TP53, BRCA1) to highlight them as gold stars
- **Box-select on the plot** — draw a rectangle to populate the table below
- **Export CSV** — download the current significant gene list

---

## Tab 2 — MA Plot

### What is it?
A companion to the Volcano plot used for **quality control**. The MA plot shows whether fold changes are consistent across the range of expression levels.

Each dot = one gene.

### How to read it

**X-axis — log₁₀(mean expression)** — how highly expressed the gene is on average
- Left = lowly expressed genes (often less reliable)
- Right = highly expressed genes (more reliable measurements)

**Y-axis — log₂ Fold Change** — same as volcano X-axis

**Horizontal lines**
- Middle line at y=0 = no change
- Dashed lines at y=±1 = 2× fold change threshold

**Colour coding** — Red = up, Blue = down, Grey = not significant

**What to look for:**
- Fold changes should be roughly symmetric around y=0 for non-significant genes
- If low-expression genes (left side) show large scatter, that is expected noise — not biology
- Large fold changes at high expression (right side) are the most trustworthy hits

---

## Tab 3 — PCA (Principal Component Analysis)

### What is it?
PCA answers: **"How similar are your samples to each other overall?"**

It reduces thousands of gene measurements down to two axes (PC1 and PC2) that capture the most variation. Each dot = one sample.

### How to read it

**PC1 (X-axis)** — the biggest source of variation in your data
**PC2 (Y-axis)** — the second biggest source of variation

The % in each axis label tells you how much of the total variation that component explains.

**What to look for:**
- Samples from the same condition (e.g. replicates) should cluster together
- The two conditions should separate along PC1 (the main axis)
- An outlier sample far from its group may indicate a technical problem

> **Note:** In this app, PCA is derived from the DE results matrix. For true sample-level PCA, provide a count matrix (genes × samples) as a second upload.

---

## Tab 4 — Sample Correlation

### What is it?
Shows the **Pearson correlation coefficient (r)** between every pair of samples, based on their gene expression profiles.

Each cell = pair of samples. The value shown is r (ranges from -1 to +1).

### How to read it

- **r ≈ 1.0** (dark colour) = samples are nearly identical in expression
- **r ≈ 0.9–0.95** = typical biological replicates — good correlation
- **r < 0.8** = unexpected — possible batch effect or sample swap
- Diagonal is always 1.0 (a sample perfectly correlates with itself)

**What to look for:**
- Replicates within the same condition should form a tightly correlated block
- Different conditions should show lower correlation between blocks
- Any sample that correlates poorly with ALL others is an outlier — investigate it

---

## Tab 5 — Heatmap

### What is it?
A heatmap of the **top 50 most variable genes** across all samples.

Each row = one gene. Each column = one sample. Colour = Z-scored expression level (how far above or below average that gene is in that sample).

### How to read it

**Colour scale:**
- Red = higher expression than average (for that gene)
- Blue = lower expression than average
- White = near average

**What to look for:**
- Columns that cluster together have similar expression profiles
- Rows that cluster together are co-regulated (turn on/off together)
- A block of red genes in one condition and blue in another = a regulated gene module

Genes are Z-scored per row so you can compare patterns regardless of absolute expression level.

---

## Tab 6 — Pathway Enrichment

### What is it?
This tab answers: **"Which biological pathways are most affected in my experiment?"**

Genes work together in coordinated groups called pathways (e.g. Cell Cycle, DNA Repair, Apoptosis). If many of your significant DE genes belong to one pathway, that pathway is likely activated or suppressed in your experiment.

### How the enrichment is calculated

For each pathway:
1. Count how many of your significant genes are in that pathway
2. Calculate a **fold enrichment** = how many more hits you got vs. random chance
3. Separate UP-regulated hits (red bars) from DOWN-regulated hits (blue bars)

The chart shows the top 20 pathways sorted by their maximum enrichment score.

### How to read it

**X-axis** — pathways (KEGG / Reactome)
**Y-axis** — fold enrichment (higher = more enriched than expected by chance)

**Red bars** = upregulated genes enriched in that pathway
**Blue bars** = downregulated genes enriched in that pathway

**What to look for:**
- A pathway with tall red AND blue bars = mixed regulation (complex response)
- A pathway with only red bars = pathway being activated
- A pathway with only blue bars = pathway being suppressed

A summary table below the chart lists hit counts and scores for each pathway.

---

## Data Format — What to Upload

The upload accepts **DESeq2 or edgeR differential expression results** as CSV or TSV. Column names are flexible — common variants are automatically detected.

| Column | Accepted names | Example values |
|--------|---------------|----------------|
| Gene symbol | `gene`, `gene_id`, `geneId`, `symbol` | `TP53`, `BRCA1` |
| Log₂ fold change | `log2FoldChange`, `log2fc`, `lfc` | `2.41`, `-1.82` |
| Raw p-value | `pvalue`, `pval`, `p.value` | `0.000032`, `3.2e-10` |
| Adjusted p-value | `padj`, `p.adj`, `fdr` | `0.0012`, `4.0e-08` |
| Base mean (optional) | `baseMean`, `mean_expr`, `avgExpr` | `340.5`, `128.7` |

**Example file:**
```
gene,baseMean,log2FoldChange,pvalue,padj
TP53,340.5,2.41,0.000032,0.0012
BRCA1,128.7,-1.82,0.00036,0.014
MYC,892.1,3.10,0.0000000012,0.000000040
GAPDH,4521.0,0.02,0.91,0.99
```

A sample test file (`test_data_1000genes.csv`) is included in this folder.

**Notes:**
- First row must be a header row with column names
- Missing values (NA, NaN) in `padj` are handled automatically
- No row limit, but files over 100,000 genes may be slow to render

---

## How Is This Data Obtained?

The data visualised here is produced **before** you use this dashboard through a multi-step pipeline:

1. **Wet lab** — grow cells/tissue under two conditions, extract RNA
2. **RNA sequencing** — sequencer counts RNA fragments per gene → FASTQ files
3. **Preprocessing** — align sequences to the genome, count per gene → count matrix
4. **DESeq2 / edgeR** (R package) — normalise counts, run statistical tests → results CSV
5. **Upload that CSV here** → this dashboard visualises the results

The count matrix (Step 3) is the large file (50–200 GB per sample) requiring HPC computing. The DESeq2/edgeR results CSV (Step 4) is small (a few MB) and is what you upload here. Steps 4–5 run fine on a laptop.
