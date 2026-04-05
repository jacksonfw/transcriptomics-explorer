"""
Transcriptomics Explorer
Interactive RNA-seq / differential expression dashboard:
  - Volcano Plot    — significant DE genes
  - MA Plot         — mean expression vs fold change
  - PCA             — sample quality check
  - Sample Corr.    — sample-to-sample Pearson correlation
  - Heatmap         — top gene expression patterns
  - Pathway Enrich. — enriched biological pathways

Run:  python app.py
Open: http://localhost:8050
"""

import base64
import io

import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dash_table, dcc, html
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage

# ─────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────
DARK   = "#0f172a"
CARD   = "#1e293b"
ACCENT = "#38bdf8"

PATHWAY_DB = {
    # Cell Proliferation & Cycle
    "Cell Cycle":           ["TP53","RB1","CDKN2A","CDK4","CDK6","E2F1","CCND1","CDC20","PLK1","AURKA"],
    "G1/S Transition":      ["CCND1","CDK4","CDK6","RB1","E2F1","CDKN1A","CDKN2A","CCNE1","CDK2"],
    "DNA Replication":      ["PCNA","MCM2","MCM7","RFC1","POLA1","POLE","RPA1","GINS1","CDC6"],
    # DNA Damage & Repair
    "DNA Repair":           ["BRCA1","BRCA2","ATM","MLH1","MSH2","RAD51","PARP1","ATR","CHEK1","CHEK2"],
    "Homologous Recomb.":   ["BRCA1","BRCA2","RAD51","PALB2","RBBP8","EME1","MUS81","BLM","RPA1"],
    "Mismatch Repair":      ["MLH1","MSH2","MSH6","PMS2","PCNA","RFC1","EXO1","POLD1"],
    # Cell Death
    "Apoptosis":            ["TP53","BCL2","BAX","CASP3","CASP9","APAF1","CYCS","BID","PUMA","NOXA"],
    "p53 Signaling":        ["TP53","MDM2","CDKN1A","BAX","PUMA","FAS","GADD45A","SESN2","ZMAT3"],
    # Signaling
    "PI3K/AKT/mTOR":        ["PIK3CA","PTEN","AKT1","MTOR","TSC1","TSC2","FOXO3","S6K1","RICTOR"],
    "RAS/MAPK":             ["KRAS","BRAF","RAF1","MAP2K1","MAPK1","EGFR","GRB2","SOS1","HRAS","NRAS"],
    "WNT/Beta-catenin":     ["APC","CTNNB1","GSK3B","AXIN1","TCF7L2","LRP5","FZD1","DVL1","MYC"],
    "Notch Signaling":      ["NOTCH1","DLL3","JAG1","HES1","RBPJ","MAML1","NUMB","LFNG","PSEN1"],
    "Hedgehog Signaling":   ["SHH","PTCH1","SMO","GLI1","GLI2","SUFU","HHIP","BOC","GAS1"],
    "JAK/STAT":             ["JAK1","JAK2","STAT3","STAT5A","STAT1","IL6","IL2","IFNG","SOCS1","SOCS3"],
    # Immune & Inflammation
    "NF-kB Signaling":      ["NFKB1","RELA","IKBKB","NFKBIA","TNF","IL1B","CXCL8","TRAF2","MAP3K7"],
    "Cytokine Signaling":   ["TNF","IL6","IL1B","CXCL8","CXCL10","CCL2","IFNG","IL10","IL4","IL12A"],
    # Metabolism & Hypoxia
    "Angiogenesis":         ["VEGFA","KDR","FLT1","HIF1A","VHL","PDGFRA","NRP1","ANG","TEK","ANGPT1"],
    "HIF1A / Hypoxia":      ["HIF1A","VHL","EPAS1","ARNT","LDHA","SLC2A1","VEGFA","EPO","BNIP3","CA9"],
    "Oxidative Phosph.":    ["ATP5F1A","NDUFS1","SDHB","UQCRC1","COX5A","TFAM","CYCS","SOD2","CAT"],
    # Structural
    "EMT":                  ["CDH1","VIM","SNAI1","SNAI2","TWIST1","ZEB1","FN1","MMP2","MMP9","TGFB1"],
    "Focal Adhesion":       ["ITGA5","ITGB1","PTK2","SRC","VCL","TLN1","PXN","BCAR1","LIMS1"],
    # Protein Homeostasis
    "Ubiquitin Proteasome": ["PSMD1","UBA1","UBE2C","FBXW7","SKP2","BTRC","VHL","MDM2","HERC2"],
}


# ─────────────────────────────────────────────
# DEMO DATA
# ─────────────────────────────────────────────

def make_demo_de(n=800, seed=42):
    rng = np.random.default_rng(seed)
    known = ["TP53","BRCA1","BRCA2","EGFR","MYC","VEGFA","TNF","IL6","CDKN2A","PTEN",
             "RB1","APC","KRAS","BRAF","PIK3CA","ERBB2","CDH1","VHL","MLH1","ATM",
             "STAT3","JAK2","BCL2","BAX","CASP3","MTOR","AKT1","FOXO3","CDK4","CDK6",
             "CCND1","E2F1","RAD51","PARP1","MSH2","APAF1","CYCS","TSC1","TSC2","NRP1",
             "KDR","FLT1","HIF1A","PDGFRA","CXCL8","NFKB1","IL1B","RAF1","GRB2","MEK1"]
    generic   = [f"GENE_{i:04d}" for i in range(1, n - len(known) + 1)]
    genes     = (known + generic)[:n]
    log2fc    = rng.normal(0, 1.2, n)
    base_mean = rng.lognormal(5, 1.5, n)
    noise     = rng.exponential(1, n)
    pval_raw  = stats.chi2.sf((log2fc**2 / 0.8 + noise * 0.3), df=1)
    rank      = np.argsort(np.argsort(pval_raw)) + 1
    pval_adj  = np.clip(pval_raw * n / rank, 0, 1)
    df = pd.DataFrame({
        "gene": genes, "log2FoldChange": log2fc, "pvalue": pval_raw,
        "padj": pval_adj, "baseMean": base_mean,
    })
    return _annotate_de(df)


def make_demo_counts(de_df, seed=7):
    rng     = np.random.default_rng(seed)
    samples = [f"Ctrl_{i}" for i in range(1,4)] + [f"Treat_{i}" for i in range(1,4)]
    genes   = de_df["gene"].values
    base    = de_df["baseMean"].values
    lfc     = de_df["log2FoldChange"].values
    mat     = np.zeros((len(genes), len(samples)))
    for j, s in enumerate(samples):
        fc    = (2 ** lfc) if "Treat" in s else np.ones(len(genes))
        mat[:,j] = rng.negative_binomial(
            np.clip(base * 0.3, 1, None).astype(int),
            0.3 * np.ones(len(genes)),
        ) * fc
    return pd.DataFrame(np.log2(mat + 1), index=genes, columns=samples)


# ─────────────────────────────────────────────
# DATA PROCESSING
# ─────────────────────────────────────────────

def _annotate_de(df):
    df = df.copy()
    for col in ["padj","pvalue","log2FoldChange"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["padj"]           = df["padj"].fillna(1.0)
    df["pvalue"]         = df["pvalue"].fillna(1.0)
    df["log2FoldChange"] = df["log2FoldChange"].fillna(0.0)
    df["baseMean"]       = pd.to_numeric(
        df.get("baseMean", pd.Series(np.ones(len(df)) * 100)),
        errors="coerce").fillna(100.0)
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["significant"] = (df["padj"] < 0.05) & (np.abs(df["log2FoldChange"]) > 1)
    df["regulation"]  = "NS"
    df.loc[(df["padj"]<0.05)&(df["log2FoldChange"]> 1),"regulation"] = "Up"
    df.loc[(df["padj"]<0.05)&(df["log2FoldChange"]<-1),"regulation"] = "Down"
    return df


def _pathway_enrichment(df):
    n_total    = len(df)
    up_genes   = set(df.loc[df["regulation"]=="Up",   "gene"].str.upper())
    down_genes = set(df.loc[df["regulation"]=="Down", "gene"].str.upper())
    bg_up   = max(len(up_genes)   / n_total, 0.001)
    bg_down = max(len(down_genes) / n_total, 0.001)
    rows = []
    for pathway, genes in PATHWAY_DB.items():
        gs        = set(g.upper() for g in genes)
        n_pw      = len(gs)
        hits_up   = len(up_genes   & gs)
        hits_down = len(down_genes & gs)
        score_up   = round((hits_up   / n_pw) / bg_up,   3) if hits_up   > 0 else 0.0
        score_down = round((hits_down / n_pw) / bg_down, 3) if hits_down > 0 else 0.0
        rows.append({
            "pathway": pathway, "n_genes": n_pw,
            "hits_up": hits_up,   "score_up":   score_up,
            "hits_down": hits_down, "score_down": score_down,
            "max_score": max(score_up, score_down),
        })
    return pd.DataFrame(rows).sort_values("max_score", ascending=False)


# ─────────────────────────────────────────────
# DEMO GLOBALS
# ─────────────────────────────────────────────
DEMO_DE     = make_demo_de()
DEMO_COUNTS = make_demo_counts(DEMO_DE)


# ─────────────────────────────────────────────
# PLOT BUILDERS
# ─────────────────────────────────────────────

def build_volcano(df, fc_thresh=1.0, pval_thresh=0.05, highlight_genes=None):
    df = df.copy()
    df["regulation"] = "NS"
    df.loc[(df["padj"]<pval_thresh)&(df["log2FoldChange"]> fc_thresh),"regulation"] = "Up"
    df.loc[(df["padj"]<pval_thresh)&(df["log2FoldChange"]<-fc_thresh),"regulation"] = "Down"
    cmap = {"Up":"#e63946","Down":"#457b9d","NS":"#adb5bd"}
    smap = {"Up":7,"Down":7,"NS":4}
    omap = {"Up":0.85,"Down":0.85,"NS":0.35}
    fig  = go.Figure()
    for reg in ["NS","Up","Down"]:
        sub = df[df["regulation"]==reg]
        fig.add_trace(go.Scatter(
            x=sub["log2FoldChange"], y=sub["-log10padj"], mode="markers",
            name=f"{reg} ({len(sub)})", text=sub["gene"],
            hovertemplate="<b>%{text}</b><br>log2FC: %{x:.3f}<br>-log10(padj): %{y:.3f}<extra></extra>",
            marker=dict(color=cmap[reg], size=smap[reg], opacity=omap[reg], line=dict(width=0)),
        ))
    if highlight_genes:
        hl = df[df["gene"].isin(highlight_genes)]
        if not hl.empty:
            fig.add_trace(go.Scatter(
                x=hl["log2FoldChange"], y=hl["-log10padj"],
                mode="markers+text", text=hl["gene"], textposition="top center",
                textfont=dict(color="#fbbf24", size=10),
                marker=dict(color="#fbbf24", size=12, symbol="star",
                            line=dict(width=1, color="#fff")),
                name="Highlighted", hoverinfo="skip",
            ))
    fig.add_hline(y=-np.log10(pval_thresh), line_dash="dash", line_color="#6c757d",
                  line_width=1, annotation_text=f"padj={pval_thresh}", annotation_position="right")
    fig.add_vline(x= fc_thresh, line_dash="dash", line_color="#6c757d", line_width=1)
    fig.add_vline(x=-fc_thresh, line_dash="dash", line_color="#6c757d", line_width=1)
    fig.update_layout(
        title=dict(text="Volcano Plot — Differential Expression", font=dict(size=16)),
        xaxis_title="log₂ Fold Change", yaxis_title="-log₁₀ (adjusted p-value)",
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(size=11)), hovermode="closest",
        xaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        yaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        margin=dict(l=60,r=20,t=50,b=50),
    )
    return fig


def build_ma(df):
    df  = df.copy()
    df["log10mean"] = np.log10(df["baseMean"].clip(lower=0.01))
    cmap = {"Up":"#e63946","Down":"#457b9d","NS":"#adb5bd"}
    smap = {"Up":7,"Down":7,"NS":4}
    omap = {"Up":0.85,"Down":0.85,"NS":0.25}
    fig  = go.Figure()
    for reg in ["NS","Up","Down"]:
        sub = df[df["regulation"]==reg]
        fig.add_trace(go.Scatter(
            x=sub["log10mean"], y=sub["log2FoldChange"], mode="markers",
            name=f"{reg} ({len(sub)})", text=sub["gene"],
            hovertemplate="<b>%{text}</b><br>log₁₀(mean): %{x:.3f}<br>log₂FC: %{y:.3f}<extra></extra>",
            marker=dict(color=cmap[reg], size=smap[reg], opacity=omap[reg], line=dict(width=0)),
        ))
    fig.add_hline(y=0,  line_color="#6c757d", line_width=1)
    fig.add_hline(y= 1, line_dash="dash", line_color="#6c757d", line_width=1,
                  annotation_text="FC=2", annotation_position="right")
    fig.add_hline(y=-1, line_dash="dash", line_color="#6c757d", line_width=1)
    fig.update_layout(
        title=dict(text="MA Plot — Mean Expression vs Fold Change", font=dict(size=16)),
        xaxis_title="log₁₀ (mean expression)", yaxis_title="log₂ Fold Change",
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(size=11)), hovermode="closest",
        xaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        yaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        margin=dict(l=60,r=20,t=50,b=50),
    )
    return fig


def build_pca(counts_df):
    mat      = counts_df.values.T
    mat_c    = mat - mat.mean(axis=0)
    U, S, _  = np.linalg.svd(mat_c, full_matrices=False)
    pcs      = U * S
    var      = (S**2) / (S**2).sum() * 100
    samples  = counts_df.columns.tolist()
    groups   = ["Control" if "Ctrl" in s else "Treated" for s in samples]
    fig = go.Figure()
    for grp, color in [("Control","#457b9d"),("Treated","#e63946")]:
        idx = [i for i,g in enumerate(groups) if g==grp]
        fig.add_trace(go.Scatter(
            x=pcs[idx,0], y=pcs[idx,1], mode="markers+text",
            text=[samples[i] for i in idx], textposition="top center",
            textfont=dict(size=9, color="#e2e8f0"),
            marker=dict(size=14, color=color, opacity=0.9,
                        line=dict(width=1.5, color="#fff")),
            name=grp,
            hovertemplate="<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
        ))
    fig.update_layout(
        title=dict(text="PCA — Sample Quality Check", font=dict(size=16)),
        xaxis_title=f"PC1 ({var[0]:.1f}% variance)",
        yaxis_title=f"PC2 ({var[1]:.1f}% variance)",
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
        xaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        yaxis=dict(gridcolor="#1e293b", zerolinecolor="#334155"),
        margin=dict(l=60,r=20,t=50,b=50),
    )
    return fig, float(var[0]), float(var[1])


def build_sample_corr(counts_df):
    mat     = counts_df.values.T        # samples × genes
    corr    = np.corrcoef(mat)          # samples × samples
    samples = counts_df.columns.tolist()
    annot   = np.round(corr, 2)
    fig = go.Figure(go.Heatmap(
        z=corr, x=samples, y=samples,
        colorscale="RdBu_r", zmin=-1, zmax=1, zmid=0,
        colorbar=dict(title="Pearson r", tickfont=dict(color="#e2e8f0"),
                      titlefont=dict(color="#e2e8f0")),
        text=annot, texttemplate="%{text}", textfont=dict(size=11),
        hovertemplate="<b>%{x}</b> vs <b>%{y}</b><br>Pearson r: %{z:.4f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text="Sample Correlation — Pearson r (log₂ counts)", font=dict(size=16)),
        xaxis=dict(tickangle=-35, tickfont=dict(size=10), side="bottom"),
        yaxis=dict(tickfont=dict(size=10), autorange="reversed"),
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        margin=dict(l=120,r=20,t=50,b=100), height=520,
    )
    return fig


def build_heatmap(de_df, counts_df, n_top=40):
    top_genes = de_df.nsmallest(n_top,"padj")["gene"].tolist()
    available = [g for g in top_genes if g in counts_df.index]
    if not available:
        available = counts_df.index[:n_top].tolist()
    sub      = counts_df.loc[available]
    row_mean = sub.mean(axis=1)
    row_std  = sub.std(axis=1).replace(0, 1)
    z        = sub.subtract(row_mean, axis=0).divide(row_std, axis=0)
    if len(z) > 2:
        order = leaves_list(linkage(z.values, method="ward"))
        z     = z.iloc[order]
    fig = go.Figure(go.Heatmap(
        z=z.values, x=z.columns.tolist(), y=z.index.tolist(),
        colorscale="RdBu_r", zmid=0,
        colorbar=dict(title="Z-score", tickfont=dict(color="#e2e8f0"),
                      titlefont=dict(color="#e2e8f0")),
        hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=f"Heatmap — Top {len(z)} Significant Genes (Z-scored)", font=dict(size=16)),
        xaxis=dict(tickangle=-35, tickfont=dict(size=10)),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        margin=dict(l=120,r=20,t=50,b=80), height=600,
    )
    return fig


def build_pathway_bar(results):
    top = results.head(20).sort_values("max_score")  # ascending so top is at top of chart
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="Upregulated genes", y=top["pathway"], x=top["score_up"],
        orientation="h", marker_color="#e63946", opacity=0.85,
        customdata=np.column_stack([top["hits_up"], top["n_genes"]]),
        hovertemplate="<b>%{y}</b><br>Enrichment: %{x:.2f}×<br>Hits: %{customdata[0]}/%{customdata[1]} genes<extra></extra>",
    ))
    fig.add_trace(go.Bar(
        name="Downregulated genes", y=top["pathway"], x=top["score_down"],
        orientation="h", marker_color="#457b9d", opacity=0.85,
        customdata=np.column_stack([top["hits_down"], top["n_genes"]]),
        hovertemplate="<b>%{y}</b><br>Enrichment: %{x:.2f}×<br>Hits: %{customdata[0]}/%{customdata[1]} genes<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text="Pathway Enrichment — Fold Enrichment over Background", font=dict(size=16)),
        xaxis=dict(title="Fold Enrichment (vs background rate)", gridcolor="#1e293b"),
        yaxis=dict(tickfont=dict(size=10)),
        barmode="group",
        plot_bgcolor=DARK, paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#e2e8f0"),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
        height=620, margin=dict(l=190,r=40,t=50,b=60),
    )
    return fig


def _pathway_table(results):
    disp = results[["pathway","n_genes","hits_up","score_up","hits_down","score_down"]].copy()
    disp.columns = ["Pathway","Pathway Genes","Up Hits","Up Score","Down Hits","Down Score"]
    disp["Up Score"]   = disp["Up Score"].round(2)
    disp["Down Score"] = disp["Down Score"].round(2)
    return dash_table.DataTable(
        data=disp.to_dict("records"),
        columns=[{"name":c,"id":c} for c in disp.columns],
        sort_action="native", page_size=10,
        style_table={"overflowX":"auto","marginTop":"12px"},
        style_cell={"background":CARD,"color":"#e2e8f0","border":"1px solid #334155",
                    "fontSize":"0.8rem","padding":"6px 10px","textAlign":"left"},
        style_header={"background":"#0f172a","fontWeight":"700","color":ACCENT,
                      "border":"1px solid #334155"},
        style_data_conditional=[
            {"if":{"row_index":"odd"},"background":"#162032"},
            {"if":{"filter_query":"{Up Hits} > 0","column_id":"Up Hits"},
             "color":"#e63946","fontWeight":"700"},
            {"if":{"filter_query":"{Down Hits} > 0","column_id":"Down Hits"},
             "color":"#457b9d","fontWeight":"700"},
        ],
    )


# ─────────────────────────────────────────────
# APP
# ─────────────────────────────────────────────
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.SLATE, dbc.icons.BOOTSTRAP],
    suppress_callback_exceptions=True,
    title="Transcriptomics Explorer",
)


def stat_card(icon, value, label, color):
    return dbc.Col(dbc.Card([dbc.CardBody([
        html.Div([
            html.I(className=f"bi bi-{icon} me-2", style={"fontSize":"1.4rem","color":color}),
            html.Span(value, style={"fontSize":"1.6rem","fontWeight":"700","color":color}),
        ], className="d-flex align-items-center"),
        html.P(label, className="mb-0 text-muted small mt-1"),
    ])], style={"background":CARD,"border":f"1px solid {color}22"}), xs=6, md=3)


app.layout = dbc.Container(fluid=True, style={"background":DARK,"minHeight":"100vh"}, children=[

    # HEADER
    dbc.Row(dbc.Col(html.Div([
        html.Div([
            html.I(className="bi bi-activity me-3", style={"fontSize":"2rem","color":ACCENT}),
            html.Div([
                html.H2("Transcriptomics Explorer", className="mb-0",
                        style={"color":ACCENT,"fontWeight":"700"}),
                html.P("Volcano · MA Plot · PCA · Sample Correlation · Heatmap · Pathway Enrichment",
                       className="mb-0 text-muted small"),
            ]),
        ], className="d-flex align-items-center"),
    ], className="py-3 px-2")), className="border-bottom mb-3"),

    # UPLOAD + STATS
    dbc.Row([
        dbc.Col([
            dcc.Upload(id="upload-data", multiple=False,
                children=html.Div([
                    html.I(className="bi bi-cloud-upload me-2",
                           style={"fontSize":"1.3rem","color":ACCENT}),
                    html.Span("Drop DESeq2 / edgeR CSV or TSV or "),
                    html.A("Browse", style={"color":ACCENT,"cursor":"pointer"}),
                ]),
                style={"borderWidth":"1px","borderStyle":"dashed","borderRadius":"8px",
                       "borderColor":"#334155","textAlign":"center","padding":"14px",
                       "background":CARD,"color":"#94a3b8","cursor":"pointer"},
            ),
            html.Div("Currently showing: 800 demo genes. Upload your own DESeq2/edgeR results to replace.",
                     className="text-muted small mt-1", id="upload-status"),
        ], md=4),
        dbc.Col(dbc.Row(id="stat-cards"), md=8),
    ], className="mb-3 g-3"),

    # TABS
    dbc.Tabs(id="main-tabs", active_tab="tab-volcano", children=[
        dbc.Tab(label="Volcano Plot",        tab_id="tab-volcano",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
        dbc.Tab(label="MA Plot",             tab_id="tab-ma",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
        dbc.Tab(label="PCA",                 tab_id="tab-pca",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
        dbc.Tab(label="Sample Correlation",  tab_id="tab-correlation",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
        dbc.Tab(label="Heatmap",             tab_id="tab-heatmap",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
        dbc.Tab(label="Pathway Enrichment",  tab_id="tab-pathway",
                label_style={"color":"#94a3b8"}, active_label_style={"color":ACCENT}),
    ], style={"borderBottom":"1px solid #1e293b"}),

    html.Div(id="tab-content", className="mt-3"),
    dcc.Store(id="store-de",     data=DEMO_DE.to_json(orient="split")),
    dcc.Store(id="store-counts", data=DEMO_COUNTS.to_json(orient="split")),
])


# ─────────────────────────────────────────────
# TAB RENDERER
# ─────────────────────────────────────────────

@app.callback(Output("tab-content","children"),
              Input("main-tabs","active_tab"),
              Input("store-de","data"),
              Input("store-counts","data"))
def render_tab(tab, de_json, counts_json):
    df     = pd.read_json(io.StringIO(de_json),     orient="split")
    counts = pd.read_json(io.StringIO(counts_json), orient="split")

    if tab == "tab-volcano":
        return dbc.Row([
            dbc.Col(dbc.Card(dbc.CardBody([
                html.H6("Filters", className="text-uppercase mb-3",
                        style={"color":ACCENT,"letterSpacing":"0.08em","fontSize":"0.75rem"}),
                html.Label("log₂FC threshold", className="small text-muted"),
                dcc.Slider(id="fc-slider", min=0, max=4, step=0.1, value=1,
                           marks={i:str(i) for i in range(5)}, tooltip={"placement":"bottom"}),
                html.Label("Adjusted p-value", className="small text-muted mt-3"),
                dcc.Dropdown(id="pval-drop",
                             options=[{"label":f"p < {v}","value":v}
                                      for v in [0.001,0.01,0.05,0.1]],
                             value=0.05, clearable=False,
                             style={"background":"#0f172a","color":"#e2e8f0"}),
                html.Label("Gene search (highlight)", className="small text-muted mt-3"),
                dcc.Input(id="gene-search", type="text", placeholder="e.g. TP53, BRCA1",
                          className="form-control form-control-sm",
                          style={"background":"#0f172a","color":"#e2e8f0","border":"1px solid #334155"}),
                html.Hr(style={"borderColor":"#334155"}),
                html.H6("Summary", className="text-uppercase mb-2",
                        style={"color":ACCENT,"letterSpacing":"0.08em","fontSize":"0.75rem"}),
                html.Div(id="volcano-summary"),
            ]), style={"background":CARD}), md=2),
            dbc.Col([
                dcc.Graph(id="volcano-plot", style={"height":"460px"},
                          config={"displayModeBar":True,"modeBarButtonsToRemove":["lasso2d"]}),
                dbc.Row([
                    dbc.Col(html.H6("Selected / Significant Genes",
                                    className="mt-1 mb-0 text-muted small"), width="auto"),
                    dbc.Col(dbc.Button(
                        [html.I(className="bi bi-download me-1"), "Export CSV"],
                        id="export-btn", size="sm", color="secondary", outline=True,
                    ), width="auto"),
                ], className="mt-3 mb-2 align-items-center"),
                dcc.Download(id="export-download"),
                html.Div(id="de-table-container"),
            ], md=10),
        ], className="g-3")

    elif tab == "tab-ma":
        n_up = int((df["regulation"]=="Up").sum())
        n_dn = int((df["regulation"]=="Down").sum())
        return dbc.Row([dbc.Col([
            dcc.Graph(figure=build_ma(df), style={"height":"480px"},
                      config={"displayModeBar":True}),
            dbc.Row([
                dbc.Col(html.Div([html.Span("Upregulated: ", className="text-muted small"),
                    html.Span(str(n_up), style={"color":"#e63946","fontWeight":"700"})]), md=3),
                dbc.Col(html.Div([html.Span("Downregulated: ", className="text-muted small"),
                    html.Span(str(n_dn), style={"color":"#457b9d","fontWeight":"700"})]), md=3),
            ], className="mt-2 px-2"),
            dbc.Alert([html.I(className="bi bi-info-circle me-2"),
                "MA plot shows mean expression (X) vs fold change (Y). "
                "Genes with high expression and large fold change (top/bottom right) "
                "are the most reliable hits. Low-expression genes with high fold change "
                "(top/bottom left) may be noise."],
                color="secondary", className="mt-3 small py-2"),
        ])])

    elif tab == "tab-pca":
        fig, var1, var2 = build_pca(counts)
        n_ctrl    = sum("Ctrl"  in c for c in counts.columns)
        n_treated = sum("Treat" in c for c in counts.columns)
        return dbc.Row([dbc.Col([
            dcc.Graph(figure=fig, style={"height":"480px"}, config={"displayModeBar":True}),
            dbc.Row([
                dbc.Col(html.Div([html.Span("PC1 variance: ", className="text-muted small"),
                    html.Span(f"{var1:.1f}%", style={"color":ACCENT,"fontWeight":"700"})]), md=3),
                dbc.Col(html.Div([html.Span("PC2 variance: ", className="text-muted small"),
                    html.Span(f"{var2:.1f}%", style={"color":ACCENT,"fontWeight":"700"})]), md=3),
                dbc.Col(html.Div([html.Span("Control: ", className="text-muted small"),
                    html.Span(str(n_ctrl), style={"color":"#457b9d","fontWeight":"700"})]), md=3),
                dbc.Col(html.Div([html.Span("Treated: ", className="text-muted small"),
                    html.Span(str(n_treated), style={"color":"#e63946","fontWeight":"700"})]), md=3),
            ], className="mt-2 px-2"),
            dbc.Alert([html.I(className="bi bi-info-circle me-2"),
                "Good QC: Control samples cluster together (blue), Treated samples cluster "
                "together (red), with clear separation between groups. Outlier samples "
                "that don't cluster with their group should be investigated. "
                "PCA here is simulated from DE results — upload a raw count matrix for real PCA."],
                color="secondary", className="mt-3 small py-2"),
        ])])

    elif tab == "tab-correlation":
        return dbc.Row([dbc.Col([
            dcc.Graph(figure=build_sample_corr(counts), style={"height":"540px"},
                      config={"displayModeBar":True}),
            dbc.Alert([html.I(className="bi bi-info-circle me-2"),
                "Values close to 1.0 = highly similar samples. "
                "Replicates within the same group (e.g. Ctrl_1 vs Ctrl_2) should have "
                "r > 0.95. Low correlation within a group may indicate a failed sample. "
                "Cross-group correlation (Control vs Treated) is expected to be lower."],
                color="secondary", className="mt-2 small py-2"),
        ])])

    elif tab == "tab-heatmap":
        n_top = min(40, max(int(df["significant"].sum()), 20))
        return dbc.Row([dbc.Col([
            dcc.Graph(figure=build_heatmap(df, counts, n_top=n_top),
                      style={"height":"650px"}, config={"displayModeBar":True}),
            dbc.Alert([html.I(className="bi bi-info-circle me-2"),
                "Top significant genes by adjusted p-value. "
                "Red = high expression, Blue = low (Z-scored per gene). "
                "Genes are clustered by expression similarity — "
                "genes in the same cluster behave similarly across samples."],
                color="secondary", className="mt-2 small py-2"),
        ])])

    elif tab == "tab-pathway":
        results   = _pathway_enrichment(df)
        sig_count = int(df["significant"].sum())
        n_enriched = int((results["max_score"] > 0).sum())
        return dbc.Row([dbc.Col([
            dcc.Graph(figure=build_pathway_bar(results), style={"height":"650px"},
                      config={"displayModeBar":True}),
            dbc.Row([
                dbc.Col(html.Div([html.Span("Significant genes tested: ", className="text-muted small"),
                    html.Span(str(sig_count), style={"color":ACCENT,"fontWeight":"700"})]), md=4),
                dbc.Col(html.Div([html.Span("Pathways with hits: ", className="text-muted small"),
                    html.Span(str(n_enriched), style={"color":"#2ecc71","fontWeight":"700"})]), md=4),
            ], className="mt-2 px-2"),
            html.H6("Enrichment Summary", className="mt-3 mb-1 text-muted small"),
            _pathway_table(results),
            dbc.Alert([html.I(className="bi bi-info-circle me-2"),
                "Enrichment score = fraction of pathway genes that are significant in your data "
                "divided by the background rate. Score > 1 means the pathway is enriched. "
                "Uses 22 representative KEGG/Reactome pathways. "
                "For full pathway analysis, run clusterProfiler in R after DESeq2."],
                color="secondary", className="mt-2 small py-2"),
        ])])


# ─────────────────────────────────────────────
# CALLBACKS
# ─────────────────────────────────────────────

@app.callback(Output("stat-cards","children"), Input("store-de","data"))
def update_stats(de_json):
    df = pd.read_json(io.StringIO(de_json), orient="split")
    return [
        stat_card("arrow-up-circle",   int((df["regulation"]=="Up").sum()),   "Upregulated",         "#e63946"),
        stat_card("arrow-down-circle", int((df["regulation"]=="Down").sum()), "Downregulated",       "#457b9d"),
        stat_card("check2-circle",     int(df["significant"].sum()),          "Significant (p<0.05)","#2ecc71"),
        stat_card("bar-chart",         f"{len(df):,}",                        "Total Genes",         ACCENT),
    ]


@app.callback(
    Output("store-de","data"),
    Output("store-counts","data"),
    Output("upload-status","children"),
    Input("upload-data","contents"),
    State("upload-data","filename"),
    prevent_initial_call=True,
)
def handle_upload(contents, filename):
    if not contents:
        return dash.no_update, dash.no_update, ""
    _, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    try:
        sep = "\t" if (filename or "").endswith(".tsv") else ","
        df  = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep=sep)
        rename = {}
        for col in df.columns:
            lc = col.lower()
            if   lc in ("gene","gene_id","geneid","symbol"):       rename[col]="gene"
            elif lc in ("log2foldchange","log2fc","lfc"):          rename[col]="log2FoldChange"
            elif lc in ("pvalue","pval","p.value","p_value"):      rename[col]="pvalue"
            elif lc in ("padj","p.adj","fdr","adjusted_pvalue"):   rename[col]="padj"
            elif lc in ("basemean","mean_expr","avgexpr"):         rename[col]="baseMean"
        df = df.rename(columns=rename)
        missing = {"gene","log2FoldChange","pvalue","padj"} - set(df.columns)
        if missing:
            return dash.no_update, dash.no_update, dbc.Alert(
                f"Missing columns: {missing}", color="danger", className="py-1")
        df     = _annotate_de(df)
        counts = make_demo_counts(df)
        return df.to_json(orient="split"), counts.to_json(orient="split"), dbc.Alert(
            f"Loaded {filename}: {len(df):,} genes", color="success", className="py-1")
    except Exception as e:
        return dash.no_update, dash.no_update, dbc.Alert(
            f"Error: {e}", color="danger", className="py-1")


@app.callback(
    Output("volcano-plot","figure"),
    Output("volcano-summary","children"),
    Input("fc-slider","value"),
    Input("pval-drop","value"),
    Input("gene-search","value"),
    Input("store-de","data"),
)
def update_volcano(fc, pval, search, de_json):
    df         = pd.read_json(io.StringIO(de_json), orient="split")
    highlights = [g.strip() for g in (search or "").split(",") if g.strip()]
    fig        = build_volcano(df, fc_thresh=fc, pval_thresh=pval, highlight_genes=highlights)
    n_up = int(((df["padj"]<pval)&(df["log2FoldChange"]> fc)).sum())
    n_dn = int(((df["padj"]<pval)&(df["log2FoldChange"]<-fc)).sum())
    summary = [
        html.Div([html.Span("▲ Up: ",   style={"color":"#e63946"}), html.Span(str(n_up))], className="small mb-1"),
        html.Div([html.Span("▼ Down: ", style={"color":"#457b9d"}), html.Span(str(n_dn))], className="small"),
    ]
    return fig, summary


@app.callback(
    Output("de-table-container","children"),
    Input("volcano-plot","selectedData"),
    Input("fc-slider","value"),
    Input("pval-drop","value"),
    Input("store-de","data"),
)
def update_de_table(selected, fc, pval, de_json):
    df = pd.read_json(io.StringIO(de_json), orient="split")
    df["regulation"] = "NS"
    df.loc[(df["padj"]<pval)&(df["log2FoldChange"]> fc),"regulation"] = "Up"
    df.loc[(df["padj"]<pval)&(df["log2FoldChange"]<-fc),"regulation"] = "Down"
    if selected and selected.get("points"):
        sub = df[df["gene"].isin([p["text"] for p in selected["points"] if "text" in p])]
    else:
        sub = df[df["regulation"]!="NS"].nlargest(50,"-log10padj")
    sub = sub[["gene","log2FoldChange","pvalue","padj","baseMean","regulation"]].copy()
    sub["log2FoldChange"] = sub["log2FoldChange"].round(3)
    sub["pvalue"]   = sub["pvalue"].apply(lambda x: f"{x:.2e}")
    sub["padj"]     = sub["padj"].apply(lambda x: f"{x:.2e}")
    sub["baseMean"] = sub["baseMean"].round(1)
    return dash_table.DataTable(
        data=sub.to_dict("records"),
        columns=[{"name":c,"id":c,"type":"text"} for c in sub.columns],
        filter_action="native", sort_action="native", page_size=12,
        style_table={"overflowX":"auto"},
        style_cell={"background":CARD,"color":"#e2e8f0","border":"1px solid #334155",
                    "fontSize":"0.8rem","padding":"6px 10px","textAlign":"left"},
        style_header={"background":"#0f172a","fontWeight":"700","color":ACCENT,
                      "border":"1px solid #334155"},
        style_data_conditional=[
            {"if":{"filter_query":'{regulation} = "Up"',   "column_id":"regulation"},"color":"#e63946","fontWeight":"700"},
            {"if":{"filter_query":'{regulation} = "Down"', "column_id":"regulation"},"color":"#457b9d","fontWeight":"700"},
            {"if":{"row_index":"odd"},"background":"#162032"},
        ],
        style_filter={"background":"#0f172a","color":"#e2e8f0","border":"1px solid #334155"},
    )


@app.callback(
    Output("export-download","data"),
    Input("export-btn","n_clicks"),
    State("store-de","data"),
    State("fc-slider","value"),
    State("pval-drop","value"),
    prevent_initial_call=True,
)
def export_csv(_, de_json, fc, pval):
    df = pd.read_json(io.StringIO(de_json), orient="split")
    df["regulation"] = "NS"
    df.loc[(df["padj"]<pval)&(df["log2FoldChange"]> fc),"regulation"] = "Up"
    df.loc[(df["padj"]<pval)&(df["log2FoldChange"]<-fc),"regulation"] = "Down"
    sig = df[df["regulation"]!="NS"].sort_values("-log10padj", ascending=False)
    return dcc.send_data_frame(
        sig[["gene","log2FoldChange","pvalue","padj","baseMean","regulation"]].to_csv,
        "significant_genes.csv", index=False,
    )


# ─────────────────────────────────────────────
if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 7860))
    app.run(debug=False, host="0.0.0.0", port=port)
