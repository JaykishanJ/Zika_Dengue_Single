#!/usr/bin/env python3
"""
STEP 08 - Publication + supplementary figures at 600 dpi PNG.

This script compiles outputs from Steps 02-07 and generates high-quality 
publication-ready figures and supplementary figures.

Outputs:
    - publication/Figure*.png
    - supplementary/FigureS*.png
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch
from matplotlib_venn import venn3
from venn import venn as venn4

warnings.filterwarnings("ignore")

# Define categorical color palette globally
CAT = {
    "Mock": "#7f8c8d", "DENV": "#c0392b", "ZIKV": "#2471a3",
    "High": "#c0392b", "Low": "#e67e22", "Bystander": "#2980b9"
}

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate publication and supplementary figures.")
    
    # Base paths
    import sys
    _root_dir = Path(__file__).resolve().parent.parent.parent
    if str(_root_dir) not in sys.path: sys.path.append(str(_root_dir))
    from src import config
    default_here = config.PROCESSED_DATA_DIR / Path(__file__).resolve().parent.name
    default_root = config.PROCESSED_DATA_DIR
    
    parser.add_argument("--root-dir", type=Path, default=default_root,
                        help="Root directory of the workspace.")
    parser.add_argument("--out-dir", type=Path, default=default_here,
                        help="Output directory for figures.")
    parser.add_argument("--dpi", type=int, default=600,
                        help="DPI for saved figures (default: 600).")
    
    # Logging
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level.")
    return parser.parse_args()

def setup_logger(level: str) -> logging.Logger:
    """Set up the logger for the script."""
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, level.upper()))
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(ch)
    return logger

def apply_plot_settings(dpi: int) -> None:
    """Configure global matplotlib rcParams for publication standard."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "axes.titleweight": "bold",
        "savefig.dpi": dpi,
        "savefig.bbox": "tight"
    })

def save_fig(fig: plt.Figure, folder: Path, name: str, dpi: int, logger: logging.Logger) -> None:
    """Save a figure and close it."""
    out_path = folder / f"{name}.png"
    try:
        fig.savefig(out_path, dpi=dpi)
        logger.info(f"Wrote {name}")
    except Exception as e:
        logger.error(f"Failed to save {name}: {e}")
    finally:
        plt.close(fig)

def ucat(ax: plt.Axes, adata: sc.AnnData, col: str, title: str) -> None:
    """Categorical scatter plot on UMAP space."""
    um = adata.obsm["X_umap"]
    v = adata.obs[col].astype(str)
    
    preset = [c for c in ["Mock", "Bystander", "Low", "High", "DENV", "ZIKV"] if c in set(v)]
    if preset:
        cats = preset
        cm = None
    else:
        cats = sorted(v.unique(), key=lambda x: (len(x), x))
        cm = plt.get_cmap("tab20")
        
    for i, c in enumerate(cats):
        m = (v == c).to_numpy()
        color = [CAT.get(c) if cm is None else cm(i % 20)]
        ax.scatter(um[m, 0], um[m, 1], s=3, c=color, label=c, linewidths=0, alpha=0.8)
        
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.legend(markerscale=3, fontsize=7, frameon=False, ncol=2 if cm else 1)

def ucont(ax: plt.Axes, adata: sc.AnnData, vals: np.ndarray, title: str, cmap: str = "magma") -> None:
    """Continuous scatter plot on UMAP space."""
    um = adata.obsm["X_umap"]
    o = np.argsort(vals)
    s = ax.scatter(um[o, 0], um[o, 1], s=3, c=np.asarray(vals)[o], cmap=cmap, linewidths=0)
    
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    plt.colorbar(s, ax=ax, shrink=0.7)

# --- Individual Figure Generation Functions ---

def generate_figure_1(adata: sc.AnnData, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig1: Single-cell landscape UMAPs."""
    logger.info("Generating Figure 1...")
    fig, ax = plt.subplots(2, 2, figsize=(10, 9))
    
    ucat(ax[0, 0], adata, "condition", "Condition")
    ucat(ax[0, 1], adata, "infection_state", "Infection state")
    ucont(ax[1, 0], adata, adata.obs["log_viral_load"].to_numpy(), "log viral load")
    
    if "leiden" in adata.obs:
        ucat(ax[1, 1], adata, "leiden", "Leiden clusters")
    else:
        ax[1, 1].axis("off")
        
    fig.suptitle("Single-cell landscape (GSE110496, Huh7)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure1_scRNA_landscape", dpi, logger)

def generate_figure_2(s3_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig2: Single-cell differential expression volcano plots."""
    logger.info("Generating Figure 2...")
    fdr, lfc = 0.05, 0.58
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    
    for a, v in zip(ax, ["DENV", "ZIKV"]):
        path = s3_dir / f"DE_High_vs_Mock_{v}.csv"
        if not path.exists():
            continue
            
        r = pd.read_csv(path).dropna(subset=["padj"])
        r["nlp"] = -np.log10(r["padj"].clip(lower=1e-300))
        
        up = (r["padj"] < fdr) & (r["log2FoldChange"] > lfc)
        dn = (r["padj"] < fdr) & (r["log2FoldChange"] < -lfc)
        
        a.scatter(r["log2FoldChange"], r["nlp"], s=5, c="#d5d8dc")
        a.scatter(r["log2FoldChange"][up], r["nlp"][up], s=8, c="#c0392b", label=f"up ({int(up.sum())})")
        a.scatter(r["log2FoldChange"][dn], r["nlp"][dn], s=8, c="#2471a3", label=f"down ({int(dn.sum())})")
        
        a.axvline(lfc, ls="--", c="k", lw=0.5)
        a.axvline(-lfc, ls="--", c="k", lw=0.5)
        a.set(xlabel="log2FC (High vs Mock)", ylabel="-log10 adj P", title=v)
        a.legend(fontsize=8)
        
    fig.suptitle("Single-cell differential expression", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure2_scRNA_volcanoes", dpi, logger)

def generate_figure_3(s4_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig3: Core signature heatmap."""
    logger.info("Generating Figure 3...")
    path = s4_dir / "convergent_core_signature.csv"
    if not path.exists(): return
    
    c = pd.read_csv(path)
    hm = c[c["concordant"] == True].set_index("symbol")[["DENV_log2FC", "ZIKV_log2FC"]].sort_values("DENV_log2FC")
    
    # Restrict to top 25 up and top 25 down genes for publication readability
    if len(hm) > 50:
        hm = pd.concat([hm.head(25), hm.tail(25)])
    
    fig, a = plt.subplots(figsize=(4, 8))
    im = a.imshow(hm.values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
    
    a.set_xticks([0, 1])
    a.set_xticklabels(["DENV", "ZIKV"], fontsize=12)
    
    a.set_yticks(range(len(hm)))
    a.set_yticklabels(hm.index, fontsize=9)
        
    a.set_title("Convergent core (Top 50)", fontsize=14, pad=10, fontweight="bold")
    
    # Adjust colorbar to be nicely proportioned
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(a)
    cax = divider.append_axes("right", size="8%", pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("log2FC", size=10)
    
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure3_convergent_core_heatmap", dpi, logger)

def generate_figure_4(adata: sc.AnnData, s5_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig4: Progression trends for top genes."""
    logger.info("Generating Figure 4...")
    path = s5_dir / "convergent_progression_signature.csv"
    if not path.exists(): return
    
    prog = pd.read_csv(path).head(9)
    symbols = adata.var["symbol"].astype(str)
    s2i = {s.upper(): i for i, s in enumerate(symbols)}
    
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    axes = axes.ravel()
    
    for ax_, (_, r) in zip(axes, prog.iterrows()):
        gid = s2i.get(str(r["symbol"]).upper())
        if gid is None: continue
        
        for v, c in [("DENV", "#c0392b"), ("ZIKV", "#2471a3")]:
            sub = adata[(adata.obs["condition"] == v).to_numpy()]
            x = sub.obs["log_viral_load"].to_numpy()
            
            e = sub[:, gid].X
            e = np.asarray(e.todense()).ravel() if hasattr(e, "todense") else np.asarray(e).ravel()
            
            o = np.argsort(x)
            d = pd.DataFrame({"x": x[o], "e": e[o]})
            w = max(20, len(d) // 15)
            
            ax_.scatter(d["x"], d["e"], s=3, alpha=0.12, c=c)
            ax_.plot(d["x"], d["e"].rolling(w, min_periods=1, center=True).mean(), c=c, lw=2, label=v)
            
        ax_.set_title(f"{r['symbol']} (ρD={r['rho_DENV']:.2f}, ρZ={r['rho_ZIKV']:.2f})", fontsize=9)
        ax_.set_xlabel("log viral load")
        ax_.legend(fontsize=6)
        
    fig.suptitle("Convergent progression genes vs viral load", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure4_progression_trends", dpi, logger)

def generate_figure_5(s6_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig5: Cross-modal evidence heatmap."""
    logger.info("Generating Figure 5...")
    path = s6_dir / "common_progression_genes_Tier1.csv"
    if not path.exists(): return
    
    t1 = pd.read_csv(path).set_index("symbol")
    dm = {"up": 1, "down": -1, "ns": 0, "absent": np.nan}
    
    M = pd.DataFrame({
        "scRNA-DENV": np.sign(t1["rho_DENV_sc"]),
        "scRNA-ZIKV": np.sign(t1["rho_ZIKV_sc"]),
        "ZIKV-macrophage": t1["ZIKV_macrophage"].map(dm),
        "ZIKV-neural": t1["ZIKV_neural"].map(dm),
        "DENV-blood": t1["DENV_blood"].map(dm)
    }).sort_values(["scRNA-DENV", "ZIKV-neural"])
    
    fig, a = plt.subplots(figsize=(6, 8))
    
    # Publication-ready discrete heatmap using seaborn
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    
    # Colormap: Blue (-1), Light Grey (0), Red (1)
    cmap = ListedColormap(["#2980b9", "#ecf0f1", "#c0392b"])
    
    # Missing values (absent) will be white natively if we use mask, but they are NaN. 
    # Seaborn handles NaN by leaving it transparent/white.
    sns.heatmap(M, cmap=cmap, vmin=-1, vmax=1, ax=a, 
                linewidths=0.5, linecolor='white',
                cbar_kws={"shrink": 0.5, "ticks": [-0.66, 0, 0.66]})
                
    a.set_xticklabels(a.get_xticklabels(), rotation=45, ha="right", fontsize=11)
    a.set_yticklabels(a.get_yticklabels(), fontsize=10)
    a.set_ylabel("")
    a.set_title("Cross-modal Evidence (Tier 1)", fontsize=14, pad=15, fontweight="bold")
    
    # Fix colorbar labels
    cbar = a.collections[0].colorbar
    cbar.set_ticklabels(["Down", "ns", "Up"])
    cbar.ax.tick_params(labelsize=11) 
    
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure5_crossmodal_evidence_heatmap", dpi, logger)

def generate_figure_6(s7_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig6: Pathway enrichment barplot."""
    logger.info("Generating Figure 6...")
    path = s7_dir / "enrichment_up.csv"
    if not path.exists(): return
    
    e = pd.read_csv(path)
    e["nlp"] = -np.log10(e["Adjusted P-value"].clip(lower=1e-300))
    e["short"] = e["Term"].str.replace(r"\(GO:\d+\)", "", regex=True) \
                          .str.replace(r"R-HSA-\d+", "", regex=True).str.strip()
                          
    pal = {
        "MSigDB_Hallmark_2020": "#c0392b",
        "GO_Biological_Process_2021": "#2471a3",
        "KEGG_2021_Human": "#27ae60",
        "Reactome_2022": "#8e44ad"
    }
    
    top = e.drop_duplicates("short").head(15).sort_values("nlp")
    fig, a = plt.subplots(figsize=(9, 6))
    a.barh(range(len(top)), top["nlp"], color=[pal.get(g, "#888") for g in top["Gene_set"]])
    
    a.set_yticks(range(len(top)))
    a.set_yticklabels(top["short"], fontsize=8)
    a.set_xlabel("-log10 adj P")
    a.set_title("Enriched pathways (up)")
    a.legend(handles=[Patch(color=c, label=l.split('_')[0]) for l, c in pal.items()],
             fontsize=7, loc="lower right")
             
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure6_pathway_enrichment", dpi, logger)

def generate_figure_7(s5b_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig7: Bulk DE volcano plots."""
    logger.info("Generating Figure 7...")
    bulk_specs = [
        ("bulk_GSE118305_ZIKV_macrophage_DE.csv", "lfc", "padj", "ZIKV macrophage (GSE118305)", "-log10 adj P", 0.58),
        ("bulk_GSE78711_ZIKV_neural_DE.csv", "lfc", "p_value", "ZIKV neural (GSE78711)", "-log10 p (Cufflinks)", 1.0),
        ("bulk_GSE279208_DENV_blood_DE.csv", "lfc", "padj", "DENV blood (GSE279208)", "-log10 adj P", 0.58),
    ]
    
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5))
    for a, (fn, lc, pc, title, ylab, cut) in zip(ax, bulk_specs):
        path = s5b_dir / fn
        if not path.exists(): continue
        
        d = pd.read_csv(path)
        if pc not in d.columns and "padj" in d.columns:
            # Fallback if 'p_value' was not generated properly by upstream scripts
            pc = "padj"
            
        d[pc] = pd.to_numeric(d.get(pc, 1.0), errors="coerce")
        d[lc] = pd.to_numeric(d.get(lc, 0.0), errors="coerce")
        d = d.dropna(subset=[lc, pc])
        d["nlp"] = -np.log10(d[pc].clip(lower=1e-300))
        
        up = d["dir"] == "up"
        dn = d["dir"] == "down"
        
        a.scatter(d[lc], d["nlp"], s=5, c="#d5d8dc", linewidths=0)
        a.scatter(d[lc][up], d["nlp"][up], s=8, c="#c0392b", linewidths=0, label=f"up ({int(up.sum())})")
        a.scatter(d[lc][dn], d["nlp"][dn], s=8, c="#2471a3", linewidths=0, label=f"down ({int(dn.sum())})")
        
        a.axvline(cut, ls="--", c="k", lw=0.5)
        a.axvline(-cut, ls="--", c="k", lw=0.5)
        a.set(xlabel="log2 fold change", ylabel=ylab, title=title)
        a.legend(fontsize=8, frameon=False)
        
    fig.suptitle("Bulk RNA-seq differential expression", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, pub_dir, "Figure7_bulk_volcanoes", dpi, logger)

def generate_figure_8(s5_dir: Path, s5b_dir: Path, s6_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig8: Venn diagram single-cell vs bulk."""
    logger.info("Generating Figure 8...")
    try:
        sc_set = set(pd.read_csv(s5_dir / "convergent_progression_signature.csv")["symbol"].str.upper())
        mac = pd.read_csv(s5b_dir / "bulk_GSE118305_ZIKV_macrophage_DE.csv")
        neu = pd.read_csv(s5b_dir / "bulk_GSE78711_ZIKV_neural_DE.csv")
        bl = pd.read_csv(s5b_dir / "bulk_GSE279208_DENV_blood_DE.csv")
        
        zikv_bulk = set(mac[mac["dir"] != "ns"]["SYM"].str.upper()) | set(neu[neu["dir"] != "ns"]["SYM"].str.upper())
        denv_bulk = set(bl[bl["dir"] != "ns"]["SYM"].str.upper())
        
        n_tier1 = int((pd.read_csv(s6_dir / "crossmodal_validated_progression_signature.csv")["tier"] == 1).sum())
        
        fig, a = plt.subplots(figsize=(9, 8))
        v = venn3([sc_set, zikv_bulk, denv_bulk],
                  set_labels=(f"Single-cell\nconvergent progression\n({len(sc_set)})", "ZIKV bulk DEGs\n(macrophage + neural)", "DENV bulk DEGs\n(blood)"), ax=a)
                  
        for t in (v.set_labels or []):
            if t:
                t.set_fontsize(10)
                t.set_fontweight("bold")
        for t in (v.subset_labels or []):
            if t:
                t.set_fontsize(11)
                
        center_count = len(sc_set & zikv_bulk & denv_bulk)
        a.set_title(f"Single-cell vs bulk overlap\n(centre = {center_count} genes in single-cell + both viruses' bulk; "
                    f"{n_tier1} are direction-concordant = the common signature)", fontsize=11, fontweight="bold")
                    
        fig.tight_layout()
        save_fig(fig, pub_dir, "Figure8_venn_singlecell_vs_bulk", dpi, logger)
    except Exception as e:
        logger.error(f"Failed to generate Figure 8: {e}")

def generate_figure_9(s3_dir: Path, s5b_dir: Path, pub_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate Fig9: Cross-modal 4-set Venn diagram."""
    logger.info("Generating Figure 9...")
    try:
        def sc_dir(virus: str, sign: str) -> set:
            r = pd.read_csv(s3_dir / f"DE_High_vs_Mock_{virus}.csv")
            s = r[(r["padj"] < 0.05) & (r["log2FoldChange"].abs() > 0.58)].dropna(subset=["symbol"])
            s = s[s["log2FoldChange"] > 0] if sign == "up" else s[s["log2FoldChange"] < 0]
            return set(s["symbol"].str.upper())
            
        def bdir(df: pd.DataFrame, sign: str) -> set:
            return set(df[df["dir"] == sign]["SYM"].str.upper())
            
        _mac = pd.read_csv(s5b_dir / "bulk_GSE118305_ZIKV_macrophage_DE.csv")
        _neu = pd.read_csv(s5b_dir / "bulk_GSE78711_ZIKV_neural_DE.csv")
        _bl = pd.read_csv(s5b_dir / "bulk_GSE279208_DENV_blood_DE.csv")
        
        zikv_up = sc_dir("ZIKV", "up") | bdir(_mac, "up") | bdir(_neu, "up")
        zikv_down = sc_dir("ZIKV", "down") | bdir(_mac, "down") | bdir(_neu, "down")
        denv_up = sc_dir("DENV", "up") | bdir(_bl, "up")
        denv_down = sc_dir("DENV", "down") | bdir(_bl, "down")
        
        data4 = {
            "Upregulated DEGs in Zika": zikv_up,
            "Downregulated DEGs in Zika": zikv_down,
            "Upregulated DEGs in Dengue": denv_up,
            "Downregulated DEGs in Dengue": denv_down
        }
        
        fig, a = plt.subplots(figsize=(12, 8))
        venn4(data4, ax=a, fontsize=9, legend_loc="upper right")
        
        a.set_title(f"Cross-Modal Gene Signature Overlap  (single-cell + bulk)\n"
                    f"convergent UP = {len(zikv_up & denv_up)}  |  "
                    f"convergent DOWN = {len(zikv_down & denv_down)}", fontsize=12, fontweight="bold")
                    
        fig.tight_layout()
        save_fig(fig, pub_dir, "Figure9_CrossModal_Signature_Venn", dpi, logger)
    except Exception as e:
        logger.error(f"Failed to generate Figure 9: {e}")

def generate_supp_figure_1(adata: sc.AnnData, sup_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate FigS1: QC Violin plots."""
    logger.info("Generating Figure S1...")
    fig, ax = plt.subplots(1, 3, figsize=(11, 4))
    
    metrics = [
        ("n_genes", "Genes/cell", [2000, 8000]),
        ("total_counts", "Total counts", None),
        ("pct_mito", "% mito", [15])
    ]
    
    for a, col, lab, thr in zip(ax, *zip(*metrics)):
        if col in adata.obs:
            vp = a.violinplot(adata.obs[col].dropna().to_numpy(), showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor("#5dade2")
                b.set_alpha(0.7)
            a.set_ylabel(lab)
            a.set_xticks([])
            for t in (thr or []):
                a.axhline(t, color="r", ls="--", lw=0.7)
                
    fig.suptitle("Supplementary: single-cell QC", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, sup_dir, "FigureS1_QC_violin", dpi, logger)

def generate_supp_figure_2(s5_dir: Path, sup_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate FigS2: Progression concordance scatter."""
    logger.info("Generating Figure S2...")
    path = s5_dir / "convergent_progression_signature.csv"
    if not path.exists(): return
    
    p = pd.read_csv(path)
    c_arr = np.where(p["direction"].str.startswith("up"), "#c0392b", "#2471a3")
    
    fig, a = plt.subplots(figsize=(6, 6))
    a.scatter(p["rho_DENV"], p["rho_ZIKV"], s=10, c=c_arr, alpha=0.6, linewidths=0)
    a.axhline(0, c="k", lw=0.5)
    a.axvline(0, c="k", lw=0.5)
    a.set(xlabel="ρ DENV", ylabel="ρ ZIKV", title=f"Convergent progression genes (n={len(p)})")
    
    fig.tight_layout()
    save_fig(fig, sup_dir, "FigureS2_progression_concordance", dpi, logger)

def generate_supp_figure_3(s3_dir: Path, s5b_dir: Path, sup_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate FigS3: DEG counts per dataset barplot."""
    logger.info("Generating Figure S3...")
    counts = {}
    
    for v in ["DENV", "ZIKV"]:
        path = s3_dir / f"DE_High_vs_Mock_{v}.csv"
        if path.exists():
            r = pd.read_csv(path).dropna(subset=["padj"])
            counts[f"scRNA {v}"] = int(((r["padj"] < 0.05) & (r["log2FoldChange"].abs() > 0.58)).sum())
            
    try:
        counts["bulk ZIKV\n(macrophage)"] = int((pd.read_csv(s5b_dir / "bulk_GSE118305_ZIKV_macrophage_DE.csv")["dir"] != "ns").sum())
        counts["bulk ZIKV\n(neural)"] = int((pd.read_csv(s5b_dir / "bulk_GSE78711_ZIKV_neural_DE.csv")["dir"] != "ns").sum())
        counts["bulk DENV\n(blood)"] = int((pd.read_csv(s5b_dir / "bulk_GSE279208_DENV_blood_DE.csv")["dir"] != "ns").sum())
    except Exception as e:
        logger.warning(f"Skipping some bulk DEG counts due to error: {e}")
        
    if not counts: return
    
    fig, a = plt.subplots(figsize=(8, 4.5))
    colors = ["#c0392b", "#2471a3", "#2980b9", "#5499c7", "#e74c3c"][:len(counts)]
    a.bar(range(len(counts)), list(counts.values()), color=colors)
    a.set_xticks(range(len(counts)))
    a.set_xticklabels(list(counts.keys()), fontsize=8)
    
    for i, vv in enumerate(counts.values()):
        a.text(i, vv, str(vv), ha="center", va="bottom", fontsize=8)
        
    a.set_ylabel("# DE genes")
    a.set_title("Supplementary: DE gene counts per dataset")
    fig.tight_layout()
    save_fig(fig, sup_dir, "FigureS3_DEG_counts_per_dataset", dpi, logger)

def generate_supp_figure_4(s6_dir: Path, sup_dir: Path, dpi: int, logger: logging.Logger) -> None:
    """Generate FigS4: Cross-modal tiers barplot."""
    logger.info("Generating Figure S4...")
    path = s6_dir / "crossmodal_validated_progression_signature.csv"
    if not path.exists(): return
    
    cm = pd.read_csv(path)
    tc = cm["tier"].value_counts().sort_index()
    
    lab = {1: "Tier1\nboth viruses", 2: "Tier2\none virus", 3: "Tier3\nconflict", 4: "Tier4\nSC only"}
    fig, a = plt.subplots(figsize=(7, 4.5))
    
    colors = ["#1a7f37", "#7bb661", "#d29d00", "#b0b0b0"][:len(tc)]
    labels = [lab.get(t, str(t)) for t in tc.index]
    
    a.bar(labels, tc.values, color=colors)
    for i, vv in enumerate(tc.values):
        a.text(i, vv, str(vv), ha="center", va="bottom", fontsize=8)
        
    a.set_ylabel("# genes")
    a.set_title("Supplementary: cross-modal tiers")
    fig.tight_layout()
    save_fig(fig, sup_dir, "FigureS4_crossmodal_tiers", dpi, logger)

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 08: Figures Generation")
    apply_plot_settings(args.dpi)
    
    root = args.root_dir
    s2 = root / "Step02_clustering"
    s3 = root / "Step03_differential_expression"
    s4 = root / "Step04_core_signature"
    s5 = root / "Step05_progression"
    s5b = root / "Step05b_bulk_differential_expression"
    s6 = root / "Step06_crossmodal_validation"
    s7 = root / "Step07_enrichment"
    
    pub_dir = args.out_dir / "publication"
    sup_dir = args.out_dir / "supplementary"
    pub_dir.mkdir(parents=True, exist_ok=True)
    sup_dir.mkdir(parents=True, exist_ok=True)
    
    # Pre-load AnnData if available
    adata_path = s2 / "adata_processed.h5ad"
    if adata_path.exists():
        logger.info(f"Loading AnnData from {adata_path}...")
        adata = sc.read_h5ad(adata_path)
        if "log_viral_load" not in adata.obs:
            adata.obs["log_viral_load"] = np.log1p(adata.obs["viral_load"].astype(float))
            
        generate_figure_1(adata, pub_dir, args.dpi, logger)
        generate_figure_4(adata, s5, pub_dir, args.dpi, logger)
        generate_supp_figure_1(adata, sup_dir, args.dpi, logger)
    else:
        logger.warning(f"AnnData not found at {adata_path}. Skipping figures 1, 4, S1.")
        
    generate_figure_2(s3, pub_dir, args.dpi, logger)
    generate_figure_3(s4, pub_dir, args.dpi, logger)
    generate_figure_5(s6, pub_dir, args.dpi, logger)
    generate_figure_6(s7, pub_dir, args.dpi, logger)
    generate_figure_7(s5b, pub_dir, args.dpi, logger)
    generate_figure_8(s5, s5b, s6, pub_dir, args.dpi, logger)
    generate_figure_9(s3, s5b, pub_dir, args.dpi, logger)
    
    generate_supp_figure_2(s5, sup_dir, args.dpi, logger)
    generate_supp_figure_3(s3, s5b, sup_dir, args.dpi, logger)
    generate_supp_figure_4(s6, sup_dir, args.dpi, logger)

    logger.info("Step 08 completed successfully.")

if __name__ == "__main__":
    main()
