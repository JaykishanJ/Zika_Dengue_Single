#!/usr/bin/env python3
"""
STEP 05b - Bulk RNA-seq differential expression (computed in-house from raw data).

This script calculates differential expression for three independent bulk cohorts.
Each DE is computed and saved as a table with a direction call (up/down/ns) per gene.
These results feed into Step 06 (cross-modal validation).

Datasets:
  * GSE118305: ZIKV-infected(4G2+) vs mock, human macrophages -> Welch t-test on log2(FPKM+1)
  * GSE78711: ZIKV vs mock, human neural progenitors (hNPC) -> Cufflinks log2FC
  * GSE279208: dengue vs healthy, whole blood -> Welch t-test on VST matrix

Outputs:
    - bulk_GSE118305_ZIKV_macrophage_DE.csv
    - bulk_GSE78711_ZIKV_neural_DE.csv
    - bulk_GSE279208_DENV_blood_DE.csv
"""

import argparse
import gzip
import logging
import re
import sys
import warnings
from pathlib import Path
from typing import List

import mygene
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Compute bulk RNA-seq differential expression.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_raw = default_root / "raw_data"
    
    default_g118305 = default_raw / "bulk" / "GSE118305" / "RNAseq_HMDM_FPKM.txt.gz"
    default_g78711 = default_raw / "bulk" / "GSE78711" / "gene_exp.txt.gz"
    default_vst = default_raw / "bulk" / "GSE279208" / "GSE279208_VST_Normalized_Matrix.csv"
    default_pheno = default_raw / "bulk" / "GSE279208" / "GSE279208_Final_Phenotype.csv"

    # Input Files
    parser.add_argument("--g118305", type=Path, default=default_g118305, help="Path to GSE118305 FPKM data.")
    parser.add_argument("--g78711", type=Path, default=default_g78711, help="Path to GSE78711 gene expression data.")
    parser.add_argument("--vst", type=Path, default=default_vst, help="Path to GSE279208 VST matrix.")
    parser.add_argument("--pheno", type=Path, default=default_pheno, help="Path to GSE279208 phenotype data.")

    # Output Dir
    parser.add_argument("--out-dir", type=Path, default=default_here, help="Output directory for CSVs.")

    # Thresholds
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold.")
    parser.add_argument("--lfc-rep", type=float, default=0.58, help="Log2 Fold Change threshold for t-tests.")
    parser.add_argument("--lfc-cuff", type=float, default=1.0, help="Log2 Fold Change threshold for Cufflinks data.")
    
    # Logging
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level.")
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

def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.empty(n)
    
    adj[o] = p[o] * n / (np.arange(n) + 1)
    run = 1.0
    for i in o[::-1]:
        run = min(run, adj[i])
        adj[i] = run
    return np.clip(adj, 0, 1)

def assign_direction(lfc: float, padj: float, cut: float, fdr: float, use_padj: bool = True) -> str:
    """Assign differential expression direction based on log2FC and padj."""
    if use_padj and not (padj < fdr):
        return "ns"
    if lfc > cut:
        return "up"
    elif lfc < -cut:
        return "down"
    return "ns"

def process_gse118305(in_file: Path, out_dir: Path, fdr: float, lfc_cut: float, logger: logging.Logger) -> None:
    """Process GSE118305: ZIKV macrophage data (infected vs mock)."""
    logger.info("Processing GSE118305 (ZIKV macrophage)...")
    if not in_file.exists():
        logger.error(f"File not found: {in_file}")
        return

    try:
        with gzip.open(in_file, "rt") as f:
            header = f.readline()
            
        match = re.search(r'-d\s+(.*?)\)', header)
        if match:
            samp = match.group(1).split()
        else:
            raise ValueError("Failed to parse samples from header.")
            
        df = pd.read_csv(in_file, sep="\t", header=0)
        df.columns = ["ID", "chr", "start", "end", "strand", "Length", "gene_name", "annot", "gene_id", "gene_type"] + samp
        
        inf = [s for s in samp if "4g2pos" in s.lower()]
        mock = [s for s in samp if "mock" in s.lower()]
        
        expr = df.groupby("gene_name")[inf + mock].sum()
        L = np.log2(expr + 1.0)
        
        lfc = L[inf].to_numpy().mean(1) - L[mock].to_numpy().mean(1)
        _, p = ttest_ind(L[inf], L[mock], axis=1, equal_var=False)
        
        mac = pd.DataFrame({"SYM": expr.index.str.upper(), "lfc": lfc, "padj": benjamini_hochberg(np.nan_to_num(p, nan=1.0))})
        mac["dir"] = [assign_direction(l, q, lfc_cut, fdr) for l, q in zip(mac["lfc"], mac["padj"])]
        mac = mac.drop_duplicates("SYM")
        
        out_csv = out_dir / "bulk_GSE118305_ZIKV_macrophage_DE.csv"
        mac.to_csv(out_csv, index=False)
        
        n_up = (mac["dir"] == "up").sum()
        n_dn = (mac["dir"] == "down").sum()
        logger.info(f"GSE118305 summary: {len(inf)} inf vs {len(mock)} mock -> {n_up} up, {n_dn} down")
        
    except Exception as e:
        logger.error(f"Failed to process GSE118305: {e}")

def process_gse78711(in_file: Path, out_dir: Path, lfc_cut: float, logger: logging.Logger) -> None:
    """Process GSE78711: ZIKV neural data (Cufflinks)."""
    logger.info("Processing GSE78711 (ZIKV neural)...")
    if not in_file.exists():
        logger.error(f"File not found: {in_file}")
        return

    try:
        z = pd.read_csv(in_file, sep="\t")
        z.columns = [c.strip() for c in z.columns]
        
        z["lfc"] = pd.to_numeric(z["log2.fold_change."], errors="coerce")
        z["SYM"] = z["gene"].astype(str).str.upper()
        
        z = z.replace([np.inf, -np.inf], np.nan).dropna(subset=["lfc"]).drop_duplicates("SYM")
        z["dir"] = [assign_direction(l, 0.0, lfc_cut, 1.0, use_padj=False) for l in z["lfc"]]
        
        out_csv = out_dir / "bulk_GSE78711_ZIKV_neural_DE.csv"
        z.to_csv(out_csv, index=False)
        
        n_up = (z["dir"] == "up").sum()
        n_dn = (z["dir"] == "down").sum()
        logger.info(f"GSE78711 summary: {n_up} up, {n_dn} down")
        
    except Exception as e:
        logger.error(f"Failed to process GSE78711: {e}")

def process_gse279208(vst_file: Path, pheno_file: Path, out_dir: Path, fdr: float, lfc_cut: float, logger: logging.Logger) -> None:
    """Process GSE279208: DENV blood data."""
    logger.info("Processing GSE279208 (DENV blood)...")
    if not vst_file.exists() or not pheno_file.exists():
        logger.error(f"Missing GSE279208 files.")
        return

    try:
        vst = pd.read_csv(vst_file, index_col=0)
        ph = pd.read_csv(pheno_file)
        
        grp = ph.set_index("SampleID")["Group"].to_dict()
        
        dis = [c for c in vst.columns if grp.get(c) in ("Dengue", "Severe_Dengue")]
        hlt = [c for c in vst.columns if grp.get(c) == "healthy"]
        
        A = vst[dis].to_numpy()
        B = vst[hlt].to_numpy()
        
        lfc = A.mean(1) - B.mean(1)
        _, p = ttest_ind(A, B, axis=1, equal_var=False)
        
        den = pd.DataFrame({"entrez": vst.index.astype(str), "lfc": lfc, "padj": benjamini_hochberg(np.nan_to_num(p, nan=1.0))})
        
        logger.info("Querying MyGene for symbol mapping (requires internet)...")
        mg = mygene.MyGeneInfo()
        res = mg.querymany(den["entrez"].tolist(), scopes="entrezgene", fields="symbol", species="human", as_dataframe=True, verbose=False)
        
        den["SYM"] = den["entrez"].map(res[~res.index.duplicated()]["symbol"]).str.upper()
        den = den.dropna(subset=["SYM"]).drop_duplicates("SYM")
        den["dir"] = [assign_direction(l, q, lfc_cut, fdr) for l, q in zip(den["lfc"], den["padj"])]
        
        out_csv = out_dir / "bulk_GSE279208_DENV_blood_DE.csv"
        den.to_csv(out_csv, index=False)
        
        n_up = (den["dir"] == "up").sum()
        n_dn = (den["dir"] == "down").sum()
        logger.info(f"GSE279208 summary: {len(dis)} dengue vs {len(hlt)} healthy -> {n_up} up, {n_dn} down")
        
    except Exception as e:
        logger.error(f"Failed to process GSE279208: {e}")

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 05b: Bulk Differential Expression")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    
    process_gse118305(args.g118305, args.out_dir, args.fdr, args.lfc_rep, logger)
    process_gse78711(args.g78711, args.out_dir, args.lfc_cuff, logger)
    process_gse279208(args.vst, args.pheno, args.out_dir, args.fdr, args.lfc_rep, logger)
    
    logger.info("Step 05b completed successfully.")

if __name__ == "__main__":
    main()
