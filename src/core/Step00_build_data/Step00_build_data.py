#!/usr/bin/env python3
"""
STEP 00 - Build the raw AnnData from GEO single-cell data.

This script reads GSE110496.zip (2,260 count files + SOFT metadata) from the workspace's raw_data/ directory,
parses the metadata, assembles the counts matrix, labels infection states, and maps gene symbols.
The output is saved as an AnnData object in h5ad format.
"""

import argparse
import io
import logging
import re
import sys
import zipfile
from pathlib import Path
from typing import List, Dict, Any, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import warnings

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Build raw AnnData from GEO single-cell data.")
    
    # Define default paths relative to this script
    import sys
    _root_dir = Path(__file__).resolve().parent.parent.parent
    if str(_root_dir) not in sys.path: sys.path.append(str(_root_dir))
    from src import config
    default_here = config.PROCESSED_DATA_DIR / Path(__file__).resolve().parent.name
    default_root = config.PROCESSED_DATA_DIR
    default_raw = config.RAW_DATA_DIR
    default_zip = default_raw / "single_cell_GSE110496" / "GSE110496.zip"
    default_symbol_cache = default_raw / "reference" / "ensembl_to_symbol.csv"
    default_out = default_here / "adata_raw.h5ad"

    parser.add_argument("--zip-file", type=Path, default=default_zip,
                        help="Path to the GSE110496.zip file.")
    parser.add_argument("--symbol-cache", type=Path, default=default_symbol_cache,
                        help="Path to the Ensembl to Symbol mapping CSV.")
    parser.add_argument("--out-file", type=Path, default=default_out,
                        help="Path to save the output AnnData object (.h5ad).")
    parser.add_argument("--log-level", type=str, default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level.")
    return parser.parse_args()

def setup_logger(level: str) -> logging.Logger:
    """Set up the logger for the script."""
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, level.upper()))
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(ch)
    return logger

def parse_soft_metadata(zf: zipfile.ZipFile, soft_name: str, logger: logging.Logger) -> pd.DataFrame:
    """Parse metadata from the SOFT file."""
    logger.info("Parsing SOFT metadata...")
    records: Dict[str, Dict[str, str]] = {}
    cur = None
    
    try:
        with zf.open(soft_name) as raw_file:
            for line in io.TextIOWrapper(raw_file, encoding="utf-8", errors="replace"):
                line = line.rstrip("\n")
                if line.startswith("^SAMPLE"):
                    cur = line.split("=", 1)[1].strip()
                    records[cur] = {"gsm": cur}
                elif cur is None:
                    continue
                elif line.startswith("!Sample_characteristics_ch1"):
                    val = line.split("=", 1)[1].strip()
                    if ":" in val:
                        k, v = val.split(":", 1)
                        records[cur][k.strip()] = v.strip()
    except Exception as e:
        logger.error(f"Failed to read SOFT file {soft_name}: {e}")
        sys.exit(1)

    meta = pd.DataFrame.from_dict(records, orient="index")
    
    # Convert numeric columns
    numeric_cols = ["moi", "time[h]", "n_dengue_molecules", "n_zika_molecules"]
    for col in numeric_cols:
        if col in meta.columns:
            meta[col] = pd.to_numeric(meta[col], errors="coerce")
            
    # Calculate viral load
    meta["viral_load"] = np.where(
        meta["virus"].eq("dengue"),
        meta.get("n_dengue_molecules"),
        meta.get("n_zika_molecules")
    ).astype(float)
    meta["viral_load"] = meta["viral_load"].fillna(0.0)
    
    # Rename columns for consistency
    meta = meta.rename(columns={"time[h]": "time_h", "cell line": "cell_line"})
    return meta

def build_counts_matrix(zf: zipfile.ZipFile, tsv_files: List[str], logger: logging.Logger) -> Tuple[np.ndarray, List[str], List[str], List[str], List[str], np.ndarray]:
    """Assemble the counts matrix from TSV files."""
    logger.info("Building counts matrix...")
    try:
        with zf.open(tsv_files[0]) as f:
            first = pd.read_csv(f, sep="\t")
    except Exception as e:
        logger.error(f"Failed to read the first TSV file {tsv_files[0]}: {e}")
        sys.exit(1)

    features = first["feature"].tolist()
    gene_mask = ~first["feature"].str.startswith("__").to_numpy()
    htseq_rows = [f for f in features if f.startswith("__")]
    genes = [f for f in features if not f.startswith("__")]
    
    X = np.zeros((len(tsv_files), len(features)), dtype=np.float32)
    gsms = []
    
    for i, file_name in enumerate(tsv_files):
        try:
            with zf.open(file_name) as f:
                df = pd.read_csv(f, sep="\t").set_index("feature").reindex(features).fillna(0)
            X[i] = df["count"].to_numpy()
            gsms.append(Path(file_name).name.split("_")[0])
        except Exception as e:
            logger.warning(f"Error processing file {file_name}: {e}")
            
        if (i + 1) % 500 == 0:
            logger.info(f"Loaded {i+1}/{len(tsv_files)} files...")
            
    return X, gsms, htseq_rows, genes, features, gene_mask

def label_infection_states(adata: sc.AnnData) -> None:
    """Label cells with their infection states."""
    vl = adata.obs["viral_load"].to_numpy()
    moi = adata.obs["moi"].to_numpy()
    
    state = np.full(adata.n_obs, "Bystander", dtype=object)
    state[vl > 10] = "High"
    state[(vl >= 1) & (vl <= 10)] = "Low"
    state[moi == 0] = "Mock"
    
    adata.obs["infection_state"] = pd.Categorical(state, categories=["Mock", "Bystander", "Low", "High"])
    adata.obs["condition"] = np.where(
        adata.obs["infection_state"] == "Mock", 
        "Mock",
        adata.obs["virus"].astype(str).str.upper().map({"DENGUE": "DENV", "ZIKA": "ZIKV"})
    )
    adata.obs["condition"] = pd.Categorical(adata.obs["condition"])
    adata.obs["time_h"] = adata.obs["time_h"].astype("category")

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 00: Data Building")
    
    if not args.zip_file.exists():
        logger.error(f"Zip file not found: {args.zip_file}")
        sys.exit(1)
        
    try:
        zf = zipfile.ZipFile(args.zip_file)
    except zipfile.BadZipFile as e:
        logger.error(f"Invalid zip file {args.zip_file}: {e}")
        sys.exit(1)
        
    tsv_files = sorted([n for n in zf.namelist() if n.endswith("_counts.tsv") and "GSM" in n])
    soft_name_candidates = [n for n in zf.namelist() if n.endswith("family.soft")]
    
    if not soft_name_candidates:
        logger.error("No SOFT file found in zip.")
        sys.exit(1)
    soft_name = soft_name_candidates[0]
    
    logger.info(f"Found {len(tsv_files)} count files and SOFT file in zip.")
    
    # Parse Metadata
    meta = parse_soft_metadata(zf, soft_name, logger)
    
    # Parse counts
    X, gsms, htseq_rows, genes, features, gene_mask = build_counts_matrix(zf, tsv_files, logger)
    
    htseq_df = pd.DataFrame({r: X[:, features.index(r)] for r in htseq_rows}, index=gsms)
    Xg = X[:, gene_mask]
    
    obs = meta.reindex(gsms).copy()
    obs.index = gsms
    obs = obs.join(htseq_df)
    
    logger.info("Creating AnnData object...")
    adata = sc.AnnData(X=sparse.csr_matrix(Xg), obs=obs)
    adata.var_names = genes
    
    adata.var["is_ercc"] = adata.var_names.str.startswith("ERCC-")
    adata.var["is_spikein"] = adata.var["is_ercc"] | adata.var_names.str.contains("NIST", case=False)
    
    # Infection states
    logger.info("Labeling infection states...")
    label_infection_states(adata)
    
    # Map gene symbols
    if args.symbol_cache.exists():
        logger.info(f"Mapping gene symbols using {args.symbol_cache}...")
        gmap = pd.read_csv(args.symbol_cache, index_col=0)["symbol"].to_dict()
        mapped = adata.var_names.map(gmap)
        adata.var["symbol"] = pd.Index(mapped).where(pd.notna(mapped), adata.var_names).astype(str)
        adata.var["mito"] = adata.var["symbol"].str.startswith("MT-")
    else:
        logger.warning(f"Symbol cache not found at {args.symbol_cache}. Skipping symbol mapping.")
        adata.var["symbol"] = adata.var_names.astype(str)
        adata.var["mito"] = False

    logger.info(f"Saving AnnData object to {args.out_file}...")
    try:
        args.out_file.parent.mkdir(parents=True, exist_ok=True)
        adata.write(args.out_file)
        logger.info(f"Successfully saved AnnData. Shape: {adata.shape}")
        
        condition_counts = dict(adata.obs['condition'].value_counts())
        state_counts = dict(adata.obs['infection_state'].value_counts())
        logger.info(f"Condition counts: {condition_counts}")
        logger.info(f"Infection states: {state_counts}")
    except Exception as e:
        logger.error(f"Failed to save output file: {e}")
        sys.exit(1)

    logger.info("Step 00 completed successfully.")

if __name__ == "__main__":
    main()
