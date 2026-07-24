#!/usr/bin/env python3
"""
STEP 06 - Symmetric cross-modal validation vs 3 bulk cohorts.

This script takes the single-cell convergent progression signature and validates it
against 3 independent bulk RNA-seq cohorts. It categorizes genes into tiers based
on their consistency across both single-cell and bulk datasets.

Inputs:
    - convergent_progression_signature.csv (from Step05_progression)
    - bulk_GSE118305_ZIKV_macrophage_DE.csv (from Step05b)
    - bulk_GSE78711_ZIKV_neural_DE.csv (from Step05b)
    - bulk_GSE279208_DENV_blood_DE.csv (from Step05b)

Outputs:
    - crossmodal_validated_progression_signature.csv: Full validated signature with tiers
    - common_progression_genes_Tier1.csv: Subset of Tier 1 genes
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Cross-modal validation against bulk cohorts.")
    
    # Default paths
    default_here = Path(__file__).resolve().parent
    default_root = default_here.parent
    default_prog = default_root / "Step05_progression" / "convergent_progression_signature.csv"
    default_bulk_dir = default_root / "Step05b_bulk_differential_expression"
    default_out_full = default_here / "crossmodal_validated_progression_signature.csv"
    default_out_tier1 = default_here / "common_progression_genes_Tier1.csv"

    # Input Files
    parser.add_argument("--prog-file", type=Path, default=default_prog, help="Path to convergent progression signature CSV.")
    parser.add_argument("--bulk-dir", type=Path, default=default_bulk_dir, help="Directory containing bulk DE CSVs.")
    
    # Outputs
    parser.add_argument("--out-full", type=Path, default=default_out_full, help="Path to save full validated signature CSV.")
    parser.add_argument("--out-tier1", type=Path, default=default_out_tier1, help="Path to save Tier 1 genes CSV.")

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

def get_bulk_direction(df: pd.DataFrame, gene: str) -> str:
    """Retrieve direction for a gene from a bulk dataframe."""
    if gene in df.index:
        return str(df.loc[gene, "dir"])
    return "absent"

def load_bulk_data(bulk_dir: Path, logger: logging.Logger) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the 3 bulk DE tables."""
    def _load(filename: str, desc: str) -> pd.DataFrame:
        path = bulk_dir / filename
        if not path.exists():
            logger.error(f"Missing bulk DE file: {path}")
            raise FileNotFoundError(f"Missing {path}")
        df = pd.read_csv(path).drop_duplicates("SYM").set_index("SYM")
        n_up = (df["dir"] == "up").sum()
        n_dn = (df["dir"] == "down").sum()
        logger.info(f"Loaded {desc}: {n_up} up / {n_dn} down")
        return df
        
    mac = _load("bulk_GSE118305_ZIKV_macrophage_DE.csv", "macrophage")
    neu = _load("bulk_GSE78711_ZIKV_neural_DE.csv", "neural")
    blood = _load("bulk_GSE279208_DENV_blood_DE.csv", "blood")
    
    return mac, neu, blood

def evaluate_evidence(row: pd.Series) -> pd.Series:
    """Format evidence strings for the output table."""
    d = row["sc_direction"]
    e = ["scRNA-DENV", "scRNA-ZIKV"]
    cf = []
    
    mapping = [
        ("ZIKV-macrophage", "ZIKV_macrophage"),
        ("ZIKV-neural", "ZIKV_neural"),
        ("DENV-blood", "DENV_blood")
    ]
    
    for nm, cl in mapping:
        if row[cl] == d:
            e.append(nm)
        elif row[cl] in ("up", "down"):
            cf.append(nm)
            
    return pd.Series({
        "n_datasets_confirm": len(e),
        "confirmed_in": "; ".join(e),
        "conflicts_in": "; ".join(cf) if cf else "none"
    })

def perform_validation(sc: pd.DataFrame, mac: pd.DataFrame, neu: pd.DataFrame, blood: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Validate progression genes against bulk data and assign tiers."""
    logger.info("Validating genes against bulk cohorts...")
    rows: List[Dict[str, Any]] = []
    
    for _, r in sc.iterrows():
        g = r["SYM"]
        d = r["sc_dir"]
        
        zmac = get_bulk_direction(mac, g)
        zneu = get_bulk_direction(neu, g)
        db = get_bulk_direction(blood, g)
        
        zmatch = sum(x == d for x in [zmac, zneu])
        zconf = any(x in ("up", "down") and x != d for x in [zmac, zneu])
        
        dmatch = int(db == d)
        dconf = db in ("up", "down") and db != d
        
        zok = (zmatch >= 1) and not zconf
        dok = (dmatch == 1) and not dconf
        
        nconf = zmatch + dmatch
        anyc = zconf or dconf
        
        if zok and dok:
            tier = 1
        elif (zok or dok) and not anyc:
            tier = 2
        elif anyc:
            tier = 3
        else:
            tier = 4
            
        gold = (zmac == d) and (zneu == d) and (db == d)
        
        rows.append({
            "symbol": r["symbol"],
            "sc_direction": d,
            "rho_DENV_sc": round(r.get("rho_DENV", 0), 3),
            "rho_ZIKV_sc": round(r.get("rho_ZIKV", 0), 3),
            "ZIKV_macrophage": zmac,
            "ZIKV_neural": zneu,
            "DENV_blood": db,
            "n_bulk_confirm": nconf,
            "all3_agree": gold,
            "conflict": anyc,
            "tier": tier,
            "in_DE_core": bool(r.get("in_DE_core", False))
        })
        
    out = pd.DataFrame(rows).sort_values(
        ["tier", "n_bulk_confirm", "sc_direction", "symbol"],
        ascending=[True, False, True, True]
    )
    
    # Append formatted evidence columns
    out = pd.concat([out.reset_index(drop=True), out.reset_index(drop=True).apply(evaluate_evidence, axis=1)], axis=1)
    
    return out

def main():
    """Main execution entry point."""
    args = parse_args()
    logger = setup_logger(args.log_level)
    
    logger.info("Starting Step 06: Cross-modal Validation")
    
    if not args.prog_file.exists():
        logger.error(f"Progression file not found: {args.prog_file}")
        sys.exit(1)
        
    logger.info(f"Loading progression signature from {args.prog_file}...")
    try:
        sc = pd.read_csv(args.prog_file)
    except Exception as e:
        logger.error(f"Failed to read progression file: {e}")
        sys.exit(1)
        
    sc["sc_dir"] = np.where(sc["direction"].str.startswith("up"), "up", "down")
    sc["SYM"] = sc["symbol"].str.upper()
    logger.info(f"Single-cell convergent progression genes: {len(sc)}")
    
    try:
        mac, neu, blood = load_bulk_data(args.bulk_dir, logger)
    except Exception as e:
        logger.error(f"Failed to load bulk data: {e}")
        sys.exit(1)
        
    out = perform_validation(sc, mac, neu, blood, logger)
    
    try:
        args.out_full.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_full, index=False)
        out[out["tier"] == 1].to_csv(args.out_tier1, index=False)
    except Exception as e:
        logger.error(f"Failed to save output files: {e}")
        sys.exit(1)
        
    n_gold = int(out["all3_agree"].sum())
    n_t1 = int((out["tier"] == 1).sum())
    n_t2 = int((out["tier"] == 2).sum())
    
    logger.info(f"Summary: GOLD(all3)={n_gold}, Tier1={n_t1}, Tier2={n_t2}")
    
    if n_t1 > 0:
        t1_genes = out[out["tier"] == 1]["symbol"].tolist()
        logger.info(f"Tier1 genes: {', '.join(t1_genes)}")
        
    logger.info("Step 06 completed successfully.")

if __name__ == "__main__":
    main()
