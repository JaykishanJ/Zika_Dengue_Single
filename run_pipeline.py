#!/usr/bin/env python3
"""
Master Orchestration Script for Zika/Dengue Publication Pipeline

This script automatically runs the full analysis pipeline sequentially:
1. Core Single Cell & Bulk RNA-seq Processing (src/core/Step00 -> Step09)
2. Visualization and Figure Generation (src/viz/Step_*)
"""

import subprocess
import sys
from pathlib import Path
import logging

def setup_logger():
    logger = logging.getLogger("PipelineRunner")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def main():
    logger = setup_logger()
    ROOT = Path(__file__).resolve().parent

    # Core processing scripts that must run in order
    core_pipeline = [
        "src/core/Step00_build_data/Step00_build_data.py",
        "src/core/Step01_qc/Step01_qc.py",
        "src/core/Step02_clustering/Step02_clustering.py",
        "src/core/Step03_differential_expression/Step03_differential_expression.py",
        "src/core/Step04_core_signature/Step04_core_signature.py",
        "src/core/Step05_progression/Step05_progression.py",
        "src/core/Step05b_bulk_differential_expression/Step05b_bulk_DE.py",
        "src/core/Step06_crossmodal_validation/Step06_crossmodal_validation.py",
        "src/core/Step07_enrichment/Step07_enrichment.py",
        "src/core/Step08_figures/Step08_figures.py",
        "src/core/Step09_marker_analysis/Step09_marker_analysis.py"
    ]

    logger.info("=======================================")
    logger.info("STARTING PIPELINE EXECUTION")
    logger.info("=======================================")
    
    if len(sys.argv) > 1 and sys.argv[1] == '--viz-only':
        logger.info("\nSkipping Part 1 (Core Data Processing) due to --viz-only flag.")
    else:
        logger.info("\n>>> PART 1: CORE DATA PROCESSING <<<")
        for script in core_pipeline:
            script_path = ROOT / script
            logger.info(f"\n--- Executing {script_path.name} ---")
            
            if not script_path.exists():
                logger.error(f"Cannot find script at {script_path}")
                sys.exit(1)
                
            result = subprocess.run([sys.executable, str(script_path)], cwd=script_path.parent)
            if result.returncode != 0:
                logger.error(f"Pipeline halted: {script_path.name} failed with exit code {result.returncode}")
                sys.exit(1)

    logger.info("\n>>> PART 2: VISUALIZATION & FIGURES <<<")
    # Execute all visualization scripts
    viz_scripts = sorted(list((ROOT / "src/viz").rglob("*.py")))
    
    success_count = 0
    fail_count = 0
    
    for script_path in viz_scripts:
        logger.info(f"\n--- Generating Figure using {script_path.name} ---")
        result = subprocess.run([sys.executable, str(script_path)], cwd=script_path.parent)
        
        if result.returncode != 0:
            logger.warning(f"Error executing {script_path.name}. (Continuing to next script...)")
            fail_count += 1
        else:
            success_count += 1

    logger.info("=======================================")
    logger.info("PIPELINE EXECUTION COMPLETE")
    logger.info(f"Visualization Summary: {success_count} succeeded, {fail_count} failed.")
    logger.info("=======================================")

if __name__ == "__main__":
    main()
