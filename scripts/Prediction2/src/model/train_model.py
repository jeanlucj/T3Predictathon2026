#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import joblib
import yaml

from src.model.models import fit_model



# Helpers

def load_config(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("trial", type=str, help="Predictathon trial name")
    parser.add_argument("--config", type=str, default="config.yaml")
    return parser.parse_args()


def subset_pheno_for_trial(unified_pheno: pd.DataFrame, trial: str) -> pd.DataFrame:
    """
    Subset unified_training_pheno_mapped.csv to rows belonging to this trial.
    """
    trial_cols = [
        c for c in unified_pheno.columns
        if "trial" in c.lower() or "study" in c.lower()
    ]
    if not trial_cols:
        raise ValueError("Could not infer trial column in unified phenotype file.")

    trial_col = trial_cols[0]
    return unified_pheno[unified_pheno[trial_col] == trial].copy()



# Main

def main():
    args = parse_args()
    trial = args.trial

    # Load config
    config = load_config(Path(args.config))
    paths = config["paths"]

    predictathon_root = Path(paths["predictathon_root"])
    trained_models_root = Path(paths["trained_models_root"])
    unified_pheno_path = Path(config["phenotypes"]["unified_training"])

    # Trial directories
    trial_raw_dir = predictathon_root / trial
    trial_proc_dir = trial_raw_dir / "processed"
    outdir = trained_models_root / trial
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\nTraining model for {trial}\n")


    # Load unified phenotype file
    if not unified_pheno_path.exists():
        raise SystemExit(f"Unified phenotype file not found: {unified_pheno_path}")

    unified_pheno = pd.read_csv(unified_pheno_path, low_memory=False)

    # Subset to this trial
    pheno = subset_pheno_for_trial(unified_pheno, trial)

    if pheno.empty:
        print(f"[train_model] No phenotype rows found for {trial}.")
        print("[train_model] Model will be genotype-only.")
    else:
        if "germplasmName_mapped" not in pheno.columns:
            raise SystemExit("Unified phenotype file must contain 'germplasmName_mapped'.")

        if "value" not in pheno.columns:
            raise SystemExit("Unified phenotype file must contain a 'value' column.")

        pheno["value"] = pd.to_numeric(pheno["value"], errors="coerce")
        bad = pheno["value"].isna().sum()
        if bad > 0:
            print(f"[train_model] Warning: Dropping {bad} rows with non-numeric phenotype values.")
            pheno = pheno.dropna(subset=["value"])

        print(f"✓ Loaded phenotype: {pheno.shape[0]} rows, "
              f"{pheno['germplasmName_mapped'].nunique()} unique mapped lines")


    # Load genotype (trial-specific, used only for prediction)
    geno_numeric_path = trial_proc_dir / "geno_numeric.npy"
    geno_lines_path = trial_proc_dir / "geno_lines.npy"

    if not geno_numeric_path.exists():
        raise SystemExit(f"Genotype numeric matrix not found: {geno_numeric_path}")
    if not geno_lines_path.exists():
        raise SystemExit(f"Genotype line list not found: {geno_lines_path}")

    geno_numeric = np.load(geno_numeric_path)
    geno_lines = np.load(geno_lines_path, allow_pickle=True).tolist()

    print(f"✓ Loaded genotype: {len(geno_lines)} lines × {geno_numeric.shape[1]} markers")


    # Load GLOBAL GRM + sample list
    grm_dir = Path(paths["global_grm_root"])
    global_grm_path = grm_dir / "GRM_global_union.npy"
    global_samples_path = grm_dir / "G_global_union_samples.txt"

    if not global_grm_path.exists():
        raise SystemExit(f"Global GRM not found: {global_grm_path}")
    if not global_samples_path.exists():
        raise SystemExit(f"Global GRM sample list not found: {global_samples_path}")

    G = np.load(global_grm_path)
    with open(global_samples_path) as f:
        global_samples = [x.strip() for x in f]

    print(f"✓ Loaded GLOBAL GRM: {G.shape}")


    # Save GRM + GRM lines for CV0/CV00
    out_grm_path = outdir / "GRM.npy"
    np.save(out_grm_path, G)

    out_lines_path = outdir / "GRM_lines.txt"
    with open(out_lines_path, "w") as f:
        for g in global_samples:
            f.write(g + "\n")

    print(f"✓ Saved GLOBAL GRM → {out_grm_path}")
    print(f"✓ Saved GLOBAL GRM lines → {out_lines_path}")


    # Fit model (using GLOBAL samples, not trial samples)
    print("\nFitting model...")

    model = fit_model(
        train_pheno=pheno if not pheno.empty else None,
        geno_numeric=geno_numeric,
        geno_lines=global_samples,   # aligned to global GRM
        G=G,
        model_type="gblup"
    )


    # Save model
    model_path = outdir / "final_model.joblib"
    joblib.dump(model, model_path)

    print(f"✓ Saved model → {model_path}")
    print("✓ Training complete.\n")


if __name__ == "__main__":
    main()