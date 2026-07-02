#!/usr/bin/env python3

import sys
from pathlib import Path


# Add repo root to PYTHONPATH BEFORE importing src.*
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from joblib import dump
from src.model.models import gblup_fit


def normalize(x):
    """Robust normalization for matching."""
    return str(x).strip().upper()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("trial", type=str)
    args = parser.parse_args()

    trial = args.trial
    repo = Path(__file__).resolve().parents[2]


    # Load unified phenotype file
    pheno = pd.read_csv(
        repo / "data" / "processed" / "unified_training_pheno_mapped.csv"
    )
    pheno["germplasmName_mapped_norm"] = pheno["germplasmName_mapped"].apply(normalize)

    # Filter to this trial
    pheno_trial = pheno[pheno["trial"] == trial].copy()

    # Original + normalized trial lines
    trial_lines_raw = (
        pheno_trial["germplasmName_mapped"]
        .dropna()
        .unique()
        .tolist()
    )
    trial_lines = [normalize(x) for x in trial_lines_raw]


    # Load global GRM + samples
    grm_dir = repo / "data" / "processed" / "global_union"
    GRM = np.load(grm_dir / "GRM_global_union.npy")

    with open(grm_dir / "G_global_union_samples.txt") as f:
        global_samples_raw = [x.strip() for x in f]

    global_samples = [normalize(x) for x in global_samples_raw]


    # Filter out lines not present in the global GRM
    missing = [l for l in trial_lines if l not in global_samples]
    trial_lines = [l for l in trial_lines if l in global_samples]

    if missing:
        print(f"[train_global] Warning: {len(missing)} lines missing from global GRM and skipped.")


    # Slice GRM to trial lines
    idx = [global_samples.index(l) for l in trial_lines]
    GRM_sub = GRM[np.ix_(idx, idx)]


    # Extract phenotype vector (line means)
    pheno_map = (
        pheno_trial.groupby("germplasmName_mapped_norm")["value"]
        .mean()
        .to_dict()
    )

    y = np.array([pheno_map[l] for l in trial_lines], dtype=float)


    # Fit GBLUP
    model = gblup_fit(GRM_sub, y)


    # Save model + GRM slice
    outdir = repo / "trained_models" / trial
    outdir.mkdir(parents=True, exist_ok=True)

    dump(model, outdir / "final_model.joblib")
    np.save(outdir / "GRM.npy", GRM_sub)

    # Save original (unnormalized) names for readability
    with open(outdir / "GRM_lines.txt", "w") as f:
        for raw, norm in zip(trial_lines_raw, trial_lines):
            if normalize(raw) == norm:
                f.write(f"{raw}\n")

    print(f"✓ Trained global model for {trial}")


if __name__ == "__main__":
    main()