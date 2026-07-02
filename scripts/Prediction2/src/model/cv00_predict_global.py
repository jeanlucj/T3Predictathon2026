#!/usr/bin/env python3

import sys
from pathlib import Path


# Add repo root to PYTHONPATH BEFORE importing src.*
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import yaml
from src.model.models import gblup_fit, gblup_predict


def load_config(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def normalize(x):
    """Robust string normalization for matching."""
    return str(x).strip().upper()


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial", type=str, required=True)
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()

    trial = args.trial
    config = load_config(Path(args.config))

    repo = Path(__file__).resolve().parents[2]


    # Load phenotype (historical)
    pheno_path = Path(config["phenotypes"]["unified_training"])
    pheno = pd.read_csv(pheno_path)
    pheno["germplasmName_mapped_norm"] = pheno["germplasmName_mapped"].apply(normalize)


    # Load global GRM + samples
    grm_dir = Path(config["paths"]["global_grm_root"])
    GRM = np.load(grm_dir / "GRM_global_union.npy")

    with open(grm_dir / "G_global_union_samples.txt") as f:
        global_samples_raw = [x.strip() for x in f]
    global_samples = [normalize(x) for x in global_samples_raw]


    # Load full accession list for this trial
    acc_path = repo / "data" / "raw" / "accession_lists" / f"{trial}.txt"
    with open(acc_path) as f:
        acc_raw = [x.strip() for x in f]
    acc_norm = [normalize(x) for x in acc_raw]

    # Lines in this trial that are also in the global GRM
    trial_norm_in_global = [a for a in acc_norm if a in global_samples]

  
    # CV00 masking: remove ALL phenotypes for these lines
    pheno_cv = pheno[
        ~pheno["germplasmName_mapped_norm"].isin(trial_norm_in_global)
    ].copy()


    # Build phenotype vector aligned to global samples
    pheno_map = (
        pheno_cv.groupby("germplasmName_mapped_norm")["value"]
        .mean()
        .to_dict()
    )

    y = []
    mask = []
    for line in global_samples:
        if line in pheno_map:
            y.append(pheno_map[line])
            mask.append(True)
        else:
            y.append(0.0)
            mask.append(False)

    y = np.array(y)
    mask = np.array(mask)


    # Fit model with CV00 masking
    GRM_train = GRM[np.ix_(mask, mask)]
    y_train = y[mask]
    model = gblup_fit(GRM_train, y_train)


    # Predict for ALL accessions in the list
    #   - genotyped: use GRM
    #   - non-genotyped: use model.mu
    genotyped = [a for a in acc_norm if a in global_samples]
    missing = [a for a in acc_norm if a not in global_samples]

    # Genotyped predictions
    idx = [global_samples.index(a) for a in genotyped]
    GRM_pred = GRM[np.ix_(idx, mask)]
    preds_geno = gblup_predict(model, GRM_pred)

    # Non-genotyped → population mean
    preds_missing = [model.mu] * len(missing)

    # Map normalized → raw names
    norm_to_raw = {normalize(r): r for r in acc_raw}

    # Build prediction map
    pred_map = {}
    for a, p in zip(genotyped, preds_geno):
        pred_map[a] = p
    for a, p in zip(missing, preds_missing):
        pred_map[a] = p

    out_lines = [norm_to_raw[a] for a in acc_norm]
    out_preds = [pred_map[a] for a in acc_norm]

    out_df = pd.DataFrame({
        "line_name": out_lines,
        "prediction": out_preds
    })
    out_df.to_csv(args.out, index=False)

    print(f"✓ CV00 predictions complete for {trial}")


if __name__ == "__main__":
    main()