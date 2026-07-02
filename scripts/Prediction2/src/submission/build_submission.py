#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd


# Resolve repo root
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


# Input: results folders
CV0_DIR = f"{REPO_ROOT}/results/cv0_predictions"
CV00_DIR = f"{REPO_ROOT}/results/cv00_predictions"

# Training metadata
TRAINING_META = f"{REPO_ROOT}/data/processed/unified_training_pheno_cleaned.csv"

# Accession lists
ACC_ROOT = f"{REPO_ROOT}/data/raw/accession_lists"

# Genotyped accessions live here
GENO_ROOT = f"{REPO_ROOT}/data/predictathon"


# Output: final submission folder
OUTDIR = f"{REPO_ROOT}/submission"
os.makedirs(OUTDIR, exist_ok=True)

print("======================================")
print(" Building Final Predictathon Submission")
print("======================================\n")


# Load training metadata
train = pd.read_csv(TRAINING_META)
train["norm"] = train["germplasmName_mapped"].astype(str).str.strip().str.upper()

# All historical trials
ALL_HISTORICAL_TRIALS = sorted(train["trial"].dropna().unique())


# Helper: write correct Trials.csv files
def write_trials_files(trial, cv_dir, cv_type):
    """
    Writes CV0_Trials.csv or CV00_Trials.csv depending on cv_type.
    """

    # Load accession list for this Predictathon trial
    acc_path = f"{ACC_ROOT}/{trial}.txt"
    with open(acc_path) as f:
        focal_acc = [x.strip().upper() for x in f]

    if cv_type == "CV0":
        # CV0 uses ALL historical trials
        trials_used = ALL_HISTORICAL_TRIALS

    elif cv_type == "CV00":
        # CV00 excludes any historical trial containing a focal-trial line
        trials_with_focal_lines = (
            train[train["norm"].isin(focal_acc)]["trial"].unique().tolist()
        )
        trials_used = sorted(
            [t for t in ALL_HISTORICAL_TRIALS if t not in trials_with_focal_lines]
        )

    else:
        raise ValueError("cv_type must be CV0 or CV00")

    # Write Trials.csv
    pd.DataFrame({"studyName": trials_used}).to_csv(
        os.path.join(cv_dir, f"{cv_type}_Trials.csv"), index=False
    )


# Helper: write Accessions.csv files
def write_accessions_file(cv_dir, geno_norm, cv_type):
    pd.DataFrame({"germplasmName": geno_norm}).to_csv(
        os.path.join(cv_dir, f"{cv_type}_Accessions.csv"), index=False
    )


# Process each trial
trials = sorted([f.replace(".csv", "") for f in os.listdir(CV0_DIR) if f.endswith(".csv")])

for trial in trials:
    print(f"Processing {trial}...")

    # Load genotyped accessions
    geno_lines_path = f"{GENO_ROOT}/{trial}/processed/geno_lines.npy"
    geno_lines = np.load(geno_lines_path, allow_pickle=True).tolist()
    geno_norm = [str(x).strip() for x in geno_lines]

    # Create trial folder
    trial_out = os.path.join(OUTDIR, trial)
    os.makedirs(trial_out, exist_ok=True)


    # CV0
    cv0_pred_file = f"{CV0_DIR}/{trial}.csv"
    df0 = pd.read_csv(cv0_pred_file)
    df0 = df0.rename(columns={"line_name": "germplasmName"})
    df0 = df0[df0["germplasmName"].isin(geno_norm)]

    cv0_dir = os.path.join(trial_out, "CV0")
    os.makedirs(cv0_dir, exist_ok=True)
    df0.to_csv(os.path.join(cv0_dir, "CV0_Predictions.csv"), index=False)

    write_trials_files(trial, cv0_dir, "CV0")
    write_accessions_file(cv0_dir, geno_norm, "CV0")


    # CV00
    cv00_pred_file = f"{CV00_DIR}/{trial}.csv"
    df00 = pd.read_csv(cv00_pred_file)
    df00 = df00.rename(columns={"line_name": "germplasmName"})
    df00 = df00[df00["germplasmName"].isin(geno_norm)]

    cv00_dir = os.path.join(trial_out, "CV00")
    os.makedirs(cv00_dir, exist_ok=True)
    df00.to_csv(os.path.join(cv00_dir, "CV00_Predictions.csv"), index=False)

    write_trials_files(trial, cv00_dir, "CV00")
    write_accessions_file(cv00_dir, geno_norm, "CV00")

print(f"\n[submission] Final submission folder built → {OUTDIR}")