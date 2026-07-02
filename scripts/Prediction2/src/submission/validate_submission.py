#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd


# Paths
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SUBMISSION = f"{REPO_ROOT}/submission"
GENO_ROOT = f"{REPO_ROOT}/data/predictathon"
ACC_ROOT = f"{REPO_ROOT}/data/raw/accession_lists"
TRAINING_META = f"{REPO_ROOT}/data/processed/unified_training_pheno_cleaned.csv"

print("======================================")
print(" Validating Predictathon Submission")
print("======================================\n")


# Load historical phenotype metadata
train = pd.read_csv(TRAINING_META)
train["norm"] = train["germplasmName_mapped"].astype(str).str.strip().str.upper()

ALL_HIST_TRIALS = set(train["trial"].dropna().unique())

# Track errors
errors = []

def normalize(x):
    return str(x).strip().upper()


# Validate each trial folder
for trial in sorted(os.listdir(SUBMISSION)):
    trial_path = os.path.join(SUBMISSION, trial)
    if not os.path.isdir(trial_path):
        continue

    print(f"Checking {trial}...")


    # Validate that trial exists in predictathon data
    geno_dir = os.path.join(GENO_ROOT, trial)
    if not os.path.isdir(geno_dir):
        errors.append(f"{trial}: No matching folder in data/predictathon/")
        continue

    # Load genotyped accessions
    geno_lines_path = os.path.join(geno_dir, "processed", "geno_lines.npy")
    geno_lines = np.load(geno_lines_path, allow_pickle=True).tolist()
    geno_norm = {normalize(x) for x in geno_lines}

    # Load focal trial accession list
    acc_path = os.path.join(ACC_ROOT, f"{trial}.txt")
    focal_acc = {normalize(x) for x in open(acc_path)}

    # Historical trials containing focal-trial lines
    trials_with_focal = set(
        train[train["norm"].isin(focal_acc)]["trial"].unique().tolist()
    )
    expected_cv00_trials = ALL_HIST_TRIALS - trials_with_focal


    # Validate CV0 and CV00 folders
    for cv_type in ["CV0", "CV00"]:
        cv_dir = os.path.join(trial_path, cv_type)
        if not os.path.isdir(cv_dir):
            errors.append(f"{trial}: Missing {cv_type}/ folder")
            continue

        # Required files
        pred_file = os.path.join(cv_dir, f"{cv_type}_Predictions.csv")
        trials_file = os.path.join(cv_dir, f"{cv_type}_Trials.csv")
        acc_file = os.path.join(cv_dir, f"{cv_type}_Accessions.csv")

        for f in [pred_file, trials_file, acc_file]:
            if not os.path.exists(f):
                errors.append(f"{trial}: Missing file {f}")
                continue


        # Validate Predictions.csv
        df = pd.read_csv(pred_file)

        # Column check
        if list(df.columns) != ["germplasmName", "prediction"]:
            errors.append(f"{trial}/{cv_type}: Incorrect columns in Predictions.csv")

        # Genotyped-only check
        preds_set = {normalize(x) for x in df["germplasmName"]}
        if not preds_set.issubset(geno_norm):
            errors.append(f"{trial}/{cv_type}: Predictions include ungenotyped accessions")

        # Missing or non-finite predictions
        if df["prediction"].isna().any():
            errors.append(f"{trial}/{cv_type}: Missing prediction values")
        if not np.isfinite(df["prediction"]).all():
            errors.append(f"{trial}/{cv_type}: Non-finite prediction values")


        # Validate Trials.csv
        trials_df = pd.read_csv(trials_file)
        submitted_trials = set(trials_df["studyName"])

        if cv_type == "CV0":
            # CV0 must include ALL historical trials
            if submitted_trials != ALL_HIST_TRIALS:
                errors.append(f"{trial}/{cv_type}: CV0 Trials.csv does not match all historical trials")

        elif cv_type == "CV00":
            # CV00 must exclude trials containing focal lines
            if submitted_trials != expected_cv00_trials:
                errors.append(f"{trial}/{cv_type}: CV00 Trials.csv incorrect")

            # CV00 must be subset of CV0
            if not submitted_trials.issubset(ALL_HIST_TRIALS):
                errors.append(f"{trial}/{cv_type}: CV00 Trials.csv contains unexpected trials")


        # Validate Accessions.csv
        acc_df = pd.read_csv(acc_file)
        submitted_acc = {normalize(x) for x in acc_df["germplasmName"]}

        if submitted_acc != geno_norm:
            errors.append(f"{trial}/{cv_type}: Accessions.csv does not match genotyped accessions")

    print(f"  ✓ {trial} passed basic checks\n")


# Final report
print("======================================")
if errors:
    print(" Submission has issues:")
    for e in errors:
        print(" -", e)
else:
    print(" All submission files validated successfully")
print("======================================")