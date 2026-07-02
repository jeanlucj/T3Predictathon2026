#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np


# Resolve repo root
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


# Input locations
HIST_PHENO = f"{REPO_ROOT}/data/raw/pheno_processed.csv"
PREDICTATHON_DIR = f"{REPO_ROOT}/data/predictathon"


# Output locations
OUTDIR = f"{REPO_ROOT}/data/processed"
os.makedirs(OUTDIR, exist_ok=True)

OUTFILE = f"{OUTDIR}/unified_training_pheno_cleaned.csv"
TRIAL_STATS_OUT = f"{OUTDIR}/trial_yield_stats_raw.csv"


# Load historical phenotype
print("[unified_pheno] Loading historical phenotype...")

hist = pd.read_csv(
    HIST_PHENO,
    encoding="latin1",
    engine="python",
    dtype=str,
    on_bad_lines="skip"
)


# Extract grain yield column
yield_cols = [
    "Grain yield - kg/ha|CO_321:0001218",
    "Grain yield - g/plot|CO_321:0001221",
    "Grain yield - main tillers - kg/ha|CO_321:0501088",
]

available = [c for c in yield_cols if c in hist.columns]
if not available:
    raise SystemExit("No grain yield columns found in historical phenotype file.")

yield_col = (
    "Grain yield - kg/ha|CO_321:0001218"
    if "Grain yield - kg/ha|CO_321:0001218" in available
    else available[0]
)

hist["value"] = pd.to_numeric(hist[yield_col], errors="coerce")
hist = hist.dropna(subset=["value"]).copy()


# Build trial column
if "studyName" not in hist.columns:
    raise SystemExit("Historical phenotype missing 'studyName' column.")
hist["trial"] = hist["studyName"]


# Ensure germplasm mapping column exists
if "germplasmName_mapped" not in hist.columns:
    hist["germplasmName_mapped"] = hist["germplasmName"]

hist = hist[["germplasmName", "germplasmName_mapped", "value", "trial"]]


# Load Predictathon training phenotypes
print("[unified_pheno] Loading Predictathon training phenotypes...")

predictathon_frames = []
for trial in sorted(os.listdir(PREDICTATHON_DIR)):
    trial_dir = f"{PREDICTATHON_DIR}/{trial}"
    pheno_path = f"{trial_dir}/training_pheno_merged.csv"

    if not os.path.exists(pheno_path):
        continue

    df = pd.read_csv(
        pheno_path,
        encoding="latin1",
        engine="python",
        dtype=str,
        on_bad_lines="skip"
    )

    df["trial"] = trial

    # Normalize germplasm column name
    if "germplasmName" not in df.columns:
        if "germplasm_name" in df.columns:
            df = df.rename(columns={"germplasm_name": "germplasmName"})
        else:
            raise SystemExit(f"Predictathon file {pheno_path} missing germplasm name column.")

    # Ensure mapping column exists
    if "germplasmName_mapped" not in df.columns:
        df["germplasmName_mapped"] = df["germplasmName"]

    # Ensure numeric phenotype
    if "value" not in df.columns:
        raise SystemExit(f"Predictathon file {pheno_path} missing 'value' column.")

    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    predictathon_frames.append(df)

if not predictathon_frames:
    raise SystemExit("No Predictathon training phenotype files found.")

pred = pd.concat(predictathon_frames, ignore_index=True)
pred = pred[["germplasmName", "germplasmName_mapped", "value", "trial"]]


# Combine historical + Predictathon
print("[unified_pheno] Combining historical + Predictathon phenotypes...")
full = pd.concat([hist, pred], ignore_index=True)


# Drop missing or zero values
full = full.dropna(subset=["value"])
full = full[full["value"] != 0]


# Save raw per‑trial mean and std BEFORE normalization
trial_stats = (
    full.groupby("trial")["value"]
        .agg(raw_mean="mean", raw_std="std")
        .reset_index()
)

trial_stats.to_csv(TRIAL_STATS_OUT, index=False)
print(f"[unified_pheno] Wrote trial stats → {TRIAL_STATS_OUT}")


# Normalize phenotype within each trial (z‑score)
def normalize(x):
    if x.std(ddof=0) == 0:
        return np.zeros(len(x))
    return (x - x.mean()) / x.std(ddof=0)

full["value"] = full.groupby("trial")["value"].transform(normalize)


# Build unified ID column for GRM matching
full["id_for_grm"] = full["germplasmName_mapped"].fillna(full["germplasmName"])


# Save output
full.to_csv(OUTFILE, index=False)

print(f"[unified_pheno] Wrote → {OUTFILE}")
print(f"[unified_pheno] Total rows: {len(full):,}")
print(f"[unified_pheno] Unique lines: {full['id_for_grm'].nunique():,}")
print("[unified_pheno] Done.")