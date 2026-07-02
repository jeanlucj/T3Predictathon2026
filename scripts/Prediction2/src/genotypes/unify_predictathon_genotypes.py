#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

ROOT = Path("data/processed")
OUT = Path("data/processed/genotypes_predictathon.csv")

def main():
    geno_files = list(ROOT.glob("*/geno_matrix.csv"))
    print(f"Found {len(geno_files)} genotype files")

    dfs = []
    for f in geno_files:
        print(f"Loading {f}")

        # Load normally
        df = pd.read_csv(f, header=0)

        # Rename the first column (blank header) to germplasmName
        first_col = df.columns[0]
        if first_col == "" or first_col.startswith("Unnamed"):
            df = df.rename(columns={first_col: "germplasmName"})
        else:
            df = df.rename(columns={first_col: "germplasmName"})

        dfs.append(df)

    # Concatenate all trials
    full = pd.concat(dfs, ignore_index=True)

    # Drop duplicate germplasm lines
    full = full.drop_duplicates(subset=["germplasmName"])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(OUT, index=False)

    print(f"Wrote unified genotype matrix with {full.shape[0]} lines and {full.shape[1]} columns to {OUT}")

if __name__ == "__main__":
    main()