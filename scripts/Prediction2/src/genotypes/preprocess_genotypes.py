#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from cyvcf2 import VCF

from src.model.models import build_grm_from_geno


def choose_best_vcf(vcf_dir: Path) -> Path:
    """
    Choose the VCF with the largest number of samples.
    Only considers .vcf or .vcf.gz files.
    """
    if not vcf_dir.exists():
        raise SystemExit(f"Genotype directory not found: {vcf_dir}")

    vcf_files = [
        p for p in vcf_dir.iterdir()
        if p.suffix == ".vcf" or p.suffixes[-2:] == [".vcf", ".gz"]
    ]

    if not vcf_files:
        raise SystemExit(f"No VCF files found in {vcf_dir}")

    best_path = None
    best_n = -1

    for path in vcf_files:
        try:
            v = VCF(str(path))
            n = len(v.samples)
            if n > best_n:
                best_n = n
                best_path = path
        except Exception:
            continue

    if best_path is None:
        raise SystemExit(f"No valid VCF could be read in {vcf_dir}")

    print(f"[preprocess_genotypes] Selected VCF with {best_n} samples → {best_path}")
    return best_path


def vcf_to_matrix(vcf_path: Path, keep_samples: set) -> pd.DataFrame:
    """
    Convert VCF to a sample × marker dosage matrix.
    Only keeps samples in keep_samples.
    """
    vcf = VCF(str(vcf_path))
    samples = [s for s in vcf.samples if s in keep_samples]

    if not samples:
        raise SystemExit(
            f"No overlapping samples between VCF and accession list for {vcf_path}"
        )

    geno_rows = []
    marker_ids = []

    for variant in vcf:
        g = variant.genotypes

        # Convert diploid genotypes to dosages (0,1,2)
        dosages = []
        for (a1, a2, phased) in g:
            if a1 is None or a2 is None or a1 < 0 or a2 < 0:
                dosages.append(np.nan)
            else:
                dosages.append(a1 + a2)

        # Keep only the samples we want
        dosages = [dosages[vcf.samples.index(s)] for s in samples]

        # Marker ID
        if variant.ID not in (None, ".", ""):
            marker_id = variant.ID
        else:
            marker_id = f"{variant.CHROM}_{variant.POS}"

        geno_rows.append(dosages)
        marker_ids.append(marker_id)

    M = pd.DataFrame(
        np.array(geno_rows).T,
        index=samples,
        columns=marker_ids,
    )

    return M


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python -m src.genotypes.preprocess_genotypes <TRIAL>")

    trial = sys.argv[1]

    # Paths
    vcf_dir = Path("data/predictathon") / trial / "genotypes"
    acc_path = Path("data/raw/accession_lists") / f"{trial}.txt"

    if not acc_path.exists():
        raise SystemExit(f"Accession list not found: {acc_path}")

    # Load accession list
    acc_list = set(x.strip() for x in open(acc_path))

    # Choose VCF
    vcf_path = choose_best_vcf(vcf_dir)

    print(f"[preprocess_genotypes] Loading VCF for {trial}...")
    M = vcf_to_matrix(vcf_path, keep_samples=acc_list)

    # Output directory
    outdir = Path("data/predictathon") / trial / "processed"
    outdir.mkdir(parents=True, exist_ok=True)

    # Save genotype matrix
    M.to_csv(outdir / "geno_matrix.csv")
    np.save(outdir / "geno_numeric.npy", M.to_numpy())
    np.save(outdir / "geno_lines.npy", M.index.to_numpy())

    # Build GRM
    print(f"[preprocess_genotypes] Building GRM for {trial}...")
    G = build_grm_from_geno(M.to_numpy())
    np.save(outdir / "GRM.npy", G)

    print(
        f"[preprocess_genotypes] {trial}: "
        f"{M.shape[0]} samples × {M.shape[1]} markers, GRM saved."
    )


if __name__ == "__main__":
    main()