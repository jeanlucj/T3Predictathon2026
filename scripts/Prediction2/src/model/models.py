#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd


# Build grm from genotype matrix
def build_grm_from_geno(geno_numeric):
    """Build a VanRaden-like GRM from a genotype dosage matrix."""
    if hasattr(geno_numeric, "to_numpy"):
        X = geno_numeric.to_numpy(dtype=float)
    else:
        X = np.asarray(geno_numeric, dtype=float)

    n_lines, n_markers = X.shape
    if n_markers == 0:
        raise ValueError("Genotype matrix has zero markers.")

    # Impute missing values
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    if inds[0].size > 0:
        X[inds] = np.take(col_means, inds[1])

    # Remove monomorphic markers
    std = X.std(axis=0)
    keep = std > 0
    if keep.sum() == 0:
        raise ValueError("All markers are monomorphic.")
    X = X[:, keep]

    # Center markers
    X_centered = X - X.mean(axis=0, keepdims=True)
    m = X_centered.shape[1]

    return (X_centered @ X_centered.T) / m


# Legacy gblup model
class LegacyGBLUPModel:
    """Legacy GBLUP container used by fit_model() and predict_for_trial()."""
    def __init__(self, u_hat=None, lines=None):
        self.u_hat = None if u_hat is None else np.asarray(u_hat, dtype=float)
        self.lines = [] if lines is None else list(lines)

    def predict(self, G):
        if self.u_hat is None:
            raise ValueError("u_hat is not set on this LegacyGBLUPModel.")
        return self.u_hat


# Fit legacy gblup model
def fit_model(train_pheno, geno_numeric, geno_lines, G, model_type="gblup"):
    """
    Legacy GBLUP:
        u = G_all_m (G_mm + λI)^(-1) y_m
    """
    n_lines = len(geno_lines)

    if train_pheno is None or len(train_pheno) == 0:
        return LegacyGBLUPModel(u_hat=np.zeros(n_lines), lines=geno_lines)

    ph = train_pheno.copy()
    if "germplasmName" not in ph.columns:
        raise ValueError("train_pheno must contain 'germplasmName'.")

    ph = ph[ph["germplasmName"].isin(geno_lines)]
    pheno_map = dict(zip(ph["germplasmName"], ph["value"]))

    y = np.array([pheno_map.get(line, np.nan) for line in geno_lines])
    mask = ~np.isnan(y)

    if mask.sum() == 0:
        return LegacyGBLUPModel(u_hat=np.zeros(n_lines), lines=geno_lines)

    y_m = y[mask]
    G_mm = G[np.ix_(mask, mask)]
    G_all_m = G[:, mask]

    lam = 1e-5
    A = G_mm + lam * np.eye(G_mm.shape[0])
    alpha = np.linalg.solve(A, y_m)

    u_hat = G_all_m @ alpha
    return LegacyGBLUPModel(u_hat=u_hat, lines=geno_lines)


# Predict using legacy model
def predict_for_trial(model, focal_trial, test_accessions, geno_numeric, geno_lines, env, G, model_type="gblup"):
    """Return u_hat for requested accessions."""
    if isinstance(model, LegacyGBLUPModel):
        u_hat = model.u_hat
        idx = {ln: i for i, ln in enumerate(model.lines)}
        preds = [u_hat[idx.get(acc, -1)] if acc in idx else 0.0 for acc in test_accessions]
        return pd.DataFrame({"germplasmName": test_accessions, "pred": preds})

    raise TypeError("Model must be LegacyGBLUPModel for legacy prediction.")


# Simple cross validation
def cross_validate_model(train_pheno, geno_numeric, geno_lines, env, G, model_type="gblup", n_folds=5):
    if "germplasmName" not in train_pheno.columns or "value" not in train_pheno.columns:
        raise ValueError("train_pheno must contain germplasmName and value.")

    line_means = train_pheno.groupby("germplasmName")["value"].mean().reset_index()
    lines = line_means["germplasmName"].tolist()

    rng = np.random.default_rng(42)
    rng.shuffle(lines)
    folds = np.array_split(lines, n_folds)

    results = []

    for fold_id, test_lines in enumerate(folds, start=1):
        train_lines = [ln for ln in lines if ln not in test_lines]

        train_df = line_means[line_means["germplasmName"].isin(train_lines)]
        test_df = line_means[line_means["germplasmName"].isin(test_lines)]

        model = fit_model(train_df, geno_numeric, geno_lines, G)

        preds = predict_for_trial(model, None, test_lines, geno_numeric, geno_lines, None, G)
        merged = test_df.merge(preds, on="germplasmName", how="left")
        merged["fold"] = fold_id

        results.append(merged)

    return pd.concat(results, ignore_index=True)


# Global grm gblup model
class GBLUPModel:
    """New global-GRM GBLUP model."""
    def __init__(self, mu, alpha):
        self.mu = float(mu)
        self.alpha = np.asarray(alpha, dtype=float)


# Fit global grm gblup
def gblup_fit(K, y, lam=1e-1):
    """Stable ridge-regularized GBLUP."""
    mu = np.mean(y)
    y_centered = y - mu

    n = K.shape[0]
    K_reg = K + lam * np.eye(n)

    alpha, *_ = np.linalg.lstsq(K_reg, y_centered, rcond=1e-6)
    return GBLUPModel(mu=mu, alpha=alpha)


# Predict using global grm gblup
def gblup_predict(model, K_pred):
    return model.mu + K_pred @ model.alpha


# Cli: build grm for a trial
def _cli_build_grm(trial):
    ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    geno_path = f"{ROOT}/data/predictathon/{trial}/processed/geno_matrix.csv"
    outdir = f"{ROOT}/data/predictathon/{trial}/processed"
    os.makedirs(outdir, exist_ok=True)

    if not os.path.exists(geno_path):
        raise SystemExit(f"Genotype matrix not found: {geno_path}")

    print(f"[models] Loading genotype matrix for {trial}...")
    geno_df = pd.read_csv(geno_path, index_col=0)

    print(f"[models] Building GRM for {trial}...")
    G = build_grm_from_geno(geno_df)

    outpath = f"{outdir}/GRM.npy"
    np.save(outpath, G)

    print(f"[models] {trial}: GRM saved → {outpath}")
    print(f"[models] GRM shape: {G.shape}")


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "build_grm":
        _cli_build_grm(sys.argv[2])
    else:
        raise SystemExit("Usage: python -m src.model.models build_grm <TRIAL>")


if __name__ == "__main__":
    main()