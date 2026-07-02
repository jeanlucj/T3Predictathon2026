#!/usr/bin/env python3
"""T3 Prediction Challenge — Genomic Prediction Pipeline.

Predicts wheat grain yield (kg/ha) for 9 test trials using
GRM-based GBLUP with G×E environmental weighting and per-trial
adaptive blending.

Workflow
--------
1.  Fetch focal-trial accession lists and location via BrAPI.
2.  Load per-study training phenotypes (grain yield).
3.  Fetch historical weather; compute per-trial environmental
    similarity weights to focal trial (G×E).
4.  Compute environmentally-weighted BLUEs (trial-adjusted means
    where observations from similar environments get higher weight).
5.  Parse VCF genotype data; apply MAF and missing-rate filters.
6.  Build GRM (Genomic Relationship Matrix = XX'/p) — uses ALL
    markers without subsampling, robust when n << p.
7.  For CV0 / CV00, solve GBLUP with LOO-CV-optimised lambda.
8.  Apply per-trial adaptive blend weights (more weight on BLUE
    when training set is small).
9.  Write submission CSVs per the competition specification.

Usage
-----
    python predict.py
"""

from __future__ import annotations

import csv
import json
import logging
import os
import shutil
import time
import zipfile
import warnings
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.impute import SimpleImputer

warnings.filterwarnings("ignore")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).resolve().parent
TRAINING_DIR = BASE_DIR / "unzipped"
GENOTYPE_DIR = BASE_DIR / "Genotype Data for Prediction Trials"
OUTPUT_DIR = BASE_DIR / "submission"
YIELD_TRAIT = "Grain yield - kg/ha|CO_321:0001218"
BRAPI_BASE = "https://wheat.triticeaetoolbox.org/brapi/v1"

TRIAL_CONFIG: dict[str, dict[str, Any]] = {
    "2025_AYT_Aurora": {
        "study_id": 10673,
        "vcf": "2025-12-16_15_59_59_fileP7Mi.vcf",
    },
    "24Crk_AY2-3": {
        "study_id": 10674,
        "vcf": "2025-12-11_18_36_07_fileY6KQ.vcf",
    },
    "25_Big6_SVREC_SVREC": {
        "study_id": 10675,
        "vcf": "2025-12-12_19_42_54_file8tAQ.vcf",
    },
    "CornellMaster_2025_McGowan": {
        "study_id": 10676,
        "vcf": "2025-12-16_15_37_13_filegz0x.vcf",
    },
    "YT_Urb_25": {
        "study_id": 10677,
        "vcf": "2026-01-12_15_02_24_file7kBP.vcf",
    },
    "AWY1_DVPWA_2024": {
        "study_id": 10678,
        "vcf": "2026-02-25_15_45_12_fileGxN4.vcf",
    },
    "OHRWW_2025_SPO": {
        "study_id": 10679,
        "vcf": "2025-12-16_14_44_02_fileoUdn.vcf",
    },
    "TCAP_2025_MANKS": {
        "study_id": 10680,
        "vcf": "2025-12-15_16_43_21_fileY4ZA.vcf",
    },
    "STP1_2025_MCG": {
        "study_id": 10681,
        "vcf": "2025-12-16_19_25_20_filek_yE.vcf",
    },
}

WEATHER_CACHE_DIR = BASE_DIR / "weather_cache"
OPEN_METEO_URL = "https://archive-api.open-meteo.com/v1/archive"
WHEAT_GDD_BASE = 0.0
MIN_ENV_TRIALS = 10
MIN_ENV_R2 = 0.05


# ──────────────────────────────────────────────────────────────────
# BrAPI helpers
# ──────────────────────────────────────────────────────────────────


def fetch_focal_trial_accessions(study_id: int) -> set[str]:
    """Return the set of germplasm names in the focal trial via BrAPI.

    Args:
        study_id: T3 study database ID.

    Returns:
        Set of germplasm names.  Empty on network/parse failure.
    """
    accessions: set[str] = set()
    page = 0
    page_size = 1000
    while True:
        url = (
            f"{BRAPI_BASE}/studies/{study_id}/germplasm"
            f"?pageSize={page_size}&page={page}"
        )
        try:
            with urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read().decode())
        except (URLError, json.JSONDecodeError, OSError) as exc:
            logger.warning("BrAPI call failed for study %s: %s", study_id, exc)
            break

        entries = data.get("result", {}).get("data", [])
        if not entries:
            break
        for entry in entries:
            name = entry.get("germplasmName", "")
            if name:
                accessions.add(name)

        pagination = data.get("metadata", {}).get("pagination", {})
        total_pages = pagination.get("totalPages", 1)
        page += 1
        if page >= total_pages:
            break

    logger.info(
        "  BrAPI: fetched %d accessions for study %d", len(accessions), study_id
    )
    return accessions


# ──────────────────────────────────────────────────────────────────
# Phenotype loading
# ──────────────────────────────────────────────────────────────────


def load_training_phenotypes(study_dir: Path) -> pd.DataFrame:
    """Load and concatenate all phenotype CSVs in *study_dir*.

    Filters for the yield trait and returns a DataFrame with columns:
        study_db_id, germplasm_name, value   (value = yield kg/ha)

    Args:
        study_dir: Directory containing phenotype CSV files.

    Returns:
        DataFrame of yield observations.
    """
    frames: list[pd.DataFrame] = []
    for csv_path in sorted(study_dir.glob("*.csv")):
        if "metadata" in csv_path.name or "vcf_names" in csv_path.name:
            continue
        try:
            df = pd.read_csv(csv_path, dtype=str)
        except Exception:
            continue

        if "observation_variable_name" not in df.columns:
            continue

        yield_df = df[df["observation_variable_name"] == YIELD_TRAIT].copy()
        if yield_df.empty:
            continue

        yield_df["value"] = pd.to_numeric(yield_df["value"], errors="coerce")
        yield_df = yield_df.dropna(subset=["value"])
        if yield_df.empty:
            continue

        frames.append(
            yield_df[["study_db_id", "germplasm_name", "value"]].copy()
        )

    if not frames:
        return pd.DataFrame(columns=["study_db_id", "germplasm_name", "value"])

    result = pd.concat(frames, ignore_index=True)
    logger.info("  Loaded %d yield observations from %s", len(result), study_dir.name)
    return result


def load_training_metadata(study_dir: Path) -> pd.DataFrame:
    """Load training trial metadata CSV.

    Args:
        study_dir: Directory containing the metadata CSV.

    Returns:
        DataFrame of training trial metadata.
    """
    for csv_path in study_dir.glob("training_trial_metadata_*.csv"):
        try:
            return pd.read_csv(csv_path, dtype=str)
        except Exception:
            continue
    return pd.DataFrame()


def load_training_vcf_urls(study_dir: Path) -> list[dict[str, str]]:
    """Load training VCF download URLs.

    Args:
        study_dir: Directory containing the training_vcf_names CSV.

    Returns:
        List of dicts with project_id, file_name, download_url.
    """
    for csv_path in study_dir.glob("training_vcf_names_*.csv"):
        try:
            df = pd.read_csv(csv_path, dtype=str)
            return df.to_dict("records")
        except Exception:
            continue
    return []


def download_training_vcf(
    url: str,
    dest_path: Path,
    timeout: int = 30,
    max_bytes: int = 200_000_000,
    max_seconds: float = 90.0,
) -> bool:
    """Download a VCF from T3 to a local path with streaming and guards.

    Args:
        url: Download URL.
        dest_path: Where to save.
        timeout: Socket-level timeout in seconds.
        max_bytes: Abort if response exceeds this size.
        max_seconds: Abort if total download exceeds this duration.

    Returns:
        True on success.
    """
    if dest_path.exists() and dest_path.stat().st_size > 100:
        return True
    t0 = time.time()
    try:
        logger.info("    Downloading training VCF → %s …", dest_path.name)
        with urlopen(url, timeout=timeout) as resp:
            chunks: list[bytes] = []
            total = 0
            while True:
                chunk = resp.read(1_048_576)  # 1 MB
                if not chunk:
                    break
                total += len(chunk)
                if total > max_bytes:
                    logger.warning("    File too large (>%d MB) — skip", max_bytes // 1_000_000)
                    return False
                if time.time() - t0 > max_seconds:
                    logger.warning("    Download too slow (>%.0fs) — skip", max_seconds)
                    return False
                chunks.append(chunk)
            dest_path.write_bytes(b"".join(chunks))
        elapsed = time.time() - t0
        logger.info("    Downloaded %.1f MB in %.0fs", total / 1e6, elapsed)
        return True
    except Exception as exc:
        logger.warning("    Download failed (%.0fs): %s", time.time() - t0, exc)
        return False


def merge_vcfs(
    primary_samples: list[str],
    primary_markers: list[str],
    primary_geno: np.ndarray,
    secondary_samples: list[str],
    secondary_markers: list[str],
    secondary_geno: np.ndarray,
) -> tuple[list[str], list[str], np.ndarray]:
    """Merge two VCF genotype matrices on shared markers.

    The primary VCF samples come first, followed by secondary samples
    NOT already present in the primary.

    Args:
        primary_samples: Sample names from primary VCF.
        primary_markers: Marker IDs from primary VCF.
        primary_geno: Genotype matrix (n_primary × n_markers_primary).
        secondary_samples: Sample names from secondary VCF.
        secondary_markers: Marker IDs from secondary VCF.
        secondary_geno: Genotype matrix (n_secondary × n_markers_secondary).

    Returns:
        Merged (samples, markers, geno_matrix) on common markers.
    """
    primary_marker_set = set(primary_markers)
    common = [m for m in secondary_markers if m in primary_marker_set]
    if len(common) < 50:
        logger.warning("    Only %d common markers — skip merge", len(common))
        return primary_samples, primary_markers, primary_geno

    primary_idx = {m: i for i, m in enumerate(primary_markers)}
    secondary_idx = {m: i for i, m in enumerate(secondary_markers)}
    cols_primary = [primary_idx[m] for m in common]
    cols_secondary = [secondary_idx[m] for m in common]

    new_primary_set = set(primary_samples)
    extra_samples = [s for s in secondary_samples if s not in new_primary_set]
    if not extra_samples:
        return primary_samples, common, primary_geno[:, cols_primary]

    extra_idx = [secondary_samples.index(s) for s in extra_samples]

    merged_geno = np.vstack([
        primary_geno[:, cols_primary],
        secondary_geno[np.array(extra_idx)][:, cols_secondary],
    ])
    merged_samples = primary_samples + extra_samples

    logger.info(
        "    Merged: %d common markers, %d extra samples → total %d",
        len(common),
        len(extra_samples),
        len(merged_samples),
    )
    return merged_samples, common, merged_geno


# ──────────────────────────────────────────────────────────────────
# VCF parsing
# ──────────────────────────────────────────────────────────────────

_GT_MAP = {
    "0/0": 0.0, "0/1": 1.0, "1/0": 1.0, "1/1": 2.0, "./.": np.nan,
    "0|0": 0.0, "0|1": 1.0, "1|0": 1.0, "1|1": 2.0, ".|.": np.nan,
}


def parse_vcf(vcf_path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Parse a VCF file into a numeric genotype matrix.

    Encoding: 0/0 → 0, 0/1 → 1, 1/1 → 2, ./. → NaN.
    Handles both unphased (/) and phased (|) genotypes,
    as well as rows with fewer fields than expected.

    Args:
        vcf_path: Path to the VCF file.

    Returns:
        Tuple of (sample_names, marker_ids, genotype_matrix)
        where genotype_matrix has shape (n_samples, n_markers).
    """
    sample_names: list[str] = []
    marker_ids: list[str] = []
    genotype_rows: list[np.ndarray] = []
    n_samples = 0

    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                sample_names = cols[9:]
                n_samples = len(sample_names)
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, snp_id = parts[0], parts[1], parts[2]
            marker_id = snp_id if snp_id != "." else f"{chrom}_{pos}"
            marker_ids.append(marker_id)

            fmt_fields = parts[8].split(":")
            gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0

            row = np.full(n_samples, np.nan, dtype=np.float32)
            gt_fields = parts[9:]
            for i, sample_gt in enumerate(gt_fields):
                if i >= n_samples:
                    break
                gt_str = sample_gt.split(":")[gt_idx]
                row[i] = _GT_MAP.get(gt_str, np.nan)
            genotype_rows.append(row)

    geno_matrix = np.stack(genotype_rows, axis=0).T  # (n_samples, n_markers)
    logger.info(
        "  VCF parsed: %d samples × %d markers from %s",
        geno_matrix.shape[0],
        geno_matrix.shape[1],
        vcf_path.name,
    )
    return sample_names, marker_ids, geno_matrix


# ──────────────────────────────────────────────────────────────────
# Weather integration (Open-Meteo Historical API)
# ──────────────────────────────────────────────────────────────────


def _estimate_season_end(start_str: str) -> str:
    """Estimate harvest date from planting date for wheat.

    Falls back to today − 5 days if the estimated end is in the future.

    Args:
        start_str: Planting/start date as YYYY-MM-DD string.

    Returns:
        Estimated end date as YYYY-MM-DD string.
    """
    dt = datetime.strptime(start_str[:10], "%Y-%m-%d")
    if dt.month >= 8:
        end = dt + timedelta(days=280)
    else:
        end = dt + timedelta(days=120)
    today = datetime.now() - timedelta(days=5)
    if end > today:
        end = today
    return end.strftime("%Y-%m-%d")


def fetch_weather(
    lat: float,
    lon: float,
    start_date: str,
    end_date: str,
) -> pd.DataFrame | None:
    """Fetch daily weather from Open-Meteo Historical API.

    Results are cached to disk under ``WEATHER_CACHE_DIR`` so that
    subsequent runs avoid redundant network calls.

    Args:
        lat: Latitude in decimal degrees.
        lon: Longitude in decimal degrees.
        start_date: Start date (YYYY-MM-DD).
        end_date: End date (YYYY-MM-DD).

    Returns:
        DataFrame with daily weather columns, or None on failure.
    """
    WEATHER_CACHE_DIR.mkdir(exist_ok=True)
    cache_key = f"{lat:.4f}_{lon:.4f}_{start_date}_{end_date}".replace("-", "")
    cache_path = WEATHER_CACHE_DIR / f"{cache_key}.json"

    if cache_path.exists():
        try:
            with open(cache_path) as f:
                data = json.load(f)
            daily = data.get("daily", {})
            return pd.DataFrame(daily) if daily else None
        except Exception:
            pass

    params = (
        f"latitude={lat}&longitude={lon}"
        f"&start_date={start_date}&end_date={end_date}"
        f"&daily=temperature_2m_mean,temperature_2m_max,"
        f"temperature_2m_min,precipitation_sum"
        f"&timezone=auto"
    )
    url = f"{OPEN_METEO_URL}?{params}"

    for attempt in range(3):
        try:
            time.sleep(0.25 * (attempt + 1))
            with urlopen(url, timeout=20) as resp:
                data = json.loads(resp.read().decode())
            with open(cache_path, "w") as f:
                json.dump(data, f)
            daily = data.get("daily", {})
            if not daily or "time" not in daily:
                return None
            return pd.DataFrame(daily)
        except Exception as exc:
            if attempt == 2:
                logger.debug(
                    "Weather fetch failed (%.2f, %.2f): %s", lat, lon, exc
                )
    return None


def compute_weather_features(daily_df: pd.DataFrame) -> dict[str, float]:
    """Derive agronomic weather features from daily observations.

    Features computed:
        temp_mean, temp_std  – growing-season temperature statistics.
        gdd                  – Growing Degree Days (base 0 °C).
        heat_days            – days with T_max > 32 °C.
        frost_days           – days with T_min < −2 °C.
        precip_total         – total precipitation (mm).
        precip_days          – days with > 1 mm precipitation.

    Args:
        daily_df: DataFrame from Open-Meteo with daily weather columns.

    Returns:
        Dict of feature name → numeric value.
    """
    features: dict[str, float] = {}

    t_mean = pd.to_numeric(
        daily_df.get("temperature_2m_mean"), errors="coerce"
    ).dropna()
    t_max = pd.to_numeric(
        daily_df.get("temperature_2m_max"), errors="coerce"
    ).dropna()
    t_min = pd.to_numeric(
        daily_df.get("temperature_2m_min"), errors="coerce"
    ).dropna()
    precip = pd.to_numeric(
        daily_df.get("precipitation_sum"), errors="coerce"
    ).dropna()

    if not t_mean.empty:
        features["temp_mean"] = float(t_mean.mean())
        features["temp_std"] = float(t_mean.std())
        features["gdd"] = float(t_mean.clip(lower=WHEAT_GDD_BASE).sum())

    if not t_max.empty:
        features["heat_days"] = float((t_max > 32).sum())

    if not t_min.empty:
        features["frost_days"] = float((t_min < -2).sum())

    if not precip.empty:
        features["precip_total"] = float(precip.sum())
        features["precip_days"] = float((precip > 1.0).sum())

    return features


def fetch_focal_study_location(study_id: int) -> dict[str, Any]:
    """Retrieve focal trial location and dates via BrAPI.

    Args:
        study_id: T3 study database ID.

    Returns:
        Dict with latitude, longitude, start_date, end_date keys.
        Values are None when the field is unavailable.
    """
    url = f"{BRAPI_BASE}/studies/{study_id}"
    for attempt in range(3):
        try:
            time.sleep(0.5 * attempt)
            with urlopen(url, timeout=20) as resp:
                data = json.loads(resp.read().decode())
            result = data.get("result", {})
            location = result.get("location", {})
            lat = location.get("latitude") or result.get("latitude")
            lon = location.get("longitude") or result.get("longitude")
            start = result.get("startDate")
            end = result.get("endDate")
            logger.info(
                "  Focal location: (%.4f, %.4f), start=%s",
                float(lat) if lat else 0,
                float(lon) if lon else 0,
                start or "NA",
            )
            return {
                "latitude": float(lat) if lat is not None else None,
                "longitude": float(lon) if lon is not None else None,
                "start_date": start,
                "end_date": end,
            }
        except Exception as exc:
            if attempt == 2:
                logger.warning(
                    "  BrAPI study details failed for %d: %s", study_id, exc
                )
    return {}


def build_env_weights(
    meta_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    focal_study_id: int,
    focal_lat: float | None,
    focal_lon: float | None,
    focal_start: str | None,
    focal_end: str | None,
) -> dict[int, float]:
    """Compute per-trial environmental similarity weights for G×E.

    Fetches weather for each training trial and the focal trial,
    then computes Gaussian similarity in the weather feature space.
    Trials with environments similar to the focal trial receive
    higher weight in downstream BLUE computation.

    Args:
        meta_df: Training trial metadata with lat/lon/dates.
        pheno_df: Training phenotype observations.
        focal_study_id: Focal trial study DB ID.
        focal_lat: Focal trial latitude.
        focal_lon: Focal trial longitude.
        focal_start: Focal trial start date.
        focal_end: Focal trial end date.

    Returns:
        Dict mapping study_db_id → similarity weight (0, 1].
        Empty dict when weather data is insufficient.
    """
    if focal_lat is None or focal_lon is None or focal_start is None:
        logger.info("    No focal trial location — skip G×E weighting")
        return {}

    if meta_df.empty:
        return {}

    sid_list: list[int] = []
    feature_rows: list[dict[str, float]] = []
    seen_keys: set[str] = set()
    fetch_count = 0

    logger.info("    Fetching weather for G×E environmental weights …")

    for _, row in meta_df.iterrows():
        sid = row.get("study_db_id")
        if sid is None:
            continue
        try:
            sid_int = int(float(sid))
        except (ValueError, TypeError):
            continue
        if sid_int == focal_study_id:
            continue

        lat = pd.to_numeric(row.get("latitude"), errors="coerce")
        lon = pd.to_numeric(row.get("longitude"), errors="coerce")
        if pd.isna(lat) or pd.isna(lon):
            continue

        start = str(row.get("start_date", ""))
        if start in ("", "NA", "nan", "None", "NaT"):
            continue
        start_clean = start[:10]

        end = str(row.get("end_date", ""))
        if end in ("", "NA", "nan", "None", "NaT"):
            end_clean = _estimate_season_end(start)
        else:
            end_clean = end[:10]

        loc_key = f"{lat:.3f}_{lon:.3f}_{start_clean[:7]}"
        if loc_key in seen_keys:
            continue
        seen_keys.add(loc_key)

        daily = fetch_weather(float(lat), float(lon), start_clean, end_clean)
        fetch_count += 1
        if fetch_count % 20 == 0:
            logger.info("    … fetched %d locations so far", fetch_count)
        if daily is None or daily.empty:
            continue

        feats = compute_weather_features(daily)
        if len(feats) < 5:
            continue

        sid_list.append(sid_int)
        feature_rows.append(feats)

    logger.info(
        "    Weather collected for %d / %d unique trial locations",
        len(feature_rows),
        len(seen_keys),
    )

    if len(feature_rows) < MIN_ENV_TRIALS:
        logger.info("    Too few trials (%d) for G×E weighting", len(feature_rows))
        return {}

    # ── Focal trial weather features ──
    focal_start_clean = focal_start[:10] if focal_start else ""
    if not focal_start_clean:
        return {}

    if focal_end and focal_end not in ("NA", "nan", "", "None"):
        focal_end_clean = focal_end[:10]
    else:
        focal_end_clean = _estimate_season_end(focal_start)

    focal_daily = fetch_weather(
        focal_lat, focal_lon, focal_start_clean, focal_end_clean
    )
    if focal_daily is None or focal_daily.empty:
        logger.info("    No weather for focal trial — skip G×E weighting")
        return {}

    focal_feats = compute_weather_features(focal_daily)
    if len(focal_feats) < 5:
        return {}

    # ── Build feature matrix and compute similarity ──
    feat_df = pd.DataFrame(feature_rows).fillna(0)
    focal_row = pd.DataFrame([focal_feats]).reindex(
        columns=feat_df.columns, fill_value=0
    )

    X_all = np.vstack([feat_df.values, focal_row.values])
    col_std = X_all.std(axis=0)
    col_std[col_std < 1e-8] = 1.0
    X_scaled = X_all / col_std

    X_train_scaled = X_scaled[:-1]
    X_focal_scaled = X_scaled[-1:].reshape(1, -1)

    dists = cdist(X_focal_scaled, X_train_scaled, metric="euclidean").flatten()

    median_dist = float(np.median(dists)) if len(dists) > 0 else 1.0
    sigma = max(median_dist, 1e-4)
    weights = np.exp(-0.5 * (dists / sigma) ** 2)

    env_w: dict[int, float] = {}
    for sid_int, w in zip(sid_list, weights):
        env_w[sid_int] = float(w)

    n_high = sum(1 for w in weights if w > 0.5)
    logger.info(
        "    G×E weights: %d trials, %d with w>0.5, median_dist=%.2f",
        len(weights), n_high, median_dist,
    )
    return env_w


# ──────────────────────────────────────────────────────────────────
# Prediction model — GRM-based GBLUP
# ──────────────────────────────────────────────────────────────────

MIN_TRAINING_SIZE = 3
MIN_MAF = 0.01
MAX_MISSING_RATE = 0.5

TRIAL_PARAMS: dict[str, dict[str, float]] = {
    "2025_AYT_Aurora":          {"blend_obs": 0.95},
    "CornellMaster_2025_McGowan": {"blend_obs": 0.85},
    "AWY1_DVPWA_2024":          {"blend_obs": 0.95},
    "OHRWW_2025_SPO":           {"blend_obs": 0.85},
    "STP1_2025_MCG":            {"blend_obs": 0.80},
    "YT_Urb_25":                {"blend_obs": 0.50},
    "TCAP_2025_MANKS":          {"blend_obs": 0.60},
    "25_Big6_SVREC_SVREC":      {"blend_obs": 0.65},
    "24Crk_AY2-3":              {"blend_obs": 0.70},
}


def _optimal_lambda_loo(G_tt: np.ndarray, y_c: np.ndarray) -> float:
    """Choose GBLUP regularisation lambda via leave-one-out CV.

    Uses the hat-matrix identity: LOO residual_i = r_i / (1 - h_ii).

    Args:
        G_tt: GRM sub-matrix for training samples (n × n).
        y_c: Centred training phenotypes.

    Returns:
        Best lambda from a log-spaced grid.
    """
    n = G_tt.shape[0]
    lambdas = np.logspace(-4, 4, 30)
    best_mse = np.inf
    best_lam = 1.0
    for lam in lambdas:
        try:
            L_inv = np.linalg.inv(G_tt + lam * np.eye(n))
            H = G_tt @ L_inv
            h_diag = np.diag(H)
            residuals = y_c - H @ y_c
            denom = 1.0 - h_diag
            denom = np.where(np.abs(denom) < 1e-10, 1e-10, denom)
            loo_mse = float(np.mean((residuals / denom) ** 2))
            if loo_mse < best_mse:
                best_mse = loo_mse
                best_lam = lam
        except np.linalg.LinAlgError:
            continue
    return best_lam


def build_and_predict(
    geno_matrix: np.ndarray,
    sample_names: list[str],
    training_yields: dict[str, float],
    all_training_yields: dict[str, float] | None = None,
    blend_obs: float = 0.70,
) -> dict[str, float]:
    """Fit GBLUP via the Genomic Relationship Matrix and predict all samples.

    Uses ALL markers through the n×n GRM (no marker subsampling),
    which is standard in plant breeding and robust when n << p.

    Args:
        geno_matrix: Genotype matrix (n_samples × n_markers).
        sample_names: Sample names aligned with geno_matrix rows.
        training_yields: Accession name → BLUE yield (VCF-overlapping).
        all_training_yields: Full training-yield dict for global mean.
        blend_obs: Weight on observed BLUE for training accessions.

    Returns:
        Dict mapping every sample name → predicted yield.
    """
    if all_training_yields is None:
        all_training_yields = training_yields

    global_mean = (
        float(np.mean(list(all_training_yields.values())))
        if all_training_yields
        else 4000.0
    )

    name_to_idx = {name: i for i, name in enumerate(sample_names)}
    train_names = [n for n in sample_names if n in training_yields]

    if len(train_names) < MIN_TRAINING_SIZE:
        logger.warning(
            "    Only %d training accessions in VCF — mean fallback (%.1f)",
            len(train_names),
            global_mean,
        )
        return {name: global_mean for name in sample_names}

    train_idx = np.array([name_to_idx[n] for n in train_names])
    y_train = np.array([training_yields[n] for n in train_names])

    X_full = geno_matrix.copy()

    # ── MAF filtering ──
    allele_freq = np.nanmean(X_full, axis=0) / 2.0
    maf = np.minimum(allele_freq, 1.0 - allele_freq)
    maf_valid = np.isfinite(maf)
    maf_pass = maf_valid & (maf > MIN_MAF)
    if maf_pass.sum() >= 10:
        X_full = X_full[:, maf_pass]

    # ── Missing-rate filtering ──
    missing_rate = np.isnan(X_full).mean(axis=0)
    miss_pass = missing_rate < MAX_MISSING_RATE
    if miss_pass.sum() >= 10:
        X_full = X_full[:, miss_pass]

    # ── Mean imputation ──
    imputer = SimpleImputer(strategy="mean")
    X_full = imputer.fit_transform(X_full)

    # ── Variance filtering ──
    col_var = np.var(X_full, axis=0)
    keep_mask = col_var > 1e-8
    if keep_mask.sum() < 10:
        logger.warning("    Too few variable markers — mean fallback (%.1f)", global_mean)
        return {name: global_mean for name in sample_names}
    X_full = X_full[:, keep_mask]

    n_samples, n_markers = X_full.shape

    # ── Centre genotypes (column means) — standard for GRM ──
    col_means = X_full.mean(axis=0)
    X_centered = X_full - col_means

    # ── Build Genomic Relationship Matrix (GRM = XX'/p) ──
    G = X_centered @ X_centered.T / n_markers
    G += np.eye(n_samples) * 1e-4  # numerical stability

    G_tt = G[np.ix_(train_idx, train_idx)]
    G_all_t = G[:, train_idx]

    mu = float(y_train.mean())
    y_c = y_train - mu

    # ── Optimise lambda via LOO-CV ──
    lam = _optimal_lambda_loo(G_tt, y_c)

    # ── Solve GBLUP ──
    n_train = len(train_idx)
    L = G_tt + lam * np.eye(n_train)
    try:
        alpha = np.linalg.solve(L, y_c)
    except np.linalg.LinAlgError:
        alpha = np.linalg.lstsq(L, y_c, rcond=None)[0]

    y_pred_all = mu + G_all_t @ alpha

    logger.info(
        "    GBLUP: λ=%.4f, n_train=%d, p=%d, n_samples=%d",
        lam, n_train, n_markers, n_samples,
    )

    predictions: dict[str, float] = {}
    for i, name in enumerate(sample_names):
        predictions[name] = float(y_pred_all[i])

    # ── Blend with observed BLUE (adaptive per trial) ──
    w = blend_obs
    for name in train_names:
        obs = training_yields[name]
        pred = predictions[name]
        predictions[name] = w * obs + (1.0 - w) * pred

    return predictions


# ──────────────────────────────────────────────────────────────────
# Main prediction loop
# ──────────────────────────────────────────────────────────────────


def compute_accession_blues(
    pheno_df: pd.DataFrame,
    focal_study_id: int,
    focal_accessions: set[str],
    cv_mode: str,
    env_weights: dict[int, float] | None = None,
) -> dict[str, float]:
    """Compute BLUEs (trial-adjusted means) per accession.

    Removes trial effects before averaging.  When *env_weights* are
    provided (G×E environmental similarity), observations from
    environmentally similar trials receive higher weight.

        adjusted_ij = y_ij − trial_mean_j + grand_mean
        BLUE_i       = weighted_mean(adjusted_ij, w_j)

    Args:
        pheno_df: All training phenotype observations.
        focal_study_id: Study DB ID of the focal trial.
        focal_accessions: Accession names in the focal trial.
        cv_mode: "CV0" or "CV00".
        env_weights: Optional per-study-id environmental similarity
            weights (higher = more similar to focal trial).

    Returns:
        Dict mapping accession name → BLUE yield.
    """
    df = pheno_df.copy()
    df = df[df["study_db_id"].astype(int) != focal_study_id]

    if cv_mode == "CV00":
        df = df[~df["germplasm_name"].isin(focal_accessions)]

    if df.empty:
        return {}

    grand_mean = df["value"].mean()
    trial_effects = df.groupby("study_db_id")["value"].transform("mean")
    df["adjusted"] = df["value"] - trial_effects + grand_mean

    if env_weights:
        df["env_w"] = df["study_db_id"].astype(int).map(env_weights).fillna(1.0)
        df["env_w"] = df["env_w"].clip(lower=0.1)

        def _weighted_mean(g: pd.DataFrame) -> float:
            w = g["env_w"].values
            v = g["adjusted"].values
            return float(np.average(v, weights=w))

        blues = df.groupby("germplasm_name").apply(_weighted_mean).to_dict()
    else:
        blues = df.groupby("germplasm_name")["adjusted"].mean().to_dict()

    return blues


def run_trial(
    trial_name: str,
    config: dict[str, Any],
    output_base: Path,
) -> None:
    """Run CV0 and CV00 predictions for a single trial.

    Args:
        trial_name: Name of the prediction trial.
        config: Configuration dict with study_id and vcf filename.
        output_base: Base output directory.
    """
    study_id = config["study_id"]
    vcf_file = config["vcf"]
    study_dir = TRAINING_DIR / f"study{study_id}"
    vcf_path = GENOTYPE_DIR / vcf_file

    logger.info("=" * 60)
    logger.info("Trial: %s  (study %d)", trial_name, study_id)
    logger.info("=" * 60)

    # ── 1. Fetch focal-trial accessions via BrAPI ──
    focal_accessions = fetch_focal_trial_accessions(study_id)

    # ── 2. Load training phenotypes ──
    pheno_df = load_training_phenotypes(study_dir)
    if pheno_df.empty:
        logger.error("  No training phenotype data for %s — skipping", trial_name)
        return

    # ── 3. Load training metadata (for output files) ──
    meta_df = load_training_metadata(study_dir)
    training_study_names: list[str] = []
    if not meta_df.empty and "study_name" in meta_df.columns:
        training_study_names = meta_df["study_name"].dropna().unique().tolist()

    # ── 3b. Fetch focal trial location for G×E environmental weighting ──
    focal_loc = fetch_focal_study_location(study_id)

    # ── 3c. Compute per-trial environmental similarity weights (G×E) ──
    env_weights = build_env_weights(
        meta_df, pheno_df, study_id,
        focal_loc.get("latitude"),
        focal_loc.get("longitude"),
        focal_loc.get("start_date"),
        focal_loc.get("end_date"),
    )

    # ── Per-trial blend weight ──
    trial_blend = TRIAL_PARAMS.get(trial_name, {}).get("blend_obs", 0.70)

    # ── 4. Parse VCF ──
    if not vcf_path.exists():
        logger.error("  VCF not found: %s — skipping", vcf_path)
        return
    sample_names, marker_ids, geno_matrix = parse_vcf(vcf_path)

    # ── 4b. Merge ALL training VCFs to maximise genotyped training set ──
    all_pheno_accessions = set(pheno_df["germplasm_name"].unique())
    vcf_set = set(sample_names)
    overlap_count = len(all_pheno_accessions & vcf_set)
    logger.info("  Phenotype-VCF overlap (focal VCF only): %d accessions", overlap_count)

    vcf_refs = load_training_vcf_urls(study_dir)
    MERGE_TIME_LIMIT = 300  # 5 minutes max for VCF downloads per trial
    OVERLAP_TARGET = 500
    CONSEC_FAIL_LIMIT = 4   # stop after N consecutive incompatible-platform VCFs
    if vcf_refs:
        logger.info(
            "  Merging training VCFs to expand genotyped training set (%d available) …",
            len(vcf_refs),
        )
        training_vcf_cache = BASE_DIR / "training_vcf_cache"
        training_vcf_cache.mkdir(exist_ok=True)
        merge_start = time.time()
        merged_count = 0
        consec_platform_fail = 0
        prev_n_samples = len(sample_names)

        for idx, ref in enumerate(vcf_refs):
            elapsed = time.time() - merge_start
            if elapsed > MERGE_TIME_LIMIT:
                logger.info("    Time limit (%.0fs) — stopping VCF merge", elapsed)
                break
            if consec_platform_fail >= CONSEC_FAIL_LIMIT:
                logger.info(
                    "    %d consecutive platform mismatches — stopping VCF merge",
                    consec_platform_fail,
                )
                break
            cur_overlap = len(all_pheno_accessions & set(sample_names))
            if cur_overlap >= OVERLAP_TARGET:
                logger.info("    Overlap target (%d) reached — stopping", cur_overlap)
                break

            url = ref.get("download_url", "")
            fname = ref.get("file_name", "unknown")
            safe_fname = fname.replace(":", "_").replace("/", "_") + ".vcf"
            dest = training_vcf_cache / safe_fname
            if download_training_vcf(url, dest):
                try:
                    s2, m2, g2 = parse_vcf(dest)
                    old_n = len(sample_names)
                    sample_names, marker_ids, geno_matrix = merge_vcfs(
                        sample_names, marker_ids, geno_matrix, s2, m2, g2,
                    )
                    if len(sample_names) > old_n:
                        merged_count += 1
                        consec_platform_fail = 0
                    else:
                        consec_platform_fail += 1
                except Exception as exc:
                    logger.warning("    Failed to merge VCF %s: %s", fname, exc)
                    consec_platform_fail += 1
            else:
                consec_platform_fail += 1

            if (idx + 1) % 5 == 0:
                cur_overlap = len(all_pheno_accessions & set(sample_names))
                logger.info(
                    "    Progress: %d/%d attempted, %d merged, overlap=%d (%.0fs)",
                    idx + 1, len(vcf_refs), merged_count, cur_overlap,
                    time.time() - merge_start,
                )

        final_overlap = len(all_pheno_accessions & set(sample_names))
        logger.info(
            "  VCF merge done: overlap %d → %d, samples %d → %d (%d VCFs merged, %.0fs)",
            overlap_count, final_overlap, prev_n_samples, len(sample_names),
            merged_count, time.time() - merge_start,
        )

    # ── 5. Filter to genotyped focal-trial accessions ──
    vcf_set = set(sample_names)
    if focal_accessions:
        genotyped_focal = focal_accessions & vcf_set
        logger.info(
            "  Focal accessions: %d total, %d genotyped",
            len(focal_accessions),
            len(genotyped_focal),
        )
    else:
        genotyped_focal = vcf_set
        logger.warning("  BrAPI returned no accessions — predicting ALL VCF samples")

    # ── 6. Build predictions for CV0 and CV00 ──
    trial_out = output_base / trial_name
    trial_out.mkdir(parents=True, exist_ok=True)

    for cv_mode in ("CV0", "CV00"):
        logger.info("  --- %s ---", cv_mode)

        training_yields = compute_accession_blues(
            pheno_df, study_id, focal_accessions, cv_mode,
            env_weights=env_weights,
        )
        logger.info("    Training accessions with yield: %d", len(training_yields))

        training_yields_in_vcf = {
            k: v for k, v in training_yields.items() if k in vcf_set
        }
        logger.info(
            "    Training accessions in VCF: %d", len(training_yields_in_vcf)
        )

        predictions = build_and_predict(
            geno_matrix, sample_names, training_yields_in_vcf, training_yields,
            blend_obs=trial_blend,
        )

        # Filter predictions to genotyped focal accessions
        focal_predictions = {
            name: predictions.get(name, np.nan)
            for name in sorted(genotyped_focal)
            if name in predictions
        }
        logger.info("    Final predictions: %d accessions", len(focal_predictions))

        # ── Write output CSVs ──
        pred_path = trial_out / f"{cv_mode}_Predictions.csv"
        with open(pred_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["germplasmName", "prediction"])
            for name in sorted(focal_predictions):
                writer.writerow([name, f"{focal_predictions[name]:.4f}"])

        # Determine training data that contributed to predictions.
        # For Ridge models, this is the VCF-overlapping accessions.
        # For mean fallback, it is ALL accessions whose yields
        # informed the global mean.
        filtered_df = pheno_df.copy()
        filtered_df = filtered_df[
            filtered_df["study_db_id"].astype(int) != study_id
        ]
        if cv_mode == "CV00":
            filtered_df = filtered_df[
                ~filtered_df["germplasm_name"].isin(focal_accessions)
            ]

        used_study_ids = set(filtered_df["study_db_id"].astype(int).unique())
        all_training_accessions = sorted(
            filtered_df["germplasm_name"].dropna().unique()
        )

        # Resolve study IDs → study names via metadata
        used_trial_names: list[str] = []
        if not meta_df.empty and "study_db_id" in meta_df.columns:
            meta_ids = pd.to_numeric(
                meta_df["study_db_id"], errors="coerce"
            )
            for idx, sid in meta_ids.items():
                if pd.notna(sid) and int(sid) in used_study_ids:
                    sname = meta_df.at[idx, "study_name"]
                    if sname and sname not in used_trial_names:
                        used_trial_names.append(sname)

        if not used_trial_names:
            used_trial_names = [str(sid) for sid in sorted(used_study_ids)]

        trials_path = trial_out / f"{cv_mode}_Trials.csv"
        with open(trials_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["studyName"])
            for tname in sorted(used_trial_names):
                writer.writerow([tname])

        acc_path = trial_out / f"{cv_mode}_Accessions.csv"
        with open(acc_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["germplasmName"])
            for acc in all_training_accessions:
                writer.writerow([acc])

    logger.info("  Output written to %s", trial_out)


# ──────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────


def main() -> None:
    """Orchestrate the full prediction pipeline."""
    logger.info("T3 Prediction Challenge — Genomic Prediction Pipeline")
    logger.info("Base directory: %s", BASE_DIR)

    # Clean & create output directory
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    for trial_name, config in TRIAL_CONFIG.items():
        try:
            run_trial(trial_name, config, OUTPUT_DIR)
        except Exception:
            logger.exception("Failed to process trial %s", trial_name)

    # Copy methods description (version-controlled at repo root) into output
    methods_src = BASE_DIR / "methods_description.txt"
    methods_dest = OUTPUT_DIR / "methods_description.txt"
    if methods_src.exists():
        shutil.copy2(methods_src, methods_dest)
        logger.info("Copied methods_description.txt → %s", methods_dest)
    else:
        logger.warning("methods_description.txt not found at %s", methods_src)

    # Create zip for submission
    logger.info("=" * 60)
    logger.info("Creating submission zip …")
    zip_out = BASE_DIR / "submission.zip"
    with zipfile.ZipFile(zip_out, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in OUTPUT_DIR.rglob("*"):
            if path.is_file():
                arc = path.relative_to(OUTPUT_DIR)
                zf.write(path, arc)
    logger.info("Submission archive: %s", zip_out)
    logger.info("Done.")


if __name__ == "__main__":
    main()
