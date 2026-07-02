# T3 Prediction Challenge — Genomic Yield Prediction

Predicting wheat grain yield (kg/ha) for 9 test trials in the
[T3/Wheat Prediction Challenge](https://wheat.triticeaetoolbox.org/guides/t3-prediction-challenge)
using **GRM-based GBLUP**, **G×E-informed weighted BLUEs** (public weather data),
and **per-trial blending** of observed BLUE with genomic prediction.

## Approach

```
Training phenotypes → trial-adjusted values → G×E-weighted BLUE per accession
                                                      │
VCF (+ optional merged training VCFs) → QC → GRM → GBLUP (λ via LOO-CV) ──→ blend ──→ predictions
```

**Two-step model**

1. **Phenotype layer** — For each CV0/CV00 split, exclude the focal trial (and for CV00, exclude focal accessions elsewhere). Remove trial means, then aggregate to accession-level BLUEs with **weights** from environmental similarity to the focal trial (Gaussian kernel on scaled weather features). *No* uniform additive “environment correction” on final yields (that would not change Pearson correlation).

2. **Genomic layer** — After MAF / missing-rate filters and mean imputation, build **G = XX′/p** on all retained SNPs. Solve **GBLUP** with ridge λ on the diagonal (LOO-CV grid). Blend training accessions: `w × BLUE + (1−w) × GBLUP` with per-trial `w`.

### Weather features (for G×E weights)

Daily data from [Open-Meteo Historical API](https://open-meteo.com/); scalars include `temp_mean`, `temp_std`, `gdd`, `heat_days`, `frost_days`, `precip_total`, `precip_days`, plus `elevation` from metadata when available.

### Prediction trials

| Trial | Study ID |
|-------|----------|
| 2025_AYT_Aurora | 10673 |
| 24Crk_AY2-3 | 10674 |
| 25_Big6_SVREC_SVREC | 10675 |
| CornellMaster_2025_McGowan | 10676 |
| YT_Urb_25 | 10677 |
| AWY1_DVPWA_2024 | 10678 |
| OHRWW_2025_SPO | 10679 |
| TCAP_2025_MANKS | 10680 |
| STP1_2025_MCG | 10681 |

## Repository structure

```
.
├── predict.py                 # Full pipeline
├── methods_description.txt    # Methods text for challenge submission (source of truth)
├── submission/                # Generated CSVs + copy of methods_description.txt
├── submission.zip             # Created by predict.py (not in git; see .gitignore)
├── .gitignore
└── README.md
```

Large inputs are not tracked: `unzipped/`, `Genotype Data for Prediction Trials/`,
`weather_cache/`, `training_vcf_cache/`.

## Quick start

### Prerequisites

```bash
pip install numpy pandas scipy scikit-learn
```

### Data setup

1. Download the 9 training zips from the challenge page and unzip into `unzipped/` (e.g. `unzipped/study10673/`, …).
2. Place the 9 focal VCF files in `Genotype Data for Prediction Trials/`.

### Run

```bash
python3 predict.py
```

The pipeline writes CSVs under `submission/`, copies `methods_description.txt` into that folder, and creates **`submission.zip`** for upload.

First run can take a long time (weather API + optional training VCF downloads); cached weather and VCFs speed up later runs.

## Key design choices

- **GBLUP in GRM form** — Works in the *n×n* genomic relationship space; uses all post-QC SNPs without random subsampling.
- **G×E** — Encoded as **weights on training trials** when building BLUEs, not as a single shift of all predictions.
- **VCF merging** — Optional download/merge of training VCFs (same marker set) to increase overlapping training accessions; stops early on time limits or repeated platform mismatches.
- **Fallback** — If fewer than three training accessions overlap the VCF, predictions use the global training mean.

## Software

- Python 3.10+
- NumPy, pandas, **SciPy**, scikit-learn

