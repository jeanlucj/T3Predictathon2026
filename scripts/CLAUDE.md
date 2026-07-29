# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this subdirectory.

## Project Overview

This subdirectory contains multiple competing approaches to the T3 Predictathon 2026 — a genomic prediction challenge to forecast wheat grain yield across 9 test trials using genotypic data from the Triticeae Toolbox (T3/Wheat database). Each subdirectory (`Prediction1`–`Prediction5`) represents a different team's or approach's implementation.

The challenge involves:
- **9 focal (test) trials**: AWY1_DVPWA_2024, TCAP_2025_MANKS, 25_Big6_SVREC_SVREC, OHRWW_2025_SPO, CornellMaster_2025_McGowan, 24Crk_AY2-3, 2025_AYT_Aurora, YT_Urb_25, STP1_2025_MCG
- **Two cross-validation scenarios**: CV0 (use all historical data) and CV00 (exclude focal-trial genotypes from training)
- **Submission format**: Per-trial predictions with accession IDs and predicted grain yield (kg/ha)

## Prediction Directory Summaries

| Directory | Approach | Language | Status |
|-----------|----------|----------|--------|
| **Prediction1** | GBLUP (Genomic Best Linear Unbiased Prediction) | R | Production-ready; fallback to mean imputation when GRM cannot be constructed |
| **Prediction2** | End-to-end GBLUP pipeline with global GRM | Python (+ YAML config) | Modular pipeline with caching; supports clean rebuild |
| **Prediction3** | Multi-trial similarity matrices (genomic + environmental) | R (workflowr) | Comprehensive environmental & genomic similarity analysis; Prediction Models/ contains GBLUP, MegaLMM, and Reaction Norm models |
| **Prediction4** | GRM-based GBLUP with G×E environmental weighting | Python | Environmentally-weighted BLUEs; adaptive trial-specific blending |
| **Prediction5** | Multi-method comparison (GBLUP vs RKHS) with environmental covariates | R (workflowr) | Dual-model approach; CV0 uses GBLUP, CV00 uses RKHS; includes weather/environmental data integration |

## Common Commands

Each algorithm is run by sourcing its scripts in filename order (`Prediction1/00_… → 03_…`)
or, for the Python ones, `bash Prediction2/run_pipeline.sh` / `python Prediction4/predict.py`.
Execution-order dependencies are in "Data Dependencies & Execution Order" below.

Two invocations that are NOT guessable from the file listing:

```bash
cd Prediction2 && bash run_pipeline.sh --clean   # force rebuild: drops cached GRMs/models/predictions
cd Prediction5 && Rscript -e "workflowr::wflow_build()"   # render the whole workflowr project at once
```

## Architecture Patterns

### Prediction1 & Prediction3: R-based Genomic Analysis

**Core steps**:
1. **Phenotype acquisition**: Download grain yield observations from T3/Wheat via BrAPI
2. **Data cleaning**: Outlier removal (|z| > 3), BLUE estimation per trial (y ~ germplasm + rep + block)
3. **Genotype processing**: VCF parsing, MAF ≥ 0.05, missing rate ≤ 20%, impute missing as rounded per-SNP mean
4. **GRM construction**: VanRaden method (XX'/p), ridge term for numerical stability
5. **GBLUP model**: Mixed-model using rrBLUP, conditional expectation for breeding values
6. **Fallback**: When <10 training accessions in GRM, use global mean yield

**Key R packages**: `vcfR`, `rrBLUP`, `AGHmatrix`, `BrAPI`, `tidyverse`

**Prediction3 extensions**: Builds trial-level similarity matrices (genomic and environmental) for multi-environment modeling.

### Prediction2: Python Modular Pipeline

`src/` splits into `genotypes/` (VCF preprocessing, GRM construction), `model/` (training, CV0,
CV00, expected accuracy) and `utils/`; `config.yaml` holds the focal-trial list and model settings.

**Pipeline caching**: steps are skipped when their outputs already exist -- use `--clean` to force
a rebuild. This is the behaviour that surprises people; nothing else about the layout is.

### Prediction4: Environmentally-Weighted GBLUP

**Novelty**: G×E (genotype × environment) integration
- Fetches historical weather data, computes environmental similarity to focal trial
- Weights training phenotypes by environmental similarity (environmentally-weighted BLUEs)
- Adaptive per-trial blending: upweight BLUE when training set is small

### Prediction5: Multi-Method Genomic Prediction with Environmental Covariates

**Pipeline order** (`analysis/`, workflowr): `phenotypic_data` (BrAPI phenotypes) →
`pullGenoTrial` (focal VCF → imputed numeric matrix) → `pullGenoAssociatedTrials` (related-trial
genotypes) → `combine_GRM` → `environmental_covariates` (NASA POWER weather) → `genom_pred`
(CV0/CV00 predictions). `analyze_genom_pred.Rmd` is model comparison and did NOT feed the
submission. `code/` holds the VCF helpers (`format_curate_vcf`) and a bundled Beagle jar.

**Key Features**:
- **Dual-model approach**: 
  - CV0: GBLUP with fixed effects (location + year)
  - CV00: RKHS (Reproducing Kernel Hilbert Space) with genomic + environmental effects (BRR for location/year, RKHS for genomic)
- **Trial filtering**: 
  - Selects phenotypic training data from same breeding program
  - If >20 related trials, filters to same location or year as focal trial
  - Alternative: aggregated phenotypic data from Jean-Luc Jannink's provided folder
- **Environment filtering**: Only trials with ≥30 genotypes included in cross-validation
- **Environmental integration**: Fetches NASA POWER weather variables for all trials; used as covariates
- **GRM combination**: Uses `T3BrapiHelpers::covariance_combiner()` to merge GRMs from multiple genotyping projects

**Data Processing**:
1. Phenotypes: BrAPI queries → filter to genotyped accessions
2. Genotypes: VCF → imputation (custom `format_curate_vcf`) → numeric matrix
3. Associated trial genotypes: BrAPI → identify best genotyping project → download VCF → process
4. GRM: Combine all marker matrices, compute relationship matrix per `T3BrapiHelpers`
5. Environmental data: NASA POWER queries (T2M, PRECTOTCORR, ALLSKY_SFC_SW_DWN)

**Model Selection Rationale**:
- GBLUP for CV0: simpler, proven across multiple trials
- RKHS for CV00: non-linear genomic effects + environmental covariates; added complexity justified by stricter masking

## Key Data Flows

### Input Data (Must Exist Locally)

- **T3 Phenotypes**: Grain yield (T3 variable ID: 84527) downloaded via BrAPI
- **VCF Genotypes**: From T3/Wheat archive; not stored in Git (large, gitignored)
- **Trial Metadata**: Downloaded from T3/Wheat (cultivar names, locations, harvest dates)
- **Environmental Data**: Weather variables (typically EnvRtype-derived in Prediction3)

### Intermediate Data

- **Training phenotypes**: BLUEs per accession per trial (merged across studies)
- **SNP matrices**: Encoded as -1 (homozygous ref), 0 (het), 1 (homozygous alt)
- **GRM**: Genomic Relationship Matrix (n_accessions × n_accessions); stored as `.npy` (Prediction2) or `.rds` (Prediction1/3)

### Output Format (Predictathon-compliant)

```
submission/
  <TRIAL>/
    CV0/
      predictions.csv      # accession_id, predicted_yield_kg_ha
      accessions_used.txt
      trials_used.txt
    CV00/
      predictions.csv
      accessions_used.txt
      trials_used.txt
```

## Important Considerations

### Data Dependencies & Execution Order

- **Prediction1**: Sequential execution critical (00 → 01 → 02 → 03)
- **Prediction2**: `run_pipeline.sh` handles dependencies; can run individually with manual setup
- **Prediction3**: Data/ scripts must run before Analysis/; Analysis/ must complete before Prediction_Models/
- **Prediction4**: Standalone; fetches data on-the-fly via BrAPI

### Genomic Data Quirks

- **Platform mismatch**: Some focal trials use genotyping platforms with zero SNP overlap to historical training accessions. In these cases:
  - Prediction1: Falls back to mean imputation
  - Prediction2: Training accessions not incorporated into GRM (infeasible prediction)
  - Prediction4: No historical accessions available for training
- **Accession harmonization**: Accession names may differ across VCF files; manual mapping may be required
- **Sample filtering**: Genotypes with >20% missing rate are excluded

### Cross-Validation Masking

- **CV0**: All accessions from training trials included; no accessions from focal trial genotype file
- **CV00**: Same as CV0 + exclude any accessions appearing in focal trial genotype file (even if in other historical trials)

Both scenarios submitted to evaluate robustness.

### Large Files (Gitignored)

VCFs, GRMs and trained models are gitignored (see `.gitignore` for the patterns) and must be
regenerated locally -- a fresh clone cannot reproduce results without re-running the genotype steps.

### Prediction5-Specific Notes

- **Trial filtering logic**: Phenotypic training set selection prioritizes same breeding program; within-program filtering by location/year if >20 trials available
- **RKHS model**: Implemented with `sommer` package; allows non-linear genomic effects and explicit environmental effect fitting
- **Environment-specific validation**: Only trials with ≥30 genotypes cross-validated; predictive ability calculated within and across environments
- **Beagle integration**: VCF imputation tool (`beagle.27Feb25.75f.jar`) available in `code/` for genotype phasing if needed
- **NASA POWER weather**: Used as environmental covariates; integration differs from Prediction3's EnvRtype approach

## References & External Dependencies

- **T3/Wheat Database**: wheat.triticeaetoolbox.org (BrAPI interface)
- **Non-CRAN R packages** (not installable without the source): `BrAPI`
  (GitHub: TriticeaeToolbox/BrAPI.R), `EnvRtype` (GitHub: allogamous/EnvRtype)
- **Publications**: Jarquin et al. 2017 (cross-validation scenarios), Endelman 2011 (rrBLUP)
