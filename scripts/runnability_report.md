# Static Audit — Runnability Report

Audit date: 2026-05-20 (revised 2026-05-21 after PII removal pass)
Auditor: Claude Code (Opus 4.7) — static analysis only, no scripts were executed.

Per-folder findings: required dependencies, hard-coded paths, expected inputs,
what's actually present on disk, and what would block a fresh run on this machine.

---

## Prediction1 — R / GBLUP

**Entry points** (sequential): `00_download_training_phenotypes.R` →
`02_prepare_phenotypes_final_1.R` → `01_vcf_to_snp_matrix_final.R` →
`03_genomic_prediction.R`

(Note: the file names suggest 00→01→02→03 but script 00 reads a file produced
by 02. Script 02 downloads `Training_Trial_Info.rds` from GitHub; so the actual
order is 02 → 00 → 01 → 03.)

**R packages**
- CRAN: `tidyverse`, `vcfR`, `rrBLUP`, `AGHmatrix`, `data.table`, `janitor`, `MASS`
- GitHub: `BrAPI` (IntegratedBreedingPlatform/BrAPI),
  `T3BrapiHelpers` (jeanlucj/T3BrapiHelpers)

**Hard-coded paths — previously a blocker, now resolved**
- All four scripts previously used Windows-only absolute paths under one
  contributor's drive layout. They have been refactored to use `here::here()`
  rooted at `scripts/Prediction1/data/...`, so they run from any platform
  once the data subfolders are populated.

**Inputs required**
- `Training_Trial_Info.rds` — downloaded from GitHub by script 02 (network)
- T3/Wheat phenotypes — downloaded via BrAPI by script 00 (network + interactive login)
- VCF files (9 focal + many historical) — must be placed in
  `scripts/Prediction1/data/Genotype/` manually; **none present in this repo**

**Present locally:** None (folder contains only the 4 scripts + README + LICENSE)

**Verdict:** With paths now portable, the script set needs only
(a) BrAPI credentials and (b) all VCFs downloaded and placed under
`scripts/Prediction1/data/Genotype/` before a fresh run is possible.

---

## Prediction2 — Python pipeline

**Entry point:** `bash run_pipeline.sh` (the `Snakefile` is alternative and
references different VCF filenames; only `run_pipeline.sh` is current)

**Python packages**
- `numpy`, `pandas`, `pyyaml`, `cyvcf2`, `joblib`
- Plus standard lib (argparse, gzip, pathlib, etc.)

**Hard-coded paths:** None problematic — uses `config.yaml` and relative paths
from script dir. Portable.

**Inputs required (per `data/predictathon/<trial>/`)**
- `accessions.csv` ✅ present for all 9
- `training_pheno_merged.csv` ✅ present for all 9
- `genotypes/*.vcf` or `*.vcf.gz` — **only STP1_2025_MCG present**;
  8 of 9 trials missing
- `data/processed/unified_training_pheno_mapped.csv` ✅ present
- `data/processed/unified_training_pheno_cleaned.csv` ✅ present

**Present locally:** Phenotype layer is fully populated. Genotype VCFs missing
for 8/9 trials.

**Verdict:** This is the closest to runnable. Install ~5 Python packages and
drop in the 8 missing VCFs and you should be able to run end-to-end. Could
even do a single-trial dry run on STP1_2025_MCG today to confirm the code path.

---

## Prediction3 — R workflowr

**Entry points (in order):**
1. `Data/T3_trial_import_and_filtering.Rmd`
2. `Data/Phenotypic_data_import.Rmd`, `Data/Genomic_data_import.Rmd`,
   `Data/Environmental_data_import.Rmd` (any order, after #1)
3. `Analysis/genomic_analysis.Rmd`, `Analysis/enviromic_analysis.Rmd`
4. `Prediction Models/GBLUP_predictions.Rmd`

**R packages (many)**
- CRAN: `tidyverse`, `vcfR`, `rrBLUP`, `sommer`, `AGHmatrix`, `janitor`, `here`,
  `kableExtra`, `devtools`, `furrr`, `progressr`, `jsonlite`, `httr`, `progress`
- GitHub: `BrAPI` (TriticeaeToolbox/BrAPI.R), `CovCombR` (cran/CovCombR),
  `EnvRtype` (allogamous/EnvRtype)
- Scripts call `install_github(...)` on every run — annoying but harmless.

**Hard-coded paths:** None — uses `here::here()` consistently. Portable.

**Inputs required**
- Phenotype/metadata CSVs in `Data/` ✅ present (`genotyped_phenotyped_trials.csv`,
  `genotyped_trial_phenotypes.csv`, EnvRtype outputs, test trial metadata, etc.)
- `Analysis/trial_Gmat.rds` ✅ present
- `Analysis/Trial_W_mat.csv` ✅ present
- `Analysis/comb_grm_mat.rds` — **NOT present** (gitignored; required by
  `GBLUP_predictions.Rmd`)
- VCF files in `data/vcf_genotyped_data/` — **NOT present**; needed by
  `genomic_analysis.Rmd` and `Genomic_data_import.Rmd`

**Present locally:** All non-genotypic intermediate outputs, plus the final
`Predictions/` (CV0/CV00 already written).

**Verdict:** Re-running the final prediction step needs only `comb_grm_mat.rds`
regenerated, which itself needs the VCFs. Re-running stage 1–2 from scratch
requires VCFs + BrAPI access. The final submission outputs are already on
disk — useful for assessment without re-running.

---

## Prediction4 — Python single-script (G×E weighted GBLUP)

**Entry point:** `python predict.py` (1272-line monolith)

**Python packages**
- `numpy`, `pandas`, `scipy`, `scikit-learn`
- Standard lib (urllib, csv, json, zipfile, etc.)

**Hard-coded paths:** None problematic — uses `Path(__file__).resolve().parent`.
Portable.

**Inputs required**
- `unzipped/study<ID>/` for all 9 study IDs (phenotypes from the challenge zips)
  — **directory does not exist**
- `Genotype Data for Prediction Trials/<trial>.vcf` for all 9 — **directory
  does not exist**
- `weather_cache/` — populated on demand from Open-Meteo API (network)

**Present locally:** Only the script, README, methods description, and prior
`submission/` folder.

**Verdict:** Cleanest, most self-contained code, but requires hand-fetching all
9 training zips and all 9 focal VCFs first. Once data is in place, a single
command runs the whole pipeline.

---

## Prediction5 — R workflowr + RKHS

**Entry points (in order):**
`phenotypic_data.Rmd` → `pullGenoTrial.Rmd` → `pullGenoAssociatedTrials.Rmd` →
`combine_GRM.Rmd` → `environmental_covariates.Rmd` → `genom_pred.Rmd`

**R packages (many)**
- CRAN: `tidyverse`, `data.table`, `apercu`, `dplyr`, `ggplot2`, `paletteer`,
  `purrr`, `bigstatsr`, `ggbeeswarm`, `janitor`, `rrBLUP`, `BGLR`, `sommer`
  (implied), `foreach`, `snow`, `doSNOW`, `vcfR`, `DT`, `progress`, `here`,
  `tidygeocoder`, `nasapower`, `daymetr`, `ggVennDiagram`, `ggpubr`,
  `gtsummary`, `lme4`, `kableExtra`
- GitHub: `BrAPI` (TriticeaeToolbox/BrAPI.R),
  `T3BrapiHelpers` (jeanlucj/T3BrapiHelpers)
- External tool: Beagle 5.x (Java jar at `code/beagle.27Feb25.75f.jar` ✅ present)

**Hard-coded paths:** None — uses `here::i_am()` and `here::here()`. Portable.

**Inputs required**
- Phenotype CSVs in `data/phenotypic/<trial>/` ✅ present for all 9
- Genotype VCFs in `data/genotypic/<trial>/` — **only .csi indices and manifest
  files present**; actual `.vcf.gz` files are gitignored and missing
- NASA POWER weather (fetched on demand; network)
- BrAPI access (network)

**Present locally:** Phenotypes ✅, intermediate model RDS/TSV in `output/` ✅,
final `submission/` folders ✅. VCFs missing.

**Verdict:** Like Prediction3, the final outputs are already on disk. Re-running
from scratch needs VCFs and network access. Easiest re-run target would be
`genom_pred.Rmd` if a combined GRM `.rds` is regenerable from the existing
material — but `combine_GRM.Rmd` needs the missing VCFs.

---

## Cross-cutting observations

|                                | Pred1      | Pred2 | Pred3 | Pred4 | Pred5         |
| ------------------------------ | ---------- | ----- | ----- | ----- | ------------- |
| Portable paths                 | ✅         | ✅    | ✅    | ✅    | ✅            |
| Pheno inputs present           | ❌         | ✅    | ✅    | ❌    | ✅            |
| VCF inputs present             | ❌         | 1/9   | ❌    | ❌    | ❌            |
| Final outputs already in repo  | ❌         | ❌    | ✅    | ✅ (prior submission) | ✅ |
| Needs network/BrAPI to start   | ✅         | ❌    | ✅    | ✅    | ✅            |
| External non-package deps      | —          | —     | —     | —     | Java (Beagle) |

**The dominant blocker** across four of five folders is the absence of VCF
genotype files. For runnability, you'd need to retrieve them from T3/Wheat
once and place them in each folder's expected layout. That's the single
highest-leverage step.

**Easiest folder to actually run** is Prediction2 — it has the most complete
local data, uses portable paths, and has the fewest Python dependencies.

**Folders where the final submission already exists** in the repo (so
"assessment" may not require re-running them) are Prediction3, Prediction4
(in `submission/`), and Prediction5.

---

## Suggested next moves

1. **Inventory VCFs.** The various scripts expect specific VCF filenames per
   trial. Extract the full canonical list (which file maps to which trial
   across all 5 folders) so the files can be fetched once and reused.
2. **Dry run Prediction2 on STP1_2025_MCG** (the one trial with a VCF in place)
   to confirm the Python toolchain works end-to-end, without needing other VCFs.
3. **Skip running entirely** for Prediction3/4/5 and just assess their
   already-written submission CSVs — since "different prediction algorithms
   for different characteristics" may be measurable from outputs alone.
4. ~~Tackle Prediction1's path problem~~ — done; the four R scripts now use
   `here::here("scripts", "Prediction1", "data", ...)` and create the
   required subdirectories on first run.
