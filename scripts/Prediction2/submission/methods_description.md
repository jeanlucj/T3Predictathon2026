# T3 Predictathon — Prediction2 Methods

## Data sources

Phenotypes, genotypes, and trial metadata were obtained from the Triticeae
Toolbox (T3/Wheat) using the graphical interface. Genotypes were provided as
trial-specific VCF files, and phenotypes as CSV files.

Across the nine Predictathon trials, accession lists from
`data/raw/accession_lists/*.txt` yielded:

- 1,837 total accessions
- 1,448 genotyped
- 959 with historical phenotypes
- 866 with both genotype and phenotype
- 582 genotyped only
- 93 phenotyped only
- 296 with neither (excluded)

Counts were obtained by intersecting accession lists with the unified
phenotype file and the global genotype sample list.

## Phenotype processing

Historical phenotypes (`pheno_processed.csv`) were cleaned and merged into a
unified table for downstream analysis.

Processing steps:

- **Trait selection:** Grain yield was taken from the appropriate CO_321
  trait column, preferring kg/ha when available.
- **Cleaning:** Non-numeric, missing, and zero values were removed.
- **Trial labeling:** Assigned from `studyName` in the historical dataset.
- **Name harmonization:** `germplasmName_mapped` was used when available;
  otherwise `germplasmName`.
- **Combination:** All historical phenotypes were merged into a table
  containing `germplasmName`, `germplasmName_mapped`, `value`, and `trial`.

Per-trial raw means and standard deviations were saved to
`trial_yield_stats_raw.csv`.

**Normalization:** Phenotypes were standardized within each trial using:

$$
z_{i,t} = \frac{y_{i,t} - \bar{y}_t}{\sigma_t}
$$

Trials with zero variance were assigned zeros. A GRM-matching identifier was
created as:

```
id_for_grm = germplasmName_mapped (if present) else germplasmName
```

The final cleaned phenotype file was saved as
`unified_training_pheno_cleaned.csv`.

## Genotype processing and GRM construction

Genotypes were processed per trial using
`src/genotypes/preprocess_genotypes.py`.

**VCF selection:** For each trial, the VCF with the largest number of samples
in `data/predictathon/<TRIAL>/genotypes/` was selected.

**VCF → dosage matrix:**

- Only samples in the trial's accession list were retained.
- Diploid genotypes were converted to dosages (0/1/2), with invalid calls set
  to NaN.
- Markers were identified by ID or `CHROM_POS`.
- Outputs (`geno_matrix.csv`, `geno_numeric.npy`, `geno_lines.npy`) were
  saved to `data/predictathon/<TRIAL>/processed/`.

**Per-trial GRM:** A VanRaden-like GRM was constructed by mean-imputing
missing genotypes, removing monomorphic markers, centering markers, and
computing:

$$
G = \frac{Z Z^\top}{m}
$$

The per-trial GRM was saved as `GRM.npy`.

**Global GRM:** A global union GRM (`GRM_global_union.npy`) and sample list
(`G_global_union_samples.txt`) were built across all genotyped accessions.
This GRM was used for model training and prediction.

## Model training

Model training used a ridge-regularized GBLUP on the global GRM, implemented
in `src/model/models.py` and called by `train_global.py`.

**Name normalization:** All names were normalized using
`str(x).strip().upper()`. Trial accessions were matched to the global GRM
sample list; unmatched lines were skipped.

**Phenotype vector:** Standardized phenotype values were averaged per line
to form $y$, aligned to the GRM.

**GRM slicing:** The global GRM was sliced to the phenotyped lines:

$$
K = G_{\text{trial},\text{trial}}
$$

**GBLUP fitting:**

$$
\mu = \mathrm{mean}(y), \quad y_c = y - \mu, \quad
\alpha = (K + \lambda I)^{-1} y_c
$$

Outputs for each trial (`trained_models/<TRIAL>/`) included:

- `final_model.joblib`
- `GRM.npy` (trial-specific slice)
- `GRM_lines.txt`

## Prediction and cross-validation

All predictions for CV0 and CV00 used the global GRM and the ridge-GBLUP
model (`gblup_fit` + `gblup_predict`).

Common components:

- Global GRM and sample list
- Standardized unified phenotype file
- Normalized accession names
- Trial accession list from `data/raw/accession_lists/<TRIAL>.txt`

For both CV0 and CV00:

- A phenotype map was built from the (possibly masked) phenotype table.
- A vector $y$ and mask were constructed over global GRM samples.
- Training GRM and phenotype vector were:

$$
K_{\text{train}} = G[\text{mask},\text{mask}], \quad
y_{\text{train}} = y[\text{mask}]
$$

- A GBLUP model was fit with `gblup_fit(K_train, y_train)`.

Predictions:

- **Genotyped accessions:**

$$
\hat{y} = \mu + K_{\text{pred}} \alpha
$$

- **Non-genotyped accessions:** assigned the population mean $\mu$.

Outputs contained `line_name` and `prediction`.

### CV0 — new environments

- No phenotypes were masked.
- The model was trained on all phenotyped lines.
- Predictions were generated for all accessions in the focal trial.

### CV00 — new lines in new environments

- All phenotypes for lines present in the focal trial were removed.
- The model was trained on the remaining phenotyped lines.
- Predictions were generated for all trial accessions.

## Submission construction

For each trial and scheme (CV0, CV00), prediction files were assembled into
the Predictathon submission format. Predictions remained on the standardized
(z-score) scale.

## Workflow management

A Snakemake workflow (`run_pipeline.sh`) orchestrated phenotype cleaning,
genotype processing, GRM construction, model training, CV0/CV00 prediction,
and submission assembly.

The pipeline can be run with:

```bash
bash run_pipeline.sh
```

and rebuilt from scratch with:

```bash
bash run_pipeline.sh --clean
```

## Code availability

[redacted]
