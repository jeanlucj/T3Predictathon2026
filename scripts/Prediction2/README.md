# Predictathon 2025 Genomic Prediction Pipeline
### Project: GBLUP Pipeline for Predictathon 2025

## Overview
This repository contains a complete, end‑to‑end pipeline for generating CV0 and CV00 predictions for all nine Predictathon trials. The workflow performs phenotype cleaning, genotype preprocessing, GRM construction, mixed‑model training, and Predictathon‑compliant output generation.

All genotype and phenotype inputs originate from T3/Wheat and must be downloaded locally before running the pipeline.

## Key Features
- Unified phenotype cleaning and accession harmonization

- Genotype merging, platform harmonization, and GRM construction

- Mixed‑model GBLUP training using the VanRaden relationship matrix

- Challenge‑compliant CV0 and CV00 masking

- Per‑trial prediction and expected‑accuracy estimation

- Automatic generation of Predictathon submission folders


## Repository Structure
```
data/
  raw/                     # Accession lists and global yield stats
  processed/               # Unified phenotypes and mapped ids
  predictathon/            # Per-trial genotype and phenotype folders

src/
  genotypes/               # Vcf preprocessing and grm construction
  model/                   # Training, cv0, cv00, expected accuracy
  utils/                   # Shared helpers

trained_models/            # Saved grms and trained models
results/                   # Cv0, cv00, and expected accuracy outputs
submission/                # Predictathon-compliant final submission

config.yaml                # Trial, trait, and model settings
run_pipeline.sh            # End-to-end pipeline runner
```

## Challenge Trials

AWY1_DVPWA_2024

TCAP_2025_MANKS

25_Big6_SVREC_SVREC

OHRWW_2025_SPO

CornellMaster_2025_McGowan

24Crk_AY2-3

2025_AYT_Aurora

YT_Urb_25

STP1_2025_MCG

## Required Inputs (Must Exist Locally)
These files must be present for the pipeline to run.
Large genotype files are not stored in GitHub and must be downloaded manually.

1. Phenotype Inputs
- data/processed/unified_training_pheno_mapped.csv: Unified phenotype table with mapped accessions and a trial column
- data/predictathon//training_pheno_merged.csv: Per‑trial merged phenotype file (optional but recommended)
3. Genotype Inputs (LFS‑sized, not stored in GitHub)
- data/predictathon//genotypes/*.vcf.gz: Raw genotype VCFs downloaded from T3/Wheat
These files are typically hundreds of MB to several GB and must remain local.

4. Metadata / Mapping Files
- Metadata files for each predictathon trial downloaded from T3/Wheat
- Any accession‑mapping tables: Required for phenotype/genotype alignment

## Running the Pipeline
```
bash run_pipeline.sh
```
This will:

- Preprocess genotypes (if GRM not already built)
- Train the GBLUP model
- Generate CV0 predictions
- Generate CV00 predictions
- Compute expected accuracy
- Write outputs to results/ and submission/

## Force a Clean Rebuild
```
bash run_pipeline.sh --clean
```
Removes cached GRMs, models, and predictions before rerunning.

## Outputs
All final Predictathon‑formatted outputs are written to:

```
submission/
  <TRIAL>/
    CV0/
    CV00/
```
Each folder contains:
- Predictions
- Accessions used
- Trials used for training

