# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a wheat genomic prediction project for the T3 Predictathon competition. The pipeline integrates genomic, phenotypic, and environmental data from the Triticeae Toolbox (T3) database to build predictive models for wheat yield across different trials and environments.

## Technology Stack

- **Language**: R (RStudio project)
- **Database Interface**: BrAPI to connect to wheat.triticeaetoolbox.org
- **Key R Packages**:
  - `BrAPI` (from github.com/TriticeaeToolbox/BrAPI.R) - T3 database connection
  - `EnvRtype` (from github.com/allogamous/EnvRtype) - environmental typing
  - `vcfR` - genomic VCF file processing
  - `rrBLUP`, `sommer`, `AGHmatrix` - genomic prediction models
  - `CovCombR` (from github.com/cran/CovCombR) - covariance matrix combination
  - Standard tidyverse packages (dplyr, readr, tidyr, stringr)
  - `janitor`, `here`, `kableExtra` - utilities

## Project Architecture

### Data Flow Pipeline

The project follows a strict sequential data processing pipeline with three main stages:

**Stage 1: Data Import** (`Data/` directory)
1. **First**: Run `T3_trial_import_and_filtering.Rmd` to filter wheat trials meeting criteria:
   - Wheat crop, harvested 2015-2025
   - >25 entries, >50% genotyped
   - Grain yield data present
2. **Then** (any order):
   - `Phenotypic_data_import.Rmd` - pulls phenotypic observations
   - `Genomic_data_import.Rmd` - downloads VCF files to `data/vcf_genotyped_data/`
   - `Environmental_data_import.Rmd` - pulls weather data via EnvRtype

**Stage 2: Similarity Matrix Construction** (`Analysis/` directory)
These can run in parallel but **must both complete before modeling**:
- `genomic_analysis.Rmd` - processes VCF files → trial-wise genomic similarity matrix (`trial_Gmat.rds`)
- `enviromic_analysis.Rmd` - processes weather data → environmental similarity matrix (`Trial_W_mat.csv`)

**Stage 3: Prediction Models** (`Prediction Models/` directory)
Models use similarity matrices from Stage 2:
- `GBLUP/` - Genomic Best Linear Unbiased Prediction
- `MegaLMM/` - Multi-Environment Genomic Linear Mixed Model
- `Reaction Norm/` - Reaction norm models

### Key Data Files

**Training Set Data**:
- `Data/genotyped_phenotyped_trials.csv` - filtered trial metadata
- `Data/genotyped_trial_phenotypes.csv` - phenotypic observations
- `Data/EnvRtype_output.csv` - raw weather data
- `Data/processed_weather_data.csv` - processed weather variables

**Test Set Data** (prediction targets):
- `Data/test_trial_metadata.csv` - test trial information
- `Data/testset_EnvRtype_output.csv` - test set weather data
- `Data/testset_processed_weather_data.csv` - processed test weather

**Intermediate Outputs**:
- `Analysis/trial_Gmat.rds` - trial × trial genomic similarity matrix
- `Analysis/Trial_W_mat.csv` - trial × trial environmental similarity matrix
- `Analysis/comb_grm_mat.rds` - combined genomic relationship matrix (gitignored, large)

## Running R Markdown Files

All analysis scripts are R Markdown (.Rmd) files that can be run in RStudio:
- Open in RStudio and click "Knit" to run entire document
- Or run individual chunks interactively with Cmd+Shift+Enter (Mac) / Ctrl+Shift+Enter (Windows)
- Outputs are generated in the same directory unless specified with `here::here()`

## Key Patterns and Conventions

### BrAPI Connection
Most scripts establish a connection to T3:
```r
conn <- createBrAPIConnection("wheat.triticeaetoolbox.org", is_breedbase = TRUE)
```

### File Path Management
All scripts use `here::here()` for relative paths from project root:
```r
readr::read_csv(here::here("data", "genotyped_phenotyped_trials.csv"))
```
This ensures scripts work regardless of working directory.

### Case Conventions
- CSV column names are often cleaned with `janitor::clean_names()` to snake_case
- T3 database fields use camelCase (e.g., `studyName`, `studyDbId`)
- Be aware of case sensitivity when joining datasets

### Genomic Data Processing
- VCF files are stored in `data/vcf_genotyped_data/` (gitignored due to size)
- VCF → dosage matrix conversion is done via custom functions in `genomic_analysis.Rmd`
- Some genotyping projects may be excluded (e.g., project 12856 is filtered out)

### Environmental Data
- Weather data uses EnvRtype package conventions
- Growth stages are defined by accumulated Growing Degree Days (GDD) with base temp 0°C
- Environmental variables are calculated across crop growth stages

### Modeling Approach
- Models use trial-level similarity matrices rather than individual-level matrices
- Training set selection for each test trial uses top N most similar trials based on W and G matrices
- The `select_training_data()` function pattern is used across modeling scripts

## Important Notes

- **Sequential execution is critical**: Data scripts must run in order, analysis scripts must complete before modeling
- **Large files are gitignored**: VCF files and combined GRM matrices are regenerated locally
- **BrAPI package installation**: Must be installed from GitHub (TriticeaeToolbox/BrAPI.R)
- **EnvRtype package**: Must be installed from GitHub (allogamous/EnvRtype)
- **CovCombR package**: Install from GitHub (cran/CovCombR)
- **Data dependencies**: All modeling scripts expect similarity matrices from Analysis/ to exist
- **Working directory**: Open `T3_predictathon.Rproj` in RStudio to set correct working directory
