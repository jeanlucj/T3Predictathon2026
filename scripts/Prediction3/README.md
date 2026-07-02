# T3_predictathon

## [Data](Data)
### *R Script*
1) **T3_trial_import_and_filtering.Rmd:** Code to import all trials from T3 and filter trials to include only ones meeting certain criteria
2) **Phenotypic_data_import.Rmd:** Code to pull phenotypic observations from all selected trials
3) **Genomic_data_import.Rmd:** Code to pull variant call files to gather genotypic data for accessions grown in selected trials
4) **Environmental_data_import.Rmd:** Code to import environmental data from EnvRtype for selected trials

### *Data Output Files*
1) **drone_traits_lookup.csv** Data frame of drone trait names and identifiers
2) **EnvRtype_output.csv** Data frame of weather data for training set candidate trials
3) **genotyped_phenotyped_trials.csv** Data frame of metadata for candidate training set trials from the trial filtering script
4) **genotyped_trial_phenotypes.csv** Data frame of phenotypes for training set data from the phenotypic import script
5) **location_soil_data.csv** Data frame of soil data from SoilGrids for training set trial locations
6) **processed_weather_data.csv** Data frame of processed weather data for training set trials
7) **selected_geno_projects.csv** Data frame of T3 genotyping_project_ids for projects containing accessions in training or test set trials
8) **test_location_soil_data.csv** Data frame of soil data from locations of test set trials
9) **test_trial_metadata.csv** Data frame of metadata for test set trials
10) **testset_EnvRtype_output.csv** Data frame of weather data for test set trials
11) **testset_processed_weather_data.csv** Data frame of processed weather data for test set trials
12) **trial_geno_search_parallel_safe_output.csv** Data frame matching accessions from all training set trials to a genotyping protocol or project

## [Analysis](Analysis)
### *R Script*
1) **enviromic_analysis.Rmd** Code to process and summarize weather data, cluster environments, and create the Enviromental Similarity matrix between trials
2) **genomic_analysis.Rmd** Code to process marker vcf files, construct a combined GRM of all trial accessions, cluster trials genetically, and create the Genomic Similarity matrix between trials

### *Data Output Files*
1) **Trial_W_mat.csv** Trial-wise environmental similarity matrix of test and training set trials
2) **Trial_Gmat.rds** Trial-wise genomic similarity matrix of test and training set trials
3) **trial_env_distances** Matrix of distances of trials in the environmental PCA space


## [figures](figures)
1) **W_rel_heatmap.png**
2) **env_pca.png**
3) **env_var_corr_plot.png**

## [Prediction Models](Prediction Models)
1) **GBLUP_predictions.Rmd** Code for training set data selection, base GBLUP prediction model, and the multi-kernel BLUP prediction model that was ultimately used to generate submitted predictions under CV0 and CV00 scenarios. Leave-One-Out Cross Validation also performed.

## [Predictions](Predictions)
• This folder contains the submitted predictions with the appropriate file structure for the T3 predictathon submission form.
