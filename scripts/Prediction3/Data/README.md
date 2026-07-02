This folder contains all the code necessary to pull data needed for subsequent analysis. **Data must be read-in sequentially**

**Code for Filtering T3 Trials** 
1) [Code to import and filter T3 trials](T3_trial_import_and_filtering.Rmd)

The compiled list of T3 trials that met the following filtering criteria can be viewed here: [genotyped_phenotyped_trials.csv](genotyped_phenotyped_trials.csv)
* Wheat 
* Harvested between 2015-2025 
* Greater than 25 entries 
* Greater than 50% of entries genotyped 
* Grain yield data present 

**Code for Importing Data for Filtered Trials** *(can run in any order)*
1) [Code to import environmental data from EnvRtype](Environmental_data_import.Rmd) 
2) [Code to import genomic data from T3](Genomic_data_import.Rmd)
3) [Code to import phenotypic data from T3](Phenotypic_data_import.Rmd)

**Output Data Files** *(need to read in for analyses)*
1) [drone_traits_lookup.csv](drone_traits_lookup.csv) : List of Crop Ontology trait identifiers and names for drone-derived traits, gathered from T3
2) [EnvRtype_output.csv](EnvRtype_output.csv) : Weather data for candidate trainning set trials pulled from EnvRtype
3) [genotyped_phenotyped_trials.csv](genotyped_phenotyped_trials.csv) : Compiled dataframe of trial information for T3 trials with genomic data and grain yield phenotypes
4) [genotyped_trial_phenotypes.csv](genotyped_drone_trial_phenotypes.csv) : Compiled dataframe of field trial phenotypic observations from candidate training set trials
5) [location_soil_data.csv](location_soil_data.csv) : Soil data for locations of candidate training set trials pulled from EnvRtype
6) [processed_weather_data.csv](processed_weather_data.csv) : Output from the EnvRtype::processWTH() function for candidate training set trial weather data
7) [selected_geno_projects.csv](selected_geno_projects.csv) : List of T3 ids for genotyping projects selected to pull genomic data for accessions grown in selected trials
8) [test_location_soil_data.csv](test_location_soil_data.csv) : Soil data for locations of trials from the prediction set
9) [test_trial_metadata.csv](test_trial_metadata.csv) : Compiled dataframe of trial information for prediction set trials
10) [testset_EnvRtype_output.csv](testset_EnvRtype_output.csv) : Weather data for prediction set trials pulled from EnvRtype
11) [testset_processed_weather_data.csv](testset_processed_weather_data.csv) : Output from the EnvRtype::processWTH() function for prediction set trial weather data
12) [trial_geno_search_parallel_safe_output.csv](trial_geno_search_parallel_safe_output.csv) : Output from trial_geno_search_parallel_safe() function in the T3_trial_import_and_filtering.Rmd script matching all accessions to a genotyping protocol or project in T3


