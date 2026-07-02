
# Methods description T3 Predictathon

The architecture of the repository is following the workflowr recommendations, all the scripts used are contained in the “analysis” folder.

# Data analysis pipeline description

## Data formatting

- “phenotypic_data.Rmd” was used to retrieve phenotypic data from T3 using BrAPI requests. All observations from the ‘focal accessions’ were requested, as well as observations from the trials tested under the same breeding program. When the number of trials was superior to 20, filtering was done to keep the trials from the same location or from the same year as the trial to be predicted. Alternatively, phenotypic data was aggregated from the folder provided by Jean-Luc Jannink. Later, a filtering was done to keep the accessions with genotypic data.

- “pullGenoTrial.Rmd” was used to convert the downloaded vcf files made available on the T3 Predictathon website into an imputed numeric matrix. The conversion was done using a custom function “format_curate_vcf”, sourced from the “code/useful_functions.R” file.

- “pullGenoAssociatedTrials.Rmd” was used to retrieve more genotypic data from T3. Briefly, given the list of trials from the “phenotypic_data.Rmd” script, we found the best genotyping project and download the vcf(s) with BrAPI requests.

- “combine_GRM.Rmd” was used to convert the other genotypic data coming from “pullGenoAssociatedTrials.Rmd” into a numeric marker matrix, using the same function “format_curate_vcf”. Then, the genomic relationship matrices (GRM) were computed from all available marker matrices, using the function “covariance_combiner” from the package “T3BrapiHelpers”.

- “environmental_covariates.Rmd” was used to get weather variables (T2M, PRECTOTCORR, ALLSKY_SFC_SW_DWN) for each trial with phenotypic data from NASA POWER.

## Predictions and output

- “genom_pred.Rmd” was used to run the genomic prediction from the formatted outputs. For both CV0 and CV00 scenarios, GBLUP and RKHS methods were implemented. For each breeding program, a subset of genomic and phenotypic data was made of the common genotypes. Environments with more the 30 genotypes were selected to run the cross-validation. For each environment, data from the validation set was masked according to the prediction scenario. The predictive ability was calculated within each environment and across-environments for each breeding program.

- “analyze_genom_pred.Rmd” was meant to compare the genomic prediction methods and scenarios (not used).
Prediction model training


The genomic relationship matrix construction was done in the script “combine_GRM.Rmd’ and described above.
CV0 and CV00 cross-validation were performed across environments for both GBLUP and RKHS. 
Sommer GBLUP with GxE and environment effects were tested but not implemented. The GBLUP model was selected for CV0, with a fixed effect of the location and year (factor). The RKHS model was selected for CV00 with location and year effects fitted with BRR and genomic effect fitted with RKHS. For both scenarios, we used a similar training set with all phenotype data available for the lines genotyped under this trial with a one-step model.


