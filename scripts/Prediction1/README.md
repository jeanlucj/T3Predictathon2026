# production--prediction
Code for T3/Wheat Predicathon


1. Overview
We implemented a Genomic Best Linear Unbiased Prediction (GBLUP) approach to predict grain yield for nine focal trials under two cross-validation scenarios (CV0 and CV00). All analyses were conducted in R (v4.3.3).

2. Training Data Selection
Training trials for each focal trial were identified using the Training_Trial_Info.rds file provided by the Predictathon organizers (retrieved from the official GitHub repository: jeanlucj/T3_predictathon_find_training_trials). This file was automatically downloaded at the start of each pipeline run to ensure the most up-to-date trial selection was used.
Training trials were categorized into two tiers:
•	Primary trials: historical trials that evaluated at least n accessions in common with the focal trial (threshold varied per focal trial: 3-6 accessions).
•	Secondary trials: trials that evaluated accessions in common with the primary pool, applying a higher overlap threshold (4-20 accessions depending on the focal trial).
Phenotypic data for all primary and secondary training trials were downloaded from the T3/Wheat database via BrAPI, using the BrAPI and T3BrapiHelpers R packages. Only grain yield observations (T3 variable ID: 84527) were retained. In total, phenotypic records from 1,415 historical training studies encompassing 17,936 unique accessions were used.

3. Phenotypic Data Processing and BLUE Estimation
Raw yield records were processed separately for each training study. The following steps were applied:
1.	Outlier removal: observations with a standardized residual |z| > 3 (relative to the trial mean and standard deviation) were excluded.
2.	BLUE estimation: Best Linear Unbiased Estimates (BLUEs) for each accession were obtained by fitting a fixed-effects linear model: yield ~ germplasm + rep + block. Terms were included only when the corresponding design factor had more than one level. When the linear model failed to converge (e.g., due to severely unbalanced designs), the per-accession mean yield was used as a fallback.
3.	Cross-study aggregation: For each focal trial, BLUEs from all training studies were pooled. When an accession appeared in multiple training studies, its BLUEs were averaged across studies before being used as the training phenotype.

4. Genotypic Data Processing
Genotype data were obtained as VCF files downloaded from the T3/Wheat archive. The following quality control steps were applied to each VCF:
4.	Duplicate SNP removal: SNPs with identical CHROM:POS keys were deduplicated, retaining one record per genomic position.
5.	Genotype encoding: Genotypes were encoded additively as -1 (homozygous reference), 0 (heterozygous), and 1 (homozygous alternative) using the vcfR package.
6.	Sample filtering: Samples with missing genotype rate > 20% were removed.
7.	SNP filtering: SNPs were retained if missing rate <= 20% and minor allele frequency (MAF) >= 0.05.
8.	Missing value imputation: Remaining missing genotypes were imputed with the rounded per-SNP mean (nearest integer), preserving the additive encoding.

5. Genomic Relationship Matrix (GRM)
For each focal trial, the SNP matrix of the focal trial was used as the starting point. Additional SNP matrices from other VCF files were iteratively merged by retaining the intersection of SNP positions, provided that at least 500 SNP positions were shared. Duplicate accessions were removed after merging.
The GRM was constructed using the VanRaden (2008) method implemented in the AGHmatrix R package (Gmatrix(), method = "VanRaden", MAF = 0.05). A small ridge term (1 x 10-4 x I) was added to the diagonal to ensure positive definiteness.
A key limitation was encountered for several focal trials: the genotyping platform used for the focal trial VCF had no SNP positions in common with VCF files containing historical training accessions. In these cases, training accessions could not be incorporated into the GRM, making genomic prediction infeasible (see Section 6).

6. GBLUP Model and Prediction
Genomic predictions were generated using the rrBLUP R package (mixed.solve()). For each focal trial, the following procedure was applied:
9.	Training phenotypes (BLUEs) were standardized (mean = 0, SD = 1) prior to model fitting.
10.	The GBLUP model was fitted on training accessions present in the GRM: y = mu + Zu + e, where u ~ N(0, sigma2_u * G).
11.	Breeding values for test accessions were obtained via the conditional expectation: u_test = G21 * G11^-1 * u_train, where G21 is the submatrix of the GRM between test and training accessions, and G11 is the GRM among training accessions. Matrix inversion was performed using solve(); the Moore-Penrose pseudoinverse (MASS::ginv()) was used as fallback when the matrix was singular.
12.	Predictions were back-transformed to the original yield scale.
Fallback procedure: When fewer than 10 training accessions were present in the GRM — due to the absence of shared SNP positions between the focal trial VCF and historical training VCF files — genomic prediction was not feasible. In these cases, the global mean of training BLUEs was assigned as the prediction for all test accessions.

Trial	CV0	CV00
24Crk_AY2-3	Mean imputation	Mean imputation
CornellMaster_2025_McGowan	Mean imputation	Mean imputation
STP1_2025_MCG	Mean imputation	Mean imputation
AWY1_DVPWA_2024	GBLUP	Mean imputation
OHRWW_2025_SPO	GBLUP	Mean imputation


7. Cross-Validation Scenarios
CV0: The training set included all accessions from primary and secondary training trials identified for each focal trial. No accessions from the focal trial itself were included in training.
CV00: Same as CV0, with the additional exclusion of any accessions present in the focal trial's genotype file (i.e., accessions that were genotyped as part of the focal trial were removed from training, regardless of whether they appeared in other historical trials.

8. Software and Packages

Package	Version	Purpose
R	4.3.3	Computing environment
vcfR	1.16.0	VCF parsing and genotype extraction
AGHmatrix	—	GRM construction (VanRaden method)
rrBLUP	—	GBLUP model fitting
tidyverse	2.0.0	Data processing
data.table	—	Fast CSV reading
BrAPI	—	T3/Wheat database access
T3BrapiHelpers	—	Trial and phenotype retrieval


9. References
VanRaden, P.M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414-4423.
Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. The Plant Genome, 4(3), 250-255.
Amadeu, R.R. et al. (2016). AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species. Crop Science, 56(3), 1114-1120.
