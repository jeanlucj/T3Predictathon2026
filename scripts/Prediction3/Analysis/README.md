This folder contains all the code necessary to run analyses for T3 predictathon methods **Code must be run sequentially**

**Code for estimating trial genomic and enviromic similarities** *(can run in any order but need to be run before modeling scripts*
1) [genomic_analysis.Rmd](genomic_analysis.Rmd) : Code to process marker vcf files, construct a combined GRM of all trial accessions, cluster trials genetically, and create the Genomic Similarity matrix between trials
2) [enviromic_analysis.Rmd](enviromic_analysis.Rmd) : Code to process and summarize weather data, cluster environments, and create the Enviromental Similarity matrix between trials


**Output Data Files** *(need to read in for analyses)*
1) [Trial_W_mat.csv](Trial_W_mat.csv) : Trial-wise environmental similarity matrix of test and training set trials
2) [trial_Gmat.rds](trial_Gmat.rds) : Trial-wise genomic similarity matrix of test and training set trials
3) [trial_env_distances.csv](trial_env_distances.csv) : Matrix of distances of trials from the centroid in the environmental PCA space
4) [trial_gen_distances.csv](trial_gen_distances.csv) : Matrix of distances of trials from the centroid in genomic PCA space
