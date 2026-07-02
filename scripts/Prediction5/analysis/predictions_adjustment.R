#### Predictions Adjuster #####

# Extract Training Sets
library(dplyr)
library(readr)
library(here)

##### 1. Load the target accessions (assuming 'pheno', 'grmList', and 'He' are already in your environment) #####
target_data <- read_csv(here("data","acc_for_prediction.csv"))

target_pheno <- target_data %>%
  dplyr::rename(GID = acc_name,
                LOC = location,
                YEAR = year) %>%  
  dplyr::mutate(YLD = NA) %>%
  dplyr::select(trial, ENV, GID, LOC, YEAR, YLD)

target_trials <- unique(target_pheno$trial)

# ==============================================================================
# EXTRACTOR 1: CV0 (Target trial removed, but historical line data kept)
# ==============================================================================
extract_training_cv0 <- function(tri) {
  
  target_sub <- target_pheno %>% dplyr::filter(trial == tri)
  if(!tri %in% names(grmList)){ return(NULL) }
  
  gid.geno <- rownames(grmList[[tri]])
  target_sub <- target_sub %>% dplyr::filter(GID %in% gid.geno)
  if(nrow(target_sub) == 0){ return(NULL) }
  
  # CV0 Logic
  pheno_sub <- pheno %>% dplyr::filter(GID %in% gid.geno)
  pheno_clean <- pheno_sub %>% dplyr::anti_join(target_sub, by = c("ENV", "GID"))
  
  robust_historical_envs <- pheno_clean %>%
    dplyr::group_by(ENV) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    dplyr::filter(n >= 30) %>%  
    dplyr::pull(ENV)
  
  rescue_env <- "Neoga_2022" 
  target_envs <- unique(target_sub$ENV)
  
  envs_to_keep <- unique(c(target_envs, robust_historical_envs, rescue_env))
  pheno_clean <- pheno_clean %>% dplyr::filter(ENV %in% envs_to_keep)
  
  if(tri == "YT_Urb_25"){
    pheno_clean <- pheno_clean %>% dplyr::mutate(numeric_year = as.numeric(as.character(YEAR)))
    recent_years <- sort(unique(pheno_clean$numeric_year), decreasing = TRUE)[1:2]
    pheno_clean <- pheno_clean %>%
      dplyr::filter(numeric_year %in% recent_years | ENV == rescue_env) %>%
      dplyr::select(-numeric_year)
  }
  
  # Final environment validation against the H-matrix
  missing_envs <- setdiff(target_envs, rownames(He))
  if(length(missing_envs) > 0){ return(NULL) }
  
  valid_envs <- intersect(unique(c(pheno_clean$ENV, target_envs)), rownames(He))
  pheno_clean <- pheno_clean %>% dplyr::filter(ENV %in% valid_envs)
  
  # Tag these records with the target trial they are meant to predict
  pheno_clean <- pheno_clean %>%
    dplyr::mutate(Target_Trial_Predicted = tri,
                  Scenario = "CV0") %>%
    dplyr::select(Scenario, Target_Trial_Predicted, ENV, GID, LOC, YEAR, YLD, dplyr::everything())
  
  return(pheno_clean)
}


# ==============================================================================
# EXTRACTOR 2: CV00 (Strict Blinding - Target lines completely ghosted)
# ==============================================================================
extract_training_cv00 <- function(tri) {
  
  target_sub <- target_pheno %>% dplyr::filter(trial == tri)
  if(!tri %in% names(grmList)){ return(NULL) }
  
  gid.geno <- rownames(grmList[[tri]])
  target_sub <- target_sub %>% dplyr::filter(GID %in% gid.geno)
  if(nrow(target_sub) == 0){ return(NULL) }
  
  target_gids <- unique(target_sub$GID)
  target_envs <- unique(target_sub$ENV)
  
  # CV00 STRICT Logic
  pheno_clean <- pheno %>% 
    dplyr::filter(GID %in% gid.geno) %>%
    dplyr::filter(!ENV %in% target_envs) %>%
    dplyr::filter(!GID %in% target_gids)
  
  robust_historical_envs <- pheno_clean %>%
    dplyr::group_by(ENV) %>%
    dplyr::summarize(n = n(), .groups = "drop") %>%
    dplyr::filter(n >= 30) %>%  
    dplyr::pull(ENV)
  
  rescue_env <- "Neoga_2022" 
  
  envs_to_keep <- unique(c(target_envs, robust_historical_envs, rescue_env))
  pheno_clean <- pheno_clean %>% dplyr::filter(ENV %in% envs_to_keep)
  
  if(tri == "YT_Urb_25"){
    pheno_clean <- pheno_clean %>% dplyr::mutate(numeric_year = as.numeric(as.character(YEAR)))
    recent_years <- sort(unique(pheno_clean$numeric_year), decreasing = TRUE)[1:2]
    pheno_clean <- pheno_clean %>%
      dplyr::filter(numeric_year %in% recent_years | ENV == rescue_env) %>%
      dplyr::select(-numeric_year)
  }
  
  # Final environment validation against the H-matrix
  missing_envs <- setdiff(target_envs, rownames(He))
  if(length(missing_envs) > 0){ return(NULL) }
  
  valid_envs <- intersect(unique(c(pheno_clean$ENV, target_envs)), rownames(He))
  pheno_clean <- pheno_clean %>% dplyr::filter(ENV %in% valid_envs)
  
  # Tag these records with the target trial they are meant to predict
  pheno_clean <- pheno_clean %>%
    dplyr::mutate(Target_Trial_Predicted = tri,
                  Scenario = "CV00") %>%
    dplyr::select(Scenario, Target_Trial_Predicted, ENV, GID, LOC, YEAR, YLD, dplyr::everything())
  
  return(pheno_clean)
}

# Run the sequential extraction loops
write("Extracting CV0 training sets...", stderr())
cv0_training_list <- lapply(target_trials, extract_training_cv0)

write("Extracting CV00 training sets...", stderr())
cv00_training_list <- lapply(target_trials, extract_training_cv00)

# Bind the results into master dataframes
master_cv0_training <- dplyr::bind_rows(cv0_training_list)
master_cv00_training <- dplyr::bind_rows(cv00_training_list)

# Export to CSV
write.csv(master_cv0_training, "CV0_Training_Sets_Extracted.csv", row.names = FALSE)
write.csv(master_cv00_training, "CV00_Training_Sets_Extracted.csv", row.names = FALSE)

write("Extraction complete! Check your working directory for the CSVs.", stdout())

#####
CV0_train <- read.csv(here("analysis","CV0_Training_Sets_Extracted.csv"))
CV00_train <- read.csv(here("analysis","CV00_Training_Sets_Extracted.csv"))

CV0 <- read.csv(here("analysis","CV0_Final_Yield_Predictions.csv"))
CV00 <- read.csv(here("analysis","CV00_Final_Yield_Predictions.csv"))



CV0_yld <- data.frame()
CV00_yld <- data.frame()

for(trials in target_trials){
  CV0_filt <- CV0 %>% filter(trial == trials)
  temp_gid <- CV0_filt[,3]
  temp_yld <- CV0_filt[,7]
  temp_yld <- scale(temp_yld)
  CV0_filt <- as.data.frame(cbind(temp_gid, temp_yld))
  colnames(CV0_filt)[1] = "GID"
  colnames(CV0_filt)[2] = "YLD"
  CV0_train_filt <- CV0_train %>% filter(Target_Trial_Predicted == trials)
  mean_yld <- mean(CV0_train_filt$YLD)
  CV0_filt <- CV0_filt %>%
    mutate(
      YLD = as.numeric(YLD) + mean_yld
    )
  CV0_yld <- rbind(CV0_yld, CV0_filt)
  
  CV00_filt <- CV00 %>% filter(trial == trials)
  temp_gid <- CV00_filt[,3]
  temp_yld <- CV00_filt[,7]
  temp_yld <- scale(temp_yld)
  CV00_filt <- as.data.frame(cbind(temp_gid, temp_yld))
  colnames(CV00_filt)[1] = "GID"
  colnames(CV00_filt)[2] = "YLD"
  CV00_train_filt <- CV00_train %>% filter(Target_Trial_Predicted == trials)
  mean_yld <- mean(CV00_train_filt$YLD)
  CV00_filt <- CV00_filt %>%
    mutate(
      YLD = as.numeric(YLD) + mean_yld
    )
  CV00_yld <- rbind(CV00_yld, CV00_filt)
}

CV0_lines <- CV0[,c(1,3)]
CV00_lines <- CV00[,c(1,3)]

CV0_yld <- left_join(CV0_yld, CV0_lines)
CV00_yld <- left_join(CV00_yld, CV00_lines)

write.csv(CV0_yld, here("analysis", "CV0_Predictions.csv"))
write.csv(CV00_yld, here("analysis", "CV00_Predictions.csv"))
