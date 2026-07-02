# Functions to work with VCFs

# Process a list of training population phenotyping trials training_trial_info
# coming out of `choose_training_trials`
find_useful_vcfs <- function(training_trial_info, min_germ_in_vcf=5){

  # All accessions that are relevant to prediction and training
  accessions <- c(training_trial_info$focal_germplasm,
                  training_trial_info$primary_germplasm,
                  training_trial_info$secondary_germplasm)
  accessions <- accessions[!is.na(accessions)]

  # accessions_by_genotyping_project is a named vector with how many accessions
  # were typed in a particular VCF
  useful_vcfs <- wheat_conn$filter_geno_projects(accessions)
  useful_vcfs_counts <- useful_vcfs$content$results$counts$accessions_by_genotyping_project

  return(useful_vcfs_counts[useful_vcfs_counts >= min_germ_in_vcf])
}

# Use find_useful_vcfs and then download them to a Predictathon-trial specific
# folder
download_useful_vcfs <- function(training_trial_info, min_germ_in_vcf=5,
                                 actually_save_vcfs=FALSE){

  study_id <- training_trial_info$study_id
  # Make sure the subdirectory is there
  subdir <- here::here("data", paste0(study_id))
  if (!dir.exists(subdir)) dir.create(subdir)

  # This prefix targets the download to the specific folder
  prefix <- paste0("study", study_id, "/training")

  # "Useful" means a VCF that genotypes at least 5 accessions we care about
  useful_vcfs <- find_useful_vcfs(training_trial_info,
                                  min_germ_in_vcf=min_germ_in_vcf)
  vcf_project_ids <- names(useful_vcfs)

  # Get information on the archived VCFs for those projects
  ifd <- lapply(
    vcf_project_ids,
    function(pid) wheat_conn$vcf_archived_list(genotyping_project_id=pid))

  # Sometimes vcf_archived_list comes up empty in which case remove then bind
  info_for_download <- ifd[which(sapply(ifd, function(df) nrow(df))>0)] |>
    dplyr::bind_rows() |>
    dplyr::select(project_id, file_name)

  file_info_for_download <- paste0(here::here("data", prefix),
                                   "_vcf_names", study_id, ".csv")
  readr::write_csv(info_for_download, file=file_info_for_download)

  info_for_download <- info_for_download |> dplyr::mutate(prefix=prefix)
  if (actually_save_vcfs){
    genotype_file_names <- purrr::pmap(info_for_download,
                                       purrr::safely(save_archived_vcf),
                                       .progress=TRUE)
    results <- lapply(genotype_file_names, function(re) re$result) |> unlist()
    error <- lapply(genotype_file_names, function(er) er$error) |> unlist()
    if (!is.null(error)){
      print("Errors occurred in downloading VCFs")
      print(error)
    }
  } else{
    make_vcf_name <- function(project_id, file_name, prefix){
      suffix <- stringr::str_extract(file_name, "file.*$")
      file_name_vcf <- paste0(prefix, "_genotypes_",
                              project_id, "_",
                              suffix, ".vcf")
      return(here::here("data", file_name_vcf))
    }
    genotype_file_names <- purrr::pmap(info_for_download, make_vcf_name) |>
      unlist()
  }

  return(list(genotype_file_names=genotype_file_names))
}
