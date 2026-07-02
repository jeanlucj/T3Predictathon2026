# Function to save an archived VCF to data for later use
# Assumes here::i_am has been set
# If you want the VCF saved to a subdirectory within data, use the prefix
save_archived_vcf <- function(project_id, file_name, prefix){
  suffix <- stringr::str_extract(file_name, "file.*$")
  file_name_vcf <- paste0(prefix, "_genotypes_",
                            project_id, "_",
                            suffix, ".vcf")
  file_name_vcf <- here::here("data", file_name_vcf)
  vcf_status <- wheat_conn$vcf_archived(output=file_name_vcf,
                          genotyping_project_id = project_id,
                          file_name = file_name)
  if (vcf_status$status$category != "Success"){
    cat("Unsuccessful VCF download. Project ID", project_id, "File Name", file_name, "\n")
  }
  return(file_name_vcf)
}

# Function to save the yield values from identified trials to data in a csv file
# If you want the csv saved to a subdirectory within data, use the prefix
save_yield_phenotypes <- function(study_id, prefix){

  file_name_yield <- paste0(prefix, "_phenotypes_", study_id, ".csv")
  file_name_yield <- here::here("data", file_name_yield)

  # Function to get the coordinates needed to do a spatial analysis
  extract_obs_unit_coord <- function(search_res){
    plot_column <- search_res$observationUnitPosition$positionCoordinateX %||%
      NA_integer_
    plot_row <- search_res$observationUnitPosition$positionCoordinateY %||%
      NA_integer_
    # Make a df out of the observationLevelRelationships
    obs_tags <- search_res$observationUnitPosition$observationLevelRelationships |> bind_rows()
    obs_tags <- dplyr::select(obs_tags, -levelOrder) |>
      tidyr::pivot_wider(names_from = levelName, values_from = levelCode) |>
      dplyr::mutate(row=plot_row, column=plot_column,
                    observationUnitDbId=search_res$observationUnitDbId)
    return(obs_tags)
  }

  # observationVariableDbIds=c(84527) is the yield variable
  yield_phenotypes <- wheat_conn$search("observations",
    body=list(studyDbIds=c(study_id),
              observationVariableDbIds=c(84527)))$combined_data |>
    purrr::map(function(sl){ sl$season <- sl$season[[1]]; sl}) |>
    dplyr::bind_rows()

  if (nrow(yield_phenotypes) > 0){
    yield_phenotypes <- yield_phenotypes |>
      dplyr::select(studyDbId, season,
                    germplasmDbId, germplasmName,
                    observationUnitDbId, observationUnitName,
                    observationVariableName, value)

    resp <- wheat_conn$search(
      "observationunits",
      body=list(observationUnitDbIds=yield_phenotypes$observationUnitDbId))

    plot_coords <- lapply(resp$combined_data, extract_obs_unit_coord) |>
      dplyr::bind_rows()

    yield_phenotypes <- merge(yield_phenotypes, plot_coords) |>
      dplyr::select(studyDbId, season,
                    germplasmDbId, germplasmName,
                    observationUnitDbId, observationUnitName,
                    rep, block, plot, row, column,
                    observationVariableName, value) |>
      dplyr::arrange(plot) |>
      janitor::clean_names()

    readr::write_csv(yield_phenotypes, file_name_yield)
    return(file_name_yield)
  } else{
    return(paste0("Trial ", study_id, " did not have any yield observations"))
  }
}

# Download all of the trials identified the evaluate related information
# to a Predictathon trial specific folder
download_training_trials <- function(training_trial_info){

  # Make sure the subdirectory is there
  subdir <- here::here("data", paste0("study", training_trial_info$study_id))
  if (!dir.exists(subdir)) dir.create(subdir)

  # This prefix targets the download to the specific folder
  # Distinguish between trials that phenotyped actual Predictathon accessions
  # (primary)
  prefix <- paste0("study", training_trial_info$study_id, "/", "primary")
  primary_file_names <- purrr::map(training_trial_info$primary_study_db_ids,
                                   save_yield_phenotypes, prefix=prefix,
                                   .progress=TRUE) |> unlist()

  # Versus trials that phenotyped accessions in other trials that had
  # Predictathon accessions
  prefix <- paste0("study", training_trial_info$study_id, "/", "secondary")
  secondary_file_names <- purrr::map(training_trial_info$secondary_study_db_ids,
                                     save_yield_phenotypes, prefix=prefix,
                                     .progress=TRUE) |> unlist()

  return(list(primary_file_names=primary_file_names,
              secondary_file_names=secondary_file_names))
}
