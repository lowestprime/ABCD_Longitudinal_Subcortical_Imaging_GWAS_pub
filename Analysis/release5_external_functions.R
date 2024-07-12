#### EXTERNAL FUNCTIONS FOR ABCD.Release.5.1.C4.Merged.Data ####
### Cooper Beaman, 7/9/24

# recode timepoint data
recode.timepoint <- Vectorize(function(eventname) {
  parts <- unlist(strsplit(eventname, split = '_'))
  t <- parts[1]
  y <- parts[2]

  if (t == 'baseline') {
    return(0)
  } else if (y == 'year') {
    return(as.numeric(t))
  } else {
    return(as.numeric(t) / 12)
  }
})

# recode sex data
recode.sex <- Vectorize(function(kbi_sex_assigned_at_birth) {
  if (is.na(kbi_sex_assigned_at_birth)) {
    return('NA')
  } else if (kbi_sex_assigned_at_birth == 1) {
    return('M')
  } else if (kbi_sex_assigned_at_birth == 2) {
    return('F')
  } else {
    return('NA')
  }
})

# average hemispheres
average_hemispheres <- function(df) {
  regions <- df %>%
    names() %>%
    str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>%
    na.omit() %>%
    unique()
  avg_data <- sapply(setNames(
    lapply(regions, function(x) {
      c(paste0(x, ifelse(x == "smri_vol_scs_aa", "l", "lh")),
        paste0(x, ifelse(x == "smri_vol_scs_aa", "r", "rh")))
    }),
    regions
  ), function(region) {
    if (all(region %in% names(df)))
      rowMeans(df[, region, drop = FALSE], na.rm = TRUE)
    else
      rep(NA_real_, nrow(df))
  })
  result <- cbind(df, avg_data) %>%
    dplyr::select(-matches("^smri_vol_scs_.*(lh|rh|aal|aar)$"))
  
  return(result)
}

# Helper function to get latest timepoint information
get_latest_timepoint_info <- function(df, timepoint) {
  df %>%
    filter(timepoint == !!timepoint & sex != 'NA') %>%
    dplyr::select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age, ethnicity)
}

# baseline_y2_roc
baseline_y2_roc <- function(df, roc_volumes) {
  # Compute ROCs
  roc_df <- df %>%
    group_by(src_subject_id, rel_family_id) %>%
    reframe(across(
      .cols = all_of(roc_volumes),
      .fns = list(
        ROC0_2 = ~ ((.x[timepoint == 2] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 2)
      ),
      .names = "{.col}_{.fn}"
    ))
  
  # Get latest timepoint information for timepoint 2
  latest_info <- get_latest_timepoint_info(df, 2)
  
  # Join latest timepoint information
  roc_df <- roc_df %>%
    left_join(latest_info, by = c("src_subject_id", "rel_family_id"))
  
  # Filter out rows with 'NA' in the sex column
  roc_df <- roc_df %>%
    filter(sex != 'NA')
  
  return(roc_df)
}

# all_timepoints_roc
all_timepoints_roc <- function(df, roc_volumes) {
  # Compute ROCs
  roc_df <- df %>%
    group_by(src_subject_id, rel_family_id) %>%
    reframe(across(
      .cols = all_of(roc_volumes),
      .fns = list(
        ROC0_2 = ~ ((.x[timepoint == 2] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 2),
        ROC0_4 = ~ ((.x[timepoint == 4] - .x[timepoint == 0]) / .x[timepoint == 0]) * (100 / 4),
        ROC2_4 = ~ ((.x[timepoint == 4] - .x[timepoint == 2]) / .x[timepoint == 2]) * (100 / 2)
      ),
      .names = "{.col}_{.fn}"
    ))
  
  # Get latest timepoint information for timepoint 4
  latest_info <- get_latest_timepoint_info(df, 4)
  
  # Join latest timepoint information
  roc_df <- roc_df %>%
    left_join(latest_info, by = c("src_subject_id", "rel_family_id"))
  
  # Filter out rows with 'NA' in the sex column
  roc_df <- roc_df %>%
    filter(sex != 'NA')
  
  return(roc_df)
}

# Function to pivot ROC data to long format
pivot_roc_to_long_format <- function(df, is_baseline_y2 = FALSE) {
  # Identify ROC columns
  roc_columns <- df %>% dplyr::select(starts_with("smri_vol_")) %>% colnames()
  
  # Pivot only the ROC columns
  df %>%
    pivot_longer(
      cols = all_of(roc_columns),
      names_to = if (is_baseline_y2) "roi" else c("roi", "Time_Comparison"),
      names_pattern = if (is_baseline_y2) "(.*)" else "(.*)_(ROC.*)",
      values_to = "Value"
    ) %>%
    # Keep the metadata columns intact
    dplyr::select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age, ethnicity, roi, everything())
}

# Function to pivot original data to long format
pivot_original_to_long_format <- function(df, roc_volumes) {
  df %>%
    pivot_longer(
      cols = all_of(roc_volumes),
      names_to = "volume_type",
      values_to = "volume"
    ) %>%
    dplyr::select(all_of(c("src_subject_id", "rel_family_id", "sex", "interview_age", "eventname", "timepoint", "ethnicity", "volume_type", "volume")))
}

# Function to remove non-ROI columns from smri.R5.1.baseline.y2.ROC
roi_filter <- function(df) {
  # Vector of column names to be removed
  columns_to_remove <- c(
    "smri_vol_scs_3rdventricle_ROC0_2",
    "smri_vol_scs_4thventricle_ROC0_2",
    "smri_vol_scs_bstem_ROC0_2",
    "smri_vol_scs_csf_ROC0_2",
    "smri_vol_scs_inflatvent_ROC0_2",
    "smri_vol_scs_ltventricle_ROC0_2"
  )
  
  # Remove specified columns
  df_cleaned <- df %>% dplyr::select(-all_of(columns_to_remove))
  
  # Columns to be placed in specific positions
  specific_columns <- c("sex", "mri_info_deviceserialnumber", "interview_age", "ethnicity")
  
  # Columns that start with "smri_vol_scs" and not in columns_to_remove
  smri_columns <- colnames(df_cleaned) %>% 
    .[grepl("^smri_vol_scs", .)] %>%
    sort()
  
  # Remaining columns after removing specific and smri columns
  remaining_columns <- setdiff(colnames(df_cleaned), c(specific_columns, smri_columns))
  
  # Final column order
  final_column_order <- c(
    setdiff(colnames(df_cleaned), c(specific_columns, smri_columns)),
    specific_columns,
    smri_columns
  )
  
  # Reorder columns
  df_cleaned <- df_cleaned %>% dplyr::select(all_of(final_column_order))
  
  return(df_cleaned)
}

# Define a function to rename columns from V3-V22 to P1-P20 (used in preprocessing ancestry PCs)
rename_columns <- function(col_names) {
  old_names <- paste0("V", 3:22)
  new_names <- paste0("PC", 1:20)
  renamed_cols <- col_names
  renamed_cols[match(old_names, col_names)] <- new_names
  return(renamed_cols)
}

# Function to create dummy variables
create_dummies <- function(data, var_names) {
  # Ensure specified columns are treated as factors
  data <- data %>%
    mutate(across(all_of(var_names), as.factor))
  
  # Generate dummy variables using model.matrix
  dummies <- model.matrix(~ . - 1, data = data %>% dplyr::select(all_of(var_names))) %>%
    as.data.frame()
  
  # Ensure valid column names
  colnames(dummies) <- make.names(colnames(dummies), unique = TRUE)
  
  # Append dummy columns and remove original factor columns
  data <- data %>%
    dplyr::select(-all_of(var_names)) %>%
    bind_cols(dummies)
  
  return(data)
}

# Function to save split data with directory creation
save_split_data <- function(data, ethnicity, sex, pheno_dir, covar_dir, date) {
  # Define output directories
  pheno_out_dir <- file.path(pheno_dir, ethnicity, sex)
  covar_disc_out_dir <- file.path(covar_dir, "Discrete", ethnicity, sex)
  covar_quant_out_dir <- file.path(covar_dir, "Quantitative", ethnicity, sex)
  
  # Verify the paths exist
  if (!dir.exists(pheno_out_dir)) dir.create(pheno_out_dir, recursive = TRUE)
  if (!dir.exists(covar_disc_out_dir)) dir.create(covar_disc_out_dir, recursive = TRUE)
  if (!dir.exists(covar_quant_out_dir)) dir.create(covar_quant_out_dir, recursive = TRUE)
  
  # Subset data
  subset_data <- data %>%
    filter(ethnicity == !!ethnicity, sex == !!sex) %>%
    mutate(across(starts_with("smri_vol_"), ~ rankTransPheno(.x, 0.5)))
  
  # Specify the columns for dummy variable creation
  dummy_vars <- c("mri_info_deviceserialnumber", "batch")
  
  # Create dummy variables
  subset_data <- create_dummies(subset_data, dummy_vars)
  
  num_samples <- nrow(subset_data)
  id <- sample(1:1000000, 1) # generate a random ID
  
  # Save phenotype files
  phenotypes <- colnames(subset_data)[grepl("^smri_vol", colnames(subset_data))]
  for (phenotype in phenotypes) {
    if (phenotype != "smri_vol_scs_intracranialv_ROC0_2") {
      phenotype_file_name <- file.path(pheno_out_dir, sprintf("%s_pheno_%s_%s_%s_%s_%d_%s.txt", 
                                                              id, date, ethnicity, sex, num_samples, phenotype))
      write.table(subset_data %>% dplyr::select(FID, IID, all_of(phenotype)), 
                  file = phenotype_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
    
    # Save corresponding qcovar file
    if (phenotype == "smri_vol_scs_wholeb_ROC0_2") {
      qcovar_file_name_no_icv <- file.path(covar_quant_out_dir, sprintf("qcovar_noICV_%s_%s_%s_%s_%d.txt", 
                                                                        id, date, ethnicity, sex, num_samples))
      write.table(subset_data %>% dplyr::select(FID, IID, interview_age, starts_with("PC")), 
                  file = qcovar_file_name_no_icv, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    } else {
      qcovar_file_name <- file.path(covar_quant_out_dir, sprintf("qcovar_%s_%s_%s_%s_%d.txt", 
                                                                 id, date, ethnicity, sex, num_samples))
      write.table(subset_data %>% dplyr::select(FID, IID, interview_age, smri_vol_scs_intracranialv_ROC0_2, starts_with("PC")), 
                  file = qcovar_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  
  # Save discrete covariate file with dummy variables for specified columns
  covar_discrete <- subset_data %>%
    dplyr::select(FID, IID, sex, starts_with("mri_info_"), starts_with("batch")) %>%
    mutate(sex = ifelse(sex == "M", 1, 2))
  
  covar_file_name <- file.path(covar_disc_out_dir, sprintf("covar_%s_%s_%s_%s_%d.txt", 
                                                           id, date, ethnicity, sex, num_samples))
  write.table(covar_discrete, file = covar_file_name, col.names = TRUE, 
              row.names = FALSE, sep = "\t", quote = FALSE)
}
