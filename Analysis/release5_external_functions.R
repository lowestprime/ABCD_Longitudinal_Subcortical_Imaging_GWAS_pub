#### EXTERNAL FUNCTIONS FOR ABCD.Release.5.1.C4.Merged.Data ####
### Cooper Beaman, 5/18/24

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
    select(-matches("^smri_vol_scs_.*(lh|rh|aal|aar)$"))
  
  return(result)
}

# Helper function to get latest timepoint information
get_latest_timepoint_info <- function(df, timepoint) {
  df %>%
    filter(timepoint == !!timepoint & sex != 'NA') %>%
    select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age, ethnicity)
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
  roc_columns <- df %>% select(starts_with("smri_vol_")) %>% colnames()
  
  # Pivot only the ROC columns
  df %>%
    pivot_longer(
      cols = all_of(roc_columns),
      names_to = if (is_baseline_y2) "roi" else c("roi", "Time_Comparison"),
      names_pattern = if (is_baseline_y2) "(.*)" else "(.*)_(ROC.*)",
      values_to = "Value"
    ) %>%
    # Keep the metadata columns intact
    select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age, ethnicity, roi, everything())
}

# Function to pivot original data to long format
pivot_original_to_long_format <- function(df, roc_volumes) {
  df %>%
    pivot_longer(
      cols = all_of(roc_volumes),
      names_to = "volume_type",
      values_to = "volume"
    ) %>%
    select(all_of(c("src_subject_id", "rel_family_id", "sex", "interview_age", "eventname", "timepoint", "ethnicity", "volume_type", "volume")))
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
  df_cleaned <- df %>% select(-all_of(columns_to_remove))
  
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
  df_cleaned <- df_cleaned %>% select(all_of(final_column_order))
  
  return(df_cleaned)
}

# split and save phenotype files
create_phenotype_files <- function(data, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Get the column names that start with "smri"
  smri_cols <- grep("^smri", names(data), value = TRUE)
  
  # Loop through each "smri" column
  for (col in smri_cols) {
    # Create a new data frame with the other columns and the current "smri" column
    phenotype_data <- data %>% select(IID, FID, sex, PC1:20, all_of(col))
    
    # Define the file name based on the phenotype column name
    file_name <- file.path(output_dir, paste0(col, ".txt"))
    
    # Save the data frame to a text file
    write.table(phenotype_data, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}


