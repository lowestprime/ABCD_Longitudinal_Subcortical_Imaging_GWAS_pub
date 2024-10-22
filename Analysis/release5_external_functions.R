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

# Function to create dummy variables and return the modified data
create_dummies <- function(data, var_names) {
  # Create dummy variables using fastDummies package
  data_with_dummies <- dummy_cols(data, select_columns = var_names, remove_first_dummy = FALSE, remove_selected_columns = TRUE)
  
  # Dummy variable column names
  dummy_colnames <- setdiff(colnames(data_with_dummies), colnames(data))
  
  # Print the number of columns after creating dummy variables
  print(paste("Number of columns after creating dummies:", ncol(data_with_dummies)))
  
  # Print dummy variable column names
  print("Newly created dummy variable column names:")
  print(dummy_colnames)
  
  return(data_with_dummies)
}

# Function to handle multicollinearity by checking single-class variables, sparse dummy variables, 
# highly correlated covariate pairs, and combining rare levels; optionally performs VIF analysis
handle_multicollinearity <- function(data, dummy_vars, rare_level_threshold = 5, log = FALSE) {
  log_content <- c()
  
  # Check for single-class variables
  single_class_vars <- data %>%
    dplyr::select(FID, IID, starts_with("mri_info_"), starts_with("batch")) %>%
    sapply(function(col) length(unique(col)) == 1)
  
  if (any(single_class_vars)) {
    single_class_var_names <- names(single_class_vars)[single_class_vars]
    warning_message <- paste("The following variables have only one class and will be removed:", 
                             toString(single_class_var_names))
    message(warning_message)
    log_content <- c(log_content, warning_message)
    data <- data %>%
      dplyr::select(-one_of(single_class_var_names))
  }
  
  # Check for sparse dummy variables
  sparse_threshold <- 0.05
  dummy_cols <- data %>% dplyr::select(starts_with("mri_info_"), starts_with("batch"))
  sparse_dummies <- colSums(dummy_cols == 1) / nrow(dummy_cols) < sparse_threshold
  
  if (log && any(sparse_dummies)) {
    sparse_message <- paste("Sparse dummy variables (less than", sparse_threshold * 100, 
                            "% of 1s):", toString(names(sparse_dummies)[sparse_dummies]))
    message(sparse_message)
    log_content <- c(log_content, sparse_message)
  }
  
  # Check for multicollinearity
  covariates <- data %>%
    dplyr::select(starts_with("mri_info_"), starts_with("batch"), starts_with("PC"), interview_age)
  cor_matrix <- cor(covariates)
  
  # Identify highly correlated pairs
  high_cor <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if (log && nrow(high_cor) > 0) {
    cor_message <- paste("Highly correlated covariate pairs:", 
                         paste(rownames(cor_matrix)[high_cor[,1]], 
                               colnames(cor_matrix)[high_cor[,2]], 
                               sep = " vs. ", collapse = ", "))
    message(cor_message)
    log_content <- c(log_content, cor_message)
  }
  
  # Address multicollinearity by combining levels of mri_info_deviceserialnumber
  if (any(sparse_dummies) || nrow(high_cor) > 0) {
    mri_counts <- table(data$mri_info_deviceserialnumber)
    rare_levels <- names(mri_counts[mri_counts < rare_level_threshold])
    if (log) {
      combine_message <- paste("Rare levels combined into 'Other':", toString(rare_levels))
      message(combine_message)
      log_content <- c(log_content, combine_message)
    }
    data$mri_info_deviceserialnumber[data$mri_info_deviceserialnumber %in% rare_levels] <- "Other"
    
    # Re-calculate dummy variables after combination
    data <- dummy_cols(data, select_columns = dummy_vars, 
                       remove_first_dummy = TRUE, remove_selected_columns = TRUE)
    
    # Recalculate covariates and correlation matrix
    covariates <- data %>%
      dplyr::select(starts_with("mri_info_"), starts_with("batch"), starts_with("PC"), interview_age)
    cor_matrix <- cor(covariates)
    
    if (log) {
      combine_message <- "Combined levels of mri_info_deviceserialnumber based on sparsity and correlation."
      message(combine_message)
      log_content <- c(log_content, combine_message)
    }
  }
  
  # Optional: VIF Analysis
  if (requireNamespace("car", quietly = TRUE)) {
    vif_values <- car::vif(lm(interview_age ~ ., data = covariates))
    high_vif <- vif_values > 5
    
    if (log && any(high_vif)) {
      vif_message <- paste("High VIF values (above 5):", toString(names(vif_values)[high_vif]))
      message(vif_message)
      log_content <- c(log_content, vif_message)
    }
    
    if (any(high_vif)) {
      data <- data %>% dplyr::select(-one_of(names(vif_values)[high_vif]))
    }
  }
  
  list(data = data, log_content = log_content)
}

# Function to create and save split data and apply FRGEpistasis rank inverse log normalization 
# with directory check/creation
save_split_data <- function(data, ethnicity, sex, pheno_dir, covar_dir, date) {
  
  # Subset data for the specific ethnicity and sex
  subset_data <- data %>%
    filter(ethnicity == !!ethnicity, sex == !!sex) %>%
    mutate(across(starts_with("smri_vol_"), ~ rankTransPheno(.x, 0.5)))
  
  # Define output directories
  pheno_out_dir <- file.path(pheno_dir, ethnicity, sex)
  covar_disc_out_dir <- file.path(covar_dir, "Discrete", ethnicity, sex)
  covar_quant_out_dir <- file.path(covar_dir, "Quantitative", ethnicity, sex)
  
  # Verify the paths exist, and if not, create them
  if (!dir.exists(pheno_out_dir)) dir.create(pheno_out_dir, recursive = TRUE)
  if (!dir.exists(covar_disc_out_dir)) dir.create(covar_disc_out_dir, recursive = TRUE)
  if (!dir.exists(covar_quant_out_dir)) dir.create(covar_quant_out_dir, recursive = TRUE)
  
  # Generate random ID for each file set
  id <- sample(1:1000000, 1)
  
  # Get the number of samples for the subset
  num_samples <- nrow(subset_data)
  
  # Save phenotype files and create covariate files
  phenotypes <- colnames(subset_data)[grepl("^smri_vol", colnames(subset_data))]
  for (phenotype in phenotypes) {
    if (phenotype != "smri_vol_scs_intracranialv_ROC0_2") {
      # Save phenotype files for non-ICV volumes
      phenotype_file_name <- file.path(pheno_out_dir, sprintf("%d_pheno_%s_%s_%s_%d_%s.txt",
                                                              id, date, ethnicity, sex, num_samples, phenotype))
      write.table(subset_data %>% dplyr::select(FID, IID, all_of(phenotype)),
                  file = phenotype_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
    
    # Create qcovar file based on whether the phenotype is "whole brain" or not
    if (phenotype == "smri_vol_scs_wholeb_ROC0_2") {
      # Create qcovar file that excludes ICV for whole brain phenotype
      qcovar_file_name_no_icv <- file.path(covar_quant_out_dir, sprintf("qcovar_noICV_%d_%s_%s_%s_%d.txt", 
                                                                        id, date, ethnicity, sex, num_samples))
      write.table(subset_data %>% dplyr::select(FID, IID, interview_age, starts_with("PC")), 
                  file = qcovar_file_name_no_icv, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    } else {
      # Create qcovar file that includes ICV for other phenotypes
      qcovar_file_name <- file.path(covar_quant_out_dir, sprintf("qcovar_%d_%s_%s_%s_%d.txt", 
                                                                 id, date, ethnicity, sex, num_samples))
      write.table(subset_data %>% dplyr::select(FID, IID, interview_age, starts_with("PC"), smri_vol_scs_intracranialv_ROC0_2),
                  file = qcovar_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  
  # Save discrete covariate file
  covar_discrete <- subset_data %>%
    dplyr::select(FID, IID, starts_with("mri_info_"), starts_with("batch"))
  covar_file_name <- file.path(covar_disc_out_dir, sprintf("covar_%d_%s_%s_%s_%d.txt",
                                                           id, date, ethnicity, sex, num_samples))
  write.table(covar_discrete, file = covar_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

# Function to create and save split data, convert discrete covar dummy variables, 
# and apply FRGEpistasis rank inverse log normalization with directory check/creation
# save_split_data <- function(data, ethnicity, sex, pheno_dir, covar_dir, date, dummy_vars, rare_level_threshold = 5, log = FALSE) {
#   
#   # Create dummy variables on the entire dataset using fastDummies package
#   data_with_dummies <- dummy_cols(data, select_columns = dummy_vars, 
#                                   remove_first_dummy = TRUE, remove_selected_columns = TRUE)
#   
#   # Define output directories
#   pheno_out_dir <- file.path(pheno_dir, ethnicity, sex)
#   covar_disc_out_dir <- file.path(covar_dir, "Discrete", ethnicity, sex)
#   covar_quant_out_dir <- file.path(covar_dir, "Quantitative", ethnicity, sex)
#   
#   # Verify the paths exist
#   if (!dir.exists(pheno_out_dir)) dir.create(pheno_out_dir, recursive = TRUE)
#   if (!dir.exists(covar_disc_out_dir)) dir.create(covar_disc_out_dir, recursive = TRUE)
#   if (!dir.exists(covar_quant_out_dir)) dir.create(covar_quant_out_dir, recursive = TRUE)
#   
#   # Subset data and apply FRGEpistasis rank inverse log transformation
#   subset_data <- data_with_dummies %>%
#     filter(ethnicity == !!ethnicity, sex == !!sex) %>%
#     mutate(across(starts_with("smri_vol_"), ~ rankTransPheno(.x, 0.5)))
#   
#   # Generate random ID for each file
#   id <- sample(1:1000000, 1)
#   
#   # Handle multicollinearity
#   result <- handle_multicollinearity(subset_data, dummy_vars, rare_level_threshold, log)
#   subset_data <- result$data
#   log_content <- result$log_content
#   
#   # Log the information if logging is enabled
#   if (log) {
#     original_colnames <- colnames(data)
#     new_colnames <- colnames(subset_data)
#     dummy_colnames <- setdiff(new_colnames, original_colnames)
#     
#     dummy_log_content <- c(
#       paste("Number of new dummy columns for", ethnicity, sex, ":", length(dummy_colnames)),
#       "Names of new dummy columns:",
#       toString(dummy_colnames)
#     )
#     
#     log_content <- c(log_content, dummy_log_content)
#     
#     log_file <- file.path(covar_disc_out_dir, sprintf("%d_log_%s_%s_%s.txt", id, date, ethnicity, sex))
#     writeLines(log_content, log_file)
#   }
#   
#   # Get number of samples for split subset
#   num_samples <- nrow(subset_data)
#   
#   # Save phenotype files
#   phenotypes <- colnames(subset_data)[grepl("^smri_vol", colnames(subset_data))]
#   for (phenotype in phenotypes) {
#     if (phenotype != "smri_vol_scs_intracranialv_ROC0_2") {
#       phenotype_file_name <- file.path(pheno_out_dir, sprintf("%d_pheno_%s_%s_%s_%d_%s.txt",
#                                                               id, date, ethnicity, sex, num_samples, phenotype))
#       write.table(subset_data %>% dplyr::select(FID, IID, all_of(phenotype)),
#                   file = phenotype_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#     }
#     
#     if (phenotype == "smri_vol_scs_wholeb_ROC0_2") {
#       qcovar_file_name_no_icv <- file.path(covar_quant_out_dir, sprintf("qcovar_noICV_%d_%s_%s_%s_%d.txt", 
#                                                                         id, date, ethnicity, sex, num_samples))
#       write.table(subset_data %>% dplyr::select(FID, IID, interview_age, starts_with("PC")), 
#                   file = qcovar_file_name_no_icv, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#     } else {
#       qcovar_file_name <- file.path(covar_quant_out_dir, sprintf("qcovar_%d_%s_%s_%s_%d.txt", 
#                                                                  id, date, ethnicity, sex, num_samples))
#       write.table(subset_data %>% dplyr::select(FID, IID, interview_age, starts_with("PC"), smri_vol_scs_intracranialv_ROC0_2),
#                   file = qcovar_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#     }
#   }
#   
#   # Save discrete covariate file with dummy variables for specified columns
#   covar_discrete <- subset_data %>%
#     dplyr::select(FID, IID, starts_with("mri_info_"), starts_with("batch"))
#   covar_file_name <- file.path(covar_disc_out_dir, sprintf("covar_%d_%s_%s_%s_%d.txt",
#                                                            id, date, ethnicity, sex, num_samples))
#   write.table(covar_discrete, file = covar_file_name, col.names = FALSE,
#               row.names = FALSE, sep = "\t", quote = FALSE)
# }

# Function to archive or clear files in 'F' and 'M' subdirectories
clean_old_runs <- function(base_path, archive = TRUE) {
  # Get all subdirectories recursively
  dirs <- list.dirs(path = base_path, full.names = TRUE, recursive = TRUE)
  
  # Loop through each directory
  for (dir in dirs) {
    # Only clear files in 'F' and 'M' subdirectories
    if (basename(dir) %in% c("F", "M")) {
      files <- list.files(path = dir, full.names = TRUE)
      
      # Exclude files in the archive subdirectory
      files <- files[!grepl("/archive/", files)]
      
      if (length(files) > 0) {
        if (archive) {
          # Create archive directory if it doesn't exist
          archive_dir <- file.path(dir, "archive")
          dir.create(archive_dir, showWarnings = FALSE)
          # Move files to archive directory
          file.copy(files, archive_dir, overwrite = TRUE)
        }
        # Remove the original files
        unlink(files)
      }
    }
  }
}

# Function to create visualizations (histogram and QQ plot) for normality check
plot_normality <- function(data, original_data, phenotype_cols) {
  plot_list <- list()
  for (col_name in phenotype_cols) {
    # Plot for transformed data
    p1 <- ggplot(data, aes(.data[[col_name]])) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "blue", alpha = 0.7) +
      geom_density(color = "red") +
      labs(title = paste("rankTransPheno Histogram", col_name))
    
    p2 <- ggplot(data, aes(sample = .data[[col_name]])) +
      stat_qq() +
      stat_qq_line() +
      labs(title = paste("rankTransPheno QQ Plot"))
    
    # Plot for original data
    p3 <- ggplot(original_data, aes(.data[[col_name]])) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "blue", alpha = 0.7) +
      geom_density(color = "red") +
      labs(title = paste("Original Histogram"))
    
    p4 <- ggplot(original_data, aes(sample = .data[[col_name]])) +
      stat_qq() +
      stat_qq_line() +
      labs(title = paste("Original QQ Plot"))
    
    combined_histogram <- p3 + p1
    combined_qqplot <- p4 + p2
    
    plot_list[[col_name]] <- list(histogram = combined_histogram, qq_plot = combined_qqplot)
  }
  return(plot_list)
}

# Function to apply statistical test for normality
test_normality <- function(data, phenotype_cols) {
  test_results <- tibble(
    Phenotype = phenotype_cols,
    P_value = map_dbl(phenotype_cols, ~ ad.test(data[[.x]])$p.value)
  )
  return(test_results)
}

