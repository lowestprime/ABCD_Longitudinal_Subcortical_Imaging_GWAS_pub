#### GCTA GWAS PREP ####
# 7/11/2024
# Phenotypes:
# roi_columns <- smri.R5.1.baseline.y2.ROC.filtered %>% 
#   dplyr::select(starts_with("smri_vol_")) %>% 
#   colnames()

## Discrete Covariates ##
# sex, batch, mri_info_deviceserialnumber

## Quantitative Covariates ##
# interview_age, bigsnpr top 10 PCs, smri_vol_scs_intracranialv (except for smri_vol_scs_wholeb)

## Priorities ##
# Ethnicity: EUR
# ROI: smri_vol_scs_wholeb (smri_vol_scs_intracranialv covar not needed)

# Load Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(purrr, readr, FRGEpistasis, dplyr, fastDummies, nortest, patchwork, data.table)

# Directories
base_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/'
r5.1_genetics <- '/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5.1/core/genetics'
image_dir <- paste0(base_dir, 'Images/')
shiny_dir <- paste0(base_dir, 'plots/Shinyapps_plotly_plots/')
table1_dir <- paste0(shiny_dir, 'Counts.Table')
table2_dir <- paste0(shiny_dir, 'ROC.Summary.Table')
gwas_dir <- paste0(base_dir, 'GCTA_GWAS/')
proc_dir <- paste0(gwas_dir, 'Processed_Data')
pheno_dir <- file.path(proc_dir, 'Phenotypes')
covar_dir <- file.path(proc_dir, 'Covariates')
anc_pc_dir <- paste0(gwas_dir, 'ANCESTRY_PCS/')

# Define ancestries and sexes
ethnicities <- c("AFR", "AMR", "EUR")
sexes <- c("F", "M")

# Specify the columns for dummy variable creation
dummy_vars <- c("mri_info_deviceserialnumber", "batch")

# Set Date
date <- format(Sys.Date(), "%m%d%Y")

# Load ABCD 5.1 external functions and Sruthi covar file
source(paste0(base_dir, 'release5_external_functions.R'))
covariate_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/covar_AllSubj_batch_gender_noNAs_baseline_11665.txt"

# Remove non-ROI columns from smri.R5.1.baseline.y2.ROC and reorder columns
smri.R5.1.baseline.y2.ROC.filtered <- roi_filter(smri.R5.1.baseline.y2.ROC)

# Rename pheno cols to GCTA format
gcta.pheno.scs.vol.roc <- smri.R5.1.baseline.y2.ROC.filtered %>%
  rename(FID = rel_family_id, IID = src_subject_id)

# Load and preprocess ancestry principal components data by renaming and selecting columns
ancestry_pcs <- fread(paste0(anc_pc_dir, "ABCD5_all_ancestries_all_individuals_PC20.txt")) %>%
  rename(IID = V1) %>%
  dplyr::select(-V2)

# Apply the renaming function to the ancestry_pcs data frame
ancestry_pcs <- ancestry_pcs %>% rename_with(rename_columns, .cols = starts_with("V"))

#### Remove samples with missing data and create split txts ####
# Read files without headers and assign appropriate column names for phenotype_data and covariate_data
covariate_data <- read.table(covariate_file, header = FALSE)
colnames(covariate_data) <- c("FID", "IID", "batch", "sex")

# Merge gcta.pheno.scs.vol.roc with covariate_data (genotyping batch) and ancestry_pcs (PCs 1-20), then filter by matching IIDs
filtered_covariate_data <- covariate_data %>% 
  semi_join(gcta.pheno.scs.vol.roc, by = "IID")
filtered_ancestry_pcs <- ancestry_pcs %>% 
  semi_join(gcta.pheno.scs.vol.roc, by = "IID")

# Merge phenotype, covariate, and ancestry PC data
merged_data <- gcta.pheno.scs.vol.roc %>%
  left_join(filtered_covariate_data %>% dplyr::select(IID, batch), by = "IID") %>%
  left_join(filtered_ancestry_pcs %>% dplyr::select(IID, starts_with("PC")), by = "IID")

merged_data_final <- merged_data %>%
  # Ensure FID and IID are the first two columns, then remove rows with any NA values
  dplyr::select(FID, IID, everything()) %>%
  drop_na() %>%
  # Replace FID with IID values
  mutate(FID = IID) %>%
  # Remove PC11 to PC20 and ethnicity columns
  dplyr::select(-PC11:-PC20, -ethnicity)

# Read the IIDs from each file
afr_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.AFR.2263.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]
amr_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.AMR.2019.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]
eur_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.EUR.6891.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]

# Make empty ethnicity column
merged_data_final$ethnicity <- NA

# Assign ethnicity based on IIDs
merged_data_final$ethnicity[merged_data_final$IID %in% afr_iids] <- "AFR"
merged_data_final$ethnicity[merged_data_final$IID %in% amr_iids] <- "AMR"
merged_data_final$ethnicity[merged_data_final$IID %in% eur_iids] <- "EUR"

# Remove samples with missing ethnicity and reorder columns
merged_data_final <- merged_data_final %>%
  filter(!is.na(ethnicity)) %>%
  relocate(batch, ethnicity, PC1:PC10, .after = interview_age)

# Apply save_split_data to each combination of ethnicity and sex
for (ethnicity in ethnicities) {
  for (sex in sexes) {
    save_split_data(merged_data_final, ethnicity, sex, pheno_dir, covar_dir, date, dummy_vars)
  }
}

# Generate normality check visualizations (histogram and QQ plot) and tests for each split phenotype
for (ethnicity in ethnicities) {
  for (sex in sexes) {
    # Define plot and test output directories
    plot_out_dir <- file.path(pheno_dir, "Plots", ethnicity, sex)
    test_out_dir <- file.path(pheno_dir, "Tests", ethnicity, sex)
    
    # Verify the paths exist
    if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)
    if (!dir.exists(test_out_dir)) dir.create(test_out_dir, recursive = TRUE)
    
    # Create subset_data_normalized to check for normality
    subset_data <- merged_data_final %>%
      filter(ethnicity == !!ethnicity, sex == !!sex)
    subset_data_normalized <- subset_data %>%
      mutate(across(starts_with("smri_vol_"), ~ rankTransPheno(.x, 0.5)))
    
    # Specify which columns to check for normality
    phenotype_cols <- colnames(subset_data_normalized)[grepl("^smri_vol_", colnames(subset_data_normalized))]
    
    # Create visualizations (histogram and QQ plot) for normality check
    plot_list <- plot_normality(subset_data_normalized, subset_data, phenotype_cols)
    for (col_name in names(plot_list)) {
      ggsave(filename = file.path(plot_out_dir, sprintf("%s_%s_%s_hist.png", ethnicity, sex, col_name)), 
             plot = plot_list[[col_name]]$histogram)
      ggsave(filename = file.path(plot_out_dir, sprintf("%s_%s_%s_qq.png", ethnicity, sex, col_name)), 
             plot = plot_list[[col_name]]$qq_plot)
    }
    
    # Statistical Test
    test_results_old <- test_normality(subset_data, phenotype_cols)
    print(test_results_old)
    test_results_new <- test_normality(subset_data_normalized, phenotype_cols)
    print(test_results_new)
    
    # Save the test results to a file
    write_csv(test_results_new, file.path(test_out_dir, sprintf("%s_%s_test_results_new.csv", ethnicity, sex)))
    write_csv(test_results_old, file.path(test_out_dir, sprintf("%s_%s_test_results_old.csv", ethnicity, sex)))
  }
}

# Run to delete all files in split processed_Data (proc_dir)
# clear_files_in_FM_subdirectories(proc_dir)
