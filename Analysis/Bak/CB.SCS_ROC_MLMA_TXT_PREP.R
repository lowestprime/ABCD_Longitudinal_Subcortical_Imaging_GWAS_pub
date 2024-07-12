#### GCTA GWAS PREP ####
# Bak
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
p_load(purrr, readr, FRGEpistasis, dplyr)

# Directories
base_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/'
r5.1_genetics <- '/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5.1/core/genetics'
image_dir <- paste0(base_dir, 'Images/')
shiny_dir <- paste0(base_dir, 'plots/Shinyapps_plotly_plots/')
table1_dir <- paste0(shiny_dir, 'Counts.Table')
table2_dir <- paste0(shiny_dir, 'ROC.Summary.Table')
gwas_dir <- paste0(base_dir, 'GCTA_GWAS/')
proc_dir <- paste0(gwas_dir, 'Processed_Data')
pheno_dir <- file.path(proc_dir, "Phenotypes")
covar_dir <- file.path(proc_dir, "Covariates")
anc_pc_dir <- paste0(gwas_dir, 'ANCESTRY_PCS/')

# Define ancestries and sexes
ethnicities <- c("AFR", "AMR", "EUR")
sexes <- c("F", "M")

# Set Date
date <- format(Sys.Date(), "%m%d%Y")

# Load ABCD 5.1 external functions and Sruthi covar file
source(paste0(base_dir, 'release5_external_functions.R'))
covariate_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/covar_AllSubj_batch_gender_noNAs_baseline_11665.txt"

# remove non-ROI columns from smri.R5.1.baseline.y2.ROC and reorder columns
smri.R5.1.baseline.y2.ROC.filtered <- roi_filter(smri.R5.1.baseline.y2.ROC)

# rename pheno cols to GCTA format
gcta.pheno.scs.vol.roc <- smri.R5.1.baseline.y2.ROC.filtered %>%
  rename(FID = rel_family_id, IID = src_subject_id)

# Load and preprocess ancestry principal components data by renaming and selecting columns
ancestry_pcs <- fread(paste0(anc_pc_dir, "ABCD5_all_ancestries_all_individuals_PC20.txt")) %>%
  rename(IID = V1) %>%
  dplyr::select(-V2)

# Apply the renaming function to the ancestry_pcs data frame
ancestry_pcs <- ancestry_pcs %>% rename_with(rename_columns, .cols = starts_with("V"))

#### Remove samples with missing data ####
# Read the files without headers
covariate_data <- read.table(covariate_file, header = FALSE)

# Assign appropriate column names for phenotype_data and covariate_data
colnames(covariate_data) <- c("FID", "IID", "batch", "sex")

# Merge gcta.pheno.scs.vol.roc with covariate_data genotyping batch and ancestry_pcs 1-20
# Filter covariate_data and ancestry_pcs to include only the IIDs present in gcta.pheno.scs.vol.roc
filtered_covariate_data <- covariate_data %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")
filtered_ancestry_pcs <- ancestry_pcs %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")

# Merge gcta.pheno.scs.vol.roc with filtered_covariate_data and filtered_ancestry_pcs
merged_data <- gcta.pheno.scs.vol.roc %>%
  left_join(filtered_covariate_data %>% dplyr::select(IID, batch), by = "IID") %>%
  left_join(filtered_ancestry_pcs %>% dplyr::select(IID, starts_with("PC")), by = "IID")

# drop any rows from merged_data with NAs (in batch and PCs columns)
# reorder columns
merged_data_no_na <- merged_data %>%
  dplyr::select(FID, IID, everything()) %>% # Ensures FID and IID are the first two columns
  drop_na()

# replace FID = rel_family_id, IID = src_subject_id, remove PC11-PC20 and reported ethnicity cols
merged_data_final <- merged_data_no_na %>%
  mutate(FID = IID) %>%
  dplyr::select(-PC11:-PC20, -ethnicity)

# split merged_data_final by TOPMed ethnicity
# Define the file names for each ethnicity
afr_file <- file.path(anc_pc_dir, "ABCD.ancestry_knn.AFR.2263.txt")
amr_file <- file.path(anc_pc_dir, "ABCD.ancestry_knn.AMR.2019.txt")
eur_file <- file.path(anc_pc_dir, "ABCD.ancestry_knn.EUR.6891.txt")

# Read the IIDs from each file
afr_iids <- read.table(afr_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
amr_iids <- read.table(amr_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
eur_iids <- read.table(eur_file, header = FALSE, stringsAsFactors = FALSE)[, 1]

# make empty ethnicity column
merged_data_final$ethnicity <- NA

# Assign ethnicity based on IIDs
merged_data_final$ethnicity[merged_data_final$IID %in% afr_iids] <- "AFR"
merged_data_final$ethnicity[merged_data_final$IID %in% amr_iids] <- "AMR"
merged_data_final$ethnicity[merged_data_final$IID %in% eur_iids] <- "EUR"

# Remove samples with missing ethnicity and reorder cols
merged_data_final <- merged_data_final %>%
  filter(!is.na(ethnicity)) %>%
  relocate(batch, ethnicity, PC1:PC10, .after = interview_age)

# Split pheno, covar and qcovar files by ethnicity and sex and then save
# for (ethnicity in ethnicities) {
#   for (sex in sexes) {
#     # Filter data by ethnicity and sex
#     df_subset <- filter_data(merged_data_final, ethnicity, sex)
#     
#     # Apply rank inverse normalization to phenotypes
#     df_subset <- apply_rank_inverse_norm(df_subset)
#     
#     # Save GCTA files
#     save_gcta_files(df_subset, ethnicity, sex, pheno_dir, covar_dir, date)
#   }
# }
for (ethnicity in ethnicities) {
  for (sex in sexes) {
    save_split_data(merged_data_final, ethnicity, sex, pheno_dir, covar_dir, date)
  }
}


#### Extra and QC ####
# Subset the merged_data_final dataframe by ethnicity
# afr_data <- merged_data_final[merged_data_final$IID %in% afr_iids, ]
# amr_data <- merged_data_final[merged_data_final$IID %in% amr_iids, ]
# eur_data <- merged_data_final[merged_data_final$IID %in% eur_iids, ]

# Optional: Check the number of rows in each subset to ensure correctness
# cat("Total number of individuals:", nrow(merged_data_final), "\n")
# cat("Number of AFR individuals:", nrow(afr_data), "\n")
# cat("Number of AMR individuals:", nrow(amr_data), "\n")
# cat("Number of EUR individuals:", nrow(eur_data), "\n")
# Create a new ethnicity column in merged_data_final

# Save final_data as a space-separated text file suitable for GCTA MLMA
# setwd(pheno_dir)
# write.table(merged_data_final, file = "gcta_mlma_SCS_ROC_master.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
# first_100_rows <- head(merged_data_final, 100)
# write.csv(first_100_rows, "first_100_rows.csv", row.names = FALSE)

# Filtering rows with any `NA` values
# rows_with_na <- merged_data %>% 
#   filter(if_any(everything(), is.na)) %>%
#   dplyr::select(FID, IID, batch, PC1:PC10)

# Filter out rows in `rows_with_na` where `IID` is not present in either `ancestry_pcs$IID` or `covariate_data$IID`
# This results in a dataset with `IID`s that are missing in either data frame
# all_IIDs_missing <- rows_with_na %>%
#   filter(!(IID %in% ancestry_pcs$IID) | !(IID %in% covariate_data$IID))

# Extract unique IIDs and save to a text file
# unique_IIDs <- all_IIDs_missing %>% 
#   pull(IID) %>%
#   unique()

# Write IIDs to a text file
# write_lines(unique_IIDs, "missing_IIDs.txt")

# Filter `smri.R5.1.baseline.y2` to retain rows with `src_subject_id` present in `all_IIDs_missing$IID`
# smri.R5.1.baseline.y2_missing_baseline_pcs <- smri.R5.1.baseline.y2 %>%
#   filter(src_subject_id %in% all_IIDs_missing$IID) %>%
#   distinct(src_subject_id, .keep_all = TRUE)

# Save merged_data_no_na as a space-separated text file
# write.table(smri.R5.1.baseline.y2_missing_baseline_pcs, file = "missing_IIDs_info.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Covariate file
# covariate_data <- final_data %>%
#   dplyr::select(src_subject_id, sex, batch, mri_info_deviceserialnumber) %>%
#   rename(IID = src_subject_id)

# save covariate file
# write.table(covariate_data, file = "covariate.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Quantitative Covariate file

# save phenotype files
# Assuming gcta.pheno.scs.vol.roc is your data frame and pheno_dir is the directory where you want to save the files

#### From EMMA/Sruthi ####
# Read ancestry ID files
# read_ancestry_ids <- function(file_path) {
#   read.table(file_path, header = FALSE, stringsAsFactors = FALSE)$V2
# }
# 
# ancestry_files <- list(
#   AFR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AFR.2263_no.sexmismatch_IDs.txt",
#   AMR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AMR.2019_no.sexmismatch_IDs.txt",
#   EUR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_EUR.6891_no.sexmismatch_IDs.txt"
# )
# 
# ancestry_samples_id_list <- map(ancestry_files, read_ancestry_ids)

# Split phenotype data by ancestry
# split_by_ancestry <- function(pheno_data, ancestry_samples_id_list) {
#   map2(ancestry_samples_id_list, names(ancestry_samples_id_list), function(ids, ancestry_name) {
#     pheno_data %>% filter(IID %in% ids) %>% mutate(ancestry = ancestry_name)
#   })
# }

# Split phenotype data by ancestry and gender
# split_by_ancestry_and_gender <- function(pheno_data_by_ancestry) {
#   gender <- c("M", "F")
#   pheno_data_by_ancestry_and_gender <- list()
#   
#   for (i in seq_along(pheno_data_by_ancestry)) {
#     pheno_data <- pheno_data_by_ancestry[[i]]
#     pheno_name <- tolower(names(pheno_data_by_ancestry[i]))
#     
#     for (j in seq_along(gender)) {
#       name <- paste(pheno_name, tolower(gender[j]), sep = "_")
#       pheno_data_by_ancestry_and_gender[[name]] <- pheno_data %>% filter(sex == gender[j])
#     }
#   }
#   
#   return(pheno_data_by_ancestry_and_gender)
# }

# Save phenotype tables with ancestry info
# save_pheno_tables <- function(pheno_data_list, output_dir, date) {
#   setwd(output_dir)
#   
#   walk2(pheno_data_list, names(pheno_data_list), function(dt, pheno_name) {
#     n <- nrow(dt)
#     file_name <- paste0(pheno_name, "_1yr_followup_before_srs_within_ancestry_group_noNAs_", n, "_table_", date, ".csv")
#     write.csv(dt, file_name, row.names = TRUE)
#   })
# }

# Main processing function
# process_phenotype_data <- function(pheno_data, ancestry_samples_id_list, output_dir, date) {
#   # Split by ancestry
#   pheno_by_ancestry <- split_by_ancestry(pheno_data, ancestry_samples_id_list)
#   
#   # Split by ancestry and gender
#   pheno_by_ancestry_and_gender <- split_by_ancestry_and_gender(pheno_by_ancestry)
#   
#   # Combine both lists
#   data <- c(pheno_by_ancestry, pheno_by_ancestry_and_gender)
#   
#   # Save tables
#   save_pheno_tables(data, output_dir, date)
# }

# Usage
# pheno <- gcta.pheno.scs.vol.roc  # Your phenotype data
# process_phenotype_data(pheno, ancestry_samples_id_list, pheno_qc_dir, date)
