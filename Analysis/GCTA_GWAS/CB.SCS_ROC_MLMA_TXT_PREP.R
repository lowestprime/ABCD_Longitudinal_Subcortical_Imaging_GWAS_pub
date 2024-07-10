#### GCTA GWAS PREP ####

# Load Packages
library(purrr)
library(readr)
library(FRGEpistasis)
library(dplyr)
library(data.table)

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
date <- format(Sys.Date(), "%Y%m%d")

# Load external functions
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

# Read and process covariate data
covariate_data <- read.table(covariate_file, header = TRUE)
colnames(covariate_data) <- c("FID", "IID", "batch", "sex")

# Filter covariate data to match phenotypes
filtered_covariate_data <- covariate_data %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")
filtered_ancestry_pcs <- ancestry_pcs %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")

# Merge phenotype, covariate, and ancestry PC data
merged_data <- gcta.pheno.scs.vol.roc %>%
  left_join(filtered_covariate_data %>% select(IID, batch), by = "IID") %>%
  left_join(filtered_ancestry_pcs %>% select(IID, starts_with("PC")), by = "IID")

# Remove NAs and reorder columns
merged_data_final <- merged_data %>%
  filter(!is.na(batch), !is.na(PC1), !is.na(PC10)) %>%
  select(FID, IID, everything())

# Assign ethnicity based on IIDs
afr_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.AFR.2263.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]
amr_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.AMR.2019.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]
eur_iids <- read.table(file.path(anc_pc_dir, "ABCD.ancestry_knn.EUR.6891.txt"), header = FALSE, stringsAsFactors = FALSE)[, 1]

merged_data_final$ethnicity <- NA
merged_data_final$ethnicity[merged_data_final$IID %in% afr_iids] <- "AFR"
merged_data_final$ethnicity[merged_data_final$IID %in% amr_iids] <- "AMR"
merged_data_final$ethnicity[merged_data_final$IID %in% eur_iids] <- "EUR"

# Remove samples with missing ethnicity and reorder columns
merged_data_final <- merged_data_final %>%
  filter(!is.na(ethnicity)) %>%
  relocate(batch, ethnicity, PC1:PC10, .after = interview_age)

# Apply the function to each combination of ethnicity and sex
for (ethnicity in ethnicities) {
  for (sex in sexes) {
    save_split_data(merged_data_final, ethnicity, sex, pheno_dir, covar_dir, date)
  }
}