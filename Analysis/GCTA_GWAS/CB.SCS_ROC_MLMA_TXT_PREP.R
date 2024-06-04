#### GCTA MLMA GWAS PREP ####
## phenotypes, covars, qcovars ##
# Phenotypes:
# roi_columns <- smri.R5.1.baseline.y2.ROC.filtered %>% 
#   select(starts_with("smri_vol_")) %>% 
#   colnames()

# Discrete Covariates:
# sex, genotyping batch, mri_info_deviceserialnumber

# Quantitative Covariates:
# interview_age, bigsnpr 20 PCs, smri_vol_scs_intracranialv (except for smri_vol_scs_wholeb)

## Remaining Tasks ##
# append genotyping batch ancestry pc covars
# split txts by sex ancestry and phenotype
# add rank based log transformation to normalize the wide phenos

## Priorities ##
# Ethnicity priority: EUR
# ROI priority: smri_vol_scs_wholeb (smri_vol_scs_intracranialv covar not needed)

library(pacman)
p_load(dplyr, purrr, readr)

# directories
base_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/'
image_dir <- paste0(base_dir, 'Images/')
shiny_dir <- paste0(base_dir, 'plots/Shinyapps_plotly_plots/')
table1_dir <- paste0(shiny_dir, 'Counts.Table')
table2_dir <- paste0(shiny_dir, 'ROC.Summary.Table')
gwas_dir <- paste0(base_dir, 'GCTA_GWAS/')
anc_pc_dir <- paste0(gwas_dir, 'ANCESTRY_PCS/')
pheno_dir <- paste0(gwas_dir, 'PHENOTYPE_DATA')
ancestry_id_dir <- "/path/to/ancestry/ids"
pheno_qc_dir <- "/path/to/output/directory"
date <- Sys.Date()

# ABCD 5.1 data and external functions
source(paste0(base_dir, 'release5_external_functions.R'))

# remove non-ROI columns from smri.R5.1.baseline.y2.ROC and reorder columns
smri.R5.1.baseline.y2.ROC.filtered <- roi_filter(smri.R5.1.baseline.y2.ROC)
# pivot long
smri.R5.1.baseline.y2.ROC.filtered.long <- pivot_roc_to_long_format(smri.R5.1.baseline.y2.ROC.filtered, is_baseline_y2 = TRUE)

# rename pheno cols to GCTA format
gcta.pheno.scs.vol.roc <- smri.R5.1.baseline.y2.ROC.filtered %>%
  rename(FID = rel_family_id, IID = src_subject_id)

# save phenotype files
# Assuming gcta.pheno.scs.vol.roc is your data frame and pheno_dir is the directory where you want to save the files
create_phenotype_files(gcta.pheno.scs.vol.roc, pheno_dir)

# Define paths


# Read ancestry ID files
read_ancestry_ids <- function(file_path) {
  read.table(file_path, header = FALSE, stringsAsFactors = FALSE)$V2
}

ancestry_files <- list(
  AFR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AFR.2263_no.sexmismatch_IDs.txt",
  AMR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AMR.2019_no.sexmismatch_IDs.txt",
  # EAS = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_EAS.382_no.sexmismatch_IDs.txt",
  EUR = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_EUR.6891_no.sexmismatch_IDs.txt"
  # SAS = "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_SAS.92_no.sexmismatch_IDs.txt"
)

ancestry_samples_id_list <- map(ancestry_files, read_ancestry_ids)

# Split phenotype data by ancestry
split_by_ancestry <- function(pheno_data, ancestry_samples_id_list) {
  map2(ancestry_samples_id_list, names(ancestry_samples_id_list), function(ids, ancestry_name) {
    pheno_data %>% filter(IID %in% ids) %>% mutate(ancestry = ancestry_name)
  })
}

# Split phenotype data by ancestry and gender
split_by_ancestry_and_gender <- function(pheno_data_by_ancestry) {
  gender <- c("M", "F")
  pheno_data_by_ancestry_and_gender <- list()
  
  for (i in seq_along(pheno_data_by_ancestry)) {
    pheno_data <- pheno_data_by_ancestry[[i]]
    pheno_name <- tolower(names(pheno_data_by_ancestry[i]))
    
    for (j in seq_along(gender)) {
      name <- paste(pheno_name, tolower(gender[j]), sep = "_")
      pheno_data_by_ancestry_and_gender[[name]] <- pheno_data %>% filter(sex == gender[j])
    }
  }
  
  return(pheno_data_by_ancestry_and_gender)
}

# Save phenotype tables with ancestry info
save_pheno_tables <- function(pheno_data_list, output_dir, date) {
  setwd(output_dir)
  
  walk2(pheno_data_list, names(pheno_data_list), function(dt, pheno_name) {
    n <- nrow(dt)
    file_name <- paste0(pheno_name, "_1yr_followup_before_srs_within_ancestry_group_noNAs_", n, "_table_", date, ".csv")
    write.csv(dt, file_name, row.names = TRUE)
  })
}

# Main processing function
process_phenotype_data <- function(pheno_data, ancestry_samples_id_list, output_dir, date) {
  # Split by ancestry
  pheno_by_ancestry <- split_by_ancestry(pheno_data, ancestry_samples_id_list)
  
  # Split by ancestry and gender
  pheno_by_ancestry_and_gender <- split_by_ancestry_and_gender(pheno_by_ancestry)
  
  # Combine both lists
  data <- c(pheno_by_ancestry, pheno_by_ancestry_and_gender)
  
  # Save tables
  save_pheno_tables(data, output_dir, date)
}

# Usage
pheno <- gcta.pheno.scs.vol.roc  # Your phenotype data
process_phenotype_data(pheno, ancestry_samples_id_list, pheno_qc_dir, date)

#### OLD ####
# Ensure ancestry PCs have the same FID and IID as in phenotype data
pheno_data <- gcta.pheno.scs.vol.roc %>%
  left_join(ancestry_pcs, by = c("IID"))

# Ensure the covariate file includes all necessary covariates
# Prepare covariate data by selecting relevant columns, renaming them, removing samples with NA sex, and merging with ancestry principal components
covar_data <- smri.R5.1.baseline.y2.long %>%
  select(src_subject_id, rel_family_id, sex, interview_age) %>%
  rename(IID = src_subject_id, FID = rel_family_id, Age = interview_age) %>%
  filter(sex != 'NA') %>%
  left_join(ancestry_pcs, by = "IID")

# Convert sex to numeric and select family ID, individual ID, sex, age, and principal components
covar_data <- covar_data %>%
  mutate(Sex = ifelse(sex == "M", 1, 2)) %>%
  select(FID, IID, Sex, Age, starts_with("PC"))

# Save final covariate file
write.table(covar_data, 
            paste0(pheno_dir, "/gcta.pheno.scs.vol.roc.txt"), 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

# Ensure ancestry PCs have the same FID and IID as in phenotype data
covar_data_unique <- covar_data %>%
  distinct(IID, .keep_all = TRUE)

pheno_data_unique <- pheno_data %>%
  distinct(IID, .keep_all = TRUE)

# Perform the join to add Sex and Age from covar_data to pheno_data
gcta.pheno.scs.vol.roc.covar <- pheno_data_unique %>%
  left_join(covar_data_unique %>% select(IID, FID, Sex, Age), by = "IID")

# Ensure the FID column exists in the final data and always move it together with IID, Sex, and Age upfront
gcta.pheno.scs.vol.roc.covar <- gcta.pheno.scs.vol.roc.covar %>%
  mutate(FID = coalesce(FID.x, FID.y)) %>%  # Coalesce, in case of FID coming from different sources
  select(-FID.x, -FID.y) %>%  # Drop extra FID columns
  relocate(FID, IID, Sex, Age) 