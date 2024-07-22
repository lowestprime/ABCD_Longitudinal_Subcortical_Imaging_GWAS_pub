#### Useful snippets ####

# We can list the libraries that are actually loaded doing
(.packages())

# Unload all currently loaded packages using pacman
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Check which packages are currently loaded
pacman::p_loaded()

# count and list the number of each unique in dataframe column
summary(factor(merged_data_final$mri_info_deviceserialnumber))
summary(factor(merged_data_final$batch))
summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))
combined_summary_df <- as.data.frame(as.table(summary(factor(c(merged_data_final$mri_info_deviceserialnumber, merged_data_final$batch)))))

# take me back home
setwd('~/')

# take me back to proj dir
setwd('~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis')

# list column names
paste(names(merged_data_final), collapse = ", ")

# List all subdirectories
list.dirs(path = '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data', full.names = T, recursive = T)

# sum different values in df
sum(merged_data_no_na$ethnicity %in% c("AFR", "AMR", "EUR"))

# get col names within a dataframe matching a specific naming scheme/pattern using grepl
colnames(smri.R5.1.all)[grepl("^smri_vol_scs", colnames(smri.R5.1.all))]

# get summary of specific cols basic data characteristics within a dataframe that match a naming scheme/pattern using grepl
summary(merged_data_normalized)[grepl("^smri_vol_scs", colnames(merged_data_normalized))]

# run rankTransPheno
merged_data_normalized <- merged_data_final %>%
  mutate(across(starts_with("smri_vol_"), ~ rankTransPheno(.x, 0.5)))

# launch gptstudio
gptstudio:::gptstudio_chat()

# Restart RStudio session
.rs.restartR()

#### data prep ####

# Check L/R hemisphere averages QC
grep(pattern = "smri_vol_scs_", x = colnames(smri.R5.1), value = T);
grep(pattern = "smri_vol_scs_.*(lh$|rh$|l$|r$)", x = colnames(smri.R5.1), value = T);
# Manually check the first few rows
head(smri.R5.1[, c("smri_vol_scs_amygdalalh", "smri_vol_scs_amygdalarh", "smri_vol_scs_amygdala")])
# Check correlations
cor(smri.R5.1$smri_vol_scs_amygdala, (smri.R5.1$smri_vol_scs_amygdalalh + smri.R5.1$smri_vol_scs_amygdalarh) / 2)
# Plotting to visually inspect the correctness of averages
ggplot(smri.R5.1, aes(x = smri_vol_scs_amygdalalh, y = smri_vol_scs_amygdalarh)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_point(aes(x = (smri_vol_scs_amygdalalh + smri_vol_scs_amygdalarh)/2, y = smri_vol_scs_amygdala), color = "blue")

# # Jack ROC code by region
# rate_of_change_per_brain_region <- data.frame(
#   'src_subject_id' = character(0),
#   'rate_of_change_timepoints' = character(0)
# )
# for (pheno in brain_region_phenotypes) {
#   rate_of_change_per_brain_region[[pheno]] <- numeric(0);
# }
# for (sample in unique(smri.R5.1.complete.C4$src_subject_id)) {
#   print(sample);
#   sample_data <- smri.R5.1.complete.C4 %>% filter(src_subject_id == sample);
#   sample_roc <- data.frame(
#     'sample' = rep(sample, 3),
#     'rate_of_change_timepoints' = c('Baseline_T4', 'Baseline_T2', 'T2_T4')
#   );
#   for (pheno in brain_region_phenotypes) {
#     Baseline_T4.pheno = ((sample_data[[pheno]][3] - sample_data[[pheno]][1]) / sample_data[[pheno]][1]) / (sample_data$timepoint[3] - sample_data$timepoint[1])
#     Baseline_T2.pheno = ((sample_data[[pheno]][2] - sample_data[[pheno]][1]) / sample_data[[pheno]][1]) / (sample_data$timepoint[2] - sample_data$timepoint[1])
#     T2_T4.pheno       = ((sample_data[[pheno]][3] - sample_data[[pheno]][2]) / sample_data[[pheno]][2]) / (sample_data$timepoint[3] - sample_data$timepoint[2])
#     pheno_roc <- data.frame(
#       'sample' = rep(sample, 3),
#       'rate_of_change_timepoints' = c('Baseline_T4', 'Baseline_T2', 'T2_T4'),
#       'pheno' = c(Baseline_T4.pheno, Baseline_T2.pheno, T2_T4.pheno)
#     );
#     names(pheno_roc)[3] <- paste0('roc_', pheno);
#     #print(pheno_roc);
#     sample_roc <- merge(
#       x = sample_roc,
#       y = pheno_roc,
#       by = c('sample', 'rate_of_change_timepoints')
#     );
#   }
#   rate_of_change_per_brain_region <- rbind(rate_of_change_per_brain_region, sample_roc);
# }
# 
# # Cooper ROC code by region
# roi_volumes <- grep("smri_vol_scs_", names(smri.R5.1.complete), value = TRUE)
# smri.R5.1.complete.baseline.T2 <- smri.R5.1.complete %>%
#   group_by(src_subject_id) %>%
#   mutate(across(all_of(roi_volumes), ~ (.x - lag(.x)) / (timepoint - lag(timepoint)), .names = "roc_{.col}"))

# Normality testing 
# qqplot by region
qqnorm(smri.R5.1.baseline.y2.ROC$smri_vol_scs_3rdventricle, main = "Q-Q plot for 3rd Ventricle Volume")
qqline(smri.R5.1.baseline.y2.ROC$smri_vol_scs_3rdventricle, col = "steelblue", lwd = 2)

# shapiro walk test random sample of 5000 from larger ROC dataset
set.seed(123)  # for reproducibility
sample_data <- sample(smri.R5.1.baseline.y2.ROC$smri_vol_scs_3rdventricle, 5000)
shapiro.test(sample_data)
# saveRDS
saveRDS(object = smri.R5.1.baseline.y2.ROC.long, file = "smri.R5.1.baseline.y2.ROC.long.rds")
# comparison_result
comparison_result <- all.equal(smri.R5.1, smri.R5.1t)
print(comparison_result)

#### read in example mlma txt input files ####

# Define file paths
# phenotype_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/Phen_AllSubj_abcd_pssrs01.txt_1yr_followup_ssrs_42_p_within_ancestry_group_noNAs_11023_w_fid.txt"
covariate_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/covar_AllSubj_batch_gender_noNAs_baseline_11665.txt"
# qcovar_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/qcovar_ABCD5_EUR_PC20_SOR_related_bigsnpr.txt"
# batchinfodir <- "/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5/genomics_sample03/smokescreen"
# batchinfodir_file <- "/u/project/lhernand/FROM_GANDALM/ABCD_Release_4/genomics_sample03/ABCD_release3.0_.batch_info.txt"

# # Modify IDs with ABCD ID format
# setwd(batchinfodir)
# Batch <- read.delim("ABCD_202209.updated.nodups.curated.batch.info",  header = TRUE, na.strings = c("", "NA"), sep = ",")
# #Batch<- rename(Batch, subjectkey = IID)
# Batch$IID2 <- Batch$IID

# Read the files without headers
# phenotype_data <- read.table(phenotype_file, header = FALSE)
covariate_data <- read.table(covariate_file, header = FALSE)
# qcovar_data <- read.table(qcovar_file, header = FALSE, sep = "\t")
# batchinfodir_data <- read.table(batchinfodir_file, header = TRUE, sep = "\t")

# Assign appropriate column names for phenotype_data and covariate_data
# colnames(phenotype_data) <- c("FID", "IID", "phenotype")
colnames(covariate_data) <- c("FID", "IID", "batch", "sex")

# Assign column names to the qcovar data
# num_pcs <- ncol(qcovar_data) - 3 # Subtracting 3 for FID, IID, and interview_age_years_nodecimal
# pc_names <- paste0("PC", 1:num_pcs)
# colnames(qcovar_data) <- c("FID", "IID", pc_names, "interview_age_years_nodecimal")
# 
# Ensure ancestry PCs have the same FID and IID as in phenotype data
# covariate_data_unique <- covariate_data %>%
#   distinct(IID, .keep_all = TRUE)
# 
# gcta.pheno.scs.vol.roc_unique <- gcta.pheno.scs.vol.roc %>%
#   distinct(IID, .keep_all = TRUE)
# 
# # Perform the join to add Sex and Age from covar_data to pheno_data
# gcta.pheno.scs.vol.roc.covar <- pheno_data_unique %>%
#   left_join(covar_data_unique %>% select(IID, FID, Sex, Age), by = "IID")
# 
# # Ensure the FID column exists in the final data and always move it together with IID, Sex, and Age upfront
# gcta.pheno.scs.vol.roc.covar <- gcta.pheno.scs.vol.roc.covar %>%
#   mutate(FID = coalesce(FID.x, FID.y)) %>%  # Coalesce, in case of FID coming from different sources
#   select(-FID.x, -FID.y) %>%  # Drop extra FID columns
#   relocate(FID, IID, Sex, Age) 

# Step 2: Merge gcta.pheno.scs.vol.roc with covariate_data genotyping batch and ancestry_pcs 1-20
# Filter covariate_data and ancestry_pcs to include only the IIDs present in gcta.pheno.scs.vol.roc
filtered_covariate_data <- covariate_data %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")
filtered_ancestry_pcs <- ancestry_pcs %>% semi_join(gcta.pheno.scs.vol.roc, by = "IID")

# Merge gcta.pheno.scs.vol.roc with filtered_covariate_data and filtered_ancestry_pcs
merged_data <- gcta.pheno.scs.vol.roc %>%
  left_join(filtered_covariate_data %>% select(IID, batch), by = "IID") %>%
  left_join(filtered_ancestry_pcs %>% select(IID, starts_with("PC")), by = "IID")

# drop any rows from merged_data with NAs (specifically known to be within the batch and PCs columns)
# reorder columns
merged_data_no_na <- merged_data %>%
  select(FID, IID, everything()) %>% # Ensures FID and IID are the first two columns
  drop_na()

# replace FID = rel_family_id, IID = src_subject_id
merged_data_no_na <- merged_data_no_na %>%
  mutate(FID = IID)

# Save final_data as a space-separated text file suitable for GCTA MLMA
setwd(pheno_dir)
write.table(merged_data_no_na, file = "gcta_mlma_SCS_ROC_master.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Filtering rows with any `NA` values
rows_with_na <- merged_data %>% 
  filter(if_any(everything(), is.na)) %>%
  select(FID, IID, batch, PC1:PC20)

# Filter out rows in `rows_with_na` where `IID` is not present in either `ancestry_pcs$IID` or `covariate_data$IID`
# This results in a dataset with `IID`s that are missing in either data frame
all_IIDs_missing <- rows_with_na %>% 
  filter(!(IID %in% ancestry_pcs$IID) | !(IID %in% covariate_data$IID))

# Extract unique IIDs and save to a text file
unique_IIDs <- all_IIDs_missing %>% 
  pull(IID) %>%
  unique()

# Write IIDs to a text file
write_lines(unique_IIDs, "missing_IIDs.txt")

# Filter `smri.R5.1.baseline.y2` to retain rows with `src_subject_id` present in `all_IIDs_missing$IID`
# Remove duplicate rows based on `src_subject_id` to ensure each subject ID is unique in the final dataset
smri.R5.1.baseline.y2_missing_baseline_pcs <- smri.R5.1.baseline.y2 %>%
  filter(src_subject_id %in% all_IIDs_missing$IID) %>%
  distinct(src_subject_id, .keep_all = TRUE)

# Save merged_data_no_na as a space-separated text file
write.table(smri.R5.1.baseline.y2_missing_baseline_pcs, file = "missing_IIDs_info.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Step 3: Create the phenotype split .txt files
# Phenotype files
# Assuming gcta.pheno.scs.vol.roc is your data frame and pheno_dir is the directory where you want to save the files
create_phenotype_files(merged_data_no_na, pheno_dir)

# Covariate file
covariate_data <- final_data %>%
  select(src_subject_id, sex, batch, mri_info_deviceserialnumber) %>%
  rename(IID = src_subject_id)
write.table(covariate_data, file = "covariate.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Quantitative Covariate file


#### ROC CALC MOD ####
# Helper function to get latest timepoint information
get_latest_timepoint_info <- function(df, timepoint) {
  df %>%
    filter(timepoint == !!timepoint & sex != 'NA') %>%
    select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age)
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

#### PIVOT MOD ####
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
    select(src_subject_id, rel_family_id, sex, mri_info_deviceserialnumber, interview_age, roi, everything())
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

# filtered phenotypes
roi_columns <- smri.R5.1.baseline.y2.ROC.filtered %>% 
  select(starts_with("smri_vol_")) %>% 
  colnames()

