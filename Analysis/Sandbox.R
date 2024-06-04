#### data prep ####

# # Extract unique c(src_subject_id, timepoint, ethnicity, sex)fix
# unique_regions <- names(smri.R5.1) %>%
#   str_extract("^smri_vol_scs_.*(?=lh$|rh$|l$|r$)") %>%
#   na.omit() %>%
#   unique();

# # Average across brain hemispheres for each applicable scs region (non-vectorized approach)
# # Get base scs region names omitting lh and rh or l and r suffixes
# unique_regions <- names(smri.R5.1) %>%
#   str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>%
#   na.omit() %>%
#   unique();
# # Define the pairs of columns for which to calculate the average
# column_pairs <- setNames(lapply(unique_regions, function(x) 
#   c(paste0(x, ifelse(x == "smri_vol_scs_aa", "l", "lh")),
#     paste0(x, ifelse(x == "smri_vol_scs_aa", "r", "rh")))), unique_regions);
# # Calculate the average for each pair of columns
# for (region in unique_regions) {
#   if (all(column_pairs[[region]] %in% names(smri.R5.1))) {
#     smri.R5.1[, region] <- rowMeans(smri.R5.1[, column_pairs[[region]]], na.rm = TRUE)
#   }
# }
# 
# # Average across brain hemispheres for each relevant scs region (vectorized approach)
# # Get base scs region names omitting lh and rh or l and r suffixes
# unique_regions <- names(smri.R5.1) %>%
#   str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>%
#   na.omit() %>%
#   unique();
# # Define the pairs of columns for which to calculate the average
# column_pairs <- setNames(lapply(unique_regions, function(x)
#   c(paste0(x, ifelse(x == "smri_vol_scs_aa", "l", "lh")),
#     paste0(x, ifelse(x == "smri_vol_scs_aa", "r", "rh")))), unique_regions);
# # Calculate the average for each pair of columns
# averages <- sapply(unique_regions, function(region) {
#   if (all(column_pairs[[region]] %in% names(smri.R5.1))) {
#     rowMeans(smri.R5.1[, column_pairs[[region]]], na.rm = TRUE)
#   } else {
#     rep(NA_real_, nrow(smri.R5.1))
#   }
# });
# smri.R5.1 <- cbind(smri.R5.1, averages)
# 
# # COMPACT VECTORIZED: Extracting unique regions, forming column pairs, calculating averages, and binding to original data
# smri.R5.1 <- cbind(smri.R5.1, sapply(setNames(lapply(names(smri.R5.1) %>% 
#                                                        str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>% 
#                                                        na.omit() %>% 
#                                                        unique(), function(x) c(paste0(x, ifelse(x == "smri_vol_scs_aa", "l", "lh")), paste0(x, ifelse(x == "smri_vol_scs_aa", "r", "rh")))), 
#                                               names(smri.R5.1) %>% str_extract("^smri_vol_scs_.*(?=(lh|rh|l|r)$)") %>% na.omit() %>% unique()), function(region) {
#                                                 if (all(region %in% names(smri.R5.1))) rowMeans(smri.R5.1[, region, drop = FALSE], na.rm = TRUE) else rep(NA_real_, nrow(smri.R5.1))
#                                               })) %>%
#   select(-matches("^smri_vol_scs_.*(lh|rh|aal|aar)$"));

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
phenotype_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/Phen_AllSubj_abcd_pssrs01.txt_1yr_followup_ssrs_42_p_within_ancestry_group_noNAs_11023_w_fid.txt"
covariate_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/covar_AllSubj_batch_gender_noNAs_baseline_11665.txt"
qcovar_file <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/qcovar_ABCD5_EUR_PC20_SOR_related_bigsnpr.txt"

# Read the files without headers
phenotype_data <- read.table(phenotype_file, header = FALSE, sep = " ")
covariate_data <- read.table(covariate_file, header = FALSE, sep = " ")
qcovar_data <- read.table(qcovar_file, header = FALSE, sep = "\t")

# Assign appropriate column names for phenotype_data and covariate_data
colnames(phenotype_data) <- c("FID", "IID", "phenotype")
colnames(covariate_data) <- c("FID", "IID", "batch", "sex")

# Assign column names to the qcovar data
num_pcs <- ncol(qcovar_data) - 3 # Subtracting 3 for FID, IID, and interview_age_years_nodecimal
pc_names <- paste0("PC", 1:num_pcs)
colnames(qcovar_data) <- c("FID", "IID", pc_names, "interview_age_years_nodecimal")

# Count the occurrences of each unique value in mri_info_deviceserialnumber
device_counts <- smri.R5.1 %>%
  count(mri_info_deviceserialnumber)

# device_counts <- summary(factor(smri.R5.1$mri_info_deviceserialnumber))
# Print the counts for device serial numbers
print(device_counts)

# Count the occurrences of each unique value in mrisdp_453
protocol_counts <- smri.R5.1 %>%
  count(mrisdp_453)

# Print the counts for MRI scanning protocols
print(protocol_counts)

# Step 1: Filter smri.R5.1.baseline.y2.ROC to remove rows with NA in the sex column
filtered_data <- smri.R5.1.baseline.y2.ROC %>%
  filter(!is.na(sex))

# Step 2: Merge with covariate_data and ancestry_pcs
# Assuming that IID in ancestry_pcs corresponds to src_subject_id in smri.R5.1.baseline.y2.ROC
merged_data <- filtered_data %>%
  left_join(covariate_data, by = c("src_subject_id" = "IID")) %>%
  left_join(ancestry_pcs, by = c("src_subject_id" = "IID"))

# Ensure the final dataset is consistent
final_data <- merged_data %>%
  select(src_subject_id, rel_family_id, all_of(roc_phenotypes), sex, batch, mri_info_deviceserialnumber, interview_age, starts_with("PC"))

# Step 3: Create the .txt files
# Phenotype files
for (phenotype in roc_phenotypes) {
  phenotype_data <- final_data %>%
    select(src_subject_id, !!sym(phenotype)) %>%
    rename(IID = src_subject_id)
  write.table(phenotype_data, file = paste0(phenotype, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Covariate file
covariate_data <- final_data %>%
  select(src_subject_id, sex, batch, mri_info_deviceserialnumber) %>%
  rename(IID = src_subject_id)
write.table(covariate_data, file = "covariate.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Quantitative Covariate file
qcovar_data <- final_data %>%
  select(src_subject_id, interview_age, starts_with("PC")) %>%
  rename(IID = src_subject_id)
write.table(qcovar_data, file = "qcovar.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

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

roc_columns <- smri.R5.1.baseline.y2.ROC.long %>% 
  select(starts_with("smri_vol_")) %>% 
  colnames()

