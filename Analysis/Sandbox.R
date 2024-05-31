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
