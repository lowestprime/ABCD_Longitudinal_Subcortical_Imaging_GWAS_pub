library(dplyr);
library(tidyr);
library(ggplot2);
library(rlang);
library(lme4);
library(plotrix);
library(ggrepel);
library(gt);
library(ggridges);

source('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/release5_external_functions.R');

#### LOAD DATA AND RECODE TIMEPOINT ####
base_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/';
load('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/ABCD.Release.5.1.C4.Merged.Data.Rda');

# recode timepoint data
abcdData.R5.1$timepoint <- lapply(
  X = abcdData.R5.1$eventname, 
  FUN = recode.timepoint
  );

abcdData.R5.1$sex <- lapply(
  X = abcdData.R5.1$sex,
  FUN = recode.sex
  );
abcdData.R5.1$sex <- sapply(abcdData.R5.1$sex, paste, collapse = ",");

samples <- unique(abcdData.R5.1$src_subject_id);
id.timepoint.fields <- c('src_subject_id', 'eventname', 'timepoint');
print(colnames(abcdData.R5.1)); # variable names
print(length(unique(abcdData.R5.1$src_subject_id))); # number of samples
# extract C4 data
C4.expr <- abcdData.R5.1 %>%
  #filter(!is.na(c4a_expression) & is.finite(c4a_expression)) %>%
  dplyr::select(matches('src_subject_id|eventname|timepoint|c4'), -pc4);

#### PLOT C4 EXPRESSION ####
# plot histogram of C4A expression
pdf('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/plots/C4A_expression.pdf');
plot <- ggplot(data = C4.expr, mapping = aes(x = c4a_expression)) +
  geom_histogram(bins = 10, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of C4a Expression",
       x = "C4a Expression Level",
       y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))
print(plot);
dev.off();

# plot histogram of C4B expression
pdf('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/plots/C4B_expression.pdf');
plot <- ggplot(data = C4.expr, mapping = aes(x = c4b_expression)) +
  geom_histogram(bins = 10, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of C4b Expression",
       x = "C4b Expression Level",
       y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))
print(plot);
dev.off();

#### IMAGING ANALYSIS ####
# extract samples with smri data
smri.R5.1 <- abcdData.R5.1 %>% filter(imgincl_t1w_include == 1); # filter by mri QC
smri.R5.1 <- smri.R5.1 %>% dplyr::select(matches('smri|src_subject_id|eventname|timepoint|ethnicity|pc|ses|interview_age|sex|site_id_l|rel_family_id|mri_info_deviceserialnumber|mrisdp_453|pps'));
smri.R5.1 <- smri.R5.1[!apply(smri.R5.1[, !names(smri.R5.1) %in% id.timepoint.fields], 1, function(x) all(is.na(x))), ]; # remove NA rows

# create histograms of mean thickness/surface area/volume per timepoint
for (t in unique(smri.R5.1$eventname)) {
  smri.R5.1.timepoint <- smri.R5.1 %>% filter(eventname == t);
  print(nrow(smri.R5.1.timepoint));
  
  out_file_thick <- paste0(base_dir, 'plots/smri_thickness_', t, '.pdf');
  pdf(out_file_thick);
  plot <- ggplot(data = smri.R5.1.timepoint, mapping = aes(x = smri_thick_cdk_mean)) +
    geom_histogram() +
    ggtitle(paste0('Mean Cortical Thickness ', t, ' (n = ', nrow(smri.R5.1.timepoint), ')')) +
    theme_bw();
  print(plot);
  dev.off();
  
  out_file_surface_area <- paste0(base_dir, 'plots/smri_surfaceArea_', t, '.pdf');
  pdf(out_file_surface_area);
  plot <- ggplot(data = smri.R5.1.timepoint, mapping = aes(x = smri_area_cdk_total)) +
    geom_histogram() +
    ggtitle(paste0('Total surface area ', t, ' (n = ', nrow(smri.R5.1.timepoint), ')')) +
    theme_bw();
  print(plot);
  dev.off();
  
  out_file_volume <- paste0(base_dir, 'plots/smri_volume_', t, '.pdf');
  pdf(out_file_volume);
  plot <- ggplot(data = smri.R5.1.timepoint, mapping = aes(x = smri_vol_cdk_total)) +
    geom_histogram() +
    ggtitle(paste0('Total volume ', t, ' (n = ', nrow(smri.R5.1.timepoint), ')')) +
    theme_bw();
  print(plot);
  dev.off();
}

# average across brain hemispheres for each region
brain_region_phenotypes <- c("smri_thick_cdk_bankssts", "smri_thick_cdk_cdacate", "smri_thick_cdk_cdmdfr",
             "smri_thick_cdk_cuneus", "smri_thick_cdk_ehinal", "smri_thick_cdk_fusiform",
             "smri_thick_cdk_ifpl", "smri_thick_cdk_iftm", "smri_thick_cdk_ihcate", "smri_thick_cdk_locc",
             "smri_thick_cdk_lobfr", "smri_thick_cdk_lingual", "smri_thick_cdk_mobfr", "smri_thick_cdk_mdtm",
             "smri_thick_cdk_parahpal", "smri_thick_cdk_paracn", "smri_thick_cdk_parsopc", "smri_thick_cdk_parsobis",
             "smri_thick_cdk_parstgris", "smri_thick_cdk_pericc", "smri_thick_cdk_postcn", "smri_thick_cdk_ptcate",
             "smri_thick_cdk_precn", "smri_thick_cdk_pc", "smri_thick_cdk_rracate", "smri_thick_cdk_rrmdfr",
             "smri_thick_cdk_sufr", "smri_thick_cdk_supl", "smri_thick_cdk_sutm", "smri_thick_cdk_sm",
             "smri_thick_cdk_frpole", "smri_thick_cdk_tmpole", "smri_thick_cdk_trvtm", "smri_thick_cdk_insula",
             "smri_area_cdk_bankssts", "smri_area_cdk_cdacate", "smri_area_cdk_cdmdfr", "smri_area_cdk_cuneus",
             "smri_area_cdk_ehinal", "smri_area_cdk_fusiform", "smri_area_cdk_ifpl", "smri_area_cdk_iftm",
             "smri_area_cdk_ihcate", "smri_area_cdk_locc", "smri_area_cdk_lobfr", "smri_area_cdk_lingual",
             "smri_area_cdk_mobfr", "smri_area_cdk_mdtm", "smri_area_cdk_parahpal", "smri_area_cdk_paracn",
             "smri_area_cdk_parsopc", "smri_area_cdk_parsobis", "smri_area_cdk_parstgris", "smri_area_cdk_pericc",
             "smri_area_cdk_postcn", "smri_area_cdk_ptcate", "smri_area_cdk_precn", "smri_area_cdk_pc", "smri_area_cdk_rracate",
             "smri_area_cdk_rrmdfr", "smri_area_cdk_sufr", "smri_area_cdk_supl", "smri_area_cdk_sutm", "smri_area_cdk_sm",
             "smri_area_cdk_frpole", "smri_area_cdk_tmpole", "smri_area_cdk_trvtm", "smri_area_cdk_insula", "smri_vol_scs_cbwmatter",
             "smri_vol_scs_ltventricle", "smri_vol_scs_inflatvent", "smri_vol_scs_crbwmatter", "smri_vol_scs_crbcortex", "smri_vol_scs_tp", 
             "smri_vol_scs_caudate", "smri_vol_scs_putamen", "smri_vol_scs_pallidum", "smri_vol_scs_hpus", "smri_vol_scs_amygdala",
             "smri_vol_scs_vedc");
for (region.measure in brain_region_phenotypes){
  smri.R5.1[, region.measure] <- rowMeans(smri.R5.1[,c(paste0(region.measure, "rh"), paste0(region.measure, "lh"))], na.rm=TRUE);
}

# extract samples with smri data at all 3 timepoints
num.timepoints.per.sample <- table(smri.R5.1$src_subject_id);
complete.samples <- names(num.timepoints.per.sample[num.timepoints.per.sample == length(unique(smri.R5.1$eventname))]);
smri.R5.1.complete <- smri.R5.1 %>% filter(src_subject_id %in% complete.samples);

# merge with C4A data
smri.R5.1.complete.C4 <- merge(
  x = smri.R5.1.complete,
  y = C4.expr,
  by = c('src_subject_id', 'eventname', 'timepoint')
);

# calculate annual rate of change of cortical thickness, surface area, and volume
src_subject_id <- unlist(smri.R5.1.complete.C4$src_subject_id);
timepoint <- unlist(smri.R5.1.complete.C4$timepoint);
ordered_indices <- order(src_subject_id, timepoint);
smri.R5.1.complete.C4 <- smri.R5.1.complete.C4[ordered_indices, ];

smri.R5.1.complete.C4 <- smri.R5.1.complete.C4 %>%
  ungroup() %>%
  group_by(src_subject_id) %>%
  mutate(
    smri_thick_cdk_mean = as.numeric(smri_thick_cdk_mean),
    timepoint = as.numeric(timepoint),
    smri_thick_cdk_annualROC = ((smri_thick_cdk_mean - lag(smri_thick_cdk_mean, default = first(smri_thick_cdk_mean))) /
                                    lag(smri_thick_cdk_mean, default = first(smri_thick_cdk_mean))) / 
                                    (timepoint - lag(timepoint, default = first(timepoint))),
    smri_area_cdk_annualROC = ((smri_area_cdk_total - lag(smri_area_cdk_total, default = first(smri_area_cdk_total))) /
                             lag(smri_area_cdk_total, default = first(smri_area_cdk_total))) / 
                             (timepoint - lag(timepoint, default = first(timepoint))),
    smri_vol_cdk_annualROC = ((smri_vol_cdk_total - lag(smri_vol_cdk_total, default = first(smri_vol_cdk_total))) /
                       lag(smri_vol_cdk_total, default = first(smri_vol_cdk_total))) / 
                       (timepoint - lag(timepoint, default = first(timepoint)))
    );

################################ END OF DATA PREPROCESSING #####################
# create sample size tables
complete.samples.table <- smri.R5.1.complete.C4 %>%
  dplyr::select(src_subject_id, timepoint, ethnicity, sex, c4a_expression, c4_allele2_structure) %>%
  group_by(ethnicity, sex, c4_allele2_structure) %>%
  tally() %>%
  filter(c4_allele2_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS")) %>% # filter by allele structure
  filter(ethnicity %in% c('AFR', 'AMR', 'EUR')); # filter by ancestry

complete.samples.table |> gt(groupname_col = 'sex');

# create violin/boxplots for annual rate of change for cortical thickness, surface area, and volume
phenos <- c('smri_thick_cdk', 'smri_area_cdk', 'smri_vol_cdk');
for (phenotype in phenos) {
  var <- paste0(phenotype, '_annualROC');
  # violinplot of annual rate of change
  violin.ROC <- ggplot(smri.R5.1.complete.C4, aes(x = timepoint, y = !! sym(var), group = timepoint)) +
    geom_violin() +
    labs(x = "Timepoint", y = "Annual Rate of Change", title = var) +
    theme_minimal()
  print(violin.ROC);
  
  if (unlist(strsplit(phenotype, split = '_'))[2] == 'thick') {
    var <- paste0(phenotype, '_mean');
  }
  else var <- paste0(phenotype, '_total');

  # boxplot of phenotype over time
  boxplot.pheno <- ggplot(smri.R5.1.complete.C4, aes(x = timepoint, y = !! sym(var), group = timepoint)) +
    geom_boxplot() +
    labs(x = "Time (yrs)", y = var, title = paste0(var, ' over time (yrs)')) +
    theme_minimal()
  print(boxplot.pheno);
}

#### CORRELATIONS WITH C4A ####
# thickness
ggplot(smri.R5.1.complete.C4, aes(x = c4a_expression, y = smri_thick_cdk_annualROC)) +
  geom_jitter(width = 0.5);
cor(
  x = smri.R5.1.complete.C4$smri_thick_cdk_annualROC,
  y = smri.R5.1.complete.C4$c4a_expression,
  use = 'pairwise.complete.obs'
  );

# surface area
ggplot(smri.R5.1.complete.C4, aes(x = c4a_expression, y = smri_area_cdk_annualROC)) +
  geom_jitter(width = 0.5);
cor(
  x = smri.R5.1.complete.C4$smri_area_cdk_annualROC,
  y = smri.R5.1.complete.C4$c4a_expression,
  use = 'pairwise.complete.obs'
  );

# volume
ggplot(smri.R5.1.complete.C4, aes(x = c4a_expression, y = smri_vol_cdk_annualROC)) +
  geom_jitter(width = 0.5);
cor(
  x = smri.R5.1.complete.C4$smri_vol_cdk_annualROC,
  y = smri.R5.1.complete.C4$c4a_expression,
  use = 'pairwise.complete.obs'
  );

#### 04_PheWAS.Behavior.1YearFollowup.Rmd ####
# all data
abcdData_alltimepoints_C4 <- smri.R5.1.complete.C4;
abcdData_alltimepoints_C4_filtered <- smri.R5.1.complete.C4 %>% filter(ethnicity %in% c("AFR", "AMR", "EUR"));

abcdData_alltimepoints_C4_filtered <- abcdData_alltimepoints_C4_filtered %>%
  filter(c4_allele1_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))
abcdData_alltimepoints_C4_filtered <- abcdData_alltimepoints_C4_filtered %>%
  filter(c4_allele2_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))

abcdData_alltimepoints_C4_filtered %>% group_by(ethnicity) %>% tally();

# baseline
abcdData_baseline <- smri.R5.1.complete.C4 %>% filter(eventname == "baseline_year_1_arm_1");
abcdData_baseline_C4 <- abcdData_baseline %>% filter(!is.na(c4a_total))

abcdData_baseline_C4_filtered <- abcdData_baseline_C4 %>% filter(ethnicity %in% c("AFR", "AMR", "EUR"))

abcdData_baseline_C4_filtered <- abcdData_baseline_C4_filtered %>%
  filter(c4_allele1_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))
abcdData_baseline_C4_filtered <- abcdData_baseline_C4_filtered %>%
  filter(c4_allele2_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))

abcdData_baseline_C4_filtered %>% group_by(ethnicity) %>% tally()

# y2
abcdData_Y2 <- smri.R5.1.complete.C4 %>% filter(eventname == "2_year_follow_up_y_arm_1");
abcdData_Y2_C4 <- abcdData_Y2 %>% filter(!is.na(c4a_total))

abcdData_Y2_C4_filtered <- abcdData_Y2_C4 %>% filter(ethnicity %in% c("AFR", "AMR", "EUR"))

abcdData_Y2_C4_filtered <- abcdData_Y2_C4_filtered %>%
  filter(c4_allele1_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))
abcdData_Y2_C4_filtered <- abcdData_Y2_C4_filtered %>%
  filter(c4_allele2_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))

abcdData_Y2_C4_filtered %>% group_by(ethnicity) %>% tally() 

# y4
abcdData_Y4 <- smri.R5.1.complete.C4 %>% filter(eventname == "4_year_follow_up_y_arm_1");
abcdData_Y4_C4 <- abcdData_Y4 %>% filter(!is.na(c4a_total))

abcdData_Y4_C4_filtered <- abcdData_Y4_C4 %>% filter(ethnicity %in% c("AFR", "AMR", "EUR"))

abcdData_Y4_C4_filtered <- abcdData_Y4_C4_filtered %>%
  filter(c4_allele1_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))
abcdData_Y4_C4_filtered <- abcdData_Y4_C4_filtered %>%
  filter(c4_allele2_structure %in% c("AL","AL-AL","AL-BL", "AL-BS", "BS"))

abcdData_Y4_C4_filtered %>% group_by(ethnicity) %>% tally()

#### 06_GREx.Brain.Rmd ####

# SPLIT BY ANCESTRY
abcdData_alltimepoints_C4_filtered.split_by_ancestry <- list();
ancestries <- c('AFR', 'AMR', 'EUR');
for (ancestry in ancestries) {
  alltimepoints.ancestry <- filter(abcdData_alltimepoints_C4_filtered, ethnicity == ancestry);
  abcdData_alltimepoints_C4_filtered.split_by_ancestry[[ancestry]] <- alltimepoints.ancestry;
}

# SPLIT BY SEX
abcdData_alltimepoints_C4_filtered.split_by_sex <- list();
sexes <- c('Female', 'Male');
for (sex in sexes) {
  sex.abbr <- substr(sex, 1, 1);
  print(sex.abbr);
  alltimepoints.sex <- filter(abcdData_alltimepoints_C4_filtered, sex == sex.abbr);
  abcdData_alltimepoints_C4_filtered.split_by_sex[[sex]] <- alltimepoints.sex;
}

#### BRAIN PHENOTYPE ASSOCIATIONS ####

# Association by ancestry
datalist <- list(
  'ALL' = abcdData_alltimepoints_C4_filtered,
  'AFR' = abcdData_alltimepoints_C4_filtered.split_by_ancestry['AFR']$AFR,
  'AMR' = abcdData_alltimepoints_C4_filtered.split_by_ancestry['AMR']$AMR,
  'EUR' = abcdData_alltimepoints_C4_filtered.split_by_ancestry['EUR']$EUR
  );

# Association by sex
datalist <- list(
  'ALL' = abcdData_alltimepoints_C4_filtered,
  'Female' = abcdData_alltimepoints_C4_filtered.split_by_sex['Female']$Female,
  'Male' = abcdData_alltimepoints_C4_filtered.split_by_sex['Male']$Male
  );

# Association by timepoint
datalist <- list(
  'BASELINE' = abcdData_baseline_C4_filtered,
  'Y2' = abcdData_Y2_C4_filtered,
  'Y4' = abcdData_Y4_C4_filtered
  );

phenos <- c(brain_region_phenotypes,
            'smri_thick_cdk_mean', 'smri_area_cdk_total', 'smri_vol_cdk_total',
            'smri_thick_cdk_annualROC', 'smri_area_cdk_annualROC', 'smri_vol_cdk_annualROC'
            );
for (j in names(datalist)) {
  df <- data.frame();
  print(j);
  for (i in phenos) {
    print(i);
    pheno_idx <- which(colnames(abcdData_alltimepoints_C4_filtered) == i);
    this_var = datalist[[j]][, pheno_idx];
    if (j %in% c('Female', 'Male')) {
      covariates <- paste0('ses_AvgIncomeParentEdu + pc1 + pc2 + pc3 + pc4 + (1|mri_info_deviceserialnumber) + (1|rel_family_id)');
    }
    else {
      covariates <- paste0('ses_AvgIncomeParentEdu + pc1 + pc2 + pc3 + pc4 + sex + (1|mri_info_deviceserialnumber) + (1|rel_family_id)');
    }
    mod0 <- lmer(
      formula = as.formula(paste0('as.matrix(this_var) ~ c4b_expression +', covariates)),
      data = datalist[[j]], REML = FALSE);
    mod1 <- lmer(
      formula = as.formula(paste0('as.matrix(this_var) ~ c4a_expression +', covariates)),
      data = datalist[[j]], REML = FALSE);
    mod2 <- lmer(
      formula = as.formula(paste0('as.matrix(this_var) ~ c4a_expression + c4b_expression +', covariates)),
      data = datalist[[j]], REML = FALSE);
    p1 <- anova(mod0, mod2);
    p2 <- anova(mod1, mod2);
    tmp_df <- data.frame(
      data = j,
      ID = i, 
      N = nobs(mod2),
      beta_C4A = summary(mod2)$coefficients[2, 1],
      beta_C4B = summary(mod2)$coefficients[3, 1],
      stderr_C4A = summary(mod2)$coefficients[2, 2],
      stderr_C4B = summary(mod2)$coefficients[3, 2],
      P_C4A = p1[2, 8],
      P_C4B = p2[2, 8]
      );
    if (nrow(tmp_df) != 0) {
      df <- rbind(df, tmp_df);
    }
  }
  df$fdr_C4A = p.adjust(df$P_C4A,'fdr');
  df$fdr_C4B = p.adjust(df$P_C4B,'fdr');
  pheno_category <- character(length(df$ID));
  for (i in seq_along(df$ID)) {
    if (grepl("vol", df$ID[i])) {
      pheno_category[i] <- "Volume"
    } else if (grepl("thick", df$ID[i])) {
      pheno_category[i] <- "Thickness"
    } else if (grepl("area", df$ID[i])) {
      pheno_category[i] <- "Area"
    } else {
      pheno_category[i] <- "NA"
    }
  }
  df$Category <- pheno_category;
  df$Direction <- ifelse(df$beta_C4A >= 0, "Positive", "Negative");
  print(df);
  pdf(file = paste0(base_dir, 'plots/C4Xsmri_phenotype_regression_', j, '.pdf'), width = 12, height = 8);
  plot <- plot.pheWAS(df);
  print(plot + ylim(0,10));
  dev.off();
  # df <- gt(df);
  # print(df);
}

# Sex differences in significant measures

#### BRAIN PHENOTYPE RATE OF CHANGE ASSOCIATIONS ####
phenos <- c('smri_thick_cdk_annualROC', 'smri_area_cdk_annualROC', 'smri_vol_cdk_annualROC');

all.ROC <- filter(smri.R5.1.complete.annualROC.C4A, eventname != 'baseline_year_1_arm_1');

afr.ROC <- filter(all.ROC, ethnicity == "AFR");
latinx.ROC <- filter(all.ROC, ethnicity == "AMR");
eur.ROC <- filter(all.ROC, ethnicity == "EUR");

male.ROC <- filter(all.ROC, sex == "M");
female.ROC <- filter(all.ROC, sex == "F");

# Association of annual rate of change by ancestry
datalist <- list(all.ROC, afr.ROC, latinx.ROC, eur.ROC);
names(datalist) <- c("ALL", "AFR", "AMR", "EUR");

# Association of annual rate of change by sex
datalist <- list(all.ROC, male.ROC, female.ROC);
names(datalist) <- c("ALL", "MALE", "FEMALE");

# Association of annual rate of change by timepoint
datalist <- list(all.ROC, abcdData_Y2_C4_filtered, abcdData_Y4_C4_filtered);
names(datalist) <- c('ALL', 'Y2', 'Y4');

for (j in 1:length(datalist)) {
  df <- data.frame() 
  print(names(datalist)[j])
  for (i in phenos) {
    pheno_idx <- which(colnames(abcdData_baseline_C4_filtered) == i);
    this_var = datalist[[j]][, pheno_idx];
    mod0 <- lmer(as.matrix(this_var) ~ c4b_expression + 
                   interview_age + ses_AvgIncomeParentEdu + pc1 + pc2 + pc3 + pc4 + 
                   + (1 | mri_info_deviceserialnumber) + (1|rel_family_id),
                 data = datalist[[j]], REML = FALSE)
    mod1 <- lmer(as.matrix(this_var) ~ c4a_expression + 
                   interview_age + ses_AvgIncomeParentEdu + pc1 + pc2 + pc3 + pc4 + 
                   + (1 | mri_info_deviceserialnumber) + (1|rel_family_id),
                 data = datalist[[j]], REML = FALSE)
    mod2 <- lmer(as.matrix(this_var) ~ c4a_expression + c4b_expression + 
                   interview_age + ses_AvgIncomeParentEdu + pc1 + pc2 + pc3 + pc4 + 
                   + (1 | mri_info_deviceserialnumber) + (1|rel_family_id),
                 data = datalist[[j]], REML = FALSE)
    p1 <- anova(mod0, mod2)
    p2 <- anova(mod1, mod2)
    df <- rbind(df, data.frame(data = names(datalist)[j],
                               ID = i, 
                               N = nobs(mod2),
                               beta_C4A = summary(mod2)$coefficients[2, 1],
                               beta_C4B = summary(mod2)$coefficients[3, 1],
                               stderr_C4A = summary(mod2)$coefficients[2, 2],
                               stderr_C4B = summary(mod2)$coefficients[3, 2],
                               P_C4A = p1[2, 8],
                               P_C4B = p2[2, 8]))
  }
  df$fdr_C4A = p.adjust(df$P_C4A,'fdr')
  df$fdr_C4B = p.adjust(df$P_C4B,'fdr')
  print(df);
  plot <- plot.pheWAS(df);
  print(plot + ylim(0,5));
}


#### CLOSER EXPLORATION OF BRAIN PHENOTYPES SIGNIFICANTLY ASSOC WITH C4A IN REGRESSIONS ####

# baseline
abcdData_baseline_C4_filtered_brain <- 
  dplyr::select(abcdData_baseline_C4_filtered, src_subject_id, eventname, timepoint, sex, ethnicity, signif.regions$ID)

abcdData_Y2_C4_filtered_brain <- 
  dplyr::select(abcdData_Y2_C4_filtered, src_subject_id, eventname, timepoint, sex, ethnicity, signif.regions$ID)

abcdData_Y4_C4_filtered_brain <- 
  dplyr::select(abcdData_Y4_C4_filtered, src_subject_id, eventname, timepoint, sex, ethnicity, signif.regions$ID)

#### Longtitudinal entorhinal cortex vs psychosis symptoms ####
library(plotrix);
library(rstatix);

## Year 2 follow-up
for (phenotype in signif.regions$ID){
  print(phenotype);
  df.year2 <- create.severity.number.dataframe(timepoint.data = abcdData_Y2_C4_filtered, pheno = phenotype);
  df.year4 <- create.severity.number.dataframe(timepoint.data = abcdData_Y4_C4_filtered, pheno = phenotype);
  
  df.year2$year <- "2-Year Follow-Up"
  df.year4$year <- "4-Year Follow-Up"
  df.year2.year4 <- rbind(df.year2, df.year4)
  df.number.of.events <- df.year2.year4[c(1:4, 9:12),]
  df.severity.of.events <- df.year2.year4[c(5:8, 13:16),]
  
  # Plot 
  my.plot <- ggplot(df.number.of.events, aes(x=Size, y=mean, fill=Size)) + 
    scale_y_continuous(breaks = c(0, 1, 2)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) 
  my.plot1 <- my.plot + theme_classic() + labs(x=phenotype, y = "Psychosis Symptom Number", fill = "SA Quartile") +
    scale_fill_manual(values=c('#999999','white')) + facet_grid(Sex~year) +
    geom_text(aes(label = format(mean, digits = 2), y = 0.7)) + theme(legend.position="left")
  pdf(file = paste0('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis', phenotype, '_x_psychosisNumEvents.longit.pdf'));
  print(my.plot1)
  dev.off();

  my.plot <- ggplot(df.severity.of.events, aes(x=Size, y=mean, fill=Size)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))
  my.plot1 <- my.plot + theme_classic() + labs(x=phenotype, y = "Psychosis Symptoms Severity", fill = "SA Quartile") +
    scale_fill_manual(values=c('#999999','white')) + facet_grid(Sex~year) +
    geom_text(aes(label = format(mean, digits = 2), y = 0.7)) + theme(legend.position="left")
  my.plot1
  pdf(file = paste0('/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/individual_brain_regions/', phenotype, '_x_psychosisSeverity.longit.pdf'));
  print(my.plot1)
  dev.off();
}

