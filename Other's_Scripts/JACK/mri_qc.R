library(dplyr);
library(VennDiagram);

base_dir <- '/u/project/lhernand/';
ABCD_dir <- paste0(base_dir, 'shared/GenomicDatasets/ABCD_Release_5.1/core/');
out_dir <- paste0(base_dir, 'cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/');

mri_qc_incl <- read.delim(
  file = paste0(ABCD_dir, 'imaging/mri_y_qc_incl.csv'),
  header = TRUE,
  sep = ','
  );
# print(head(mri_qc_incl));
# print(unique(mri_qc_incl$eventname));

num_samples <- length(unique(mri_qc_incl$src_subject_id));
samples_baseline <- mri_qc_incl[mri_qc_incl$eventname == 'baseline_year_1_arm_1',];
samples_2yr <- mri_qc_incl[mri_qc_incl$eventname == '2_year_follow_up_y_arm_1',];
samples_4yr <- mri_qc_incl[mri_qc_incl$eventname == '4_year_follow_up_y_arm_1',];
print(paste0('Number of unique samples: ', num_samples));
print(paste0('Number of subjects with baseline: ', nrow(samples_baseline)));
print(paste0('Number of subjects with 2yr follow up: ', nrow(samples_2yr)));
print(paste0('Number of subjects with 4yr follow up: ', nrow(samples_4yr)));

# mri_scans_all <- c('imgincl_t1w_include', 'imgincl_t2w_include',
#              'imgincl_dmri_include', 'imgincl_rsfmri_include',
#              'imgincl_mid_include', 'imgincl_nback_include',
#              'imgincl_sst_include');
mri_scans_structural <- c('imgincl_t1w_include', 'imgincl_dmri_include', 'imgincl_rsfmri_include');
samples_baseline_allpass <- samples_baseline %>% filter(rowSums(select(., all_of(mri_scans_structural)) == 1) == length(mri_scans_structural));
samples_2yr_allpass <- samples_2yr %>% filter(rowSums(select(., all_of(mri_scans_structural)) == 1) == length(mri_scans_structural));
samples_4yr_allpass <- samples_4yr %>% filter(rowSums(select(., all_of(mri_scans_structural)) == 1) == length(mri_scans_structural));
nrow(samples_baseline_allpass);
nrow(samples_2yr_allpass);
nrow(samples_4yr_allpass);

samples_alltimepoints_allpass <- intersect(
  x = intersect(
    x = samples_baseline_allpass$src_subject_id,
    y = samples_2yr_allpass$src_subject_id),
  y = samples_4yr_allpass$src_subject_id
  );
length(samples_alltimepoints_allpass);

qc <- 'ALLPASS_structural';
venn.diagram(
  x = list(samples_baseline_allpass$src_subject_id, samples_2yr_allpass$src_subject_id, samples_4yr_allpass$src_subject_id),
  category.names = c("BASELINE" , "2YR" , "4YR"),
  filename = paste0(out_dir, 'mri_qc_', qc, '.tiff'),
  diable_logging = TRUE
  );

for (scan in mri_scans_structural) {
  qc <- unlist(strsplit(scan, split = '_'))[2];
  print(scan);
  samples_baseline_pass <- samples_baseline[samples_baseline[, scan] == 1,]$src_subject_id;
  samples_2yr_pass <- samples_2yr[samples_2yr[, scan] == 1,]$src_subject_id;
  samples_4yr_pass <- samples_4yr[samples_4yr[, scan] == 1,]$src_subject_id;
  venn.diagram(
    x = list(samples_baseline_pass, samples_2yr_pass, samples_4yr_pass),
    category.names = c("BASELINE" , "2YR" , "4YR"),
    filename = paste0(out_dir, 'mri_qc_', qc, '.tiff'),
    diable_logging = TRUE
    );
}