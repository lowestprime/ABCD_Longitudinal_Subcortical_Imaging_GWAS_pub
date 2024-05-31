#### Merge ABCD Release 5.1 longitudinal data with C4 GREx/ABCD baseline data ####
# Jack Dodson
# 3/14/24

library(tidyverse);

base_dir <- '/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5.1';
core_dir <- paste(base_dir, 'core', sep = '/');
substudies_dir <- paste(base_dir, 'substudies', sep = '/');

print(base_dir);
print(core_dir);
print(substudies_dir);

# load baseline merged data and extract variable names
load('/u/project/lhernand/jpdodson/ABCD/ABCD.Release.3.0.C4.Merged.Data.PhewasVars.Rda');

# create data frame with release 4 variable names and types
release_4_var_types <- sapply(abcdData, class);
release_4_vars <- data.frame(
    var_name = colnames(abcdData),
    class = unlist(unname(release_4_var_types))[-1],
    stringsAsFactors = FALSE
    );
extra_vars <- c('demo_comb_income_v2', 'demo_prnt_ed_v2', 'kbi_sex_assigned_at_birth', 'mri_info_deviceserialnumber');
extra_vars_types <- c('numeric', 'numeric', 'numeric', 'character');
release_4_vars <- rbind(
    release_4_vars,
    data.frame(var_name = extra_vars, class = extra_vars_types, stringsAsFactors = FALSE)
    );

# load release 5.1 data dictionary -- this will allow lookup of new tables
ABCD_Release5.1_data_dict <- read.delim(
    file = '/u/project/lhernand/jpdodson/ABCD/ABCD_Release_5.1_data_dictionary.csv',
    header = TRUE,
    sep = ','
    );

# merge release 4 and release 5 tables (DTI not included)
table_dict <- merge(release_4_vars, ABCD_Release5.1_data_dict, by = 'var_name');
table_dict <- table_dict %>% select(table_name_nda, table_name, var_name, class);
colnames(table_dict) <- c('table_name_r4', 'table_name_r5.1', 'var_name', 'class');
print(table_dict %>% filter(table_name_r5.1 == 'gish_y_gi'));
print(nrow(table_dict));

# get unique table names for release 5.1 and move longitudinal data to the front
release_5_tables <- unique(table_dict$table_name_r5.1);
longitudinal_table <- 'abcd_y_lt';
release_5_tables <- c(longitudinal_table, release_5_tables[release_5_tables != longitudinal_table]);
print(release_5_tables);

for (table in release_5_tables) {
    rel_path <- list.files(
        path = core_dir,
        pattern = table,
        recursive = TRUE
        );
    full_path <- paste(core_dir, rel_path, sep = '/');
    print(full_path);

    data <- read.delim(
        file = full_path,
        header = TRUE,
        na.strings = c('', 'NA'),
        sep = ','
        );
    
    vars_to_pull <- table_dict[which(table_dict$table_name_r5.1 == table), ]$var_name;
    if (table == 'mri_y_smr_thk_dsk') {
        vars_to_pull <- c(vars_to_pull, 'smri_thick_cdk_mean');
    }
    if (table == 'mri_y_smr_vol_dsk') {
        vars_to_pull <- c(vars_to_pull, 'smri_vol_cdk_total');
    }
    if (table == 'mri_y_smr_area_dsk') {
        vars_to_pull <- c(vars_to_pull, 'smri_area_cdk_total');
    }
    print(table_dict %>% filter(table_name_r5.1 == table));

    data <- data %>% select(src_subject_id, eventname, vars_to_pull);

    if ( !exists('abcdData.R5.1') ) {
        abcdData.R5.1 <- data;
    }
    else {
        abcdData.R5.1 <- merge(abcdData.R5.1, data, by = c('src_subject_id', 'eventname'), all.x = TRUE);
    }
}
print(colnames(abcdData.R5.1));
print(dim(abcdData.R5.1));
save(abcdData.R5.1, file = '/u/project/lhernand/jpdodson/ABCD/ABCD.Release.5.1.Merged.Data.Rda');

load('/u/project/lhernand/jpdodson/ABCD/ABCD.Release.5.1.Merged.Data.Rda');

# release_5_var_types <- sapply(abcdData.R5.1, class);
# release_5_vars <- data.frame(
#     var_name = colnames(abcdData.R5.1),
#     class = unlist(unname(release_5_var_types)),
#     stringsAsFactors = FALSE
#     );
# release_5_vars <- release_5_vars %>% filter(!var_name %in% c('src_subject_id', 'eventname'));

# can address later -- some variables are numeric in r4 and integer in r5. does this matter?
# for (var in release_5_vars$var_name) {
#     r4 <- release_4_vars %>% filter(var_name == var);
#     r5 <- release_5_vars %>% filter(var_name == var);
#     print(r4);
#     print(r5);
#     print(r4$class == r5$class);
# }

# Load C4 data and merge
ABCD.C4.imputed.multiethnic <- read.delim(
    file = '/u/project/lhernand/jpdodson/ABCD/ABCD.C4.imputed.multiethnic.txt',
    header = TRUE,
    na.strings = c('', 'NA')
    );
colnames(ABCD.C4.imputed.multiethnic)[1] <- 'src_subject_id'; 
ABCD.C4.imputed.multiethnic$eventname <- 'baseline_year_1_arm_1';
abcdData.R5.1 <- merge(
    x = abcdData.R5.1,
    y = ABCD.C4.imputed.multiethnic,
    by = c('src_subject_id', 'eventname'), 
    all.x = TRUE
    );

# Load and merge ancestry data
ancestry.pcs <- read.delim(
    file = '/u/project/lhernand/jpdodson/ABCD/ABCD_ancestry_knn.txt',
    header = TRUE,
    na.strings = c('', 'NA')
    );
colnames(ancestry.pcs) <- c('src_subject_id', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', 'pc9', 'pc10', 'pc11', 'pc12', 'pc13', 'pc14', 'pc15', 'pc16', 'pc17', 'pc18', 'pc19', 'pc20', 'data', 'ethnicity');
ancestry.pcs$src_subject_id <- lapply(
    X = ancestry.pcs$src_subject_id,
    FUN = function(x) paste0('NDAR_', unlist(strsplit(x, split = '_'))[3]));
ancestry.pcs$eventname <- 'baseline_year_1_arm_1';
abcdData.R5.1 <- merge(
    x = abcdData.R5.1,
    y = ancestry.pcs,
    by = c('src_subject_id', 'eventname'),
    all.x = TRUE
    );

abcdData.R5.1 <- abcdData.R5.1 %>%
  select(-rel_family_id.y, -kbi_sex_assigned_at_birth.y) %>%
  rename(rel_family_id = rel_family_id.x, sex = kbi_sex_assigned_at_birth.x);

abcdData.R5.1$ses_AvgIncomeParentEdu <- ((as.numeric(abcdData.R5.1$demo_comb_income_v2) + as.numeric(abcdData.R5.1$demo_prnt_ed_v2))/2);
abcdData.R5.1$pps_y_ss_number_plus1 <- (1 + abcdData.R5.1$pps_y_ss_number);
abcdData.R5.1$pps_y_ss_severity_score_plus1 <- (1 + abcdData.R5.1$pps_y_ss_severity_score);

# fill variables only measured at baseline
fill_var <- c('rel_family_id', 'race_ethnicity', 'sex', 'ses_AvgIncomeParentEdu', 'ssrs_p_ss_sum', 'c4a_total', 'c4b_total', 'c4herv_total', 'c4a_dosage', 'c4b_dosage', 'avg_dosage', 'c4_allele1_structure', 'c4_allele2_structure', 'c4a_expression', 'c4b_expression', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', 'pc9', 'pc10', 'pc11', 'pc12', 'pc13', 'pc14', 'pc15', 'pc16', 'pc17', 'pc18', 'pc19', 'pc20', 'data', 'ethnicity')
for(i in fill_var){
abcdData.R5.1 <- abcdData.R5.1 %>%
  dplyr::group_by(src_subject_id) %>%
  fill(i, .direction = 'downup') %>%
  dplyr::ungroup()
}

print(colnames(abcdData.R5.1));
print(dim(abcdData.R5.1));

# save(abcdData.R5.1, file = '/u/project/lhernand/jpdodson/ABCD/ABCD.Release.5.1.C4.Merged.Data.Rda');

#### specific column extraction (smri_vol_scs_intracranialv) ####
# Identifying the table containing `smri_vol_scs_intracranialv`
suppressPackageStartupMessages(suppressWarnings({
  p_load(tidyverse)
}))

# Define directories
base_dir <- '/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5.1'
core_dir <- paste(base_dir, 'core', sep = '/')

# Load existing merged data
load('/u/project/lhernand/jpdodson/ABCD/ABCD.Release.5.1.C4.Merged.Data.Rda')

# Find the relevant table containing 'smri_vol_scs_intracranialv'
table_dict <- read.delim(
  file = '/u/project/lhernand/jpdodson/ABCD/ABCD_Release_5.1_data_dictionary.csv',
  header = TRUE,
  sep = ','
)

# Identify the table containing the desired variable
smri_table <- table_dict %>%
  filter(var_name == 'smri_vol_scs_intracranialv') %>%
  pull(table_name) %>%
  unique()

# Construct the file path for the relevant table
rel_path <- list.files(
  path = core_dir,
  pattern = smri_table,
  recursive = TRUE
)
full_path <- paste(core_dir, rel_path, sep = '/')

# Load the data from the relevant table
smri_data <- read.delim(
  file = full_path,
  header = TRUE,
  na.strings = c('', 'NA'),
  sep = ','
) %>% 
  select(src_subject_id, eventname, smri_vol_scs_intracranialv)

# Merge with the existing dataframe
abcdData.R5.1 <- abcdData.R5.1 %>%
  left_join(smri_data, by = c('src_subject_id', 'eventname'))

# Rename and save the updated dataframe
abcdData.R5.1.ICV <- abcdData.R5.1
saveRDS(abcdData.R5.1.ICV, file = '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/ABCD.Release.5.1.Merged.Data.ICV.rds')

# Check the result
print(head(abcdData.R5.1.ICV))
print(dim(abcdData.R5.1.ICV))