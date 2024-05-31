#### Merge phenotypes from ABCD Release 5.1 longitudinal data  ####
# Cooper Beaman
# 5/30/24

# Specific column extraction (smri_vol_scs_intracranialv)
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
saveRDS(abcdData.R5.1, file = '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/ABCD.Release.5.1.Merged.Data.ICV.rds')

# Check the result
print(head(abcdData.R5.1))
if ('smri_vol_scs_intracranialv' %in% colnames(abcdData.R5.1)) {
  # Count the number of blank rows (NA or empty string)
  blank_count <- sum(is.na(abcdData.R5.1[['smri_vol_scs_intracranialv']]) | abcdData.R5.1[['smri_vol_scs_intracranialv']] == "")
  print(paste("Number of blank rows in column", 'smri_vol_scs_intracranialv', ":", blank_count))
} else {
  stop(paste("Column", 'smri_vol_scs_intracranialv', "does not exist in the dataframe."))
}