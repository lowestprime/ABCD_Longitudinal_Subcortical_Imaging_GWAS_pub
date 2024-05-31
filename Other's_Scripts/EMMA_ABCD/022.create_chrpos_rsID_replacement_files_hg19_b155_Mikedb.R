#### Create chrpos_rsID replacement files for genotyped data ####
rm(list = ls())

# libraries
suppressPackageStartupMessages(suppressWarnings({
    library(tidyverse)
    library(glue)
    library(janitor)
    library(arrow)
    library(tictoc)
}))
mikedbdir <- "/u/project/gandalm/margolis/alexandria/refs/dbsnp/v155"
source("/u/project/gandalm/emmamk/tools/MichaelDB_dbSNPquery/zz-functions.R")


#### directory ####
qcdir <- "/u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/02.QCplink"


#### Read in target .bim file ####
# format a data frame for the target file
file <- "/u/project/gandalm/shared/GenomicDatasets/ABCD_Release_4/genomics_sample03/genotype_QCed/ABCD_release_3.0_QCed.bim"
dt <- read.delim(file, header = FALSE, na.strings = c("", "NA"))
dt <- dt[-3]
names(dt) <- c("chr", "id", "pos", "alt", "ref")
dt$chrpos_refalt <- paste(dt$chr, dt$pos, dt$ref, dt$alt, sep = ":")
dt$chrpos_altref <- paste(dt$chr, dt$pos, dt$alt, dt$ref, sep = ":")

#. query for the original bim file ----
# build the query table
dt_chrpos_altref <- dt %>% select(chr, pos, ref, alt)

build <- "GRCh37"
dir <- paste(mikedbdir, "03-parquet-db", build, sep = "/")
db_SNP <- open_dataset(dir)

# query the dbSNP database
dbsnp_res_chrpos_altref <- query_for_rsid(db_SNP, dt_chrpos_altref)

dbsnp_res_chrpos_altref$checked_allele_order <- "alt.alt_ref.ref"
dbsnp_res_chrpos_altref$estimated_allele_order <- "original"
dbsnp_res_chrpos_altref$chrpos_refalt_original <- paste(dbsnp_res_chrpos_altref$chr, dbsnp_res_chrpos_altref$pos, dbsnp_res_chrpos_altref$ref, dbsnp_res_chrpos_altref$alt, sep = ":")
dbsnp_res_chrpos_altref <- dplyr::select(dbsnp_res_chrpos_altref, chrpos_refalt_original, chr, pos, ref, alt, id, checked_allele_order, estimated_allele_order)
names(dbsnp_res_chrpos_altref)[1] <- "chrpos_refalt"

nrow(dbsnp_res_chrpos_altref)
# 395884

# save the results
setwd(qcdir)
write.table(dbsnp_res_chrpos_altref, "221117_ABCDr4_originalbim_variants_rsIDs_chrposrefalt_originalorder_mikedb.txt",
    quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)

# check duplicated SNPs
dbsnp_res_chrpos_altref %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    nrow()
# 69    

dbsnp_res_chrpos_altref %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    head()

dt %>% filter(chr == "1", pos == "20806974")

dbsnp_res_chrpos_altref_dupSNPs <- dbsnp_res_chrpos_altref %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    dplyr::select(chrpos_refalt) %>%
    unique()
nrow(dbsnp_res_chrpos_altref_dupSNPs)
head(dbsnp_res_chrpos_altref_dupSNPs, 3)

# use original bim info for these SNPs
dt$estimated_allele_order[dt$chrpos_refalt %in% dbsnp_res_chrpos_altref_dupSNPs$chrpos_refalt] <- "mapped_original_dup_use_original_bim"
dt$checked_allele_order[dt$chrpos_refalt %in% dbsnp_res_chrpos_altref_dupSNPs$chrpos_refalt] <- "alt.alt_ref.ref"

# remove duplicated SNPs
dbsnp_res_chrpos_altref_nodup <- dbsnp_res_chrpos_altref %>% filter(!(chrpos_refalt %in% dbsnp_res_chrpos_altref_dupSNPs$chrpos_refalt))
nrow(dbsnp_res_chrpos_altref_nodup)

nrow(dbsnp_res_chrpos_altref)
dbsnp_res_chrpos_altref %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    nrow()

395884 - 69


#. query for the flipped alleles for the bim file ----
dt_chrpos_refalt <- dt %>% select(chr, pos, ref, alt)
names(dt_chrpos_refalt)[3:4] <- c("alt", "ref") # assuming these are A1/A2 flipped variants

build <- "GRCh37"
dir <- paste(mikedbdir, "03-parquet-db", build, sep = "/")
db_SNP <- open_dataset(dir)

# query the dbSNP database
dbsnp_res_chrpos_refalt <- query_for_rsid(db_SNP, dt_chrpos_refalt)

# check the results
dbsnp_res_chrpos_refalt$checked_allele_order <- "alt.ref_ref.alt"
dbsnp_res_chrpos_refalt$estimated_allele_order <- "flipped"
dbsnp_res_chrpos_refalt$chrpos_refalt_flipped <- paste(dbsnp_res_chrpos_refalt$chr, dbsnp_res_chrpos_refalt$pos, dbsnp_res_chrpos_refalt$ref, dbsnp_res_chrpos_refalt$alt, sep = ":")
dbsnp_res_chrpos_refalt <- dplyr::select(dbsnp_res_chrpos_refalt, chrpos_refalt_flipped, chr, pos, ref, alt, id, checked_allele_order, estimated_allele_order)
names(dbsnp_res_chrpos_refalt)[1] <- "chrpos_refalt"

nrow(dbsnp_res_chrpos_refalt)
head(dbsnp_res_chrpos_refalt, 3)

# check duplicated SNPs
dbsnp_res_chrpos_refalt %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    nrow()
dbsnp_res_chrpos_refalt %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    head()

dt %>% filter(chr == "8", pos == "99053590")

dbsnp_res_chrpos_refalt_dupSNPs <- dbsnp_res_chrpos_refalt %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    dplyr::select(chrpos_refalt, alt, ref) %>%
    unique()
nrow(dbsnp_res_chrpos_refalt_dupSNPs)
head(dbsnp_res_chrpos_refalt_dupSNPs, 3)

dbsnp_res_chrpos_refalt_dupSNPs

# use original rsID from the bim file for these SNPs
for (i in 1:nrow(dbsnp_res_chrpos_refalt_dupSNPs)) {
    chrpos_refalt <- dbsnp_res_chrpos_refalt_dupSNPs[i, ]$chrpos_refalt
    dt[dt$chrpos_altref == chrpos_refalt, ]$estimated_allele_order <- "mapped_flipped_dup_use_flipped_alleles"
    dt[dt$chrpos_altref == chrpos_refalt, ]$checked_allele_order <- "alt.ref_ref.alt"
    dt[dt$chrpos_altref == chrpos_refalt, ]$alt <- dbsnp_res_chrpos_refalt_dupSNPs[dbsnp_res_chrpos_refalt_dupSNPs$chrpos_refalt == chrpos_refalt, ]$alt
    dt[dt$chrpos_altref == chrpos_refalt, ]$ref <- dbsnp_res_chrpos_refalt_dupSNPs[dbsnp_res_chrpos_refalt_dupSNPs$chrpos_refalt == chrpos_refalt, ]$ref
}

# remove duplicated SNPs
dbsnp_res_chrpos_refalt_nodup <- dbsnp_res_chrpos_refalt %>% filter(!(chrpos_refalt %in% dbsnp_res_chrpos_refalt_dupSNPs$chrpos_refalt))
nrow(dbsnp_res_chrpos_refalt_nodup)
nrow(dbsnp_res_chrpos_refalt)
dbsnp_res_chrpos_refalt %>%
    group_by(chrpos_refalt) %>%
    filter(n() > 1) %>%
    nrow()

94107 - 18

# save the results
setwd(qcdir)
write.table(dbsnp_res_chrpos_refalt, "221117_ABCDr4_originalbim_variants_rsIDs_chrposrefalt_flippedorder_mikedb.txt",
    quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


#. get unmapped SNPs ---- 
dt_unmapped <- dt %>% filter(!(chrpos_refalt %in% dbsnp_res_chrpos_altref_nodup$chrpos_refalt) &
    !(chrpos_altref %in% dbsnp_res_chrpos_refalt_nodup$chrpos_refalt))
nrow(dt_unmapped)
head(dt_unmapped)

# mark "unmapped"
dt_unmapped$estimated_allele_order[is.na(dt_unmapped$estimated_allele_order)] <- "unmapped_use_originalbim"
dt_unmapped$checked_allele_order[is.na(dt_unmapped$checked_allele_order)] <- "unmapped"
head(dt_unmapped, 3)

summary(factor(dt_unmapped$estimated_allele_order))
# mapped_flipped_dup_use_flipped_alleles 6
# mapped_original_dup_use_original_bim 22
# unmapped_use_originalbim 26666

dt_unmapped_sub <- dt_unmapped %>% dplyr::select(names(dbsnp_res_chrpos_altref_nodup))


#### Save the results for resetting flipped alleles ####
dt_merged <- do.call("rbind", list(dbsnp_res_chrpos_altref_nodup, dbsnp_res_chrpos_refalt_nodup, dt_unmapped_sub))
summary(factor(dt_merged$estimated_allele_order))
# flipped 94089
# mapped_flipped_dup_use_flipped_alleles 6
# mapped_original_dup_use_original_bim 22
# original 395815
# unmapped_use_originalbim26666

# merge query table with the original bim file
dt_merged_comp_original <- merge(dt_merged, dt[, c(1, 3, 2)], by = c("chr", "pos"), all = TRUE)
head(dt_merged_comp_original)
names(dt_merged_comp_original)[c(6, 9)] <- c("id_queryied_and_unmapped_original", "id_original_bim")

setwd(qcdir)
write.table(dt_merged_comp_original, "221117_ABCDr4_originalbim_variants_rsIDs_original.flipped.unmapped.merged_mikedb.txt",
    quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


#### Save the results for variant ID replacement ####
dt_merged_comp_original_idmismatch <- dt_merged_comp_original %>% filter(id_queryied_and_unmapped_original != id_original_bim)
n <- nrow(dt_merged_comp_original_idmismatch)
filename <- paste0("ABCDr4_SNPid_replacement_mikedb155_", n, ".txt")
data <- dt_merged_comp_original_idmismatch[c("id_original_bim", "id_queryied_and_unmapped_original")]

setwd(qcdir)
write.table(data, filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE)

