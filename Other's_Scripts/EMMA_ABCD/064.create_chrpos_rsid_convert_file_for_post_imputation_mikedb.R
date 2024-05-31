#### Create chr:pos -> rsID replacement file for post-imputation process ####

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


#### directories and files ####
# dir <- "/Users/Emma/Hoffman"
dir <- "/u/project/gandalm"
filedir <- paste(dir, "shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/vcf_imputed_filtered_plink", sep = "/")
post_imputation_dir <- paste(dir, "emmamk/ASD_GWAS/ABCDr4/06.post_imputation", sep = "/")
outdir <- paste(post_imputation_dir, "chrpos_rsid_biallelic", sep = "/")
if(!dir.exists(outdir)){dir.create(outdir)}
list.files(filedir, pattern = ".bim")


#### Query Michael's Database ####
### CHANGE GENOME BUILD ####
build <- "GRCh38"
dir <- paste(mikedbdir, "03-parquet-db", build, sep = "/")
db_SNP <- open_dataset(dir)

for (i in 1:22) {
    print(paste("#### chr", i, "####"))
    print("Loading in the file ...")
    # load in the file
    file <- paste0("chr", i, ".imp.flt.r208.bim")
    file <- paste(filedir, file, sep = "/")

    data <- read.delim(file, header = FALSE, na.strings = c("", "NA"))
    names(data)[c(1:2,4:6)] <- c("chr", "chrpos_refalt", "pos", "alt", "ref")


    # Subset to necessary columns
    query_tbl <- data %>% dplyr::select(chr, pos, ref, alt)


    # Query Michael's Database
    print("Started the query ...")
    dbsnp_res <- query_for_rsid(db_SNP, query_tbl)

    print(paste("nrow(query_tbl):", nrow(query_tbl)))
    print(paste("nrow(dbsnp_res):", nrow(dbsnp_res)))
    diff <- nrow(query_tbl) - nrow(dbsnp_res)
    print(paste("nrow(query_tbl) - nrow(dbsnp_res):", diff))
    

    # Create necessary columns and subset
    dbsnp_res <- dbsnp_res %>% dplyr::select(chr, pos, ref, alt, id)
    dbsnp_res$chrpos <- paste(paste0("chr", dbsnp_res$chr), dbsnp_res$pos, sep = ":")
    dbsnp_res$chrpos_refalt <- paste(paste0("chr", dbsnp_res$chr), dbsnp_res$pos, dbsnp_res$ref, dbsnp_res$alt, sep = ":")
    dbsnp_res$rsid_num <- str_replace(dbsnp_res$id, "rs", "")


    # Merge datasets
    # Keep SNPs with rsIDs
    data_merged <- merge(data, dbsnp_res, by = "chrpos_refalt", all.y= TRUE)

    # check if any dup SNPs - chr:pos:ref:alt == chr:pos:ref:alt
    dup_SNPs <- data_merged %>% dplyr::select(chrpos_refalt, id) %>% 
          group_by(chrpos_refalt) %>% 
          filter(n() > 1)
    print(paste("dup SNPs: ", nrow(dup_SNPs)))

    # remove dup SNPs (chr:pos:ref:alt)
    data_merged <- data_merged %>% filter(!(id %in% dup_SNPs$id) & id != "NA") %>%
          dplyr::select(chrpos, id, rsid_num)

    # if chr:pos dups; keep the 1st rsID (younger one)
    data_nodup <- data_merged %>% group_by(chrpos) %>% 
        arrange(rsid_num) %>% 
        slice_head(n = 1) %>%
        ungroup() %>% 
        dplyr::select("chrpos", "id")

    uniq_nrow <- nrow(data_nodup)
    print(paste("unique nrow():", uniq_nrow))
    num_removed_variants <- nrow(data_merged) - uniq_nrow
    print(paste("number of removed dup variants:", num_removed_variants))

    # Save the table
    print(paste("saving the file ..."))
    n <- nrow(data_nodup)
    file <- paste0("chrpos_rsid_biallelic_chr", i, "_", n, ".txt")
    file <- paste(outdir, file, sep = "/")
    print(file)
    write.table(data_nodup, file, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = FALSE)
    cat("\n\n")
}

# Concatenated files using shell script
# cd /u/project/gandalm/emmamk/ASD_GWAS/ABCDr4/06.post_imputation/chrpos_rsid_biallelic
# cat chrpos_rsid_biallelic_chr{1..22}_*txt > chrpos_rsid_biallelic_chr1-22.txt
