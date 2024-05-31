#### Format phenotype and covariate files ####
rm(list=ls())

#module load gcc/10.2.0
#module load R/4.3.0
#module load libxml2/2.10.4

# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("FRGEpistasis")

# libraries
suppressPackageStartupMessages(suppressWarnings({
    library(tidyverse)
    library(FRGEpistasis)
    library(patchwork)
}))

#### Directories and Files ####
date <- "03122024"
batchinfodir <- "/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_5/genomics_sample03/smokescreen"
qcdir <- "/u/project/lhernand/sganesh/gwas_srs/basic_qc_plink/ABCD_Release_5_QC_plinkfiles"
ancestry_id_dir <- "/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed"
pcdir <- "/u/project/lhernand/sganesh/gwas_srs/pc_unrelated"
phenodir <- "/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_4/phenotypes"
pheno_qc_dir <- "/u/project/lhernand/sganesh/gwas_srs/phenotypes/unrelated"
plotdir <- paste(pheno_qc_dir, "plots", sep = "/")
if (!dir.exists(plotdir)) {dir.create(plotdir)}
scores <- list(srs = "ssrs_42_p")

#### Read in phenotype files ####
setwd(phenodir)
# check file names
list.files(pattern = "ssrs")

# read in the files
files <- list.files(pattern = "abcd_pssrs01.txt")
names(files) <- c("ssrs")

pheno <- list()
for (file in files) {
    if (file == "abcd_pssrs01.txt") {
        score = scores[["srs"]]
    }
    print(file)
    print(score)
    dt <- read.delim(file, header = TRUE, na.strings = c("", "NA"))
    dt <- dt[-c(1),]
    dt <- dplyr::select(dt, subjectkey, eventname, all_of(score), interview_age, sex)

    age_col <- "interview_age"        
    dt[,c(score, age_col)] <- lapply(dt[,c(score, age_col)], as.numeric)
    pheno[[file]] <- dt
}

# check the data
lapply(pheno, head, n = 3)
lapply(pheno, nrow)

#. quick check on the pheno data ----
for (i in seq_along(pheno)) {
    dt <- pheno[[i]]
    score <- scores[[i]]
    
    # https://stackoverflow.com/questions/63403373/r-how-to-summarize-and-group-by-variables-as-column-names
    print(
        dt %>% group_by(eventname, sex) %>% 
            summarise(n = n(),
                      across(interview_age, mean, na.rm = TRUE), 
                      across(all_of(score), mean, na.rm = TRUE), .groups = 'drop')
        )
}
# SRS has only for 1_year_follow_up_y_arm_1
# Age are the same across measures in ABCD study as long as using the same timepoint/eventname.
# For ASD_GWAS project, we will use `1_year_follow_up_y_arm_1` SRS.

for (i in seq_along(pheno)) {
    print(names(pheno[i]))
    dt <- pheno[[i]]
    print(paste("nrow(dt_original):", nrow(dt)))
    dt <- dt %>% filter(eventname == "1_year_follow_up_y_arm_1")
    print(paste("nrow(dt_1yr_followup):", nrow(dt)))
    pheno[[i]] <- dt
    cat("\n")
}

#### Tidy up phenodata ####
for (i in seq_along(pheno)){
    
    dt <- pheno[[i]]
    name <- names(pheno[i])
    samplesize_original <- nrow(dt)
    
    cat("########## ", toupper(name), " ##########\n")
    
    ### Remove samples have NAs in the target score
    cat("### Removing NAs in the target score column ###\n")
    score <- ifelse(str_detect(name, "ssrs"), scores[[1]], NA)
    nrow_before <- nrow(dt)
    print(paste0("samplesize before: ", nrow_before))
    
    print("NAs in the target score:")
    print(table(is.na(dt[,score])))
    cat("\n")

    # remove NAs
    keep <- !is.na(dt[,score])
    dt <- dt[keep,]
    
    nrow_after <- nrow(dt)
    diff <- nrow_before - nrow_after
    
    print(paste0("samplesize after: ", nrow_after))
    print(paste0(diff, " samples removed (NAs in the score column)"))
    cat("\n")
    
    
    ### Remove samples have NAs in the age column
    cat("### Removing NAs in the age column ###\n")
    age_col <- "interview_age"
    nrow_before <- nrow(dt)
    print(paste0("samplesize before: ", nrow_before))
    
    print("NAs in the target score:")
    print(table(is.na(dt[,age_col])))
    cat("\n")

    # remove NAs
    keep <- !is.na(dt[,age_col])
    dt <- dt[keep,]
    
    nrow_after <- nrow(dt)
    diff <- nrow_before - nrow_after
    
    print(paste0("samplesize after: ", nrow_after))
    print(paste0(diff, " samples removed (NAs in the age column)"))
    cat("\n")    
      
    
    ### All: get necessary columns
    cat("### Subsetting to necessary columns ###\n")
    print(paste0("ncol() before: ", ncol(dt)))
    dt <- dt %>% dplyr::select(subjectkey, all_of(age_col), sex, all_of(score))
    print(paste0("ncol() after: ", ncol(dt)))
    cat("\n")
    
    
    ### print final dimention of the data
    cat("### Sample size ###\n")
    samplesize_processed <- nrow(dt)
    diff <- samplesize_original - samplesize_processed
    print(paste0("original: ", samplesize_original))
    print(paste0("processed: ", nrow(dt)))
    print(paste0("unique samples: ", length(dt$subjectkey)))
    print(paste0(diff, " samples removed through the process"))
    cat("\n\n")
    
    pheno[[i]] <- dt
}

# nrow() after tiding up the data
lapply(pheno, nrow)


#### Modify IDs with ABCD ID format ####
setwd(batchinfodir)
Batch <- read.delim("ABCD_202209.updated.nodups.curated.batch.info",  header = TRUE, na.strings = c("", "NA"), sep = ",")
#Batch<- rename(Batch, subjectkey = IID)
Batch$IID2 <- Batch$IID



#### Split by ancestry and gender ####
setwd(ancestry_id_dir)
list.files(pattern = "IDs")


# IDs for each ancestry group
afr <- read.table("/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AFR.2263_no.sexmismatch_IDs.txt")
amr <- read.table("/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_AMR.2019_no.sexmismatch_IDs.txt")
eas <- read.table("/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_EAS.382_no.sexmismatch_IDs.txt")
eur <- read.table("/u/project/lhernand/sganesh/gwas_srs/TOPMed_imputed/ABCDr5_EUR.6891_no.sexmismatch_IDs.txt")
#sas <- read.table("ABCDr5_SAS.92_no.sexmismatch_IDs.txt")

ancestry_samples_id_list <- list(AFR = afr$V2, AMR = amr$V2, EAS = eas$V2, EUR = eur$V2, SAS = sas$V2)

# before splitting into ancestry groups
lapply(ancestry_samples_id_list, head, n = 3)
lapply(pheno, dim)

# split by ancestry
pheno_by_ancestry <- list()
pheno_data <- pheno[[1]]
pheno_name <- tolower(names(pheno[1]))
for (j in seq_along(ancestry_samples_id_list)) {
    ids <- ancestry_samples_id_list[[j]]
#    print(head(ids),3)
    ancestry_name <- tolower(names(ancestry_samples_id_list[j]))
#    print(head(ancestry_name))
    name <- paste(pheno_name, ancestry_name, sep = "_")
#    print(head(name))
    pheno_by_ancestry[[name]] <- pheno_data %>% filter(subjectkey %in% ids)  
}

# after splitting into ancestry groups
lapply(pheno_by_ancestry, dim)

# split phenodata to ancestry groups and by gender
seq_along(pheno_by_ancestry)
pheno_by_ancestry_and_gender <- list()
gender <- c("M", "F")
for (i in seq_along(pheno_by_ancestry)) {
    pheno_data <- pheno_by_ancestry[[i]]
    pheno_name <- tolower(names(pheno_by_ancestry[i]))
    
    for (j in seq_along(gender)) {
        name <- paste(pheno_name, tolower(gender[j]), sep = "_")
        pheno_by_ancestry_and_gender[[name]] <- pheno_data %>% filter(sex == gender[j])
    }
}

# after splitting into ancestry groups and gender
lapply(pheno_by_ancestry_and_gender, dim)
length(pheno_by_ancestry_and_gender) #10

data <- c(pheno_by_ancestry, pheno_by_ancestry_and_gender)
length(data) #15
names(pheno)
head(pheno_by_ancestry[[1]])

#. save tables for each phenotype containing ancestry info ----
pheno_by_ancestry_before_srs <- list()
n <- length(pheno_by_ancestry) # 25 (5 pheno * 5 ancestry groups)
for (i in seq_along(pheno)){
    pheno_name <- names(pheno[i])
    keyword <- paste0(pheno_name, "_")
    pheno_dt_list <- data[1:n][str_detect(names(data[1:n]), keyword)]
    pheno_by_ancestry_before_srs[[pheno_name]] <- do.call("rbind", pheno_dt_list)
}
lapply(pheno_by_ancestry_before_srs, head, n = 3)
lapply(pheno_by_ancestry_before_srs, nrow) #11023

# check if all samples are unique
lapply(pheno_by_ancestry_before_srs, function(dt){
    print(head(dt$subjectkey),3)
    return(length(unique(dt$subjectkey)) == nrow(dt))
}) #true
 
setwd(pheno_qc_dir)
tables <- pheno_by_ancestry_before_srs
for (i in seq_along(tables)) {
    pheno_name <- names(tables[i])
    print(head(pheno_name))
    dt <- tables[[i]]
    print(head(dt))
    n <- nrow(dt)
    file_name <- paste0(pheno_name, "_1yr_followup_before.srs_within_ancestry_group_noNAs_", n, "_table_", date, ".csv")
    print(file_name)
    write.csv(dt, file_name, row.names = TRUE)
}

#### Rank transform phenotypes ####
for (i in seq_along(data)){
    dt <- data[[i]]
    name <- names(data[i])    
    score <- ifelse(str_detect(name, "ssrs"), scores[[1]],NA)
    score_srs <- paste(score, "srs", sep = "_")
    age_col <- "interview_age"
    age_col_years <- "interview_age_years"
    age_col_years_nodecimal <- "interview_age_years_nodecimal"
    n <- nrow(dt)
    
    dt[,c(score, age_col)] <- lapply(dt[,c(score, age_col)], as.numeric)
    dt[,age_col_years] <- round(dt[,age_col]/12, 2)
    dt[,age_col_years_nodecimal] <- round(dt[,age_col]/12, 0)
    dt[,score_srs] <- rankTransPheno(dt[,score], 0.5)
    data[[i]] <- dt
}

#. histogram ----
summary(data[[1]]$interview_age_years_nodecimal)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  10.00   10.00   11.00   10.91   11.00   12.00 

summary(factor(data[[1]]$interview_age_years_nodecimal))
# 10:681 11:866 12:495 

# histograms
setwd(plotdir)
for (i in seq_along(data)){
    dt <- data[[i]]
    name <- names(data[i])    
    score <- ifelse(str_detect(name, "ssrs"), scores[[1]], NA)
    score_srs <- paste(score, "srs", sep = "_")
    age_col <- "interview_age_years"
    n <- nrow(dt)

    png(paste("hist", name, score, n, "png", sep = "."), width = 800, height = 800)
    par(mfrow = c(2,2))
    # score
    hist(dt[,score], main = paste0("Histogram: " , toupper(name), " (n=", n, ")", "\n", score), xlab = score, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    hist(dt[,score_srs], main = paste0("Histogram: " , toupper(name), " (n=", n, ")", "\n", score_srs), xlab = score_srs, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    # age in months
    hist(dt[,age_col], main = paste0("Histogram: " , toupper(name), " (n=", n, ")", "\n", age_col), xlab = age_col, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    dev.off()    

    cat("\n\n")
}

# histograms for each ancestry by gender ----
data <- lapply(data, function(dt){
    dt$sex <- factor(dt$sex, levels = c("M", "F"))
    return(dt)
})

for (i in seq_along(pheno_by_ancestry)){
    dt <- data[[i]]
    name <- names(data[i])    
    score <- ifelse(str_detect(name, "srs"), scores[[1]], NA)
    
    score_srs <- paste(score, "srs", sep = "_")
    n <- nrow(dt)
    
    p1 <- ggplot(dt, aes(x = dt[,score], fill = sex, color = sex)) +
        geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
        theme_bw() +
        labs(title = paste("Histogram:", toupper(name)),
             subtitle = paste0(score, " (n=", n, ")"), 
             x        = "Scores",
             y        = "Freqency") +
        theme(text = element_text(size = 16), legend.position='none', panel.grid = element_blank()) +
        scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
        scale_color_manual(values = c("#00BFC4", "#F8766D")) +
        guides(fill = guide_legend("Gender"), color = guide_legend("Gender"))

    p2 <- ggplot(dt, aes(x = dt[,score_srs], fill = sex, color = sex)) +
        geom_histogram(position="identity", alpha=0.5, binwidth = 0.5) +
        theme_bw() +
        labs(title = paste("Histogram:", toupper(name)),
             subtitle = paste0(score_srs, " (n= ", n, ")"), 
             x        = "Scores",
             y        = "") +
        guides(fill = guide_legend("Gender"), color = guide_legend("Gender")) +
        theme(text = element_text(size = 16), legend.position='none', panel.grid = element_blank()) +
        scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
        scale_color_manual(values = c("#00BFC4", "#F8766D")) +
        guides(fill = guide_legend("Gender"), color = guide_legend("Gender"))

#     suppressWarnings(print(p1 + p2 + plot_layout(guides = "collect") & 
#                 theme(legend.position = "right", legend.text = element_text(size = 14))))

    png(paste("hist", name, score, "by_gender", n, "png", sep = "."), width = 800, height = 400)    
    suppressWarnings(print(p1 + p2 + plot_layout(guides = "collect") & 
            theme(legend.position = "right", legend.text = element_text(size = 14))))
    dev.off()
}

# Histograms for each ancestry by age ----
data <- lapply(data, function(dt){
    age_median <- median(dt$interview_age_years)
    dt$age_median <- ifelse(dt$interview_age_years > age_median, "top50%", "bottom50%")
    dt$age_median <- factor(dt$age_median, levels = c("top50%", "bottom50%"))
    return(dt)
})

#. histograms for each ancestry by age ----
data <- lapply(data, function(dt){
    age_median <- median(dt$interview_age_years)
    dt$age_median <- ifelse(dt$interview_age_years > age_median, "top50%", "bottom50%")
    dt$age_median <- factor(dt$age_median, levels = c("top50%", "bottom50%"))
    return(dt)
})

for (i in seq_along(pheno_by_ancestry)){
    dt <- data[[i]]
    name <- names(data[i])    
    score <- ifelse(str_detect(name, "srs"), scores[[1]], NA)
    
    score_srs <- paste(score, "srs", sep = "_")
    n <- nrow(dt)
    
    p1 <- ggplot(dt, aes(x = dt[,score], fill = age_median, color = age_median)) +
        geom_histogram(position="identity", alpha=0.5, binwidth = 5) +
        theme_bw() +
        labs(title = paste("Histogram:", toupper(name)),
             subtitle = paste0(score, " (n=", n, ")"), 
             x        = "Scores",
             y        = "Freqency") +
        theme(text = element_text(size = 16), legend.position='none', panel.grid = element_blank()) +
        scale_fill_manual(values = c("orange", "skyblue1")) +
        scale_color_manual(values = c("orange", "skyblue1")) +
        guides(fill = guide_legend("Age by Median"), color = guide_legend("Age by Median"))

    p2 <- ggplot(dt, aes(x = dt[,score_srs], fill = age_median, color = age_median)) +
        geom_histogram(position="identity", alpha=0.5, binwidth = 0.5) +
        theme_bw() +
        labs(title = paste("Histogram:", toupper(name)),
             subtitle = paste0(score_srs, " (n= ", n, ")"), 
             x        = "Scores",
             y        = "") +
        guides(fill = guide_legend("Age by Median"), color = guide_legend("Age by Median")) +
        theme(text = element_text(size = 16), legend.position='none', panel.grid = element_blank()) +
        scale_fill_manual(values = c("orange", "skyblue1")) +
        scale_color_manual(values = c("orange", "skyblue1")) +
        guides(fill = guide_legend("Age by Median"), color = guide_legend("Age by Median"))

#     suppressWarnings(print(p1 + p2 + plot_layout(guides = "collect") & 
#                 theme(legend.position = "right", legend.text = element_text(size = 14))))

    png(paste("hist", name, score, "by_age_median", n, "png", sep = "."), width = 800, height = 400)    
    suppressWarnings(print(p1 + p2 + plot_layout(guides = "collect") & 
            theme(legend.position = "right", legend.text = element_text(size = 14))))
    dev.off()
}

#### Combine phenodata ####
#. by ancestry ----
pheno_by_ancestry_srs <- list()
n <- length(pheno_by_ancestry) 
for (i in seq_along(pheno)){
    pheno_name <- names(pheno[i])
    keyword <- paste0(pheno_name, "_")
    pheno_dt_list <- data[1:n][str_detect(names(data[1:n]), keyword)]
    pheno_by_ancestry_srs[[pheno_name]] <- do.call("rbind", pheno_dt_list)
}
lapply(pheno_by_ancestry_srs, nrow) #85
lapply(pheno_by_ancestry_srs, head, n = 3)


#. check if all samples are unique ----
lapply(pheno_by_ancestry_srs, function(dt){
    return(length(unique(dt$subjectkey)) == nrow(dt))
}) # TRUE

# Age Summary
summary(dt$interview_age_years)
# $srs
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   9.92   10.33   11.00   10.95   11.50   12.00 

setwd(pheno_qc_dir)
for (i in seq_along(pheno_by_ancestry_srs)) {
    pheno_name <- names(pheno_by_ancestry_srs[i])
    dt <- pheno_by_ancestry_srs[[i]]
    n <- nrow(dt)
    file_name <- paste0(pheno_name, "_1yr_followup_after.srs_within_ancestry_group_noNAs_age.as.is_", n, "_table_", date, ".csv")
    print(file_name)
    write.csv(dt, file_name, row.names = TRUE)
}

#. save phenotype files ----
scores <- list(srs = "ssrs_42_p")
for (i in seq_along(pheno_by_ancestry_srs)) {
    pheno_name <- names(pheno_by_ancestry_srs[i])
    dt <- pheno_by_ancestry_srs[[i]]
    score <- ifelse(str_detect(pheno_name, "ssrs"), scores[[1]], NA)
    n <- nrow(dt)
    file_name <- paste0("Phen_AllSubj_", pheno_name, "_1yr_followup_", score, "_within_ancestry_group_noNAs_", n, ".txt")
    print(file_name)

    dt <- dt %>% dplyr::select(subjectkey, subjectkey, all_of(score))
    print(head(dt, 3))
    write.table(dt, file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    cat("\n")
}

#### Create qcovar files ####
# read in ancestry PCs
setwd(pcdir)
files <- list.files(pattern = "eigenvec")
files

#### IDs for each ancestry group ####
pcs <- list()
for (i in seq_along(files)) {
    file <- files[i]
    ancestry <- tolower(str_extract(file, "AFR|AMR|EAS|EUR|SAS"))
    pcs[[ancestry]] <- read.table(file, na.strings = c("", "NA"))
}
lapply(pcs, head, n = 3)

pcs <- do.call("rbind", pcs)
dim(pcs)

names(pcs) <- c("fid", "iid", paste0(rep("PC", 10), 1:10))
head(pcs, 3)
tail(pcs, 3)

#. save quantitative covariate files ----
setwd(pheno_qc_dir)
for (i in seq_along(pheno_by_ancestry_srs)) {
    par(mfrow = c(1,2))
    pheno_name <- names(pheno_by_ancestry_srs[i])
    print(pheno_name)

    dt <- pheno_by_ancestry_srs[[i]]
    dt$ancestry <- str_extract(rownames(dt), "afr|amr|eas|eur|sas")
    dt_pcs <- merge(dt, pcs, by.x = "subjectkey", by.y = "iid")
    print(paste("nrow(dt_pcs):", nrow(dt_pcs)))
    
    # PCs1-10
    dt_pcs10 <- dt_pcs %>% dplyr::select(subjectkey,subjectkey, interview_age_years_nodecimal, paste0(rep("PC", 10), 1:10))
    print(head(dt_pcs10, 3))
    n <- nrow(dt_pcs10)
    print(n)

    file_name <- paste0("qcovar_AllSubj_", pheno_name, "_1yr_followup_age.as.is_PCs1-10_within_ancestry_group_noNAs_", n, ".txt")
    print(file_name)
    write.table(dt_pcs10, file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    cat("\n")
}

#### Create discrete covariate files ####
head(Batch, 3)
#names(Batch)[4] <- "batch"

#. merge with gender info ----
setwd(phenodir)
dt_tracking <- read.delim("abcd_lt01.txt", na.strings = c("", "NA"))
dt_tracking <- dt_tracking[-1,]
head(dt_tracking, 3)

# get baseline data
dt_tracking <- dt_tracking %>% filter(eventname == "baseline_year_1_arm_1") %>% dplyr::select(subjectkey, interview_age, sex)
covar <- merge(Batch, dt_tracking, by.x = "IID", by.y = "subjectkey")
nrow(covar)
head(covar, 3)

# format covariate file
covar <- covar %>% dplyr::select(IID, subjectkey, BATCH, sex)
head(covar, 3)

summary(factor(covar$sex))
#   F    M 
# 5585 6080 

summary(factor(covar$batch))
#BATCH_1_Saliva     BATCH_1_WB BATCH_2_Saliva     BATCH_2_WB BATCH_3_Saliva 
#           617             11           4358            194           4020 
#    BATCH_3_WB BATCH_4_Saliva        BATCH_5 BATCH_5_Saliva 
#            84           1860            318            203 

#. save discrete covariate file ----
setwd(pheno_qc_dir)

# batch, gender
batch_gender <- covar
n <- nrow(batch_gender)
filename <- paste0("covar_AllSubj_batch_gender_noNAs_baseline_", n, ".txt")
write.table(batch_gender, filename, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# Write per ancestry phenotype file for GCTA
pheno <- read.table("Phen_AllSubj_abcd_pssrs01.txt_1yr_followup_ssrs_42_p_within_ancestry_group_noNAs_11023.txt", sep = "\t", header=FALSE)
setwd(ancestry_id_dir)
sas <- read.table("ABCDr5_SAS.92_no.sexmismatch_IDs.txt")
sas_pheno <- merge(sas,pheno, by="V1")
sas_pheno <- sas_pheno[,c("V1","V2.y")]
setwd(pheno_qc_dir)
write.table(sas_pheno, "Phen_SAS_abcd_pssrs01_1yr_followup_ssrs_42_p.txt", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)

covar$V4 <- covar$V1
covar <- covar[,c(4,1,2,3)]
write.table(covar, "covar_AllSubj_batch_gender_noNAs_baseline_11665.txt", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
