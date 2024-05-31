#### Make PCA bi-plots for ancestry estimation results ####
rm(list=ls())

# options
options(scipen=100, digits=3)
options(stringsAsFactors = T)

# libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(class)
    library(tidyverse)
})

#### directories and files ####
dir <- "/u/project/gandalm"
ancestry_dir <- paste(dir, "emmamk/ASD_GWAS/ABCDr4/03.ancestry_estimation", sep = "/")
study <- "ABCD"
study_name <- "ABCD_release_3.0_QCed.NDARids.2idfixed.plate461rm.2badsubjects.rm.chr1.22.maf0.01.geno0.05.hwe1e6.mind0.1.refallele.reset_RSidOnly"
ref_name <- "ALL.autosomes.phase3"
setwd(ancestry_dir)

#### read in the eigenvectors, produced through ancestry estimation pipeline ####
eigenvec_name <- paste(study_name, ref_name, "eigenvec", sep = ".")
eigenvec <- read.table(eigenvec_name, header = TRUE, sep = '\t', comment.char = "")
names(eigenvec)[1] <- "FID"
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]

# get 1KG sample IDs and Population information
pedfile <- paste(dir, "shared/refGenomes/1000genomes/chrs/20130606_g1k.ped", sep = "/")
PED <- read.table(pedfile, header = TRUE, sep = '\t', comment.char = "")
PED <- select(PED, Individual.ID, Population)

# set Population as "study" for study sample IDs
study_tempdf <- data.frame(Individual.ID = setdiff(rownames(eigenvec), PED$Individual.ID), Population = study)

# combine 1KG and study sample
KG_study <- rbind(PED, study_tempdf)

# make sure all sample IDs are from the eigenvec file
KG_study <- KG_study[which(KG_study$Individual.ID %in% rownames(eigenvec)),]
KG_study <- KG_study[match(rownames(eigenvec), KG_study$Individual.ID),]
all(KG_study$Individual.ID == rownames(eigenvec)) == TRUE

# make KG_study$Population as factor
KG_study$Population <- factor(KG_study$Population, levels=c(
    study,
    "ACB","ASW","ESN","GWD","LWK","MSL","YRI", # AFR
    "CLM","MXL","PEL","PUR", # AMR
    "CDX","CHB","CHS","JPT","KHV", # EAS
    "CEU","FIN","GBR","IBS","TSI", # EUR
    "BEB","GIH","ITU","PJL","STU")) # SAS

# set colors for each population
col <- colorRampPalette(c(
    "purple",
    "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
    "forestgreen","forestgreen","forestgreen","forestgreen",
    "grey","grey","grey","grey","grey",
    "royalblue","royalblue","royalblue","royalblue","royalblue",
    "black","black","black","black","black"))(length(unique(KG_study$Population)))[factor(KG_study$Population)]

# set color for study samples
col_rgb <-  rgb(matrix(col2rgb("purple")/255, ncol=3), alpha = 0.2, maxColorValue = 1)

# generate PCA bi-plots
par(mfrow = c(1,2))
project.pca <- eigenvec

#. PC1 vs. PC2 ----
study_samples <- rownames(project.pca) %in% study_tempdf$Individual.ID
plot(project.pca[,1], project.pca[,2],
     type = 'n',
     main = paste(study, "and 1KG PC1 vs. PC2"),
     adj = 0.5,
     xlab = 'PC1',
     ylab = 'PC2',
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(project.pca[study_samples,1], project.pca[study_samples,2], col = col_rgb, pch = 20, cex = 1.2)
points(project.pca[!study_samples,1], project.pca[!study_samples,2], col = col, pch = 20, cex = 1.2)
legend('topleft',
       bty = 'n',
       title = '',
       c(study, 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
       fill = c('purple', 'yellow', 'forestgreen', 'grey', 'royalblue', 'black'))


#. PC1 vs. PC3 ----
plot(project.pca[,1], project.pca[,3],
     type = 'n',
     main = paste(study, "and 1KG PC1 vs. PC3"),
     adj = 0.5,
     xlab = 'PC1',
     ylab = 'PC3',
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(project.pca[study_samples,1], project.pca[study_samples,3], col = col_rgb, pch = 20, cex = 1.2)
points(project.pca[!study_samples,1], project.pca[!study_samples,3], col = col, pch = 20, cex = 1.2)
legend('topleft',
       bty = 'n',
       title = '',
       c(study, 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'),
       fill = c('purple', 'yellow', 'forestgreen', 'grey', 'royalblue', 'black'))



#### k-nearest neighbor (kNN) ####
file <- paste(study_name, ref_name, sep = ".")

# prepare a data frame with PCs and population info
pca <- read.table(paste0(file, ".eigenvec"), header = TRUE, sep = "\t", comment.char = ""); names(pca)[1] <- "FID"
pop <- read.table(paste0(file, ".pop"), header = TRUE)
df <- merge(pca, pop)
df <- data.frame(df, SuperFam = df$POP)
df$POP <- factor(df$POP, levels = c(
    "SET",
    "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI",
    "CLM", "MXL", "PEL", "PUR",
    "CDX", "CHB", "CHS", "JPT", "KHV",
    "CEU", "FIN", "GBR", "IBS", "TSI",
    "BEB", "GIH", "ITU", "PJL", "STU"
))
SF <- factor(c(
    study,
    "AFR", "AFR", "AFR", "AFR", "AFR", "AFR", "AFR",
    "AMR", "AMR", "AMR", "AMR",
    "EAS", "EAS", "EAS", "EAS", "EAS",
    "EUR", "EUR", "EUR", "EUR", "EUR",
    "SAS", "SAS", "SAS", "SAS", "SAS"
))

# get Super population labels
for (i in 1:length(levels(df$POP))) {
    df$SuperFam <- gsub(levels(df$POP)[i], SF[i], df$SuperFam)
}
df$SuperFam <- factor(df$SuperFam, levels = c(study, "AFR", "AMR", "EAS", "EUR", "SAS"))

# apply kNN
pr <- knn(df[df$POP != "SET", 3:7], df[df$POP == "SET", 3:7], cl = df$SuperFam[df$POP != "SET"], k = 13)
df$SuperFam[df$POP == "SET"] <- pr
df <- df[df$POP == "SET",]

# create ID file for each population
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
for (pop in pops) {
    df_pop <- df[df$SuperFam == pop,]
    num <- nrow(df_pop)
    write.table(df_pop, paste(study, "ancestry_knn", pop, num, "txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

