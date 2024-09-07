#### GCTA GWAS PLOTTING ####
# 9/7/2024

# Load Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ggplot2, fastman, manhplot, hudson, ggmanh, devtools)

# Define Directories
mlma_dir <- '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run'
plot_dir <- '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/plots'