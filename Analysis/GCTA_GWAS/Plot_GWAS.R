#### GCTA GWAS PLOTTING ####
# 9/7/2024

# Load Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ggplot2, fastman, manhplot, hudson, ggmanh, CMplot, qqman, pheatmap, ComplexHeatmap, devtools, data.table, TrumpetPlots, locuszoomr)

# Define Directories
mlma_dir <- '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/test_run'
plot_dir <- '~/project-lhernand/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/plots'

# Ensure directories exist
if (!dir.exists(mlma_dir)) dir.create(mlma_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Process all .mlma files in mlma_dir
mlma_files <- list.files(mlma_dir, pattern = "\\.mlma$", full.names = TRUE)

# Helper function to rename and clean columns
process_gwas_results <- function(file) {
  gwasResults <- fread(file, header = TRUE)
  # Rename the necessary columns
  colnames(gwasResults)[which(names(gwasResults) == "Chr")] <- "CHR"
  colnames(gwasResults)[which(names(gwasResults) == "bp")] <- "BP"
  colnames(gwasResults)[which(names(gwasResults) == "p")] <- "P"
  
  # Remove NA and non-finite p-values
  gwasResults <- gwasResults[!is.na(gwasResults$P) & is.finite(gwasResults$P), ]
  return(gwasResults)
}

# List all .mlma files in the directory
mlma_files <- list.files(mlma_dir, pattern = "\\.mlma$", full.names = TRUE)

#### SECTION: CMplot Circular Manhattan Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  CMplot(gwasResults, plot.type = "c", LOG10 = TRUE, threshold = 5e-8,
         highlight = gwasResults$SNP[gwasResults$P < 5e-8],
         highlight.text = gwasResults$SNP[gwasResults$P < 5e-8],
         file.output = TRUE, file = "svg",
         file.name = file.path(plot_dir, paste0(plot_prefix, "_CMplot")))
}

#### SECTION: qqman Manhattan Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  svg(file.path(plot_dir, paste0(plot_prefix, "_qqman_manhattan.svg")))
  manhattan(gwasResults, chr = "CHR", bp = "BP", p = "P", col = c("blue4", "orange3"))
  dev.off()
}

#### SECTION: fastman Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  fastman_plot <- fastman(gwasResults, chr = "CHR", bp = "BP", p = "P", col = "matlab")
  ggsave(file.path(plot_dir, paste0(plot_prefix, "_fastman.svg")), plot = fastman_plot)
}

#### SECTION: Manhattan++ Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  manhplot(gwasResults, chr = "CHR", bp = "BP", p = "P",
           highlight = gwasResults$SNP[gwasResults$P < 5e-8],
           file = file.path(plot_dir, paste0(plot_prefix, "_ManhattanPlusPlus.svg")))
}

#### SECTION: Hudson Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  hudson_manhattan(gwasResults, chr = "CHR", bp = "BP", p = "P",
                   file_name = file.path(plot_dir, paste0(plot_prefix, "_Hudson.svg")))
}

#### SECTION: ggmanh Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  gwasData <- manhattan_data_preprocess(gwasResults, chr = "CHR", bp = "BP", p = "P")
  ggmanh_plot <- ggmanh(gwasData, chr = "CHR", bp = "BP", p = "P")
  ggsave(file.path(plot_dir, paste0(plot_prefix, "_ggmanh.svg")), plot = ggmanh_plot)
}

#### SECTION: Trumpet Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  TrumpetPlots(gwasResults, chr = "CHR", bp = "BP", p = "P", out = file.path(plot_dir, paste0(plot_prefix, "_TrumpetPlot.svg")))
}

#### SECTION: Modern QQ Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  svg(file.path(plot_dir, paste0(plot_prefix, "_qq_plot.svg")))
  qq(gwasResults$P)
  dev.off()
}

#### SECTION: Phenotype Distribution Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  svg(file.path(plot_dir, paste0(plot_prefix, "_PhenotypeDistribution.svg")))
  ggplot(gwasResults, aes(x = Genotype, y = BrainPhenotype)) + 
    geom_boxplot() +
    labs(title = "Brain Phenotype by Genotype", x = "Genotype", y = "Brain Phenotype")
  dev.off()
}

#### SECTION: pheatmap ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  # Example heatmap data - replace with relevant SNP effect data
  pheatmap(matrix(rnorm(100), nrow = 10),
           file = file.path(plot_dir, paste0(plot_prefix, "_pheatmap.svg")))
}

#### SECTION: ComplexHeatmap ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  ComplexHeatmap::Heatmap(matrix(rnorm(100), nrow = 10),
                          file = file.path(plot_dir, paste0(plot_prefix, "_ComplexHeatmap.svg")))
}

#### SECTION: LocusZoom Plot ####
for (file in mlma_files) {
  gwasResults <- process_gwas_results(file)
  plot_prefix <- gsub("\\.mlma$", "", basename(file))
  
  locuszoomr::locuszoom(gwasResults, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                        output = file.path(plot_dir, paste0(plot_prefix, "_LocusZoom.svg")))
}
