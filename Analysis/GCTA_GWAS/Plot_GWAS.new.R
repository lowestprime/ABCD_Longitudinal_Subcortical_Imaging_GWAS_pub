#### Comprehensive GWAS Analysis and Visualization Script ####
# Date: 9/7/2024

# Load Required Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  ggplot2, data.table, qqman, CMplot, pheatmap, ComplexHeatmap, devtools,
  TrumpetPlots, ensembldb, AnnotationFilter, GenomicRanges, biovizBase, ggbio
)

# Install and load packages from GitHub
if (!require("manhplot", quietly = TRUE)) {
  devtools::install_github("cgrace1978/manhplot", dependencies = TRUE, force = TRUE)
  library(manhplot)
}

if (!require("hudson", quietly = TRUE)) {
  devtools::install_github("anastasia-lucas/hudson")
  library(hudson)
}

if (!require("ggmanh", quietly = TRUE)) {
  BiocManager::install("ggmanh")
  library(ggmanh)
}

if (!require("locuszoomr", quietly = TRUE)) {
  devtools::install_github("myles-lewis/locuszoomr")
  library(locuszoomr)
}

if (!require("fastman", quietly = TRUE)) {
  devtools::install_github("slowkow/fastman")
  library(fastman)
}

# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ensembldb", "AnnotationFilter", "GenomicRanges", "biovizBase", "ggbio", "ComplexHeatmap"))

# Load additional libraries
library(GenomicRanges)
library(biovizBase)
library(ggbio)

# Define Parent Directory Containing .mlma Files
mlma_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/'

# Check if Directory Exists
if (!dir.exists(mlma_dir)) stop("The specified mlma_dir does not exist.")

# Recursively List All .mlma Files in Parent Directory
mlma_files <- list.files(mlma_dir, pattern = "\\.mlma$", full.names = TRUE, recursive = TRUE)

if (length(mlma_files) == 0) stop("No .mlma files found in the specified directory.")

# Helper Function to Process GWAS Results
process_gwas_results <- function(file) {
  gwasResults <- fread(file, header = TRUE)
  
  # Rename Necessary Columns
  setnames(gwasResults, old = c("Chr", "bp", "p", "SNP", "beta", "se", "freq"), 
           new = c("CHR", "BP", "P", "SNP", "Beta", "SE", "Freq"), skip_absent = TRUE)
  
  # Remove NA and Non-finite P-values
  gwasResults <- gwasResults[!is.na(P) & is.finite(P)]
  
  # Ensure CHR and BP are Numeric
  gwasResults[, CHR := as.numeric(CHR)]
  gwasResults[, BP := as.numeric(BP)]
  
  # Remove SNPs with Missing CHR or BP
  gwasResults <- gwasResults[!is.na(CHR) & !is.na(BP)]
  
  # Keep Only Autosomal Chromosomes 1-22
  gwasResults <- gwasResults[CHR %in% 1:22]
  
  return(gwasResults)
}

# Function to Generate CMplot Circular Manhattan Plot
generate_cmplot <- function(gwasResults, plot_dir, plot_prefix) {
  CMplot(
    gwasResults, plot.type = "c", LOG10 = TRUE, threshold = 5e-8,
    highlight = gwasResults$SNP[gwasResults$P < 5e-8],
    highlight.text = gwasResults$SNP[gwasResults$P < 5e-8],
    file.output = TRUE, file = "svg", memo = plot_prefix,
    dpi = 300, file.name = file.path(plot_dir, paste0(plot_prefix, "_CMplot")), verbose = FALSE
  )
}

# Function to Generate qqman Manhattan Plot
generate_qqman_manhattan <- function(gwasResults, plot_dir, plot_prefix) {
  svg(file.path(plot_dir, paste0(plot_prefix, "_qqman_manhattan.svg")))
  manhattan(
    gwasResults, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
    col = c("blue4", "orange3"), genomewideline = -log10(5e-8),
    suggestiveline = -log10(1e-5), chrlabs = as.character(1:22)
  )
  dev.off()
}

# Function to Generate fastman Plot
generate_fastman_plot <- function(gwasResults, plot_dir, plot_prefix) {
  fastman_plot <- fastman(gwasResults, chr = "CHR", bp = "BP", p = "P", col = "matlab")
  ggsave(file.path(plot_dir, paste0(plot_prefix, "_fastman.svg")), plot = fastman_plot, width = 10, height = 6)
}

# Function to Generate Manhattan++ Plot using manhplot package
generate_manhplot <- function(gwasResults, plot_dir, plot_prefix) {
  # Generate Manhattan++ Plot
  manhplot(
    data = gwasResults,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    annotate_p = 5e-8,
    annotate_top = TRUE,
    file_name = file.path(plot_dir, paste0(plot_prefix, "_ManhattanPlusPlus.svg")),
    file_type = "svg"
  )
}

# Function to Generate ggmanh Plot
generate_ggmanh_plot <- function(gwasResults, plot_dir, plot_prefix) {
  # Adjust column names if necessary
  gwasData <- copy(gwasResults)
  setnames(gwasData, old = c("CHR", "BP", "P"), new = c("chrom", "pos", "pvalue"), skip_absent = TRUE)
  
  # Generate ggmanh Plot
  ggmanh_plot <- ggmanh(
    data = gwasData,
    chrom = "chrom",
    pos = "pos",
    pvalue = "pvalue"
  )
  
  # Save plot
  ggsave(file.path(plot_dir, paste0(plot_prefix, "_ggmanh.svg")), plot = ggmanh_plot, width = 10, height = 6)
}

# Function to Generate Trumpet Plot
generate_trumpet_plot <- function(gwasResults, plot_dir, plot_prefix) {
  # Ensure required columns are present
  required_cols <- c("SNP", "Freq", "Beta")
  if (all(required_cols %in% names(gwasResults))) {
    # Prepare dataset for TrumpetPlots
    # Adding placeholder columns if necessary
    if (!"Analysis" %in% names(gwasResults)) gwasResults$Analysis <- "GWAS"
    if (!"Gene" %in% names(gwasResults)) gwasResults$Gene <- ""
    
    # Generate Trumpet Plot
    svg(file.path(plot_dir, paste0(plot_prefix, "_TrumpetPlot.svg")), width = 8, height = 6)
    plot_trumpets(
      dataset = gwasResults,
      rsID = "SNP",
      freq = "Freq",
      A1_beta = "Beta",
      Analysis = "Analysis",
      Gene = "Gene",
      calculate_power = TRUE,
      show_power_curves = TRUE,
      threshold = c(0.7, 0.9),
      N = max(gwasResults$N, na.rm = TRUE),  # Assuming 'N' column is available
      alpha = 5e-08,
      Nfreq = 500,
      power_color_palette = c("purple", "deeppink"),
      analysis_color_palette = c("#018571", "#a6611a")
    )
    dev.off()
  } else {
    warning("Required columns (SNP, Freq, Beta) not found in file:", plot_prefix, "Skipping Trumpet Plot.")
  }
}

# Function to Generate Modern QQ Plot
generate_qq_plot <- function(gwasResults, plot_dir, plot_prefix) {
  svg(file.path(plot_dir, paste0(plot_prefix, "_qq_plot.svg")))
  qq(gwasResults$P, main = paste0("QQ Plot for ", plot_prefix))
  dev.off()
}

# Function to Generate Phenotype Distribution Plot
generate_phenotype_distribution_plot <- function(gwasResults, plot_dir, plot_prefix) {
  if ("Genotype" %in% names(gwasResults) & "BrainPhenotype" %in% names(gwasResults)) {
    svg(file.path(plot_dir, paste0(plot_prefix, "_PhenotypeDistribution.svg")))
    ggplot(gwasResults, aes(x = as.factor(Genotype), y = BrainPhenotype)) + 
      geom_boxplot() +
      labs(title = "Brain Phenotype by Genotype", x = "Genotype", y = "Brain Phenotype") +
      theme_minimal()
    dev.off()
  } else {
    warning("Columns 'Genotype' and 'BrainPhenotype' not found in file:", plot_prefix, "Skipping Phenotype Distribution Plot.")
  }
}

# Function to Generate pheatmap Heatmap
generate_pheatmap <- function(gwasResults, plot_dir, plot_prefix) {
  # Example heatmap data - replace with relevant SNP effect data if available
  heatmap_data <- matrix(rnorm(100), nrow = 10)
  rownames(heatmap_data) <- paste0("Gene", 1:10)
  colnames(heatmap_data) <- paste0("Sample", 1:10)
  
  svg(file.path(plot_dir, paste0(plot_prefix, "_pheatmap.svg")), width = 8, height = 6)
  pheatmap(heatmap_data, main = paste0("Heatmap for ", plot_prefix))
  dev.off()
}

# Function to Generate ComplexHeatmap
generate_complex_heatmap <- function(gwasResults, plot_dir, plot_prefix) {
  # Example heatmap data - replace with relevant data if available
  heatmap_data <- matrix(rnorm(100), nrow = 10)
  rownames(heatmap_data) <- paste0("Gene", 1:10)
  colnames(heatmap_data) <- paste0("Sample", 1:10)
  
  svg(file.path(plot_dir, paste0(plot_prefix, "_ComplexHeatmap.svg")), width = 8, height = 6)
  Heatmap(heatmap_data,
          name = "Expression",
          show_row_names = TRUE,
          show_column_names = TRUE,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          column_title = paste0("Complex Heatmap for ", plot_prefix),
          row_title = "Genes")
  dev.off()
}

# Function to Generate LocusZoom-like Regional Association Plot
generate_regional_plot <- function(gwasResults, plot_dir, plot_prefix) {
  # Identify the top SNP
  top_snp <- gwasResults[which.min(P), ]
  
  # Define the region around the top SNP (±500kb)
  region_start <- max(0, top_snp$BP - 500000)
  region_end <- top_snp$BP + 500000
  region_chr <- top_snp$CHR
  
  # Subset data for the region
  region_data <- gwasResults[CHR == region_chr & BP >= region_start & BP <= region_end]
  
  if (nrow(region_data) > 0) {
    # Create GRanges object for plotting
    gr <- GRanges(
      seqnames = Rle(paste0("chr", region_data$CHR)),
      ranges = IRanges(start = region_data$BP, width = 1),
      SNP = region_data$SNP,
      P = region_data$P
    )
    
    # Fetch gene annotations for the region
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    genes <- genes(txdb, columns = c("gene_id"), filter = GRanges(seqnames = paste0("chr", region_chr), 
                                                                  ranges = IRanges(region_start, region_end)))
    
    # Plot regional association with gene annotations
    svg(file.path(plot_dir, paste0(plot_prefix, "_RegionalPlot.svg")), width = 10, height = 6)
    
    # Association plot
    p1 <- ggplot(as.data.frame(gr), aes(x = start, y = -log10(P))) +
      geom_point(color = "blue", size = 1.5) +
      theme_bw() +
      labs(title = paste0("Regional Association Plot for ", top_snp$SNP),
           x = paste0("Chromosome ", region_chr, " Position (bp)"),
           y = "-log10(P-value)") +
      geom_vline(xintercept = top_snp$BP, linetype = "dashed", color = "red") +
      geom_text_repel(aes(label = ifelse(-log10(P) > 7, as.character(SNP), "")), size = 3)
    
    # Gene track
    p2 <- autoplot(genes, fill = "lightblue") +
      theme_bw() +
      labs(y = "Genes") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    # Combine plots
    gridExtra::grid.arrange(p1, p2, heights = c(3, 1), ncol = 1)
    
    dev.off()
  } else {
    warning("No data available for regional plot in file:", plot_prefix)
  }
}

# Function to Generate LocusZoom Plot using locuszoomr
generate_locuszoom_plot <- function(gwasResults, plot_dir, plot_prefix) {
  # Identify the top SNP
  top_snp <- gwasResults[which.min(P), ]
  
  # Generate LocusZoom Plot
  locuszoom(
    data = gwasResults,
    snp = top_snp$SNP,
    build = "hg19",  # or 'hg38' if your data uses hg38 coordinates
    pval_col = "P",
    chr_col = "CHR",
    pos_col = "BP",
    flank = 500000,  # 500kb on each side
    out = file.path(plot_dir, paste0(plot_prefix, "_LocusZoom"))
  )
}

# Function to Generate Hudson Plot (requires two datasets)
generate_hudson_plot <- function(gwasResults1, gwasResults2, plot_dir, plot_prefix) {
  # Ensure both datasets have the necessary columns
  required_cols <- c("CHR", "BP", "P")
  if (all(required_cols %in% names(gwasResults1)) & all(required_cols %in% names(gwasResults2))) {
    # Generate Hudson Plot
    hudson_manhattan(
      data1 = gwasResults1,
      data2 = gwasResults2,
      chr_col1 = "CHR",
      pos_col1 = "BP",
      p_col1 = "P",
      chr_col2 = "CHR",
      pos_col2 = "BP",
      p_col2 = "P",
      highlight_p = 5e-8,
      file_name = file.path(plot_dir, paste0(plot_prefix, "_Hudson.svg")),
      file_format = "svg"
    )
  } else {
    warning("Required columns not found in one or both datasets. Skipping Hudson Plot.")
  }
}

#### Main Processing Loop ####
# For Hudson plots, we need to specify pairs of files to compare
# Example: hudson_pairs <- list(c(file1, file2), c(file3, file4))
hudson_pairs <- list()  # Define your pairs here

for (file in mlma_files) {
  cat("Processing file:", file, "\n")
  
  # Process GWAS Results
  gwasResults <- process_gwas_results(file)
  
  if (nrow(gwasResults) == 0) {
    warning("No valid data in file:", file)
    next
  }
  
  # Define Plot Directory and Create It If It Doesn't Exist
  file_dir <- dirname(file)
  plot_dir <- file.path(file_dir, 'plots')
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define Plot Prefix
  plot_prefix <- sub("\\.mlma$", "", basename(file))
  
  # Generate CMplot Circular Manhattan Plot
  generate_cmplot(gwasResults, plot_dir, plot_prefix)
  
  # Generate qqman Manhattan Plot
  generate_qqman_manhattan(gwasResults, plot_dir, plot_prefix)
  
  # Generate fastman Plot
  generate_fastman_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate Manhattan++ Plot using manhplot
  generate_manhplot(gwasResults, plot_dir, plot_prefix)
  
  # Generate ggmanh Plot
  generate_ggmanh_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate Trumpet Plot
  generate_trumpet_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate Modern QQ Plot
  generate_qq_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate Phenotype Distribution Plot
  generate_phenotype_distribution_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate pheatmap Heatmap
  generate_pheatmap(gwasResults, plot_dir, plot_prefix)
  
  # Generate ComplexHeatmap
  generate_complex_heatmap(gwasResults, plot_dir, plot_prefix)
  
  # Generate LocusZoom-like Regional Association Plot
  generate_regional_plot(gwasResults, plot_dir, plot_prefix)
  
  # Generate LocusZoom Plot using locuszoomr
  generate_locuszoom_plot(gwasResults, plot_dir, plot_prefix)
}

# Generate Hudson Plots for specified pairs
for (pair in hudson_pairs) {
  file1 <- pair[1]
  file2 <- pair[2]
  cat("Generating Hudson Plot for files:", file1, "and", file2, "\n")
  
  # Process GWAS Results for both files
  gwasResults1 <- process_gwas_results(file1)
  gwasResults2 <- process_gwas_results(file2)
  
  if (nrow(gwasResults1) == 0 | nrow(gwasResults2) == 0) {
    warning("One or both files have no valid data. Skipping Hudson Plot for this pair.")
    next
  }
  
  # Define Plot Directory and Create It If It Doesn't Exist
  file_dir <- dirname(file1)
  plot_dir <- file.path(file_dir, 'plots')
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define Plot Prefix
  plot_prefix <- paste0(sub("\\.mlma$", "", basename(file1)), "_vs_", sub("\\.mlma$", "", basename(file2)))
  
  # Generate Hudson Plot
  generate_hudson_plot(gwasResults1, gwasResults2, plot_dir, plot_prefix)
}