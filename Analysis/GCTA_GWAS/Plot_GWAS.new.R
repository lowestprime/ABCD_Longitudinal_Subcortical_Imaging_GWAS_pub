#### Comprehensive GWAS Analysis and Visualization Script ####
# Date: 9/7/2024

# This script allows you to generate GWAS plots manually, one at a time,
# by running individual code blocks in RStudio.

# -----------------------------
# Load Required Packages
# -----------------------------
# Ensure required installers are available
for (p in c("devtools", "BiocManager")) if (!require(p, quietly = TRUE)) install.packages(p)

# List of packages and their installation sources
packages <- list(
  "manhplot"    = "cgrace1978/manhplot",
  "hudson"      = "anastasia-lucas/hudson",
  "locuszoomr"  = "myles-lewis/locuszoomr",
  "fastman"     = "slowkow/fastman",
  "ggmanh"      = "bioc",
  "biovizBase"  = "github",  # biovizBase needs to be installed from GitHub
  "ggbio"       = "bioc"     # ggbio is a Bioconductor package
)

# Install and load packages
for (pkg in names(packages)) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (packages[[pkg]] == "bioc") {
      # Install development version of Bioconductor packages
      BiocManager::install(pkg, version = "devel")
    } else if (packages[[pkg]] == "github") {
      # Special case for biovizBase
      devtools::install_github("lawremi/biovizBase")
    } else {
      # Install GitHub packages
      devtools::install_github(packages[[pkg]], dependencies = TRUE, force = TRUE)
    }
  }
}

# Load Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  ggplot2, data.table, qqman, CMplot, pheatmap, ComplexHeatmap, devtools,
  TrumpetPlots, ensembldb, AnnotationFilter, GenomicRanges, biovizBase, ggbio,
  dplyr, manhplot, hudson, locuszoomr, fastman, ggmanh
)

# -----------------------------
# List .mlma Files Organized by Directory
# -----------------------------
# Define Parent Directory Containing .mlma Files
mlma_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/'

# Check if Directory Exists
if (!dir.exists(mlma_dir)) stop("The specified mlma_dir does not exist.")

# Recursively List All .mlma Files in Parent Directory
mlma_files <- list.files(mlma_dir, pattern = "\\.mlma$", full.names = TRUE, recursive = TRUE)

if (length(mlma_files) == 0) stop("No .mlma files found in the specified directory.")

# Process the file paths to extract directory components
paths_split <- strsplit(mlma_files, split = "/")

# Create a data frame with the directory components
df <- data.frame(
  full_path = mlma_files,
  stringsAsFactors = FALSE
)

# Extract last three components of the path
df$last_three <- sapply(paths_split, function(x) paste(x[(length(x)-2):length(x)], collapse = "/"))
components <- strsplit(df$last_three, split = "/")
df$dir1 <- sapply(components, function(x) x[1])
df$dir2 <- sapply(components, function(x) x[2])
df$file <- sapply(components, function(x) x[3])

# Arrange the data frame
df <- df %>% arrange(dir1, dir2, file)

# Group by directories and number the files
df <- df %>% group_by(dir1, dir2) %>% mutate(file_num = row_number())

# Print the organized list
prev_dir <- ""
dir_count <- 0
file_count <- 0

for (i in seq_len(nrow(df))) {
  dir <- paste0(df$dir1[i], "/", df$dir2[i])
  file_num <- df$file_num[i]
  file_name <- df$file[i]
  
  if (dir != prev_dir) {
    cat(dir, "\n")
    prev_dir <- dir
    dir_count <- dir_count + 1
  }
  cat(file_num, ". ", file_name, "\n", sep = "")
  file_count <- file_count + 1
}

# Print total counts
cat(dir_count, "dirs, ", file_count, " files\n", sep = "")

# -----------------------------
# Helper Functions
# -----------------------------
# Function to Process GWAS Results
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

# -----------------------------
# Load and Process GWAS Results from a Specific File
# -----------------------------
# Select a specific .mlma file to process
# For example, choose the first file in the list
file_to_process <- mlma_files[1]

# Alternatively, specify the file path directly
# file_to_process <- "/path/to/your/file.mlma"

# Process the selected GWAS results file
gwasResults <- process_gwas_results(file_to_process)

# Check if data is loaded correctly
if (nrow(gwasResults) == 0) {
  warning("No valid data in file:", file_to_process)
} else {
  cat("Data loaded successfully for file:", file_to_process, "\n")
}

# Define Plot Directory and Create It If It Doesn't Exist
file_dir <- dirname(file_to_process)
plot_dir <- file.path(file_dir, 'plots')
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Define Plot Prefix
plot_prefix <- sub("\\.mlma$", "", basename(file_to_process))

# -----------------------------
# Generate Plots Manually
# -----------------------------
# Now you can generate plots manually by running individual code blocks below.

# --- Example: Generate QQ Plot ---
# Run this block to generate a QQ plot
svg(file.path(plot_dir, paste0(plot_prefix, "_qq_plot.svg")))
qq(gwasResults$P, main = paste0("QQ Plot for ", plot_prefix))
dev.off()

# --- Example: Generate Manhattan Plot using qqman ---
# Run this block to generate a Manhattan plot
svg(file.path(plot_dir, paste0(plot_prefix, "_qqman_manhattan.svg")))
manhattan(
  gwasResults, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  col = c("blue4", "orange3"), genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5), chrlabs = as.character(1:22),
  main = paste0("Manhattan Plot for ", plot_prefix)
)
dev.off()

# --- Example: Generate CMplot Circular Manhattan Plot ---
# Run this block to generate a circular Manhattan plot
CMplot(
  gwasResults, plot.type = "c", LOG10 = TRUE, threshold = 5e-8,
  highlight = gwasResults$SNP[gwasResults$P < 5e-8],
  highlight.text = gwasResults$SNP[gwasResults$P < 5e-8],
  file.output = TRUE, file = "svg", memo = plot_prefix,
  dpi = 300, file.name = file.path(plot_dir, paste0(plot_prefix, "_CMplot")), verbose = FALSE
)

# --- Example: Generate Trumpet Plot ---
# Ensure required columns are present
if (all(c("SNP", "Freq", "Beta") %in% names(gwasResults))) {
  # Prepare dataset for TrumpetPlots
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
  warning("Required columns (SNP, Freq, Beta) not found. Skipping Trumpet Plot.")
}

# --- Example: Generate LocusZoom Plot using locuszoomr ---
# Run this block to generate a LocusZoom plot
top_snp <- gwasResults[which.min(P), ]
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

# --- Example: Generate Manhattan++ Plot using manhplot ---
# Run this block to generate a Manhattan++ plot
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

# --- Example: Generate Hudson Plot ---
# Note: Hudson plots require two datasets
# Select a second file to compare
file_to_compare <- mlma_files[2]  # For example, the second file in the list

# Process the second GWAS results file
gwasResults2 <- process_gwas_results(file_to_compare)

# Check if data is loaded correctly
if (nrow(gwasResults2) == 0) {
  warning("No valid data in file:", file_to_compare)
} else {
  cat("Data loaded successfully for file:", file_to_compare, "\n")
  
  # Generate Hudson Plot
  hudson_manhattan(
    data1 = gwasResults,
    data2 = gwasResults2,
    chr_col1 = "CHR",
    pos_col1 = "BP",
    p_col1 = "P",
    chr_col2 = "CHR",
    pos_col2 = "BP",
    p_col2 = "P",
    highlight_p = 5e-8,
    file_name = file.path(plot_dir, paste0(plot_prefix, "_vs_", sub("\\.mlma$", "", basename(file_to_compare)), "_Hudson.svg")),
    file_format = "svg"
  )
}

# --- Additional Plotting Functions ---
# You can include additional plotting functions here and run them as needed.

# -----------------------------
# Notes:
# -----------------------------
# - To process a different .mlma file, change the 'file_to_process' variable to point to your desired file.
# - Ensure that the necessary columns are present in your data before running a plot function.
# - You can run each code block individually in RStudio to generate the corresponding plot.
# - The generated plots will be saved in the 'plots' directory within the directory of your .mlma file.

# End of Script
