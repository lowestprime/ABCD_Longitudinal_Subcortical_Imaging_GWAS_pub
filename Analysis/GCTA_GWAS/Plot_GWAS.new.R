#### Comprehensive GWAS Analysis and Visualization Script ####
# Date: 9/7/2024

# This script allows you to generate GWAS plots manually, one at a time,
# by running individual code blocks in RStudio.

# -----------------------------
# Load Required Packages
# -----------------------------
# Ensure required installers are available
for (p in c("devtools", "BiocManager", "conflicted")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

# Resolve conflicts
conflicts_prefer(
  dplyr::filter, dplyr::select, dplyr::rename, dplyr::mutate, data.table::first,BiocGenerics::combine,
  dplyr::recode, dplyr::slice, dplyr::setdiff, fs::path, BiocGenerics::combine, ggplot2::stat_qq_line
)


# List of packages and their installation sources
packages <- list(
  "manhplot"    = "cgrace1978/manhplot",
  "hudson"      = "anastasia-lucas/hudson",
  "locuszoomr"  = "myles-lewis/locuszoomr",
  "fastman"     = "kaustubhad/fastman",
  "ggrepel"     = "slowkow/ggrepel",
  "ggmanh"      = "bioc",
  "biovizBase"  = "github",
  "ggbio"       = "bioc",
  "scattermore" = "exaexa/scattermore",
  "ggfastman"   = "roman-tremmel/ggfastman",
  "fastqq"      = "gumeo/fastqq"
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
      devtools::install_github(packages[[pkg]], dependencies = T, force = T)
    }
  }
}

# Load Packages
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  ggplot2, data.table, qqplotr, CMplot, pheatmap, ComplexHeatmap, devtools,
  TrumpetPlots, ensembldb, AnnotationFilter, GenomicRanges, biovizBase, ggbio, svglite,
  manhplot, hudson, locuszoomr, fastman, fs, ggmanh, ggrepel, tidyr, dplyr, purrr,
  ggfastman, ggforce, ggrastr, scattermore, fastqq, hexbin
)


# -----------------------------
# List .mlma Files Organized by Directory
# -----------------------------
# Define Parent Directory Containing .mlma Files
mlma_dir <- '/u/project/lhernand/cobeaman/ABCD_Longitudinal_Subcortical_Imaging_GWAS/Analysis/GCTA_GWAS/Processed_Data/Results/'

# Load and process files
mlma_files <- dir_ls(mlma_dir, glob = "*.mlma", recurse = TRUE)  # Assign to mlma_files
if (!length(mlma_files)) stop("No .mlma files found.")

# # Process file paths into a tibble with directory and file info
# df <- tibble(full_path = mlma_files) %>%
#   mutate(path = map_chr(full_path, ~ paste(tail(strsplit(.x, "/")[[1]], 3), collapse = "/"))) %>%
#   separate(path, into = c("dir1", "dir2", "file"), sep = "/") %>%
#   arrange(dir1, dir2, file) %>%
#   group_by(dir1, dir2) %>%
#   mutate(file_num = row_number()) %>%
#   ungroup()
# 
# # Group by directories and collapse file names, print each directory once
# df %>%
#   group_by(dir1, dir2) %>%
#   reframe(
#     dir_path = paste(dir1, dir2, sep = "/"),
#     files = paste0(file_num, ". ", file, collapse = "\n")
#   ) %>%
#   distinct(dir_path, files) %>%  # Ensure distinct groups
#   pwalk(function(dir_path, files) {
#     cat(dir_path, "\n", files, "\n\n")
#   })
# 
# # Print the final summary of total directories and files
# cat("\n", n_distinct(df$dir1, df$dir2), "dirs,", nrow(df), "files\n")

# -----------------------------
# Helper Functions
# -----------------------------
# Function to Process GWAS Results
process_gwas_results <- function(file) {
  # Read the data and perform the operations step by step
  gwasResults <- fread(file, header = TRUE)[
    # Rename columns where necessary
    , setnames(.SD, old = c("Chr", "bp", "p", "SNP", "beta", "se", "freq"), 
               new = c("CHR", "BP", "P", "SNP", "Beta", "SE", "Freq"), 
               skip_absent = TRUE)][
                 # Filter for valid P values and numeric CHR and BP columns
                 !is.na(P) & is.finite(P) & !is.na(CHR) & !is.na(BP) & CHR %in% 1:22
               ]
  
  # Convert CHR and BP to numeric if not already
  gwasResults[, `:=`(CHR = as.numeric(CHR), BP = as.numeric(BP))]
  
  return(gwasResults)
}

# -----------------------------
# Load and Process GWAS Results from a Specific File
# -----------------------------
# Select a specific .mlma file to process
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
  
  # Create plot directory and prefix
  plot_dir <- file.path(dirname(file_to_process), 'plots')
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_prefix <- sub("\\.mlma$", "", basename(file_to_process))
}

# -----------------------------
# Generate Plots
# -----------------------------

# --- Generate Optimized QQ Plot with ggfastman ---
generate_qq_plots <- function(gwasResults, plot_dir, plot_prefix) {
  # Prepare data for ggfastman
  data_fastman <- gwasResults
  colnames(data_fastman)[colnames(data_fastman) == "P"] <- "pval"
  
  # Filter out invalid p-values and set a reasonable lower threshold
  data_fastman <- data_fastman[data_fastman$pval > 1e-300 & data_fastman$pval < 1, ]
  
  # Prioritize the most significant p-values by keeping the smallest ones
  # Retain the top 5% smallest p-values and downsample the remaining data
  significant_threshold <- quantile(data_fastman$pval, 0.05)  # Top 5% of p-values
  significant_points <- data_fastman[data_fastman$pval <= significant_threshold, ]
  non_significant_points <- data_fastman[data_fastman$pval > significant_threshold, ]
  
  # Downsample the non-significant points to 50,000
  if (nrow(non_significant_points) > 50000) {
    set.seed(123)  # Ensure reproducibility
    non_significant_points <- non_significant_points[sample(1:nrow(non_significant_points), 50000), ]
  }
  
  # Combine significant and downsampled non-significant points
  final_data <- rbind(significant_points, non_significant_points)
  
  # Generate QQ plot using ggfastman
  qq_plot <- ggfastman::fast_qq(
    final_data$pval,  # Use filtered and downsampled p-values
    speed = "slow",   # Slow mode to keep points vectorized
    pointsize = 1.5,  # Adjust point size for better visibility
    linecolor = "deeppink",
    ci_color = "steelblue",
    ci_alpha = 0.3,
    log10 = TRUE,
    inflation_method = "median",
    title = paste0("QQ Plot - ", plot_prefix)
  ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save the plot as an optimized SVG
  svg_file_path <- file.path(plot_dir, paste0(plot_prefix, "_qqplot_significant.svg"))
  svglite(svg_file_path, width = 7, height = 7)
  print(qq_plot)
  dev.off()
  
  cat("QQ plot with prioritized significant points saved successfully as SVG in", plot_dir, "\n")
}

# Generate QQ plots
generate_qq_plots(gwasResults, plot_dir, plot_prefix)

# --- Generate Manhattan Plot using ggfastman::fast_manhattan ---
# Run this block to generate Manhattan plot
svglite(file.path(plot_dir, paste0(plot_prefix, "_ggfastman_manhattan.svg")),
        system_fonts = list(sans = "Arial"),
        fix_text_size = FALSE)

# Prepare data
df <- gwasResults %>%
  filter(!is.na(P), P > 0, P <= 1) %>%
  mutate(
    logP = -log10(P),
    chr = as.factor(CHR),
    pos = BP
  )

# Define significance thresholds
genomewide_line <- -log10(5e-8)
suggestive_line <- -log10(1e-5)

# Identify top SNPs for annotation
top_snps <- df %>%
  group_by(chr) %>%
  top_n(-1, P) %>%  # SNP with lowest P-value in each chromosome
  ungroup() %>%
  filter(P <= 5e-8)  # Genome-wide significant SNPs

# Generate Manhattan plot with annotations
manhattan_plot <- ggfastman::fast_manhattan(df,
                                            build = 'hg19',  # Adjust genome build if needed
                                            col1 = "deepskyblue3",
                                            col2 = "coral2",
                                            ymax = max(df$logP, genomewide_line + 1),
                                            title = paste0("Manhattan Plot - ", plot_prefix)) +
  geom_hline(yintercept = genomewide_line, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = suggestive_line, linetype = "dashed", color = "grey") +
  ggrepel::geom_text_repel(data = top_snps,
                           aes(x = pos, y = logP, label = SNP),
                           size = 3,
                           color = "black",
                           max.overlaps = 10) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(manhattan_plot)
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
