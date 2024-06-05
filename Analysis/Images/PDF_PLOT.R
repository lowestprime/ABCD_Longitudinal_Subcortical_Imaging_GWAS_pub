# Load necessary libraries
library(pdftools)
library(magick)

# Define the path to the PDF file
pdf_path <- "path/to/your/ROC_sMRI_Vol_SCS_All.pdf"

# Extract images from the PDF
pdf_info <- pdf_info(pdf_path)
pdf_images <- pdf_convert(pdf_path, format = 'png', pages = 1:pdf_info$pages)

# Create a function to split images by region
split_image_by_region <- function(image_path, rows, cols) {
  image <- image_read(image_path)
  image_info <- image_info(image)
  
  # Calculate the width and height of each region
  region_width <- image_info$width / cols
  region_height <- image_info$height / rows
  
  regions <- list()
  for (i in 0:(rows-1)) {
    for (j in 0:(cols-1)) {
      region <- image_crop(image, geometry_area(region_width, region_height, j * region_width, i * region_height))
      regions <- append(regions, list(region))
    }
  }
  return(regions)
}

# Define the number of rows and columns to split the images into
rows <- 3
cols <- 3

# Split each extracted image into regions
all_regions <- lapply(pdf_images, split_image_by_region, rows, cols)

# Combine all regions into a summary image
combined_image <- image_blank(width = 2000, height = 2000, color = "white")
current_x <- 0
current_y <- 0
for (regions in all_regions) {
  for (region in regions) {
    combined_image <- image_composite(combined_image, region, offset = geometry_point(current_x, current_y))
    current_x <- current_x + region$width
    if (current_x >= combined_image$width) {
      current_x <- 0
      current_y <- current_y + region$height
    }
  }
}

# Save the combined summary image
image_write(combined_image, path = "summary_image.png")
