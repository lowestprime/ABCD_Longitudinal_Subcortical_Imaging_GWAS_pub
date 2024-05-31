# Function to add heading level
add_heading <- function(level) {
  library(rstudioapi)
  # Get the current document context
  context <- rstudioapi::getActiveDocumentContext()

  # Get the line number where the cursor is
  line_num <- context$selection[[1]]$range$start[[1]]

  # Get the current line text
  line <- context$contents[line_num]

  # Remove existing # symbols at the start and end
  line <- sub("^#+\\s*|\\s*#+$", "", line)

  # Create the new heading
  new_line <- paste0(strrep("#", level), " ", line, " ", strrep("#", level))

  # Define the range for the entire line
  range <- rstudioapi::document_range(
    rstudioapi::document_position(line_num, 1),
    rstudioapi::document_position(line_num, nchar(context$contents[line_num]) + 1)
  )

  # Modify the document
  rstudioapi::modifyRange(location = range, text = new_line)
}

# Wrapper functions for each heading level
add_heading_1 <- function() add_heading(1)
add_heading_2 <- function() add_heading(2)
add_heading_3 <- function() add_heading(3)
add_heading_4 <- function() add_heading(4)
add_heading_5 <- function() add_heading(5)
add_heading_6 <- function() add_heading(6)
