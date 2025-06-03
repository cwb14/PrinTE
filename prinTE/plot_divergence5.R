#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(viridis)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide input files as command line arguments.", call. = FALSE)
}

# Parse optional output ID
output_id <- NULL
if (length(args) > 1) {
  output_id <- args[length(args)]
  input_files <- args[1:(length(args)-1)]
} else {
  input_files <- args
}

# Determine PDF filename
pdfname <- if (!is.null(output_id)) {
  paste0(output_id, "_divergence_density.pdf")
} else {
  "divergence_density.pdf"
}

# Open PDF device
pdf(pdfname, width = 8, height = 6)

# Initialize combined data frame
combined_data <- data.frame()

# Process each file
for (i in seq_along(input_files)) {
  file <- input_files[i]
  filename <- basename(file)

  # Try reading first line
  first_line <- readLines(file, n = 1)

  # Attempt to parse header
  possible_header <- sub("^#", "", first_line)
  header <- unlist(strsplit(possible_header, "\\s+"))

  # Assume Identity is in column 2 unless found in header
  identity_col <- 2
  has_identity <- "Identity" %in% header

  if (has_identity) {
    data <- read.table(file, skip = 1, col.names = header,
                       stringsAsFactors = FALSE, comment.char = "", fill = TRUE)
    identity_name <- "Identity"
  } else {
    data <- read.table(file, header = FALSE, stringsAsFactors = FALSE,
                       comment.char = "", fill = TRUE)
    if (ncol(data) < identity_col) {
      warning("File '", file, "' has fewer than ", identity_col,
              " columns. Skipping this file.")
      next
    }
    identity_name <- paste0("V", identity_col)
    names(data)[identity_col] <- identity_name
  }

  # Extract and convert Identity
  identity <- suppressWarnings(as.numeric(data[[identity_name]]))

  # Warn about non-numeric values
  if (any(is.na(identity))) {
    warning("Non-numeric values found in Identity column in file: ", file)
  }

  # Compute divergence
  divergence <- 1 - identity
  valid_divergence <- divergence[!is.na(divergence) & divergence >= 0 & divergence <= 1]

  # Group assignment
  if (i == 1) {
    group_label <- "Reference"
  } else {
    generation_match <- str_match(filename, "^gen(\\d+)_")
    if (ncol(generation_match) < 2 || is.na(generation_match[1, 2])) {
      warning("Filename does not match expected pattern: gen\\d+_*. File: ", file)
      group_label <- paste0("Unknown_", i)
    } else {
      generation <- as.numeric(generation_match[1, 2])
      group_label <- paste("Generation", generation)
    }
  }

  # Append valid data
  if (length(valid_divergence) > 0) {
    combined_data <- rbind(combined_data,
                           data.frame(
                             divergence = valid_divergence,
                             group = group_label
                           ))
  } else {
    warning("No valid divergence values found in file: ", file)
  }
}

# Generate color palette
groups <- unique(combined_data$group)
ref_color <- c("Reference" = "red")
other_groups <- groups[groups != "Reference"]

# Sort generations numerically
if (length(other_groups) > 0) {
  gen_numbers <- as.numeric(sub("Generation ", "", other_groups))
  order_idx <- order(gen_numbers)
  sorted_groups <- other_groups[order_idx]

  n_colors <- length(sorted_groups)
  palette <- viridis_pal(option = "plasma")(n_colors)
  names(palette) <- sorted_groups
  all_colors <- c(ref_color, palette)
} else {
  all_colors <- ref_color
}

# Create plot
plot <- ggplot(combined_data, aes(x = divergence, group = group, color = group)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = all_colors, name = "Group") +
  labs(x = "Divergence (1 - Identity)", y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

# Save and close
print(plot)
dev.off()
