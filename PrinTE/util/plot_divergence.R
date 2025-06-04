#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_divergence7.R <input1> [<input2> ...] <output_id> [use_identity]", call. = FALSE)
}

# Check if last arg is TRUE or FALSE
last_arg <- args[length(args)]
use_identity <- FALSE
if (last_arg == "TRUE" || last_arg == "FALSE") {
  use_identity <- as.logical(last_arg)
  output_id    <- args[length(args) - 1]
  input_files  <- args[1:(length(args) - 2)]
} else {
  output_id   <- last_arg
  input_files <- args[1:(length(args) - 1)]
}
pdfname <- paste0(output_id, "_divergence_density.pdf")

combined_data <- data.frame()

for (file in input_files) {
  fn <- basename(file)
  cat("Processing:", fn, "\n")
  
  # Detect header style
  first_line <- readLines(file, n = 1)
  if (startsWith(first_line, "#")) {
    header_fields <- strsplit(sub("^#", "", first_line), "\\s+")[[1]]
    df <- tryCatch({
      read.table(file, skip = 1, col.names = header_fields,
                 stringsAsFactors = FALSE, comment.char = "", fill = TRUE)
    }, error = function(e) {
      warning("Could not read file: ", fn, "\n", e)
      return(NULL)
    })
  } else {
    df <- tryCatch({
      read.table(file, header = TRUE,
                 stringsAsFactors = FALSE, comment.char = "", fill = TRUE)
    }, error = function(e) {
      warning("Could not read file: ", fn, "\n", e)
      return(NULL)
    })
  }
  if (is.null(df)) next
  
  # Decide which field to use
  divergence <- NULL
  divergence_source <- NULL

  if (use_identity && "Identity" %in% names(df)) {
    divergence <- 1 - as.numeric(df$Identity)
    divergence_source <- "1 - Identity"
  } else if (!use_identity && "Insertion_Time" %in% names(df)) {
    divergence <- as.numeric(df$Insertion_Time)
    divergence_source <- "Insertion_Time"
  } else if ("K2P_d" %in% names(df)) {
    divergence <- as.numeric(df$K2P_d)
    divergence_source <- "K2P_d"
  } else if ("Identity" %in% names(df)) {
    divergence <- 1 - as.numeric(df$Identity)
    divergence_source <- "1 - Identity"
  } else {
    warning("No valid divergence field found in: ", fn)
    next
  }

  # Filter valid
  valid <- !is.na(divergence)
  if (divergence_source != "Insertion_Time") {
    valid <- valid & divergence >= 0 & divergence <= 1
  }
  divergence <- divergence[valid]
  if (length(divergence) == 0) {
    warning("No valid divergence values in: ", fn)
    next
  }

  # Assign group
  if (grepl("^gen(\\d+)_", fn, ignore.case = TRUE)) {
    gen_num    <- as.numeric(sub("^gen(\\d+)_.*", "\\1", fn, ignore.case = TRUE))
    group_name <- paste("Generation", gen_num)
  } else {
    group_name <- "Reference"
  }

  combined_data <- rbind(
    combined_data,
    data.frame(divergence = divergence, group = group_name, source = divergence_source)
  )
}

if (nrow(combined_data) == 0) {
  stop("No valid data to plot.")
}

# Color palette
groups       <- unique(combined_data$group)
ref_color    <- c(Reference = "red")
other_groups <- setdiff(groups, "Reference")

if (length(other_groups) > 0) {
  gen_nums <- as.numeric(sub("Generation ", "", other_groups))
  ord      <- order(gen_nums)
  sorted_g <- other_groups[ord]
  pal      <- viridis_pal(option = "plasma")(length(sorted_g))
  names(pal) <- sorted_g
  all_colors <- c(ref_color, pal)
} else {
  all_colors <- ref_color
}

x_label <- unique(combined_data$source)
if (length(x_label) > 1) x_label <- "Divergence / Time"

# Plotting
pdf(pdfname, width = 8, height = 6)
p <- ggplot(combined_data, aes(x = divergence, color = group)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = all_colors, name = "Group") +
  labs(x = x_label, y = "Density") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position    = "right")
print(p)
dev.off()
