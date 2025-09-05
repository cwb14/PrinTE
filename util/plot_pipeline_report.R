#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(optparse)
})

# ----------------------------
# Command-line option parsing
# ----------------------------
option_list <- list(
  make_option(
    c("-c", "--collapse"),
    action = "store_true",
    default = FALSE,
    help = "Collapse TE_inserts into a single total value"
  )
)

opt_parser <- OptionParser(usage = "%prog [options] input_file", option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = TRUE)

if (length(opt$args) == 0) {
  stop("You must provide the input file as a positional argument.")
}

input_file <- opt$args[1]
collapse <- opt$options$collapse

# ----------------------------
# Load and preprocess data
# ----------------------------
raw_data <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE)

# Process generation and TE insertions
te_parsed <- raw_data |>
  mutate(
    generation_m = Generation / 1e6,
    total_insertions = as.integer(str_extract(TE_inserts.nest.nonnest., "^\\d+")),
    nest_insertions = if (!collapse) as.integer(str_match(TE_inserts.nest.nonnest., "\\((\\d+)/")[, 2]) else NA_integer_,
    nonnest_insertions = if (!collapse) as.integer(str_match(TE_inserts.nest.nonnest., "/(\\d+)\\)")[, 2]) else NA_integer_,
    calculated_deletions_k = Calculated_TE_deletions / 1e3,
    actual_deletions_k = Actual_TE_deletions / 1e3,
    total_insertions_k = total_insertions / 1e3,
    nest_insertions_k = nest_insertions / 1e3,
    nonnest_insertions_k = nonnest_insertions / 1e3
  )

# Choose columns to pivot
if (collapse) {
  plot_data <- te_parsed |>
    select(generation_m,
           calculated_deletions_k,
           actual_deletions_k,
           total_insertions_k) |>
    pivot_longer(
      cols = -generation_m,
      names_to = "Type",
      values_to = "Count_k"
    ) |>
    mutate(
      Type = factor(
        Type,
        levels = c("calculated_deletions_k", "actual_deletions_k", "total_insertions_k"),
        labels = c("Calculated deletions", "Actual deletions", "TE insertions")
      )
    )
} else {
  plot_data <- te_parsed |>
    select(generation_m,
           calculated_deletions_k,
           actual_deletions_k,
           nest_insertions_k,
           nonnest_insertions_k) |>
    pivot_longer(
      cols = -generation_m,
      names_to = "Type",
      values_to = "Count_k"
    ) |>
    mutate(
      Type = factor(
        Type,
        levels = c("calculated_deletions_k", "actual_deletions_k",
                   "nest_insertions_k", "nonnest_insertions_k"),
        labels = c("Calculated deletions", "Actual deletions",
                   "Nest TE insertions", "Non-nest TE insertions")
      )
    )
}

# ----------------------------
# Plot aesthetics
# ----------------------------
line_colors <- c(
  "Calculated deletions" = "#1b9e77",
  "Actual deletions" = "#d95f02",
  "Nest TE insertions" = "#7570b3",
  "Non-nest TE insertions" = "#e7298a",
  "TE insertions" = "#7570b3"
)

line_styles <- c(
  "Calculated deletions" = "solid",
  "Actual deletions" = "dashed",
  "Nest TE insertions" = "solid",
  "Non-nest TE insertions" = "solid",
  "TE insertions" = "solid"
)

# ----------------------------
# Generate plot
# ----------------------------
plot <- ggplot(plot_data, aes(x = generation_m, y = Count_k,
                              color = Type, linetype = Type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_styles) +
  scale_x_continuous(
    name = "Generation (millions)",
    expand = c(0, 0),
    limits = range(te_parsed$generation_m)
  ) +
  scale_y_continuous(
    name = "Count (thousands)",
    expand = c(0, 0),
    limits = c(0, max(plot_data$Count_k, na.rm = TRUE))
  ) +
  labs(
    title = "TE Insertions and Deletions Over Generations"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# ----------------------------
# Save output
# ----------------------------
ggsave("te_dynamics_plot.pdf", plot = plot, width = 10, height = 7, dpi = 300)
message("Saved plot to 'pipeline_report.pdf'")
