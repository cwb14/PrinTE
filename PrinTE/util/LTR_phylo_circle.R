#!/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
outfig   <- if (length(args) >= 2 && !grepl("^[0-9.]+$", args[2])) args[2] else sub("\\.nwk$", ".pdf", treefile)
custom_scale <- if (length(args) >= 2 && grepl("^[0-9.]+$", args[2])) as.numeric(args[2]) else
                if (length(args) >= 3 && grepl("^[0-9.]+$", args[3])) as.numeric(args[3]) else NA

# Load libraries
library(ggplot2)
library(ggtree)
library(ape)
library(dplyr)

# Read the tree
tree <- read.tree(treefile)
if (is.null(tree)) stop("Tree could not be read.")

# Extract tip labels and assign groups
labels <- tree$tip.label
group  <- sapply(labels, function(x) {
  tag <- sub("#.*", "", x)
  if      (grepl("^uniq",   tag)) "uniq"
  else if (grepl("^shared", tag)) "shared"
  else                             NA
})

# Prune tips not in uniq/shared
keep <- labels[!is.na(group)]
if (length(keep) == 0) stop("No uniq/shared tips to plot.")
tree <- drop.tip(tree, setdiff(labels, keep))

# Recompute labels & groups after pruning
labels <- tree$tip.label
group  <- sapply(labels, function(x) {
  tag <- sub("#.*", "", x)
  if (grepl("^uniq", tag)) "uniq" else "shared"
})
label_df <- data.frame(label = labels,
                       group = factor(group, levels = c("uniq", "shared")))

# Group for branch coloring
grp_list     <- split(label_df$label, label_df$group)
tree_grouped <- groupOTU(tree, grp_list, group_name = "group")

# Compute geometry for placing a horizontal scale bar
depths        <- node.depth.edgelength(tree)
max_depth     <- max(depths[1:length(tree$tip.label)])
n_tips        <- length(tree$tip.label)
radial_offset <- max_depth * 1.05
y_middle      <- (n_tips + 1) / 2

# Build the plot
p <- ggtree(tree_grouped, aes(color = group), layout = "circular") +
  scale_color_manual(values = c("uniq" = "#e31a1c", "shared" = "#1f78b4")) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# Add the scale bar
if (!is.na(custom_scale)) {
  p <- p + geom_treescale(
    x        = radial_offset,
    y        = y_middle,
    width    = custom_scale,
    fontsize = 4,
    linesize = 0.7,
    offset   = 12
  )
} else {
  p <- p + geom_treescale(
    x        = radial_offset,
    y        = y_middle,
    fontsize = 4,
    linesize = 0.7,
    offset   = 12
  )
}

# Save to PDF
ggsave(outfig, plot = p, width = 13.5, height = 8.4, units = "in")
