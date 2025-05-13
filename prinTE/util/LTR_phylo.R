#!/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
outfig  <- if (length(args) >= 2) args[2] else sub("\\.nwk$", ".pdf", treefile)

# Load required libraries
library(ggplot2)
library(ggtree)
library(treeio)
library(ape)
library(dplyr)

# Read the tree
tree <- read.tree(treefile)
if (is.null(tree)) stop("Tree could not be read.")

labels <- tree$tip.label

# Determine group from label prefix (before the '#')
group <- sapply(labels, function(x) {
  tag <- sub("#.*", "", x)
  if      (grepl("^uniq",   tag)) "uniq"
  else if (grepl("^shared", tag)) "shared"
  else                             NA
})

# Filter out non-matching tips
keep_tips <- labels[!is.na(group)]
if (length(keep_tips) == 0) stop("No matching tips found for 'uniq' or 'shared'. Tree is empty after filtering.")
tree <- drop.tip(tree, setdiff(labels, keep_tips))
if (is.null(tree)) stop("All tips were dropped; no tree left to plot.")

# Update labels and groups after pruning
labels <- tree$tip.label
group  <- sapply(labels, function(x) {
  tag <- sub("#.*", "", x)
  if (grepl("^uniq", tag)) "uniq" else "shared"
})
label_df <- data.frame(label = labels,
                       group = factor(group, levels = c("uniq","shared")))

# Group tree tips by category for branch coloring
grp_list     <- split(label_df$label, label_df$group)
tree_grouped <- groupOTU(tree, grp_list, group_name = "group")

# Plot the tree with scale bar using rectangular layout
p <- ggtree(tree_grouped, aes(color = group), layout = "rectangular") +
  scale_color_manual(
    values = c("uniq" = "#e31a1c", "shared" = "#1f78b4")
  ) +
  geom_treescale(fontsize = 4, linesize = 0.7) +  # autoscaled distance bar
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save to PDF
ggsave(outfig, plot = p, width = 13.5, height = 8.4, units = "in")
