#!/usr/bin/env Rscript

# Plot an LTR-RT tree with tips coloured by whether their label is prefixed
# 'uniq' or 'shared'. Tips matching neither prefix are pruned.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(dplyr)
  library(optparse)
})

option_list <- list(
  make_option(c("-l", "--layout"), type = "character", default = "rectangular",
              help = "Tree layout: rectangular or circular [default %default]"),
  make_option(c("-s", "--scale"), type = "double", default = NA,
              help = "Scale-bar width in tree distance units. Circular layout only; omit to autoscale."),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output PDF [default: the tree file with .nwk swapped for .pdf]")
)

opt_parser <- OptionParser(usage = "%prog [options] tree.nwk [out.pdf]",
                           option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = TRUE)

if (length(opt$args) == 0) stop("You must provide a Newick tree as a positional argument.")

layout <- match.arg(opt$options$layout, c("rectangular", "circular"))
treefile <- opt$args[1]
outfig <- if (!is.null(opt$options$out)) opt$options$out
          else if (length(opt$args) >= 2) opt$args[2]
          else sub("\\.nwk$", ".pdf", treefile)

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

p <- ggtree(tree_grouped, aes(color = group), layout = layout) +
  scale_color_manual(
    values = c("uniq" = "#e31a1c", "shared" = "#1f78b4")
  )

if (layout == "circular") {
  # A circular layout gives geom_treescale nowhere sensible to autoplace itself,
  # so put the bar just outside the tips at the vertical midpoint.
  depths        <- node.depth.edgelength(tree)
  max_depth     <- max(depths[1:length(tree$tip.label)])
  radial_offset <- max_depth * 1.05
  y_middle      <- (length(tree$tip.label) + 1) / 2

  scale_args <- list(x = radial_offset, y = y_middle,
                     fontsize = 4, linesize = 0.7, offset = 12)
  if (!is.na(opt$options$scale)) scale_args$width <- opt$options$scale

  p <- p + do.call(geom_treescale, scale_args) +
    theme(
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
} else {
  p <- p + geom_treescale(fontsize = 4, linesize = 0.7) +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.margin = margin(10, 10, 10, 10)
    )
}

ggsave(outfig, plot = p, width = 13.5, height = 8.4, units = "in")
