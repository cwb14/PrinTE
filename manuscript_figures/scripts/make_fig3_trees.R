#!/usr/bin/env Rscript
# Fig. 3D/E helper: text-free circular Tekay RT trees as vector PDFs (cairo), for placement by make_fig3.py.
# Same tree processing and colours as the published LTR_phylo7.R (ggtree circular, groupOTU colouring by
# tip prefix: "shared" #1f78b4, "uniq" #e31a1c). Both trees are drawn on ONE radial scale so that a single
# scale bar is valid for both: branch length `radius` maps to 0.4 x page width (ggplot2 coord_polar),
# i.e. page_mm = mm_per_unit * radius / 0.4. No titles, legends, scale bars or other text.
#
# usage: Rscript make_fig3_trees.R <tree1.nwk> <out1.pdf> <tree2.nwk> <out2.pdf> <mm_per_unit> <radius> <linewidth_pt>
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("usage: make_fig3_trees.R tree1 out1 tree2 out2 mm_per_unit radius linewidth_pt")
mm_per_unit <- as.numeric(args[5])
radius      <- as.numeric(args[6])
lw          <- as.numeric(args[7]) / .pt   # ggplot linewidth units (mm-ish) from points
page_mm     <- mm_per_unit * radius / 0.4

process_tree <- function(treefile) {  # verbatim logic of LTR_phylo7.R::process_tree
  tree <- read.tree(treefile)
  if (is.null(tree)) stop(paste("Tree could not be read:", treefile))
  labels <- tree$tip.label
  group  <- sapply(labels, function(x) {
    tag <- sub("#.*", "", x)
    if      (grepl("^uniq",   tag)) "uniq"
    else if (grepl("^shared", tag)) "shared"
    else                             NA
  })
  keep <- labels[!is.na(group)]
  if (length(keep) == 0) stop(paste("No uniq/shared tips to plot in", treefile))
  tree <- drop.tip(tree, setdiff(labels, keep))
  labels <- tree$tip.label
  group  <- sapply(labels, function(x) {
    tag <- sub("#.*", "", x)
    if (grepl("^uniq", tag)) "uniq" else "shared"
  })
  label_df <- data.frame(label = labels, group = factor(group, levels = c("uniq", "shared")))
  grp_list <- split(label_df$label, label_df$group)
  groupOTU(tree, grp_list, group_name = "group")
}

draw <- function(treefile, out) {
  tr <- process_tree(treefile)
  depth <- max(node.depth.edgelength(tr))
  if (depth > radius) stop(sprintf("%s: max root-to-tip %.3f exceeds radius %.3f", treefile, depth, radius))
  p <- suppressMessages(
    ggtree(tr, aes(color = group), layout = "circular", linewidth = lw) +
      scale_color_manual(values = c("uniq" = "#e31a1c", "shared" = "#1f78b4")) +
      scale_x_continuous(limits = c(0, radius), expand = c(0, 0)) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  )
  ggsave(out, plot = p, width = page_mm, height = page_mm, units = "mm", device = cairo_pdf, bg = "transparent")
  cat(sprintf("%s\ttips=%d\tmax_depth=%.4f\tpage_mm=%.3f\n", out, length(tr$tip.label), depth, page_mm))
}

draw(args[1], args[2])
draw(args[3], args[4])
