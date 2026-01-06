#!/usr/bin/env Rscript

# ancestral_reconstruction.R
# Estimate ancestral genome sizes and print back-transformed values (bp)
# with 95% confidence intervals and representative descendant tips.
# Also outputs a contour map phylogeny (contMap) of genome size (bp).

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

## 1. Read tree -------------------------------------------------------------

#tree_file <- "iqtree_rerooted.nwk"
tree_file <- "iqtree_rerooted.dated.nwk"
tree <- read.tree(tree_file)

cat("Read tree with", length(tree$tip.label), "tips and",
    tree$Nnode, "internal nodes.\n\n")

## 2. Tip genome sizes (bp) -------------------------------------------------
## Use full species names

# Mapping from abbreviations (if present) to full names
abbr_to_full <- c(
  Aaren = "Arabidopsis arenosa",
  Ahall = "Arabidopsis halleri",
  Alyra = "Arabidopsis lyrata",
  Chisp = "Camelina hispida",
  Athal = "Arabidopsis thaliana"
)

# If the tree uses abbreviations, convert to full names
if (all(tree$tip.label %in% names(abbr_to_full))) {
  tree$tip.label <- abbr_to_full[tree$tip.label]
  cat("Converted abbreviated tip labels to full species names.\n")
}

cat("Tip labels in tree:\n")
print(tree$tip.label)
cat("\n")

# Genome sizes in bp, named by full species names
genome_bp <- c(
  "Arabidopsis thaliana"   = 119668634,
  "Camelina hispida"       = 283323681,
  "Arabidopsis arenosa"    = 149659193,
  "Arabidopsis halleri"    = 227207528,
  "Arabidopsis lyrata"     = 187501846
)

# Cummulative intact LTR-RT lengths.
#genome_bp <- c(
#  "Arabidopsis thaliana" = 5692421,
#  "Camelina hispida" = 63180233,
#  "Arabidopsis arenosa" = 9610561,
#  "Arabidopsis halleri" = 36203795,
#  "Arabidopsis lyrata" = 20200266
#)

# Reorder genome_bp to match tree tip order
if (!all(tree$tip.label %in% names(genome_bp))) {
  stop("Some tip labels in the tree do not match the expected species names.")
}
genome_bp <- genome_bp[tree$tip.label]

cat("Genome sizes (bp) in tree order:\n")
print(genome_bp)
cat("\n")

## 3. Log-transform and reconstruct ----------------------------------------

# log10-transform genome sizes for ancestral reconstruction
genome_log <- log10(genome_bp)

anc <- fastAnc(tree, genome_log, vars = TRUE, CI = TRUE)

cat("Ancestral states on log10 scale:\n")
print(anc$ace)
cat("\nVariances of ancestral estimates:\n")
print(anc$var)
cat("\n95% CI (log10 scale):\n")
print(anc$CI95)
cat("\n")

## 4. Back-transform to bp --------------------------------------------------

node_ids <- as.integer(names(anc$ace))

res <- data.frame(
  node               = node_ids,
  representative_tips = NA_character_,  # filled later
  log10_est          = anc$ace,
  log10_CI_lower     = anc$CI95[, 1],
  log10_CI_upper     = anc$CI95[, 2],
  bp_est             = 10^anc$ace,
  bp_CI_lower        = 10^anc$CI95[, 1],
  bp_CI_upper        = 10^anc$CI95[, 2],
  row.names          = NULL
)

res <- res[order(res$node), ]

## 5. Add representative descendant tips -----------------------------------

# Helper: get all descendant tips of a node
get_desc_tip_names <- function(tree, node) {
  all_desc <- phytools::getDescendants(tree, node)
  # getDescendants returns indices of *all* descendants (nodes + tips)
  tip_idx <- all_desc[all_desc <= length(tree$tip.label)]
  tree$tip.label[tip_idx]
}

rep_tips <- character(nrow(res))

for (i in seq_len(nrow(res))) {
  n <- res$node[i]
  tips <- get_desc_tip_names(tree, n)
  
  if (length(tips) >= 2) {
    rep_tips[i] <- paste0("(", tips[1], ", ", tips[2], ")")
  } else if (length(tips) == 1) {
    rep_tips[i] <- paste0("(", tips[1], ")")
  } else {
    rep_tips[i] <- "(?)"
  }
}

res$representative_tips <- rep_tips

# Order columns
res <- res[, c("node", "representative_tips",
               "log10_est", "log10_CI_lower", "log10_CI_upper",
               "bp_est",    "bp_CI_lower",    "bp_CI_upper")]

cat("Back-transformed ancestral genome sizes (bp) with 95% CI:\n")
print(res)
cat("\n")

## 6. Optional: Print full descendant sets ----------------------------------
# Uncomment if desired:
#
# cat("\nNode -> ALL descendant tips:\n")
# for (i in seq_len(nrow(res))) {
#   n <- res$node[i]
#   tips <- get_desc_tip_names(tree, n)
#   cat("Node", n, ":", paste(tips, collapse = ", "), "\n")
# }
# cat("\n")

## 7. Contour map phylogeny (contMap) ---------------------------------------
## Color scale = genome size in bp (not log10)

cat("Generating contour map phylogeny (contMap) for genome size (bp)...\n")

# Back-transform ancestral states to bp for coloring
anc_bp <- 10^anc$ace

# Build contMap object without plotting first, using bp as the trait
cont_obj <- contMap(
  tree,
  x          = genome_bp,  # tip values in bp
  anc.states = anc_bp,     # ancestral values in bp
  plot       = FALSE,
  res        = 200
)

# (Optional) customize colors:
# cont_obj$cols <- setMap(c("blue", "yellow", "red"))

pdf_file <- "ancestral_genome_size_contMap.pdf"
pdf(pdf_file, width = 7, height = 4)

plot(
  cont_obj,
  legend  = TRUE,   # include color legend (in bp)
  fsize   = 0.8,    # tip label font size
  outline = FALSE   # cleaner look
)

mtext("Genome size (bp)", side = 3, line = 1)

dev.off()

cat("Contour map phylogeny written to", pdf_file, "\n")
