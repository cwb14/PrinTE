#!/usr/bin/env Rscript

# fig7_contmap.R — Fig. 7A tree for the Nature Communications final figure.
# Copy of /data2/chris/fungi/PrinTE/grid/ancestral_reconstruction_gs.R (2026-08-19)
# with ONLY output-format changes (the reconstruction is untouched):
#   * cairo_pdf(family = --family, pointsize = --pointsize) instead of pdf(),
#     so text is set in Arial and embedded as TrueType. Arial is resolved by
#     fontconfig: run with FONTCONFIG_FILE pointing at a fonts.conf that adds
#     the mscorefonts directory (see make_fig7.py).
#   * --tree_lwd / --bar_lwd / --circle_cex so line weights and node circles
#     can be matched to the final print size (the PDF is drawn at 1:1 size).
#   * also writes <out>.tip_genome_sizes.tsv (tip, genome_bp) for Source Data.
#
# ancestral_reconstruction_flexible.R
## Ancestral genome size (remove --label_nodes for prettier tree. 
# Rscript PrinTE/grid/ancestral_reconstruction_gs.R --suffix .fa --newick subset.nwk --abrev abrev.tsv --out ancestral_genome_size --label_nodes
#
## Ancestral LTR-RT size.
# Rscript PrinTE/grid/ancestral_reconstruction_gs.R --suffix _ltr.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --newick subset.nwk --abrev abrev.tsv --out ancestral_ltrrt_size --label_nodes
#
## Peak at the PDF to figure out with nodes youre interested in. Im interested in 20.
## Ancestral LTR-RT or genome size of node 20.
# cat ancestral_genome_size.ancestral_genome_sizes.tsv | awk '$1 == 20' | cut -f 6
#
# Flexible ancestral genome size reconstruction from FASTA files + Newick tree.
#
# Key behavior:
#   - Taxon IDs (prefix/abbreviations) are inferred directly from the Newick tip labels.
#   - FASTA files are looked for in the SAME directory as the --newick file,
#     with filenames constructed as:  <newick_dir>/<tip_label><suffix>
#
# Options:
#   --suffix        Filename suffix appended to each tip label to locate FASTA files (required)
#   --newick        Newick tree file (required)
#   --abrev         OPTIONAL TSV mapping abbreviations -> full names (2 columns)
#   --out           OPTIONAL output stem for files (default: "ancestral_genome_size")
#   --res           OPTIONAL contMap resolution (default: 200)
#   --label_nodes   OPTIONAL flag: add internal node numbers to the output PDF (default: off)
#   --node_cex      OPTIONAL: size of node-number labels (default: 0.65)
#   --node_adj      OPTIONAL: comma-separated adj for node labels (default: 1.2,-0.2)
#
# Notes:
# - Genome size is computed from a samtools faidx index (.fai): the sum of the
#   sequence-length column (col 2). If a "<fasta>.fai" already exists it is
#   reused; otherwise it is created once with `samtools faidx`. This is far
#   faster than scanning the FASTA in R. Requires `samtools` on PATH.
# - If --abrev is provided:
#     * FASTA reading uses the original tree tip labels (abbreviations) to locate files.
#     * Reconstruction/plotting uses full names (tree tips are relabeled after sizes are read).
# - Tip labels in the tree must be unique.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(optparse)
})

# ----------------------------- CLI -----------------------------------------

option_list <- list(
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Suffix appended to each Newick tip label to form FASTA filename (required). Example: .genome.fa",
              metavar = "SUFFIX"),
  make_option(c("--newick"), type = "character", default = NULL,
              help = "Newick tree file path (required). Example: subset.nwk",
              metavar = "FILE"),
  make_option(c("--abrev"), type = "character", default = NULL,
              help = "Optional TSV mapping: abbreviation<TAB>full_name (2 columns).",
              metavar = "FILE"),
  make_option(c("--out"), type = "character", default = "ancestral_genome_size",
              help = "Output stem for files (default: ancestral_genome_size).",
              metavar = "STEM"),
  make_option(c("--res"), type = "integer", default = 200,
              help = "contMap resolution (default: 200).",
              metavar = "INT"),
  make_option(c("--label_nodes"), action = "store_true", default = FALSE,
              help = "If set, draw internal node numbers on the contMap PDF."),
  make_option(c("--node_cex"), type = "double", default = 0.65,
              help = "Size (cex) of node-number labels when --label_nodes is set (default: 0.65).",
              metavar = "FLOAT"),
  make_option(c("--node_adj"), type = "character", default = "1.2,-0.2",
              help = "Comma-separated adj for node labels (x,y). Default: 1.2,-0.2",
              metavar = "X,Y"),
  make_option(c("--pdf_width"), type = "double", default = 7.0,
              help = "Output PDF width in inches (default: 7).",
              metavar = "FLOAT"),
  make_option(c("--pdf_height"), type = "double", default = 2.25,
              help = "Output PDF height in inches (default: 2.25). Increase to fit a taller embed cell without letterboxing.",
              metavar = "FLOAT"),
  make_option(c("--family"), type = "character", default = "Arial",
              help = "cairo_pdf font family (default: Arial).", metavar = "NAME"),
  make_option(c("--pointsize"), type = "double", default = 12,
              help = "cairo_pdf base pointsize; text = pointsize * cex (default: 12).",
              metavar = "FLOAT"),
  make_option(c("--tree_lwd"), type = "double", default = 4,
              help = "Branch line width (phytools default 4).", metavar = "FLOAT"),
  make_option(c("--bar_lwd"), type = "double", default = 3,
              help = "Colour-bar line width (default 3).", metavar = "FLOAT"),
  make_option(c("--circle_cex"), type = "double", default = 1.3,
              help = "cex of the highlighted-node circles (default 1.3).", metavar = "FLOAT")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$suffix) || is.null(opt$newick)) {
  cat("\nERROR: --suffix and --newick are required.\n\n")
  cat("Example:\n")
  cat("  Rscript ancestral_reconstruction_flexible.R --suffix .genome.fa --newick subset.nwk --abrev abrev.tsv\n\n")
  quit(status = 1)
}

if (!file.exists(opt$newick)) stop(paste("Newick file not found:", opt$newick))

# ------------------------- helpers -----------------------------------------

read_abbrev_tsv <- function(path) {
  x <- read.table(path, sep = "\t", header = FALSE, quote = "", comment.char = "",
                  stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(x) < 2) stop("--abrev file must have at least 2 tab-separated columns: abrev<TAB>full_name")
  ab <- trimws(x[[1]])
  full <- trimws(x[[2]])
  if (any(ab == "" | full == "")) stop("--abrev has empty abbreviation or full-name entries.")
  stats::setNames(full, ab)
}

# Genome size in bp via samtools faidx: sum of the .fai sequence-length column.
# Reuses an existing "<fasta>.fai" if present; otherwise builds it once.
calc_genome_size_bp <- function(fasta_file) {
  if (!file.exists(fasta_file)) stop(paste("FASTA file not found:", fasta_file))

  fai_file <- paste0(fasta_file, ".fai")

  if (!file.exists(fai_file)) {
    out <- suppressWarnings(system2("samtools",
                                    args   = c("faidx", shQuote(fasta_file)),
                                    stdout = TRUE, stderr = TRUE))
    st <- attr(out, "status")
    if (!is.null(st) && st != 0) {
      stop(paste0("samtools faidx failed for ", fasta_file, ":\n",
                  paste(out, collapse = "\n")))
    }
    if (!file.exists(fai_file)) {
      stop(paste("samtools faidx did not produce expected index:", fai_file))
    }
  }

  fai <- read.table(fai_file, sep = "\t", header = FALSE, quote = "",
                    comment.char = "", stringsAsFactors = FALSE)
  if (ncol(fai) < 2) stop(paste("Malformed .fai (need >=2 columns):", fai_file))
  if (nrow(fai) == 0) stop(paste("Empty .fai index:", fai_file))

  # as.numeric (double) avoids 32-bit integer overflow on large genomes.
  sum(as.numeric(fai[[2]]))
}

# Get descendant tip names for a node
get_desc_tip_names <- function(tree, node) {
  all_desc <- phytools::getDescendants(tree, node)
  tip_idx <- all_desc[all_desc <= length(tree$tip.label)]
  tree$tip.label[tip_idx]
}

# Parse "x,y" into numeric length-2 vector
parse_adj <- function(adj_string) {
  parts <- strsplit(adj_string, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  if (length(parts) != 2) stop("--node_adj must be two comma-separated numbers, e.g. 1.2,-0.2")
  out <- suppressWarnings(as.numeric(parts))
  if (any(is.na(out))) stop("--node_adj must be numeric, e.g. 1.2,-0.2")
  out
}

# -------------------------- read tree --------------------------------------

tree <- read.tree(opt$newick)
cat("Read tree:", opt$newick, "\n")
cat("Tips:", length(tree$tip.label), " | Internal nodes:", tree$Nnode, "\n\n")

if (any(duplicated(tree$tip.label))) {
  dups <- unique(tree$tip.label[duplicated(tree$tip.label)])
  cat("ERROR: Duplicate tip labels detected:\n")
  print(dups)
  stop("Tip labels must be unique so genome sizes can be mapped unambiguously.")
}

cat("Tip labels in tree:\n")
print(tree$tip.label)
cat("\n")

# Keep original labels for file lookup (even if we later relabel to full names)
tip_labels_for_files <- tree$tip.label

# --------------------- compute genome sizes --------------------------------

cat("Computing genome sizes from FASTA files...\n")

# FASTA files are expected in the same directory as the Newick tree file.
fasta_dir <- dirname(opt$newick)
cat("Looking for FASTA files in: ", fasta_dir, "\n", sep = "")

# Fail fast if samtools is missing (used to build/read .fai indices).
if (Sys.which("samtools") == "") {
  stop(paste("samtools not found on PATH. samtools is required to index",
             "FASTA files (e.g. `mamba install -c bioconda samtools`)."))
}

genome_bp <- numeric(length(tip_labels_for_files))
names(genome_bp) <- tip_labels_for_files

for (p in tip_labels_for_files) {
  fasta_file <- file.path(fasta_dir, paste0(p, opt$suffix))
  reused <- file.exists(paste0(fasta_file, ".fai"))
  bp <- calc_genome_size_bp(fasta_file)
  genome_bp[p] <- bp
  cat("  ", p, " -> ", fasta_file,
      if (reused) " [reused .fai]" else " [built .fai]",
      " : ", bp, " bp\n", sep = "")
}
cat("\n")

# --------------------- optional abbreviation mapping -----------------------

if (!is.null(opt$abrev)) {
  if (!file.exists(opt$abrev)) stop(paste("Abbrev TSV not found:", opt$abrev))
  abbr_to_full <- read_abbrev_tsv(opt$abrev)

  tree_tips <- tree$tip.label

  if (!all(tree_tips %in% names(abbr_to_full))) {
    missing <- setdiff(tree_tips, names(abbr_to_full))
    cat("ERROR: Some tree tips are not present in the abbreviation column of --abrev:\n")
    print(missing)
    stop("abrev.tsv must include a mapping for every tip label in the Newick.")
  }

  cat("Relabeling tree tips using --abrev (abbrev -> full name).\n")

  tree$tip.label <- unname(abbr_to_full[tree$tip.label])
  names(genome_bp) <- unname(abbr_to_full[names(genome_bp)])

  cat("\nTip labels after relabeling:\n")
  print(tree$tip.label)
  cat("\n")

  if (any(duplicated(tree$tip.label))) {
    dups <- unique(tree$tip.label[duplicated(tree$tip.label)])
    cat("ERROR: Duplicate full names after relabeling:\n")
    print(dups)
    stop("Full names must be unique after mapping; otherwise tips collide.")
  }
}

# --------------------- match genome sizes to (possibly relabeled) tree tips -

if (!all(tree$tip.label %in% names(genome_bp))) {
  missing <- setdiff(tree$tip.label, names(genome_bp))
  cat("ERROR: Some tip labels in the tree do not have genome sizes.\n")
  cat("Missing tips:\n")
  print(missing)
  cat("\nAvailable genome_bp names:\n")
  print(names(genome_bp))
  cat("\n")
  stop("Tip labels in the tree do not match genome-size names.")
}

genome_bp <- genome_bp[tree$tip.label]

cat("Genome sizes (bp) in tree order:\n")
print(genome_bp)
cat("\n")

# --------------------- reconstruct on log10 scale --------------------------

genome_log <- log10(genome_bp)
anc <- fastAnc(tree, genome_log, vars = TRUE, CI = TRUE)

cat("Ancestral states (log10):\n")
print(anc$ace)
cat("\nVariances:\n")
print(anc$var)
cat("\n95% CI (log10):\n")
print(anc$CI95)
cat("\n")

# --------------------- back-transform + representative tips ----------------

node_ids <- as.integer(names(anc$ace))

res <- data.frame(
  node                = node_ids,
  representative_tips = NA_character_,
  log10_est           = anc$ace,
  log10_CI_lower      = anc$CI95[, 1],
  log10_CI_upper      = anc$CI95[, 2],
  bp_est              = 10^anc$ace,
  bp_CI_lower         = 10^anc$CI95[, 1],
  bp_CI_upper         = 10^anc$CI95[, 2],
  row.names           = NULL
)

res <- res[order(res$node), ]

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

res <- res[, c("node", "representative_tips",
               "log10_est", "log10_CI_lower", "log10_CI_upper",
               "bp_est", "bp_CI_lower", "bp_CI_upper")]

cat("Back-transformed ancestral genome sizes (bp) with 95% CI:\n")
print(res)
cat("\n")

# Write TSV
tsv_file <- paste0(opt$out, ".ancestral_genome_sizes.tsv")
write.table(res, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Table written to ", tsv_file, "\n\n", sep = "")

# --------------------- contMap on bp scale ---------------------------------

cat("Generating contMap for genome size (bp)...\n")
anc_bp <- 10^anc$ace

cont_obj <- contMap(
  tree,
  x          = genome_bp,
  anc.states = anc_bp,
  plot       = FALSE,
  res        = opt$res
)

pdf_file <- paste0(opt$out, ".contMap.pdf")
# 7:2.25 aspect (shorter than the old 7:3) so the embedded banner in
# plots/panel_plot.py takes less vertical space. The Python --phylo-height-
# ratio default (1.2896) is paired to this ratio so the tree still fills the
# full figure width with no letterboxing; change both together.
tip_file <- paste0(opt$out, ".tip_genome_sizes.tsv")
write.table(data.frame(tip = names(genome_bp), genome_bp = as.numeric(genome_bp)),
            file = tip_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Tip sizes written to ", tip_file, "\n", sep = "")

cairo_pdf(pdf_file, width = opt$pdf_width, height = opt$pdf_height,
          family = opt$family, pointsize = opt$pointsize)

# longest root-to-tip distance (tree "height") in branch-length units;
# the tree is time-calibrated in millions of generations, so this is also the
# bar length in millions of generations.
tree_height <- max(phytools::nodeHeights(cont_obj$tree)[, 2])
tree_height_legend <- as.integer(round(tree_height))  # nearest whole number

# Plot the tree only. We suppress phytools' built-in colour bar (legend =
# FALSE) and draw our own publication-quality dual legend below it, because
# the built-in path force-formats the end labels as raw bp (round(lims,sig))
# and offers no subtitle for a rightwards bar. `underscore = TRUE` keeps the
# real tip labels (e.g. Melme_tre1); the phytools default rewrites "_"->" ".
# ylim mirrors exactly what plot.densityMap reserves for its own legend
# (c(1 - 0.12*(N-1), N)), so the figure geometry is unchanged.
n_tip <- length(cont_obj$tree$tip.label)
plot(
  cont_obj,
  legend     = FALSE,
  fsize      = 0.85,
  lwd        = opt$tree_lwd,
  outline    = FALSE,
  underscore = TRUE,
  ylim       = c(1 - 0.12 * (n_tip - 1), n_tip)
)

# ---- custom dual legend -----------------------------------------------------
# Colour encodes genome size; the bar's drawn length spans tree_height_legend
# branch-length units, i.e. that many million generations of evolutionary time on the tree's
# x-axis. End labels report the genome-size range in Mb, rounded and
# thousands-separated (e.g. 113 ... 1,018) instead of raw bp (1018398822).
mb_fmt <- function(v) formatC(as.integer(round(v / 1e6)),
                              big.mark = ",", format = "d")
lo_lab <- mb_fmt(cont_obj$lims[1])
hi_lab <- mb_fmt(cont_obj$lims[2])

# --- legend layout tunables --------------------------------------------------
# BAR_DROP : push the whole scale bar (and its labels) this many tree-y units
#            below the phytools default position. Larger = lower.
# LAB_OFF  : vertical gap of the title + Mb end values ABOVE the bar.
# SUB_OFF  : vertical gap of the "Bar length = N million generations"
#            subtitle BELOW the bar.
#            LAB_OFF/SUB_OFF are R text() `offset` values (fractions of a
#            character height for pos = 3/1); the phytools/text() default is
#            0.5, so values < 0.5 sit the words closer to the bar.
BAR_DROP <- 0.30
LAB_OFF  <- 0.32
SUB_OFF  <- 0.32

bar_x <- 0
bar_y <- 1 - 0.08 * (n_tip - 1) - BAR_DROP   # phytools default, nudged down

# Draw the colour bar ONLY. Passing "" (not NULL) for title/subtitle makes
# add.color.bar's internal text() calls render empty strings (nothing), and
# lims = NULL suppresses its raw-bp end labels — so we control every word's
# distance from the bar ourselves, below, via LAB_OFF / SUB_OFF.
add.color.bar(
  leg      = tree_height_legend,
  cols     = cont_obj$cols,
  title    = "",
  subtitle = "",
  lims     = NULL,
  prompt   = FALSE,
  x        = bar_x,
  y        = bar_y,
  lwd      = opt$bar_lwd,
  fsize    = 0.85,
  outline  = FALSE
)

# Title + Mb end values, snug ABOVE the bar (pos = 3, tight offset).
text(x = bar_x + tree_height_legend / 2, y = bar_y, "Genome Size (Mb)",
     pos = 3, offset = LAB_OFF, cex = 0.85)
text(x = bar_x,                      y = bar_y, lo_lab,
     pos = 3, offset = LAB_OFF, cex = 0.85)
text(x = bar_x + tree_height_legend, y = bar_y, hi_lab,
     pos = 3, offset = LAB_OFF, cex = 0.85)

# Subtitle, snug BELOW the bar (pos = 1, tight offset).
text(x = bar_x + tree_height_legend / 2, y = bar_y,
     paste0("Bar length = ", tree_height_legend, " million generations"),
     pos = 1, offset = SUB_OFF, cex = 0.85)

# ---- Highlight selected divergence nodes -----------------------------------
# Mark the MRCA of each tip pair so the reader can cross-reference these
# splits with downstream panels. The fill colour of each circle matches the
# species colour used in the corresponding row of plots/panel_plot.py:
#   Crori1/Corcom2          -> green  (#228833, Cronartium ribicola row)
#   PuccoSD80_1/PuccoNC29_1 -> blue   (#4477AA, Puccinia coronata row)
#   myrtle_rust/PuccoNC29_1 -> purple (#AA3377, Austropuccinia psidii row)
highlight_pairs <- list(
  list(tips = c("Crori1",      "Corcom2"),     col = "#228833"),
  list(tips = c("PuccoSD80_1", "PuccoNC29_1"), col = "#4477AA"),
  list(tips = c("myrtle_rust", "PuccoNC29_1"), col = "#AA3377")
)
highlight_nodes <- integer(0)
highlight_cols  <- character(0)
for (entry in highlight_pairs) {
  pair <- entry$tips
  if (all(pair %in% cont_obj$tree$tip.label)) {
    n_mrca <- ape::getMRCA(cont_obj$tree, pair)
    if (!is.null(n_mrca)) {
      highlight_nodes <- c(highlight_nodes, n_mrca)
      highlight_cols  <- c(highlight_cols,  entry$col)
    }
  } else {
    missing <- setdiff(pair, cont_obj$tree$tip.label)
    cat("WARN: highlight pair tip(s) not in tree, skipping: ",
        paste(missing, collapse = ", "), "\n", sep = "")
  }
}
if (length(highlight_nodes) > 0) {
  nodelabels(
    node = highlight_nodes,
    pch  = 21,
    cex  = opt$circle_cex,
    bg   = highlight_cols,
    col  = "black",
    lwd  = 1.0
  )
}

# ---- OPTIONAL: node numbers on the plot ----
if (isTRUE(opt$label_nodes)) {
  adj_xy <- parse_adj(opt$node_adj)

  n_tips <- length(cont_obj$tree$tip.label)
  n_nodes <- cont_obj$tree$Nnode
  internal_nodes <- (n_tips + 1):(n_tips + n_nodes)

  # Draw internal node numbers
  nodelabels(
    text  = internal_nodes,
    node  = internal_nodes,
    frame = "none",
    cex   = opt$node_cex,
    adj   = adj_xy
  )

  # Also add a small note so you remember the scheme
  mtext("Internal node numbers shown", side = 1, line = 2, cex = 0.7)
}

dev.off()

cat("Contour map phylogeny written to ", pdf_file, "\n", sep = "")
cat("Done.\n")

