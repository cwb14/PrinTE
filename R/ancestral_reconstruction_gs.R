#!/usr/bin/env Rscript

# ancestral_reconstruction_gs.R
## Ancestral genome size (remove --label_nodes for prettier tree. 
# Rscript R/ancestral_reconstruction_gs.R --suffix .fa --newick subset.nwk --abrev abrev.tsv --out ancestral_genome_size --label_nodes
#
## Ancestral LTR-RT size.
# Rscript R/ancestral_reconstruction_gs.R --suffix _ltr.ltrharvest.full_length.dedup.fa.rexdb-plant.cls.lib.fa --newick subset.nwk --abrev abrev.tsv --out ancestral_ltrrt_size --label_nodes
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
              metavar = "X,Y")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$suffix) || is.null(opt$newick)) {
  cat("\nERROR: --suffix and --newick are required.\n\n")
  cat("Example:\n")
  cat("  Rscript ancestral_reconstruction_gs.R --suffix .genome.fa --newick subset.nwk --abrev abrev.tsv\n\n")
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
pdf(pdf_file, width = 7, height = 4)

# longest root-to-tip distance (tree "height") in branch-length units
tree_height <- max(phytools::nodeHeights(cont_obj$tree)[, 2])
tree_height_legend <- as.integer(round(tree_height))  # nearest whole number

plot(
  cont_obj,
  legend  = tree_height_legend,
  fsize   = 0.8,
  outline = FALSE
)

mtext("Genome size (bp)", side = 3, line = 1)

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
