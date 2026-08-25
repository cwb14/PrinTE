#!/usr/bin/env Rscript

# ancestral_reconstruction_ltr_age.R
#
# For each tip taxon, run summary_stats.py on column 7 from:
#   <taxon><suffix>
# Then reconstruct ancestral totals + per-bin counts at internal nodes via fastAnc.
#
# Key fix vs your current version:
# - “representative tips” is NOT unique/reliable (nested clades can share the same first 1–2 tips).
# - This script instead writes:
#     (1) a node map TSV that includes the FULL descendant-tip set for each node (unique identifier),
#     (2) a PDF of the tree with internal node numbers plotted on the tree,
#     (3) includes clade_size + clade_signature + desc_tips in every output row.
#
# Usage:
#   Rscript ancestral_reconstruction_ltr_age.R --newick subset.nwk --suffix _ltr_kmer2ltr_dedup
#
# Peak at the PDF to figure out which nodes you need. I need 20.
#   cat ancestral_summary_bins.tsv | awk '$2 == 20'  | cut -f 7-9 | awk -F'\t' 'BEGIN{OFS="\t"} { $3 = sprintf("%.0f", $3); print }' > crori.ancestral.LTR.bins
#
# Convert to frequency.
#   awk -v OFS="\t" '{a[NR]=$0; v[NR]=$3; s+=$3} END{for(i=1;i<=NR;i++) print a[i], v[i]/s}' crori.ancestral.LTR.bins | cut -f 1,2,4 > crori.ancestral.LTR.bins.freq
#
# Optional:
#   --summary_path PATH   path to summary_stats.py (default: run it as a module)
#   --abrev FILE          two-column TSV: <abbrev> <full_name>
#   --bins INT
#   --bin_max NUM
#   --python CMD
#   --tree_pdf FILE       default: tree_with_node_ids.pdf
#   --node_map FILE       default: node_map.tsv
#   --max_tips_sig INT    number of tips to show in signature (default 6)

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

## ------------------------ 0. Defaults --------------------------------------

default_newick      <- "iqtree_rerooted.dated.nwk"
default_suffix      <- ".LTRs.alns.results"
default_bins        <- 50
default_bin_max     <- 0.15
default_python_cmd  <- "python"

default_tree_pdf    <- "tree_with_node_ids.pdf"
default_node_map    <- "node_map.tsv"
default_max_tips_sig <- 6

totals_outfile <- "ancestral_summary_totals.tsv"
bins_outfile   <- "ancestral_summary_bins.tsv"

## ------------------------ 0b. CLI parsing ----------------------------------

print_usage_and_exit <- function(exit_code = 0) {
  cat(
"Usage:
  ancestral_reconstruction_ltr_age.R [--newick FILE] [--suffix SUFFIX]
                                     [--summary_path PATH] [--abrev FILE]
                                     [--bins INT] [--bin_max NUM]
                                     [--python CMD]
                                     [--tree_pdf FILE] [--node_map FILE]
                                     [--max_tips_sig INT]

Required inputs:
  --newick FILE
      Newick tree (tips must match result-file prefixes unless --abrev is used).
      Default: iqtree_rerooted.dated.nwk

  --suffix SUFFIX
      Results-file suffix appended to each tip label (or abbreviation).
      Default: .LTRs.alns.results

Optional:
  --summary_path PATH
      Path to a summary_stats.py file. If not provided, it is run from the
      installed printe package as: python -m printe.grid.summary_stats

  --abrev FILE
      Two-column TSV mapping: <abbrev> <full_name>
      If provided, tree tip labels that match <abbrev> are converted to <full_name>
      for reconstruction/output, but file lookup still uses <abbrev>.

  --bins INT
      Number of bins for summary_stats.py. Default: 50

  --bin_max NUM
      Maximum value for binning (passed to summary_stats.py). Default: 0.15

  --python CMD
      Python executable (python or python3). Default: python

  --tree_pdf FILE
      Output PDF of the tree with internal node IDs labeled.
      Default: tree_with_node_ids.pdf

  --node_map FILE
      Output TSV mapping node ID -> descendant tips (unique clade identifier).
      Default: node_map.tsv

  --max_tips_sig INT
      How many tips to include in the short clade signature.
      Default: 6
", sep = ""
  )
  quit(status = exit_code)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && any(args %in% c("-h", "--help"))) print_usage_and_exit(0)

get_arg_value <- function(flag, args, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[i + 1]
}

newick_file    <- get_arg_value("--newick",       args, default_newick)
suffix         <- get_arg_value("--suffix",       args, default_suffix)
summary_path   <- get_arg_value("--summary_path", args, NA_character_)
abrev_file     <- get_arg_value("--abrev",        args, NA_character_)
bins           <- as.integer(get_arg_value("--bins",    args, as.character(default_bins)))
bin_max        <- as.numeric(get_arg_value("--bin_max", args, as.character(default_bin_max)))
python_cmd     <- get_arg_value("--python",       args, default_python_cmd)

tree_pdf_file  <- get_arg_value("--tree_pdf",     args, default_tree_pdf)
node_map_file  <- get_arg_value("--node_map",     args, default_node_map)
max_tips_sig   <- as.integer(get_arg_value("--max_tips_sig", args, as.character(default_max_tips_sig)))

if (is.na(bins) || bins <= 0) stop("--bins must be a positive integer.")
if (is.na(bin_max) || bin_max <= 0) stop("--bin_max must be a positive number.")
if (is.na(max_tips_sig) || max_tips_sig <= 0) stop("--max_tips_sig must be a positive integer.")
if (!file.exists(newick_file)) stop("Newick file not found: ", newick_file)

## ------------------------ 0c. Resolve script directory ----------------------

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE))
  }
  return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

script_dir <- get_script_dir()

## ------------------------ 0d. Resolve summary_stats -------------------------

# summary_stats ships as part of the printe package, so the default is to run it as
# a module. --summary_path still takes a plain file path for anyone working from a
# loose checkout. Held as a pre-quoted fragment because the module form must not be
# shQuote'd - "-m printe.grid.summary_stats" is two arguments, not one.
if (is.na(summary_path) || summary_path == "") {
  summary_invocation <- "-m printe.grid.summary_stats"
  summary_label      <- "python -m printe.grid.summary_stats"
} else {
  if (!file.exists(summary_path)) {
    stop(
      "summary_stats.py not found at: ", summary_path, "\n",
      "Omit --summary_path to run it from the installed printe package instead."
    )
  }
  summary_invocation <- shQuote(summary_path)
  summary_label      <- summary_path
}

## ------------------------ 1. Read tree -------------------------------------

tree <- read.tree(newick_file)

cat("Tree file:", newick_file, "\n")
cat("Read tree with", length(tree$tip.label), "tips and", tree$Nnode, "internal nodes.\n\n")
cat("Original tip labels in tree:\n")
print(tree$tip.label)
cat("\n")

## ------------------------ 1b. Optional abbreviation -> full names ----------

# By default: identity mapping (tip label is the file prefix).
# Names are analysis labels; values are file-prefix abbreviations used to find files.
tip_to_abbr <- setNames(tree$tip.label, tree$tip.label)

if (!is.na(abrev_file) && abrev_file != "") {
  if (!file.exists(abrev_file)) stop("Abbreviation mapping file not found: ", abrev_file)

  ab <- tryCatch(
    read.table(abrev_file, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE),
    error = function(e) stop("Error reading --abrev file: ", conditionMessage(e))
  )
  if (ncol(ab) < 2) stop("--abrev file must have at least 2 columns: <abbrev> <full_name>")

  abbr_to_full <- setNames(ab[[2]], ab[[1]])

  if (all(tree$tip.label %in% names(abbr_to_full))) {
    old_tips <- tree$tip.label
    tree$tip.label <- unname(abbr_to_full[tree$tip.label])

    full_to_abbr <- setNames(names(abbr_to_full), abbr_to_full)
    tip_to_abbr <- full_to_abbr[tree$tip.label]

    cat("Converted abbreviated tip labels to full species names using:", abrev_file, "\n\n")
    cat("Updated tip labels in tree (used for analysis/output):\n")
    print(tree$tip.label)
    cat("\n")
  } else {
    stop(
      "With --abrev provided, ALL tree tip labels must be abbreviations found in the first column of --abrev.\n",
      "Tree tips did not all match abbreviations in: ", abrev_file
    )
  }
}

tip_species <- tree$tip.label  # analysis labels

## ------------------------ 1c. Force internal-node numbering visibility ------

# IMPORTANT: fastAnc returns ace names as internal node numbers:
#   (Ntip+1) ... (Ntip+Nnode)
Ntip  <- length(tree$tip.label)
Nnode <- tree$Nnode
node_ids <- (Ntip + 1):(Ntip + Nnode)

## ------------------------ 2. Helper: run summary_stats.py -------------------

run_summary_stats_for_values <- function(values,
                                         bins,
                                         bin_max,
                                         python_cmd,
                                         summary_invocation) {
  tmpfile <- tempfile(pattern = "summary_input_", fileext = ".txt")
  on.exit(unlink(tmpfile), add = TRUE)

  writeLines(format(values, digits = 10, scientific = FALSE), con = tmpfile)

  cmd <- sprintf(
    '%s %s --bins %d --bin_max %g --input %s',
    shQuote(python_cmd),
    summary_invocation,
    bins,
    bin_max,
    shQuote(tmpfile)
  )

  out <- tryCatch(
    system(cmd, intern = TRUE),
    error = function(e) stop("Error running summary_stats.py: ", conditionMessage(e))
  )

  if (length(out) == 0) stop("summary_stats.py produced no output.")

  total_read <- NA_real_
  total_used <- NA_real_

  bins_df <- data.frame(
    bin_index = integer(),
    bin_start = numeric(),
    bin_end   = numeric(),
    count     = numeric(),
    stringsAsFactors = FALSE
  )

  for (line in out) {
    line <- trimws(line)
    if (line == "") next

    if (startsWith(line, "#")) {
      if (grepl("^#\\s*total_values_read", line)) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) total_read <- as.numeric(parts[2])
      } else if (grepl("^#\\s*total_values_used_", line)) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) total_used <- as.numeric(parts[2])
      }
    } else {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) < 4) next
      bins_df <- rbind(
        bins_df,
        data.frame(
          bin_index = as.integer(parts[1]),
          bin_start = as.numeric(parts[2]),
          bin_end   = as.numeric(parts[3]),
          count     = as.numeric(parts[4]),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (is.na(total_read) || is.na(total_used)) {
    stop("Could not parse total_values_read or total_values_used_ from summary_stats.py output.")
  }
  if (nrow(bins_df) == 0) stop("No bin rows parsed from summary_stats.py output.")

  list(total_read = total_read, total_used = total_used, bins = bins_df)
}

## ------------------------ 3. Run summary_stats.py for each tip --------------

stats_by_species <- list()

for (sp in tip_species) {
  abbr <- tip_to_abbr[[sp]]
  res_file <- paste0(abbr, suffix)

  if (!file.exists(res_file)) {
    stop("Results file not found for species ", sp, " (prefix ", abbr, "): ", res_file)
  }

  cat("Reading values for", sp, "from", res_file, "...\n")

  df <- tryCatch({
    # Result files carry a '#'-prefixed header line (e.g. "#name<TAB>LTR_len<TAB>...").
    # We cannot use comment.char="#" to skip it because the name column contains a
    # '#' mid-field (e.g. "CP094635.1:4318-9747#LTR/Copia/Ale"), which would be
    # truncated. Instead drop only lines that START with '#', then parse the rest.
    lines <- readLines(res_file)
    lines <- lines[!startsWith(lines, "#")]
    read.table(text = lines, header = FALSE, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE)
  },
    error = function(e) stop("Error reading ", res_file, ": ", conditionMessage(e))
  )

  if (ncol(df) < 7) stop("File ", res_file, " has fewer than 7 columns.")

  values <- as.numeric(df[[7]])
  values <- values[is.finite(values)]
  if (length(values) == 0) stop("No numeric values in column 7 of ", res_file)

  cat("  -> ", length(values), " values read from column 7.\n", sep = "")

  stats <- run_summary_stats_for_values(
    values       = values,
    bins         = bins,
    bin_max      = bin_max,
    python_cmd   = python_cmd,
    summary_invocation = summary_invocation
  )

  cat("  total_values_read =", stats$total_read,
      ", total_values_used_ =", stats$total_used, "\n\n")

  stats_by_species[[sp]] <- stats
}

## Confirm bins are consistent across species
species1 <- tip_species[1]
bins_ref <- stats_by_species[[species1]]$bins[, c("bin_index", "bin_start", "bin_end")]

for (sp in tip_species[-1]) {
  tmp <- stats_by_species[[sp]]$bins[, c("bin_index", "bin_start", "bin_end")]
  if (!isTRUE(all.equal(bins_ref, tmp, tolerance = 1e-12))) {
    stop("Bin definitions differ between species: ", species1, " and ", sp)
  }
}

cat("All species have consistent bin definitions.\n\n")

## ------------------------ 4. Prepare trait vectors -------------------------

tot_read_vec <- setNames(numeric(length(tip_species)), tip_species)
tot_used_vec <- setNames(numeric(length(tip_species)), tip_species)

for (sp in tip_species) {
  s <- stats_by_species[[sp]]
  tot_read_vec[sp] <- s$total_read
  tot_used_vec[sp] <- s$total_used
}

bin_indices <- bins_ref$bin_index
bin_starts  <- bins_ref$bin_start
bin_ends    <- bins_ref$bin_end

## ------------------------ 5. Node identity: descendant-tip sets + plot ------

# Full, unique node identity = the set of descendant tips under that node.
# We also compute a short "signature" for readability (NOT used for identity).
get_desc_tip_names <- function(tree, node) {
  all_desc <- phytools::getDescendants(tree, node)
  tip_idx <- all_desc[all_desc <= length(tree$tip.label)]
  tree$tip.label[tip_idx]
}

make_signature <- function(tips, max_show = 6) {
  tips <- sort(unique(tips))
  n <- length(tips)
  if (n == 0) return("(?)")
  show <- head(tips, max_show)
  if (n <= max_show) {
    paste0(n, " tips: ", paste(show, collapse = ", "))
  } else {
    paste0(n, " tips: ", paste(show, collapse = ", "), ", ...")
  }
}

# Precompute node -> descendant tips (sorted) + signature + size
node_desc_tips <- vector("list", length(node_ids))
names(node_desc_tips) <- as.character(node_ids)

node_sig_map  <- setNames(character(length(node_ids)), as.character(node_ids))
node_size_map <- setNames(integer(length(node_ids)),   as.character(node_ids))
node_tips_str <- setNames(character(length(node_ids)), as.character(node_ids))

for (node in node_ids) {
  tips <- sort(unique(get_desc_tip_names(tree, node)))
  node_desc_tips[[as.character(node)]] <- tips
  node_size_map[as.character(node)] <- length(tips)
  node_sig_map[as.character(node)]  <- make_signature(tips, max_tips_sig)
  node_tips_str[as.character(node)] <- paste(tips, collapse = ";")
}

# Write a node map file for easy lookup in spreadsheets/scripts
node_map_df <- data.frame(
  node_id       = node_ids,
  clade_size    = as.integer(node_size_map[as.character(node_ids)]),
  clade_signature = unname(node_sig_map[as.character(node_ids)]),
  descendant_tips = unname(node_tips_str[as.character(node_ids)]),
  stringsAsFactors = FALSE
)

cat("Writing node map to", node_map_file, "...\n")
write.table(
  node_map_df,
  file      = node_map_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Also generate a PDF labeling internal node IDs directly on the tree
cat("Writing labeled tree PDF to", tree_pdf_file, "...\n")
pdf(tree_pdf_file, width = 10, height = 8)
plot(tree, cex = 0.8, no.margin = TRUE)
nodelabels(text = node_ids, node = node_ids, frame = "none", cex = 0.7)
tiplabels(frame = "none", cex = 0.7)
title(main = "Tree with internal node IDs (matches fastAnc node numbering)")
dev.off()

cat("\nNode identity outputs created:\n")
cat("  Node map TSV :", node_map_file, "\n")
cat("  Tree PDF     :", tree_pdf_file, "\n\n")

## ------------------------ 6. Helper: run fastAnc on a trait ----------------

reconstruct_trait <- function(tree, trait_vec, stat_type,
                              bin_index = NA_integer_,
                              bin_start = NA_real_,
                              bin_end   = NA_real_,
                              node_sig_map,
                              node_size_map,
                              node_tips_str) {

  anc <- fastAnc(tree, trait_vec, vars = TRUE, CI = TRUE)
  node_vec <- as.integer(names(anc$ace))

  data.frame(
    stat_type        = stat_type,
    node             = node_vec,
    clade_size       = as.integer(node_size_map[as.character(node_vec)]),
    clade_signature  = unname(node_sig_map[as.character(node_vec)]),
    descendant_tips  = unname(node_tips_str[as.character(node_vec)]),
    bin_index        = rep(bin_index, length(node_vec)),
    bin_start        = rep(bin_start, length(node_vec)),
    bin_end          = rep(bin_end, length(node_vec)),
    est              = anc$ace,
    CI_lower         = anc$CI95[, 1],
    CI_upper         = anc$CI95[, 2],
    stringsAsFactors = FALSE
  )
}

## ------------------------ 7. Reconstruct totals ----------------------------

cat("Reconstructing ancestral total_values_read...\n")
totals_read_res <- reconstruct_trait(
  tree         = tree,
  trait_vec    = tot_read_vec,
  stat_type    = "total_values_read",
  node_sig_map = node_sig_map,
  node_size_map = node_size_map,
  node_tips_str = node_tips_str
)

cat("Reconstructing ancestral total_values_used...\n")
totals_used_res <- reconstruct_trait(
  tree         = tree,
  trait_vec    = tot_used_vec,
  stat_type    = "total_values_used",
  node_sig_map = node_sig_map,
  node_size_map = node_size_map,
  node_tips_str = node_tips_str
)

totals_res <- rbind(totals_read_res, totals_used_res)

## ------------------------ 8. Reconstruct per-bin counts --------------------

cat("Reconstructing ancestral bin counts for", length(bin_indices), "bins...\n")

bin_results_list <- list()

for (i in seq_along(bin_indices)) {
  idx     <- bin_indices[i]
  start_i <- bin_starts[i]
  end_i   <- bin_ends[i]

  cat("  Bin", idx, " [", start_i, ",", end_i, "]: reconstructing...\n")

  counts_vec <- setNames(numeric(length(tip_species)), tip_species)
  for (sp in tip_species) {
    s <- stats_by_species[[sp]]$bins
    row_i <- which(s$bin_index == idx)
    if (length(row_i) != 1) {
      stop("Unexpected bin_index lookup for species ", sp, " and bin ", idx)
    }
    counts_vec[sp] <- s$count[row_i]
  }

  bin_res <- reconstruct_trait(
    tree          = tree,
    trait_vec     = counts_vec,
    stat_type     = "bin_count",
    bin_index     = idx,
    bin_start     = start_i,
    bin_end       = end_i,
    node_sig_map  = node_sig_map,
    node_size_map = node_size_map,
    node_tips_str = node_tips_str
  )

  bin_results_list[[length(bin_results_list) + 1]] <- bin_res
}

bins_res <- do.call(rbind, bin_results_list)

## ------------------------ 9. Write outputs ---------------------------------

cat("\nWriting ancestral totals to", totals_outfile, "...\n")
write.table(
  totals_res,
  file      = totals_outfile,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("Writing ancestral bin summaries to", bins_outfile, "...\n")
write.table(
  bins_res,
  file      = bins_outfile,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("\nDone.\n")
cat("  Tree           :", newick_file, "\n")
cat("  Suffix         :", suffix, "\n")
cat("  summary_stats  :", summary_label, "\n")
if (!is.na(abrev_file) && abrev_file != "") cat("  Abbrev mapping :", abrev_file, "\n")
cat("  Node map TSV   :", node_map_file, "\n")
cat("  Tree PDF       :", tree_pdf_file, "\n")
cat("  Totals file    :", totals_outfile, "\n")
cat("  Bins file      :", bins_outfile, "\n")
