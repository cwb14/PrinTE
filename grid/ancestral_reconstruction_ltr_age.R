#!/usr/bin/env Rscript

# ancestral_reconstruction.R
# For each tip taxon, run summary_stats.py on column 7 from:
#   <taxon><suffix>
# Then reconstruct ancestral totals + per-bin counts at internal nodes via fastAnc.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

## ------------------------ 0. Defaults --------------------------------------

# Defaults (override via CLI args)
default_newick      <- "iqtree_rerooted.dated.nwk"
default_suffix      <- ".LTRs.alns.results"
default_bins        <- 50
default_bin_max     <- 0.15
default_python_cmd  <- "python"

totals_outfile <- "ancestral_summary_totals.tsv"
bins_outfile   <- "ancestral_summary_bins.tsv"

## ------------------------ 0b. CLI parsing ----------------------------------

print_usage_and_exit <- function(exit_code = 0) {
  cat(
"Usage:
  ancestral_reconstruction.R [--newick FILE] [--suffix SUFFIX]
                             [--summary_path PATH] [--abrev FILE]
                             [--bins INT] [--bin_max NUM]
                             [--python CMD]

Required inputs:
  --newick FILE
      Newick tree (tips must match result-file prefixes unless --abrev is used).
      Default: iqtree_rerooted.dated.nwk

  --suffix SUFFIX
      Results-file suffix appended to each tip label (or abbreviation).
      Example: .LTRs.alns.results
      Default: .LTRs.alns.results

Optional:
  --summary_path PATH
      Path to summary_stats.py. If not provided, the script will look for
      summary_stats.py in the same directory as this R script.

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

Examples:
  ancestral_reconstruction.R --suffix .LTRs.alns.results --newick subset.nwk

  ancestral_reconstruction.R --suffix .LTRs.alns.results --newick subset.nwk \\
    --summary_path ../path/to/summary_stats.py

  ancestral_reconstruction.R --newick subset.nwk --suffix .LTRs.alns.results \\
    --abrev abrev.tsv
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

newick_file   <- get_arg_value("--newick",       args, default_newick)
suffix        <- get_arg_value("--suffix",       args, default_suffix)
summary_path  <- get_arg_value("--summary_path", args, NA_character_)
abrev_file    <- get_arg_value("--abrev",        args, NA_character_)
bins          <- as.integer(get_arg_value("--bins",   args, as.character(default_bins)))
bin_max       <- as.numeric(get_arg_value("--bin_max",args, as.character(default_bin_max)))
python_cmd    <- get_arg_value("--python",       args, default_python_cmd)

if (is.na(bins) || bins <= 0) stop("--bins must be a positive integer.")
if (is.na(bin_max) || bin_max <= 0) stop("--bin_max must be a positive number.")
if (!file.exists(newick_file)) stop("Newick file not found: ", newick_file)

## ------------------------ 0c. Resolve script directory ----------------------

# Find where THIS R script lives (best effort). If invoked in weird ways, fall back to getwd().
get_script_dir <- function() {
  # Common, robust approach for Rscript:
  # Look for the --file=... argument in commandArgs()
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE))
  }
  # Fallback
  return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

script_dir <- get_script_dir()

## ------------------------ 0d. Resolve summary_stats.py ----------------------

if (is.na(summary_path) || summary_path == "") {
  summary_path <- file.path(script_dir, "summary_stats.py")
}
if (!file.exists(summary_path)) {
  stop(
    "summary_stats.py not found at: ", summary_path, "\n",
    "Either place summary_stats.py in the same directory as this R script (", script_dir, "),\n",
    "or pass --summary_path /path/to/summary_stats.py"
  )
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
tip_to_abbr <- setNames(tree$tip.label, tree$tip.label)  # names are "analysis label", values are file-prefix abbr

if (!is.na(abrev_file) && abrev_file != "") {
  if (!file.exists(abrev_file)) stop("Abbreviation mapping file not found: ", abrev_file)

  ab <- tryCatch(
    read.table(abrev_file, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE),
    error = function(e) stop("Error reading --abrev file: ", conditionMessage(e))
  )
  if (ncol(ab) < 2) stop("--abrev file must have at least 2 columns: <abbrev> <full_name>")

  abbr_to_full <- setNames(ab[[2]], ab[[1]])

  # If tree tips are abbreviations present in mapping, rename tips to full names for analysis/output
  if (all(tree$tip.label %in% names(abbr_to_full))) {
    old_tips <- tree$tip.label
    tree$tip.label <- unname(abbr_to_full[tree$tip.label])

    # Now we need a mapping from analysis label (full) -> file prefix (abbrev)
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

## ------------------------ 2. Helper: run summary_stats.py -------------------

run_summary_stats_for_values <- function(values,
                                         bins,
                                         bin_max,
                                         python_cmd,
                                         summary_path) {
  tmpfile <- tempfile(pattern = "summary_input_", fileext = ".txt")
  on.exit(unlink(tmpfile), add = TRUE)

  writeLines(format(values, digits = 10, scientific = FALSE), con = tmpfile)

  cmd <- sprintf(
    '%s %s --bins %d --bin_max %g --input %s',
    shQuote(python_cmd),
    shQuote(summary_path),
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

  df <- tryCatch(
    read.table(res_file, header = FALSE, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE),
    error = function(e) stop("Error reading ", res_file, ": ", conditionMessage(e))
  )

  if (ncol(df) < 7) stop("File ", res_file, " has fewer than 7 columns.")

  values <- as.numeric(df[[7]])
  values <- values[is.finite(values)]
  if (length(values) == 0) stop("No numeric values in column 7 of ", res_file)

  cat("  -> ", length(values), " values read from column 7.\n", sep = "")

  stats <- run_summary_stats_for_values(
    values      = values,
    bins        = bins,
    bin_max     = bin_max,
    python_cmd  = python_cmd,
    summary_path = summary_path
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

## ------------------------ 5. Representative tips per node ------------------

get_desc_tip_names <- function(tree, node) {
  all_desc <- phytools::getDescendants(tree, node)
  tip_idx <- all_desc[all_desc <= length(tree$tip.label)]
  tree$tip.label[tip_idx]
}

Ntip  <- length(tree$tip.label)
Nnode <- tree$Nnode
node_ids <- (Ntip + 1):(Ntip + Nnode)

rep_tip_map <- setNames(character(length(node_ids)), as.character(node_ids))
for (node in node_ids) {
  tips <- get_desc_tip_names(tree, node)
  if (length(tips) >= 2) {
    rep_tip_map[as.character(node)] <- paste0("(", tips[1], ", ", tips[2], ")")
  } else if (length(tips) == 1) {
    rep_tip_map[as.character(node)] <- paste0("(", tips[1], ")")
  } else {
    rep_tip_map[as.character(node)] <- "(?)"
  }
}

## ------------------------ 6. Helper: run fastAnc on a trait ----------------

reconstruct_trait <- function(tree, trait_vec, stat_type,
                              bin_index = NA_integer_,
                              bin_start = NA_real_,
                              bin_end   = NA_real_,
                              rep_tip_map) {

  anc <- fastAnc(tree, trait_vec, vars = TRUE, CI = TRUE)
  node_vec <- as.integer(names(anc$ace))

  data.frame(
    stat_type           = stat_type,
    node                = node_vec,
    representative_tips = rep_tip_map[as.character(node_vec)],
    bin_index           = rep(bin_index, length(node_vec)),
    bin_start           = rep(bin_start, length(node_vec)),
    bin_end             = rep(bin_end, length(node_vec)),
    est                 = anc$ace,
    CI_lower            = anc$CI95[, 1],
    CI_upper            = anc$CI95[, 2],
    stringsAsFactors    = FALSE
  )
}

## ------------------------ 7. Reconstruct totals ----------------------------

cat("Reconstructing ancestral total_values_read...\n")
totals_read_res <- reconstruct_trait(
  tree        = tree,
  trait_vec   = tot_read_vec,
  stat_type   = "total_values_read",
  rep_tip_map = rep_tip_map
)

cat("Reconstructing ancestral total_values_used...\n")
totals_used_res <- reconstruct_trait(
  tree        = tree,
  trait_vec   = tot_used_vec,
  stat_type   = "total_values_used",
  rep_tip_map = rep_tip_map
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
    tree        = tree,
    trait_vec   = counts_vec,
    stat_type   = "bin_count",
    bin_index   = idx,
    bin_start   = start_i,
    bin_end     = end_i,
    rep_tip_map = rep_tip_map
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
cat("  summary_stats  :", summary_path, "\n")
if (!is.na(abrev_file) && abrev_file != "") cat("  Abbrev mapping :", abrev_file, "\n")
cat("  Totals file    :", totals_outfile, "\n")
cat("  Bins file      :", bins_outfile, "\n")
