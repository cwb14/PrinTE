#!/usr/bin/env Rscript

# ancestral_reconstruction.R
# Run summary_stats3.py for each tip, then estimate ancestral summary stats
# (totals and per-bin counts) at each internal node using fastAnc.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

## ------------------------ 0. Config ----------------------------------------

# Tree file (dated)
tree_file <- "iqtree_rerooted.dated.nwk"

# Python command and summary stats script
python_cmd      <- "python"             # change to "python3" if needed
summary_script  <- "summary_stats.py"

# Parameters for summary_stats3.py
bins    <- 50
bin_max <- 0.15

# Mapping from abbreviations to full species names
abbr_to_full <- c(
  Aaren = "Arabidopsis arenosa",
  Ahall = "Arabidopsis halleri",
  Alyra = "Arabidopsis lyrata",
  Chisp = "Camelina hispida",
  Athal = "Arabidopsis thaliana"
)

# File name pattern for LTR alignment results (per abbreviation)
# e.g. "Athal.LTRs.alns.results"
results_pattern <- "%s.LTRs.alns.results"

# Output files
totals_outfile <- "ancestral_summary_totals.tsv"
bins_outfile   <- "ancestral_summary_bins.tsv"

## ------------------------ 1. Read tree -------------------------------------

tree <- read.tree(tree_file)

cat("Read tree with", length(tree$tip.label), "tips and",
    tree$Nnode, "internal nodes.\n\n")

## Handle abbreviated tip labels -> full species names, if needed
if (all(tree$tip.label %in% names(abbr_to_full))) {
  tree$tip.label <- abbr_to_full[tree$tip.label]
  cat("Converted abbreviated tip labels to full species names.\n")
}

cat("Tip labels in tree:\n")
print(tree$tip.label)
cat("\n")

# Make reverse mapping: full -> abbrev
full_to_abbr <- setNames(names(abbr_to_full), abbr_to_full)

if (!all(tree$tip.label %in% names(full_to_abbr))) {
  stop("Some tip labels in the tree do not have known abbreviations in full_to_abbr.")
}

## ------------------------ 2. Helper: run summary_stats3.py -----------------

run_summary_stats_for_values <- function(values,
                                         bins,
                                         bin_max,
                                         python_cmd,
                                         summary_script) {
  # Write values (one per line) to a temporary file to emulate:
  #   --input <(cat file | cut -f 7)
  tmpfile <- tempfile(pattern = "summary_input_", fileext = ".txt")
  on.exit(unlink(tmpfile), add = TRUE)

  # Write numeric values as plain text, one per line
  writeLines(format(values, digits = 10, scientific = FALSE), con = tmpfile)

  # Build command
  cmd <- sprintf(
    '%s %s --bins %d --bin_max %g --input %s',
    shQuote(python_cmd),
    shQuote(summary_script),
    bins,
    bin_max,
    shQuote(tmpfile)
  )

  # Run and capture stdout
  out <- tryCatch(
    system(cmd, intern = TRUE),
    error = function(e) {
      stop("Error running summary_stats3.py: ", conditionMessage(e))
    }
  )

  if (length(out) == 0) {
    stop("summary_stats3.py produced no output.")
  }

  # Parse output
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
      # header lines
      if (grepl("^#\\s*total_values_read", line)) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) {
          total_read <- as.numeric(parts[2])
        }
      } else if (grepl("^#\\s*total_values_used_", line)) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 2) {
          total_used <- as.numeric(parts[2])
        }
      }
    } else {
      # data lines: bin_index, bin_start, bin_end, count
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
    stop("Could not parse total_values_read or total_values_used_ from summary_stats3 output.")
  }

  if (nrow(bins_df) == 0) {
    stop("No bin rows parsed from summary_stats3 output.")
  }

  list(
    total_read = total_read,
    total_used = total_used,
    bins       = bins_df
  )
}

## ------------------------ 3. Run summary_stats3.py for each tip ------------

tip_species <- tree$tip.label

stats_by_species <- list()

for (sp in tip_species) {
  abbr <- full_to_abbr[[sp]]
  res_file <- sprintf(results_pattern, abbr)

  if (!file.exists(res_file)) {
    stop("Results file not found for species ", sp, ": ", res_file)
  }

  cat("Reading values for", sp, "from", res_file, "...\n")

  # Read the .results file; assume tab-delimited, >= 7 columns, no header
  df <- tryCatch(
    read.table(res_file, header = FALSE, sep = "\t", comment.char = ""),
    error = function(e) {
      stop("Error reading ", res_file, ": ", conditionMessage(e))
    }
  )

  if (ncol(df) < 7) {
    stop("File ", res_file, " has fewer than 7 columns.")
  }

  values <- df[[7]]
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    stop("No numeric values in column 7 of ", res_file)
  }

  cat("  -> ", length(values), "values read from column 7.\n", sep = "")

  stats <- run_summary_stats_for_values(
    values       = values,
    bins         = bins,
    bin_max      = bin_max,
    python_cmd   = python_cmd,
    summary_script = summary_script
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
  if (!all.equal(bins_ref, tmp, tolerance = 1e-12)) {
    stop("Bin definitions differ between species: ", species1, " and ", sp)
  }
}

cat("All species have consistent bin definitions.\n\n")

## ------------------------ 4. Prepare trait vectors -------------------------

# total_values_read and total_values_used for each tip
tot_read_vec <- setNames(numeric(length(tip_species)), tip_species)
tot_used_vec <- setNames(numeric(length(tip_species)), tip_species)

for (sp in tip_species) {
  s <- stats_by_species[[sp]]
  tot_read_vec[sp] <- s$total_read
  tot_used_vec[sp] <- s$total_used
}

# Per-bin counts: we will reconstruct separate traits for each bin_index
bin_indices <- bins_ref$bin_index
bin_starts  <- bins_ref$bin_start
bin_ends    <- bins_ref$bin_end

## ------------------------ 5. Helper: representative tips per node ----------

get_desc_tip_names <- function(tree, node) {
  all_desc <- phytools::getDescendants(tree, node)
  # tips are indices <= Ntip
  tip_idx <- all_desc[all_desc <= length(tree$tip.label)]
  tree$tip.label[tip_idx]
}

# Precompute representative tips for all internal nodes
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
    stat_type          = stat_type,
    node               = node_vec,
    representative_tips = rep_tip_map[as.character(node_vec)],
    bin_index          = rep(bin_index, length(node_vec)),
    bin_start          = rep(bin_start, length(node_vec)),
    bin_end            = rep(bin_end, length(node_vec)),
    est                = anc$ace,
    CI_lower           = anc$CI95[, 1],
    CI_upper           = anc$CI95[, 2],
    stringsAsFactors   = FALSE
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
  idx <- bin_indices[i]
  start_i <- bin_starts[i]
  end_i   <- bin_ends[i]

  cat("  Bin", idx, " [", start_i, ",", end_i, "]: reconstructing...\n")

  # Trait vector of counts for this bin across species
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
cat("  Totals file:", totals_outfile, "\n")
cat("  Bins file  :", bins_outfile, "\n")
