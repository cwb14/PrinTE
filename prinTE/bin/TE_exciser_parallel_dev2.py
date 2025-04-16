#!/usr/bin/env python3
"""
Script to simulate TE excisions from a genome FASTA based on a BED file annotation.
It randomly removes features (TEs only) and adjusts the BED and genome sequence accordingly.
See the help text below for details.

New functionality: If a TE chosen for excision is one of:
   NEST_TE_IN_TE, NEST_TE_IN_GENE, or NON_NEST_GROUP_TE,
and its feature_ID contains "LTRlen" (and does NOT contain "_SOLO"),
then with a probability given by --soloLTR_freq (in percent) the excision is done partially.
In a partial excision the sequence from (START + LTRlen) to END is removed, and the BED
entry is updated so that its end becomes (START+LTRlen) and "_SOLO" is appended to the feature_ID,
while preserving supplemental info fields.

Additional filtering is applied to ensure that non-INTACT_LTR candidates are excluded:
  1) If the first supplemental attribute in NAME contains "CUT_BY", the candidate is not considered.
  2) For candidates with supplemental info fields, if any other entry within 100 lines (based on BED order)
     has an identical TSD, identical strand, and a NAME that is an exact or prefix match,
     the candidate is filtered out.
     
Note: Consolidation of nest groups is reserved for full-length excisions.
  
NEW EDIT:
For the non-fixed rate model (specified with "--rate"), the number of TE excisions is still determined
by: rate * generations * (sum of disrupted gene weights). However, when selecting excision events,
the sampling weight is now modified so that gene-disrupting TEs (those classified as NEST_TE_IN_GENE)
are further favored. Their base weight (based on feature length and parameter k) is multiplied by:
    (1 + sel_coeff * (gene_weight) * generations)
where gene_weight is obtained from gene_selection.tsv (or defaults to 1), and sel_coeff is provided by --sel_coeff.
This adjustment reflects selection acting against individuals that fail to delete gene-disrupting TEs.
  
Additionally, if a fixed excision rate is provided with --fix_ex, the gene weights and selection coefficient are bypassed.
  
The script also produces publication-quality figures for the log-normal distribution and the weighted candidate selection curve.

NEW CHROMATIN BIAS:
In rate mode, the likelihood of a TE excision is further modulated by its chromatin context.
A new parameter --euch_het_buffer specifies a buffer (in bp) around gene features to be treated as euchromatin.
Any TE whose midpoint falls in one of these regions will have its base weight multiplied
by --euch_het_bias (e.g. 1.1 increases the weight by 10%). Regions outside these merged gene intervals
are assumed to be heterochromatin and receive no additional boost.

MULTIPROCESSING EDIT:
A new option -m [int] specifies the maximum number of chromosomes to process in parallel.
For example, '-m 10' processes up to 10 chromosomes concurrently.
Each chromosome gets a unique deterministic seed (global seed + index) so that results are reproducible.
Multiprocessing.Pool.map() is used to parallelize the genome update step.
"""

import argparse
import sys
import random
import copy
import os
import math
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

# =============================================================================
# Data Structures and Helpers
# =============================================================================

class BedEntry:
    def __init__(self, chrom, start, end, name, strand, tsd, lineno):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.strand = strand
        self.tsd = tsd
        self.lineno = lineno   # original order line number (for debugging)
        self.group = None      # if part of a nest group, reference to the group id
        self.subtype = None    # classification (see below)
        # Derived attributes:
        # feature_ID: first attribute before ';'
        self.feature_id = name.split(';')[0]
        # supplemental info (list of attributes after first)
        self.supp = name.split(';')[1:] if ';' in name else []
    
    def __str__(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), self.name, self.strand, self.tsd])
    
    def length(self):
        return self.end - self.start

def parse_bed(bed_file):
    """Parse the BED file and return a list of BedEntry objects."""
    entries = []
    with open(bed_file) as f:
        for lineno, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                sys.exit(f"Error: Line {lineno+1} in BED file does not have 6 columns.")
            entry = BedEntry(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], lineno)
            entries.append(entry)
    return entries
    
def _select_removals_per_chrom(task):
    """
    Wrapper for per-chromosome select_removals.
    task: tuple of (chrom, entries, nest_groups, num_excision, seed, k, gene_weights, generations,
                    sel_coeff, euchromatin_intervals, euch_het_bias)
    Returns: set of selected BedEntry objects for this chromosome
    """
    chrom, entries, nest_groups, num_excision, seed, k, gene_weights, generations, sel_coeff, euch_intervals, euch_het_bias = task
    return select_removals(
        entries, nest_groups, num_excision, seed, k, gene_weights, generations, sel_coeff,
        euch_intervals=euch_intervals, euch_het_bias=euch_het_bias
    )

def parse_fasta(fasta_file):
    """Parse FASTA file and return a dict: {chrom: SeqRecord}."""
    recs = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        recs[rec.id] = rec
    return recs

def build_euchromatin_intervals(entries, buffer):
    """
    Build a dictionary of euchromatin intervals from gene entries.
    For each gene (feature_ID starts with "gene"), extend its interval by the specified buffer.
    Overlapping intervals on the same chromosome are merged.
    Returns a dict: {chrom: [(start, end), ...]}.
    """
    intervals_by_chrom = defaultdict(list)
    for e in entries:
        if e.feature_id.startswith("gene"):
            start = max(0, e.start - buffer)
            end = e.end + buffer
            intervals_by_chrom[e.chrom].append((start, end))
    merged = {}
    for chrom, intervals in intervals_by_chrom.items():
        if not intervals:
            continue
        intervals.sort(key=lambda x: x[0])
        merged_list = []
        current_start, current_end = intervals[0]
        for s, e_val in intervals[1:]:
            if s <= current_end:  # overlapping intervals
                current_end = max(current_end, e_val)
            else:
                merged_list.append((current_start, current_end))
                current_start, current_end = s, e_val
        merged_list.append((current_start, current_end))
        merged[chrom] = merged_list
    return merged

# =============================================================================
# Classification of BED Entries
# =============================================================================

def classify_entries(entries):
    """
    Classify entries into groups.
      - Genes are those with feature_ID starting with 'gene'
      - A gene is 'disrupted' if its NAME contains a semicolon.
    Then identify nest groups.
    
    Nest groups occur when a BED entry (the middle one) has a supplemental attribute
    starting with 'NESTED_IN:' AND is immediately flanked (previous and next entries)
    by lines that share identical feature_ID (only the first part, ignoring supplemental info),
    strand and TSD.
    
    If the flanking feature_ID begins with 'gene', the nest group is a NEST_GROUP_GENE;
    otherwise it is a NEST_GROUP_TE.
    
    For nest groups, assign the following subtypes:
      - For TE nest groups:
           CUT_PAIR_TE_1, NEST_TE_IN_TE, CUT_PAIR_TE_2.
      - For gene nest groups:
           CUT_PAIR_GENE_1, NEST_TE_IN_GENE, CUT_PAIR_GENE_2.
    All other entries with feature_ID starting with 'gene' become NON_NEST_GROUP_GENES.
    All remaining entries become NON_NEST_GROUP_TE.
    """
    n = len(entries)
    i = 0
    used = [False]*n
    nest_groups = {}
    group_id = 0
    
    while i < n:
        if i > 0 and i < n-1 and not used[i-1] and not used[i] and not used[i+1]:
            middle = entries[i]
            if middle.supp and middle.supp[0].startswith("NESTED_IN:"):
                prev = entries[i-1]
                nxt = entries[i+1]
                # Compare based on feature_ID (first attribute) only.
                if (prev.feature_id == nxt.feature_id and
                    prev.strand == nxt.strand and
                    prev.tsd == nxt.tsd):
                    group_id += 1
                    if prev.feature_id.startswith("gene"):
                        prev.subtype = "CUT_PAIR_GENE_1"
                        middle.subtype = "NEST_TE_IN_GENE"
                        nxt.subtype = "CUT_PAIR_GENE_2"
                    else:
                        prev.subtype = "CUT_PAIR_TE_1"
                        middle.subtype = "NEST_TE_IN_TE"
                        nxt.subtype = "CUT_PAIR_TE_2"
                    for ent in (prev, middle, nxt):
                        ent.group = group_id
                        used[ent.lineno] = True
                    nest_groups[group_id] = (prev, middle, nxt)
                    i += 2
                    continue
        i += 1

    for ent in entries:
        if ent.subtype is None:
            if ent.feature_id.startswith("gene"):
                ent.subtype = "NON_NEST_GROUP_GENES"
            else:
                ent.subtype = "NON_NEST_GROUP_TE"
    return entries, nest_groups

# =============================================================================
# Calculate number of TE excisions to make
# =============================================================================

def calculate_excision_count(entries, rate, generations, gene_weights):
    """
    For disrupted genes, use the disrupted_gene_weight from gene_weights instead of counting each as 1.
    A gene is considered disrupted if its feature_ID starts with 'gene' and its NAME contains a semicolon.
    Sum the weights of all unique disrupted genes and calculate:
         lambda_val = rate * generations * (sum of disrupted gene weights)
    Then sample the number of TE excisions from a Poisson distribution with lambda_val.
    """
    disrupted_genes = {e.feature_id for e in entries if e.feature_id.startswith("gene") and (';' in e.name)}
    weight_sum = sum(gene_weights.get(gene, 1) for gene in disrupted_genes)
    print(f"Disrupted genes (unique): {len(disrupted_genes)}")
    print(f"Sum of disrupted gene weights: {weight_sum:.4f}")
    
    lambda_val = rate * generations * weight_sum
    excision_count = np.random.poisson(lambda_val)
    print(f"Lambda for Poisson (rate mode): {lambda_val:.4f}")
    print(f"Number of TE excisions to simulate: {excision_count}")
    return excision_count

# =============================================================================
# Select Removal Events (Weighted by length, chromatin state and selection for gene-disrupting TEs)
# =============================================================================

from collections import defaultdict
from multiprocessing import Pool

def select_removals(entries, nest_groups, num_excision, seed, k, gene_weights, generations, sel_coeff,
                    euch_intervals=None, euch_het_bias=1.0):
    """
    Eligible for removal are:
       - NON_NEST_GROUP_TE entries,
       - In nest groups of type TE: the middle entry (NEST_TE_IN_TE) is eligible,
         and flanking entries (CUT_PAIR_TE_1 and CUT_PAIR_TE_2) are eligible.
       - In nest groups of type GENE: only the middle (NEST_TE_IN_GENE) is eligible.
       
    Genes (NON_NEST_GROUP_GENES, CUT_PAIR_GENE_1, CUT_PAIR_GENE_2) are immune.
    
    This version selects from eligible features using a weighted probability distribution based on length.
    For each eligible entry, let L = end - start, and Lmax be the maximum L among eligible entries.
    Then assign the base excision weight as:
         weight_base = exp( - k * (1 - (L / Lmax)) )
    Additionally, if euch_intervals is provided (rate mode), and the candidate’s midpoint falls within any
    euchromatin interval, multiply the weight by euch_het_bias.
    
    NEW EDIT: For candidates that are gene-disrupting (subtype NEST_TE_IN_GENE), their weight is multiplied
    by an additional factor: (1 + sel_coeff * gene_weight * generations), where gene_weight is obtained from gene_weights
    using the gene ID from the flanking entry in the nest group.
    
    Removal events are then sampled without replacement using these weights.
    
    For nest groups, if a candidate is chosen from a group, then other members from that group are skipped.
    """
    random.seed(seed)
    # Build eligible candidate list
    eligible = []
    for e in entries:
        if e.subtype in ["NON_NEST_GROUP_TE", "CUT_PAIR_TE_1", "CUT_PAIR_TE_2",
                         "NEST_TE_IN_TE", "NEST_TE_IN_GENE"]:
            eligible.append(e)
    
    if not eligible:
        return set()
    
    # Determine Lmax among eligible entries.
    Lmax = max(e.length() for e in eligible)
    
    # Pre-calculate weights for each eligible candidate.
    weights = {}
    for e in eligible:
        base_weight = math.exp(-k * (1 - (e.length() / Lmax)))
        # Apply chromatin state bias if euch_intervals is provided.
        if euch_intervals is not None:
            # Use the midpoint of the TE candidate.
            mid = (e.start + e.end) // 2
            if e.chrom in euch_intervals:
                for (istart, iend) in euch_intervals[e.chrom]:
                    if istart <= mid <= iend:
                        base_weight *= euch_het_bias
                        break
        # If the candidate is gene-disrupting, further weight it.
        if e.subtype == "NEST_TE_IN_GENE" and e.group is not None:
            # For nest groups with gene disruption, use the gene ID from one of the flanking entries.
            group = nest_groups.get(e.group)
            if group:
                # Assume the first flanking entry (CUT_PAIR_GENE_1) has the gene ID.
                gene_id = group[0].feature_id
                gene_weight = gene_weights.get(gene_id, 1)
                multiplier = 1 + sel_coeff * gene_weight * generations
                base_weight *= multiplier
        weights[id(e)] = base_weight

    removals = set()
    group_removed = {}
    
    # Perform weighted sampling without replacement.
    # Continue until we have num_excision events or no eligible candidates remain.
    while len(removals) < num_excision and eligible:
        total_weight = sum(weights[id(e)] for e in eligible)
        # Pick a random threshold in [0, total_weight)
        r = random.uniform(0, total_weight)
        cumulative = 0.0
        chosen = None
        for e in eligible:
            cumulative += weights[id(e)]
            if cumulative >= r:
                chosen = e
                break
        if chosen is None:
            chosen = eligible[-1]
        
        # If this candidate belongs to a nest group, check if one from the group is already selected.
        if chosen.group is not None:
            if chosen.group in group_removed:
                eligible = [e for e in eligible if e is not chosen]
                continue
            else:
                if chosen.subtype in ["NEST_TE_IN_TE", "NEST_TE_IN_GENE"]:
                    group_removed[chosen.group] = "middle"
                else:
                    group_removed[chosen.group] = "flank"
                # Remove all candidates from this group from eligible.
                eligible = [e for e in eligible if e.group != chosen.group]
        else:
            eligible = [e for e in eligible if e is not chosen]
        
        removals.add(chosen)
    print(f"Selected {len(removals)} removal events.")
    return removals

def parallel_select_removals(entries, nest_groups, rate, generations, gene_weights,
                              sel_coeff, euch_intervals, euch_het_bias, seed, k, max_proc,
                              fix_ex=None, genome_size=None):
    """
    Parallel version of select_removals, operating per chromosome.
    Supports both variable-rate (default) and fixed-rate (--fix_ex) modes.

    Parameters:
        - fix_ex: If provided, use fixed excision rate model.
        - genome_size: Required if fix_ex is provided.
    """
    # Partition entries and nest groups by chromosome
    entries_by_chrom = defaultdict(list)
    nest_groups_by_chrom = defaultdict(dict)
    for e in entries:
        entries_by_chrom[e.chrom].append(e)
    for gid, group in nest_groups.items():
        chrom = group[1].chrom
        nest_groups_by_chrom[chrom][gid] = group

    excision_counts = {}

    if fix_ex is not None:
        # Fixed-rate mode: use genome size to determine number of excisions
        if genome_size is None:
            raise ValueError("genome_size must be provided when using --fix_ex.")
        lambda_val = fix_ex * genome_size * generations
        total_excision_count = np.random.poisson(lambda_val)
        print(f"[Fixed-rate mode] Global lambda: {lambda_val:.4f}, total TE excisions: {total_excision_count}")

        # Distribute excisions across chromosomes proportional to chromosome size
        chrom_lengths = {chrom: sum(e.length() for e in ents)
                         for chrom, ents in entries_by_chrom.items()}
        total_length = sum(chrom_lengths.values())
        for chrom in entries_by_chrom:
            prop = chrom_lengths[chrom] / total_length if total_length > 0 else 0
            chrom_lambda = prop * total_excision_count
            excision_counts[chrom] = np.random.poisson(chrom_lambda)
    else:
        # Rate-based mode: use disrupted gene weights
        for chrom, ents in entries_by_chrom.items():
            disrupted = {e.feature_id for e in ents
                         if e.feature_id.startswith("gene") and ";" in e.name}
            weight_sum = sum(gene_weights.get(g, 1) for g in disrupted)
            lambda_val = rate * generations * weight_sum
            excision_counts[chrom] = np.random.poisson(lambda_val)

    # Build task list
    tasks = []
    for idx, chrom in enumerate(entries_by_chrom):
        ents = entries_by_chrom[chrom]
        ngroups = nest_groups_by_chrom.get(chrom, {})
        num_exc = excision_counts.get(chrom, 0)
        chrom_seed = seed + idx
        euch = euch_intervals.get(chrom) if euch_intervals else None

        sc = sel_coeff if fix_ex is None else 0.0  # ignore selection coefficient in fixed mode

        tasks.append((chrom, ents, ngroups, num_exc, chrom_seed,
                      k, gene_weights, generations, sc, euch, euch_het_bias))

    # Run in parallel
    if max_proc > 1:
        with Pool(processes=max_proc) as pool:
            results = pool.map(_select_removals_per_chrom, tasks)
    else:
        results = map(_select_removals_per_chrom, tasks)

    removals = set()
    for r in results:
        removals.update(r)
    print(f"Total removal events selected across all chromosomes: {len(removals)}")
    return removals

# =============================================================================
# Helper for multiprocessing: Genome update per chromosome
# =============================================================================

def update_chromosome(task):
    """
    Process the genome sequence update for one chromosome.
    task: (chrom, rec, events, seed_val)
    Sets the random seed for reproducibility.
    Returns (chrom, updated SeqRecord, total_shift) where total_shift is the sum of removed bases.
    """
    chrom, rec, events, seed_val = task
    random.seed(seed_val)
    np.random.seed(seed_val)
    seq = list(str(rec.seq))
    total_shift = 0
    events_sorted = sorted(events, key=lambda x: x[0])
    for (rstart, rend, rlen) in events_sorted:
        adj_start = rstart - total_shift
        adj_end = rend - total_shift
        del seq[adj_start:adj_end]
        total_shift += rlen
    updated_seq = "".join(seq)
    new_rec = SeqRecord(Seq(updated_seq), id=rec.id, description="")
    return chrom, new_rec, total_shift

# =============================================================================
# Simulation: Remove sequences and adjust bed coordinates
# =============================================================================

def simulate_excision(genome_records, entries, nest_groups, removals, soloLTR_freq, max_proc, global_seed):
    """
    This function modifies the genome sequences and bed entries.
    
    There are two main cases:
    
    1) Simple removal (NON_NEST_GROUP_TE, or flanking entries in a nest group that get excised)
       - For a given bed entry to remove, remove from the genome the sequence from start to (end - TSD_length)
         (if TSD is not 'NA'; if TSD=='NA', remove the whole region).
       - Remove the bed entry.
       - For subsequent coordinates on that chromosome, shift by the removed length.
       
    2) Nested removal: if a nest-group middle (NEST_TE_IN_TE or NEST_TE_IN_GENE) is selected
       and it is NOT a partial (solo) excision, then:
       - Remove the sequence corresponding to the middle entry.
       - Then consolidate the flanking entries into one bed entry.
         The new bed entry will have:
             start = CUT_PAIR_x_1.start (unchanged)
             end = CUT_PAIR_x_2.end - (middle.length())
         (i.e. the gap from removal is subtracted)
       - The consolidated entry replaces the three original ones.
       
    NEW: For a removal event that qualifies as an INTACT_LTR candidate (its feature_ID contains 'LTRlen'
         and does not contain '_SOLO') and its subtype is one of [NEST_TE_IN_TE, NEST_TE_IN_GENE, NON_NEST_GROUP_TE],
         we further filter candidates:
           a) If the first supplemental attribute in the NAME contains "CUT_BY", it is not considered.
           b) If the candidate has supplemental info fields and any other entry within 100 lines (before or after)
              has identical TSD and strand and a NAME that is an exact or prefix match, it is not considered.
         For those candidates passing these filters, with probability soloLTR_freq (in percent) a partial excision is performed.
         In a partial excision, the sequence from (start + LTRlen) to end is removed, and the BED entry is updated so that
         end becomes (start+LTRlen) and '_SOLO' is appended to the feature_ID while preserving supplemental info.
         Note: No consolidation is performed after a partial excision.
         
    MULTIPROCESSING: Genome sequences are processed per chromosome concurrently using a Pool.
    The global seed is used to assign each chromosome a unique seed for reproducibility.
    """
    def qualifies_as_intact_ltr(e, all_entries):
        # Condition (a): if the first supplemental attribute contains "CUT_BY", reject.
        if e.supp and "CUT_BY" in e.supp[0]:
            return False
        # Condition (b): if there is at least one supplemental field, check nearby entries.
        if e.supp:
            for other in all_entries:
                if other is e:
                    continue
                if abs(other.lineno - e.lineno) <= 100:
                    if other.tsd == e.tsd and other.strand == e.strand:
                        # Check if the NAMEs have an exact or prefix match.
                        if e.name.startswith(other.name) or other.name.startswith(e.name):
                            return False
        return True

    # Iterate over removals in a deterministic order.
    sorted_removals = sorted(removals, key=lambda x: (x.chrom, x.start, x.end, x.lineno))
    
    # Determine which removal events will be partial excisions.
    # Map entry id -> LTR offset (extracted from "LTRlen:XXX") for those chosen for partial excision.
    partial_info = {}
    for e in sorted_removals:
        if e.subtype in ["NEST_TE_IN_TE", "NEST_TE_IN_GENE", "NON_NEST_GROUP_TE"]:
            if ("LTRlen" in e.feature_id) and ("_SOLO" not in e.feature_id):
                if not qualifies_as_intact_ltr(e, entries):
                    continue
                m = re.search(r"LTRlen:(\d+)", e.feature_id)
                if m:
                    ltr_val = int(m.group(1))
                    if random.random() < (soloLTR_freq / 100.0):
                        partial_info[id(e)] = ltr_val

    removals_by_chrom = defaultdict(list)
    new_entries = []
    to_remove = set(sorted_removals)
    nest_middle_removed = set()
    for e in sorted_removals:
        if e.subtype in ["NEST_TE_IN_TE", "NEST_TE_IN_GENE"]:
            if id(e) not in partial_info:
                nest_middle_removed.add(e.group)
    
    # Process nest groups.
    for gid, group in nest_groups.items():
        if gid in nest_middle_removed:
            rem_len = group[1].length()
            chrom = group[1].chrom
            removals_by_chrom[chrom].append( (group[1].start, group[1].end, rem_len) )
            print(f"Excision (nest group consolidation): Group {gid} - Removing middle element {group[1].chrom}:{group[1].start}-{group[1].end} (length {rem_len})")
            new_start = group[0].start
            new_end = group[2].end - rem_len
            new_name = group[0].feature_id  
            new_strand = group[0].strand
            new_tsd = group[0].tsd
            new_entry = BedEntry(chrom, new_start, new_end, new_name, new_strand, new_tsd, -1)
            new_entry.subtype = "CONSOLIDATED_NEST"
            new_entries.append(new_entry)
        else:
            for e in group:
                if e in to_remove:
                    if id(e) in partial_info:
                        ltr_val = partial_info[id(e)]
                        rem_len = e.end - (e.start + ltr_val)
                        removals_by_chrom[e.chrom].append( (e.start + ltr_val, e.end, rem_len) )
                        print(f"Partial excision in group {gid}: {e.chrom}:{e.start}-{e.end} reduced to {e.start}-{e.start+ltr_val} (removed {rem_len} bases)")
                        e.end = e.start + ltr_val
                        new_feature_id = e.feature_id + "_SOLO"
                        e.feature_id = new_feature_id
                        if e.supp:
                            e.name = new_feature_id + ";" + ";".join(e.supp)
                        else:
                            e.name = new_feature_id
                        new_entries.append(e)
                    else:
                        tsd_len = len(e.tsd) if e.tsd != "NA" else 0
                        rem_len = e.length() - tsd_len
                        removals_by_chrom[e.chrom].append( (e.start, e.end - tsd_len, rem_len) )
                        print(f"Full excision in nest group {gid}: {e.chrom}:{e.start}-{e.end} (removed {rem_len} bases)")
                else:
                    new_entries.append(e)
    # Process non-group entries.
    for e in entries:
        if e.group is None:
            if e in to_remove:
                if id(e) in partial_info:
                    ltr_val = partial_info[id(e)]
                    rem_len = e.end - (e.start + ltr_val)
                    removals_by_chrom[e.chrom].append( (e.start + ltr_val, e.end, rem_len) )
                    print(f"Partial excision: {e.chrom}:{e.start}-{e.end} reduced to {e.chrom}:{e.start}-{e.start+ltr_val} (removed {rem_len} bases)")
                    e.end = e.start + ltr_val
                    new_feature_id = e.feature_id + "_SOLO"
                    e.feature_id = new_feature_id
                    if e.supp:
                        e.name = new_feature_id + ";" + ";".join(e.supp)
                    else:
                        e.name = new_feature_id
                    new_entries.append(e)
                else:
                    tsd_len = len(e.tsd) if e.tsd != "NA" else 0
                    rem_len = e.length() - tsd_len
                    removals_by_chrom[e.chrom].append( (e.start, e.end - tsd_len, rem_len) )
                    print(f"Full excision: {e.chrom}:{e.start}-{e.end} (removed {rem_len} bases)")
            else:
                new_entries.append(e)
    
    # -----------------------------
    # Parallel genome sequence update per chromosome:
    # -----------------------------
    # Prepare tasks: for each chromosome, get its removals events and assign a unique seed.
    tasks = []
    for idx, (chrom, rec) in enumerate(sorted(genome_records.items())):
        events = removals_by_chrom.get(chrom, [])
        seed_val = global_seed + idx  # unique and deterministic seed for each chromosome
        tasks.append((chrom, rec, events, seed_val))
    
    updated_genome = {}
    # Use multiprocessing if more than one process is requested.
    if max_proc > 1:
        with Pool(processes=max_proc) as pool:
            results = pool.map(update_chromosome, tasks)
    else:
        results = list(map(update_chromosome, tasks))
    
    for chrom, new_rec, _ in results:
        updated_genome[chrom] = new_rec

    # Adjust bed coordinates based on removed lengths.
    for entry in new_entries:
        shift = 0
        events = sorted(removals_by_chrom.get(entry.chrom, []), key=lambda x: x[0])
        for (rstart, rend, rlen) in events:
            if rstart < entry.start:
                shift += rlen
        entry.start -= shift
        entry.end -= shift

    new_entries.sort(key=lambda x: (x.chrom, x.start))
    return updated_genome, new_entries

# =============================================================================
# Failsafe consolidation for adjacent gene entries
# =============================================================================

def fail_safe_consolidation(bed_entries):
    """
    After simulation, if two or more gene entries with the same feature_ID, strand, and TSD
    are directly adjacent in the sorted BED, consolidate them into a single entry.
    In the consolidated entry, only the feature_ID is kept as NAME (i.e. supplemental attributes are removed).
    """
    consolidated = []
    # Group entries by chromosome.
    entries_by_chrom = defaultdict(list)
    for e in bed_entries:
        entries_by_chrom[e.chrom].append(e)
    
    for chrom, entries in entries_by_chrom.items():
        entries.sort(key=lambda x: x.start)
        i = 0
        while i < len(entries):
            curr = entries[i]
            # Only consolidate genes (feature_ID starting with "gene")
            if curr.feature_id.startswith("gene"):
                j = i + 1
                new_start = curr.start
                new_end = curr.end
                # Consolidate consecutive entries with same feature_ID, strand, and TSD.
                while (j < len(entries) and 
                       entries[j].feature_id.startswith("gene") and
                       entries[j].feature_id == curr.feature_id and
                       entries[j].strand == curr.strand and
                       entries[j].tsd == curr.tsd):
                    new_end = entries[j].end
                    j += 1
                if j > i + 1:
                    # Create a new consolidated entry.
                    new_entry = BedEntry(chrom, new_start, new_end, curr.feature_id, curr.strand, curr.tsd, curr.lineno)
                    new_entry.subtype = "CONSOLIDATED_FAILSAFE"
                    consolidated.append(new_entry)
                    i = j
                    continue
            consolidated.append(curr)
            i += 1
    # Sort the consolidated list before returning.
    consolidated.sort(key=lambda x: (x.chrom, x.start))
    return consolidated

# =============================================================================
# Write outputs
# =============================================================================

def write_fasta(genome_records, output_prefix):
    out_file = output_prefix + ".fasta"
    SeqIO.write(list(genome_records.values()), out_file, "fasta")
    print(f"Updated FASTA written to {out_file}")

def write_bed(bed_entries, output_prefix):
    out_file = output_prefix + ".bed"
    with open(out_file, "w") as f:
        for e in bed_entries:
            f.write(str(e) + "\n")
    print(f"Updated BED written to {out_file}")

# =============================================================================
# Plotting Functions
# =============================================================================

def plot_lognormal(sigma, outname):
    """
    Plot the log-normal gene selection curve.
    The lognormal distribution used has parameters: mu = sigma^2 and sigma.
    """
    mu = sigma ** 2
    # Define x range: avoid 0 and go up to a high quantile
    x_max = math.exp(mu + 3*sigma)
    x = np.linspace(0.001, x_max, 500)
    pdf = (1/(x * sigma * np.sqrt(2*math.pi)) *
           np.exp(- (np.log(x) - mu)**2 / (2*sigma**2)))
    plt.figure(figsize=(6,4))
    plt.plot(x, pdf, lw=2)
    plt.xlabel("Gene Weight")
    plt.ylabel("Probability Density")
    plt.title("Log-normal Distribution (μ = {:.2f}, σ = {:.2f})".format(mu, sigma))
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    print(f"Log-normal distribution figure saved as {outname}")

def plot_weighted_candidate_curve(k, Lmax, outname):
    """
    Plot weighted candidate selection curve.
    The weight is defined as: weight = exp( - k * (1 - (L / Lmax)) ) for L in [0, Lmax].
    """
    L = np.linspace(0, Lmax, 500)
    weight = np.exp(-k * (1 - (L / Lmax)))
    plt.figure(figsize=(6,4))
    plt.plot(L, weight, lw=2)
    plt.xlabel("Feature Length (L)")
    plt.ylabel("Excision Weight")
    plt.title("Weighted Candidate Selection Curve (k = {:.2f})".format(k))
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()
    print(f"Weighted candidate selection curve figure saved as {outname}")

# =============================================================================
# Main function
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Simulate TE excisions from a genome based on a BED file.\n"
                    "Randomly remove TE features (with special handling of nest groups) and "
                    "update the genome FASTA and BED coordinates accordingly.\n\n"
                    "New: For INTACT_LTR elements (TE entries with 'LTRlen' in feature_ID and not '_SOLO'), "
                    "partial excision may be performed. Use --soloLTR_freq to set the percentage (e.g., 10 for 10%).\n"
                    "Additional filtering excludes candidates if the first supplemental attribute contains 'CUT_BY' "
                    "or if a nearby entry (within 100 lines) has an identical TSD, strand, and a NAME that is an exact or prefix match.\n\n"
                    "For disrupted genes, a weight is applied. By default, the number of TE excisions is calculated as:\n"
                    "    rate * generations * (sum of disrupted gene weights).\n"
                    "NEW: Use --fix_ex to bypass this calculation and specify a fixed excision rate (e.g., 1e-6).\n"
                    "      In this case, the number of excisions is calculated as: --fix_ex * genome_size * generations.\n\n"
                    "EDIT: For the non-fixed rate model (--rate), candidate excisions are now further weighted so that\n"
                    "gene-disrupting TEs (NEST_TE_IN_GENE) are favored. Their base weight (from feature length and --k)\n"
                    "is multiplied by (1 + sel_coeff * gene_weight * generations), where gene_weight comes from gene_selection.tsv.\n\n"
                    "CHROMATIN BIAS (rate mode only): Define euchromatin regions as any region within a buffer around a gene.\n"
                    "Candidates whose midpoints lie within these merged intervals get their weight multiplied by --euch_het_bias.\n"
                    "For example, '--euch_het_buffer 10000 --euch_het_bias 1.1' treats 10kb upstream and downstream of each gene\n"
                    "as euchromatin and gives candidate TE excisions in these regions a 10% weight boost.\n\n"
                    "MULTIPROCESSING: Use -m [int] to process up to that many chromosomes in parallel.\n"
                    "Publication-quality figures are generated for the log-normal distribution and the weighted candidate selection curve.\n"
                    "Use --no_fig to disable figure generation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--rate", type=float, default=1e-4, help="Rate of TE deletion per generation per disrupted gene weight.")
    parser.add_argument("--generations", type=int, default=1, help="Number of generations to simulate.")
    parser.add_argument("--bed", required=True, help="Existing BED file with TE/gene coordinates.")
    parser.add_argument("--output", required=True, help="Output prefix (for .bed and .fasta).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--soloLTR_freq", type=float, default=0, help="Percentage of INTACT_LTR excisions to perform partially (e.g., 10 for 10%).")
    parser.add_argument("--sigma", type=float, default=1.0, help="Sigma for log-normal distribution of disrupted_gene_weight (used if gene_selection.tsv is not present).")
    parser.add_argument("--k", type=float, default=1.0, help="Decay rate for weighted excision selection based on sequence length.")
    parser.add_argument("--no_fig", action="store_true", help="Disable generating PDF figures.")
    # New argument: fixed excision rate. When provided, ignore --rate, --generations, gene weights, and selection coefficient.
    parser.add_argument("--fix_ex", type=float, help="If provided, use this fixed excision rate to calculate the number of TE excisions as: fix_ex * genome_size * generations.")
    # New argument: selection coefficient for gene-disrupting TEs.
    parser.add_argument("--sel_coeff", type=float, default=0.0, help="Selection coefficient to weight excision of gene-disrupting TEs.")
    # New arguments for chromatin bias in rate mode.
    parser.add_argument("--euch_het_buffer", type=int, default=0,
                        help="Buffer (in bp) around gene features to be considered euchromatin. Only interpreted in rate mode.")
    parser.add_argument("--euch_het_bias", type=float, default=1.0,
                        help="Bias factor to increase probability of TE excision in euchromatin. Only interpreted in rate mode.")
    # New argument for multiprocessing: max number of chromosomes processed in parallel.
    parser.add_argument("-m", type=int, default=1, help="Maximum number of chromosomes to process in parallel.")
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    print("Parsing genome FASTA...")
    genome_records = parse_fasta(args.genome)
    # Calculate genome_size: total number of bases in the genome.
    genome_size = sum(len(rec.seq) for rec in genome_records.values())
    print(f"Genome size (total bases): {genome_size}")

    print("Parsing BED file...")
    entries = parse_bed(args.bed)

    # Extract unique geneIDs (feature_IDs beginning with 'gene') from the BED entries.
    unique_genes = {e.feature_id for e in entries if e.feature_id.startswith("gene")}
    gene_sel_file = "./gene_selection.tsv"
    gene_weights = {}
    if os.path.exists(gene_sel_file):
        print(f"Loading gene selection weights from {gene_sel_file} ...")
        with open(gene_sel_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    gene_weights[parts[0]] = float(parts[1])
    else:
        # Only generate synthetic gene weights if fixed excision rate is not used
        if args.fix_ex is None:
            print(f"{gene_sel_file} not found. Generating synthetic gene selection weights ...")
            mu = args.sigma ** 2
            # Use sorted list of unique genes to ensure deterministic ordering.
            for gene in sorted(unique_genes):
                gene_weights[gene] = random.lognormvariate(mu, args.sigma)
            with open(gene_sel_file, "w") as f:
                for gene, weight in gene_weights.items():
                    f.write(f"{gene}\t{weight:.4f}\n")
            print(f"Synthetic gene selection weights written to {gene_sel_file}")
        else:
            print("Fixed excision rate provided (--fix_ex), so gene weights generation is bypassed.")

    # Optionally generate log-normal distribution figure.
    if not args.no_fig and args.fix_ex is None:
        plot_lognormal(args.sigma, "lognormal_distribution.pdf")

    print("Classifying BED entries...")
    entries, nest_groups = classify_entries(entries)

    # For plotting the weighted candidate selection curve, determine Lmax from eligible entries.
    eligible_entries = [e for e in entries if e.subtype in ["NON_NEST_GROUP_TE", "CUT_PAIR_TE_1", "CUT_PAIR_TE_2",
                                                           "NEST_TE_IN_TE", "NEST_TE_IN_GENE"]]
    if eligible_entries:
        Lmax = max(e.length() for e in eligible_entries)
        if not args.no_fig:
            plot_weighted_candidate_curve(args.k, Lmax, "weighted_candidate_selection.pdf")
    else:
        print("No eligible entries for weighted candidate selection curve plot.")

    # Determine the number of excisions.
    if args.fix_ex is not None:
        lambda_val = args.fix_ex * genome_size * args.generations
        excision_count = np.random.poisson(lambda_val)
        print(f"Using fixed excision rate: {args.fix_ex}")
        print(f"Lambda for Poisson (fix_ex mode): {lambda_val:.4f}")
        print(f"Calculated number of TE excisions: {excision_count}")
    else:
        excision_count = calculate_excision_count(entries, args.rate, args.generations, gene_weights)
    
    # Build euchromatin intervals if in rate mode and a buffer is specified.
    euch_intervals = None
    if args.fix_ex is None and args.euch_het_buffer > 0:
        euch_intervals = build_euchromatin_intervals(entries, args.euch_het_buffer)
    
    removals = parallel_select_removals(
        entries=entries,
        nest_groups=nest_groups,
        rate=args.rate,
        generations=args.generations,
        gene_weights=gene_weights,
        sel_coeff=args.sel_coeff,
        euch_intervals=euch_intervals,
        euch_het_bias=args.euch_het_bias,
        seed=args.seed,
        k=args.k,
        max_proc=args.m,
        fix_ex=args.fix_ex,
        genome_size=genome_size
    )
    # Print removal events in a sorted order for determinism.
    print("Selected removal events:")
    for e in sorted(removals, key=lambda x: (x.chrom, x.start, x.end, x.lineno)):
        print(f" - {e.chrom}:{e.start}-{e.end}, {e.feature_id}, subtype: {e.subtype}, length: {e.length()}")

    print("Simulating excisions and updating coordinates...")
    updated_genome, updated_bed = simulate_excision(genome_records, entries, nest_groups, removals, args.soloLTR_freq, args.m, args.seed)

    # Failsafe: Consolidate adjacent gene entries that share the same feature_ID.
    updated_bed = fail_safe_consolidation(updated_bed)

    write_fasta(updated_genome, args.output)
    write_bed(updated_bed, args.output)

if __name__ == "__main__":
    main()
