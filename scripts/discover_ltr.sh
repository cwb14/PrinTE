#!/usr/bin/env bash
# Annotate the LTR-RTs in every simulation directory of a grid search, so the results
# can be scored. Run this from the directory holding the per-parameter simulation
# directories, between the sweep and the scoring step.
#
# This is a template. The annotator is synLTR, which is not part of PrinTE, so point
# SYNLTR at your clone. Everything else has a sensible default you can override:
#
#   SYNLTR=~/synLTR ./scripts/discover_ltr.sh
#   SYNLTR=~/synLTR GENERATION=5100000 THREADS=64 ./scripts/discover_ltr.sh
#
# See docs/ltr-annotation.md for what synLTR produces and docs/nextflow.md for where
# this fits in a sweep.

set -euo pipefail

SYNLTR="${SYNLTR:?Set SYNLTR to your synLTR checkout, e.g. SYNLTR=~/synLTR $0}"
GENERATION="${GENERATION:-}"
THREADS="${THREADS:-16}"
PROTEINS="${PROTEINS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/TAIR10.pep.fa.gz}"

harvest="${SYNLTR}/module2/ltrharvest.py"
[[ -f "$harvest" ]] || { echo "Not found: $harvest" >&2; exit 1; }
[[ -f "$PROTEINS" ]] || { echo "Not found: $PROTEINS" >&2; exit 1; }

for d in */; do
  [[ -d "$d" ]] || continue

  # Default to the last generation this directory actually produced, so the same
  # command works across runs with different --generation_end.
  if [[ -n "$GENERATION" ]]; then
    fasta="gen${GENERATION}_final.fasta"
  else
    fasta=$(ls "$d"gen*_final.fasta 2>/dev/null | sort -V | tail -1 | xargs -r basename)
  fi
  if [[ -z "$fasta" || ! -f "$d$fasta" ]]; then
    echo "  $d: no gen*_final.fasta, skipping."
    continue
  fi

  prefix="${fasta%_final.fasta}"
  if [[ -s "$d${prefix}_r1_ltr.tsv" ]]; then
    echo "  $d: already annotated, skipping."
    continue
  fi

  echo ">>> $d ($fasta)"
  (
    cd "$d"
    python "$harvest" \
      --genome "$fasta" \
      --proteins "$PROTEINS" \
      --threads "$THREADS" \
      --scn-min-ltr-len 100 \
      --scn-min-ret-len 800 \
      --scn-max-ret-len 15000 \
      --scn-min-int-len 500 \
      --scn-max-int-len 12000 \
      --out-prefix "$prefix"
  )
done
