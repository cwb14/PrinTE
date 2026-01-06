#!/usr/bin/env bash
set -euo pipefail

fasta="gen5100000_final.fasta"
out="gen5100000_ltrharvest_kmer2ltr_dedup"

for d in */; do
  [[ -d "$d" ]] || continue
  echo ">>> Entering directory: $d"

  pushd "$d" > /dev/null

  # 1) Skip if the fasta file doesn't exist
  if [[ ! -f "$fasta" ]]; then
    echo "    - $fasta not found, skipping."
    popd > /dev/null
    continue
  fi

  # 2) Skip if the output file exists and is non-empty
  if [[ -s "$out" ]]; then
    echo "    - $out exists and is non-empty, skipping."
    popd > /dev/null
    continue
  fi

  # 3) Otherwise run
  echo "    - Running ltrharvest.py in $(pwd)"
  python ../ltrharvest.py \
    --genome "$fasta" \
    --proteins ../PrinTE/data/TAIR10.pep.fa.gz \
    --threads 200 \
    --scn-min-ltr-len 100 \
    --scn-min-ret-len 800 \
    --scn-max-ret-len 15000 \
    --scn-min-int-len 500 \
    --scn-max-int-len 12000 \
    --out-prefix gen5100000_ltrharvest 

  popd > /dev/null
done
