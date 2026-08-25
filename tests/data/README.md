# Test fixtures

Small synthetic files, about 120 KB in total, so the suite runs in CI without carrying
real genomes. Everything here is derived from the bundled reference libraries.

| file | what it is |
|---|---|
| `tiny_TE.lib` | 21 TE sequences, 35 kb. 13 LTR-RTs carrying `~LTRlen:` tags plus 8 non-LTR elements, drawn from a mix of superfamilies so no `ratios.tsv` row goes unexercised. |
| `tiny_cds.fa` | 100 CDS, 300-1200 bp each, from TAIR10. |
| `tiny_ratios.tsv` | Weights covering exactly the superfamilies present in `tiny_TE.lib`, with eight families flagged `cutpaste` so the Type 2 path is exercised. |
| `sample.gff` | Five hand-written TEgenomeSimulator-style records for the `gff_to_bed` converter, including two partial elements and both `TSD_5`/`TSD_3` cases. |

## Rebuilding them

`tiny_TE.lib` and `tiny_cds.fa` are subsets of the shipped libraries. To regenerate:

```bash
# Build a cleaned library once, which is what PrinTE writes as lib_clean.fa
printe --burnin_only -sz 200kb -cn 1 -P 4 -itp 10 -t 4 -s 1

python - <<'PY'
import random, re
random.seed(11)

def read_fa(path):
    recs, hdr, seq = [], None, []
    for line in open(path):
        if line.startswith('>'):
            if hdr: recs.append((hdr, ''.join(seq)))
            hdr, seq = line.strip(), []
        else: seq.append(line.strip())
    if hdr: recs.append((hdr, ''.join(seq)))
    return recs

lib = read_fa('lib_clean.fa')
ltr = [r for r in lib if 'LTRlen' in r[0] and len(r[1]) <= 4000]
non = [r for r in lib if 'LTRlen' not in r[0] and len(r[1]) <= 2500]
pick = random.sample(ltr, 13) + random.sample(non, 8)
random.shuffle(pick)
with open('tiny_TE.lib', 'w') as o:
    for h, s in pick:
        o.write(h + '\n' + s + '\n')
PY
```

`tiny_ratios.tsv` must then be regenerated to match, or the fixture library will contain
superfamilies with no declared weight. `tests/unit/test_ratios.py` checks exactly that.

`sample.gff` is hand-written; edit it directly.
