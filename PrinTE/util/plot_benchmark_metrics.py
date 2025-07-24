#!/usr/bin/env python3
"""
Script to run bedtools.py across three replicates in parallel and plot
detection_rate, precision, FDR, and sensitivity.
Generates a grouped bar chart saved as metrics.pdf in the current directory.

I can run:
bash ../TESS/TESS.sh PrinTE --burnin_only --cds_percent 0 --TE_num 800 --chr_number 1 --size 1Gb -s 21 --TE_lib ltr-db.fa --threads 20
bash ../TESS/TESS.sh PrinTE --burnin_only --cds_percent 0 --TE_num 800 --chr_number 1 --size 1Gb -s 32 --TE_lib ltr-db.fa --threads 20
bash ../TESS/TESS.sh PrinTE --burnin_only --cds_percent 0 --TE_num 800 --chr_number 1 --size 1Gb -s 44 --TE_lib ltr-db.fa --threads 20

Then LTR detection using any method...

Then feed SCN or pass listed formmatted data into this script. 
"""

import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed

# Define datasets and SCN file templates
datasets = [
    ('harvest',           'rep{rep}/burnin/burnin.fa.harvest.combine.scn'),
    ('finder',            'rep{rep}/burnin/burnin.fa.finder.combine.scn'),
    ('merged',            'rep{rep}/burnin/burnin_merged.scn'),
    ('pass_trunc_clean',  'rep{rep}/burnin/burnin.fa.defalse.clean'),
]

bed_template   = 'rep{rep}/burnin.bed'
r_threshold    = '0.96'

# Prepare data structures, now including 'sensitivity'
metrics = {
    name: {
        'detection_rate': [],
        'precision': [],
        'FDR': [],
        'sensitivity': []
    }
    for name, _ in datasets
}

# Regex patterns to parse bedtools.py output
pattern_overlap    = re.compile(r'Overlapping entries: \d+ \((\d+) unique\)')
pattern_unique_scn = re.compile(r'Entries unique to SCN/PASS file: (\d+)')
pattern_unique_bed = re.compile(r'Entries unique to BED file: (\d+)')

def run_analysis(rep, name, scn_template):
    """Run bedtools.py for one replicate & dataset, parse output, compute metrics."""
    scn_path = scn_template.format(rep=rep)
    bed_path = bed_template.format(rep=rep)
    cmd = [
        'python', 'bedtools.py',
        '-pass_scn', scn_path,
        '-bed',      bed_path,
        '-r',        r_threshold
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd,
                            capture_output=True,
                            text=True,
                            check=True)
    out = result.stdout

    # Extract counts
    tp = int(pattern_overlap.search(out).group(1))       # true positives
    fp = int(pattern_unique_scn.search(out).group(1))    # false positives
    fn = int(pattern_unique_bed.search(out).group(1))    # false negatives

    # Compute metrics
    detection_rate = tp / (tp + fn)         # same as sensitivity
    precision      = tp / (tp + fp)
    fdr            = fp / (tp + fp)
    sensitivity    = tp / (tp + fn)         # explicitly labeled

    return name, detection_rate, precision, fdr, sensitivity

# Run analyses in parallel
futures = []
with ThreadPoolExecutor() as executor:
    for rep in [1, 2, 3]:
        for name, template in datasets:
            futures.append(
                executor.submit(run_analysis, rep, name, template)
            )

    for future in as_completed(futures):
        name, det, prec, fdr, sens = future.result()
        metrics[name]['detection_rate'].append(det)
        metrics[name]['precision'].append(prec)
        metrics[name]['FDR'].append(fdr)
        metrics[name]['sensitivity'].append(sens)

# Compute means and standard deviations
names = [n for n, _ in datasets]
metric_keys = ['detection_rate', 'precision', 'FDR', 'sensitivity']

means = {
    m: [np.mean(metrics[n][m]) for n in names]
    for m in metric_keys
}
stds = {
    m: [np.std(metrics[n][m], ddof=1) for n in names]
    for m in metric_keys
}

# Plot grouped bar chart
x = np.arange(len(names))
width = 0.2
# Center the 4 bars around each x: offsets = [-1.5, -0.5, +0.5, +1.5] * width
offsets = [(i - (len(metric_keys)-1)/2) * width for i in range(len(metric_keys))]

fig, ax = plt.subplots(figsize=(8, 4))
for i, metric in enumerate(metric_keys):
    ax.bar(
        x + offsets[i],
        means[metric],
        width,
        yerr=stds[metric],
        label=metric.replace('_', ' ').title()
    )

ax.set_xticks(x)
ax.set_xticklabels(names, rotation=45, ha='right')
ax.set_ylabel('Value')
ax.set_title('Metrics across datasets')
ax.legend()
plt.tight_layout()
plt.savefig('metrics.pdf')
print("Saved plot to metrics.pdf")
