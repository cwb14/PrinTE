#!/usr/bin/env python3
"""Shared helpers for conserved cut-and-paste TE modeling.

Dependency-free (stdlib only) so it is trivially unit-testable and importable
by both nest_inserter_parallel.py and TE_exciser_parallel.py (siblings in
bin/, which is on sys.path[0] when those scripts run by path).

ratios.tsv 5th column ('transposition'): 'cutpaste' (case-insensitive) marks a
(class, superfamily) family as conserved cut-and-paste. Anything else (missing
column / blank / unrecognized) is copy-and-paste (legacy behavior).
"""
import os
import re

_FAMILY_RE = re.compile(r"[^#]+#([^/]+)/([^~;]+)")


def load_cutpaste_set(ratio_path):
    """Return {(class, superfamily), ...} flagged 'cutpaste' in the file.

    Tolerates comment (#) / blank lines and rows with <5 fields. A None or
    missing path yields an empty set (=> all copy-and-paste).
    """
    result = set()
    if not ratio_path or not os.path.exists(ratio_path):
        return result
    with open(ratio_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            if parts[4].lower() == "cutpaste":
                result.add((parts[0], parts[1]))
    return result


def parse_family(feature_id):
    """Extract (class, superfamily) from a feature_id like 'name#CLASS/SUPER'.

    Suffixes (~LTRlen:N, _FRAG, _SOLO) live outside the class/superfamily
    capture, so they do not affect the result. Returns ('unknown','unknown')
    when the pattern does not match (e.g. gene features).
    """
    m = _FAMILY_RE.match(feature_id)
    if m:
        # groups may carry leading/trailing space only for malformed ids;
        # the strip() is intentional defensive cleanup -- do not remove.
        return (m.group(1).strip(), m.group(2).strip())
    return ("unknown", "unknown")


def write_debt(path, tally):
    """Write {(class,superfamily): count} as 'class\\tsuperfamily\\tcount'.

    Always (re)creates the file, truncating. Only counts > 0 are written;
    an all-zero/empty tally yields a valid empty file. Raises OSError if
    *path* is not writable (fail-fast; callers handle).
    """
    with open(path, "w") as fh:
        for (cls, sf), cnt in sorted(tally.items()):
            if cnt > 0:
                fh.write(f"{cls}\t{sf}\t{int(cnt)}\n")


def read_debt(path):
    """Parse a debt file into {(class,superfamily): count}.

    Missing file => {}. Malformed lines (wrong arity / non-int count) are
    skipped silently.
    """
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                continue
            try:
                cnt = int(parts[2])
            except ValueError:
                continue
            if cnt > 0:
                result[(parts[0], parts[1])] = cnt
    return result
