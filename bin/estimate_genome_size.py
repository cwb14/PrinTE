#!/usr/bin/env python3
"""Back-of-napkin terminal genome-size / event-count estimator for PrinTE.

OBSERVE-ONLY. Predicts, before a fixed-rate run, the per-step genome-size
trajectory (bp and FASTA bytes) and cumulative insertion/deletion event
counts, and (in --check mode) compares the prediction to the actual
per-step trajectory. Makes NO decisions. See the design spec.
"""
import argparse
import gzip
import io
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

BYTES_PER_BP = 61.0 / 60.0  # Biopython SeqIO "fasta" wraps at 60 chars (+1 newline)

_LTRLEN_RE = re.compile(r"LTRlen:(\d+)")

_EST_COLS = ["ins", "del", "sr", "k", "ge", "step", "iterations",
             "g0_bp", "e_lins", "e_lrem", "r", "bytes_factor",
             "iter", "pred_size_bp", "pred_size_bytes", "pred_cum_nins", "pred_cum_ndel"]


def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")


def parse_ratios(path) -> Dict[Tuple[str, str], float]:
    """Read ratios.tsv -> {(class, superfamily): intact_frequency}. Skips
    blank lines and lines starting with '#'."""
    out = {}
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            cls, sf, freq = parts[0].strip(), parts[1].strip(), parts[2].strip()
            try:
                out[(cls, sf)] = float(freq)
            except ValueError:
                continue
    return out


def library_mean_lengths(path) -> Dict[Tuple[str, str], float]:
    """Parse a TE library FASTA (headers '>name#CLASS/SUPERFAMILY') and return
    {(class, superfamily): mean sequence length in bp}."""
    sums: Dict[Tuple[str, str], int] = {}
    counts: Dict[Tuple[str, str], int] = {}
    cur_key = None
    cur_len = 0

    def _flush():
        nonlocal cur_key, cur_len
        if cur_key is not None:
            sums[cur_key] = sums.get(cur_key, 0) + cur_len
            counts[cur_key] = counts.get(cur_key, 0) + 1
        cur_key, cur_len = None, 0

    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                header = line[1:].strip()
                cur_key = _key_from_classsf(header.split("#", 1)[1]) if "#" in header else None
                cur_len = 0
            else:
                cur_len += len(line.strip())
    _flush()
    return {k: sums[k] / counts[k] for k in sums}


def _key_from_classsf(classsf: str) -> Tuple[str, str]:
    """'LTR/Copia' -> ('LTR', 'Copia'); 'LTR/Copia~LTRlen:200' -> ('LTR','Copia')."""
    classsf = classsf.split("~", 1)[0].strip()
    if "/" in classsf:
        cls, sf = classsf.split("/", 1)
        return (cls.strip(), sf.strip())
    return (classsf.strip(), "unknown")


def expected_tsd_length(te_class: str, te_superfamily: str) -> float:
    """Expected TSD length (bp) added per insertion. Mirrors get_tsd_length in
    nest_inserter_parallel.py with the *means* of its random ranges."""
    c = (te_class or "").strip().lower()
    s = (te_superfamily or "").strip().lower()
    const5, l1, cacta, hat, mule, tourist = 5.0, 12.5, 3.0, 6.5, 8.5, 6.0
    table = {
        ("ltr", "copia"): const5, ("ltr", "gypsy"): const5, ("ltr", "solo"): const5,
        ("ltr", "ty3"): const5, ("ltr", "unknown"): const5, ("ltr", "crm"): const5,
        ("ltr", "trim"): const5,
        ("dna", "harbinger"): 3.0, ("dna", "mariner"): 2.0,
        ("dnaauto", "helitron"): 0.0, ("dna", "helitron"): 0.0, ("dnanona", "helitron"): 0.0,
        ("mite", "stow"): tourist, ("mite", "tourist"): tourist,
        ("sine", "trna"): l1, ("sine", "unknown"): l1,
        ("line", "l1"): l1, ("line", "unknown"): l1, ("line", "rte"): l1,
        ("dnaauto", "cacta"): cacta, ("dna", "cacta"): cacta, ("dnanona", "cacta"): cacta,
        ("dna", "hat"): hat, ("dnanona", "hat"): hat, ("dnaauto", "hat"): hat,
        ("dnaauto", "mudr"): mule, ("dnanona", "mule"): mule, ("dna", "mudr"): mule,
        ("dnaauto", "mule"): mule,
        ("dna", "dta"): hat, ("mite", "dta"): hat,
        ("dna", "dtc"): cacta, ("mite", "dtc"): cacta,
        ("dna", "dth"): tourist, ("mite", "dth"): tourist,
        ("dna", "dtm"): mule, ("mite", "dtm"): mule,
        ("dna", "dtt"): tourist, ("mite", "dtt"): tourist,
        ("dnaauto", "cactg"): cacta, ("dnaauto", "mle"): tourist,
        ("dnaauto", "pile"): tourist, ("dnaauto", "pole"): tourist,
        ("dnanona", "cactg"): cacta, ("dnanona", "mle"): tourist, ("dnanona", "muletir"): mule,
        ("dnanona", "pile"): tourist, ("dnanona", "pole"): tourist,
        ("dnanona", "tourist"): tourist, ("dnanona", "unknown"): tourist,
        ("dna", "mutator"): mule, ("dna", "tc1_mariner"): tourist, ("dna", "unknown"): tourist,
    }
    return table.get((c, s), 5.0)


def expected_insertion_length(ratios: Dict[Tuple[str, str], float],
                              lib_means: Dict[Tuple[str, str], float]) -> float:
    """E[L_ins] = sum over categories present in BOTH ratios and library of
    w_cat * (mean_lib_len + expected_tsd), with weights renormalized over the
    present categories."""
    shared = [(c, sf) for (c, sf) in ratios if (c, sf) in lib_means]
    total_w = sum(ratios[k] for k in shared)
    if total_w <= 0:
        return 0.0
    out = 0.0
    for (c, sf) in shared:
        w = ratios[(c, sf)] / total_w
        out += w * (lib_means[(c, sf)] + expected_tsd_length(c, sf))
    return out


def parse_bed_elements(path) -> List[dict]:
    """Read a PrinTE BED. Keep TE features (feature_id contains '#'); drop genes.
    Returns dicts with: length, ltrlen, tsd_len, solo_qualifies.

    solo_qualifies mirrors TE_exciser_parallel.py:526 — feature_id contains
    'LTRlen' and not '_SOLO' and not '_FRAG'.
    """
    out = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 6:
                continue
            start, end, fid, tsd = int(f[1]), int(f[2]), f[3], f[4]
            if "#" not in fid:  # gene or non-TE
                continue
            primary = fid.split(";", 1)[0]  # element annotation before nesting markers
            m = _LTRLEN_RE.search(primary)
            ltrlen = int(m.group(1)) if m else 0
            tsd_len = 0 if tsd.strip().upper() == "NA" else len(tsd.strip())
            solo_q = ("LTRlen" in fid) and ("_SOLO" not in fid) and ("_FRAG" not in fid)
            out.append({
                "length": end - start,
                "ltrlen": ltrlen,
                "tsd_len": tsd_len,
                "solo_qualifies": solo_q,
            })
    return out


def expected_removal_length(elements: List[dict], k: float, solo_rate: float) -> float:
    """E[L_rem]: mean bp removed per excision, with length-bias k and solo
    conversion. Mirrors TE_exciser_parallel.py weighting (:393-394) and the
    solo/TSD bp accounting (:616-672).

    weight_i = exp(-k * (1 - min(L_i, L95)/L95))   # longer favored as k grows
    L_sel    = sum(w_i * L_i) / sum(w_i)            # full element length removed
    p_solo   = (solo_rate/100) * (weighted fraction that solo-qualify)
    saving   = p_solo * E[LTRlen | qualifying] + (1 - p_solo) * E[TSD]
    E[L_rem] = L_sel - saving
    """
    if not elements:
        return 0.0
    lengths = np.array([e["length"] for e in elements], dtype=float)
    L95 = float(np.percentile(lengths, 95))
    if L95 <= 0:
        L95 = float(lengths.max()) or 1.0
    L_eff = np.minimum(lengths, L95)
    w = np.exp(-k * (1.0 - L_eff / L95))
    W = float(w.sum())
    if W <= 0:
        return 0.0
    L_sel = float((w * lengths).sum() / W)

    qual = np.array([1.0 if e["solo_qualifies"] else 0.0 for e in elements])
    ltrlen = np.array([e["ltrlen"] for e in elements], dtype=float)
    tsd = np.array([e["tsd_len"] for e in elements], dtype=float)

    Wq = float((w * qual).sum())
    f_intactLTR = Wq / W
    p_solo = (solo_rate / 100.0) * f_intactLTR
    E_ltrlen = float((w * qual * ltrlen).sum() / Wq) if Wq > 0 else 0.0
    E_tsd = float((w * tsd).sum() / W)
    saving = p_solo * E_ltrlen + (1.0 - p_solo) * E_tsd
    return L_sel - saving


def genome_size_bp(fasta_path) -> int:
    """Total sequence length (bp). Prefers a FASTA `.fai` index (instant); if
    absent, creates one with pysam so repeat runs on the same genome are fast,
    then reads it. Falls back to a streaming scan if indexing is unavailable
    (e.g. pysam missing, or a gzipped FASTA that samtools cannot faidx)."""
    fasta_path = str(fasta_path)
    fai = fasta_path + ".fai"
    if not Path(fai).exists() and not fasta_path.endswith(".gz"):
        try:
            import pysam
            pysam.faidx(fasta_path)  # writes <fasta>.fai
        except Exception:
            pass  # fall through to streaming
    if Path(fai).exists():
        total = 0
        with open(fai) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) >= 2:
                    total += int(parts[1])
        return total
    return _stream_genome_size_bp(fasta_path)


def _stream_genome_size_bp(fasta_path) -> int:
    """Fallback: sum sequence lengths by streaming the FASTA (low memory)."""
    total = 0
    with _open_text(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def bp_to_bytes(bp: float) -> int:
    """Approximate FASTA file bytes for `bp` of sequence. Biopython 'fasta'
    wraps at 60 chars (1 newline / 60 bp); header bytes are negligible at the
    genome scale and are calibrated empirically by --check."""
    return int(round(bp * BYTES_PER_BP))


def project_trajectory(g0_bp: float, ins: float, dele: float, step: float,
                       iterations: int, e_lins: float, e_lrem: float):
    """Iterate G_i = G_{i-1}*(1+r), r = step*(ins*E[L_ins] - del*E[L_rem]).
    Returns (rows, r) where each row has iter, size_bp, cum_nins, cum_ndel.
    Event counts use the size at the START of each step (as the simulator does).
    Size is floored at 1 bp to avoid runaway-negative projections."""
    r = step * (ins * e_lins - dele * e_lrem)
    g = float(g0_bp)
    cum_nins = 0.0
    cum_ndel = 0.0
    rows = [{"iter": 0, "size_bp": g, "cum_nins": 0.0, "cum_ndel": 0.0}]
    for i in range(1, iterations + 1):
        cum_nins += ins * g * step
        cum_ndel += dele * g * step
        g = max(1.0, g * (1.0 + r))
        rows.append({"iter": i, "size_bp": g, "cum_nins": cum_nins, "cum_ndel": cum_ndel})
    return rows, r


def write_estimate(path, params, g0_bp, e_lins, e_lrem, r, rows):
    """Long format: one row per iteration, scalars repeated (easy to concat)."""
    iterations = max(row["iter"] for row in rows)
    with open(path, "w") as fh:
        fh.write("\t".join(_EST_COLS) + "\n")
        for row in rows:
            rec = [
                params["ins"], params["del"], params["sr"], params["k"],
                params["ge"], params["step"], iterations,
                f"{g0_bp:.0f}", f"{e_lins:.4f}", f"{e_lrem:.4f}", f"{r:.8g}",
                f"{BYTES_PER_BP:.6f}",
                row["iter"], f"{row['size_bp']:.0f}", bp_to_bytes(row["size_bp"]),
                f"{row['cum_nins']:.2f}", f"{row['cum_ndel']:.2f}",
            ]
            fh.write("\t".join(str(x) for x in rec) + "\n")


def read_actual_trajectory(path) -> Dict[int, int]:
    """PrinTE.sh genome_size_trajectory.tsv -> {iteration: bytes}."""
    out = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                out[int(parts[0])] = int(parts[1])
    return out


def run_check(estimate_path, trajectory_path, out_path) -> List[dict]:
    """Join predicted (per iter, bytes) with actual per-step bytes; write
    napkin_check.tsv with pct_err. Returns the rows (also for printing)."""
    pred = {}
    scalars = {}
    with open(estimate_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            pred[int(d["iter"])] = int(d["pred_size_bytes"])
            scalars = d  # last row's scalars are fine (repeated)
    actual = read_actual_trajectory(trajectory_path) if Path(trajectory_path).exists() else {}

    cols = ["ins", "del", "sr", "k", "step", "iter",
            "pred_size_bytes", "actual_size_bytes", "pct_err"]
    out_rows = []
    for it in sorted(set(pred) | set(actual)):
        p = pred.get(it)
        a = actual.get(it)
        pct = "" if (p is None or a in (None, 0)) else f"{(p - a) / a * 100:.2f}"
        out_rows.append({
            "ins": scalars.get("ins", ""), "del": scalars.get("del", ""),
            "sr": scalars.get("sr", ""), "k": scalars.get("k", ""),
            "step": scalars.get("step", ""), "iter": it,
            "pred_size_bytes": "" if p is None else p,
            "actual_size_bytes": "" if a is None else a, "pct_err": pct,
        })
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return out_rows


def main():
    p = argparse.ArgumentParser(description="PrinTE back-of-napkin genome-size estimator (observe-only)")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--check", action="store_true", help="Comparison mode: predicted vs actual trajectory")
    # estimate-mode inputs
    p.add_argument("--fasta"); p.add_argument("--bed"); p.add_argument("--lib"); p.add_argument("--ratios")
    p.add_argument("--ins", type=float); p.add_argument("--del", dest="dele", type=float)
    p.add_argument("--step", type=int); p.add_argument("--ge", type=int)
    p.add_argument("--solo-rate", type=float, default=95.0); p.add_argument("--k", type=float, default=10.0)
    p.add_argument("--out", default="napkin_estimate.tsv",
                   help="estimate mode: where to WRITE the prediction")
    # check-mode inputs
    p.add_argument("--estimate-in", default="napkin_estimate.tsv",
                   help="check mode: the prediction file to READ")
    p.add_argument("--trajectory", default="genome_size_trajectory.tsv")
    p.add_argument("--check-out", default="napkin_check.tsv")
    args = p.parse_args()

    if args.check:
        rows = run_check(args.estimate_in, args.trajectory, args.check_out)
        print("Napkin check (predicted vs actual genome bytes):")
        for r in rows:
            print(f"  iter {r['iter']}: pred={r['pred_size_bytes']} "
                  f"actual={r['actual_size_bytes']} err={r['pct_err']}%")
        print(f"Wrote {args.check_out}")
        return 0

    for req in ("fasta", "bed", "lib", "ratios", "ins", "dele", "step", "ge"):
        if getattr(args, req) is None:
            p.error(f"--{'del' if req == 'dele' else req} is required in estimate mode")

    g0 = genome_size_bp(args.fasta)
    ratios = parse_ratios(args.ratios)
    libm = library_mean_lengths(args.lib)
    e_lins = expected_insertion_length(ratios, libm)
    e_lrem = expected_removal_length(parse_bed_elements(args.bed), args.k, args.solo_rate)
    iterations = args.ge // args.step
    rows, r = project_trajectory(g0, args.ins, args.dele, args.step, iterations, e_lins, e_lrem)
    write_estimate(args.out, {"ins": args.ins, "del": args.dele, "sr": args.solo_rate,
                              "k": args.k, "ge": args.ge, "step": args.step},
                   g0, e_lins, e_lrem, r, rows)

    term = rows[-1]
    print(f"Napkin estimate (fixed-rate): E[L_ins]={e_lins:.1f} bp, E[L_rem]={e_lrem:.1f} bp, "
          f"per-step r={r:.4g}")
    print(f"  G0={g0:,} bp  ->  predicted terminal={term['size_bp']:,.0f} bp "
          f"({bp_to_bytes(term['size_bp']):,} bytes) over {iterations} steps")
    print(f"  predicted cumulative events: insertions={term['cum_nins']:,.0f}, "
          f"deletions={term['cum_ndel']:,.0f}")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
