"""Pure planning + rendering for PrinTE sbatch generation.

Takes a ClusterProfile (from slurm_introspect) plus an SbatchRequest of user
inputs and produces an SbatchPlan (sizing, packing, chunking, throttle), then
renders the sbatch script, packing runner, and optional chunk submitter. No
subprocess calls here -- fully unit-testable without a cluster.
"""

import math
from dataclasses import dataclass, field
from pathlib import Path

from .slurm_introspect import ClusterProfile, _time_to_minutes


@dataclass
class SbatchRequest:
    threads: int
    n_combos: int
    job_name: str
    outdir: str
    partition: str | None
    account: str | None
    qos: str | None
    cpus_per_task: int | None
    mem: str | None
    time: str
    pack: str                       # "auto" | "on" | "off"
    combos_per_node: int | None
    max_concurrent: int | None
    use_logical_cpus: bool
    no_probe: bool
    probe_timeout: int


@dataclass
class SbatchPlan:
    mode: str                       # "single" | "packed"
    combos_per_element: int
    n_combos: int
    n_elements: int
    array_throttle: int | None
    cpus_per_task: int
    exclusive: bool
    mem_directive: str
    time_str: str
    partition: str | None
    account: str | None
    qos: str | None
    chunks: list[tuple[int, int]]
    directives: list[str]
    header_comments: list[str]
    degraded: bool
    warnings_out: list[str] = field(default_factory=list)


def _minutes_to_hms(total_min: int) -> str:
    h, m = divmod(max(1, total_min), 60)
    return f"{h:02d}:{m:02d}:00"


def plan_sbatch(profile: ClusterProfile, req: SbatchRequest) -> SbatchPlan:
    if req.n_combos <= 0:
        raise ValueError("plan_sbatch requires n_combos >= 1")
    warnings = list(profile.warnings)
    degraded = not profile.slurm_available and not req.no_probe

    # --- cpus-per-task ---
    cpus = req.cpus_per_task or req.threads
    cap = profile.node_cpus or cpus
    if profile.max_cpus_per_node:
        cap = min(cap, profile.max_cpus_per_node)
    cpus = max(1, min(cpus, cap))
    if cpus < req.threads:
        warnings.append(
            f"cpus-per-task clamped from {req.threads} to {cpus} (node CPU limit)")

    # --- packing decision ---
    if req.pack == "on":
        do_pack = True
    elif req.pack == "off":
        do_pack = False
    else:
        do_pack = (profile.exclusivity == "exclusive")

    if do_pack:
        usable = profile.max_cpus_per_node or profile.node_cpus or cpus
        if (profile.threads_per_core and profile.threads_per_core > 1
                and not req.use_logical_cpus and profile.node_cores):
            usable = profile.node_cores
        K = req.combos_per_node or max(1, usable // cpus)
        K = max(1, min(K, req.n_combos))
    else:
        K = 1
    if do_pack and K == 1:
        warnings.append("Packing requested but only 1 combo fits per node; single mode.")
    n_elements = math.ceil(req.n_combos / K)
    mode = "packed" if (do_pack and K > 1) else "single"
    exclusive = (mode == "packed")

    # --- memory ---
    if req.mem:
        mem_directive = f"--mem={req.mem}"
    elif mode == "packed":
        mem_directive = "--mem=0"                       # all node RAM, paired w/ --exclusive
    elif profile.node_mem_mib and profile.node_cpus:
        share = math.floor(cpus * 0.92 * profile.node_mem_mib / profile.node_cpus)
        mem_directive = f"--mem={max(1, share)}"
    else:
        mem_directive = "--mem-per-cpu=4000"

    # --- time ---
    req_min, _ = _time_to_minutes(req.time)
    if req_min is None:
        req_min = 24 * 60
    if profile.max_time_min and not profile.max_time_unlimited:
        eff = min(req_min, profile.max_time_min)
        if eff < req_min:
            warnings.append(f"--time capped to partition MaxTime ({_minutes_to_hms(eff)})")
    else:
        eff = req_min
    time_str = _minutes_to_hms(eff)

    # --- array chunking + throttle ---
    mas = profile.max_array_size if profile.max_array_size is not None else 1001
    chunk = 1 if mas == 0 else mas
    if profile.max_submit_jobs:
        chunk = min(chunk, profile.max_submit_jobs)
    chunk = max(1, chunk)
    chunks = [(b, min(chunk, n_elements - b)) for b in range(0, n_elements, chunk)]

    if req.max_concurrent:
        throttle = req.max_concurrent
    elif profile.max_submit_jobs:
        throttle = min(profile.max_submit_jobs, n_elements)
    else:
        throttle = min(n_elements, 16)

    # --- partition (emit only if known) ---
    partition = profile.partition

    # --- directives (first chunk span) ---
    first_count = chunks[0][1] if chunks else n_elements
    array_spec = f"0-{first_count - 1}%{throttle}"
    directives = [
        f"#SBATCH --job-name={req.job_name}",
        f"#SBATCH --array={array_spec}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
    ]
    if exclusive:
        directives.append("#SBATCH --exclusive")
    directives.append(f"#SBATCH {mem_directive}")
    directives.append(f"#SBATCH --time={time_str}")
    directives.append(f"#SBATCH --output={req.outdir}/%x_%A_%a.out")
    directives.append(f"#SBATCH --error={req.outdir}/%x_%A_%a.err")
    if partition:
        directives.append(f"#SBATCH --partition={partition}")
    if req.account or profile.account:
        directives.append(f"#SBATCH --account={req.account or profile.account}")
    if req.qos or profile.qos:
        directives.append(f"#SBATCH --qos={req.qos or profile.qos}")

    # --- transparency header ---
    overrides = [name for name, val in [
        ("partition", req.partition), ("account", req.account), ("qos", req.qos),
        ("cpus", req.cpus_per_task), ("mem", req.mem),
        ("combos-per-node", req.combos_per_node),
        ("max-concurrent", req.max_concurrent), ("pack", req.pack if req.pack != "auto" else None),
    ] if val]
    header = [
        "Auto-generated by PrinTE guided_search (slurm_sbatch).",
        f"probe: {'OK' if profile.slurm_available else ('disabled (--slurm-no-probe)' if req.no_probe else 'FAILED -> portable fallback')}; "
        f"partition={partition or '(admin default)'} ({profile.partition_source})",
        f"exclusivity: {profile.exclusivity} [{profile.exclusivity_signal}]",
        f"node(min): cpus={profile.node_cpus} cores={profile.node_cores} "
        f"mem={profile.node_mem_mib}MiB threads/core={profile.threads_per_core}",
        f"mode={mode} K(combos/element)={K} elements={n_elements} throttle=%{throttle}",
        f"cpus-per-task={cpus} mem={mem_directive} time={time_str} exclusive={exclusive}",
        f"chunks={len(chunks)}" + ("  (>1: use submit_chunks.sh)" if len(chunks) > 1 else ""),
        f"overrides: {', '.join(overrides) if overrides else '(none)'}",
    ]
    for w in warnings:
        header.append(f"WARN: {w}")

    return SbatchPlan(
        mode=mode, combos_per_element=K, n_combos=req.n_combos, n_elements=n_elements,
        array_throttle=throttle, cpus_per_task=cpus, exclusive=exclusive,
        mem_directive=mem_directive, time_str=time_str, partition=partition,
        account=req.account or profile.account, qos=req.qos or profile.qos,
        chunks=chunks, directives=directives, header_comments=header,
        degraded=degraded, warnings_out=warnings,
    )


def _runner_passthrough(args) -> str:
    """The PrinTE pass-through arg tail (everything except --combo-index/--threads)."""
    clean_lib = getattr(args, "clean_lib", None)
    lib_arg = (f'--clean-lib "{clean_lib}"' if clean_lib
               else f'--te-lib "{args.te_lib}"')
    parts = [
        f'--printe-script "{args.printe_script}"',
        f'--ge "{args.ge}"', f'--st "{args.st}"', f'--mut "{args.mut}"',
        f'--mxgs "{args.mxgs}"', f'--mngs "{args.mngs}"', f'--tstv "{args.tstv}"',
        f'--bed "{args.bed}"', f'--fasta "{args.fasta}"',
        lib_arg, f'--ratios "{args.ratios}"',
    ]
    return " \\\n    ".join(parts)


def _header_block(plan: "SbatchPlan") -> str:
    return "\n".join(f"# {ln}" for ln in plan.header_comments)


def render_sbatch_script(plan: "SbatchPlan", args, runner_script: Path,
                         combos_tsv: Path) -> str:
    """Render the full sbatch script text for `plan` (single or packed body)."""
    runner = str(runner_script)
    tsv = str(combos_tsv)
    passthrough = _runner_passthrough(args)
    directives = "\n".join(plan.directives)
    header = _header_block(plan)

    if plan.mode == "single":
        body = f"""set -euo pipefail

echo "[INFO] Host: $(hostname)"
echo "[INFO] Date: $(date)"
echo "[INFO] PWD:  $(pwd)"
GI=$(( ${{SLURM_ARRAY_TASK_ID:-0}} + ${{BASE:-0}} ))
echo "[INFO] Combo index: $GI"
THREADS="${{SLURM_CPUS_PER_TASK:-{plan.cpus_per_task}}}"

python "{runner}" \\
    --run-one \\
    --combo-file "{tsv}" \\
    --combo-index "$GI" \\
    {passthrough} \\
    --threads "$THREADS"
"""
    else:
        body = f"""set -uo pipefail            # NOT -e: collect every combo's status

LOGDIR="{plan_outdir(plan)}/${{SLURM_ARRAY_JOB_ID:-job}}_${{SLURM_ARRAY_TASK_ID:-0}}"
mkdir -p "$LOGDIR"
THREADS="${{SLURM_CPUS_PER_TASK:-{plan.cpus_per_task}}}"
CPUS=$( (command -v python3 >/dev/null && \\
        python3 -c 'import os;print(len(os.sched_getaffinity(0)))') 2>/dev/null || nproc)
E=$(( ${{SLURM_ARRAY_TASK_ID:-0}} + ${{BASE:-0}} ))
K={plan.combos_per_element}
N={plan.n_combos}
JOBS=$(( CPUS/THREADS > 0 ? CPUS/THREADS : 1 ))
echo "[INFO] Host $(hostname): element $E packing up to $K combo(s), $THREADS thread(s) each"

run_one() {{
  python "{runner}" \\
    --run-one \\
    --combo-file "{tsv}" \\
    --combo-index "$1" \\
    {passthrough} \\
    --threads "$THREADS"
}}
export -f run_one
export THREADS

# Fast path: GNU parallel (joblog + halt-on-fail) if available
if command -v parallel >/dev/null 2>&1; then
  for ((j=0; j<K; j++)); do gi=$((E*K+j)); [ "$gi" -lt "$N" ] && echo "$gi"; done | \\
    parallel --will-cite -j "$JOBS" --joblog "$LOGDIR/joblog.tsv" --halt soon,fail=1 \\
      'run_one {{}} > '"$LOGDIR"'/combo_{{}}.log 2>&1'
  exit $?
fi

# Portable fallback: background & + per-PID wait + status aggregation
declare -a PIDS GIL
rc=0
trap 'kill -- -$$ 2>/dev/null; exit 143' TERM INT
for ((j=0; j<K; j++)); do
  gi=$((E*K+j)); [ "$gi" -ge "$N" ] && break
  run_one "$gi" > "$LOGDIR/combo_${{gi}}.log" 2>&1 &
  PIDS+=($!); GIL+=("$gi")
done
for i in "${{!PIDS[@]}}"; do
  if wait "${{PIDS[$i]}}"; then
    echo "OK   ${{GIL[$i]}}" >> "$LOGDIR/status.txt"
  else
    crc=$?; rc=1; echo "FAIL ${{GIL[$i]}} rc=$crc" >> "$LOGDIR/status.txt"
  fi
done
echo "[INFO] element $E done rc=$rc"
exit $rc
"""

    return f"""#!/bin/bash
{header}
{directives}

{body}"""


def plan_outdir(plan: "SbatchPlan") -> str:
    """Recover the outdir from the --output directive for packed log subdirs."""
    for d in plan.directives:
        if d.startswith("#SBATCH --output="):
            val = d.split("=", 1)[1]
            return val.rsplit("/", 1)[0]
    return "slurm_logs"


def summarize_plan(plan: "SbatchPlan") -> str:
    """Plain-text summary for stdout (same content as the script header block)."""
    lines = ["[SLURM] " + plan.header_comments[0]]
    lines += ["[SLURM] " + ln for ln in plan.header_comments[1:]]
    return "\n".join(lines)


def render_chunk_submitter(plan: "SbatchPlan", sbatch_path: Path) -> str | None:
    """Return a bash submitter looping sbatch per chunk, or None if single chunk."""
    if len(plan.chunks) <= 1:
        return None
    lines = [
        "#!/bin/bash",
        "# Submit the array in chunks (n_elements exceeds MaxArraySize).",
        "# Each chunk passes BASE so combo indices stay globally correct.",
        "set -euo pipefail",
        "",
    ]
    for base, count in plan.chunks:
        lines.append(
            f"sbatch --export=ALL,BASE={base} "
            f"--array=0-{count - 1}%{plan.array_throttle} {sbatch_path}")
    lines.append("")
    return "\n".join(lines)


def build_request_from_args(args, n_combos: int) -> SbatchRequest:
    """Map an argparse Namespace (from either caller) -> SbatchRequest.

    getattr defaults make this tolerant of callers that lack the newer flags.
    None for cpus/mem means 'auto-size'.
    """
    return SbatchRequest(
        threads=int(args.threads),
        n_combos=n_combos,
        job_name=getattr(args, "slurm_job_name", None) or "printe_grid",
        outdir=getattr(args, "slurm_outdir", None) or "slurm_logs",
        partition=getattr(args, "slurm_partition", None),
        account=getattr(args, "slurm_account", None),
        qos=getattr(args, "slurm_qos", None),
        cpus_per_task=getattr(args, "slurm_cpus", None),
        mem=getattr(args, "slurm_mem", None),
        time=getattr(args, "slurm_time", None) or "24:00:00",
        pack=getattr(args, "slurm_pack", None) or "auto",
        combos_per_node=getattr(args, "slurm_combos_per_node", None),
        max_concurrent=getattr(args, "slurm_max_concurrent", None),
        use_logical_cpus=bool(getattr(args, "slurm_use_logical_cpus", False)),
        no_probe=bool(getattr(args, "slurm_no_probe", False)),
        probe_timeout=int(getattr(args, "slurm_probe_timeout", None) or 8),
    )
