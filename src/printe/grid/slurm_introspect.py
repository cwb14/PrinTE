"""Generation-time SLURM cluster introspection for sbatch generation.

Impure discovery layer: runs sinfo/scontrol/sacctmgr with hard timeouts and
strict validation (exit 0 AND non-empty AND parseable), returning a
ClusterProfile. NEVER raises: any failed probe degrades that field to None and
appends a warning. Parses machine-readable output only (scontrol --oneliner
key=value, sinfo -o '%...'), never human column layouts.
"""

import os
import shutil
import subprocess
from dataclasses import dataclass, field


@dataclass
class ClusterProfile:
    slurm_available: bool = False
    partition: str | None = None
    partition_source: str = "none"          # "user" | "default-detected" | "none"
    exclusivity: str = "unknown"            # "exclusive" | "shared" | "forced" | "unknown"
    exclusivity_signal: str = ""
    node_cpus: int | None = None
    node_cores: int | None = None
    threads_per_core: int | None = None
    node_mem_mib: int | None = None
    max_cpus_per_node: int | None = None
    max_time_min: int | None = None
    max_time_unlimited: bool = False
    default_time_min: int | None = None
    max_array_size: int | None = None
    max_submit_jobs: int | None = None
    account: str | None = None
    qos: str | None = None
    select_type: str | None = None
    warnings: list[str] = field(default_factory=list)
    provenance: dict[str, str] = field(default_factory=dict)


def _run(cmd: list[str], timeout_s: int) -> str | None:
    """Return stripped stdout iff exit 0 AND non-empty. Else None. Never raises."""
    if shutil.which(cmd[0]) is None:
        return None
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    except (subprocess.TimeoutExpired, OSError):
        return None
    if res.returncode != 0:
        return None
    out = (res.stdout or "").strip()
    return out or None


def _parse_oneliner(line: str) -> dict[str, str]:
    """Parse whitespace-separated Key=Value tokens (scontrol --oneliner)."""
    out: dict[str, str] = {}
    for tok in line.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def _grep_config(cfg: str, key: str) -> str | None:
    """Find 'Key = Value' in `scontrol show config` output (case-insensitive)."""
    kl = key.lower()
    for ln in cfg.splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            if k.strip().lower() == kl:
                return v.strip()
    return None


def _strip_plus(s: str) -> str:
    return s.rstrip("+")


def _to_int(s: str | None) -> int | None:
    if s is None:
        return None
    s = _strip_plus(s.strip())
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


def _time_to_minutes(s: str | None) -> tuple[int | None, bool]:
    """Return (minutes, unlimited). minutes None if unparseable/unlimited/none."""
    if s is None:
        return None, False
    t = s.strip()
    if t.upper() in {"UNLIMITED", "INFINITE"}:
        return None, True
    if t.upper() in {"NONE", "N/A", ""}:
        return None, False
    try:
        days = 0
        if "-" in t:
            d, t = t.split("-", 1)
            days = int(d)
        parts = [int(p) for p in t.split(":")]
        if len(parts) == 3:
            h, m, sec = parts
        elif len(parts) == 2:
            h, m, sec = 0, parts[0], parts[1]
        elif len(parts) == 1:
            h, m, sec = 0, parts[0], 0
        else:
            return None, False
        return days * 1440 + h * 60 + m + (1 if sec else 0), False
    except (ValueError, TypeError):
        return None, False


def classify_partition(oversub: str | None, exclusive_field: str | None,
                       select_type: str | None) -> tuple[str, str]:
    """Map SLURM signals -> (class, deciding-signal-explanation).

    class in {"exclusive","shared","forced","unknown"}. Precedence (verified):
    FORCE > OverSubscribe=EXCLUSIVE / Exclusive=NODE|USER|TOPO > SelectType.
    Unknown SelectType => "shared" (safe: under-pack beats over-pack).
    """
    os_ = (oversub or "").upper()
    exf = (exclusive_field or "").upper()
    sel = (select_type or "").lower()
    if os_.startswith("FORCE"):
        return "forced", "OverSubscribe=FORCE (force-shared; never pack)"
    if os_ == "EXCLUSIVE":
        return "exclusive", "OverSubscribe=EXCLUSIVE"
    if exf in {"NODE", "USER", "TOPO"}:
        return "exclusive", f"Exclusive={exf}"
    if not sel:
        return "shared", "SelectType unknown -> safe default (shared)"
    if sel == "select/linear" and os_ in {"NO", "YES", ""}:
        return "exclusive", "SelectType=select/linear (whole node per job)"
    if sel in {"select/cons_tres", "select/cons_res"}:
        return "shared", f"SelectType={sel} (cores co-scheduled)"
    return "unknown", "unrecognized SelectType/OverSubscribe signals"


def _detect_default_partition(timeout_s: int) -> str | None:
    out = _run(["sinfo", "-h", "-o", "%P"], timeout_s)
    if not out:
        return None
    names = [ln.strip() for ln in out.splitlines() if ln.strip()]
    for n in names:
        if n.endswith("*"):
            return n.rstrip("*")
    return None  # no default marked; caller leaves partition unset


def _discover_accounts(profile: "ClusterProfile", timeout_s: int) -> None:
    user = os.environ.get("USER", "")
    if not user:
        return
    out = _run(["sacctmgr", "-nP", "show", "assoc", f"user={user}",
                "format=Account,QOS,DefaultQOS,MaxSubmitJobs"], timeout_s)
    if not out:
        return
    first = out.splitlines()[0].split("|")
    # format=Account,QOS,DefaultQOS,MaxSubmitJobs
    if len(first) >= 1 and first[0].strip():
        profile.account = first[0].strip()
        profile.provenance["account"] = "sacctmgr show assoc"
    if len(first) >= 3 and first[2].strip():
        profile.qos = first[2].strip()
    elif len(first) >= 2 and first[1].strip():
        profile.qos = first[1].split(",")[0].strip()
    if len(first) >= 4:
        msj = _to_int(first[3])
        if msj:
            profile.max_submit_jobs = msj
            profile.provenance["max_submit_jobs"] = "sacctmgr show assoc"


def probe_cluster(requested_partition: str | None, *, timeout_s: int = 8,
                  do_probe: bool = True) -> ClusterProfile:
    """Probe the cluster (login node) and return a ClusterProfile. Never raises."""
    p = ClusterProfile()

    if not do_probe:
        p.warnings.append("Probing disabled (--slurm-no-probe); portable defaults.")
        if requested_partition:
            p.partition, p.partition_source = requested_partition, "user"
        return p

    # 1. Effective partition
    if requested_partition:
        p.partition, p.partition_source = requested_partition, "user"
    else:
        default = _detect_default_partition(timeout_s)
        if default:
            p.partition, p.partition_source = default, "default-detected"
            p.provenance["partition"] = "sinfo -h -o '%P' (default '*')"

    # 2. Cluster config: SelectType, MaxArraySize
    cfg = _run(["scontrol", "show", "config"], timeout_s)
    if cfg:
        p.slurm_available = True
        sel = _grep_config(cfg, "SelectType")
        if sel:
            p.select_type = sel
            p.provenance["select_type"] = "scontrol show config"
        mas = _grep_config(cfg, "MaxArraySize")
        mas_i = _to_int(mas) if mas else None
        if mas_i is not None:
            p.max_array_size = mas_i
            p.provenance["max_array_size"] = "scontrol show config"

    # 3. Partition policy + 4. smallest node — only if we know the partition
    if p.partition:
        line = _run(["scontrol", "show", "partition", p.partition, "--oneliner"], timeout_s)
        oversub = exclusive_field = None
        if line:
            p.slurm_available = True
            d = _parse_oneliner(line)
            oversub = d.get("OverSubscribe") or d.get("Shared")
            exclusive_field = d.get("Exclusive")
            mt, mt_unl = _time_to_minutes(d.get("MaxTime"))
            p.max_time_min, p.max_time_unlimited = mt, mt_unl
            p.default_time_min = _time_to_minutes(d.get("DefaultTime"))[0]
            mcpn = d.get("MaxCPUsPerNode")
            if mcpn and mcpn.upper() not in {"UNLIMITED", "INFINITE"}:
                p.max_cpus_per_node = _to_int(mcpn)
            p.provenance["exclusivity"] = "scontrol show partition --oneliner"
        else:
            h = _run(["sinfo", "-h", "-p", p.partition, "-o", "%h"], timeout_s)
            if h:
                p.slurm_available = True
                oversub = h.splitlines()[0].strip()
                p.provenance["exclusivity"] = "sinfo -o '%h'"

        cls, sig = classify_partition(oversub, exclusive_field, p.select_type)
        p.exclusivity, p.exclusivity_signal = cls, sig

        nodes = _run(["sinfo", "-e", "-h", "-p", p.partition, "-o", "%c %m %X %Y %Z"], timeout_s)
        if nodes:
            p.slurm_available = True
            min_cpus = min_mem = cores = tpc = None
            for ln in nodes.splitlines():
                f = ln.split()
                if len(f) < 5:
                    continue
                c, m = _to_int(f[0]), _to_int(f[1])
                x, y, z = _to_int(f[2]), _to_int(f[3]), _to_int(f[4])
                if None in (c, m, x, y, z):
                    continue
                if min_cpus is None or c < min_cpus:
                    min_cpus, cores, tpc = c, x * y, z
                if min_mem is None or m < min_mem:
                    min_mem = m
            p.node_cpus, p.node_cores, p.threads_per_core, p.node_mem_mib = \
                min_cpus, cores, tpc, min_mem
            p.provenance["node_sizing"] = "sinfo -e -o '%c %m %X %Y %Z' (smallest)"

    # 5. Accounts / QOS / submit limits (best-effort)
    _discover_accounts(p, timeout_s)

    if not p.slurm_available:
        p.warnings.append("No SLURM introspection succeeded; using portable defaults.")
    return p
