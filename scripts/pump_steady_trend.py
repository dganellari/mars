#!/usr/bin/env python3
"""Compact the per-step pump diagnostics into the few trends that decide
whether the run is reaching steady state.

The raw log is far too big to move around, and the interesting part is always
the tail. This prints one row per report plus the three derived quantities:

  gain      u_rms / Q_out -- plateaus if the startup transient is converging,
                             keeps climbing if it is not. Q_out tracks the ramp
                             exactly, so this divides out the forcing.
  divRCmax  stabilized divergence max; the one diagnostic that did NOT simply
            scale with the ramp in the 300-step probe.
  max/p999  flat => the peak is a coherent flow feature, climbing => a local
            spike is forming. At >1 rank the percentiles are MPI_MAX-reduced,
            so treat the level as a lower bound and read the trend.

usage: pump_steady_trend.py LOG [--every N]
"""
import re, sys

def num(pat, s):
    m = re.search(pat, s)
    return float(m.group(1)) if m else None

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    every = 1
    for a in sys.argv[1:]:
        if a.startswith("--every="):
            every = int(a.split("=")[1])
    if not args:
        sys.exit(__doc__)

    rows, pend, cur = [], {}, None
    for line in open(args[0], errors="replace"):
        if "[u-profile]" in line:
            pend["p999"] = num(r"p999=([-\d.eE+]+)", line)
            pend["umax_p"] = num(r"max=([-\d.eE+]+)", line)
            pend["ratio"] = num(r"max/p999=([-\d.eE+]+)", line)
            pend["hratio"] = num(r"h_peak/h_med=([-\d.eE+]+)", line)
            m = re.search(r"peak-node=(\w+)", line)
            pend["peak"] = m.group(1) if m else "?"
        elif "[outlet-continuity]" in line:
            pend["ocrms"] = num(r"rms=([-\d.eE+]+)", line)
            pend["oc_n"] = num(r"corrections=(\d+)", line)
        elif line.startswith("Step"):
            cur = dict(pend)
            pend = {}
            cur["step"] = int(line.split()[1])
            for k, p in (("u_rms", r"u_rms=([-\d.eE+]+)"), ("u_max", r"u_max=([-\d.eE+]+)"),
                         ("dur", r"d\(u_rms\)=([-\d.eE+]+)"), ("divRC", r"divRC\*L/U=([-\d.eE+]+)"),
                         ("divRCrms", r"divRCrms\*L/U=([-\d.eE+]+)"), ("cg_p", r"cg_p=(-?\d+)")):
                cur[k] = num(p, line)
            rows.append(cur)
        elif cur is not None and "[bc-sanity]" in line:
            cur["Q_out"] = num(r"Q_out=([-\d.eE+]+)", line)
        elif cur is not None and "[mass-balanceRC]" in line:
            cur["imbRC"] = num(r"imbalance=([-\d.eE+]+)", line)

    if not rows:
        sys.exit("no Step records found")

    hdr = ("step", "u_rms", "u_max", "gain", "d(u_rms)", "divRCmax", "divRCrms",
           "max/p999", "peak", "h_p/h_m", "imbRC%", "cg_p")
    print("".join(f"{h:>10}" for h in hdr))
    for r in rows[::every]:
        g = (r["u_rms"] / r["Q_out"]) if r.get("Q_out") else None
        def f(v, spec="10.3g"):
            return format(v, spec) if isinstance(v, float) else f"{str(v):>10}"
        print(f"{r['step']:10d}{f(r['u_rms'])}{f(r['u_max'])}{f(g)}{f(r['dur'])}"
              f"{f(r['divRC'])}{f(r['divRCrms'])}{f(r.get('ratio'))}"
              f"{str(r.get('peak','?')):>10}{f(r.get('hratio'))}"
              f"{f(r.get('imbRC'))}{f(r['cg_p'])}")

    a, b = rows[0], rows[-1]
    if a.get("Q_out") and b.get("Q_out"):
        ga, gb = a["u_rms"]/a["Q_out"], b["u_rms"]/b["Q_out"]
        print(f"\ngain {ga:.1f} -> {gb:.1f} over steps {a['step']}..{b['step']}"
              f"   ({100*(gb-ga)/ga:+.1f}%)")
    print(f"divRCmax {a['divRC']} -> {b['divRC']};  "
          f"max/p999 {a.get('ratio')} -> {b.get('ratio')}")

main()
