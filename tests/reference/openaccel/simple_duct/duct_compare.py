#!/usr/bin/env python3
"""Rectangular-duct SIMPLE comparator (contract: README.md, "Mathematical contract").

  run    PREFIX --mesh MESH.json --ranks R
                                   one SIMPLE result (PREFIX-fields.csv or the parts listed by
                                   PREFIX-fields.json, PREFIX-metrics.csv, PREFIX.log, and the exit
                                   status in PREFIX.exit or PREFIX.run.json) against the analytic
                                   duct solution; MESH.json must sit next to the Exodus file it
                                   describes
  study  RUNDIR --levels 8,16,32 --ranks 1,2,4
                                   every RUNDIR/duct-<cells>-<ranks> result (mesh RUNDIR/duct-<cells>.json):
                                   each run, rank parity per level and refinement across levels

Per run: complete evidence (Exodus file matching its SHA-256, fields, metrics, a log with the
controls and CONVERGED, a recorded exit status of 0); the run's identity (the log's rank count
equals the label, its advection scheme equals --advection, and it ends at the metrics' final
iteration, with one header and one final line); converged; the comparison window starts beyond the a-priori entrance
length and ends before the outlet's upstream influence; two a-posteriori indicators hold: the
window sections agree and the section-mean pressure is linear, each to 10% of the measured error.
The indicators detect entrance or outlet transients that vary across the window; they are not a
bound on contamination (a transient nearly uniform over the window would pass them).
Reported: profile errors of the window's reference section against the Boussinesq series (L2 and max), the least-squares pressure gradient against G = 3 mu U/(a^2 K),
transverse velocity, wall slip, flow rate, and the measured development and outlet lengths.
Study: P-rank fields equal the one-rank fields to --rank-tol (velocity/U, pressure/(rho U^2));
errors decrease monotonically under refinement with observed order >= 0.7 on the finest pair;
the analytic G and profile lie inside the grid-convergence-index band (Roache) that the three
finest levels estimate from the solutions alone. Exit 0
only on PASS. Standard library only, Python 3.6 compatible.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duct_analytic as da  # noqa: E402
from duct_mesh import Lattice  # noqa: E402

FIELD_COLUMNS = ["node", "x", "y", "z", "u", "v", "w", "p"]
METRIC_COLUMNS = ["iteration", "momentum", "continuity", "mass_balance", "du", "dp", "dflux", "cancellation",
                  "inlet_kg_s", "outlet_kg_s", "umax_m_s", "closed_faces", "changed_faces"]
# Dunavant degree-4 rule: (weight, barycentric) on the reference triangle, weights sum to 1.
_A, _B, _C, _D = 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.816847572980459
QUADRATURE = [(0.223381589678011, (_B, _A, _A)), (0.223381589678011, (_A, _B, _A)), (0.223381589678011, (_A, _A, _B)),
              (0.109951743655322, (_D, _C, _C)), (0.109951743655322, (_C, _D, _C)), (0.109951743655322, (_C, _C, _D))]
STOKES_FACTOR = 1e4      # a-priori margin: slowest Stokes mode decays by 1e4 before the window
INDICATOR = 0.1          # window indicators allowed, relative to the measured error
ORDER_MIN = 0.7          # observed order on the finest pair (first-order upwind, pre-asymptotic)
GCI_SAFETY = 1.25        # Roache's factor for three levels with the observed order
ORDER_RANGE = (0.5, 3.0) # three-level observed orders accepted for a GCI band
FLOOR = 1e-9             # errors below this count as zero (exact synthetic fields)


ADVECTION = ("upwind", "high-resolution")


class Controls(object):
    def __init__(self, rho=1.0, mu=0.1, inlet_velocity=0.1, residual_tol=1e-8, mass_tol=1e-8,
                 change_tol=1e-8, rank_tol=1e-6, advection="upwind"):
        self.rho, self.mu, self.U = rho, mu, inlet_velocity
        self.residual_tol, self.mass_tol, self.change_tol, self.rank_tol = residual_tol, mass_tol, change_tol, rank_tol
        self.advection = advection
        if not all(math.isfinite(v) and v > 0 for v in (rho, mu, inlet_velocity, residual_tol, mass_tol, change_tol, rank_tol)):
            raise ValueError("controls and tolerances must be finite and positive")
        if advection not in ADVECTION:
            raise ValueError("advection must be one of %s" % ", ".join(ADVECTION))


# ------------------------------------------------------------------ inputs
def load_mesh(path):
    with open(path) as f:
        info = json.load(f)
    if info.get("format") != "mars-simple-duct-v1":
        raise ValueError("%s: not a duct_mesh.py description" % path)
    lattice = Lattice(info["cells"], info["length"], info["width"], info["height"], info["stretch"])
    if (lattice.nx, lattice.ny, lattice.nz, lattice.nodes) != (info["nx"], info["ny"], info["nz"], info["nodes"]):
        raise ValueError("%s: lattice sizes disagree with its parameters" % path)
    # The runs read the Exodus file; without it the description proves nothing about their mesh.
    name = info.get("exodus")
    if not name or os.path.basename(name) != name:
        raise ValueError("%s: no Exodus file named" % path)
    exo = os.path.join(os.path.dirname(path), name)
    if not os.path.isfile(exo):
        raise ValueError("%s: Exodus file %s missing; its SHA-256 cannot be verified" % (path, exo))
    h = hashlib.sha256()
    with open(exo, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    if h.hexdigest() != info.get("sha256"):
        raise ValueError("%s: SHA-256 differs from %s" % (exo, path))
    return lattice


def field_files(prefix, lattice):
    """PREFIX-fields.csv (--field-output gathered), else the per-rank parts listed by
    PREFIX-fields.json (--field-output distributed, format mars-simple-fields-v1)."""
    if os.path.isfile(prefix + "-fields.csv"):
        return [prefix + "-fields.csv"]
    manifest = prefix + "-fields.json"
    if not os.path.isfile(manifest):
        raise ValueError("%s: no -fields.csv or -fields.json" % prefix)
    with open(manifest) as f:
        info = json.load(f)
    parts = info.get("parts")
    if info.get("format") != "mars-simple-fields-v1" or info.get("nodes") != lattice.nodes or not parts \
            or any(not isinstance(q, str) or os.path.basename(q) != q for q in parts):
        raise ValueError("%s: not a distributed field manifest for %d nodes" % (manifest, lattice.nodes))
    return [os.path.join(os.path.dirname(manifest), q) for q in parts]


def load_fields(prefix, lattice):
    """u, v, w, p by lattice node from every field file of PREFIX; each node exactly once, and
    coordinates must be the lattice's."""
    n = lattice.nodes
    u, v, w, p = [None] * n, [0.0] * n, [0.0] * n, [0.0] * n
    scale = max(lattice.length, lattice.width, lattice.height)
    worst = 0.0
    paths = field_files(prefix, lattice)
    for path in paths:
        with open(path, newline="") as f:
            reader = csv.reader(f)
            if next(reader, None) != FIELD_COLUMNS:
                raise ValueError("%s: expected columns %s" % (path, ",".join(FIELD_COLUMNS)))
            for row in reader:
                if len(row) != 8:
                    raise ValueError("%s: malformed row" % path)
                g = int(float(row[0]))
                if g < 0 or g >= n or float(row[0]) != g or u[g] is not None:
                    raise ValueError("%s: invalid or duplicate node %s" % (path, row[0]))
                values = [float(x) for x in row[1:]]
                if not all(math.isfinite(x) for x in values):
                    raise ValueError("%s: nonfinite value at node %d" % (path, g))
                i, j, k = lattice.ijk(g)
                worst = max(worst, abs(values[0] - lattice.x(i)), abs(values[1] - lattice.y(j)), abs(values[2] - lattice.z(k)))
                u[g], v[g], w[g], p[g] = values[3:]
    if any(x is None for x in u):
        raise ValueError("%s: %d lattice nodes missing" % (prefix, sum(x is None for x in u)))
    if worst > 1e-12 * scale:
        raise ValueError("%s: coordinates differ from the duct lattice by %.3e" % (prefix, worst))
    return u, v, w, p


def load_metrics(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        if next(reader, None) != METRIC_COLUMNS:
            raise ValueError("%s: expected columns %s" % (path, ",".join(METRIC_COLUMNS)))
        last = None
        for row in reader:
            if len(row) != len(METRIC_COLUMNS):
                raise ValueError("%s: malformed row" % path)
            last = row
    if last is None:
        raise ValueError("%s: no iterations" % path)
    out = dict(zip(METRIC_COLUMNS, (float(x) for x in last)))
    if not all(math.isfinite(x) for x in out.values()):
        raise ValueError("%s: nonfinite final metrics" % path)
    return out


HEADER = re.compile(r"^SIMPLE Tet4\b.*?\b(\d+) ranks\b.*?,\s*(upwind|high-resolution),\s*laminar")
FINAL = re.compile(r"^(CONVERGED|NOT CONVERGED)\b.*?\biterations=(\d+) ranks=(\d+)")


def read_log(path):
    """What the solver says about its own run: {'verdict', 'controls', 'headers', 'finals'}.

    verdict is 'CONVERGED', 'NOT CONVERGED' (also on ERROR) or None; controls come from the
    'rho=... mu=... inlet_speed=...' line; headers lists (ranks, advection) of every 'SIMPLE Tet4'
    line and finals (verdict, iterations, ranks) of every final line, so concatenated or
    truncated logs are visible to the caller."""
    out = {"verdict": None, "controls": {}, "headers": [], "finals": []}
    with open(path, errors="replace") as f:
        for line in f:
            header, final = HEADER.match(line), FINAL.match(line)
            if header:
                out["headers"].append((int(header.group(1)), header.group(2)))
            elif final:
                out["finals"].append((final.group(1), int(final.group(2)), int(final.group(3))))
                out["verdict"] = final.group(1)
            elif line.startswith("NOT CONVERGED") or line.startswith("ERROR"):
                out["verdict"] = "NOT CONVERGED"
            elif line.startswith("CONVERGED"):
                out["verdict"] = "CONVERGED"
            elif line.startswith("rho="):
                out["controls"] = dict((k, float(v)) for k, v in re.findall(r"(\w+)=([-+0-9.eE]+)", line))
    if any(v == "NOT CONVERGED" for v, _, _ in out["finals"]):
        out["verdict"] = "NOT CONVERGED"
    return out


def read_exit(prefix):
    """Recorded exit statuses of the run: PREFIX.exit (a single integer, written by the GPU
    recipe) and/or the 'exit' of PREFIX.run.json (run_host_study.py)."""
    statuses = []
    if os.path.isfile(prefix + ".exit"):
        with open(prefix + ".exit") as f:
            text = f.read().strip()
        if not re.match(r"^-?\d+$", text):
            raise ValueError("%s.exit: expected one integer, found %r" % (prefix, text[:40]))
        statuses.append(int(text))
    if os.path.isfile(prefix + ".run.json"):
        with open(prefix + ".run.json") as f:
            status = json.load(f).get("exit")
        if not isinstance(status, int):
            raise ValueError("%s.run.json: no integer exit status" % prefix)
        statuses.append(status)
    return statuses


# ------------------------------------------------------------------ section geometry
class Section(object):
    """Triangulated cross-section (the Kuhn hex-face diagonals) and the analytic profile on it."""
    _cache = {}

    def __init__(self, lattice, duct, gradient):
        ny, nz = lattice.ny, lattice.nz
        self.count = (ny + 1) * (nz + 1)
        self.ys = [lattice.y(j) for j in range(ny + 1)]
        self.zs = [lattice.z(k) for k in range(nz + 1)]
        self.center = ny // 2 + (ny + 1) * (nz // 2)
        self.wall = [j in (0, ny) or k in (0, nz) for k in range(nz + 1) for j in range(ny + 1)]
        self.exact = [duct.velocity(y, z, gradient) for z in self.zs for y in self.ys]
        self.triangles = []   # (node0, node1, node2, area, analytic at quadrature points)
        for k in range(nz):
            for j in range(ny):
                a, b, c, d = j + (ny + 1) * k, j + 1 + (ny + 1) * k, j + 1 + (ny + 1) * (k + 1), j + (ny + 1) * (k + 1)
                for tri in ((a, b, c), (a, d, c)):
                    pts = [(self.ys[t % (ny + 1)], self.zs[t // (ny + 1)]) for t in tri]
                    area = 0.5 * abs((pts[1][0] - pts[0][0]) * (pts[2][1] - pts[0][1]) - (pts[2][0] - pts[0][0]) * (pts[1][1] - pts[0][1]))
                    q = [duct.velocity(sum(l * pt[0] for l, pt in zip(lam, pts)), sum(l * pt[1] for l, pt in zip(lam, pts)), gradient)
                         for _, lam in QUADRATURE]
                    self.triangles.append((tri, area, q))
        self.area = sum(t[1] for t in self.triangles)
        self.exact_l2 = math.sqrt(sum(area * sum(wq * uq * uq for (wq, _), uq in zip(QUADRATURE, q))
                                      for _, area, q in self.triangles))

    @classmethod
    def get(cls, lattice, duct, gradient):
        key = (lattice.nx, lattice.ny, lattice.nz, lattice.length, lattice.width, lattice.height, duct.viscosity, gradient)
        if key not in cls._cache:
            cls._cache[key] = cls(lattice, duct, gradient)
        return cls._cache[key]

    def mean(self, values):
        return sum(area * (values[a] + values[b] + values[c]) / 3 for (a, b, c), area, _ in self.triangles) / self.area

    def l2_error(self, values):
        total = 0.0
        for (a, b, c), area, q in self.triangles:
            for (wq, lam), uq in zip(QUADRATURE, q):
                e = lam[0] * values[a] + lam[1] * values[b] + lam[2] * values[c] - uq
                total += area * wq * e * e
        return math.sqrt(total) / self.exact_l2


def plane(values, lattice, i):
    return values[i::lattice.nx + 1]   # node i + (nx+1) t, t = j + (ny+1) k


# ------------------------------------------------------------------ one run
def default_window(lattice, duct, controls):
    """A-priori window: Stokes decay by STOKES_FACTOR and twice the Durst length from the inlet,
    Stokes decay from the outlet; rounded inwards to multiples of H/2 so every level shares it."""
    reynolds = controls.rho * controls.U * duct.hydraulic_diameter / controls.mu
    start = max(da.stokes_decay_length(duct, STOKES_FACTOR), 2 * da.durst_entrance_length(reynolds, duct.hydraulic_diameter))
    end = lattice.length - da.stokes_decay_length(duct, STOKES_FACTOR)
    step = lattice.height / 2
    return math.ceil(start / step - 1e-9) * step, math.floor(end / step + 1e-9) * step, start, end, reynolds


def analyze(prefix, lattice, controls, window=None, ranks=None):
    """Summary dict of one run; 'failures' lists every violated per-run condition. ranks is the
    rank count the run is labelled with; its log must say the same."""
    duct = da.Duct(lattice.width, lattice.height, controls.mu)
    G = duct.pressure_gradient(controls.U)
    uc = duct.centerline(G)
    out = {"prefix": prefix, "cells": lattice.nz, "h": lattice.hz, "failures": [],
           "analytic": {"G": G, "K": duct.shape_factor(), "u_centerline": uc, "fRe": duct.fanning_re(),
                        "G_planar": da.planar_pressure_gradient(lattice.height, controls.mu, controls.U)}}
    fail = out["failures"].append
    try:
        u, v, w, p = load_fields(prefix, lattice)
        metrics = load_metrics(prefix + "-metrics.csv")
    except (OSError, ValueError, StopIteration) as e:
        fail("input: %s" % e)
        return out
    # The log is evidence: the solver's own verdict, identity and the physical controls it ran with.
    try:
        record = read_log(prefix + ".log")
    except OSError as e:
        record = {"verdict": None, "controls": {}, "headers": [], "finals": []}
        fail("log: %s" % e)
    log, used = record["verdict"], record["controls"]
    out["iterations"] = int(metrics["iteration"])
    out["final_metrics"] = metrics
    out["log"] = log
    if log != "CONVERGED":
        fail("log: no CONVERGED line")
    # Run identity: one run per log, the labelled rank count, the compared scheme, and the
    # iteration the metrics end at.
    if len(record["headers"]) != 1 or len(record["finals"]) != 1:
        fail("identity: the log holds %d run headers and %d final lines; expected one of each"
             % (len(record["headers"]), len(record["finals"])))
    else:
        (header_ranks, scheme), (_, iterations, final_ranks) = record["headers"][0], record["finals"][0]
        out["log_identity"] = {"ranks": final_ranks, "advection": scheme, "iterations": iterations}
        if header_ranks != final_ranks or (ranks is not None and final_ranks != ranks):
            fail("identity: labelled %s ranks, the log says %d (header) and %d (final)" % (ranks, header_ranks, final_ranks))
        if scheme != controls.advection:
            fail("identity: the run used %s advection, the comparison expects %s" % (scheme, controls.advection))
        if iterations != out["iterations"]:
            fail("identity: the log ends at iteration %d, the metrics at %d" % (iterations, out["iterations"]))
    # The launcher's exit status: recorded by the GPU recipe (.exit) or run_host_study (.run.json).
    try:
        statuses = read_exit(prefix)
    except (OSError, ValueError) as e:
        statuses = None
        fail("exit: %s" % e)
    if statuses is not None:
        out["exit"] = statuses
        if not statuses:
            fail("exit: no exit status recorded (%s.exit or %s.run.json)" % (prefix, prefix))
        elif any(s != 0 for s in statuses):
            fail("exit: the run exited %s" % ", ".join(str(s) for s in statuses))
    # The analytic solution is for these controls: the run must state them and have used them.
    for key, value in (("rho", controls.rho), ("mu", controls.mu), ("inlet_speed", controls.U)):
        if key not in used:
            fail("controls: the log states no %s (expected a 'rho=... mu=... inlet_speed=...' line)" % key)
        elif not abs(used[key] - value) <= 1e-12 * value:
            fail("controls: the run used %s=%g, the comparison assumes %g" % (key, used[key], value))
    inflow = controls.rho * controls.U * lattice.width * lattice.height
    if not abs(abs(metrics["inlet_kg_s"]) - inflow) <= 1e-9 * inflow:
        fail("inlet mass flux %.12g kg/s, expected rho U W H = %.12g" % (abs(metrics["inlet_kg_s"]), inflow))
    for key, tol in (("momentum", controls.residual_tol), ("continuity", controls.residual_tol),
                     ("mass_balance", controls.mass_tol), ("du", controls.change_tol), ("dp", controls.change_tol),
                     ("dflux", controls.change_tol), ("cancellation", 1e-10)):
        if not metrics[key] <= tol:
            fail("not converged: final %s %.3e > %.1e" % (key, metrics[key], tol))
    if metrics["closed_faces"] != 0 or metrics["changed_faces"] != 0:
        fail("outlet reversal: %d closed, %d changed faces" % (metrics["closed_faces"], metrics["changed_faces"]))

    start, end, start_min, end_max, reynolds = default_window(lattice, duct, controls)
    if window:
        start, end = window
    planes = [i for i in range(lattice.nx + 1) if start - 1e-9 <= lattice.x(i) <= end + 1e-9]
    out["window"] = {"start": start, "end": end, "a_priori_start": start_min, "a_priori_end": end_max,
                     "reynolds_dh": reynolds, "planes": len(planes),
                     "durst_length": da.durst_entrance_length(reynolds, duct.hydraulic_diameter),
                     "stokes_length": da.stokes_decay_length(duct, STOKES_FACTOR)}
    if start < start_min - 1e-9 or end > end_max + 1e-9:
        fail("window [%g, %g] inside the a-priori entrance or outlet region [%.3f, %.3f]" % (start, end, start_min, end_max))
    if len(planes) < 3:
        fail("window holds %d node planes; need 3" % len(planes))
        return out
    ref = min(planes, key=lambda i: abs(lattice.x(i) - 0.5 * (start + end)))
    sec = Section.get(lattice, duct, G)

    # Reference section: profile, transverse velocity, wall slip, flow rate.
    ur, vr, wr = plane(u, lattice, ref), plane(v, lattice, ref), plane(w, lattice, ref)
    err = [abs(a - b) / uc for a, b in zip(ur, sec.exact)]
    out["x_ref"] = lattice.x(ref)
    out["profile_l2"] = sec.l2_error(ur)
    out["profile_max"] = max(err)
    out["profile_max_interior"] = max(e for e, wall in zip(err, sec.wall) if not wall)
    out["wall_slip"] = max(math.sqrt(a * a + b * b + c * c) / uc for a, b, c, wall in zip(ur, vr, wr, sec.wall) if wall)
    out["transverse_max"] = max(math.hypot(b, c) for b, c in zip(vr, wr)) / controls.U
    out["centerline_ratio"] = ur[sec.center] / uc
    out["flow_rate_ratio"] = sec.mean(ur) / controls.U

    # Pressure gradient: least squares of the section-mean pressure over the window.
    xs = [lattice.x(i) for i in planes]
    pm = [sec.mean(plane(p, lattice, i)) for i in planes]
    xbar, pbar = sum(xs) / len(xs), sum(pm) / len(pm)
    slope = sum((x - xbar) * (q - pbar) for x, q in zip(xs, pm)) / sum((x - xbar) ** 2 for x in xs)
    out["G"] = -slope
    out["G_error"] = (-slope - G) / G
    out["pressure_linearity"] = max(abs(q - pbar - slope * (x - xbar)) for x, q in zip(xs, pm)) / (G * (end - start))

    # Development: deviation of every section from the reference section (velocity vector / u_c).
    dev = []
    for i in range(lattice.nx + 1):
        a, b, c = plane(u, lattice, i), plane(v, lattice, i), plane(w, lattice, i)
        dev.append(max(math.sqrt((a[t] - ur[t]) ** 2 + (b[t] - vr[t]) ** 2 + (c[t] - wr[t]) ** 2) for t in range(sec.count)) / uc)
    center = [plane(u, lattice, i)[sec.center] for i in range(lattice.nx + 1)]

    def upstream(eps):
        i = ref
        while i > 0 and dev[i - 1] <= eps:
            i -= 1
        return lattice.x(i)

    def downstream(eps):
        i = ref
        while i < lattice.nx and dev[i + 1] <= eps:
            i += 1
        return lattice.length - lattice.x(i)
    i = ref
    while i > 0 and abs(center[i - 1] / center[ref] - 1) <= 0.01:
        i -= 1
    out["development"] = {"profile_1e-2": upstream(1e-2), "profile_1e-3": upstream(1e-3),
                          "centerline_99": lattice.x(i), "outlet_1e-2": downstream(1e-2), "outlet_1e-3": downstream(1e-3)}
    out["window_section_indicator"] = max(dev[i] for i in planes)
    out["profile_max_window"] = max(max(abs(a - b) / uc for a, b in zip(plane(u, lattice, i), sec.exact)) for i in planes)
    out["profile_l2_window"] = max(sec.l2_error(plane(u, lattice, i)) for i in planes)
    # Indicators, not bounds: they see transients that vary across the window.
    allowed = max(INDICATOR * out["profile_max"], 1e-6)
    if out["window_section_indicator"] > allowed:
        fail("window indicator: sections differ by %.3e u_c > %.3e (entrance or outlet transient in the window)"
             % (out["window_section_indicator"], allowed))
    allowed = max(INDICATOR * abs(out["G_error"]), 1e-6)
    if out["pressure_linearity"] > allowed:
        fail("window indicator: section-mean pressure not linear: %.3e > %.3e" % (out["pressure_linearity"], allowed))
    out["sections"] = [{"x": lattice.x(i), "deviation": dev[i], "u_center": center[i],
                        "p_mean": sec.mean(plane(p, lattice, i))} for i in range(lattice.nx + 1)]
    return out


# ------------------------------------------------------------------ rank parity
def field_difference(a_prefix, b_prefix, lattice, controls):
    a = load_fields(a_prefix, lattice)
    b = load_fields(b_prefix, lattice)
    vel = max(math.sqrt((a[0][g] - b[0][g]) ** 2 + (a[1][g] - b[1][g]) ** 2 + (a[2][g] - b[2][g]) ** 2) for g in range(lattice.nodes))
    pres = max(abs(a[3][g] - b[3][g]) for g in range(lattice.nodes))
    return vel / controls.U, pres / (controls.rho * controls.U ** 2)


# ------------------------------------------------------------------ refinement
def observed_order(coarse, fine, ratio=2.0):
    if coarse <= FLOOR or fine <= FLOOR:
        return None
    return math.log(coarse / fine) / math.log(ratio)


def gci(c, m, f):
    """Three-level (ratio 2) observed order, Roache fine-level band and extrapolated limit.

    band = F_s |f - m| / (2^p - 1) estimates |f - limit| from the solutions alone; None when the
    sequence is not monotone (no asymptotic range)."""
    d1, d2 = c - m, m - f
    if d1 == 0 or d2 == 0 or (d1 > 0) != (d2 > 0):
        return None
    p = math.log(d1 / d2) / math.log(2)
    return p, GCI_SAFETY * abs(f - m) / (2 ** p - 1), f + (f - m) / (2 ** p - 1)


def profile_gci(prefixes, lattices, sec_m, duct, G, x_ref):
    """GCI of the reference-section profile at the medium level's nodes (nested lattices).

    The order comes from the solutions alone: RMS(u_m - u_c) / RMS(u_f - u_m) over the coarse
    nodes. Returns order, band = F_s RMS(u_f - u_m) / (2^p - 1), the finest RMS error against the
    analytic profile and the extrapolated RMS error, all over the medium nodes and / u_c."""
    fields = [load_fields(pre, lat)[0] for pre, lat in zip(prefixes, lattices)]
    lc, lm, lf = lattices
    if any(abs(x_ref / lat.hx - round(x_ref / lat.hx)) > 1e-9 for lat in lattices):
        raise ValueError("reference section x=%g is not a node plane of every level" % x_ref)
    uc = duct.centerline(G)
    sc, sm, sf = [plane(u, lat, int(round(x_ref / lat.hx))) for u, lat in zip(fields, lattices)]
    coarse, fine, rows = [], [], []
    for k in range(lm.nz + 1):
        for j in range(lm.ny + 1):
            um, uf = sm[j + (lm.ny + 1) * k], sf[2 * j + (lf.ny + 1) * 2 * k]
            rows.append((um, uf, sec_m.exact[j + (lm.ny + 1) * k]))
            if j % 2 == 0 and k % 2 == 0:
                coarse.append(um - sc[j // 2 + (lc.ny + 1) * (k // 2)])
                fine.append(uf - um)
    rms = lambda xs: math.sqrt(sum(x * x for x in xs) / len(xs))  # noqa: E731
    if rms(fine) == 0 or rms(coarse) <= rms(fine):
        return None
    p = math.log(rms(coarse) / rms(fine)) / math.log(2)
    step = rms([uf - um for um, uf, _ in rows])
    return {"order": p, "band": GCI_SAFETY * step / (2 ** p - 1) / uc,
            "finest_rms": rms([uf - e for _, uf, e in rows]) / uc,
            "extrapolated_rms": rms([uf + (uf - um) / (2 ** p - 1) - e for um, uf, e in rows]) / uc}


def study(rundir, levels, ranks, controls, window=None, refinement=True):
    levels, ranks = sorted(levels), sorted(ranks)
    report = {"rundir": rundir, "levels": levels, "ranks": ranks, "runs": {}, "parity": {}, "refinement": {}, "failures": []}
    fail = report["failures"].append
    lattices = {}
    for c in levels:
        try:
            lattices[c] = load_mesh(os.path.join(rundir, "duct-%d.json" % c))
        except (OSError, ValueError, KeyError) as e:
            fail("mesh %d: %s" % (c, e))
    if report["failures"]:
        return report
    for c in levels:
        for r in ranks:
            s = analyze(os.path.join(rundir, "duct-%d-%d" % (c, r)), lattices[c], controls, window, ranks=r)
            report["runs"]["%d-%d" % (c, r)] = s
            for f in s["failures"]:
                fail("cells %d, %d ranks: %s" % (c, r, f))
    # Rank parity against the fewest-rank run of each level.
    base = ranks[0]
    for c in levels:
        for r in ranks[1:]:
            key = "%d-%d" % (c, r)
            try:
                dv, dp = field_difference(os.path.join(rundir, "duct-%d-%d" % (c, base)), os.path.join(rundir, "duct-%d-%d" % (c, r)),
                                          lattices[c], controls)
            except (OSError, ValueError) as e:
                fail("parity cells %d, %d vs %d ranks: %s" % (c, r, base, e))
                continue
            a, b = report["runs"]["%d-%d" % (c, base)], report["runs"][key]
            report["parity"][key] = {"velocity": dv, "pressure": dp,
                                     "iterations": [a.get("iterations"), b.get("iterations")]}
            if not (dv <= controls.rank_tol and dp <= controls.rank_tol):
                fail("parity cells %d: %d vs %d ranks differ by velocity/U %.3e, pressure/(rho U^2) %.3e > %.1e"
                     % (c, r, base, dv, dp, controls.rank_tol))
    if not refinement:
        return report
    # Refinement on the fewest-rank runs.
    runs = [report["runs"]["%d-%d" % (c, base)] for c in levels]
    if len(levels) < 3:
        fail("refinement needs three nested levels, got %d" % len(levels))
    if any(levels[k + 1] != 2 * levels[k] for k in range(len(levels) - 1)):
        fail("levels must double: %s" % levels)
    if report["failures"] or any("profile_l2" not in s for s in runs):
        return report
    names = ("profile_l2", "profile_max_interior", "G_error", "wall_slip", "transverse_max", "flow_rate_ratio")
    table = dict((n, [abs(s[n] - (1 if n == "flow_rate_ratio" else 0)) for s in runs]) for n in names)
    report["refinement"]["errors"] = table
    # The P1 section flow rate mixes interpolation and wall-slip terms of opposite sign: reported only.
    for n in names[:-1]:
        for k in range(len(levels) - 1):
            if table[n][k] > FLOOR and not table[n][k + 1] < table[n][k]:
                fail("refinement: %s does not decrease from %d to %d cells (%.3e -> %.3e)"
                     % (n, levels[k], levels[k + 1], table[n][k], table[n][k + 1]))
    orders = {}
    for n in ("profile_l2", "G_error"):
        orders[n] = [observed_order(table[n][k], table[n][k + 1]) for k in range(len(levels) - 1)]
        last = orders[n][-1]
        if last is not None and last < ORDER_MIN:
            fail("refinement: observed order of %s on the finest pair %.2f < %.1f" % (n, last, ORDER_MIN))
    report["refinement"]["orders"] = orders
    # Grid convergence index (Roache): the analytic G and profile must lie inside the fine-level
    # error band that the three finest solutions estimate without the analytic solution.
    duct = da.Duct(lattices[levels[0]].width, lattices[levels[0]].height, controls.mu)
    G = duct.pressure_gradient(controls.U)
    g3 = [s["G"] for s in runs[-3:]]
    band = gci(*g3)
    if band is None:
        fail("refinement: GCI G undefined: G is not monotone over the three finest levels %s" % g3)
    else:
        p, width, limit = band
        report["refinement"]["G_gci"] = {"order": p, "band": width / G, "limit": limit, "limit_error": (limit - G) / G,
                                         "finest_error": (g3[-1] - G) / G}
        if not ORDER_RANGE[0] <= p <= ORDER_RANGE[1]:
            fail("refinement: GCI G order %.2f outside [%.1f, %.1f]" % (p, ORDER_RANGE[0], ORDER_RANGE[1]))
        elif not abs(g3[-1] - G) <= width * (1 + 1e-12):
            fail("refinement: GCI G: finest error %.3e exceeds the three-level band %.3e (limit %.6g, analytic %.6g)"
                 % (abs(g3[-1] - G) / G, width / G, limit, G))
    lc, lm, lf = (lattices[c] for c in levels[-3:])
    try:
        q = profile_gci([os.path.join(rundir, "duct-%d-%d" % (c, base)) for c in levels[-3:]],
                        [lc, lm, lf], Section.get(lm, duct, G), duct, G, runs[-1]["x_ref"])
    except ValueError as e:
        fail("refinement: GCI profile: %s" % e)
        return report
    if q is None:
        fail("refinement: GCI profile undefined: the profile does not converge over the three finest levels")
        return report
    report["refinement"]["profile_gci"] = q
    if not ORDER_RANGE[0] <= q["order"] <= ORDER_RANGE[1]:
        fail("refinement: GCI profile order %.2f outside [%.1f, %.1f]" % (q["order"], ORDER_RANGE[0], ORDER_RANGE[1]))
    elif not q["finest_rms"] <= q["band"] * (1 + 1e-12):
        fail("refinement: GCI profile: finest RMS error %.3e exceeds the three-level band %.3e" % (q["finest_rms"], q["band"]))
    return report


# ------------------------------------------------------------------ reporting
def run_lines(s):
    if "profile_l2" not in s:
        return ["%s: %s" % (s["prefix"], "; ".join(s["failures"]))]
    d = s["development"]
    return ["%s: %s after %d iterations" % (s["prefix"], "PASS" if not s["failures"] else "FAIL", s.get("iterations", -1)),
            "  window [%g, %g] (a priori >= %.3f, <= %.3f; Re_Dh %.3g), reference x=%g"
            % (s["window"]["start"], s["window"]["end"], s["window"]["a_priori_start"], s["window"]["a_priori_end"],
               s["window"]["reynolds_dh"], s["x_ref"]),
            "  profile L2 %.4e  max %.4e  max interior %.4e  wall slip %.4e  transverse %.3e U"
            % (s["profile_l2"], s["profile_max"], s["profile_max_interior"], s["wall_slip"], s["transverse_max"]),
            "  G %.8g vs analytic %.8g (error %+.4e; planar parabola would be %.4g)  linearity %.2e"
            % (s["G"], s["analytic"]["G"], s["G_error"], s["analytic"]["G_planar"], s["pressure_linearity"]),
            "  centerline u/u_c %.6f  flow rate %.6f U  window section indicator %.2e u_c"
            % (s["centerline_ratio"], s["flow_rate_ratio"], s["window_section_indicator"]),
            "  development: profile 1e-2 %.3f, 1e-3 %.3f, centerline 99%% %.3f; outlet influence 1e-2 %.3f, 1e-3 %.3f"
            % (d["profile_1e-2"], d["profile_1e-3"], d["centerline_99"], d["outlet_1e-2"], d["outlet_1e-3"])] + \
           ["  FAIL: %s" % f for f in s["failures"]]


def markdown(report):
    lines = ["# Rectangular duct SIMPLE study", "", "Run directory: `%s`" % report["rundir"], ""]
    if report["runs"]:
        lines += ["| cells | ranks | iterations | profile L2 | max interior | G error | wall slip | transverse/U | window indicator | verdict |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
        for key in sorted(report["runs"], key=lambda k: tuple(int(x) for x in k.split("-"))):
            s = report["runs"][key]
            c, r = key.split("-")
            if "profile_l2" in s:
                lines.append("| %s | %s | %d | %.4e | %.4e | %+.4e | %.4e | %.3e | %.2e | %s |"
                             % (c, r, s["iterations"], s["profile_l2"], s["profile_max_interior"], s["G_error"],
                                s["wall_slip"], s["transverse_max"], s["window_section_indicator"], "FAIL" if s["failures"] else "PASS"))
            else:
                lines.append("| %s | %s | - | - | - | - | - | - | - | FAIL |" % (c, r))
    if report["parity"]:
        lines += ["", "| rank parity | max velocity/U | max pressure/(rho U^2) | iterations |", "|---|---|---|---|"]
        for key in sorted(report["parity"], key=lambda k: tuple(int(x) for x in k.split("-"))):
            q = report["parity"][key]
            lines.append("| cells %s vs %d rank | %.3e | %.3e | %s |" % (key.replace("-", ", ranks "), report["ranks"][0],
                                                                        q["velocity"], q["pressure"], q["iterations"]))
    ref = report.get("refinement", {})
    if "orders" in ref:
        lines += ["", "Observed orders (consecutive levels): profile L2 %s, G %s" % (
            ["%.2f" % o if o is not None else "-" for o in ref["orders"]["profile_l2"]],
            ["%.2f" % o if o is not None else "-" for o in ref["orders"]["G_error"]])]
    if "G_gci" in ref:
        g = ref["G_gci"]
        lines.append("GCI G (3 finest): order %.2f, finest error %+.3e within band %.3e; extrapolated G %.8g (error %+.3e)"
                     % (g["order"], g["finest_error"], g["band"], g["limit"], g["limit_error"]))
    if "profile_gci" in ref:
        q = ref["profile_gci"]
        lines.append("GCI profile (3 finest): order %.2f, finest RMS error %.3e within band %.3e; extrapolated RMS error %.3e"
                     % (q["order"], q["finest_rms"], q["band"], q["extrapolated_rms"]))
    lines += ["", "**%s**" % ("PASS" if not report["failures"] else "FAIL")] + ["- %s" % f for f in report["failures"]]
    return "\n".join(lines) + "\n"


def controls_from(o):
    return Controls(o.rho, o.mu, o.inlet_velocity, o.residual_tol, o.mass_tol, o.change_tol, o.rank_tol, o.advection)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command")
    r = sub.add_parser("run")
    r.add_argument("prefix")
    r.add_argument("--mesh", required=True)
    r.add_argument("--ranks", type=int, required=True, help="rank count the run is labelled with")
    r.add_argument("--summary")
    s = sub.add_parser("study")
    s.add_argument("rundir")
    s.add_argument("--levels", required=True)
    s.add_argument("--ranks", default="1,2,4")
    s.add_argument("--report")
    s.add_argument("--parity-only", action="store_true", help="per-run and rank checks, no refinement")
    for q in (r, s):
        q.add_argument("--rho", type=float, default=1.0)
        q.add_argument("--mu", type=float, default=0.1)
        q.add_argument("--inlet-velocity", type=float, default=0.1)
        q.add_argument("--residual-tol", type=float, default=1e-8)
        q.add_argument("--mass-tol", type=float, default=1e-8)
        q.add_argument("--change-tol", type=float, default=1e-8)
        q.add_argument("--rank-tol", type=float, default=1e-6)
        q.add_argument("--window", type=float, nargs=2, metavar=("START", "END"))
        q.add_argument("--advection", default="upwind", choices=ADVECTION, help="scheme the runs must report")
    o = p.parse_args(argv)
    try:
        controls = controls_from(o)
        if o.command == "run":
            summary = analyze(o.prefix, load_mesh(o.mesh), controls, o.window, ranks=o.ranks)
            print("\n".join(run_lines(summary)))
            if o.summary:
                with open(o.summary, "x") as f:
                    json.dump(summary, f, indent=1, sort_keys=True)
            return 0 if not summary["failures"] else 1
        if o.command == "study":
            report = study(o.rundir, [int(x) for x in o.levels.split(",")], [int(x) for x in o.ranks.split(",")], controls, o.window,
                           not o.parity_only)
            for key in sorted(report["runs"], key=lambda k: tuple(int(x) for x in k.split("-"))):
                print("\n".join(run_lines(report["runs"][key])))
            text = markdown(report)
            print(text)
            if o.report:
                with open(o.report, "x") as f:
                    f.write(text)
                with open(os.path.splitext(o.report)[0] + ".json", "x") as f:
                    json.dump(report, f, indent=1, sort_keys=True)
            return 0 if not report["failures"] else 1
        p.print_help()
        return 2
    except (OSError, ValueError, KeyError) as e:
        print("FAIL: %s" % e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
