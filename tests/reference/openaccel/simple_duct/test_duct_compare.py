#!/usr/bin/env python3
"""Comparator tests on synthetic duct fields with known defects.

A synthetic run is analytic + c h^order error terms, a decaying entrance transient from the
uniform inlet, an outlet transient and O(h) transverse velocity, written in the production
CSV formats. The clean family must PASS; each defect must FAIL for the stated reason: the
planar Poiseuille parabola, a biased or low-order pressure gradient, a flow-rate floor,
entrance or outlet transients inside the window, a window inside the a-priori entrance region,
rank mismatches, malformed or unconverged output, outlet reversal, stagnant transverse
velocity, too few or non-monotone levels, and a mesh that is not the one described.
"""
import json
import math
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import duct_analytic as da  # noqa: E402
import duct_compare as dc  # noqa: E402
from duct_mesh import Lattice  # noqa: E402

U, MU, RHO = 0.1, 0.1, 1.0


def write_mesh(d, cells, lattice):
    info = {"format": "mars-simple-duct-v1", "cells": cells, "nx": lattice.nx, "ny": lattice.ny, "nz": lattice.nz,
            "length": lattice.length, "width": lattice.width, "height": lattice.height, "stretch": lattice.stretch,
            "nodes": lattice.nodes, "exodus": "duct-%d.exo" % cells, "sha256": "0" * 64}
    with open(os.path.join(d, "duct-%d.json" % cells), "w") as f:
        json.dump(info, f)


def synthetic(d, levels=(4, 8, 16), ranks=(1,), profile=0.3, g=-0.2, order=1.0, g_bias=0.0, g_signs=None,
              entrance=0.25, outlet=0.2, outlet_amplitude=0.02, transverse=0.05, transverse_order=1.0,
              planar=False, flow_scale=1.0, edit=None, metrics_edit=None, log="CONVERGED iterations=3000",
              header="rho=1 mu=0.10000000000000001 nu=0.10000000000000001 inlet_speed=0.10000000000000001 (inward normal)"):
    """Write duct-<c>.json and duct-<c>-<r>-{fields,metrics}.csv and .log for every level and rank."""
    duct = da.Duct(2, 1, MU)
    G = duct.pressure_gradient(U)
    for n, c in enumerate(levels):
        lat = Lattice(c)
        write_mesh(d, c, lat)
        scale = lat.hz / 0.25
        e_u = profile * scale ** order
        e_g = g_bias + (g_signs[n] if g_signs else 1) * g * scale ** order
        e_t = transverse * scale ** transverse_order
        ys, zs = [lat.y(j) for j in range(lat.ny + 1)], [lat.z(k) for k in range(lat.nz + 1)]
        exact = [[duct.velocity(y, z, G) for y in ys] for z in zs]
        rows = []
        for k in range(lat.nz + 1):
            for j in range(lat.ny + 1):
                y, z = ys[j], zs[k]
                phi = math.cos(math.pi * y / lat.width) * math.cos(math.pi * z / lat.height)
                base = da.planar_velocity(z, lat.height, U) if planar else exact[k][j]
                developed = flow_scale * (base + e_u * U * phi)
                for i in range(lat.nx + 1):
                    x = lat.x(i)
                    uu = developed + (U - developed) * math.exp(-x / entrance) \
                        + outlet_amplitude * U * phi * math.exp(-(lat.length - x) / outlet)
                    Gp = da.planar_pressure_gradient(lat.height, MU, U) if planar else G
                    rows.append((lat.node(i, j, k), x, y, z, uu, e_t * U * phi * math.sin(math.pi * y / lat.width), 0.0,
                                 Gp * (1 + e_g) * (lat.length - x)))
        rows.sort()
        for r in ranks:
            prefix = os.path.join(d, "duct-%d-%d" % (c, r))
            data = [list(row) for row in rows]
            if edit:
                edit(c, r, data)
            with open(prefix + "-fields.csv", "w") as f:
                f.write("node,x,y,z,u,v,w,p\n")
                for row in data:
                    f.write(",".join(repr(v) if not isinstance(v, str) else v for v in [int(row[0])] + row[1:]) + "\n")
            m = dict((k, 0.0) for k in dc.METRIC_COLUMNS)
            m.update(iteration=3000, momentum=5e-9, continuity=1e-12, mass_balance=1e-13, du=1e-12, dp=1e-11,
                     dflux=1e-13, cancellation=1e-15, inlet_kg_s=-RHO * U * 2, outlet_kg_s=RHO * U * 2, umax_m_s=0.2)
            if metrics_edit:
                metrics_edit(c, r, m)
            with open(prefix + "-metrics.csv", "w") as f:
                f.write(",".join(dc.METRIC_COLUMNS) + "\n0" + ",0" * (len(dc.METRIC_COLUMNS) - 1) + "\n")
                f.write(",".join(repr(m[k]) for k in dc.METRIC_COLUMNS) + "\n")
            with open(prefix + ".log", "w") as f:
                f.write("SIMPLE Tet4, %d ranks\n%s outlet_pressure=0 reference_length=1\n[simple] iteration=0\n%s ranks=%d\n"
                        % (r, header, log, r))


class Study(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="duct-compare-")

    def tearDown(self):
        shutil.rmtree(self.dir)

    def study(self, levels=(4, 8, 16), ranks=(1,), window=None, **kw):
        synthetic(self.dir, levels, ranks, **kw)
        return dc.study(self.dir, list(levels), list(ranks), dc.Controls(RHO, MU, U), window)

    def assertFailsWith(self, report, text):
        self.assertTrue(report["failures"], "expected FAIL")
        self.assertTrue(any(text in f for f in report["failures"]),
                        "no failure mentions %r:\n%s" % (text, "\n".join(report["failures"])))

    # ---- the clean family
    def test_convergent_family_passes(self):
        r = self.study(ranks=(1, 2, 4))
        self.assertEqual(r["failures"], [])
        # c h phi at the nodes plus the O(h^2) P1 interpolation error of the analytic profile.
        self.assertGreater(r["refinement"]["orders"]["profile_l2"][-1], 0.75)
        self.assertLess(r["refinement"]["orders"]["profile_l2"][-1], 1.25)
        self.assertAlmostEqual(r["refinement"]["orders"]["G_error"][-1], 1.0, delta=1e-9)
        self.assertLess(abs(r["refinement"]["G_gci"]["limit_error"]), 1e-9)
        q = r["refinement"]["profile_gci"]
        self.assertLess(q["extrapolated_rms"], 0.1 * q["finest_rms"])
        self.assertLessEqual(q["finest_rms"], q["band"])
        run = r["runs"]["16-1"]
        self.assertAlmostEqual(run["G_error"], -0.2 / 4, delta=1e-12)
        self.assertLess(run["window_contamination"], 1e-4)
        self.assertLess(run["development"]["profile_1e-3"], 2.5)

    def test_second_order_family_passes(self):
        self.assertEqual(self.study(order=2.0, transverse_order=2.0)["failures"], [])

    def test_rank_differences_within_tolerance_pass(self):
        def edit(c, r, rows):
            if r == 4:
                rows[7][4] += 5e-8 * U
                for row in rows:
                    row[7] += 5e-8 * RHO * U * U
        self.assertEqual(self.study(ranks=(1, 2, 4), edit=edit)["failures"], [])

    # ---- wrong physics
    def test_planar_parabola_fails(self):
        r = self.study(planar=True, profile=0.0, g=0.0)
        self.assertFailsWith(r, "does not decrease")
        self.assertFailsWith(r, "GCI G undefined")
        self.assertGreater(abs(r["runs"]["16-1"]["G_error"]), 0.3)

    def test_converging_to_a_biased_pressure_gradient_fails(self):
        # Errors 0.44, 0.24, 0.14 decrease at orders 0.87, 0.78, but towards 1.04 G: the finest
        # error exceeds the band 1.25 |G_16 - G_8| / (2 - 1) = 0.125 G.
        r = self.study(g=0.4, g_bias=0.04)
        self.assertFailsWith(r, "GCI G: finest error")
        self.assertEqual(len(r["failures"]), 1, r["failures"])

    def test_low_order_fails(self):
        self.assertFailsWith(self.study(order=0.3, transverse_order=0.3), "observed order")

    def test_oscillating_pressure_gradient_fails(self):
        # |G error| halves each level, but G itself oscillates about the analytic value.
        self.assertFailsWith(self.study(g_signs=(1, -1, 1)), "GCI G undefined")

    def test_flow_rate_floor_fails(self):
        r = self.study(flow_scale=1.02)
        self.assertFailsWith(r, "GCI profile: finest RMS error")

    def test_stagnant_transverse_velocity_fails(self):
        self.assertFailsWith(self.study(transverse_order=0.0), "transverse_max does not decrease")

    # ---- entrance and outlet development
    def test_entrance_transient_in_window_fails(self):
        self.assertFailsWith(self.study(entrance=1.0), "window not fully developed")

    def test_outlet_transient_in_window_fails(self):
        self.assertFailsWith(self.study(outlet=1.2, outlet_amplitude=0.5), "window not fully developed")

    def test_window_inside_entrance_region_fails(self):
        self.assertFailsWith(self.study(window=(1.0, 6.0)), "a-priori entrance or outlet region")

    def test_development_lengths_are_measured(self):
        r = self.study(levels=(4, 8, 16), entrance=0.4)
        d = r["runs"]["16-1"]["development"]
        # |U - u_fd| e^{-x/0.4} / u_c <= 1e-3 once x >= 0.4 ln(1000 max|U - u_fd| / u_c).
        self.assertGreater(d["profile_1e-3"], 2.0)
        self.assertLessEqual(d["profile_1e-3"], 3.0)

    # ---- rank parity
    def test_rank_velocity_mismatch_fails(self):
        def edit(c, r, rows):
            if r == 4 and c == 8:
                rows[1234][4] += 1e-5 * U
        self.assertFailsWith(self.study(ranks=(1, 2, 4), edit=edit), "parity cells 8: 4 vs 1 ranks")

    def test_rank_pressure_offset_fails(self):
        def edit(c, r, rows):
            if r == 2:
                for row in rows:
                    row[7] += 1e-5 * RHO * U * U
        self.assertFailsWith(self.study(ranks=(1, 2), edit=edit), "pressure/(rho U^2)")

    # ---- malformed or unconverged output
    def test_missing_node_fails(self):
        self.assertFailsWith(self.study(edit=lambda c, r, rows: rows.pop(17) if c == 8 else None), "lattice nodes missing")

    def test_duplicate_node_fails(self):
        def edit(c, r, rows):
            if c == 8:
                rows[18][0] = rows[17][0]
        self.assertFailsWith(self.study(edit=edit), "duplicate node")

    def test_nonfinite_value_fails(self):
        def edit(c, r, rows):
            if c == 4:
                rows[3][6] = float("nan")
        self.assertFailsWith(self.study(edit=edit), "nonfinite")

    def test_other_mesh_fails(self):
        def edit(c, r, rows):
            if c == 16:
                rows[40][2] += 1e-6
        self.assertFailsWith(self.study(edit=edit), "coordinates differ from the duct lattice")

    def test_unconverged_metrics_fail(self):
        def metrics(c, r, m):
            if c == 16:
                m["momentum"] = 3e-7
        self.assertFailsWith(self.study(metrics_edit=metrics), "not converged: final momentum")

    def test_unconverged_log_fails(self):
        self.assertFailsWith(self.study(log="NOT CONVERGED: iteration limit iterations=30000"), "no CONVERGED line")

    def test_outlet_reversal_fails(self):
        def metrics(c, r, m):
            m["closed_faces"] = 3
        self.assertFailsWith(self.study(metrics_edit=metrics), "outlet reversal")

    def test_other_viscosity_fails(self):
        # Same profile shape, but G scales with mu: a run at mu=0.2 must not be compared at mu=0.1.
        self.assertFailsWith(self.study(header="rho=1 mu=0.2 nu=0.2 inlet_speed=0.1"), "controls: the run used mu=0.2")

    def test_other_inflow_fails(self):
        def metrics(c, r, m):
            m["inlet_kg_s"] *= 1.01
        self.assertFailsWith(self.study(metrics_edit=metrics), "inlet mass flux")

    def test_two_levels_fail(self):
        self.assertFailsWith(self.study(levels=(8, 16)), "three nested levels")

    def test_exodus_hash_mismatch_fails(self):
        synthetic(self.dir, (4, 8, 16), (1,))
        with open(os.path.join(self.dir, "duct-8.exo"), "wb") as f:
            f.write(b"not the mesh")
        r = dc.study(self.dir, [4, 8, 16], [1], dc.Controls(RHO, MU, U))
        self.assertFailsWith(r, "SHA-256")

    # ---- command line
    def test_command_line_exit_codes(self):
        synthetic(self.dir, (4, 8, 16), (1, 2))
        self.assertEqual(dc.main(["study", self.dir, "--levels", "4,8,16", "--ranks", "1,2"]), 0)
        self.assertEqual(dc.main(["run", os.path.join(self.dir, "duct-8-2"), "--mesh", os.path.join(self.dir, "duct-8.json")]), 0)
        self.assertEqual(dc.main(["run", os.path.join(self.dir, "duct-8-2"), "--mesh", os.path.join(self.dir, "duct-4.json")]), 1)
        self.assertEqual(dc.main(["study", self.dir, "--levels", "4,8,16", "--ranks", "1,2", "--rank-tol", "-1"]), 1)


if __name__ == "__main__":
    unittest.main()
