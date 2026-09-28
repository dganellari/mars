#!/usr/bin/env python3
"""Independent checks of duct_analytic.py.

The Dirichlet Poisson problem has exactly one solution, so a finite-difference residual of
mu lap(u) + G, zero wall values and the symmetries together verify the series; the flow-rate
quadrature then checks the closed-form K, and tabulated Fanning numbers (Shah & London 1978)
check both against the literature. With --cxx BINARY the C++ header (duct_analytic.hpp,
printed by `analytic_check --table`) must agree with this module to 1e-13.
"""
import cmath
import math
import os
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duct_analytic as da  # noqa: E402

CXX = None
SHAPES = [(2.0, 1.0), (1.0, 2.0), (1.0, 1.0), (3.0, 0.75)]


def gauss(n):
    """Gauss-Legendre nodes and weights on [-1, 1]."""
    xs, ws = [], []
    for i in range(n):
        t = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, t
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * t * p1 - (k - 1) * p0) / k
            dp = n * (t * p1 - p0) / (t * t - 1)
            step = p1 / dp
            t -= step
            if abs(step) < 1e-16:
                break
        xs.append(t)
        ws.append(2 / ((1 - t * t) * dp * dp))
    return xs, ws


def flow_rate(duct, gradient, panels=24, order=8):
    """Composite Gauss quadrature of u over the quarter section, graded towards the walls."""
    xs, ws = gauss(order)

    def edges(length):
        return [length * (1 - (1 - k / panels) ** 2) for k in range(panels + 1)]
    ey, ez = edges(duct.width / 2), edges(duct.height / 2)
    total = 0.0
    for p in range(panels):
        hy, cy = (ey[p + 1] - ey[p]) / 2, (ey[p + 1] + ey[p]) / 2
        for q in range(panels):
            hz, cz = (ez[q + 1] - ez[q]) / 2, (ez[q + 1] + ez[q]) / 2
            for xi, wi in zip(xs, ws):
                for xj, wj in zip(xs, ws):
                    total += wi * wj * hy * hz * duct.velocity(cy + hy * xi, cz + hz * xj, gradient)
    return 4 * total


class AnalyticDuct(unittest.TestCase):
    mu, U = 0.1, 0.1

    def each(self):
        for w, h in SHAPES:
            d = da.Duct(w, h, self.mu)
            yield d, d.pressure_gradient(self.U)

    def test_poisson_residual(self):
        for d, G in self.each():
            step = 2e-3 * d.a
            worst = 0.0
            for fy in (0.0, 0.3, -0.55, 0.8, 0.93):
                for fz in (0.0, -0.25, 0.5, 0.77, 0.9):
                    y, z = fy * d.width / 2, fz * d.height / 2

                    def second(f):
                        return (-f(-2) + 16 * f(-1) - 30 * f(0) + 16 * f(1) - f(2)) / (12 * step * step)
                    lap = (second(lambda k: d.velocity(y + k * step, z, G))
                           + second(lambda k: d.velocity(y, z + k * step, G)))
                    worst = max(worst, abs(self.mu * lap + G) / G)
            self.assertLess(worst, 1e-6, (d.width, d.height))

    def test_no_slip_limits(self):
        for d, G in self.each():
            uc, eps = d.centerline(G), 1e-9 * d.a
            for f in (-0.9, -0.4, 0.0, 0.35, 0.8):
                for y, z in ((d.width / 2 - eps, f * d.height / 2), (-d.width / 2 + eps, f * d.height / 2),
                             (f * d.width / 2, d.height / 2 - eps), (f * d.width / 2, -d.height / 2 + eps)):
                    self.assertLess(abs(d.velocity(y, z, G)) / uc, 1e-7)
            self.assertEqual(d.velocity(d.width / 2, 0.1, G), 0.0)

    def test_symmetry_and_rotation(self):
        for d, G in self.each():
            r = da.Duct(d.height, d.width, self.mu)
            for y in (0.1, 0.37):
                for z in (0.05, 0.33):
                    v = d.velocity(y * d.width, z * d.height, G)
                    self.assertEqual(v, d.velocity(-y * d.width, z * d.height, G))
                    self.assertEqual(v, d.velocity(y * d.width, -z * d.height, G))
                    self.assertEqual(v, r.velocity(z * d.height, y * d.width, G))

    def test_two_expansions_agree(self):
        for d, G in self.each():
            uc = d.centerline(G)
            for fs in (0.0, 0.3, 0.6, 0.85):
                for ft in (0.0, 0.4, 0.7, 0.9):
                    s, t = fs * d.a, ft * d.b
                    diff = abs(d.expansion(s, t, d.a, d.b, G) - d.expansion(t, s, d.b, d.a, G))
                    self.assertLess(diff / uc, 1e-13)

    def test_flow_rate_quadrature(self):
        for d, G in self.each():
            q = flow_rate(d, G)
            self.assertLess(abs(q - self.U * d.area) / (self.U * d.area), 1e-9, (d.width, d.height))

    def test_literature(self):
        # Shah & London (1978), Table 45 (Fanning f Re) and the square-duct u_max/u_mean.
        self.assertAlmostEqual(da.Duct(1, 1).fanning_re(), 14.22708, delta=5e-5)
        self.assertAlmostEqual(da.Duct(2, 1).fanning_re(), 15.54806, delta=5e-5)
        self.assertAlmostEqual(da.Duct(4, 1).fanning_re(), 18.23278, delta=5e-5)
        sq = da.Duct(1, 1, self.mu)
        self.assertAlmostEqual(sq.centerline(sq.pressure_gradient(self.U)) / self.U, 2.09624, delta=5e-5)
        plates = da.Duct(1000, 1, self.mu)
        self.assertLess(abs(plates.fanning_re() - 24), 0.05)

    def test_planar_parabola_is_not_the_duct(self):
        # The negative control the comparator must reject: same flow rate, wrong G and profile.
        d = da.Duct(2, 1, self.mu)
        G = d.pressure_gradient(self.U)
        Gp = da.planar_pressure_gradient(d.height, self.mu, self.U)
        self.assertAlmostEqual(Gp, 0.12, places=14)
        self.assertGreater((G - Gp) / G, 0.3)
        self.assertGreater(abs(da.planar_velocity(0, d.height, self.U) - d.centerline(G)) / self.U, 0.4)
        self.assertGreater(da.planar_velocity(0, d.height, self.U), 0)
        self.assertEqual(d.velocity(d.width / 2, 0, G), 0.0)   # the parabola is 1.5 U on this wall

    def test_entrance_estimates(self):
        z = complex(4.2, 2.25)   # first root of sin z + z = 0
        for _ in range(50):
            z -= (cmath.sin(z) + z) / (cmath.cos(z) + 1)
        self.assertAlmostEqual(z.real / 2, da.FADLE, places=6)
        self.assertAlmostEqual(da.durst_entrance_length(0, 1), 0.619, places=12)
        d = da.Duct(2, 1, self.mu)
        self.assertAlmostEqual(da.stokes_decay_length(d, 1e4), 0.5 * math.log(1e4) / da.FADLE, places=14)

    def test_cxx_agrees(self):
        if CXX is None:
            self.skipTest("C++ table not requested (--cxx BINARY)")
        out = subprocess.check_output([CXX, "--table"], universal_newlines=True).split("\n")
        rows = [r.split() for r in out if r.startswith("u ") or r.startswith("K ")]
        self.assertGreater(len(rows), 20)
        for r in rows:
            w, h = float(r[1]), float(r[2])
            d = da.Duct(w, h, self.mu)
            if r[0] == "K":
                self.assertLess(abs(float(r[3]) - d.shape_factor()), 1e-14)
                self.assertLess(abs(float(r[4]) - d.pressure_gradient(self.U)) / float(r[4]), 1e-14)
            else:
                G = d.pressure_gradient(self.U)
                self.assertLess(abs(float(r[5]) - d.velocity(float(r[3]), float(r[4]), G)) / d.centerline(G), 1e-13)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--cxx":
        CXX = sys.argv[2]
        del sys.argv[1:3]
    unittest.main()
