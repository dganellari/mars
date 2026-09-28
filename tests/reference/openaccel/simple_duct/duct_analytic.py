"""Fully developed laminar flow in a rectangular duct (Boussinesq series).

Cross-section |y| <= W/2, |z| <= H/2, flow along +x, G = -dp/dx > 0:

    mu (u_yy + u_zz) = -G,   u = 0 on all four walls.

With a the shorter and b the longer half-side, s the coordinate across the short side, t the
one along the long side, and lambda_i = i pi / (2a) for odd i:

    u = G/(2 mu) (a^2 - s^2)
        - 16 a^2 G / (mu pi^3) sum_i (-1)^((i-1)/2) i^-3 cos(lambda_i s) cosh(lambda_i t) / cosh(lambda_i b)

    U_mean = a^2 G K / (3 mu),   K = 1 - 192 a / (pi^5 b) sum_i tanh(i pi b / (2a)) / i^5

The first term alone is the planar Poiseuille parabola; the series is what makes the solution
vanish on the short walls |t| = b. For a finite duct the parabola is wrong everywhere (for
W:H = 2:1 its pressure gradient is 31% low at the same flow rate), so it is never used here
except as a negative control.

Numerics: the series terms decay like exp(-lambda_i (b - |t|)); the same solution written with
the sides exchanged decays like exp(-i pi (a - |s|) / (2b)). Each point uses the faster one.
cosh ratios are evaluated as exp(lambda (|t| - b)) (1 + exp(-2 lambda |t|)) / (1 + exp(-2 lambda b)),
which cannot overflow. Summation stops when the remaining terms are below 1e-18 of the leading
scale, or after 40000 terms (tail < 2e-10; reached only within ~1e-4 a of a corner).

Python 3.6 compatible, standard library only (the Alps system python).
"""
import math

PI = math.pi


class Duct(object):
    def __init__(self, width=2.0, height=1.0, viscosity=0.1):
        if not all(math.isfinite(v) and v > 0 for v in (width, height, viscosity)):
            raise ValueError("duct width, height and viscosity must be finite and positive")
        self.width, self.height, self.viscosity = float(width), float(height), float(viscosity)

    @property
    def a(self):
        return 0.5 * min(self.width, self.height)

    @property
    def b(self):
        return 0.5 * max(self.width, self.height)

    @property
    def area(self):
        return self.width * self.height

    @property
    def hydraulic_diameter(self):
        return 2 * self.width * self.height / (self.width + self.height)

    def shape_factor(self):
        """K(b/a): 1 for parallel plates, 0.42173 for the square, 0.68605 for 2:1."""
        ratio = self.b / self.a
        total = 0.0
        for i in range(1, 200000, 2):
            term = math.tanh(i * PI * ratio / 2) / float(i) ** 5
            total += term
            if term < 1e-19 * total:
                break
        return 1 - 192 / (PI ** 5 * ratio) * total

    def pressure_gradient(self, mean_velocity):
        """G = -dp/dx that carries the mean velocity U."""
        return 3 * self.viscosity * mean_velocity / (self.a ** 2 * self.shape_factor())

    def mean_velocity(self, gradient):
        return self.a ** 2 * gradient * self.shape_factor() / (3 * self.viscosity)

    def fanning_re(self):
        """Fanning f Re_Dh = G D_h^2 / (2 mu U): 14.227 square, 15.548 for 2:1, 24 for plates."""
        d = self.hydraulic_diameter
        return 3 * d * d / (2 * self.a ** 2 * self.shape_factor())

    def expansion(self, across, along, half_across, half_along, gradient):
        """Plane profile across `half_across` plus the series decaying along `along`."""
        series = 0.0
        for i in range(1, 80000, 2):
            lam = i * PI / (2 * half_across)
            decay = math.exp(lam * (along - half_along))
            ratio = decay * (1 + math.exp(-2 * lam * along)) / (1 + math.exp(-2 * lam * half_along))
            sign = 1.0 if ((i - 1) // 2) % 2 == 0 else -1.0
            series += sign * math.cos(lam * across) * ratio / float(i) ** 3
            if 2 * decay / float(i) ** 3 < 1e-18:
                break
        return gradient / self.viscosity * (0.5 * (half_across ** 2 - across ** 2)
                                            - 16 * half_across ** 2 / PI ** 3 * series)

    def velocity(self, y, z, gradient):
        """u(y, z) for pressure gradient G (centred coordinates)."""
        z_short = self.height <= self.width
        s, t = (abs(z), abs(y)) if z_short else (abs(y), abs(z))
        a, b = self.a, self.b
        if s >= a or t >= b:
            return 0.0
        if (b - t) / a >= (a - s) / b:
            return self.expansion(s, t, a, b, gradient)
        return self.expansion(t, s, b, a, gradient)

    def centerline(self, gradient):
        return self.velocity(0.0, 0.0, gradient)


def planar_velocity(z, height, mean_velocity):
    """Plane Poiseuille parabola across the height (negative control only)."""
    return 1.5 * mean_velocity * max(0.0, 1 - (2 * z / height) ** 2)


def planar_pressure_gradient(height, viscosity, mean_velocity):
    return 12 * viscosity * mean_velocity / height ** 2


def durst_entrance_length(reynolds, diameter):
    """Durst et al. (J. Fluids Eng. 127, 2005) pipe correlation applied with D = D_h.

    L/D = [0.619^1.6 + (0.0567 Re)^1.6]^(1/1.6). An order-of-magnitude estimate for a
    rectangular duct; the comparator measures the development length of the discrete solution.
    """
    return diameter * (0.619 ** 1.6 + (0.0567 * reynolds) ** 1.6) ** (1 / 1.6)


# Plane-channel estimate used only to choose the window margin, not a duct decay bound.
# FADLE = Re(z1)/2, where z1 = 4.21239 + 2.25073i solves sin z + z = 0.
FADLE = 2.1061961


def stokes_decay_length(duct, factor):
    """Window margin from the plane-channel decay estimate and requested `factor`."""
    return duct.a * math.log(factor) / FADLE
