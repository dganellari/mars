#!/usr/bin/env python3
"""SUPG host validation -- 1D convection-diffusion, the textbook stabilization test.

Validates the exact tau + streamline-diffusion algebra the GPU momentum kernel will use, BEFORE GPU.
Problem: -nu u'' + a u' = 0 on [0,1], u(0)=0, u(1)=1. Exact = (exp(Pe x)-1)/(exp(Pe)-1), Pe=a/nu,
a boundary layer at x=1. P1 Galerkin convection is central -> OSCILLATES once the element Peclet
Pe_e = a h/(2 nu) > 1. SUPG adds tau*(a.gradN_i)*(residual); on P1 the residual is a u' (viscous=0),
so the streamline term is tau*a^2*integral(N_i' N_j') = tau*a^2 * stiffness -> extra streamline
diffusion that kills the wiggles. We use the SAME tau the GPU uses:
  Sb = sum_i |a . gradN_i| (= 2a/h in 1D),  h_s = 2|a|/Sb,
  tau = 1 / sqrt( (2|a|/h_s)^2 + (4 nu/h_s^2)^2 ) = 1/sqrt(Sb^2 + (nu Sb^2/|a|^2)^2).
The wiggle signature is min(u) < 0 (the exact solution is monotone in [0,1]).
"""
import numpy as np


def solve_1d(n, a, nu, supg):
    h = 1.0 / n
    N = n + 1
    A = np.zeros((N, N)); b = np.zeros(N)
    # tau via the GPU formula (1D: |a|=a, each element has gradN = -+1/h so a.gradN_i = -+a/h)
    Sb = 2.0 * a / h                       # sum_i |a.gradN_i| over the 2 element nodes
    adv = Sb; diff = nu * Sb * Sb / (a * a)
    tau = 1.0 / np.sqrt(adv * adv + diff * diff) if supg else 0.0
    Kdiff = (nu / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])         # diffusion
    Cconv = (a / 2.0) * np.array([[-1.0, 1.0], [-1.0, 1.0]])        # central-Galerkin convection
    Ssupg = (tau * a * a / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])  # SUPG streamline diffusion (tau a^2 stiffness)
    for e in range(n):
        idx = [e, e + 1]
        A[np.ix_(idx, idx)] += Kdiff + Cconv + Ssupg
    A[0, :] = 0; A[0, 0] = 1; b[0] = 0.0       # u(0)=0
    A[-1, :] = 0; A[-1, -1] = 1; b[-1] = 1.0   # u(1)=1
    u = np.linalg.solve(A, b)
    x = np.linspace(0, 1, N)
    Pe = a / nu
    uex = (np.exp(np.clip(Pe * x, -700, 700)) - 1) / (np.exp(min(Pe, 700)) - 1)
    return x, u, uex, a * h / (2 * nu), tau


if __name__ == "__main__":
    a = 1.0
    print("1D convection-diffusion, boundary layer at x=1. Wiggle signature: min(u) < 0.\n")
    print(f"{'nu':>8} {'Pe_e':>7} | {'Galerkin err':>13} {'min(u)':>10} | {'SUPG err':>10} {'min(u)':>10} {'tau':>10}")
    for nu in [0.05, 0.02, 0.01, 0.005, 0.002]:
        n = 20
        _, ug, uex, Pee, _ = solve_1d(n, a, nu, supg=False)
        _, us, _, _, tau = solve_1d(n, a, nu, supg=True)
        eg = np.max(np.abs(ug - uex)); es = np.max(np.abs(us - uex))
        print(f"{nu:>8.3f} {Pee:>7.2f} | {eg:>13.3e} {ug.min():>10.3e} | {es:>10.3e} {us.min():>10.3e} {tau:>10.3e}")
    print("\nExpect: Galerkin min(u)<0 (oscillates) once Pe_e>1 and large error; SUPG min(u)>=~0 (smooth),"
          " bounded error. That confirms the tau + streamline term stabilize as designed.")
