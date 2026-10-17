# Mixed precision for the cell-wise p-multigrid solve (cellwise_pmg_ref.py): can the
# V-cycle run in FP32 inside the FP64 Krylov solver?
#
# The left-preconditioned BiCGStab of the solver (cellwise_pmg_ref.bicgstab, as on the
# GPU) assumes a fixed linear preconditioner; an FP32 V-cycle is deterministic but not
# linear at the 1e-7 level, which acts like a varying preconditioner. The recursive
# residual such a solver stops on can then drift from the true one, so each run also
# reports the true algebraic error against the manufactured solution and the true
# preconditioned residual ||P(b - A x)|| / ||P b|| (both with the FP64 V-cycle).
#
# Variants:
#   left64   left-preconditioned BiCGStab, FP64 V-cycle (the solver today)
#   left32   the same with an FP32 V-cycle
#   right64  right-preconditioned BiCGStab (A M y = b, x = M y), FP64 V-cycle: the
#            residual is the true assembled b - A x, and M may vary between calls
#            (the form of flexible Krylov methods for a varying preconditioner); it
#            stops on that residual, a different (here stricter) test than left's
#   right32  the same with an FP32 V-cycle
#
# Run: python3 cellwise_mixed_ref.py ne:coarse ...   e.g. 4:exact 8:exact 8:100
import copy, sys, time
import numpy as np
import cellwise_pmg_ref as ref


def fp32_copy(mg):
    """The multigrid with every matrix, metric, diagonal and the coarse solve in FP32."""
    m = copy.deepcopy(mg)
    f32 = np.float32
    for b in m.B:
        for name in ("Btil", "Dtil", "D", "W", "zeta", "xi"):
            setattr(b, name, getattr(b, name).astype(f32))
    for lev in m.levels:
        for name in ("G", "mask", "mult", "diag"):
            setattr(lev, name, getattr(lev, name).astype(f32))
    m.I = [I.astype(f32) for I in m.I]
    if m.coarse is not None:
        c = m.coarse
        Ainv = c.Ainv.astype(f32)

        def coarse32(r):
            rhs = np.zeros(c.NV, f32); np.add.at(rhs, c.gid, r)
            x = np.zeros(c.NV, f32); x[c.idx] = Ainv @ rhs[c.idx]
            return x[c.gid]
        m.coarse = coarse32
    return m


def right_bicgstab(lev, M, b, tol=1e-10, maxit=300):
    """BiCGStab on A M y = b. Residuals are assembled (DSS, Dirichlet rows zero) and
    continuous, so the weighted dot product is the assembled one; M takes an
    unassembled vector (the V-cycle's input), so an assembled r is split evenly."""
    asm = lambda y: lev.mask * lev.dss(y)
    Mv = lambda r: M(r / lev.mult)
    dot = lev.wdot
    x = np.zeros_like(b); r = asm(b); rh = r.copy(); r0 = np.sqrt(dot(r, r))
    rho0 = al = om = 1.0; v = np.zeros_like(b); pp = np.zeros_like(b)
    for it in range(1, maxit + 1):
        rho = dot(rh, r); be = (rho / rho0) * (al / om); rho0 = rho
        pp = r + be * (pp - om * v); ph = Mv(pp); v = asm(lev.A(ph)); al = rho / dot(rh, v)
        s = r - al * v; sh = Mv(s); t = asm(lev.A(sh)); om = dot(t, s) / dot(t, t)
        x = x + al * ph + om * sh; r = s - om * t
        if np.sqrt(dot(r, r)) <= tol * r0: return x, it
    return x, maxit


def run(ne, coarse):
    t = time.time()
    mg = ref.PMultigrid(ne, 0.05, [7, 3, 1], 0, 3, coarse)
    m32 = fp32_copy(mg)
    f = mg.levels[0]
    u = f.mask * np.prod(np.sin(np.pi * f.x), -1)
    b = f.A(u)
    P64 = lambda r: mg.vcycle(r)
    P32 = lambda r: m32.vcycle(r.astype(np.float32)).astype(np.float64)
    norm = lambda v: np.sqrt(f.wdot(v, v))
    Pb = norm(P64(b))
    for name, solve in (("left64", lambda: ref.bicgstab(f, P64, b)),
                        ("left32", lambda: ref.bicgstab(f, P32, b)),
                        ("right64", lambda: right_bicgstab(f, P64, b)),
                        ("right32", lambda: right_bicgstab(f, P32, b))):
        x, it = solve()
        err = norm(x - u) / norm(u)
        true_res = norm(P64(b - f.A(x))) / Pb
        print(f"ne={ne:3d} coarse={coarse:>5}  {name:8s} {it:3d} it   error {err:.1e}   "
              f"true ||P r||/||P b|| {true_res:.1e}", flush=True)
    print(f"  ({time.time() - t:.0f} s)", flush=True)


if __name__ == "__main__":
    for spec in sys.argv[1:]:
        ne, coarse = spec.split(":")
        run(int(ne), coarse)
