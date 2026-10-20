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
#   right32 rot       the FP32 V-cycle with every transfer summed in reversed axis order
#                     on half the elements, as rotated frames do on the GPU: copies of a
#                     node then differ by rounding, the residual cannot see such jumps,
#                     and they stay in x as an error
#   right32 rot cont  the same with the top prolongation made continuous as in MARS
#                     (marsir::PmgPtxF32): written on the counting copies only, then
#                     DSS; the jump of the V-cycle output must be exactly 0
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


def reversed_order(fn, flip):
    """fn (prolong or restrict) summing in reversed axis order on the flagged elements:
    the same map, other rounding."""
    up = fn is ref.prolong
    def go(x, I):
        if x.dtype != np.float32:
            return fn(x, I)
        nf, nc = I.shape
        n_in, n_out = (nc, nf) if up else (nf, nc)
        v = x.reshape(-1, n_in, n_in, n_in)
        out = np.empty((v.shape[0], n_out ** 3), np.float32)
        out[~flip] = fn(v[~flip].reshape(-1, n_in ** 3), I)
        r = fn(np.ascontiguousarray(v[flip][:, ::-1, ::-1, ::-1]).reshape(-1, n_in ** 3),
               np.ascontiguousarray(I[::-1, ::-1]))
        out[flip] = r.reshape(-1, n_out, n_out, n_out)[:, ::-1, ::-1, ::-1].reshape(-1, n_out ** 3)
        return out
    return go


def counting(lev):
    """One copy per node (the copy in the lowest element along each axis)."""
    n, ne = lev.B.n, lev.ne
    eidx = np.indices((ne, ne, ne)).reshape(3, -1).T
    l = np.indices((n, n, n)).reshape(3, -1).T
    keep = np.ones((lev.E, n ** 3), bool)
    for ax in range(3):
        keep &= (l[None, :, ax] != 0) | (eidx[:, None, ax] == 0)
    return keep.astype(np.float32)


def continuous(m, cnt):
    """m.vcycle with the top prolongation written on the counting copies, then DSS
    (pre = 0, as marsir::PmgPtxF32)."""
    plain = m.vcycle
    def vcycle(b, l=0):
        if l != 0:
            return plain(b, l)
        f = m.levels[0]
        e = ref.prolong(m.vcycle(ref.restrict(f.mask * b, m.I[0]), 1), m.I[0])
        return m.smooth(0, b, m.post, f.dss(cnt * e))
    return vcycle


def jump(lev, x):
    """max |copy - mean of the node's copies| / max |x|."""
    return float(np.max(np.abs(x - lev.dss(x / lev.mult))) / np.max(np.abs(x)))


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
    flip = np.random.default_rng(7).random(f.E) < 0.5
    plain_prolong, plain_restrict = ref.prolong, ref.restrict
    rot_prolong, rot_restrict = reversed_order(plain_prolong, flip), reversed_order(plain_restrict, flip)
    m32c = copy.copy(m32)
    m32c.vcycle = continuous(m32c, counting(f))

    def rotated(vcycle):
        def go(r):
            ref.prolong, ref.restrict = rot_prolong, rot_restrict
            try:
                return vcycle(r.astype(np.float32)).astype(np.float64)
            finally:
                ref.prolong, ref.restrict = plain_prolong, plain_restrict
        return go
    P32rot, P32cont = rotated(m32.vcycle), rotated(m32c.vcycle)
    u = f.mask * np.prod(np.sin(np.pi * f.x), -1)
    b = f.A(u)
    P64 = lambda r: mg.vcycle(r)
    P32 = lambda r: m32.vcycle(r.astype(np.float32)).astype(np.float64)
    norm = lambda v: np.sqrt(f.wdot(v, v))
    Pb = norm(P64(b))
    print(f"ne={ne:3d} coarse={coarse:>5}  jump of one V-cycle output: rot {jump(f, P32rot(b)):.1e}, "
          f"rot cont {jump(f, P32cont(b)):.1e}", flush=True)
    for name, solve in (("left64", lambda: ref.bicgstab(f, P64, b)),
                        ("left32", lambda: ref.bicgstab(f, P32, b)),
                        ("right64", lambda: right_bicgstab(f, P64, b)),
                        ("right32", lambda: right_bicgstab(f, P32, b)),
                        ("right32 rot", lambda: right_bicgstab(f, P32rot, b, tol=1e-8)),
                        ("right32 rot cont", lambda: right_bicgstab(f, P32cont, b, tol=1e-8))):
        x, it = solve()
        err = norm(x - u) / norm(u)
        true_res = norm(P64(b - f.A(x))) / Pb
        print(f"ne={ne:3d} coarse={coarse:>5}  {name:16s} {it:3d} it   error {err:.1e}   "
              f"true ||P r||/||P b|| {true_res:.1e}", flush=True)
    print(f"  ({time.time() - t:.0f} s)", flush=True)


if __name__ == "__main__":
    for spec in sys.argv[1:]:
        ne, coarse = spec.split(":")
        run(int(ne), coarse)
