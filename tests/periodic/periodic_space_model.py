"""Host model of the multi-rank periodic storage of MARS and of PeriodicNavierStokes.

Every simulated rank holds node SLOTS: owned nodes plus ghosts. A periodic
slave (on a max face) and its master (on the min face) are two slots of one
physical point. The model follows the MARS data layout:
  - elements are Morton-partitioned; each rank scatters only its owned
    elements, and the reverse halo adds ghost slots into their owners
  - the node halo refreshes ghosts that are corners of owned elements
    (fullTopo=True refreshes every ghost)
  - partner[slot] is the local slot of the final master, or -1 if absent
  - a slave owned on rank A with a master owned on rank B moves data through
    an explicit pair exchange, like crossRankPeriodicPairSum/Broadcast

The CVFEM hex operators use the formulas of the CUDA kernels on a uniform
cube: divergence, its transpose, the skew/upwind advection flux and the
cvfem_hex_diffusion_lhs element matrix. PeriodicNS is the scheme of
backend/distributed/unstructured/fem/mars_periodic_ns.hpp step for step, so
its numbers are a reference for the GPU run on the same mesh.
"""
import numpy as np

LR = np.array([(0, 1), (1, 2), (2, 3), (0, 3), (4, 5), (5, 6), (6, 7), (4, 7),
               (0, 4), (1, 5), (2, 6), (3, 7)])
OFF = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
# hexDerivConst of mars_cvfem_hex_kernel.hpp
DERIV = np.array([
    [[-1.0, -0.5, -0.5], [1.0, -0.5, -0.5], [0.0, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    [[-0.5, 0.0, 0.0], [0.5, -1.0, -0.5], [0.5, 1.0, -0.5], [-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0]],
    [[0.0, -0.5, 0.0], [0.0, -0.5, 0.0], [1.0, 0.5, -0.5], [-1.0, 0.5, -0.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.5]],
    [[-0.5, -1.0, -0.5], [0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.5, 1.0, -0.5], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5]],
    [[0.0, 0.0, -0.5], [0.0, 0.0, -0.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, -0.5, 0.5], [1.0, -0.5, 0.5], [0.0, 0.5, 0.0], [0.0, 0.5, 0.0]],
    [[0.0, 0.0, 0.0], [0.0, 0.0, -0.5], [0.0, 0.0, -0.5], [0.0, 0.0, 0.0], [-0.5, 0.0, 0.0], [0.5, -1.0, 0.5], [0.5, 1.0, 0.5], [-0.5, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -0.5], [0.0, 0.0, -0.5], [0.0, -0.5, 0.0], [0.0, -0.5, 0.0], [1.0, 0.5, 0.5], [-1.0, 0.5, 0.5]],
    [[0.0, 0.0, -0.5], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -0.5], [-0.5, -1.0, 0.5], [0.5, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.5, 1.0, 0.5]],
    [[-0.5, -0.5, -1.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [-0.5, -0.5, 1.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.5, 0.0]],
    [[-0.5, 0.0, 0.0], [0.5, -0.5, -1.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0], [-0.5, 0.0, 0.0], [0.5, -0.5, 1.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.5, 0.5, -1.0], [-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.5, 0.5, 1.0], [-0.5, 0.0, 0.0]],
    [[0.0, -0.5, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.5, 0.5, -1.0], [0.0, -0.5, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.5, 0.5, 1.0]],
])


def morton(a, b, c):
    key = 0
    for bit in range(10):
        key |= ((a >> bit) & 1) << (3 * bit) | ((b >> bit) & 1) << (3 * bit + 1) | ((c >> bit) & 1) << (3 * bit + 2)
    return key


class Mesh:
    """Uniform N^3 hex cube on [0,1]^3; the max faces keep their own node slots."""

    def __init__(self, N):
        self.N, self.h = N, 1.0 / N
        n1 = N + 1
        g = np.arange(n1 ** 3)
        self.ijk = np.stack([g % n1, (g // n1) % n1, g // (n1 * n1)], axis=1)
        self.nNodes = n1 ** 3
        self.xyz = self.ijk * self.h
        w = self.ijk % N
        self.master = w[:, 0] + n1 * (w[:, 1] + n1 * w[:, 2])
        self.isSlave = self.master != g
        e = np.array([(a, b, c) for c in range(N) for b in range(N) for a in range(N)])
        self.eijk = e
        self.conn = np.stack([(e[:, 0] + o[0]) + n1 * ((e[:, 1] + o[1]) + n1 * (e[:, 2] + o[2])) for o in OFF], axis=1)
        self.nElem = len(e)
        h = self.h
        # sub-control face area vectors of a uniform hex, from L to R, |A| = h^2/4
        self.A = np.array([(OFF[r] - OFF[l]) * h * h / 4 for l, r in LR], dtype=float)
        self.Ke = self._element_K()
        self.Ve = h ** 3

    def _element_K(self):
        coords = OFF * self.h
        K = np.zeros((8, 8))
        for ip, (l, r) in enumerate(LR):
            J = np.einsum("nj,ni->ij", DERIV[ip], coords)
            dndx = np.einsum("ji,nj->ni", np.linalg.inv(J), DERIV[ip])
            diff = -(dndx @ self.A[ip])
            K[l, :] += diff
            K[r, :] -= diff
        return K


class Rank:
    pass


def partition(mesh, R, fullTopo=False):
    keys = np.array([morton(*t) for t in mesh.eijk])
    chunks = np.array_split(np.argsort(keys, kind="stable"), R)
    elemOwner = np.empty(mesh.nElem, dtype=int)
    for r, c in enumerate(chunks):
        elemOwner[c] = r
    N = mesh.N
    cl = np.minimum(mesh.ijk, N - 1)
    nodeOwner = elemOwner[cl[:, 0] + N * (cl[:, 1] + N * cl[:, 2])]
    masterSetOfElem = [set(mesh.master[mesh.conn[e]]) for e in range(mesh.nElem)]
    elemsOfNode = [[] for _ in range(mesh.nNodes)]
    for e in range(mesh.nElem):
        for g in mesh.conn[e]:
            elemsOfNode[g].append(e)
    ranks = []
    for r in range(R):
        rk = Rank()
        rk.r = r
        owned = np.sort(chunks[r])
        ownedMasters = set().union(*(masterSetOfElem[e] for e in owned))
        # periodic cornerstone halo: every element touching an owned element,
        # across periodic faces too, plus the element star of each owned node
        halo = set(e for e in range(mesh.nElem) if masterSetOfElem[e] & ownedMasters)
        for g in np.nonzero(nodeOwner == r)[0]:
            halo |= set(elemsOfNode[g])
        halo -= set(owned.tolist())
        rk.localElems = np.concatenate([owned, np.array(sorted(halo), dtype=int)])
        rk.nodes = np.unique(mesh.conn[rk.localElems])
        rk.g2l = np.full(mesh.nNodes, -1)
        rk.g2l[rk.nodes] = np.arange(len(rk.nodes))
        rk.n = len(rk.nodes)
        rk.own = nodeOwner[rk.nodes] == r
        rk.lconn = rk.g2l[mesh.conn[owned]]
        inOwnedElem = np.zeros(rk.n, bool)
        inOwnedElem[rk.lconn.ravel()] = True
        rk.topoGhosts = np.nonzero((inOwnedElem | fullTopo) & ~rk.own)[0]
        m = rk.g2l[mesh.master[rk.nodes]]
        rk.partner = np.where(mesh.isSlave[rk.nodes], m, -1)
        ranks.append(rk)
    for rk in ranks:
        gh = rk.topoGhosts
        rk.ghostOwner = nodeOwner[rk.nodes[gh]]
        rk.ghostOwnerSlot = np.array([ranks[q].g2l[g] for q, g in zip(rk.ghostOwner, rk.nodes[gh])], dtype=int)
        ownedSlave = rk.own & mesh.isSlave[rk.nodes]
        rk.missingMaster = int(np.sum(ownedSlave & (rk.partner < 0)))
        s = np.nonzero(ownedSlave & (rk.partner >= 0))[0]
        same = rk.own[rk.partner[s]]
        rk.sameS, rk.sameM = s[same], rk.partner[s[same]]
        rk.crossS = s[~same]
        gm = rk.nodes[rk.partner[rk.crossS]]
        rk.crossQ = nodeOwner[gm]
        rk.crossM = np.array([ranks[q].g2l[g] for q, g in zip(rk.crossQ, gm)], dtype=int)
        rk.reduced = rk.own & ~mesh.isSlave[rk.nodes]
    return ranks


class Dist:
    """Distributed node fields: one array per rank, shape (n,) or (n, 3)."""

    def __init__(self, mesh, ranks):
        self.mesh, self.ranks = mesh, ranks

    def zeros(self, ncomp=None):
        return [np.zeros(rk.n if ncomp is None else (rk.n, ncomp)) for rk in self.ranks]

    def from_global(self, f):
        return [np.array(f[rk.nodes], dtype=float) for rk in self.ranks]

    def to_global(self, F, ncomp=None):
        out = np.zeros(self.mesh.nNodes if ncomp is None else (self.mesh.nNodes, ncomp))
        for rk in self.ranks:
            out[rk.nodes[rk.reduced]] = F[rk.r][rk.reduced]
        return out

    # cstone node halo
    def exchange(self, F):
        for rk in self.ranks:
            for q in np.unique(rk.ghostOwner):
                sel = rk.ghostOwner == q
                F[rk.r][rk.topoGhosts[sel]] = F[q][rk.ghostOwnerSlot[sel]]

    def reverse_add(self, F):
        sent = [(rk, np.copy(F[rk.r][rk.topoGhosts])) for rk in self.ranks]
        for rk, vals in sent:
            for q in np.unique(rk.ghostOwner):
                sel = rk.ghostOwner == q
                np.add.at(F[q], rk.ghostOwnerSlot[sel], vals[sel])

    # the reduced periodic space: P and P^T
    def fold(self, F):
        for rk in self.ranks:
            np.add.at(F[rk.r], rk.sameM, F[rk.r][rk.sameS])
            F[rk.r][rk.sameS] = 0
        for rk in self.ranks:
            for q in np.unique(rk.crossQ):
                sel = rk.crossQ == q
                np.add.at(F[q], rk.crossM[sel], F[rk.r][rk.crossS[sel]])
            F[rk.r][rk.crossS] = 0

    def bcast(self, F):
        for rk in self.ranks:
            F[rk.r][rk.sameS] = F[rk.r][rk.sameM]
            for q in np.unique(rk.crossQ):
                sel = rk.crossQ == q
                F[rk.r][rk.crossS[sel]] = F[q][rk.crossM[sel]]

    def restrict(self, F):
        self.reverse_add(F)
        self.fold(F)

    def prolong(self, F):
        self.bcast(F)
        self.exchange(F)

    def rdot(self, F, G):
        return sum(np.sum(F[rk.r][rk.reduced] * G[rk.r][rk.reduced]) for rk in self.ranks)

    def rsum(self, F):
        return sum(np.sum(F[rk.r][rk.reduced]) for rk in self.ranks)

    def rmax(self, F):
        return max(np.max(np.abs(F[rk.r][rk.reduced])) for rk in self.ranks)

    # element scatters over OWNED elements
    def scatter_div(self, U):
        out = self.zeros()
        for rk in self.ranks:
            for ip, (l, r) in enumerate(LR):
                iL, iR = rk.lconn[:, l], rk.lconn[:, r]
                flow = 0.5 * (U[rk.r][iL] + U[rk.r][iR]) @ self.mesh.A[ip]
                np.add.at(out[rk.r], iL, flow)
                np.add.at(out[rk.r], iR, -flow)
        return out

    def scatter_divT(self, P):
        out = self.zeros(3)
        for rk in self.ranks:
            for ip, (l, r) in enumerate(LR):
                iL, iR = rk.lconn[:, l], rk.lconn[:, r]
                c = np.outer(0.5 * (P[rk.r][iL] - P[rk.r][iR]), self.mesh.A[ip])
                np.add.at(out[rk.r], iL, c)
                np.add.at(out[rk.r], iR, c)
        return out

    def scatter_K(self, Q):
        out = self.zeros()
        for rk in self.ranks:
            np.add.at(out[rk.r], rk.lconn.ravel(), (Q[rk.r][rk.lconn] @ self.mesh.Ke.T).ravel())
        return out

    def scatter_Kdiag(self):
        out = self.zeros()
        for rk in self.ranks:
            np.add.at(out[rk.r], rk.lconn.ravel(), np.tile(np.diag(self.mesh.Ke), len(rk.lconn)))
        return out

    def scatter_adv(self, U, skew):
        out = self.zeros(3)
        for rk in self.ranks:
            for ip, (l, r) in enumerate(LR):
                iL, iR = rk.lconn[:, l], rk.lconn[:, r]
                mdot = 0.5 * (U[rk.r][iL] + U[rk.r][iR]) @ self.mesh.A[ip]
                qL, qR = U[rk.r][iL], U[rk.r][iR]
                if skew:
                    np.add.at(out[rk.r], iL, -0.5 * mdot[:, None] * (2 * qL + qR))
                    np.add.at(out[rk.r], iR, 0.5 * mdot[:, None] * (qL + 2 * qR))
                else:
                    flux = mdot[:, None] * np.where((mdot > 0)[:, None], qL, qR)
                    np.add.at(out[rk.r], iL, -flux)
                    np.add.at(out[rk.r], iR, flux)
        return out

    def scatter_mass(self):
        out = self.zeros()
        for rk in self.ranks:
            np.add.at(out[rk.r], rk.lconn.ravel(), self.mesh.Ve / 8)
        return out

    def scatter_pdiag(self, mass):
        out = self.zeros()
        for rk in self.ranks:
            for ip, (l, r) in enumerate(LR):
                iL, iR = rk.lconn[:, l], rk.lconn[:, r]
                val = 0.25 * self.mesh.A[ip] @ self.mesh.A[ip] * (1 / mass[rk.r][iL] + 1 / mass[rk.r][iR])
                np.add.at(out[rk.r], iL, val)
                np.add.at(out[rk.r], iR, val)
        return out


def seam_mismatch(d, F):
    """max |F[slave] - F[master]| over slot pairs an owned element can read."""
    worst = 0.0
    for rk in d.ranks:
        live = rk.own.copy()
        live[rk.topoGhosts] = True
        s = np.nonzero((rk.partner >= 0) & live)[0]
        s = s[live[rk.partner[s]]]
        if len(s):
            worst = max(worst, np.max(np.abs(F[rk.r][s] - F[rk.r][rk.partner[s]])))
    return worst


class PeriodicNS:
    """mars_periodic_ns.hpp: every field in range(P), every operator P^T A P."""

    def __init__(self, d, nu, rho, dt, skew, tol=1e-10, maxIter=1000, legacyPressureCoef=False):
        self.d, self.nu, self.rho, self.dt, self.skew = d, nu, rho, dt, skew
        self.tol, self.maxIter = tol, maxIter
        self.legacyPressureCoef = legacyPressureCoef
        self.mass = d.scatter_mass()
        d.restrict(self.mass)
        d.prolong(self.mass)
        self.diagK = d.scatter_Kdiag()
        d.restrict(self.diagK)
        self.diagP = d.scatter_pdiag(self.mass)
        d.restrict(self.diagP)
        self.steps = 0
        self.its = (0, 0, 0, 0)

    def inv_mass(self, F):
        # ghost slots outside the node halo keep mass 0; they are never read
        inv = [np.divide(1.0, m, out=np.zeros_like(m), where=m != 0) for m in self.mass]
        return [f * w[:, None] if f.ndim == 2 else f * w for w, f in zip(inv, F)]

    def remove_mean(self, F):
        mean = self.d.rsum(F) / sum(int(np.sum(rk.reduced)) for rk in self.d.ranks)
        return [f - mean for f in F]

    def gradient(self, P):
        """M^-1 G p = -grad p at DOF slots; p prolonged."""
        g = self.d.scatter_divT(P)
        self.d.restrict(g)
        return self.inv_mass(g)

    def divergence(self, U):
        out = self.d.scatter_div(U)
        self.d.restrict(out)
        return out

    def pressure_op(self, X):
        self.d.prolong(X)
        g = self.gradient(X)
        self.d.prolong(g)
        return self.divergence(g)

    def visc_op(self, c):
        def op(Q):
            self.d.prolong(Q)
            out = self.d.scatter_K(Q)
            self.d.restrict(out)
            return [self.nu * o + c * m * q for o, m, q in zip(out, self.mass, Q)]
        return op

    def pcg(self, op, b, x, diag, bRef, zeroGuess):
        d = self.d
        r = [bi.copy() for bi in b] if zeroGuess else [bi - ai for bi, ai in zip(b, op([xi.copy() for xi in x]))]

        def jacobi(r):
            return [np.where(rk.reduced, ri / np.where(rk.reduced, di, 1.0), 0.0) for rk, ri, di in zip(d.ranks, r, diag)]

        z = jacobi(r)
        p = [zi.copy() for zi in z]
        rz, rr = d.rdot(r, z), d.rdot(r, r)
        stop = self.tol * max(np.sqrt(d.rdot(b, b)), bRef)
        it = 0
        while np.sqrt(rr) > stop:
            if it == self.maxIter:
                raise RuntimeError("CG did not converge")
            Ap = op([pi.copy() for pi in p])
            alpha = rz / d.rdot(p, Ap)
            x = [xi + alpha * pi for xi, pi in zip(x, p)]
            r = [ri - alpha * ai for ri, ai in zip(r, Ap)]
            z = jacobi(r)
            rrn, rzn = d.rdot(r, r), d.rdot(r, z)
            p = [zi + (rzn / rz) * pi for zi, pi in zip(z, p)]
            rr, rz, it = rrn, rzn, it + 1
        d.prolong(x)
        return x, it

    def project(self, U, dtEff):
        """u^{n+1} = u** + dtEff/rho M^-1 G phi with A phi = -(rho/dtEff) D u**."""
        d = self.d
        coef = self.rho / (self.dt if self.legacyPressureCoef else dtEff)
        b = self.remove_mean([-coef * x for x in self.divergence(U)])
        flux = [np.linalg.norm(u, axis=1) * np.cbrt(m * m) for u, m in zip(U, self.mass)]
        bRef = coef * np.sqrt(d.rdot(flux, flux))
        phi, its = self.pcg(self.pressure_op, b, d.zeros(), self.diagP, bRef, True)
        phi = self.remove_mean(phi)
        g = self.gradient(phi)
        Unew = [u + (dtEff / self.rho) * gi for u, gi in zip(U, g)]
        d.prolong(Unew)
        return Unew, phi, its

    def start(self, U, P):
        self.d.prolong(U)
        P = self.remove_mean(P)
        self.d.prolong(P)
        self.U, self.P, self.Um1, self.advm1 = U, P, None, None
        self.steps = 0

    def step(self):
        d, dt = self.d, self.dt
        bdf2 = self.steps > 0
        dtEff = 2 * dt / 3 if bdf2 else dt
        g = self.gradient(self.P)
        adv = d.scatter_adv(self.U, self.skew)
        d.restrict(adv)
        Us = []
        for r_, rk in enumerate(d.ranks):
            invM = np.divide(1.0, self.mass[r_], out=np.zeros_like(self.mass[r_]), where=self.mass[r_] != 0)[:, None]
            if bdf2:
                us = (4 * self.U[r_] - self.Um1[r_]) / 3 + \
                    (2 * dt / 3) * ((2 * adv[r_] - self.advm1[r_]) * invM + g[r_] / self.rho)
            else:
                us = self.U[r_] + dt * (adv[r_] * invM + g[r_] / self.rho)
            Us.append(us)
        d.prolong(Us)
        c = 1.0 / dtEff
        diagV = [np.where(rk.reduced, c * m + self.nu * dk, 1.0) for rk, m, dk in zip(d.ranks, self.mass, self.diagK)]
        its = []
        comps = []
        for k in range(3):
            q0 = [u[:, k].copy() for u in Us]
            rhs = [c * m * q for m, q in zip(self.mass, q0)]
            q, it = self.pcg(self.visc_op(c), rhs, q0, diagV, 0.0, False)
            comps.append(q)
            its.append(it)
        Uss = [np.stack([comps[k][r_] for k in range(3)], axis=1) for r_ in range(len(d.ranks))]
        Unew, phi, itP = self.project(Uss, dtEff)
        P = self.remove_mean([p + f for p, f in zip(self.P, phi)])
        d.prolong(P)
        self.Um1, self.advm1 = self.U, adv
        self.U, self.P = Unew, P
        self.steps += 1
        self.its = (*its, itP)

    def kinetic_energy(self):
        return 0.5 * sum(np.sum(self.mass[rk.r][rk.reduced] * np.sum(self.U[rk.r][rk.reduced] ** 2, axis=1))
                         for rk in self.d.ranks)

    def max_divergence(self):
        return self.d.rmax(self.inv_mass(self.divergence(self.U)))


def tgv_initial(mesh, V0=1.0, rho=1.0):
    k = 2 * np.pi
    x, y, z = (mesh.xyz[:, i] * k for i in range(3))
    U = np.stack([V0 * np.sin(x) * np.cos(y) * np.cos(z),
                  -V0 * np.cos(x) * np.sin(y) * np.cos(z),
                  np.zeros_like(x)], axis=1)
    P = rho * V0 * V0 / 16 * (np.cos(2 * x) + np.cos(2 * y)) * (np.cos(2 * z) + 2)
    return U, P


def smooth_field(mesh, seed):
    rng = np.random.default_rng(seed)
    k = 2 * np.pi
    x, y, z = (mesh.xyz[:, i] * k for i in range(3))
    f = np.zeros(mesh.nNodes)
    for _ in range(4):
        a, b, c = rng.integers(0, 3, 3)
        f += rng.normal() * np.cos(a * x + b * y + c * z + rng.uniform(0, 2 * np.pi))
    return f
