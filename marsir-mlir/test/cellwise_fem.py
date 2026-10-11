# Building blocks of the element-local (cell-wise) solver references: the p = 7
# Knaus CVFEM basis and operator, the metric, the closed-form diagonal, and one
# structured block of hexahedra on the unit cube per Level, with its DSS cascade.
import numpy as np
from numpy.polynomial import legendre as L

p = 7; n = p + 1; P = p; nn = n * n; n3 = nn * n
c = np.zeros(p + 1); c[p] = 1
zeta = np.sort(np.concatenate(([-1.0], L.legroots(L.legder(c)), [1.0])))
xi, _ = L.leggauss(P)

def lag(nodes, x):
    m = len(nodes); v = np.ones(m); d = np.zeros(m)
    for j in range(m):
        others = [k for k in range(m) if k != j]
        den = np.prod([nodes[j] - nodes[k] for k in others])
        v[j] = np.prod([x - nodes[k] for k in others]) / den
        d[j] = sum(np.prod([x - nodes[k] for k in others if k != i]) for i in others) / den
    return v, d

Btil = np.array([lag(zeta, x)[0] for x in xi]); Dtil = np.array([lag(zeta, x)[1] for x in xi])
D = np.array([lag(zeta, x)[1] for x in zeta])
_pad = np.concatenate(([-1.0], xi, [1.0]))
_Winv = np.zeros((n, n))
for i in range(n):
    _, d = lag(_pad, zeta[i]); _Winv[i] = -np.cumsum(d[:n])
W = np.linalg.inv(_Winv)

def apply(u, G):
    """y = A u per element, u and y (E, n^3); G (E, 3, P, 3, n, n) = [dir][face][g0 g1 g2][s][r]."""
    u = u.reshape(-1, n, n, n); y = np.zeros_like(u)
    for d in range(3):
        U = np.moveaxis(u, 1 + d, 1).reshape(-1, n, nn)
        Y = np.moveaxis(y, 1 + d, 1)
        interp = (Btil @ U).reshape(-1, P, n, n)
        deriv = (Dtil @ U).reshape(-1, P, n, n)
        dt2 = interp @ D.T            # derivative along r
        dt1 = D @ interp              # derivative along s
        g = G[:, d]
        flux = g[:, :, 2] * deriv + g[:, :, 0] * dt2 + g[:, :, 1] * dt1
        intf = W @ (flux @ W.T)
        Y[:, :P] -= intf; Y[:, 1:] += intf
    return y.reshape(-1, n3)

t1 = (1, 0, 0); t2 = (2, 2, 1)
corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])

def shape(rf):
    f = 1 + corner_ref * rf[..., None, :]
    dN = np.empty(f.shape)
    for k in range(3):
        g = f.copy(); g[..., k] = corner_ref[:, k]
        dN[..., k] = 0.125 * g.prod(-1)
    return 0.125 * f.prod(-1), dN

def metric(C):
    G = np.empty((len(C), 3, P, 3, n, n))
    for d in range(3):
        rf = np.zeros((P, n, n, 3))
        rf[..., d] = xi[:, None, None]; rf[..., t1[d]] = zeta[None, :, None]
        rf[..., t2[d]] = zeta[None, None, :]
        J = np.einsum("eca,lsrcb->elsrab", C, shape(rf)[1])
        Ji = np.linalg.inv(J)
        gv = np.linalg.det(J)[..., None] * np.einsum("...ak,...k->...a", Ji, Ji[..., d, :])
        G[:, d, :, 0] = gv[..., t2[d]]; G[:, d, :, 1] = gv[..., t1[d]]; G[:, d, :, 2] = gv[..., d]
    return G

def diagonal(G):
    dg = np.zeros((len(G), n3))
    for node in range(n3):
        a, b, cc = node // nn, (node // n) % n, node % n
        for d in range(3):
            A = (a, b, cc)[d]; B = b if d == 0 else a; C = b if d == 2 else cc
            for l, sign in ((A - 1, 1.0), (A, -1.0)):
                if not 0 <= l < P:
                    continue
                g0, g1, g2 = G[:, d, l, 0], G[:, d, l, 1], G[:, d, l, 2]
                term = W[B, B] * W[C, C] * g2[:, B, C] * Dtil[l, A]
                term = term + W[B, B] * Btil[l, A] * (g0[:, B, :] @ (W[C] * D[:, C]))
                term = term + W[C, C] * Btil[l, A] * (g1[:, :, C] @ (W[B] * D[:, B]))
                dg[:, node] += sign * term
    return dg

def deform_vertices(ne, deform):
    lat = np.arange(ne + 1) / ne
    V = np.stack(np.meshgrid(lat, lat, lat, indexing="ij"), -1)
    return V + deform * np.prod(np.sin(np.pi * V), -1, keepdims=True)

class Level:
    """Structured block of ne^3 hexahedra on the unit cube, element e = (ex ne + ey) ne + ez."""
    def __init__(self, ne, deform):
        self.ne = ne; self.E = ne ** 3
        V = deform_vertices(ne, deform)
        eidx = np.indices((ne, ne, ne)).reshape(3, -1).T
        self.corners = np.stack([V[tuple((eidx + (s > 0)).T)] for s in corner_ref], 1)
        self.G = metric(self.corners)
        la, lb, lc = np.indices((n, n, n)).reshape(3, -1)
        def outer(e, l):
            return ((e[:, None] == 0) & (l[None, :] == 0)) | ((e[:, None] == ne - 1) & (l[None, :] == p))
        bnd = outer(eidx[:, 0], la) | outer(eidx[:, 1], lb) | outer(eidx[:, 2], lc)
        self.mask = (~bnd).astype(float)
        self.mult = self.dss(np.ones((self.E, n3)))
        self.diag = self.dss(diagonal(self.G))
        rf_nodes = np.stack([zeta[la], zeta[lb], zeta[lc]], -1)
        self.x = np.einsum("jc,eca->eja", shape(rf_nodes)[0], self.corners)

    def dss(self, fc):
        ne = self.ne
        v = fc.reshape(ne, ne, ne, n, n, n).copy()
        for ax in range(3):
            lo = [slice(None)] * 6; hi = [slice(None)] * 6
            lo[ax] = slice(0, ne - 1); hi[ax] = slice(1, ne)
            lo[3 + ax] = n - 1; hi[3 + ax] = 0
            s = v[tuple(lo)] + v[tuple(hi)]
            v[tuple(lo)] = s; v[tuple(hi)] = s
        return v.reshape(self.E, n3)

    def A(self, u): return apply(u, self.G)
    def PJ(self, r): return self.mask * self.dss(r) / self.diag
    def wdot(self, u, v): return float(np.sum(u * v / self.mult))
