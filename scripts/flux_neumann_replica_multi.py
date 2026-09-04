#!/usr/bin/env python3
"""Multi-element replica for the flux-consistent opening rows.

Extends scripts/flux_neumann_replica.py past its 1-element artefact: a 3x3x3 node
grid (Kuhn-split into 48 tets) has a genuine INTERIOR node to pin, and opening
faces made of several triangles sharing nodes -- so it exercises face-mass
assembly, which one element cannot.

Validates:   K phi = -(rho/dt) [ D(u**) + S_target ]
             S_i   = sum_j (A_f/12)(1+delta_ij) Un_j     summed over opening faces
"""
import numpy as np
from collections import defaultdict

rho, dt = 1000.0, 1e-3
N = 3
idx = lambda i, j, k: (i * N + j) * N + k
P = np.array([[i/(N-1), j/(N-1), k/(N-1)]
              for i in range(N) for j in range(N) for k in range(N)])

KUHN = [((0,0,0),(1,0,0),(1,1,0),(1,1,1)), ((0,0,0),(1,0,0),(1,0,1),(1,1,1)),
        ((0,0,0),(0,1,0),(1,1,0),(1,1,1)), ((0,0,0),(0,1,0),(0,1,1),(1,1,1)),
        ((0,0,0),(0,0,1),(1,0,1),(1,1,1)), ((0,0,0),(0,0,1),(0,1,1),(1,1,1))]
tets = []
for ci in range(N-1):
    for cj in range(N-1):
        for ck in range(N-1):
            for t in KUHN:
                e = [idx(ci+a, cj+b, ck+c) for (a,b,c) in t]
                J = np.column_stack([P[e[1]]-P[e[0]], P[e[2]]-P[e[0]], P[e[3]]-P[e[0]]])
                if np.linalg.det(J) < 0: e[1], e[2] = e[2], e[1]
                tets.append(e)
tets = np.array(tets)
nn = len(P)

# --- assemble K, lumped mass, and per-element gradients ----------------------
K, M, G = np.zeros((nn, nn)), np.zeros(nn), []
for e in tets:
    J  = np.column_stack([P[e[1]]-P[e[0]], P[e[2]]-P[e[0]], P[e[3]]-P[e[0]]])
    V  = np.linalg.det(J)/6.0
    Ji = np.linalg.inv(J)
    g  = np.zeros((4,3)); g[1],g[2],g[3] = Ji[0,:],Ji[1,:],Ji[2,:]; g[0] = -(g[1]+g[2]+g[3])
    G.append((V, g))
    K[np.ix_(e, e)] += V * (g @ g.T)
    M[e] += V/4.0

# --- boundary faces: appear in exactly one tet ------------------------------
cnt, owner = defaultdict(int), {}
for ei, e in enumerate(tets):
    for omit in range(4):
        f = tuple(sorted(e[i] for i in range(4) if i != omit))
        cnt[f] += 1; owner[f] = (ei, e[omit])
bnd = [(f, owner[f][1]) for f, c in cnt.items() if c == 1]

def face_geom(f, opp):
    a, b, c = P[f[0]], P[f[1]], P[f[2]]
    nv = np.cross(b-a, c-a); A = 0.5*np.linalg.norm(nv); n = nv/np.linalg.norm(nv)
    if n @ (P[opp] - a) > 0: n = -n
    return A, n, (A/12.0)*(np.ones((3,3))+np.eye(3))

inflow  = [(f, *face_geom(f, o)) for f, o in bnd if all(abs(P[i][2]) < 1e-12 for i in f)]
outflow = [(f, *face_geom(f, o)) for f, o in bnd if all(abs(P[i][2]-1) < 1e-12 for i in f)]
A_in  = sum(A for _, A, _, _ in inflow)
A_out = sum(A for _, A, _, _ in outflow)
Un_in  = -0.5
Un_out = -(A_in*Un_in)/A_out                        # balanced: oint U.n == 0
print(f"{len(tets)} tets, {nn} nodes;  inflow {len(inflow)} tris (A={A_in:.3f}), "
      f"outflow {len(outflow)} tris (A={A_out:.3f})")

rng   = np.random.default_rng(0)
ustar = rng.normal(scale=0.1, size=(nn, 3))

D = np.zeros(nn)
for e, (V, g) in zip(tets, G):
    us = ustar[e].sum(axis=0)
    for i in range(4): D[e[i]] += -(V/4.0)*(g[i] @ us)
S = np.zeros(nn)
for grp, Un in ((inflow, Un_in), (outflow, Un_out)):
    for f, A, n, Mf in grp: S[list(f)] += Mf @ np.full(3, Un)
print(f"compatibility  sum_i S_i = {S.sum():+.3e}   (must be ~0)")

# --- solve with the INTERIOR node pinned ------------------------------------
interior = [i for i in range(nn) if np.all((P[i] > 1e-12) & (P[i] < 1-1e-12))]
piv = interior[0]
print(f"pinning interior node {piv} at {P[piv]}")

rhs = -(rho/dt) * (D + S)
Kp, bp = K.copy(), rhs.copy()
Kp[piv,:] = 0.0; Kp[piv,piv] = 1.0; bp[piv] = 0.0
phi = np.linalg.solve(Kp, bp)

Gphi = np.zeros((nn,3))
for e, (V, g) in zip(tets, G):
    Gphi[e] += (V/4.0) * (g.T @ phi[e])
unew = ustar - (dt/rho) * (Gphi / M[:,None])

for name, grp, Un, A in (("inflow", inflow, Un_in, A_in), ("outflow", outflow, Un_out, A_out)):
    net = sum(float(np.ones(3) @ (Mf @ np.array([unew[i] @ n for i in f])))
              for f, _, n, Mf in grp)
    tgt = A*Un
    print(f"  {name:8s} net={net:+.6f}  target={tgt:+.6f}  rel.err={abs(net-tgt)/abs(tgt):.2%}")
