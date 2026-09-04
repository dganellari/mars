#!/usr/bin/env python3
"""One-element replica: flux-consistent opening rows for MARS's projection.

Requirement, in MARS's own discrete operators (mars_fem_projection.hpp):
  full weak divergence  D_full = D + S,   D_i(u) = -(V/4) gradN_i . sum_j u_j
                                          S_i(u) = oint N_i (u.n) dS
  corrector             u^{n+1} = u** - (dt/rho) M^-1 G phi
Demanding D_full(u^{n+1}) = 0 with the opening flux held at its TARGET gives
      D(M^-1 G phi) = (rho/dt) [ D(u**) + S_target ]
i.e. the operator is the Gram form D M^-1 G (for one element == -K exactly),
and the opening enters as a known target flux, NOT as a dphi/dn condition.

Two lessons this replica already caught:
  * the corrector supplies ONE constant gradient per element, so a per-node
    normal-velocity target is over-constrained -- the meaningful check is the
    NET flux through each opening face;
  * sum_i S_i = oint U.n must be ZERO or the pure-Neumann system is incompatible
    and the pin silently absorbs the imbalance (this is OpenFOAM's adjustPhi,
    already present in MARS as the opening-flux oScale rescaling).
"""
import numpy as np

rho, dt = 1000.0, 1e-3
p = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
J    = np.column_stack([p[1]-p[0], p[2]-p[0], p[3]-p[0]])
V    = np.linalg.det(J) / 6.0
Ji   = np.linalg.inv(J)
g    = np.zeros((4,3)); g[1],g[2],g[3] = Ji[0,:],Ji[1,:],Ji[2,:]; g[0] = -(g[1]+g[2]+g[3])
K    = V * (g @ g.T)
M    = np.full(4, V/4.0)

def face_data(ids, opposite):
    a, b, c = p[ids[0]], p[ids[1]], p[ids[2]]
    nv = np.cross(b-a, c-a); A = 0.5*np.linalg.norm(nv); nn = nv/np.linalg.norm(nv)
    if nn @ (p[opposite] - a) > 0: nn = -nn          # outward = away from the 4th node
    return ids, A, nn, (A/12.0)*(np.ones((3,3))+np.eye(3))

fIn  = face_data([0,1,2], 3)      # inflow face
fOut = face_data([0,1,3], 2)      # outflow face
Un_in  = -0.5                                     # into the domain
Un_out = -(fIn[1]*Un_in)/fOut[1]                  # balanced: sum of oint U.n == 0

rng = np.random.default_rng(0)
ustar = rng.normal(scale=0.1, size=(4,3))

D = np.array([-(V/4.0)*(g[i] @ ustar.sum(axis=0)) for i in range(4)])
S = np.zeros(4)
for (ids, A, nn, Mf), Un in ((fIn, Un_in), (fOut, Un_out)):
    S[ids] += Mf @ np.full(3, Un)
print(f"compatibility  sum_i S_i = {S.sum():+.3e}   (must be ~0)")

for sign in (+1.0, -1.0):
    rhs = sign * (rho/dt) * (D + S)
    Kp, bp = K.copy(), rhs.copy()
    Kp[2,:] = 0.0; Kp[2,2] = 1.0; bp[2] = 0.0     # one null-space pin
    phi = np.linalg.solve(Kp, bp)
    unew = ustar - (dt/rho) * ((V/4.0)*np.tile(g.T @ phi, (4,1))/M[:,None])
    print(f"  sign {sign:+.0f}:", end="")
    for name, (ids, A, nn, Mf), Un in (("in", fIn, Un_in), ("out", fOut, Un_out)):
        got = float(np.ones(3) @ (Mf @ np.array([unew[i] @ nn for i in ids])))
        print(f"   {name}: net={got:+.5f} target={A*Un:+.5f}", end="")
    print()
