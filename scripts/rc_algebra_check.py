# Decide, with numbers, whether K - A_gram is a genuinely different operator or another
# scalar multiple of K. Replicates MARS's median-dual SCS construction exactly
# (mars_cvfem_tet_area.hpp: LR pairs {0,1,1,2,0,2,0,3,1,3,2,3}, other-nodes = the complement,
# dual quad M -> Fa -> C -> Fb, oriented L->R).
import numpy as np
np.set_printoptions(precision=4, suppress=True, linewidth=140)

LR    = [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]
OTHER = [(2,3),(0,3),(1,3),(1,2),(0,2),(0,1)]

def quad_area_normal(p0,p1,p2,p3):
    a1,b1 = p1-p0, p2-p0
    a2,b2 = p2-p0, p3-p0
    return 0.5*(np.cross(a1,b1) + np.cross(a2,b2))

def tet_geometry(X):
    # P1 shape-function gradients and volume
    J = np.array([X[1]-X[0], X[2]-X[0], X[3]-X[0]]).T
    det = np.linalg.det(J); vol = det/6.0
    Jinv = np.linalg.inv(J)
    g = np.zeros((4,3))
    g[1:,:] = Jinv                        # grad lambda_i = row i of J^-1
    g[0,:]  = -g[1:,:].sum(axis=0)        # partition of unity
    # SCS area vectors
    C = X.mean(axis=0)
    A = np.zeros((6,3))
    for ip,(L,R) in enumerate(LR):
        a,b = OTHER[ip]
        M  = 0.5*(X[L]+X[R])
        Fa = (X[L]+X[R]+X[a])/3.0
        Fb = (X[L]+X[R]+X[b])/3.0
        v = quad_area_normal(M,Fa,C,Fb)
        if v @ (X[R]-X[L]) < 0: v = -v
        A[ip] = v
    return vol, g, A

def operators(X):
    vol, g, A = tet_geometry(X)
    # (1) the identity the review claims: B_i = -vol * grad N_i
    B = np.zeros((4,3))
    for ip,(L,R) in enumerate(LR):
        B[L] += A[ip]; B[R] -= A[ip]
    # (2) Galerkin stiffness
    K = vol * (g @ g.T)
    # (3) CVFEM divergence D: div_i(u) = sum_ip s_i * 0.5*(u_L+u_R) . A_ip
    D = np.zeros((4,4,3))
    for ip,(L,R) in enumerate(LR):
        for i,s in ((L,+1.0),(R,-1.0)):
            D[i,L] += s*0.5*A[ip]
            D[i,R] += s*0.5*A[ip]
    Vd = np.full(4, vol/4.0)               # median-dual volume per node
    Agram = np.einsum('ijd,kjd->ik', D, D / Vd[None,:,None])
    return vol, g, B, K, Agram

print("=== identity check: B_i vs -vol*gradN_i ===")
rng = np.random.default_rng(7)
for t in range(3):
    X = rng.normal(size=(4,3))
    if np.linalg.det(np.array([X[1]-X[0],X[2]-X[0],X[3]-X[0]]).T) < 0: X[[2,3]] = X[[3,2]]
    vol,g,B,K,Ag = operators(X)
    ratio = B / (-vol*g)
    print(f" tet{t}: max|B/(-vol*gradN) - 1| = {np.abs(ratio-1).max():.2e}")

print("\n=== is K - A_gram a scalar multiple of K? ===")
for t in range(3):
    X = rng.normal(size=(4,3))
    if np.linalg.det(np.array([X[1]-X[0],X[2]-X[0],X[3]-X[0]]).T) < 0: X[[2,3]] = X[[3,2]]
    vol,g,B,K,Ag = operators(X)
    Kn = K/np.abs(K).max()
    An = Ag/np.abs(K).max()
    Dm = Kn - An
    # best scalar c minimising |Dm - c*Kn|
    c = (Dm*Kn).sum()/ (Kn*Kn).sum()
    resid = np.abs(Dm - c*Kn).max() / max(np.abs(Dm).max(), 1e-300)
    print(f" tet{t}: A_gram/K best-fit scale = {(Ag*K).sum()/(K*K).sum():+.4f}"
          f"   |K-Agram| / |K| = {np.abs(Dm).max():.3e}"
          f"   residual after removing best c*K = {resid:.3e}")

print("\n=== null space + symmetry of the difference ===")
X = rng.normal(size=(4,3))
if np.linalg.det(np.array([X[1]-X[0],X[2]-X[0],X[3]-X[0]]).T) < 0: X[[2,3]] = X[[3,2]]
vol,g,B,K,Ag = operators(X)
Dm = K - Ag
print(" row sums K-Agram :", Dm.sum(axis=1))
print(" asymmetry        :", np.abs(Dm-Dm.T).max()/max(np.abs(Dm).max(),1e-300))
print(" eigenvalues      :", np.linalg.eigvalsh((Dm+Dm.T)/2))
