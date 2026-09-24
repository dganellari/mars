# Does A_gram = D M^-1 D^T fit inside K's sparsity pattern, or is it wider?
# Decides whether the corrected operator can be a values-only blend of two existing matrices.
import numpy as np, itertools
LR    = [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]
OTHER = [(2,3),(0,3),(1,3),(1,2),(0,2),(0,1)]

def quad_area_normal(p0,p1,p2,p3):
    return 0.5*(np.cross(p1-p0,p2-p0) + np.cross(p2-p0,p3-p0))

def elem_ops(X):
    J = np.array([X[1]-X[0], X[2]-X[0], X[3]-X[0]]).T
    vol = np.linalg.det(J)/6.0
    Jinv = np.linalg.inv(J)
    g = np.zeros((4,3)); g[1:,:] = Jinv; g[0,:] = -g[1:,:].sum(axis=0)
    C = X.mean(axis=0); A = np.zeros((6,3))
    for ip,(L,R) in enumerate(LR):
        a,b = OTHER[ip]
        M  = 0.5*(X[L]+X[R]); Fa=(X[L]+X[R]+X[a])/3.0; Fb=(X[L]+X[R]+X[b])/3.0
        v = quad_area_normal(M,Fa,C,Fb)
        if v @ (X[R]-X[L]) < 0: v = -v
        A[ip] = v
    return vol, g, A

# Two tets sharing the face (0,1,2): nodes 0..4
P = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0.4,0.4,-0.8]], dtype=float)
tets = [(0,1,2,3), (0,2,1,4)]
n = len(P)
K = np.zeros((n,n)); Dg = np.zeros((n,n,3)); Vd = np.zeros(n)
for t in tets:
    X = P[list(t)]
    vol,g,A = elem_ops(X)
    if vol < 0:  # keep positive orientation
        t = (t[0],t[2],t[1],t[3]); X = P[list(t)]; vol,g,A = elem_ops(X)
    for i in range(4):
        Vd[t[i]] += vol/4.0
        for j in range(4):
            K[t[i],t[j]] += vol*(g[i]@g[j])
    for ip,(L,R) in enumerate(LR):
        for loc,s in ((L,+1.0),(R,-1.0)):
            Dg[t[loc], t[L]] += s*0.5*A[ip]
            Dg[t[loc], t[R]] += s*0.5*A[ip]

Ag = np.einsum('ijd,kjd->ik', Dg, Dg/Vd[None,:,None])
tol = 1e-12
pK  = np.abs(K)  > tol
pAg = np.abs(Ag) > tol
print("K   nnz:", pK.sum(), " A_gram nnz:", pAg.sum())
print("A_gram entries OUTSIDE K's pattern:", int((pAg & ~pK).sum()))
print("K entries outside A_gram's pattern:", int((pK & ~pAg).sum()))
print("\nK pattern:\n", pK.astype(int))
print("A_gram pattern:\n", pAg.astype(int))
