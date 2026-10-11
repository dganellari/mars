#!/usr/bin/env python3
"""Reference for the element-local DSS on UNSTRUCTURED hex meshes
(backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp; design:
docs/design/cellwise_unstructured_dss.md).

Tables come from global vertex ids only: face pairs with each side's canonical face
frame, edge stars with reversed bits, vertex stars. Every copy of a node gets the sum
of all copies taken in one canonical order (copies by global element), so all copies
hold identical bits and the result does not depend on how elements are split over
ranks. Two kernel shapes are specified and must agree bit for bit: a per-copy gather
(each copy sums its partners) and the entity-centric form the GPU uses (one sum per
face/edge/vertex node, written to every copy).

Checks, on a block whose elements sit in random proper rotations of their frames (all
8 face orientations, reversed edges) and on extruded meshes with edges shared by 3 and
by 5 elements:
  * every copy of a node holds the same bits, equal to a canonical-order assembly;
  * the entity-centric form equals the per-copy gather bit for bit;
  * splitting the elements over 2, 3 or 5 ranks, each reading only the remote copies a
    single-phase exchange delivers (copies of nodes it also holds), gives the same bits.

Run: python3 test/cellwise_unstructured_ref.py
"""
import itertools, numpy as np

p = 7; n = p + 1
# MARS corner order (c_hexCornerRef): signs per local corner
CR = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
def corner_of(sign):            # local corner index with this sign vector
    return int(np.where((CR == np.array(sign)).all(1))[0][0])
def lidx(a, b, c): return (a * n + b) * n + c

# --- the 24 proper rotations of the cube, as signed axis permutations ---
ROT = []
for perm in itertools.permutations(range(3)):
    for sg in itertools.product((-1, 1), repeat=3):
        M = np.zeros((3, 3), int)
        for i in range(3): M[i, perm[i]] = sg[i]
        if round(np.linalg.det(M)) == 1: ROT.append(M)
assert len(ROT) == 24

def rotate_element(conn_e, R):
    """Same hexahedron, local frame rotated: new local corner c sits where old corner R s_c was."""
    return [conn_e[corner_of(R @ CR[c])] for c in range(8)]

# --- local entities of a hex in node coordinates ---
def face_nodes(f):
    """Face f (axis f//2, side f%2): local (i, j) -> node index, (i, j) over the other two axes."""
    k, side = f // 2, f % 2
    o = [a for a in range(3) if a != k]
    idx = np.empty((n, n), int)
    for i in range(n):
        for j in range(n):
            x = [0, 0, 0]; x[k] = side * p; x[o[0]] = i; x[o[1]] = j
            idx[i, j] = lidx(*x)
    return idx
FACE_NODES = [face_nodes(f) for f in range(6)]
def face_corner_local(f):        # local corners at face positions (0,0),(p,0),(0,p),(p,p)
    k, side = f // 2, f % 2
    o = [a for a in range(3) if a != k]
    out = []
    for (i, j) in [(0, 0), (p, 0), (0, p), (p, p)]:
        s = [0, 0, 0]; s[k] = 2 * side - 1; s[o[0]] = -1 if i == 0 else 1; s[o[1]] = -1 if j == 0 else 1
        out.append(corner_of(s))
    return out
FACE_CORNERS = [face_corner_local(f) for f in range(6)]
EDGES = []                         # (axis, fixed other coords) -> node indices along the edge, end corners
for k in range(3):
    o = [a for a in range(3) if a != k]
    for u in (0, p):
        for v in (0, p):
            nodes = []
            for t in range(n):
                x = [0, 0, 0]; x[k] = t; x[o[0]] = u; x[o[1]] = v
                nodes.append(lidx(*x))
            s0 = [0, 0, 0]; s0[k] = -1; s0[o[0]] = -1 if u == 0 else 1; s0[o[1]] = -1 if v == 0 else 1
            s1 = list(s0); s1[k] = 1
            EDGES.append((np.array(nodes), corner_of(s0), corner_of(s1)))
CORNER_NODE = [lidx(*[0 if s < 0 else p for s in CR[c]]) for c in range(8)]

def canonical_face_map(vids):
    """vids: global vertex ids at face positions (0,0),(p,0),(0,p),(p,p). Returns canon[i,j]:
    canonical position I*n+J. Origin at the smallest id; first axis toward its smaller neighbour."""
    pos = [(0, 0), (p, 0), (0, p), (p, p)]
    o = int(np.argmin(vids))
    oi, oj = pos[o]
    ni = pos.index((p - oi, oj)); nj = pos.index((oi, p - oj))     # neighbours along i and along j
    first_is_i = vids[ni] < vids[nj]
    canon = np.empty((n, n), int)
    for i in range(n):
        for j in range(n):
            di = i if oi == 0 else p - i      # distance from the origin along i
            dj = j if oj == 0 else p - j
            I, J = (di, dj) if first_is_i else (dj, di)
            canon[i, j] = I * n + J
    return canon

class Tables:
    """Per-element entity tables of a conforming hex mesh, canonical order = global element id."""
    def __init__(self, conn):
        self.conn = conn = [list(map(int, c)) for c in conn]
        E = len(conn)
        faces = {}                    # sorted vertex ids -> list of (e, f, canon map)
        for e in range(E):
            for f in range(6):
                v = [conn[e][c] for c in FACE_CORNERS[f]]
                faces.setdefault(tuple(sorted(v)), []).append((e, f, canonical_face_map(v)))
        self.partner = -np.ones((E, 6), int)       # neighbour element per face
        self.pmap = np.zeros((E, 6, n, n), int)    # own face (i, j) -> neighbour's node index
        for key, cps in faces.items():
            assert len(cps) <= 2, "non-conforming face"
            if len(cps) == 2:
                (e0, f0, c0), (e1, f1, c1) = cps
                inv1 = {c1[i, j]: FACE_NODES[f1][i, j] for i in range(n) for j in range(n)}
                inv0 = {c0[i, j]: FACE_NODES[f0][i, j] for i in range(n) for j in range(n)}
                self.partner[e0, f0], self.partner[e1, f1] = e1, e0
                self.pmap[e0, f0] = np.vectorize(inv1.get)(c0)
                self.pmap[e1, f1] = np.vectorize(inv0.get)(c1)
        edges = {}                    # (min id, max id) -> [(e, k, reversed)]
        for e in range(E):
            for k, (nodes, c0, c1) in enumerate(EDGES):
                a, b = conn[e][c0], conn[e][c1]
                edges.setdefault((min(a, b), max(a, b)), []).append((e, k, a > b))
        self.edge_star = {key: sorted(v) for key, v in edges.items()}       # canonical: by element
        self.edge_of = {}             # (e, k) -> (key, reversed)
        for key, cps in self.edge_star.items():
            for (e, k, rev) in cps: self.edge_of[(e, k)] = (key, rev)
        verts = {}
        for e in range(E):
            for c in range(8): verts.setdefault(conn[e][c], []).append((e, c))
        self.vert_star = {v: sorted(cps) for v, cps in verts.items()}

    def node_class(self):
        """For each local node: ('I'), ('F', f, i, j), ('E', k, t) or ('V', c)."""
        cls = [None] * (n ** 3)
        for f in range(6):
            for i in range(1, p):
                for j in range(1, p): cls[FACE_NODES[f][i, j]] = ('F', f, i, j)
        for k, (nodes, c0, c1) in enumerate(EDGES):
            for t in range(1, p): cls[nodes[t]] = ('E', k, t)
        for c in range(8): cls[CORNER_NODE[c]] = ('V', c)
        return [x if x is not None else ('I',) for x in cls]

def dss_gather(T, value, E_local=None):
    """value(e, l) -> copy l of element e (any element, local or ghost). Returns sums for the
    elements in E_local (all by default), every copy summed in canonical order."""
    cls = T.node_class()
    elems = range(len(T.conn)) if E_local is None else E_local
    out = {}
    for e in elems:
        y = np.empty(n ** 3)
        for l in range(n ** 3):
            c = cls[l]
            if c[0] == 'I':
                y[l] = value(e, l)
            elif c[0] == 'F':
                f, i, j = c[1:]
                q = T.partner[e, f]
                if q < 0:
                    y[l] = value(e, l)
                else:
                    pl = T.pmap[e, f, i, j]
                    lo, hi = ((e, l), (q, pl)) if e < q else ((q, pl), (e, l))
                    y[l] = value(*lo) + value(*hi)
            elif c[0] == 'E':
                k, t = c[1:]
                key, rev = T.edge_of[(e, k)]
                tc = p - t if rev else t                      # canonical position on the edge
                s = 0.0
                for (e2, k2, rev2) in T.edge_star[key]:
                    s += value(e2, EDGES[k2][0][p - tc if rev2 else tc])
                y[l] = s
            else:
                s = 0.0
                for (e2, c2) in T.vert_star[T.conn[e][c[1]]]:
                    s += value(e2, CORNER_NODE[c2])
                y[l] = s
        out[e] = y
    return out

def global_ids(T):
    """Global node identity of every (e, l), for the assembled reference."""
    cls = T.node_class(); ids = {}; gid = {}
    for e, ce in enumerate(T.conn):
        for l in range(n ** 3):
            c = cls[l]
            if c[0] == 'I': key = ('I', e, l)
            elif c[0] == 'V': key = ('V', ce[c[1]])
            elif c[0] == 'E':
                k, t = c[1:]; ek, rev = T.edge_of[(e, k)]; key = ('E', ek, p - t if rev else t)
            else:
                f, i, j = c[1:]
                v = [ce[cc] for cc in FACE_CORNERS[f]]
                key = ('F', tuple(sorted(v)), int(canonical_face_map(v)[i, j]))
            ids[(e, l)] = gid.setdefault(key, len(gid))
    return ids, len(gid)

def dss_entity(T, vals, counter=None):
    """Design B: one pass per entity type. Each face, edge and vertex node is summed ONCE
    over its copies (canonical order) and written to every copy; interior nodes keep their
    value. counter['reads'] counts value reads."""
    E = len(T.conn)
    y = vals.copy()
    reads = E * (n - 2) ** 3                      # interior nodes, read once
    for e in range(E):                            # faces: the lower element of each pair drives
        for f in range(6):
            q = T.partner[e, f]
            for i in range(1, p):
                for j in range(1, p):
                    l = FACE_NODES[f][i, j]
                    if q < 0:
                        reads += 1
                        continue
                    if e > q: continue
                    pl = T.pmap[e, f, i, j]
                    s = vals[e, l] + vals[q, pl]
                    reads += 2
                    y[e, l] = s; y[q, pl] = s
    for key, star in T.edge_star.items():
        for tc in range(1, p):
            s = 0.0
            for (e2, k2, rev2) in star:
                s += vals[e2, EDGES[k2][0][p - tc if rev2 else tc]]
            reads += len(star)
            for (e2, k2, rev2) in star:
                y[e2, EDGES[k2][0][p - tc if rev2 else tc]] = s
    for v, star in T.vert_star.items():
        s = 0.0
        for (e2, c2) in star: s += vals[e2, CORNER_NODE[c2]]
        reads += len(star)
        for (e2, c2) in star: y[e2, CORNER_NODE[c2]] = s
    if counter is not None: counter['reads'] = reads
    return y


def main():

    rng = np.random.default_rng(5)

    def block_mesh(NX, NY, NZ):
        vid = lambda i, j, k: (i * (NY + 1) + j) * (NZ + 1) + k
        conn = []
        for ex in range(NX):
            for ey in range(NY):
                for ez in range(NZ):
                    conn.append([vid(ex + (s[0] > 0), ey + (s[1] > 0), ez + (s[2] > 0)) for s in CR])
        return conn

    def extrude(quads, nv2d, layers):
        conn = []
        for k in range(layers):
            for q in quads:
                conn.append([v + nv2d * k for v in q] + [v + nv2d * (k + 1) for v in q])
        return conn

    TRI = [(0, 3, 6, 5), (3, 1, 4, 6), (6, 4, 2, 5)]                 # centre vertex 6: valence 3
    PENT = [(0, 6 + (i - 1) % 5, 1 + i, 6 + i) for i in range(5)]   # centre vertex 0: valence 5

    meshes = {
        "block 3x3x2": block_mesh(3, 3, 2),
        "triangle x3, 3 layers": extrude(TRI, 7, 3),
        "pentagon x5, 2 layers": extrude(PENT, 11, 2),
    }

    ok = True
    for name, conn in meshes.items():
        conn = [rotate_element(c, ROT[rng.integers(24)]) for c in conn]    # random local frames
        T = Tables(conn)
        E = len(conn)
        vals = rng.standard_normal((E, n ** 3))
        y = dss_gather(T, lambda e, l: vals[e, l])
        ids, G = global_ids(T)
        # 1. continuity: every copy of a global node holds the same bits
        copies = {}
        for (e, l), g in ids.items(): copies.setdefault(g, []).append((e, l))
        cont = all(len({y[e][l] for (e, l) in cps}) == 1 for cps in copies.values())
        # 2. bitwise equal to the assembled sum taken in canonical order (copies by element id)
        ref = {g: sum(vals[e, l] for (e, l) in sorted(cps)) for g, cps in copies.items()}
        bitwise = all(y[e][l] == ref[ids[(e, l)]] for (e, l) in ids)
        # 3. ... and to an ordinary scatter-add to rounding
        acc = np.zeros(G); np.add.at(acc, [ids[(e, l)] for e in range(E) for l in range(n ** 3)], vals.ravel())
        close = max(abs(y[e][l] - acc[ids[(e, l)]]) for (e, l) in ids) / np.abs(acc).max()
        ev = sorted({len(s) for s in T.edge_star.values()}); vv = sorted({len(s) for s in T.vert_star.values()})
        print(f"{name}: {E} elements, {G} global nodes; edge valences {ev}, vertex valences {vv}")
        print(f"  continuity {cont}, bit-identical to canonical-order assembly {bitwise}, vs scatter-add {close:.1e}")
        ok &= cont and bitwise and close < 1e-14

        # 4. several ranks: each rank gathers its own elements, reading remote copies ONLY if the
        #    single-phase exchange delivers them (copies of nodes this rank also holds).
        for P in (2, 3, 5):
            part = rng.integers(P, size=E)
            holders = {}
            for (e, l), g in ids.items(): holders.setdefault(g, set()).add(part[e])
            res, sent = {}, {}
            for r in range(P):
                mine = [e for e in range(E) if part[e] == r]
                def value(e, l, r=r):
                    if part[e] == r: return vals[e, l]
                    g = ids[(e, l)]
                    if r not in holders[g]:
                        raise RuntimeError("read a remote copy the exchange does not deliver")
                    sent[(part[e], r)] = sent.get((part[e], r), 0) + 1
                    return vals[e, l]
                res.update(dss_gather(T, value, mine))
            same = all(np.array_equal(res[e], y[e]) for e in range(E))
            ok &= same
            print(f"  {P} ranks: bit-identical to one domain {same}; remote reads (with repeats) {sum(sent.values())}")
        yb = dss_entity(T, vals)
        same_b = all(np.array_equal(y[e], yb[e]) for e in range(E))
        ok &= same_b
        print(f"  entity-centric form vs per-copy gather: bit-identical {same_b}")
    print("CELL-WISE UNSTRUCTURED REFERENCE:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
