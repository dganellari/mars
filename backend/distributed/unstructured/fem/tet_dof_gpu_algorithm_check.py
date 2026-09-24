"""Validate the DEVICE algorithm for tet HO DOF numbering against the host
std::map reference, before any CUDA is written.

Host (mars_ho_dof_handler_tet.hpp): key = sorted [(gid, wq)] with zero weights
dropped, deduped through std::map, ids assigned in FIRST-ENCOUNTER order.

Device plan: pack the same key into a fixed 4-slot 256-bit record, sort it with
stable LSD radix over uint64 chunks, mark unique boundaries, scan to ids, scatter.
Device ids come out in SORTED-key order, so the two numberings are a permutation
-- exactly as on the hex path. What must match is the partition: which (elem,node)
slots share a DOF.
"""
import numpy as np

Q = 4294967296.0          # 2^32, as in the header
UNUSED_GID = (1 << 31) - 1  # sorts last, so short keys pad canonically


def kuhn_mesh(ncells):
    """6 tets per cell, all sharing the cell diagonal (the header's test mesh)."""
    def vid(x, y, z): return (x * (ncells + 1) + y) * (ncells + 1) + z
    KUHN = [(0,1,3,7),(0,1,5,7),(0,2,3,7),(0,2,6,7),(0,4,5,7),(0,4,6,7)]
    corner = lambda c: ((c >> 2) & 1, (c >> 1) & 1, c & 1)
    ec = []
    for x in range(ncells):
        for y in range(ncells):
            for z in range(ncells):
                for t in KUHN:
                    ec.append([vid(x + corner(c)[0], y + corner(c)[1], z + corner(c)[2]) for c in t])
    return np.array(ec, dtype=np.int64)


def node_weights(Z):
    """Collapsed (Duffy) barycentric weights per tensor node, as in the header."""
    n = len(Z)
    out = np.zeros((n, n, n, 4))
    for ia in range(n):
        for ib in range(n):
            for ic in range(n):
                a, b, c = Z[ia], Z[ib], Z[ic]
                r, s, t = a * (1 - b) * (1 - c), b * (1 - c), c
                out[ia, ib, ic] = (1.0 - r - s - t, r, s, t)
    return out


def host_reference(ec, W, n):
    """std::map<vector<pair<int,long long>>, int>, first-encounter ids."""
    NN = n * n * n
    elemDof = np.full(len(ec) * NN, -1, dtype=np.int64)
    dofOf, nxt = {}, 0
    for e, corners in enumerate(ec):
        g = np.sort(corners)
        for ia in range(n):
            for ib in range(n):
                for ic in range(n):
                    w = W[ia, ib, ic]
                    key = tuple(sorted((int(g[c]), int(round(w[c] * Q)))
                                       for c in range(4) if w[c] > 1e-9))
                    if key not in dofOf:
                        dofOf[key] = nxt
                        nxt += 1
                    elemDof[e * NN + (ia * n + ib) * n + ic] = dofOf[key]
    return elemDof, nxt


def device_plan(ec, W, n):
    """The device algorithm: fixed 4-slot packed key -> sort -> unique -> scatter."""
    NN = n * n * n
    N = len(ec) * NN
    # 1. pack: 4 slots of (gid, wq), zero weights pushed to UNUSED and sorted last
    keys = np.empty((N, 8), dtype=np.uint64)
    for e, corners in enumerate(ec):
        g = np.sort(corners)
        for ia in range(n):
            for ib in range(n):
                for ic in range(n):
                    w = W[ia, ib, ic]
                    pairs = [(int(g[c]), int(round(w[c] * Q))) for c in range(4) if w[c] > 1e-9]
                    pairs.sort()
                    while len(pairs) < 4:
                        pairs.append((UNUSED_GID, 0))
                    row = e * NN + (ia * n + ib) * n + ic
                    keys[row] = [v for p in pairs for v in p]
    # 2. stable LSD radix over the 8 uint64 lanes, least-significant lane first
    order = np.arange(N)
    for lane in range(7, -1, -1):
        order = order[np.argsort(keys[order, lane], kind="stable")]
    srt = keys[order]
    # 3. unique boundaries -> dense ids -> scatter
    newkey = np.ones(N, dtype=np.int64)
    newkey[1:] = (srt[1:] != srt[:-1]).any(axis=1)
    ids = np.cumsum(newkey) - 1
    elemDof = np.empty(N, dtype=np.int64)
    elemDof[order] = ids
    return elemDof, int(ids[-1]) + 1


def partitions_match(a, b):
    """Same identification classes? (ids are a permutation, the partition is not)"""
    ia, ib = {}, {}
    for x, y in zip(a, b):
        if ia.setdefault(x, y) != y or ib.setdefault(y, x) != x:
            return False
    return True


for ncells, n in ((2, 2), (2, 3), (3, 3), (2, 4), (3, 4), (2, 5)):
    Z = np.linspace(0.0, 1.0, n)          # stand-in for the GLL nodes
    ec, W = kuhn_mesh(ncells), node_weights(Z)
    h, nh = host_reference(ec, W, n)
    d, nd = device_plan(ec, W, n)
    ok = (nh == nd) and partitions_match(h, d)
    print(f"ncells={ncells} n={n} nElem={len(ec):5d} | host numDof={nh:6d} device={nd:6d} | "
          f"{'PASS' if ok else 'FAIL'}")


# ---------------------------------------------------------------------------
# NODAL handler (HoTetDofHandler): integer barycentric weights, corners NOT
# pre-sorted -- the key sort alone canonicalizes. Separate path, separate risk.
# ---------------------------------------------------------------------------
def bary_nodes(P):
    """All (i,j,k) with i+j+k <= P, the Np nodal barycentric indices."""
    return [(i, j, k) for i in range(P + 1) for j in range(P + 1 - i)
            for k in range(P + 1 - i - j)]


def nodal_host(ec, P):
    bary = bary_nodes(P)
    Np = len(bary)
    elemDof = np.full(len(ec) * Np, -1, dtype=np.int64)
    dofOf, nxt = {}, 0
    for e, g in enumerate(ec):                    # NOTE: g is NOT sorted
        for m, (i, j, k) in enumerate(bary):
            w = (P - i - j - k, i, j, k)
            key = tuple(sorted((int(g[c]), w[c]) for c in range(4) if w[c] > 0))
            if key not in dofOf:
                dofOf[key] = nxt; nxt += 1
            elemDof[e * Np + m] = dofOf[key]
    return elemDof, nxt, Np


def nodal_device(ec, P):
    bary = bary_nodes(P)
    Np = len(bary)
    N = len(ec) * Np
    keys = np.empty((N, 4), dtype=np.uint64)      # 4 packed lanes: gid<<32 | w
    for e, g in enumerate(ec):
        for m, (i, j, k) in enumerate(bary):
            w = (P - i - j - k, i, j, k)
            slots = []
            for c in range(4):
                slots.append((int(g[c]) << 32 | w[c]) if w[c] > 0
                             else (((1 << 31) - 1) << 32))
            slots.sort()
            keys[e * Np + m] = slots
    order = np.arange(N)
    for lane in range(3, -1, -1):
        order = order[np.argsort(keys[order, lane], kind="stable")]
    srt = keys[order]
    newkey = np.ones(N, dtype=np.int64)
    newkey[1:] = (srt[1:] != srt[:-1]).any(axis=1)
    ids = np.cumsum(newkey) - 1
    elemDof = np.empty(N, dtype=np.int64)
    elemDof[order] = ids
    return elemDof, int(ids[-1]) + 1


print()
for ncells, P in ((2, 1), (2, 2), (3, 2), (2, 3), (3, 3), (2, 4)):
    ec = kuhn_mesh(ncells)
    h, nh, Np = nodal_host(ec, P)
    d, nd = nodal_device(ec, P)
    ok = (nh == nd) and partitions_match(h, d)
    print(f"NODAL ncells={ncells} P={P} Np={Np:3d} nElem={len(ec):5d} | "
          f"host numDof={nh:6d} device={nd:6d} | {'PASS' if ok else 'FAIL'}")
