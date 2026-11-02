#!/usr/bin/env python3
"""Check of the node-line smoother's block matrices (backend/distributed/unstructured/solvers/
mars_cellwise_lines.hpp): the element lines and neighbour tables from the C++ host builder
(tests/cellwise_lines_host_check.cpp), then build_kernel's claims, runs and band assembly emulated
step for step (own element, lateral neighbours, corner couplings at valence-3 edges, assembled
diagonal), and every band entry compared with the assembled operator restricted to the block
(scripts/bp3_smoother_study.py). On one rank the GPU blocks are these, so any difference is a bug.

Cases: the cube sheared (deform 0.05, stretch 100), unstretched (lines in mixed directions) and
with graded boxes; with --mesh, blocks of a real Hex8 Exodus mesh around --seeds.
Run from the MARS root: python3 scripts/test_node_lines.py [--mesh large_hex_only.e --seeds 1593098,2107199]
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bp3_smoother_study as S  # noqa: E402

P = 3
n = P + 1
n3 = n ** 3


def other_axes(a):
    return (1 if a == 0 else 0), (1 if a == 2 else 2)


def transverse(code, i, j):
    if code & 1:
        i, j = j, i
    if code & 2:
        i = P - i
    if code & 4:
        j = P - j
    return i, j


def idx(x):
    return (x[0] * n + x[1]) * n + x[2]


def tables(tool, corners, vid, work):
    corners.astype(np.float64).tofile(os.path.join(work, "corners.bin"))
    vid.astype(np.int64).tofile(os.path.join(work, "vid.bin"))
    subprocess.run([tool, str(len(corners)), os.path.join(work, "corners.bin"), os.path.join(work, "vid.bin"), work],
                   check=True, capture_output=True)
    rd = lambda name, t: np.fromfile(os.path.join(work, name + ".bin"), dtype=t)
    return {k: rd(k, t) for k, t in (("start", np.int32), ("elem", np.int32), ("code", np.uint8),
                                     ("face_elem", np.int32), ("face_code", np.uint8), ("edge_start", np.int32),
                                     ("edge_elem", np.int32), ("edge_code", np.uint8))}


def check(T, corners, vid):
    """Largest |emulated - assembled| over all band entries, relative to the largest entry."""
    gidx, free = S.numbering(vid, P)
    lv = S.Level(P, corners, gidx, free)
    K = lv.matrices()
    keys, vals = lv.assembled()
    diag = lv.nodes(lv.diag)
    start, elem, code = T["start"], T["elem"], T["code"]
    chains = (len(start) - 1) * n * n

    def walk(chain):
        l = chain // (n * n)
        i0, j0 = (chain % (n * n)) // n, chain % n
        for m in range(start[l], start[l + 1]):
            cd = int(code[m]); a = cd & 3; rev = (cd >> 2) & 1
            i, j = transverse(cd >> 3, i0, j0)
            k1, k2 = other_axes(a)
            for s in range(0 if m == start[l] else 1, P + 1):
                x = [0, 0, 0]; x[a] = P - s if rev else s; x[k1] = i; x[k2] = j
                yield m, s, x[a], i, j, int(elem[m]) * n3 + idx(x)

    def laterals(m, i, j):
        out = []
        on1, on2 = i in (0, P), j in (0, P)
        if on1 and on2:
            g = (i == P) + 2 * (j == P)
            for q in range(T["edge_start"][m * 4 + g], T["edge_start"][m * 4 + g + 1]):
                ec = int(T["edge_code"][q])
                out.append((int(T["edge_elem"][q]), ec & 3, (ec >> 2 & 1) * P, (ec >> 3 & 1) * P, (ec >> 4) & 1))
        elif on1 or on2:
            f = (i == P) if on1 else 2 + (j == P)
            other = int(T["face_elem"][m * 4 + f])
            if other >= 0:
                fc = int(T["face_code"][m * 4 + f]); an = fc & 3; kn = (fc >> 2) & 3; t = j if on1 else i
                x = [0, 0, 0]; x[kn] = (fc >> 4 & 1) * P; x[3 - an - kn] = P - t if (fc >> 6 & 1) else t
                o1, o2 = other_axes(an)
                out.append((other, an, x[o1], x[o2], (fc >> 5) & 1))
        return out

    def local(an, c1, c2, flip, k):
        o1, o2 = other_axes(an)
        x = [0, 0, 0]; x[o1] = c1; x[o2] = c2; x[an] = P - k if flip else k
        return idx(x)

    owner = {}
    for chain in range(chains):   # claim_kernel: the lowest chain wins a node
        for _, _, _, _, _, t in walk(chain):
            g = gidx.ravel()[t]
            if free[g] and g not in owner:
                owner[g] = chain
    err = 0.0; scale = 0.0
    for chain in range(chains):   # build_kernel
        st = dict(rows=[], band={}, blocks=[], cur=-1, prev=False, ci=0, cj=0, seg=[-1] * n, kseg=[0] * n,
                  pcur=-1, pci=0, pcj=0, pseg=[-1] * n, pkseg=[0] * n)

        def add(entry, flip):
            seg, kseg = st["seg"], st["kseg"]
            for s1 in range(P):
                if seg[s1] < 0:
                    continue
                for s2 in range(s1 + 1, n):
                    if seg[s2] < 0:
                        break
                    k, l = (P - kseg[s1], P - kseg[s2]) if flip else (kseg[s1], kseg[s2])
                    key = (seg[s2], s2 - s1)
                    st["band"][key] = st["band"].get(key, 0.0) + entry(k, l)

        def flush():
            seg, kseg, pseg, pkseg = st["seg"], st["kseg"], st["pseg"], st["pkseg"]
            cur, ci, cj = st["cur"], st["ci"], st["cj"]
            if st["pcur"] >= 0:
                before = laterals(st["pcur"], st["pci"], st["pcj"])
                for (en, an, v1, v2, fl) in laterals(cur, ci, cj):
                    for (ep, ap, w1, w2, fp) in before:
                        if ep != en:
                            continue
                        for s1 in range(P):
                            if any(pseg[q] < 0 for q in range(s1, n)):
                                continue
                            for s2 in range(1, n):
                                if seg[s2] < 0:
                                    break
                                d = seg[s2] - pseg[s1]
                                if d > P:
                                    break
                                key = (seg[s2], d)
                                st["band"][key] = st["band"].get(key, 0.0) + K[en][local(ap, w1, w2, fp, pkseg[s1]),
                                                                                   local(an, v1, v2, fl, kseg[s2])]
            e = int(elem[cur]); a = int(code[cur]) & 3; k1, k2 = other_axes(a)

            def own(k, l):
                x = [0, 0, 0]; x[k1] = ci; x[k2] = cj; x[a] = k; y = list(x); y[a] = l
                return K[e][idx(x), idx(y)]
            add(own, False)
            for (en, an, v1, v2, fl) in laterals(cur, ci, cj):
                add(lambda k, l, en=en, an=an, v1=v1, v2=v2: K[en][local(an, v1, v2, 0, k), local(an, v1, v2, 0, l)], fl)

        for m, s, k, i, j, t in walk(chain):
            if m != st["cur"]:
                if st["cur"] >= 0:
                    flush()
                    st.update(pcur=st["cur"], pci=st["ci"], pcj=st["cj"], pseg=list(st["seg"]), pkseg=list(st["kseg"]))
                st["seg"] = [-1] * n
                if st["cur"] >= 0:
                    st["seg"][0] = len(st["rows"]) - 1 if st["prev"] else -1
                    st["kseg"][0] = P if (int(code[m]) >> 2) & 1 else 0
                st.update(cur=m, ci=i, cj=j)
            g = gidx.ravel()[t]
            won = bool(free[g]) and owner.get(g) == chain
            st["kseg"][s] = k
            if won:
                if not st["prev"]:
                    st["blocks"].append(len(st["rows"]))
                st["seg"][s] = len(st["rows"]); st["rows"].append(g)
            st["prev"] = won
        if st["cur"] >= 0:
            flush()
        rows = st["rows"]
        for b0, b1 in zip(st["blocks"], st["blocks"][1:] + [len(rows)]):
            for r in range(b0, b1):
                for d in range(0, min(P, r - b0) + 1):
                    kk = rows[r] * lv.N + rows[r - d]
                    pos = min(np.searchsorted(keys, kk), len(keys) - 1)
                    exact = vals[pos] if keys[pos] == kk else 0.0
                    emul = diag[rows[r]] if d == 0 else st["band"].get((r, d), 0.0)
                    err = max(err, abs(emul - exact)); scale = max(scale, abs(exact))
    return err / scale


def block_of(mesh, seed, layers):
    import netCDF4 as nc
    f = nc.Dataset(mesh)
    X = np.stack([f.variables[c][:].data for c in ("coordx", "coordy", "coordz")], 1)
    conn = np.concatenate([f.variables[f"connect{b + 1}"][:].data for b in range(len(f.dimensions["num_el_blk"]))]) - 1
    flat = conn.ravel(); order = np.argsort(flat, kind="stable")
    start = np.searchsorted(flat[order], np.arange(len(X) + 1))
    elems = {seed}; frontier = [seed]
    for _ in range(layers):
        nodes = np.unique(conn[frontier])
        near = set((np.concatenate([order[start[v]:start[v + 1]] for v in nodes]) // 8).tolist())
        frontier = list(near - elems); elems |= near
    B = np.array(sorted(elems))
    return X[conn[B]], conn[B].astype(np.int64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", help="built tests/cellwise_lines_host_check (default: compile it)")
    ap.add_argument("--mesh")
    ap.add_argument("--seeds", default="1593098,2107199,15109394")
    ap.add_argument("--layers", type=int, default=3)
    args = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with tempfile.TemporaryDirectory() as work:
        tool = args.tool
        if not tool:
            tool = os.path.join(work, "cellwise_lines_host_check")
            subprocess.run(["c++", "-std=c++20", "-O2", "-I" + root, os.path.join(root, "tests/cellwise_lines_host_check.cpp"),
                            "-o", tool], check=True)
        cases = [("cube 6^3, deform 0.05, stretch 100", S.cube(6, 0.05, 100.0)),
                 ("cube 6^3, deform 0.05, stretch 1 (mixed directions)", S.cube(6, 0.05, 1.0)),
                 ("boxes 6^3, stretch 100, graded 10x", S.cube(6, 0.0, 100.0, 2.302585))]
        if args.mesh:
            cases += [(f"{os.path.basename(args.mesh)}, {args.layers} layers around {s}",
                       block_of(args.mesh, int(s), args.layers)) for s in args.seeds.split(",")]
        bad = 0
        for name, (corners, vid) in cases:
            rel = check(tables(tool, corners, vid, work), corners, vid)
            ok = rel <= 1e-12
            bad += not ok
            print(f"{name}: {len(corners)} elements, max |emulated - assembled| / max = {rel:.1e}  {'PASS' if ok else 'FAIL'}")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
