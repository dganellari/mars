#!/usr/bin/env python3
# Shape statistics of a Hex8 Exodus mesh: edge aspect ratio, Jacobian condition number at the
# element center, minimum scaled Jacobian over the corners, and volume.
# Usage: hex_mesh_quality.py mesh.e [per_element.npz]
import sys
import numpy as np
import netCDF4 as nc

path = sys.argv[1]
f = nc.Dataset(path)
x = np.stack([f.variables[c][:].data for c in ("coordx", "coordy", "coordz")], axis=1)
conn_var = f.variables["connect1"]
ne = conn_var.shape[0]

# Exodus (= VTK) Hex8 corner order: 0-3 bottom face, 4-7 top face.
edges = np.array([[[0, 1], [3, 2], [4, 5], [7, 6]],
                  [[0, 3], [1, 2], [4, 7], [5, 6]],
                  [[0, 4], [1, 5], [2, 6], [3, 7]]])
# For each corner, its three edges along xi, eta, zeta, oriented to give a right-handed frame.
corner_frames = np.array([[0, 1, 3, 4], [1, 2, 0, 5], [2, 3, 1, 6], [3, 0, 2, 7],
                          [4, 7, 5, 0], [5, 4, 6, 1], [6, 5, 7, 2], [7, 6, 4, 3]])

ar = np.empty(ne); cond = np.empty(ne); sj = np.empty(ne); vol = np.empty(ne)
chunk = 1 << 20
for s in range(0, ne, chunk):
    e = min(ne, s + chunk)
    p = x[conn_var[s:e, :] - 1]                       # (n, 8, 3)
    vec = p[:, edges[:, :, 1]] - p[:, edges[:, :, 0]]  # (n, 3, 4, 3)
    lengths = np.linalg.norm(vec, axis=3).mean(axis=2) # mean edge length per direction
    ar[s:e] = lengths.max(axis=1) / lengths.min(axis=1)
    J = vec.mean(axis=2).transpose(0, 2, 1)            # columns = mean edge vectors = 2 dx/dxi at center
    sv = np.linalg.svd(J, compute_uv=False)
    cond[s:e] = sv[:, 0] / sv[:, 2]
    vol[s:e] = np.abs(np.linalg.det(J))
    a = p[:, corner_frames[:, 1]] - p[:, corner_frames[:, 0]]
    b = p[:, corner_frames[:, 2]] - p[:, corner_frames[:, 0]]
    c = p[:, corner_frames[:, 3]] - p[:, corner_frames[:, 0]]
    det = np.einsum("nki,nki->nk", np.cross(a, b), c)
    norm = np.linalg.norm(a, axis=2) * np.linalg.norm(b, axis=2) * np.linalg.norm(c, axis=2)
    sj[s:e] = (det / norm).min(axis=1)

def report(name, v, fmt="%.3g"):
    q = np.percentile(v, [0, 50, 90, 99, 99.9, 100])
    print("%-22s" % name + "  ".join(("%s " + fmt) % (l, t) for l, t in
          zip(("min", "p50", "p90", "p99", "p99.9", "max"), q)))

print("%s: %d Hex8, %d nodes" % (path, ne, x.shape[0]))
report("edge aspect ratio", ar)
report("cond(J) at center", cond)
report("min scaled Jacobian", sj)
report("volume (rel. median)", vol / np.median(vol))
for t in (2, 5, 10, 20, 50, 100, 1000):
    print("cond(J) > %-5d %6.2f %%   edge AR > %-5d %6.2f %%" %
          (t, 100.0 * np.mean(cond > t), t, 100.0 * np.mean(ar > t)))
print("scaled Jacobian < 0.5: %.2f %%, < 0.2: %.2f %%, <= 0: %d elements" %
      (100.0 * np.mean(sj < 0.5), 100.0 * np.mean(sj < 0.2), int(np.sum(sj <= 0))))

names = ["".join(ch.decode() for ch in row if ch) for row in f.variables["ss_names"][:].data]
for i, n in enumerate(names):
    el = np.unique(f.variables["elem_ss%d" % (i + 1)][:].data - 1)
    print("side set %-7s %6d elements: cond(J) p50 %.3g p99 %.3g max %.3g; edge AR p50 %.3g max %.3g" %
          (n, el.size, np.median(cond[el]), np.percentile(cond[el], 99), cond[el].max(),
           np.median(ar[el]), ar[el].max()))
if len(sys.argv) > 2:
    np.savez_compressed(sys.argv[2], ar=ar.astype(np.float32), cond=cond.astype(np.float32),
                        sj=sj.astype(np.float32), vol=vol.astype(np.float32))
