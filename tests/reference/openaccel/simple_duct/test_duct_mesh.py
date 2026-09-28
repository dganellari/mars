#!/usr/bin/env python3
"""duct_mesh.py: lattice topology, orientation, boundary tagging and the Exodus/netCDF bytes.

Independent of the writer: a separate netCDF (CDF-2) parser reads the file back; every tet has
positive volume; exactly the manifold exterior faces are tagged, each on its own plane; the
generator is deterministic and refuses to overwrite. With --cxx BINARY the C++ mirror
(duct_host_run --dump-mesh) must produce the same canonical text: sizes, connectivity, side
sets and line order exactly, coordinates within Lattice.coordinate_tolerance() (both languages
evaluate the same contraction-free operations, but another toolchain may round differently).
"""
import array
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import duct_mesh as dm  # noqa: E402

CXX = None


def contents(path, mode="rb"):
    with open(path, mode) as f:
        return f.read()


def canonical_difference(a, b):
    """Largest coordinate difference between two canonical dumps. Every other token (sizes,
    connectivity, side sets, line order) must match exactly, else AssertionError."""
    la, lb = a.rstrip("\n").split("\n"), b.rstrip("\n").split("\n")
    if len(la) != len(lb):
        raise AssertionError("dumps differ in length: %d vs %d lines" % (len(la), len(lb)))
    worst = 0.0
    for n, (x, y) in enumerate(zip(la, lb)):
        tx, ty = x.split(), y.split()
        if tx[:1] == ["n"] and ty[:1] == ["n"] and len(tx) == len(ty) == 4:
            worst = max(worst, max(abs(float(p) - float(q)) for p, q in zip(tx[1:], ty[1:])))
        elif x != y:
            raise AssertionError("line %d differs: %r vs %r" % (n + 1, x, y))
    return worst


def read_netcdf(path):
    """Minimal CDF-1/CDF-2 reader: {dims}, {global attrs}, {var: (dims, type, attrs, values)}."""
    data = contents(path)
    pos = [0]

    def take(fmt):
        size = struct.calcsize(fmt)
        out = struct.unpack(">" + fmt, data[pos[0]:pos[0] + size])
        pos[0] += size
        return out

    def name():
        (n,) = take("i")
        s = data[pos[0]:pos[0] + n].decode("ascii")
        pos[0] += n + (4 - n % 4) % 4
        return s
    codes = {1: ("b", 1), 2: ("c", 1), 3: ("h", 2), 4: ("i", 4), 5: ("f", 4), 6: ("d", 8)}

    def attrs():
        tag, count = take("ii")
        assert tag in (0, 12)
        out = {}
        for _ in range(count):
            n = name()
            kind, m = take("ii")
            code, size = codes[kind]
            raw = data[pos[0]:pos[0] + m * size]
            pos[0] += m * size + (4 - m * size % 4) % 4
            out[n] = raw.decode("ascii") if kind == 2 else list(struct.unpack(">%d%s" % (m, code), raw))
        return out
    magic = data[:4]
    assert magic in (b"CDF\x01", b"CDF\x02"), magic
    pos[0] = 4
    take("i")   # numrecs
    tag, count = take("ii")
    assert tag == 10
    dims = []
    for _ in range(count):
        n = name()
        dims.append((n, take("i")[0]))
    gatts = attrs()
    tag, count = take("ii")
    assert tag == 11
    variables = {}
    for _ in range(count):
        n = name()
        (rank,) = take("i")
        ids = take("%di" % rank) if rank else ()
        a = attrs()
        kind, vsize = take("ii")
        begin = take("q" if magic == b"CDF\x02" else "i")[0]
        code, size = codes[kind]
        total = 1
        for d in ids:
            total *= dims[d][1]
        raw = data[begin:begin + total * size]
        values = raw if kind == 2 else array.array(code, raw)
        if kind != 2 and sys.byteorder == "little":
            values.byteswap()
        variables[n] = ([dims[d][0] for d in ids], kind, a, values)
    return dict(dims), gatts, variables


class Mesh(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="duct-mesh-")

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_rejects_invalid_lattices(self):
        for args in ((3,), (4, 7.0, 2.0, 1.0, 3.0), (2, 7.0, 1.5, 1.0, 2.0), (0,)):
            with self.assertRaises(ValueError):
                dm.Lattice(*args)

    def test_topology_orientation_and_tags(self):
        for cells, args in ((4, ()), (2, (3.0, 2.0, 1.0, 1.0)), (4, (3.0, 1.0, 2.0, 2.0)), (6, (7.5, 2.0, 1.0, 2.5))):
            lat = dm.Lattice(cells, *args)
            xs, ys, zs = lat.coordinates()
            conn = lat.connectivity()
            self.assertEqual(len(conn), 4 * lat.elements)
            volume = 0.0
            faces = {}
            for e in range(lat.elements):
                n = conn[4 * e:4 * e + 4]
                self.assertEqual(len(set(n)), 4)
                d = [[c[n[r]] - c[n[0]] for c in (xs, ys, zs)] for r in (1, 2, 3)]
                det = (d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) - d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0])
                       + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]))
                self.assertGreater(det, 0)
                volume += det / 6
                for o, f in enumerate(dm.FACE_NODES):
                    faces.setdefault(tuple(sorted(n[j] for j in f)), []).append((e, o))
            self.assertAlmostEqual(volume, lat.length * lat.width * lat.height, places=10)
            exterior = set(v[0] for v in faces.values() if len(v) == 1)
            self.assertTrue(all(len(v) <= 2 for v in faces.values()))
            sets = lat.side_sets()
            tagged = [f for name in dm.SIDE_SETS for f in sets[name]]
            self.assertEqual(sorted(tagged), sorted(exterior))
            self.assertEqual(dict((k, len(v)) for k, v in sets.items()), lat.boundary_faces())
            for name, planes in (("inlet", [(xs, 0.0)]), ("outlet", [(xs, lat.length)]),
                                 ("walls", [(ys, -lat.width / 2), (ys, lat.width / 2), (zs, -lat.height / 2), (zs, lat.height / 2)])):
                for e, o in sets[name]:
                    pts = [conn[4 * e + j] for j in dm.FACE_NODES[o]]
                    self.assertTrue(any(all(abs(c[g] - v) < 1e-12 for g in pts) for c, v in planes), (name, e, o))

    def test_exodus_bytes(self):
        lat = dm.Lattice(4)
        path = os.path.join(self.dir, "duct-4.exo")
        self.assertEqual(dm.main(["--cells", "4", "--output", path]), 0)
        dims, gatts, v = read_netcdf(path)
        self.assertEqual((dims["num_nodes"], dims["num_elem"], dims["num_el_blk"], dims["num_nod_per_el1"], dims["num_dim"]),
                         (lat.nodes, lat.elements, 1, 4, 3))
        self.assertEqual(gatts["floating_point_word_size"], [8])
        xs, ys, zs = lat.coordinates()
        self.assertEqual(list(v["coordx"][3]), list(xs))
        self.assertEqual(list(v["coordy"][3]), list(ys))
        self.assertEqual(list(v["coordz"][3]), list(zs))
        self.assertEqual(v["connect1"][2]["elem_type"], "TETRA4")
        self.assertEqual(list(v["connect1"][3]), [g + 1 for g in lat.connectivity()])
        names = v["ss_names"][3]
        self.assertEqual([names[33 * s:33 * s + 33].rstrip(b"\0").decode() for s in range(3)], list(dm.SIDE_SETS))
        sets = lat.side_sets()
        for s, name in enumerate(dm.SIDE_SETS):
            self.assertEqual(list(v["elem_ss%d" % (s + 1)][3]), [e + 1 for e, _ in sets[name]])
            self.assertEqual(list(v["side_ss%d" % (s + 1)][3]), [o + 1 for _, o in sets[name]])
        info = json.loads(contents(path[:-4] + ".json", "r"))
        self.assertEqual((info["nodes"], info["elements"], info["side_sets"]), (lat.nodes, lat.elements, lat.boundary_faces()))
        self.assertEqual(info["sha256"], dm.sha256(path))

    def test_deterministic_and_no_overwrite(self):
        a, b = os.path.join(self.dir, "a.exo"), os.path.join(self.dir, "b.exo")
        self.assertEqual(dm.main(["--cells", "4", "--output", a]), 0)
        self.assertEqual(dm.main(["--cells", "4", "--output", b]), 0)
        self.assertEqual(contents(a), contents(b))
        self.assertEqual(dm.main(["--cells", "4", "--output", a]), 1)
        self.assertEqual(dm.main(["--cells", "5", "--output", os.path.join(self.dir, "c.exo")]), 1)

    def test_ncdump_reads_the_file(self):
        ncdump = shutil.which("ncdump")
        if not ncdump:
            self.skipTest("ncdump not installed")
        path = os.path.join(self.dir, "duct-4.exo")
        dm.main(["--cells", "4", "--output", path])
        header = subprocess.check_output([ncdump, "-h", path], universal_newlines=True)
        self.assertIn("num_nodes = 675 ;", header)
        self.assertIn('connect1:elem_type = "TETRA4" ;', header)

    def test_cxx_mirror_is_identical(self):
        if CXX is None:
            self.skipTest("C++ mirror not requested (--cxx BINARY)")
        for cells, extra, lattice in ((4, [], dm.Lattice(4)), (6, ["--length", "7.5", "--stretch", "2.5"], dm.Lattice(6, 7.5, 2.0, 1.0, 2.5))):
            py, cx = os.path.join(self.dir, "py-%d.txt" % cells), os.path.join(self.dir, "cxx-%d.txt" % cells)
            self.assertEqual(dm.main(["--cells", str(cells), "--dump", py] + extra), 0)
            subprocess.check_call([CXX, "--cells", str(cells), "--dump-mesh", cx] + extra)
            self.assertLessEqual(canonical_difference(contents(py, "r"), contents(cx, "r")), lattice.coordinate_tolerance())

    def test_canonical_comparison_is_strict_where_it_must_be(self):
        text = dm.canonical(dm.Lattice(4))
        tol = dm.Lattice(4).coordinate_tolerance()
        self.assertEqual(canonical_difference(text, text), 0.0)
        # One ulp in a coordinate is rounding; a swapped node, another face or a moved node is not.
        lines = text.split("\n")
        n = next(i for i, l in enumerate(lines) if l.startswith("n ") and float(l.split()[2]) != 0)
        v = float(lines[n].split()[2])
        bumped = lines[:n] + ["n %s %.17g %s" % (lines[n].split()[1], v + abs(v) * sys.float_info.epsilon, lines[n].split()[3])] + lines[n + 1:]
        self.assertLessEqual(canonical_difference(text, "\n".join(bumped)), tol)
        moved = lines[:n] + ["n %s %.17g %s" % (lines[n].split()[1], v + 1e-9, lines[n].split()[3])] + lines[n + 1:]
        self.assertGreater(canonical_difference(text, "\n".join(moved)), tol)
        e = next(i for i, l in enumerate(lines) if l.startswith("e "))
        a = lines[e].split()
        swapped = lines[:e] + [" ".join([a[0], a[1], a[2], a[4], a[3]])] + lines[e + 1:]
        with self.assertRaises(AssertionError):
            canonical_difference(text, "\n".join(swapped))
        s = next(i for i, l in enumerate(lines) if l.startswith("s "))
        with self.assertRaises(AssertionError):
            canonical_difference(text, "\n".join(lines[:s] + lines[s + 1:]))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--cxx":
        CXX = sys.argv[2]
        del sys.argv[1:3]
    unittest.main()
