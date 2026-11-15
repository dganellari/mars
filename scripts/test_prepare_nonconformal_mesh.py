import importlib.util
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from scipy.io import netcdf_file

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).with_name('prepare_nonconformal_mesh.py'))
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def fixture(path, defect=None):
    with netcdf_file(path, 'w', version=2) as f:
        dims = dict(num_dim=3, num_el_blk=2, num_nodes=8, num_elem=2, four=4, one=1, num_side_sets=5, length=32)
        for name, value in dims.items():
            f.createDimension(name, value)
        x = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 1],
                      [2, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 1]], dtype=float)
        if defect == 'nonfinite':
            x[0, 0] = np.nan
        if defect == 'degenerate':
            x[0] = x[1]
        for d, axis in enumerate('xyz'):
            f.createVariable('coord' + axis, 'd', ('num_nodes',))[:] = x[:, d]
        for i in (1, 2):
            v = f.createVariable('connect' + str(i), 'i', ('one', 'four'))
            v.elem_type = b'HEX8' if defect == 'hex' and i == 2 else b'TETRA'
            values = np.arange(4 * (i - 1) + 1, 4 * i + 1)
            if defect == 'shared_nodes' and i == 2:
                values[1:] = [2, 3, 4]
            if defect == 'bad_index' and i == 2:
                values[0] = 9
            v[:] = values[None, :]
        labels = ['inlet', 'outlet', 'interface_left', 'interface_right', 'wall']
        ss = [[(1, 1)], [(2, 1)], [(1, 2)], [(2, 2)], [(1, 3), (1, 4), (2, 3), (2, 4)]]
        if defect == 'missing_face':
            ss[-1].pop()
        if defect == 'duplicate_face':
            ss[-1].append((1, 1))
        if defect == 'wrong_block':
            ss[2], ss[3] = ss[3], ss[2]
        if defect == 'bad_side':
            ss[2] = [(1, 5)]
        v = f.createVariable('ss_names', 'c', ('num_side_sets', 'length'))
        v[:] = np.asarray([list(s.ljust(32)) for s in labels], dtype='S1')
        for i, entries in enumerate(ss, 1):
            dim = 'ss' + str(i)
            f.createDimension(dim, len(entries))
            f.createVariable('elem_ss' + str(i), 'i', (dim,))[:] = [a for a, _ in entries]
            f.createVariable('side_ss' + str(i), 'i', (dim,))[:] = [b for _, b in entries]


class MeshPreparation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.mesh = Path(self.tmp.name) / 'test.exo'

    def test_source_ids_and_faces_preserved(self):
        fixture(self.mesh)
        arrays = prepare.prepare(self.mesh)
        x, cells, k, boundary, side, tags = arrays
        np.testing.assert_array_equal(cells, np.arange(8).reshape(2, 4))
        np.testing.assert_array_equal(x[1:4], x[5:8])
        np.testing.assert_array_equal(side, [0, 0, 0, 0, 1, 1, 1, 1])
        np.testing.assert_array_equal(tags, [[0, 0, 0], [1, 0, 1]])
        np.testing.assert_array_equal(boundary, [1, 1, 0, 1, 2, 2, 0, 2])
        self.assertEqual(prepare.OPPOSITE, (2, 0, 1, 3))
        out = Path(self.tmp.name) / 'mesh.mnc'
        prepare.write_mesh(out, arrays)
        content = out.read_bytes()
        self.assertEqual(content[:8], b'MARSNC01')
        self.assertEqual(struct.unpack('<QQQ', content[8:32]), (8, 2, 2))
        self.assertEqual(len(content), 32 + 32 * 8 + 24 * 2 + 12 * 2)
        with self.assertRaises(FileExistsError):
            prepare.write_mesh(out, arrays)

    def test_flow_boundary_labels(self):
        fixture(self.mesh)
        arrays = prepare.prepare(self.mesh, flow_boundaries=True)
        self.assertEqual(len(arrays), 7)
        np.testing.assert_array_equal(arrays[6], [[0, 2, 0], [1, 2, 1], [0, 1, 2], [0, 3, 2], [1, 1, 2], [1, 3, 2]])
        out = Path(self.tmp.name) / 'flow.mnc'
        prepare.write_mesh(out, arrays)
        content = out.read_bytes()
        self.assertEqual(content[:8], b'MARSNC02')
        self.assertEqual(struct.unpack('<QQQQ', content[8:40]), (8, 2, 2, 6))
        self.assertEqual(len(content), 40 + 32*8 + 24*2 + 12*2 + 12*6)

    def test_input_rejections(self):
        for defect in ['hex', 'shared_nodes', 'bad_index', 'nonfinite', 'degenerate',
                       'missing_face', 'duplicate_face', 'wrong_block', 'bad_side']:
            with self.subTest(defect=defect):
                fixture(self.mesh, defect)
                with self.assertRaises(ValueError):
                    prepare.prepare(self.mesh)

    def test_diffusivity_and_names(self):
        fixture(self.mesh)
        for k in (0, -1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                prepare.prepare(self.mesh, k_left=k)
        with self.assertRaises(ValueError):
            prepare.prepare(self.mesh, wall='inlet')
        with self.assertRaises(ValueError):
            prepare.prepare(self.mesh, interface_left='absent')
        np.testing.assert_array_equal(prepare.prepare(self.mesh, k_right=4)[2], [1, 4])


if __name__ == '__main__':
    unittest.main()
