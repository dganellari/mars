"""Transport/coverage tests; fixture numbers are hand-worked, not reference runs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from openaccel_node_check import CONTRACT, STAGES, load_nodes, pack_nodes


def fixture_rows():
    diagonal = [16, 9, -7, 3, 32, 4, -8, 5, 48]
    return [
        dict(inputs=dict(density=[2], volume=[3], pseudo_dt=[.5], mass_divergence=[-4],
                         velocity=[1, -2, 3], pressure_gradient=[2, 4, -1], force=[1, 2, 3], source=[-1, 1, 2],
                         coriolis=[0]*9),
             outputs=dict(lhs=[16, 0, 0, 0, 16, 0, 0, 0, 16], rhs=[-10, 5, 6])),
        dict(inputs=dict(lhs=[4, 9, -7, 3, 8, 4, -8, 5, 12], alpha=[.25]), outputs=dict(lhs=diagonal)),
        dict(inputs=dict(volume=[8], row_blocks=[-4, 1e8, -1e8, 1e8, -8, 1e8, 1e8, 1e8, -12]+diagonal,
                         diagonal_block=[1], consistent=[1], fractional_step=[0], transient=[0],
                         small=[2.220446049250313e-16]),
             outputs=dict(d=[.5, .25, 1/6], d_tilde=[2/3, 1/3, 2/9])),
        dict(inputs=dict(rhs=[8, -4, 0], factor=[.75]), outputs=dict(rhs=[6, -3, 0])),
    ]


def write_fixture(path, rows=None, ranks=1):
    rows = fixture_rows() if rows is None else rows
    for s, stage in enumerate(STAGES):
        for rank in range(ranks):
            head = dict(schema=1, kind='header', producer='openaccel', fixture='public_channel',
                        reference_revision=CONTRACT['reference']['revision'],
                        solver_revision=CONTRACT['reference']['solver_gitlink']['revision'],
                        stage=stage, call=1, rank=rank, ranks=ranks)
            node_rows = [dict(rows[s], kind='node', id=42)] if rank == 0 else []
            values = [head]+node_rows+[dict(kind='end', records=len(node_rows))]
            (path/(stage+'.rank'+str(rank)+'.jsonl')).write_text(''.join(json.dumps(v)+'\n' for v in values))


class NodeCaptureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.data = self.root/'nodes'
        self.data.mkdir()
        write_fixture(self.data)

    def tearDown(self):
        self.tmp.cleanup()

    def mutate(self, stage, action):
        path = self.data/(STAGES[stage]+'.rank0.jsonl')
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        action(rows)
        path.write_text(''.join(json.dumps(v)+'\n' for v in rows))

    def test_pack(self):
        pack_nodes(self.data, self.root/'inputs.txt')
        data = (self.root/'inputs.txt').read_text()
        self.assertTrue(data.startswith('MARS_PUBLIC_NODE_REPLAY_V1 4\n'))
        self.assertEqual(json.loads((self.root/'inputs.json').read_text())['records'], 4)
        with self.assertRaises(ValueError):
            pack_nodes(self.data, self.root/'inputs.txt')

    def test_empty_rank(self):
        write_fixture(self.data, ranks=2)
        self.assertEqual(len(load_nodes(self.data)[0]), 4)
        (self.data/(STAGES[0]+'.rank1.jsonl')).unlink()
        with self.assertRaisesRegex(ValueError, 'missing rank'):
            load_nodes(self.data)

    def test_rejected_mutations(self):
        cases = [
            (0, lambda r: r[0].update(reference_revision='wrong')),
            (0, lambda r: r[0].update(fixture='private')),
            (0, lambda r: r[1]['inputs']['coriolis'].__setitem__(0, 1)),
            (0, lambda r: r[1]['inputs']['density'].__setitem__(0, -1)),
            (1, lambda r: r[1]['inputs']['alpha'].__setitem__(0, 0)),
            (2, lambda r: r[1]['inputs']['fractional_step'].__setitem__(0, 1)),
            (2, lambda r: r[1]['inputs']['transient'].__setitem__(0, 1)),
            (2, lambda r: r[1]['inputs']['small'].__setitem__(0, 0)),
            (2, lambda r: r[1]['inputs']['row_blocks'].__setitem__(9, 4)),
            (2, lambda r: r[1]['inputs']['row_blocks'].pop()),
            (2, lambda r: r[1]['inputs']['diagonal_block'].__setitem__(0, 3)),
            (3, lambda r: r[1]['outputs']['rhs'].__setitem__(0, float('nan'))),
            (3, lambda r: r[1].update(id=43)),
            (3, lambda r: r[-1].update(records=5)),
        ]
        for stage, action in cases:
            with self.subTest(stage=stage, action=action):
                write_fixture(self.data)
                self.mutate(stage, action)
                with self.assertRaises(ValueError):
                    load_nodes(self.data)

    def test_duplicate_boundary_node(self):
        def duplicate(rows):
            rows.insert(2, copy.deepcopy(rows[1]))
            rows[-1]['records'] = 2
        self.mutate(3, duplicate)
        with self.assertRaisesRegex(ValueError, 'duplicate owned'):
            load_nodes(self.data)

    def test_missing_stage(self):
        (self.data/(STAGES[1]+'.rank0.jsonl')).unlink()
        with self.assertRaisesRegex(ValueError, 'missing node stages'):
            load_nodes(self.data)

    def test_influence_input_is_not_an_oracle(self):
        # Finite altered output is structurally accepted; only production replay can reject it.
        self.mutate(2, lambda r: r[1]['outputs']['d'].__setitem__(0, .125))
        self.assertEqual(load_nodes(self.data)[0][(2, 1, 42)]['outputs']['d'][0], .125)


if __name__ == '__main__':
    unittest.main()
