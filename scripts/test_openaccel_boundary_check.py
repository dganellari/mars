"""Boundary transport tests with hand-worked synthetic blocks, not reference runs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from openaccel_boundary_check import CONTRACT, STAGES, load_boundary, pack_boundary


def fixture_rows():
    rows = []
    for stage in range(6):
        wall = stage == 5
        width = 9 if wall else (4 if stage < 3 else 12)
        shape = [11/18 if s == n else 7/36 for s in range(3) for n in range(3)]
        x = dict(face_nodes=list(range(3)) if wall else [1,2,3],
                 nearest=list(range(3)) if wall else [1,2,3],
                 opposing=[0]*3, reversal=[0]*3, area=[1/6]*9, shape=shape)
        lhs, rhs = [0.]*(width*width), [0.]*width
        if stage < 3:
            x['density'] = [2]*3
        if stage in (0,2,3,5):
            x['boundary_velocity'] = ([1,0,0] if stage in (0,3) else [0,0,0])*3
        if stage in (1,3,4,5):
            x['velocity'] = [2,0,0]*(3 if wall else 4)
        if stage in (1,3,4):
            x['gradient'] = [-1,-1,-1,1,0,0,0,1,0,0,0,1]*3
        if stage in (1,3):
            x['bc_multiplier'] = [1,0,0,0]
        if stage == 1:
            x.update(pressure=[0]*4, pressure_gradient=[0]*12, influence_lhs=[2]*9, influence_rhs=[2]*9)
        if stage in (3,4):
            x.update(viscosity=[0]*3, stored_flux=[3]*3)
        if wall:
            x['wall_coefficient'] = [6]*3
        for sample in range(3):
            n = sample+1
            if stage == 0: rhs[n] = -1/3
            if stage == 1:
                rhs[n] = -2/3
                lhs[4*n] = 2
            if stage == 3: rhs[3*n] = -3
            if stage == 4:
                rhs[3*n] = -6
                for i in range(3): lhs[(3*n+i)*12+3*n+i] = 3
            if wall:
                rhs[3*sample:3*sample+3] = [-8,4,4]
                for i in range(3):
                    for f in range(3):
                        for j in range(3):
                            lhs[(3*sample+i)*9+3*f+j] = (4 if i == j else -2)*shape[3*sample+f]
        rows.append(dict(kind='face', id=21+stage%3, nodes=[20,30,40] if wall else [10,20,30,40],
                         components=1 if stage < 3 else 3, inputs=x, outputs=dict(lhs=lhs,rhs=rhs)))
    return rows


def write_fixture(path, ranks=1):
    for stage, row in enumerate(fixture_rows()):
        for rank in range(ranks):
            h = dict(kind='header', schema=1, producer='openaccel', fixture='public_channel',
                     reference_revision=CONTRACT['reference']['revision'],
                     solver_revision=CONTRACT['reference']['solver_gitlink']['revision'],
                     stage=STAGES[stage], call=1, rank=rank, ranks=ranks)
            values = [h]+([row] if rank == 0 else [])+[dict(kind='end',records=1 if rank == 0 else 0)]
            (path/(STAGES[stage]+'.rank'+str(rank)+'.jsonl')).write_text(''.join(json.dumps(v)+'\n' for v in values))


class BoundaryCaptureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.data = self.root/'boundary'
        self.data.mkdir()
        write_fixture(self.data)

    def tearDown(self):
        self.tmp.cleanup()

    def mutate(self, stage, action):
        path = self.data/(STAGES[stage]+'.rank0.jsonl')
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        action(rows)
        path.write_text(''.join(json.dumps(v)+'\n' for v in rows))

    def test_pack_and_no_overwrite(self):
        pack_boundary(self.data,self.root/'inputs.txt')
        self.assertTrue((self.root/'inputs.txt').read_text().startswith('MARS_PUBLIC_BOUNDARY_REPLAY_V1 6\n'))
        self.assertEqual(json.loads((self.root/'inputs.json').read_text())['records'],6)
        with self.assertRaises(ValueError): pack_boundary(self.data,self.root/'inputs.txt')

    def test_empty_rank_and_missing_rank(self):
        write_fixture(self.data,ranks=2)
        self.assertEqual(len(load_boundary(self.data)[0]),6)
        (self.data/(STAGES[0]+'.rank1.jsonl')).unlink()
        with self.assertRaisesRegex(ValueError,'missing boundary rank'): load_boundary(self.data)

    def test_rejected_mutations(self):
        cases = [
            (0,lambda r:r[0].update(reference_revision='wrong')),
            (0,lambda r:r[0].update(fixture='private')),
            (0,lambda r:r[1]['inputs']['density'].__setitem__(0,0)),
            (0,lambda r:r[1]['inputs']['shape'].__setitem__(0,.7)),
            (0,lambda r:r[1]['inputs']['nearest'].__setitem__(0,0)),
            (0,lambda r:r[1]['inputs']['opposing'].__setitem__(0,1)),
            (0,lambda r:r[1]['inputs']['reversal'].__setitem__(0,2)),
            (0,lambda r:r[1]['inputs']['area'].__setitem__(slice(0,3),[0,0,0])),
            (1,lambda r:r[1]['inputs']['bc_multiplier'].__setitem__(0,0)),
            (1,lambda r:r[1]['inputs']['influence_rhs'].pop()),
            (3,lambda r:r[1]['inputs']['viscosity'].__setitem__(0,-1)),
            (3,lambda r:r[1]['nodes'].__setitem__(1,99)),
            (4,lambda r:r[1]['inputs']['area'].__setitem__(0,-1/6)),
            (5,lambda r:r[1]['inputs']['face_nodes'].__setitem__(slice(None),[1,2,0])),
            (5,lambda r:r[1]['inputs']['wall_coefficient'].__setitem__(0,-1)),
            (5,lambda r:r[1]['outputs']['lhs'].__setitem__(0,float('nan'))),
            (5,lambda r:r[-1].update(records=0)),
        ]
        for stage, action in cases:
            with self.subTest(stage=stage,action=action):
                write_fixture(self.data)
                self.mutate(stage,action)
                with self.assertRaises(ValueError): load_boundary(self.data)

    def test_duplicate_owned_face(self):
        def duplicate(rows):
            rows.insert(2,copy.deepcopy(rows[1])); rows[-1]['records']=2
        self.mutate(0,duplicate)
        with self.assertRaisesRegex(ValueError,'duplicate owned'): load_boundary(self.data)

    def test_missing_stage(self):
        (self.data/(STAGES[0]+'.rank0.jsonl')).unlink()
        with self.assertRaisesRegex(ValueError,'missing boundary stages'): load_boundary(self.data)

    def test_global_ids_allow_parent_permutation(self):
        def permute(rows):
            r=rows[1];r['nodes']=[20,30,40,10]
            r['inputs'].update(face_nodes=[0,1,2],nearest=[0,1,2],opposing=[3]*3,bc_multiplier=[0,0,0,1])
        self.mutate(3,permute)
        self.assertEqual(len(load_boundary(self.data)[0]),6)

    def test_finite_wrong_output_needs_replay(self):
        self.mutate(1,lambda r:r[1]['outputs']['lhs'].__setitem__(4,0))
        self.assertEqual(len(load_boundary(self.data)[0]),6)


if __name__ == '__main__':
    unittest.main()
