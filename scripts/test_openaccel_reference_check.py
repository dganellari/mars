#!/usr/bin/env python3
"""GPT/Codex: transport/identity tests, never presented as OpenAccel results."""
import copy
from itertools import combinations
import json
from pathlib import Path
import tempfile
import unittest

import openaccel_reference_check as checker


def records(stage, nodes=(11, 31, 21, 41), reverse=False):
    components = checker.STAGES[stage]
    dofs = [(node, c) for node in nodes for c in range(components)]
    samples, scatter = [], {node: 0. for node in nodes}
    for i, left in enumerate(sorted(nodes)):
        for right in sorted(nodes)[i+1:]:
            flux = (left-right)/100
            scatter[left] -= flux
            scatter[right] += flux
            a, b, sign = (right, left, -1) if reverse else (left, right, 1)
            samples.append(dict(kind="sample", parent=91, left=a, right=b,
                                flux=sign*flux, area=[sign*0.25, sign*0.5, sign*0.75]))
    block = dict(kind="block", parent=91, nodes=list(nodes), components=components,
                 lhs=[r[0]+0.01*c[0]+0.1*r[1]+0.001*c[1] for r in dofs for c in dofs],
                 rhs=[scatter[r[0]] if components == 1 else r[0]+r[1]/10 for r in dofs])
    return samples+[block]


def input_record(stage, nodes, reverse):
    edges = list(combinations(sorted(nodes), 2))
    if reverse:
        edges = [(b, a) for a, b in reversed(edges)]
    fields = {}
    for name, (association, width) in checker.input_layout(stage).items():
        if association == "node":
            values = [node+c/10 for node in nodes for c in range(width)]
        elif association == "sample_node":
            values = [min(a,b)+max(a,b)/100+node/10000+c/100000
                      for a,b in edges for node in nodes for c in range(width)]
        elif association == "oriented_sample":
            values = [(a-b)/100 for a,b in edges]
        else:
            values = [0.]
        if name in ("force", "original_force", "mesh_velocity"):
            values = [0.]*len(values)
        fields[name] = values
    return dict(kind="inputs", parent=91, nodes=list(nodes),
                edges=[node for pair in edges for node in pair], fields=fields)


def write_dump(directory, ranks=1, owner=0, nodes=(11, 31, 21, 41), reverse=False, inputs=False):
    directory.mkdir()
    for stage in checker.STAGES:
        for rank in range(ranks):
            h = dict(kind="header", schema=2 if inputs else 1, precision="float64", coverage="local-interior-only",
                     reference_revision=checker.CONTRACT["reference"]["revision"],
                     solver_revision=checker.CONTRACT["reference"]["solver_gitlink"]["revision"],
                     fixture="unit_tet", producer="harness-test", stage=stage, call=1, rank=rank, ranks=ranks)
            body = records(stage, nodes, reverse) if rank == owner else []
            if inputs and rank == owner:
                body.append(input_record(stage, nodes, reverse))
            rows = [h]+body+[dict(kind="end", records=len(body))]
            (directory/f"{stage}.{rank}.jsonl").write_text("".join(json.dumps(r)+"\n" for r in rows))


class ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ref = self.root/"ref"
        write_dump(self.ref)

    def load(self, path=None):
        return checker.load_dump(path or self.ref, allow_test=True)

    def mutate(self, operation, stage="momentum.interior"):
        path = self.ref/f"{stage}.0.jsonl"
        rows = [json.loads(x) for x in path.read_text().splitlines()]
        operation(rows)
        path.write_text("".join(json.dumps(r)+"\n" for r in rows))

    def test_permutation_orientation_and_empty_ranks(self):
        other = self.root/"other"
        write_dump(other, ranks=4, owner=2, nodes=(41, 21, 31, 11), reverse=True)
        self.assertEqual(checker.compare(self.load(), self.load(other)), 224)

    def test_transposed_nonsymmetric_block_rejected(self):
        ref = self.load()
        other = copy.deepcopy(ref)
        block = other["blocks"][("momentum.interior", 1, 91)]
        block["lhs"] = {(col, row): value for (row, col), value in block["lhs"].items()}
        with self.assertRaisesRegex(ValueError, "comparison failed"):
            checker.compare(ref, other)

    def test_equal_total_flux_does_not_hide_sample_error(self):
        ref = self.load()
        other = copy.deepcopy(ref)
        keys = [k for k in ref["samples"] if k[0] == "momentum.interior"]
        other["samples"][keys[0]]["flux"] += 0.2
        other["samples"][keys[1]]["flux"] -= 0.2
        with self.assertRaisesRegex(ValueError, "comparison failed"):
            checker.compare(ref, other)

    def test_roundoff_tolerance(self):
        ref = self.load()
        other = copy.deepcopy(ref)
        block = other["blocks"][("momentum.interior", 1, 91)]["lhs"]
        block[next(iter(block))] += 1e-13
        checker.compare(ref, other)

    def test_pressure_scatter_rejected(self):
        self.mutate(lambda r: r[-2]["rhs"].__setitem__(0, 999), "pressure.interior")
        with self.assertRaisesRegex(ValueError, "sample scatter"):
            self.load()

    def test_missing_rank(self):
        self.mutate(lambda r: r[0].__setitem__("ranks", 4))
        with self.assertRaisesRegex(ValueError, "missing rank"):
            self.load()

    def test_missing_stage(self):
        (self.ref/"pressure.interior.0.jsonl").unlink()
        with self.assertRaisesRegex(ValueError, "both interior stages"):
            self.load()

    def test_truncated_file(self):
        self.mutate(lambda r: r.pop())
        with self.assertRaisesRegex(ValueError, "incomplete export"):
            self.load()

    def test_missing_sample(self):
        def remove(rows):
            rows.pop(1)
            rows[-1]["records"] -= 1
        self.mutate(remove)
        with self.assertRaisesRegex(ValueError, "sample coverage"):
            self.load()

    def test_duplicate_owned_block(self):
        def duplicate(rows):
            rows.insert(-1, copy.deepcopy(rows[-2]))
            rows[-1]["records"] += 1
        self.mutate(duplicate)
        with self.assertRaisesRegex(ValueError, "duplicate owned block"):
            self.load()

    def test_nonfinite(self):
        self.mutate(lambda r: r[-2]["lhs"].__setitem__(0, float('nan')))
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            self.load()

    def test_wrong_source(self):
        self.mutate(lambda r: r[0].__setitem__("reference_revision", "wrong"))
        with self.assertRaisesRegex(ValueError, "wrong source pin"):
            self.load()

    def test_bad_dimensions(self):
        self.mutate(lambda r: r[-2]["lhs"].pop())
        with self.assertRaisesRegex(ValueError, "dimensions"):
            self.load()

    def test_test_producer_cannot_certify_run(self):
        with self.assertRaisesRegex(ValueError, "cannot certify"):
            checker.load_dump(self.ref)

    def test_duplicate_json_key(self):
        path = self.ref/"momentum.interior.0.jsonl"
        path.write_text(path.read_text().replace('"call": 1', '"call": 1, "call": 1'))
        with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
            self.load()

    def frozen_dump(self):
        result = self.root/"frozen"
        write_dump(result, inputs=True)
        self.ref = result
        return self.load()

    def test_input_node_sample_permutation_and_empty_rank(self):
        ref = self.frozen_dump()
        other = self.root/"other"
        write_dump(other, ranks=4, owner=2, nodes=(41,21,31,11), reverse=True, inputs=True)
        self.assertGreater(checker.compare(ref, self.load(other)), 224)

    def test_frozen_input_mutation_with_identical_matrix_rejected(self):
        ref = self.frozen_dump()
        self.mutate(lambda r: r[-2]["fields"]["velocity"].__setitem__(0, 999))
        with self.assertRaisesRegex(ValueError, "input velocity"):
            checker.compare(ref, self.load())

    def test_frozen_capture_cannot_accept_old_binary(self):
        with self.assertRaisesRegex(ValueError, "frozen inputs are required"):
            checker.load_dump(self.ref, allow_test=True, require_inputs=True)

    def test_missing_input_block_rejected(self):
        self.frozen_dump()
        def remove(rows):
            rows.pop(-2)
            rows[-1]["records"] -= 1
        self.mutate(remove)
        with self.assertRaisesRegex(ValueError, "missing/extra frozen input blocks"):
            self.load()

    def test_wrong_input_shape_rejected(self):
        self.frozen_dump()
        self.mutate(lambda r: r[-2]["fields"]["shape_gradient"].pop())
        with self.assertRaisesRegex(ValueError, "dimensions"):
            self.load()

    def test_input_flux_and_recorded_flux_must_match(self):
        self.frozen_dump()
        self.mutate(lambda r: r[-2]["fields"]["stored_flux"].__setitem__(0, 999))
        with self.assertRaisesRegex(ValueError, "stored input flux differs"):
            self.load()

    def test_unsupported_physics_rejected(self):
        self.frozen_dump()
        self.mutate(lambda r: r[-2]["fields"]["nso"].__setitem__(0, 1))
        with self.assertRaisesRegex(ValueError, "unsupported frozen input: nso"):
            self.load()


if __name__ == "__main__":
    unittest.main()
