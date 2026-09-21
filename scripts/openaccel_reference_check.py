#!/usr/bin/env python3
"""Compare public local-interior exports by global IDs, independently of rank/order.

GPT/Codex, 2026-09-12. This is not the full OpenAccel numerical-contract gate.
"""
import argparse
from collections import defaultdict
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = Path(__file__).with_name("contract_v1.json")
if not CONTRACT_PATH.is_file():
    CONTRACT_PATH = ROOT / "tests/data/public_openaccel_reference/contract_v1.json"
CONTRACT = json.loads(CONTRACT_PATH.read_text())
STAGES = {"momentum.interior": 3, "pressure.interior": 1}


def input_layout(stage):
    layout = {name: ("node", width) for name, width in
              (("coordinates", 3), ("velocity", 3), ("density", 1))}
    layout.update({"velocity_shape": ("sample_node", 1), "coordinate_shape": ("sample_node", 1),
                   "shape_gradient": ("sample_node", 3)})
    controls = ["compressible", "velocity_shifted", "gradient_shifted"]
    if stage == "momentum.interior":
        layout.update({"viscosity": ("node", 1), "velocity_blend": ("node", 3),
                       "velocity_gradient": ("node", 9), "stored_flux": ("oriented_sample", 1)})
        controls += ["nso", "nso_fourth_factor"]
    else:
        layout.update({name: ("node", width) for name, width in
                       (("pressure", 1), ("pressure_gradient", 3), ("influence_lhs", 3),
                        ("influence_rhs", 3), ("density_blend", 1), ("density_gradient", 3),
                        ("force", 3), ("original_force", 3), ("mesh_velocity", 3))})
        controls += ["consistent", "harmonic_gradient", "mesh_moving", "frame_rotating"]
    layout.update({name: ("element", 1) for name in controls})
    return layout


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, name, minimum=1):
    require(type(value) is int and value >= minimum, f"invalid {name}: {value!r}")
    return value


def numbers(value, length, name):
    require(isinstance(value, list) and len(value) == length, f"wrong {name} dimensions")
    require(all(type(x) in (int, float) and math.isfinite(x) for x in value), f"nonfinite/invalid {name}")
    return value


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def parse_inputs(record, stage, rank):
    nodes, edges, fields = record["nodes"], record["edges"], record["fields"]
    require(isinstance(nodes, list) and len(nodes) == 4, "expected four frozen-input nodes")
    for node in nodes:
        integer(node, "input node ID")
    require(len(set(nodes)) == 4, "repeated frozen-input node")
    require(isinstance(edges, list) and len(edges) == 12, "expected six input edges")
    for node in edges:
        integer(node, "input edge node ID")
        require(node in nodes, "input edge endpoint absent from element")
    pairs = [tuple(sorted(edges[i:i+2])) for i in range(0, 12, 2)]
    require(set(pairs) == set(combinations(sorted(nodes), 2)), "invalid input edge coverage")
    layout = input_layout(stage)
    require(isinstance(fields, dict) and fields.keys() == layout.keys(), "missing/extra frozen input fields")
    mapped = {}
    for name, (association, width) in layout.items():
        if association == "node":
            keys = [(node, c) for node in nodes for c in range(width)]
        elif association == "sample_node":
            keys = [(pair, node, c) for pair in pairs for node in nodes for c in range(width)]
        elif association == "oriented_sample":
            keys = [(pair, c) for pair in pairs for c in range(width)]
        else:
            keys = list(range(width))
        values = numbers(fields[name], len(keys), "input " + name)
        if association == "oriented_sample":
            values = [value * (1 if edges[2*i] < edges[2*i+1] else -1)
                      for i, value in enumerate(values)]
        if association == "element" and name != "nso_fourth_factor":
            require(values[0] in (0, 1), "invalid input switch: " + name)
        mapped[name] = dict(zip(keys, values))
    # These omitted physics need additional frozen data before a replay is meaningful.
    unsupported = ["compressible", "nso"] if stage == "momentum.interior" else [
        "compressible", "mesh_moving", "frame_rotating", "force", "original_force", "mesh_velocity"]
    for name in unsupported:
        require(all(value == 0 for value in mapped[name].values()), "unsupported frozen input: " + name)
    return {"nodes": frozenset(nodes), "rank": rank, "fields": mapped}


def load_dump(directory, *, allow_test=False, require_inputs=False):
    files = sorted(directory.glob("*.jsonl"))
    require(files, f"no exports in {directory}")
    headers, blocks, samples, hashes, inputs = {}, {}, {}, {}, {}
    schemas = set()
    signature = None
    producers = set()
    for path in files:
        raw = path.read_bytes()
        hashes[path.name] = hashlib.sha256(raw).hexdigest()
        rows = [json.loads(line, object_pairs_hook=unique_object) for line in raw.splitlines()]
        require(all(isinstance(row, dict) for row in rows), "export records must be JSON objects")
        require(len(rows) >= 2 and rows[0].get("kind") == "header"
                and rows[-1].get("kind") == "end", f"incomplete export: {path.name}")
        h = rows[0]
        require(type(h["schema"]) is int and h["schema"] in (1, 2) and h["precision"] == "float64"
                and h["coverage"] == "local-interior-only", "unsupported export schema/coverage")
        schemas.add(h["schema"])
        ref = CONTRACT["reference"]
        require(h["reference_revision"] == ref["revision"]
                and h["solver_revision"] == ref["solver_gitlink"]["revision"], "wrong source pin")
        require(h["fixture"] in {f["id"] for f in CONTRACT["fixtures"]}, "unknown public fixture")
        require(h["producer"] in ({"openaccel", "mars", "harness-test"} if allow_test
                                  else {"openaccel", "mars"}), "test/unknown producer cannot certify a run")
        producers.add(h["producer"])
        current = (h["fixture"], h["reference_revision"], h["solver_revision"])
        require(signature in (None, current), "mixed fixtures/revisions")
        signature = current
        stage = h["stage"]
        require(stage in STAGES, "unknown stage")
        call = integer(h["call"], "call")
        rank = integer(h["rank"], "rank", 0)
        ranks = integer(h["ranks"], "ranks")
        require(rank < ranks, "rank outside communicator")
        key = (stage, call, rank)
        require(key not in headers, f"duplicate rank/stage/call: {key}")
        headers[key] = ranks
        require(type(rows[-1].get("records")) is int and rows[-1]["records"] == len(rows)-2,
                "wrong export footer count")
        for record in rows[1:-1]:
            parent = integer(record["parent"], "parent")
            block_key = (stage, call, parent)
            if record["kind"] == "block":
                require(block_key not in blocks, f"duplicate owned block: {block_key}")
                nodes = record["nodes"]
                require(isinstance(nodes, list) and len(nodes) == 4, "expected Tet4 node IDs")
                for node in nodes:
                    integer(node, "node ID")
                require(len(set(nodes)) == 4, "repeated node ID")
                components = integer(record["components"], "components")
                require(components == STAGES[stage], "wrong component count for stage")
                dofs = [(node, c) for node in nodes for c in range(components)]
                n = len(dofs)
                lhs = numbers(record["lhs"], n*n, "lhs")
                rhs = numbers(record["rhs"], n, "rhs")
                blocks[block_key] = {
                    "nodes": frozenset(nodes), "rank": rank,
                    "lhs": {(row, col): lhs[i*n+j] for i, row in enumerate(dofs)
                            for j, col in enumerate(dofs)},
                    "rhs": {row: rhs[i] for i, row in enumerate(dofs)},
                }
            elif record["kind"] == "sample":
                left, right = integer(record["left"], "left"), integer(record["right"], "right")
                require(left != right, "sample connects a node to itself")
                flux = numbers([record["flux"]], 1, "flux")[0]
                area = numbers(record["area"], 3, "area")
                sign = 1 if left < right else -1
                pair = tuple(sorted((left, right)))
                key = (*block_key, *pair)
                require(key not in samples, f"duplicate owned sample: {key}")
                samples[key] = {"rank": rank, "flux": sign*flux,
                                "area": [sign*x for x in area]}
            elif record["kind"] == "inputs":
                require(h["schema"] == 2, "frozen inputs require schema 2")
                require(block_key not in inputs, f"duplicate frozen inputs: {block_key}")
                inputs[block_key] = parse_inputs(record, stage, rank)
            else:
                raise ValueError(f"unexpected record kind: {record['kind']}")
    require(len(schemas) == 1, "mixed input coverage schemas")
    require(not require_inputs or schemas == {2}, "frozen inputs are required")
    if schemas == {2}:
        require(inputs.keys() == blocks.keys(), "missing/extra frozen input blocks")
        for key, data in inputs.items():
            require(data["rank"] == blocks[key]["rank"] and data["nodes"] == blocks[key]["nodes"],
                    "frozen input ownership/connectivity differs from block")
    require(len(producers) == 1, "mixed producers")
    require({h[0] for h in headers} == set(STAGES), "both interior stages are required")
    calls = {stage: {h[1] for h in headers if h[0] == stage} for stage in STAGES}
    require(calls["momentum.interior"] == calls["pressure.interior"], "stage call counts differ")
    for stage, indices in calls.items():
        require(indices == set(range(1, max(indices)+1)), "missing assembly call")
        for call in indices:
            entries = {rank: size for (s, c, rank), size in headers.items() if (s, c) == (stage, call)}
            require(len(set(entries.values())) == 1, "inconsistent communicator size")
            require(set(entries) == set(range(next(iter(entries.values())))), "missing rank output")
            require(any(k[:2] == (stage, call) for k in blocks), "globally empty assembly call")
    require(len(set(headers.values())) == 1, "communicator size changed within run")
    for call in calls["momentum.interior"]:
        parents = [{k[2] for k in blocks if k[:2] == (stage, call)} for stage in STAGES]
        require(parents[0] == parents[1], "momentum/pressure element coverage differs")
    by_block = defaultdict(set)
    for key in samples:
        require(key[:3] in blocks, "sample has no local block")
        by_block[key[:3]].add(key)
    for key, block in blocks.items():
        expected = {(*key, *pair) for pair in combinations(sorted(block["nodes"]), 2)}
        actual = by_block[key]
        require(actual == expected, f"wrong Tet4 sample coverage: {key}")
        require(all(samples[s]["rank"] == block["rank"] for s in actual), "split block/sample ownership")
        if key in inputs and key[0] == "momentum.interior":
            stored = inputs[key]["fields"]["stored_flux"]
            require(all(stored[(s[3:], 0)] == samples[s]["flux"] for s in actual),
                    "stored input flux differs from assembly sample")
        if key[0] == "pressure.interior":
            scatter = {node: 0.0 for node in block["nodes"]}
            for s in sorted(actual):
                scatter[s[3]] -= samples[s]["flux"]
                scatter[s[4]] += samples[s]["flux"]
            scale = max(1., *(abs(samples[s]["flux"]) for s in actual),
                        *(abs(x) for x in block["rhs"].values()))
            require(all(abs(scatter[node]-block["rhs"][(node, 0)]) <= 1e-12*scale
                        for node in scatter), f"pressure RHS differs from actual sample scatter: {key}")
    return {"signature": signature, "producer": next(iter(producers)), "blocks": blocks,
            "samples": samples, "inputs": inputs, "hashes": hashes}


def compare(reference, candidate):
    require(reference["signature"] == candidate["signature"], "fixture/source mismatch")
    failures = []
    count = 0

    def values(a, b, label):
        nonlocal count
        require(a.keys() == b.keys(), f"different IDs/components: {label}")
        scale = max(1., *(abs(x) for x in a.values()))
        for key, value in a.items():
            count += 1
            if abs(value-b[key]) > 1e-12*scale and len(failures) < 12:
                failures.append(f"{label} {key}: reference={value:.17g}, candidate={b[key]:.17g}")

    for kind in ["blocks", "samples", "inputs"]:
        require(reference[kind].keys() == candidate[kind].keys(), f"missing/extra {kind} IDs")
    for key, block in reference["blocks"].items():
        other = candidate["blocks"][key]
        require(block["nodes"] == other["nodes"], f"changed connectivity: {key}")
        for field in ["lhs", "rhs"]:
            values(block[field], other[field], f"{key} {field}")
    for key, data in reference["inputs"].items():
        other = candidate["inputs"][key]
        for name, field in data["fields"].items():
            values(field, other["fields"][name], f"{key} input {name}")
    # Scale flux and area separately; their physical units differ.
    grouped = defaultdict(list)
    for key in reference["samples"]:
        grouped[key[:3]].append(key)
    for parent, keys in grouped.items():
        for field in ["flux", "area"]:
            a, b = {}, {}
            for key in keys:
                av, bv = reference["samples"][key][field], candidate["samples"][key][field]
                for c, (x, y) in enumerate(zip([av] if field == "flux" else av,
                                              [bv] if field == "flux" else bv)):
                    a[(key[3:], c)], b[(key[3:], c)] = x, y
            values(a, b, f"{parent} {field}")
    require(not failures, "comparison failed:\n" + "\n".join(failures))
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", type=Path, help="new JSON report path, including input hashes")
    args = parser.parse_args()
    try:
        require(args.reference.resolve() != args.candidate.resolve(), "refusing self-comparison")
        ref, got = load_dump(args.reference), load_dump(args.candidate)
        require(ref["producer"] == "openaccel", "reference must come from OpenAccel")
        count = compare(ref, got)
        report = {"status": "PASS", "coverage": "local-interior-only", "values": count,
                  "reference_artifact_sha256": ref["hashes"], "candidate_artifact_sha256": got["hashes"],
                  "frozen_inputs_compared": bool(ref["inputs"]),
                  "full_contract_passed": False}
        if args.report:
            with args.report.open("x") as out:
                json.dump(report, out, indent=2)
                out.write("\n")
        print(f"PASS: {count} local interior values by global ID; "
              f"frozen inputs compared={bool(ref['inputs'])}; full contract pending")
    except (ValueError, KeyError, TypeError, OSError, OverflowError) as error:
        parser.exit(1, f"FAIL: {error}\n")


if __name__ == "__main__":
    main()
