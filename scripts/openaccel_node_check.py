#!/usr/bin/env python3
"""Validate and pack public steady node/relaxation exports without solver replicas."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys

from openaccel_reference_check import CONTRACT, integer, numbers, require, unique_object

STAGES = ('momentum.node', 'momentum.relaxation', 'momentum.influence', 'momentum.boundary_relaxation')
INPUTS = (
    dict(density=1, volume=1, pseudo_dt=1, mass_divergence=1, velocity=3,
         pressure_gradient=3, force=3, source=3, coriolis=9),
    dict(lhs=9, alpha=1),
    dict(volume=1, row_blocks=None, diagonal_block=1, consistent=1, fractional_step=1, transient=1, small=1),
    dict(rhs=3, factor=1),
)
OUTPUTS = (dict(lhs=9, rhs=3), dict(lhs=9), dict(d=3, d_tilde=3), dict(rhs=3))


def load_nodes(directory):
    records, headers, hashes = {}, {}, {}
    for path in sorted(directory.glob('*.jsonl')):
        content = path.read_bytes()
        rows = [json.loads(line, object_pairs_hook=unique_object) for line in content.splitlines()]
        require(len(rows) >= 2, 'empty/truncated node file')
        head, tail = rows[0], rows[-1]
        require(head.get('kind') == 'header' and head.get('schema') == 1, 'invalid node header')
        require(head.get('producer') == 'openaccel' and head.get('fixture') == 'public_channel', 'not a public reference')
        require(head.get('reference_revision') == CONTRACT['reference']['revision']
                and head.get('solver_revision') == CONTRACT['reference']['solver_gitlink']['revision'], 'wrong source pin')
        require(head.get('stage') in STAGES, 'unknown node stage')
        stage = STAGES.index(head['stage'])
        call, rank, ranks = integer(head['call'], 'call'), integer(head['rank'], 'rank', 0), integer(head['ranks'], 'ranks')
        require(rank < ranks, 'rank outside communicator')
        key = (stage, call, rank)
        require(key not in headers, 'duplicate node header')
        headers[key] = ranks
        require(tail.get('kind') == 'end' and type(tail.get('records')) is int
                and tail['records'] == len(rows)-2, 'invalid node footer')
        hashes[path.name] = hashlib.sha256(content).hexdigest()
        for row in rows[1:-1]:
            require(row.get('kind') == 'node', 'invalid node record')
            node = integer(row['id'], 'global node')
            key = (stage, call, node)
            require(key not in records, 'duplicate owned node (including repeated boundary relaxation)')
            for group, layout in (('inputs', INPUTS[stage]), ('outputs', OUTPUTS[stage])):
                fields = row[group]
                require(isinstance(fields, dict) and fields.keys() == layout.keys(), 'missing/extra node fields')
                for name, width in layout.items():
                    if width is None:
                        require(isinstance(fields[name], list) and len(fields[name]) > 0
                                and len(fields[name]) % 9 == 0, 'invalid block row')
                    numbers(fields[name], width if width is not None else len(fields[name]), name)
            x = row['inputs']
            if stage == 0:
                require(x['density'][0] > 0 and x['volume'][0] > 0 and x['pseudo_dt'][0] > 0, 'invalid steady node scales')
                require(all(v == 0 for v in x['coriolis']), 'rotating frame not supported')
            elif stage == 1:
                require(0 < x['alpha'][0] <= 1, 'invalid momentum relaxation')
            elif stage == 2:
                require(x['volume'][0] > 0 and x['fractional_step'] == [0] and x['transient'] == [0], 'requires steady SIMPLE/SIMPLEC')
                require(x['consistent'][0] in (0, 1), 'invalid consistent flag')
                require(x['small'] == [sys.float_info.epsilon], 'reference SMALL differs')
                d = x['diagonal_block'][0]
                require(d == int(d) and 0 <= d < len(x['row_blocks'])//9, 'invalid diagonal block')
            else:
                require(0 < x['factor'][0] <= 1, 'invalid boundary relaxation')
            records[key] = dict(row, rank=rank)
    require(headers and set(s for s, c, r in headers) == set(range(4)), 'missing node stages')
    require(len(set(headers.values())) == 1, 'communicator changed between stages')
    calls = {s: {c for a, c, r in headers if a == s} for s in range(4)}
    require(all(calls[s] == calls[0] for s in range(4)), 'node stage call counts differ')
    require(calls[0] == set(range(1, max(calls[0])+1)), 'missing node call')
    for s in range(4):
        for c in calls[0]:
            ranks = {r for a, b, r in headers if (a, b) == (s, c)}
            require(ranks == set(range(next(iter(headers.values())))), 'missing rank file, including empty rank')
    for c in calls[0]:
        sets = [{n for s, call, n in records if (s, call) == (stage, c)} for stage in range(4)]
        require(sets[0] and sets[0] == sets[1] == sets[2], 'node/relaxation/influence coverage differs')
        require(sets[3] <= sets[0], 'boundary node outside domain')
        for n in sets[3]:
            require(records[(3, c, n)]['rank'] == records[(0, c, n)]['rank'], 'boundary ownership changed')
        for n in sets[0]:
            require(len({records[(s, c, n)]['rank'] for s in range(3)}) == 1, 'node ownership changed')
            # Influence must consume the actual relaxed diagonal, not an independently rescaled one.
            relaxed = records[(1, c, n)]['outputs']['lhs']
            influence = records[(2, c, n)]['inputs']
            offset = 9*int(influence['diagonal_block'][0])
            require(relaxed == influence['row_blocks'][offset:offset+9], 'influence uses a different diagonal stage')
    return records, hashes


def pack_nodes(directory, output):
    records, hashes = load_nodes(directory)
    require(not output.exists() and not output.with_suffix('.json').exists(), 'output/provenance already exists')
    with io.StringIO() as stream:
        stream.write('MARS_PUBLIC_NODE_REPLAY_V1 {}\n'.format(len(records)))
        for (stage, call, node), row in sorted(records.items()):
            x, expected = row['inputs'], row['outputs']
            stream.write('{} {} {}\n'.format(stage, call, node))
            if stage == 0:
                names = ['density', 'volume', 'pseudo_dt', 'mass_divergence', 'velocity', 'pressure_gradient', 'force', 'source']
                values = [v for name in names for v in x[name]]
            elif stage == 1:
                values = x['alpha']+x['lhs']
            elif stage == 2:
                values = [len(x['row_blocks'])//9, int(x['diagonal_block'][0]), int(x['consistent'][0])]+x['volume']+x['row_blocks']
            else:
                values = x['factor']+x['rhs']
            stream.write(' '.join(format(v, '.17g') for v in values)+'\n')
            expected_values = expected.get('lhs', [0]*9)+expected.get('rhs', [0]*3)+expected.get('d', [0]*3)+expected.get('d_tilde', [0]*3)
            stream.write(' '.join(format(v, '.17g') for v in expected_values)+'\n')
        packed = stream.getvalue()
    with output.open('x') as stream:
        stream.write(packed)
    output.with_suffix('.json').write_text(json.dumps(dict(source_sha256=hashes, records=len(records),
        packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(), coverage='steady-node-and-postassembly'), indent=2)+'\n')
    print('Prepared {} public node records: {}'.format(len(records), output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exports', type=Path)
    parser.add_argument('--pack', type=Path)
    args = parser.parse_args()
    try:
        if args.pack:
            pack_nodes(args.exports, args.pack)
        else:
            records, hashes = load_nodes(args.exports)
            print('PASS: node capture coverage/identity: {} records, {} files; numerical replay pending'.format(len(records), len(hashes)))
    except (ValueError, KeyError, TypeError, OSError, OverflowError) as error:
        parser.exit(1, 'ERROR: '+str(error)+'\n')
