"""Validate public boundary captures and pack native local blocks for GPU replay."""
import argparse
import hashlib
import json
from pathlib import Path

from openaccel_reference_check import CONTRACT, integer, numbers, require, unique_object

STAGES = ('pressure.inlet', 'pressure.outlet', 'pressure.wall', 'momentum.inlet', 'momentum.outlet', 'momentum.wall')
ARRAYS = dict(area=9, shape=9, gradient=36, velocity=12, boundary_velocity=9, viscosity=3,
              density=3, pressure=4, pressure_gradient=12, influence_lhs=9, influence_rhs=9,
              bc_multiplier=4, stored_flux=3, wall_coefficient=3)
MAPS = ('face_nodes', 'nearest', 'opposing', 'reversal')
COMMON = set(MAPS) | {'area', 'shape'}
FIELDS = [COMMON | set(v.split()) for v in (
    'density boundary_velocity',
    'density velocity pressure gradient pressure_gradient influence_lhs influence_rhs bc_multiplier',
    'density boundary_velocity',
    'velocity gradient viscosity stored_flux boundary_velocity bc_multiplier',
    'velocity gradient viscosity stored_flux',
    'velocity boundary_velocity wall_coefficient')]


def load_boundary(directory):
    records, headers, hashes = {}, {}, {}
    for path in sorted(directory.glob('*.jsonl')):
        content = path.read_bytes()
        rows = [json.loads(line, object_pairs_hook=unique_object) for line in content.splitlines()]
        require(len(rows) >= 2, 'truncated boundary export')
        h, tail = rows[0], rows[-1]
        require(h.get('kind') == 'header' and h.get('schema') == 1 and h.get('producer') == 'openaccel'
                and h.get('fixture') == 'public_channel', 'invalid public boundary header')
        require(h.get('reference_revision') == CONTRACT['reference']['revision'] and
                h.get('solver_revision') == CONTRACT['reference']['solver_gitlink']['revision'], 'wrong source pin')
        require(h.get('stage') in STAGES, 'invalid boundary stage')
        stage = STAGES.index(h['stage'])
        call, rank, ranks = integer(h['call'], 'call'), integer(h['rank'], 'rank', 0), integer(h['ranks'], 'ranks')
        require(rank < ranks and (stage, call, rank) not in headers, 'invalid/duplicate boundary rank')
        headers[(stage, call, rank)] = ranks
        require(tail.get('kind') == 'end' and type(tail.get('records')) is int and tail['records'] == len(rows)-2, 'invalid boundary footer')
        hashes[path.name] = hashlib.sha256(content).hexdigest()
        for row in rows[1:-1]:
            require(row.get('kind') == 'face', 'invalid face record')
            face = integer(row['id'], 'face ID')
            key = (stage, call, face)
            require(key not in records, 'duplicate owned face')
            count, components = (3 if stage == 5 else 4), (1 if stage < 3 else 3)
            ids = row['nodes']
            require(len(ids) == count and len(set(ids)) == count, 'invalid connected node IDs')
            for node in ids: integer(node, 'node ID')
            require(row['components'] == components, 'wrong component count')
            x, y = row['inputs'], row['outputs']
            require(x.keys() == FIELDS[stage], 'missing/extra boundary inputs')
            require(y.keys() == {'lhs','rhs'}, 'missing/extra boundary outputs')
            for name, values in x.items():
                width = 3 if name in MAPS else (9 if stage == 5 and name == 'velocity' else ARRAYS[name])
                numbers(values, width, name)
            for name in MAPS:
                require(all(v == int(v) for v in x[name]), 'noninteger boundary map')
            require(set(x['reversal']) <= {0,1}, 'invalid reversal flags')
            require(len(set(x['face_nodes'])) == 3 and all(0 <= n < count for n in x['face_nodes']), 'invalid face map')
            require(set(x['nearest']) == set(x['face_nodes']), 'invalid nearest-node map')
            if stage != 5:
                require(all(0 <= n < 4 and n not in x['face_nodes'] for n in x['opposing']), 'invalid opposite-node map')
            else:
                require(x['face_nodes'] == [0,1,2], 'wall block must use native face ordering')
            require(all(sum(v*v for v in x['area'][3*i:3*i+3]) > 0 for i in range(3)), 'zero face sample area')
            for i in range(3):
                weights=x['shape'][3*i:3*i+3]
                require(all(v >= 0 for v in weights) and abs(sum(weights)-1) < 1e-12, 'invalid interpolation weights')
            if 'density' in x: require(all(v > 0 for v in x['density']), 'invalid density')
            if 'viscosity' in x: require(all(v >= 0 for v in x['viscosity']), 'invalid viscosity')
            if 'wall_coefficient' in x: require(all(v >= 0 for v in x['wall_coefficient']), 'invalid wall coefficient')
            if 'bc_multiplier' in x:
                require(x['bc_multiplier'] == [0 if n in x['face_nodes'] else 1 for n in range(4)], 'incorrect boundary column mask')
            numbers(y['lhs'], (count*components)**2, 'lhs'); numbers(y['rhs'], count*components, 'rhs')
            records[key] = dict(row, rank=rank)
    require(headers and {s for s,c,r in headers} == set(range(6)), 'missing boundary stages')
    require(len(set(headers.values())) == 1, 'communicator changed')
    calls={c for s,c,r in headers}
    require(calls == set(range(1,max(calls)+1)), 'missing boundary call')
    ranks=next(iter(headers.values()))
    require(set(headers) == {(s,c,r) for s in range(6) for c in calls for r in range(ranks)}, 'missing boundary rank/call file')
    for c in calls:
        for s in range(3):
            p={f for st,call,f in records if (st,call)==(s,c)}
            m={f for st,call,f in records if (st,call)==(s+3,c)}
            require(p and p==m, 'pressure/momentum boundary coverage differs')
            for f in p:
                a,b=records[(s,c,f)],records[(s+3,c,f)]
                require(a['rank']==b['rank'], 'face ownership changed')
                def face_ids(row): return {row['nodes'][int(i)] for i in row['inputs']['face_nodes']}
                require(face_ids(a)==face_ids(b), 'face connectivity changed')
                # Compare samples by their nearest global node, allowing local permutations.
                def samples(row):
                    x=row['inputs']; return {row['nodes'][int(n)]:x['area'][3*i:3*i+3] for i,n in enumerate(x['nearest'])}
                require(samples(a)==samples(b), 'boundary area/orientation changed between stages')
    return records,hashes


def pack_boundary(directory, output):
    records, hashes = load_boundary(directory)
    require(not output.exists() and not output.with_suffix('.json').exists(), 'output/provenance exists')
    lines=['MARS_PUBLIC_BOUNDARY_REPLAY_V1 '+str(len(records))]
    for (stage,call,face), row in sorted(records.items()):
        x=row['inputs']; lines.append('{} {} {}'.format(stage,call,face))
        values=[int(v) for name in MAPS for v in x[name]]
        for name,width in ARRAYS.items():
            a=x.get(name,[]); values += a+[0]*(width-len(a))
        lines.append(' '.join(format(v,'.17g') for v in values))
        for name,width in (('lhs',144),('rhs',12)):
            a=row['outputs'][name]; lines.append(' '.join(format(v,'.17g') for v in a+[0]*(width-len(a))))
    with output.open('x') as f: f.write('\n'.join(lines)+'\n')
    output.with_suffix('.json').write_text(json.dumps(dict(records=len(records),source_sha256=hashes,
        packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),coverage='frozen-boundary-blocks'),indent=2)+'\n')
    print('Prepared {} public boundary records: {}'.format(len(records),output))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exports',type=Path); parser.add_argument('--pack',type=Path)
    args=parser.parse_args()
    try:
        if args.pack: pack_boundary(args.exports,args.pack)
        else:
            records,hashes=load_boundary(args.exports)
            print('PASS: boundary capture coverage/identity: {} blocks; numerical replay pending'.format(len(records)))
    except (ValueError,KeyError,TypeError,OSError,OverflowError) as error:
        parser.exit(1,'ERROR: '+str(error)+'\n')
