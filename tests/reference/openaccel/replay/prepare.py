#!/usr/bin/env python3
"""Pack validated public captures for the standalone MARS host/CUDA gate."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if not (HERE / 'openaccel_reference_check.py').is_file():
    sys.path.insert(0, str(HERE.parents[3] / 'scripts'))
from openaccel_reference_check import load_dump, unique_object

FIELDS = [('coordinates', 12), ('velocity', 12), ('density', 4),
          ('velocity_shape', 24), ('coordinate_shape', 24), ('shape_gradient', 72),
          ('area', 18), ('viscosity', 4), ('velocity_blend', 12),
          ('velocity_gradient', 36), ('stored_flux', 6), ('pressure', 4),
          ('pressure_gradient', 12), ('influence_lhs', 12), ('influence_rhs', 12),
          ('density_blend', 4), ('density_gradient', 12)]


def prepare(exports, output):
    data = load_dump(exports, require_inputs=True)
    if data['producer'] != 'openaccel' or data['signature'][0] != 'public_channel':
        raise ValueError('only the named public OpenAccel channel is permitted')
    if output.exists():
        raise ValueError('output exists; preserve it and choose a fresh path')
    with io.StringIO() as out:
        out.write('MARS_PUBLIC_TET_REPLAY_V1 {}\n'.format(len(data['inputs'])))
        for path in sorted(exports.glob('*.jsonl')):
            rows = [json.loads(line, object_pairs_hook=unique_object) for line in path.read_text().splitlines()]
            stage, call = rows[0]['stage'], rows[0]['call']
            for row in rows:
                if row['kind'] != 'inputs':
                    continue
                f = row['fields']
                if stage == 'pressure.interior' and f['harmonic_gradient'] != [0]:
                    raise ValueError('harmonic reconstruction is outside this incompressible replay')
                nodes, edges = row['nodes'], row['edges']
                key = (stage, call, row['parent'])
                block = data['blocks'][key]
                components = 3 if stage == 'momentum.interior' else 1
                out.write('{} {} {}\n'.format(int(components == 3), call, row['parent']))
                out.write(' '.join(str(nodes.index(n)) for n in edges)+'\n')
                areas, fluxes = [], []
                for i in range(0, 12, 2):
                    left, right = edges[i:i+2]
                    sample = data['samples'][key+tuple(sorted((left, right)))]
                    sign = 1 if left < right else -1
                    areas.extend(sign*x for x in sample['area'])
                    fluxes.append(sign*sample['flux'])
                for name, width in FIELDS:
                    values = areas if name == 'area' else f.get(name, [0]*width)
                    if len(values) != width:
                        raise ValueError('wrong field width: '+name)
                    out.write(' '.join(format(x, '.17g') for x in values)+'\n')
                dofs = [(n, c) for n in nodes for c in range(components)]
                lhs = [block['lhs'][(a, b)] for a in dofs for b in dofs]
                rhs = [block['rhs'][dof] for dof in dofs]
                for values, size in ((lhs, 144), (rhs, 12), (fluxes, 6)):
                    out.write(' '.join(format(x, '.17g') for x in values+[0]*(size-len(values)))+'\n')
        packed = out.getvalue()
    with output.open('x') as out:
        out.write(packed)
    record = {'author': 'GPT/Codex', 'scope': 'public frozen interior assembly only',
              'source_sha256': data['hashes'], 'blocks': len(data['inputs']),
              'packed_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Prepared {} public blocks: {}'.format(len(data['inputs']), output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exports', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    try:
        prepare(args.exports, args.output)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, 'ERROR: '+str(error)+'\n')
