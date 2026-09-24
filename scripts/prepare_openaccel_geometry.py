"""Pack the existing public channel captures for native geometry/gradient checks."""
import argparse
from collections import defaultdict
import hashlib
import io
import json
from pathlib import Path

from openaccel_reference_check import load_dump, require, unique_object
from openaccel_boundary_check import load_boundary

FACES = ((0,1,3), (1,2,3), (0,3,2), (0,2,1))
EDGES = ((0,1), (1,2), (0,2), (0,3), (1,3), (2,3))
# Fix the interpolation/history/BC profile, not just the mesh's friendly name.
DECK_SHA256 = 'dc2d6be5ca0d2d4845891cef62c3eb2bc12b857f7a72d669365ecfab957452f5'
MESH_SHA256 = '3b3a2d8101291955d496f9c696ce250cb5acc5c2e1e34eae05c9db5ae5f861fd'


def verify_run(directory, hashes, hash_key):
    run = json.loads((directory.parent / 'run.json').read_text(), object_pairs_hook=unique_object)
    require(run['returncode'] == 0 and run['status'].endswith('_capture_completed'), 'capture did not complete')
    require(run['fixture'] == 'public_channel' and run['bundle']['sha256']['input.i'] == DECK_SHA256
            and run['bundle']['sha256']['channel.exo'] == MESH_SHA256, 'unsupported public mesh/deck profile')
    require(run[hash_key] == hashes, 'exports differ from completed capture hashes')
    return run['binary_sha256']


def read_frames(directory):
    frames = {}
    for path in sorted(directory.glob('*.jsonl')):
        rows = [json.loads(line, object_pairs_hook=unique_object) for line in path.read_text().splitlines()]
        h = rows[0]
        require(h['ranks'] == 1 and h['rank'] == 0, 'this integration gate requires one reference rank')
        key = h['stage'], h['call']
        frames[key] = {row['parent']: row for row in rows if row['kind'] == 'inputs'}
    require(set(frames) == {(s,c) for s in ('momentum.interior','pressure.interior') for c in (1,2)},
            'need both captured fields at both iterations')
    return frames


def build_data(interior, boundary):
    capture = load_dump(interior, require_inputs=True)
    require(capture['producer'] == 'openaccel' and capture['signature'][0] == 'public_channel', 'not a public reference')
    binary = verify_run(interior, capture['hashes'], 'export_sha256')
    boundary_records, boundary_hashes = load_boundary(boundary)
    # Boundary exports are under run/exports/boundary; run metadata is two levels above.
    boundary_run = boundary.parent.parent / 'run.json'
    record = json.loads(boundary_run.read_text(), object_pairs_hook=unique_object)
    require(record['returncode'] == 0 and record['status'] == 'boundary_capture_completed'
            and record['bundle']['sha256']['input.i'] == DECK_SHA256
            and record['bundle']['sha256']['channel.exo'] == MESH_SHA256
            and record['boundary_sha256'] == boundary_hashes, 'boundary capture provenance mismatch')
    frames = read_frames(interior)
    base = frames[('pressure.interior',1)]
    coordinates, faces = {}, defaultdict(list)
    for parent, row in sorted(base.items()):
        nodes, xyz = row['nodes'], row['fields']['coordinates']
        require(row['edges'] == [nodes[n] for pair in EDGES for n in pair], 'unsupported Tet4 ordering')
        for i,n in enumerate(nodes):
            point = xyz[3*i:3*i+3]
            require(n not in coordinates or coordinates[n] == point, 'inconsistent node coordinates')
            coordinates[n] = point
        for ordinal, face in enumerate(FACES):
            ids = tuple(nodes[n] for n in face)
            faces[tuple(sorted(ids))].append((parent, ordinal, ids))
    require(all(len(v) <= 2 for v in faces.values()), 'nonmanifold reference mesh')
    exterior = {key: value[0] for key,value in faces.items() if len(value) == 1}
    require((len(coordinates),len(base),len(exterior)) == (425,1536,576), 'incomplete public topology')
    boundaries = {}
    for (stage,call,face), row in boundary_records.items():
        if stage >= 3 or call != 1: continue
        require(row['rank'] == 0, 'distributed boundary capture is not supported here')
        x = row['inputs']; ids = [row['nodes'][int(n)] for n in x['face_nodes']]
        key = tuple(sorted(ids))
        require(key in exterior and key not in boundaries, 'missing/duplicate exterior boundary')
        parent,ordinal,ordered = exterior[key]
        nearest = [row['nodes'][int(n)] for n in x['nearest']]
        area, shape = [], []
        for n in ordered:
            s = nearest.index(n)
            area.extend(x['area'][3*s:3*s+3])
            shape.extend(x['shape'][3*s+ids.index(column)] for column in ordered)
        boundaries[key] = parent,ordinal,area,shape
    require(boundaries.keys() == exterior.keys(), 'boundary capture does not close the public mesh')
    nodal_frames = []
    for (stage,call), frame in sorted(frames.items()):
        require(frame.keys() == base.keys(), 'element coverage changed')
        components = 1 if stage == 'pressure.interior' else 3
        values, gradients = {}, {}
        field_name, gradient_name = ('pressure','pressure_gradient') if components == 1 else ('velocity','velocity_gradient')
        for parent,row in frame.items():
            f = row['fields']; original = base[parent]
            require(row['nodes'] == original['nodes'] and row['edges'] == original['edges']
                    and f['coordinates'] == original['fields']['coordinates'], 'mesh changed between calls')
            require(f['velocity_shifted'] == [0] and f['gradient_shifted'] == [1]
                    and f['compressible'] == [0], 'unsupported interpolation/physics')
            for name in ('velocity_shape','coordinate_shape','shape_gradient'):
                require(f[name] == original['fields'][name], 'geometry changed between captures')
            for i,n in enumerate(row['nodes']):
                for target,name,width in ((values,field_name,components),(gradients,gradient_name,3*components)):
                    item = f[name][i*width:(i+1)*width]
                    require(n not in target or target[n] == item, 'inconsistent shared nodal field')
                    target[n] = item
        nodal_frames.append((components,call,values,gradients))
    return dict(base=base, coordinates=coordinates, boundaries=boundaries, frames=nodal_frames,
                samples=capture['samples'], interior_sha256=capture['hashes'], boundary_sha256=boundary_hashes,
                interior_binary=binary, boundary_binary=record['binary_sha256'])


def prepare(interior, boundary, output):
    require(not output.exists() and not output.with_suffix('.json').exists(), 'output/provenance exists')
    data = build_data(interior,boundary)
    nodes = sorted(data['coordinates']); parents = sorted(data['base'])
    node_index = {n:i for i,n in enumerate(nodes)}; element_index = {n:i for i,n in enumerate(parents)}
    text = io.StringIO()
    def line(values): text.write(' '.join(format(v,'.17g') if isinstance(v,float) else str(v) for v in values)+'\n')
    text.write('MARS_PUBLIC_GEOMETRY_V1 425 1536 576 4\n')
    for n in nodes: line([n]+data['coordinates'][n])
    for p in parents:
        row = data['base'][p]; f = row['fields']; ids = row['nodes']
        line([p]+[node_index[n] for n in ids])
        areas = []
        for left,right in EDGES:
            a,b = ids[left],ids[right]
            sample = data['samples'][('pressure.interior',1,p,min(a,b),max(a,b))]
            areas.extend((1 if a < b else -1)*v for v in sample['area'])
        line(f['shape_gradient']+areas+f['velocity_shape']+f['coordinate_shape'])
    for key,(parent,ordinal,area,shape) in sorted(data['boundaries'].items()):
        line([element_index[parent],ordinal]); line(area+shape)
    for components,call,values,gradients in data['frames']:
        line([components,call])
        for n in nodes: line(values[n]+gradients[n])
    with output.open('x') as handle: handle.write(text.getvalue())
    metadata = {k:v for k,v in data.items() if k.endswith('sha256') or k.endswith('binary')}
    metadata.update(scope='public native Tet4 geometry and assembled nodal gradients; no flow solve',
                    packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
    output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Prepared public geometry: 425 nodes, 1536 tetrahedra, 576 boundary faces, 4 gradient states')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('interior',type=Path); parser.add_argument('--boundary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    try: prepare(args.interior,args.boundary,args.output)
    except (ValueError,KeyError,TypeError,OSError) as error: parser.exit(1,'ERROR: '+str(error)+'\n')
