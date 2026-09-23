"""Pack topology, boundary tags and independent two-iteration SIMPLE state oracles."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from openaccel_reference_check import require
from openaccel_boundary_check import load_boundary
from openaccel_update_check import load_updates
from prepare_openaccel_geometry import build_data, verify_run, FACES, EDGES


def prepare(interior, boundary, updates, output, mesh_only=False):
    require(not output.exists() and not output.with_suffix('.json').exists(), 'output/provenance exists')
    geometry = build_data(interior, boundary)
    b, _ = load_boundary(boundary)
    ids = sorted(geometry['coordinates']); parents = sorted(geometry['base'])
    ni = {n:i for i,n in enumerate(ids)}; ei = {n:i for i,n in enumerate(parents)}
    face_info = {}
    for (stage,call,identity), row in b.items():
        if stage >= 3 or call != 1: continue
        x = row['inputs']; local = [row['nodes'][int(j)] for j in x['face_nodes']]
        key = tuple(sorted(local)); parent, ordinal, _, _ = geometry['boundaries'][key]
        ordered = [geometry['base'][parent]['nodes'][j] for j in FACES[ordinal]]
        nearest = [row['nodes'][int(j)] for j in x['nearest']]
        face_info[key] = (parent,ordinal,stage,identity,[nearest.index(n) for n in ordered])
    faces = [face_info[k] for k in sorted(face_info)]
    require(len(faces)==576, 'incomplete boundary tags')
    out = io.StringIO()
    def line(values):
        out.write(' '.join(format(v,'.17g') if isinstance(v,float) else str(v) for v in values)+'\n')
    out.write('MARS_PUBLIC_SIMPLE_MESH_V1 425 1536 576\n' if mesh_only else 'MARS_PUBLIC_SIMPLE_V1 425 1536 576 2\n')
    for n in ids: line(geometry['coordinates'][n])
    for p in parents: line([ni[n] for n in geometry['base'][p]['nodes']])
    for p,f,kind,identity,mapping in faces: line([ei[p],f,kind])
    if mesh_only:
        metadata = {k:v for k,v in geometry.items() if k.endswith('sha256') or k.endswith('binary')}
        with output.open('x') as handle: handle.write(out.getvalue())
        metadata.update(packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                        scope='public channel topology and boundary tags only; no state oracles',
                        node_global_ids=ids)
        output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
        print('Prepared public mesh-only SIMPLE input: '+str(output))
        return
    require(updates is not None, '--updates is required unless --mesh-only is used')
    u, hashes = load_updates(updates)
    binary = verify_run(updates.parent, hashes, 'update_sha256')
    require({it for s,it,n,k in u} == {1,2}, 'requires exactly two captured iterations')
    for it in (1,2):
        for stage in ('momentum','increment','pressure','velocity','influence','trace','interior_flux','boundary_flux','mass_divergence'):
            values = []
            if stage in ('momentum','velocity','influence'):
                for n in ids:
                    r = u[2,it,n,0]
                    values.extend(r['outputs'] if stage=='velocity' else r['inputs'][3:6] if stage=='influence' else r['inputs'][:3])
            elif stage in ('increment','pressure'):
                values = [u[1 if stage=='increment' else 0,it,n,0]['outputs'][0] for n in ids]
            elif stage == 'trace':
                for p,f,kind,identity,mapping in faces:
                    values.extend([u[8,it,identity,s]['outputs'][0] for s in mapping] if kind==1 else [0.,0.,0.])
            elif stage == 'interior_flux':
                values = [u[3,it,p,s]['outputs'][0] for p in parents for s in range(6)]
            elif stage == 'boundary_flux':
                for p,f,kind,identity,mapping in faces:
                    if kind==0: values.extend(u[4,it,identity,s]['outputs'][0] for s in mapping)
                    elif kind==1:
                        reversal = u[6,it,identity,0]['outputs']
                        require(reversal[3:]==[0,0,0], 'native reversal selection is outside this gate')
                        values.extend(reversal[s] for s in mapping)
                    else: values.extend([0.,0.,0.])
            else:
                values = [0.]*len(ids)
                for p in parents:
                    nodes = geometry['base'][p]['nodes']
                    for s,(left,right) in enumerate(EDGES):
                        q = u[3,it,p,s]['outputs'][0]
                        values[ni[nodes[left]]] += q; values[ni[nodes[right]]] -= q
                for p,f,kind,identity,mapping in faces:
                    if kind==2: continue
                    nodes = geometry['base'][p]['nodes']
                    for j,s in enumerate(mapping):
                        q = u[4,it,identity,s]['outputs'][0] if kind==0 else u[6,it,identity,0]['outputs'][s]
                        values[ni[nodes[FACES[f][j]]]] += q
            line(values)
    with output.open('x') as handle: handle.write(out.getvalue())
    metadata = {k:v for k,v in geometry.items() if k.endswith('sha256') or k.endswith('binary')}
    metadata.update(update_sha256=hashes,update_binary=binary,
                    packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                    scope='two native single-rank public SIMPLE iterations; no reference state injected')
    output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Prepared public SIMPLE topology and two independent state oracles: '+str(output))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('interior',type=Path)
    parser.add_argument('--boundary',type=Path,required=True)
    parser.add_argument('--updates',type=Path)
    parser.add_argument('--mesh-only',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try: prepare(args.interior,args.boundary,args.updates,args.output,args.mesh_only)
    except (ValueError,KeyError,TypeError,OSError) as error: parser.exit(1,'ERROR: '+str(error)+'\n')
