"""Pack the pinned public channel's four complete frozen assembly states."""
import argparse
import hashlib
import io
import json
from pathlib import Path

from openaccel_reference_check import load_dump, require
from openaccel_boundary_check import load_boundary, ARRAYS, MAPS
from openaccel_node_check import load_nodes
from prepare_openaccel_geometry import build_data, read_frames, verify_run, FACES

INTERIOR_FIELDS = ('coordinates','velocity','density','viscosity','velocity_blend',
                   'stored_flux','pressure','density_blend','density_gradient')
BOUNDARY_FIELDS = ('velocity','boundary_velocity','viscosity','density','pressure',
                   'stored_flux','wall_coefficient')
NODE_FIELDS = ('density','volume','pseudo_dt','mass_divergence','velocity','pressure_gradient','force','source')


def near(actual, expected, message):
    require(len(actual) == len(expected), message+' size')
    require(all(abs(a-b) <= 1e-12*max(1,abs(b)) for a,b in zip(actual,expected)), message)


def prepare(interior, boundary, nodes, output):
    require(not output.exists() and not output.with_suffix('.json').exists(), 'output/provenance exists')
    geometry = build_data(interior,boundary)
    reference = load_dump(interior,require_inputs=True)
    boundaries, _ = load_boundary(boundary)
    nodal, node_hashes = load_nodes(nodes)
    node_binary = verify_run(nodes.parent,node_hashes,'node_sha256')
    require({c for s,c,n in nodal} == {1,2} and all(r['rank'] == 0 for r in nodal.values()),
            'requires two single-rank node states')
    frames = read_frames(interior)
    ids = sorted(geometry['coordinates']); parents = sorted(geometry['base'])
    ni = {n:i for i,n in enumerate(ids)}; ei = {n:i for i,n in enumerate(parents)}
    connectivity = [geometry['base'][p]['nodes'] for p in parents]
    neighbors = {n:{n} for n in ids}
    for cell in connectivity:
        for n in cell: neighbors[n].update(cell)
    positions = {}; offsets = [0]; columns = []
    for n in ids:
        for col in sorted(neighbors[n]):
            positions[(n,col)] = len(columns); columns.append(ni[col])
        offsets.append(len(columns))
    exterior = {}
    for parent,row in geometry['base'].items():
        for face,local in enumerate(FACES):
            key = tuple(sorted(row['nodes'][n] for n in local))
            if key in geometry['boundaries']: exterior[key] = (ei[parent],face)
    boundary_nodes = {n for key in exterior for n in key}
    out = io.StringIO()
    def line(values): out.write(' '.join(format(v,'.17g') if isinstance(v,float) else str(v) for v in values)+'\n')
    out.write('MARS_PUBLIC_ASSEMBLY_V1 425 1536 576 4 {}\n'.format(len(columns)))
    for n in ids: line(geometry['coordinates'][n])
    for cell in connectivity: line([ni[n] for n in cell])
    for key in sorted(exterior): line(exterior[key])
    line(offsets); line(columns)
    for call in (1,2):
        pressure_state = {}
        pressure_gradient = {}
        for row in frames[('pressure.interior',call)].values():
            for j,n in enumerate(row['nodes']):
                pressure_state[n] = row['fields']['pressure'][j]
                pressure_gradient[n] = row['fields']['pressure_gradient'][3*j:3*j+3]
        # Pressure is unchanged between this iteration's momentum and pressure assembly.
        for n in ids:
            near(nodal[(0,call,n)]['inputs']['pressure_gradient'],pressure_gradient[n],
                 'momentum/pressure reconstruction stage differs')
        for components,stage in ((3,'momentum.interior'),(1,'pressure.interior')):
            frame = frames[(stage,call)]
            field = {}
            for row in frame.values():
                f = row['fields']
                for j,n in enumerate(row['nodes']):
                    values = f['velocity'][3*j:3*j+3]+[pressure_state[n]]
                    require(n not in field or field[n] == values, 'inconsistent shared stage state')
                    field[n] = values
            line([components,call])
            for n in ids: line(field[n])
            lhs = [0.]*(len(columns)*components*components); rhs = [0.]*(len(ids)*components)
            def scatter(cell, block_lhs, block_rhs):
                width = len(cell)*components
                for r,n in enumerate(cell):
                    for i in range(components):
                        rhs[ni[n]*components+i] += block_rhs[r*components+i]
                        for c,col in enumerate(cell):
                            k = positions[(n,col)]*components*components+i*components
                            for j in range(components): lhs[k+j] += block_lhs[(r*components+i)*width+c*components+j]
            for parent in parents:
                row = frame[parent]; f = row['fields']; cell = row['nodes']
                for name in INTERIOR_FIELDS:
                    width = {'viscosity':4,'velocity_blend':12,'stored_flux':6,'pressure':4,
                             'density_blend':4,'density_gradient':12}.get(name)
                    values = f.get(name,[0.]*width) if width is not None else f[name]
                    line(values)
                block = reference['blocks'][(stage,call,parent)]
                dofs = [(n,i) for n in cell for i in range(components)]
                scatter(cell,[block['lhs'][(a,b)] for a in dofs for b in dofs],[block['rhs'][a] for a in dofs])
            selected = [(s,r) for (s,c,f),r in sorted(boundaries.items()) if c == call and (s >= 3) == (components == 3)]
            require(len(selected) == 576, 'incomplete boundary stage')
            for s,row in selected:
                x = row['inputs']; cell = row['nodes']; local = [cell[int(i)] for i in x['face_nodes']]
                element,face = exterior[tuple(sorted(local))]
                line([s,element,face]+[ni[n] for n in cell]+[-1]*(4-len(cell)))
                for name in MAPS: line(x[name])
                for name in BOUNDARY_FIELDS:
                    values = x.get(name,[])
                    line(values+[0.]*(ARRAYS[name]-len(values)))
                scatter(cell,row['outputs']['lhs'],row['outputs']['rhs'])
            if components == 3:
                require({n for s,c,n in nodal if s == 3 and c == call} == boundary_nodes,
                        'captured boundary relaxation differs from exterior-node union')
                for n in ids:
                    row = nodal[(0,call,n)]; x = row['inputs']
                    near(x['velocity'],field[n][:3],'node/interior velocity mismatch')
                    for name in NODE_FIELDS:
                        line([0.] if name == 'volume' else [0.]*3 if name == 'pressure_gradient' else x[name])
                    scatter([n],row['outputs']['lhs'],row['outputs']['rhs'])
                for n in ids:
                    alpha = nodal[(1,call,n)]['inputs']['alpha'][0]
                    factor = .75 if n in boundary_nodes else 1.
                    require(nodal[(2,call,n)]['inputs']['consistent'] == [0], 'this profile requires SIMPLE')
                    line([alpha,factor])
                    pos = 9*positions[(n,n)]
                    near(lhs[pos:pos+9],nodal[(1,call,n)]['inputs']['lhs'],'assembled diagonal differs from actual reference')
                    for j in range(3): lhs[pos+3*j+j] *= 1/alpha
                    if n in boundary_nodes:
                        near(rhs[3*ni[n]:3*ni[n]+3],nodal[(3,call,n)]['inputs']['rhs'],
                             'assembled boundary RHS differs from actual reference')
                    for j in range(3): rhs[3*ni[n]+j] *= factor
            line(lhs); line(rhs)
            if components == 3:
                for n in ids: line(nodal[(2,call,n)]['outputs']['d'])
    with output.open('x') as handle: handle.write(out.getvalue())
    metadata = {k:v for k,v in geometry.items() if k.endswith('sha256') or k.endswith('binary')}
    metadata.update(node_sha256=node_hashes,node_binary=node_binary,
                    packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                    blocks=len(columns),scope='four public frozen assembled systems; no linear solves or SIMPLE iteration')
    output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Prepared four public assembled systems: {} node blocks; reference diagonal/RHS joins passed'.format(len(columns)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('interior',type=Path)
    parser.add_argument('--boundary',type=Path,required=True)
    parser.add_argument('--nodes',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    try: prepare(args.interior,args.boundary,args.nodes,args.output)
    except (ValueError,KeyError,TypeError,OSError) as error: parser.exit(1,'ERROR: '+str(error)+'\n')
