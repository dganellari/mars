#!/usr/bin/env python3
"""GPT/Codex, 2026-09-10: public channel fixture, independent of any case mesh.

Regeneration needs numpy and netCDF4; running the committed fixture needs neither.
"""
import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from generate_tet_cube import KUHN_TETS


def geometry(opening_width=1.):
    nx, ny, nz = 16, 4, 4
    points = np.array([(4*i/nx, j/ny, k/nz)
                       for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)])
    plane = (nx+1)*(ny+1)
    offsets = np.array([0, 1, nx+2, nx+1, plane, plane+1, plane+nx+2, plane+nx+1])
    bases = np.array([k*plane+j*(nx+1)+i
                      for k in range(nz) for j in range(ny) for i in range(nx)])
    tets = (bases[:, None, None]+offsets[KUHN_TETS][None, :, :]).reshape(-1, 4)
    # Exodus TET4 side numbers; checked against geometry, not just the reader's table.
    faces = ((0, 1, 3), (1, 2, 3), (0, 3, 2), (0, 2, 1))
    incidence = {}
    patches = {name: [] for name in ('inlet', 'outlet', 'walls')}
    for element, tet in enumerate(tets):
        for side, face in enumerate(faces, 1):
            nodes = tet[list(face)]
            incidence.setdefault(tuple(sorted(nodes)), []).append((element, side, nodes))
    for entries in incidence.values():
        if len(entries) == 2:
            continue
        assert len(entries) == 1
        element, side, nodes = entries[0]
        x = points[nodes, 0]
        opening = np.all(points[nodes, 1:] <= opening_width)
        name = 'inlet' if opening and np.all(x == 0) else 'outlet' if opening and np.all(x == 4) else 'walls'
        patches[name].append((element+1, side))
    return points, tets, faces, patches


def verify(path, opening_width=1.):
    # Read back the serialized mesh, not the arrays used to write it.
    with Dataset(path) as mesh:
        xyz = np.column_stack([mesh.variables['coord'+d][:] for d in 'xyz'])
        tets = np.asarray(mesh.variables['connect1'][:], dtype=int)-1
        p = xyz[tets]
        det = np.linalg.det(np.stack([p[:, k]-p[:, 0] for k in (1, 2, 3)], axis=2))
        assert np.all(det > 0) and abs(det.sum()/6-4) < 1e-12
        faces = ((0, 1, 3), (1, 2, 3), (0, 3, 2), (0, 2, 1))
        occurrences = {}
        for tet in tets:
            for face in faces:
                key = tuple(sorted(tet[list(face)]))
                occurrences[key] = occurrences.get(key, 0)+1
        assert all(count in (1, 2) for count in occurrences.values())
        exterior = {key for key, count in occurrences.items() if count == 1}
        seen, total_area = set(), np.zeros(3)
        opening_area = opening_width**2
        for patch, expected_area in enumerate((opening_area, opening_area, 18.-2*opening_area), 1):
            area = 0.
            for element, side in zip(mesh.variables[f'elem_ss{patch}'][:],
                                     mesh.variables[f'side_ss{patch}'][:]):
                tet = tets[int(element)-1]
                local = faces[int(side)-1]
                nodes = tet[list(local)]
                key = tuple(sorted(nodes))
                assert key in exterior and key not in seen
                seen.add(key)
                a, b, c = xyz[nodes]
                normal = .5*np.cross(b-a, c-a)
                opposite = next(k for k in range(4) if k not in local)
                assert np.dot(normal, xyz[tet[opposite]]-a) < 0
                if patch in (1, 2):
                    assert np.all(xyz[nodes, 1:] <= opening_width)
                if patch == 1:
                    assert np.all(xyz[nodes, 0] == 0) and normal[0] < 0
                if patch == 2:
                    assert np.all(xyz[nodes, 0] == 4) and normal[0] > 0
                area += np.linalg.norm(normal)
                total_area += normal
            assert abs(area-expected_area) < 1e-12
        assert seen == exterior and np.linalg.norm(total_area) < 1e-12
    print('PASS: public channel read-back, positive tets, volume, boundary closure, areas and winding')


def write(path, opening_width=1.):
    if opening_width not in (1., .25):
        raise ValueError("supported public opening widths are 1 and 0.25")
    xyz, tets, _, patches = geometry(opening_width)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, 'w', format='NETCDF3_CLASSIC') as mesh:
        mesh.title = ('PUBLIC_SYNTHETIC_OUTLET_CHANNEL_V1' if opening_width == 1.
                      else 'PUBLIC_SYNTHETIC_CORNER_OPENINGS_V1')
        mesh.api_version = np.float32(7.22)
        mesh.version = np.float32(7.22)
        mesh.floating_point_word_size = np.int32(8)
        for name, size in dict(num_dim=3, num_nodes=len(xyz), num_elem=len(tets),
                               num_el_blk=1, num_el_in_blk1=len(tets), num_nod_per_el1=4,
                               num_side_sets=3, len_name=33).items():
            mesh.createDimension(name, size)
        for d, axis in enumerate('xyz'):
            mesh.createVariable('coord'+axis, 'f8', ('num_nodes',))[:] = xyz[:, d]
        conn = mesh.createVariable('connect1', 'i4', ('num_el_in_blk1', 'num_nod_per_el1'))
        conn.elem_type = 'TETRA4'
        conn[:] = tets+1
        for prefix, dimension, ids in (('eb', 'num_el_blk', [1]), ('ss', 'num_side_sets', [1,2,3])):
            prop = mesh.createVariable(prefix+'_prop1', 'i4', (dimension,))
            prop.setncattr('name', 'ID')
            prop[:] = ids
            mesh.createVariable(prefix+'_status', 'i4', (dimension,))[:] = 1
        names = np.zeros((3, 33), dtype='S1')
        for i, (name, entries) in enumerate(patches.items(), 1):
            names[i-1, :len(name)] = np.frombuffer(name.encode('ascii'), dtype='S1')
            dim = f'num_side_ss{i}'
            mesh.createDimension(dim, len(entries))
            for j, kind in enumerate(('elem', 'side')):
                mesh.createVariable(f'{kind}_ss{i}', 'i4', (dim,))[:] = np.array(entries)[:, j]
        mesh.createVariable('ss_names', 'S1', ('num_side_sets', 'len_name'))[:] = names
    verify(path, opening_width)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--opening-width', type=float, choices=(1., .25), default=1.)
    args = parser.parse_args()
    write(args.output, args.opening_width)
