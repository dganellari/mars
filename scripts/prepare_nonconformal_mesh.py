#!/usr/bin/env python3
"""Read two disjoint Tet4 Exodus blocks for the stationary interface examples."""
import argparse
from collections import defaultdict
from pathlib import Path
import struct
import numpy as np

# Exodus Tet4 side number minus one -> opposite local corner.
OPPOSITE = (2, 0, 1, 3)
MAGIC = b'MARSNC01'
MAX_ELEMENTS = 1000000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_exodus(path):
    try:
        from netCDF4 import Dataset
    except ImportError:
        from scipy.io import netcdf_file
        return netcdf_file(path, 'r', mmap=False)
    result = Dataset(path, 'r')
    result.set_auto_mask(False)
    return result


def names(variable):
    return [b''.join(row).decode('utf-8').rstrip('\x00 ').strip() for row in np.asarray(variable[:], dtype='S1')]


def prepare(path, interface_left='interface_left', interface_right='interface_right',
            inlet='inlet', outlet='outlet', wall='wall', k_left=1., k_right=1., flow_boundaries=False):
    require(all(np.isfinite(k) and k > 0 for k in (k_left, k_right)), 'diffusivities must be positive and finite')
    selected = [inlet, outlet, interface_left, interface_right, wall]
    require(len(set(selected)) == 5, 'the five side-set names must be distinct')
    with read_exodus(path) as f:
        dim = lambda key: int(f.dimensions[key]) if isinstance(f.dimensions[key], int) else len(f.dimensions[key])
        require(dim('num_dim') == 3 and dim('num_el_blk') == 2, 'requires exactly two 3D Tet4 blocks')
        n, e = dim('num_nodes'), dim('num_elem')
        require(0 < e <= MAX_ELEMENTS and 0 < n <= 4 * e, 'mesh exceeds this small example input limit')
        x = (np.column_stack([f.variables['coord' + c][:] for c in 'xyz'])
             if 'coordx' in f.variables else np.asarray(f.variables['coord'][:]).T).astype('<f8')
        require(x.shape == (n, 3) and np.isfinite(x).all(), 'invalid coordinates')
        blocks = []
        for i in (1, 2):
            v = f.variables['connect' + str(i)]
            kind = v.elem_type
            if isinstance(kind, bytes):
                kind = kind.decode()
            require(kind.strip().upper() in ('TETRA', 'TETRA4', 'TET4'), 'only Tet4 is supported; mesh is not converted')
            data = np.asarray(v[:])
            require(data.dtype.kind in 'iu' and data.ndim == 2 and data.shape[1] == 4 and len(data) > 0,
                    'invalid Tet4 connectivity')
            require(data.min() >= 1 and data.max() <= n, 'connectivity outside Exodus node range')
            blocks.append(data.astype('<i4') - 1)
        cells = np.concatenate(blocks)
        require(len(cells) == e, 'element count mismatch')
        labels = names(f.variables['ss_names'])
        require(len(set(labels)) == len(labels) and set(labels) == set(selected),
                'every side set must be mapped exactly once to inlet, outlet, two interfaces or wall')
        sets = {}
        for i, name in enumerate(labels, 1):
            elem = np.asarray(f.variables['elem_ss' + str(i)][:])
            side = np.asarray(f.variables['side_ss' + str(i)][:])
            require(elem.dtype.kind in 'iu' and side.dtype.kind in 'iu' and elem.ndim == side.ndim == 1
                    and elem.shape == side.shape and len(elem) > 0, 'invalid or empty side set')
            require(elem.min() >= 1 and elem.max() <= e and side.min() >= 1 and side.max() <= 4,
                    'side-set element or face outside valid range')
            sets[name] = [(int(a) - 1, OPPOSITE[int(b) - 1]) for a, b in zip(elem, side)]
    node_side = np.full(n, -1, dtype='<i4')
    for side, block in enumerate(blocks):
        ids = np.unique(block)
        require(np.all(node_side[ids] == -1), 'blocks must have separate node IDs, including coincident interface nodes')
        node_side[ids] = side
    require(np.all(node_side >= 0), 'unused nodes are not supported')
    require(np.all(np.diff(np.sort(cells, axis=1), axis=1) > 0), 'repeated corner in a tetrahedron')
    xyz = x[cells]
    edge = xyz[:, 1:] - xyz[:, :1]
    scale = np.max(np.linalg.norm(xyz[:, :, None] - xyz[:, None, :], axis=-1), axis=(1, 2))
    require(np.all(scale > 0), 'collapsed tetrahedron')
    determinant = np.linalg.det(edge / scale[:, None, None])
    require(np.all(np.isfinite(determinant) & (np.abs(determinant) > 64 * np.finfo(float).eps)),
            'degenerate tetrahedron')
    # This is file preflight. Runtime graph construction and integration run on device.
    faces = defaultdict(list)
    for elem, cell in enumerate(cells):
        for opposite in range(4):
            faces[tuple(sorted(int(cell[a]) for a in range(4) if a != opposite))].append((elem, opposite))
    require(all(len(v) <= 2 for v in faces.values()), 'nonmanifold face')
    exterior = {v[0] for v in faces.values() if len(v) == 1}
    declared = [entry for entries in sets.values() for entry in entries]
    require(len(set(declared)) == len(declared), 'duplicate boundary face')
    require(set(declared) == exterior, 'side sets must cover every exterior face exactly once and no internal faces')
    boundary = np.zeros(n, dtype='<i4')
    for name, code in ((inlet, 1), (outlet, 2)):
        for elem, opposite in sets[name]:
            ids = [cells[elem, a] for a in range(4) if a != opposite]
            require(np.all((boundary[ids] == 0) | (boundary[ids] == code)), 'conflicting inlet and outlet values')
            boundary[ids] = code
    tags = []
    for side, name in enumerate((interface_left, interface_right)):
        for elem, opposite in sets[name]:
            require(np.all(node_side[cells[elem]] == side), 'interface side does not belong to its declared block')
            tags.append((elem, opposite, side))
    # Reject disconnected islands rather than leave an unanchored scalar nullspace.
    parent = list(range(e))
    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for entries in faces.values():
        if len(entries) == 2:
            parent[root(entries[0][0])] = root(entries[1][0])
    require(len({root(i) for i in range(e)}) == 2, 'each block must be face-connected')
    for block in blocks:
        require(np.any(boundary[np.unique(block)] != 0), 'each block needs an inlet or outlet anchor')
    k = np.concatenate((np.full(len(blocks[0]), k_left), np.full(len(blocks[1]), k_right))).astype('<f8')
    arrays = x, cells, k, boundary, node_side, np.asarray(tags, dtype='<i4')
    if flow_boundaries:
        outer = [(elem, opposite, kind) for kind, name in enumerate((inlet, outlet, wall))
                 for elem, opposite in sets[name]]
        arrays += (np.asarray(outer, dtype='<i4'),)
    return arrays


def write_mesh(output, arrays):
    x, cells, k, boundary, side, tags = arrays[:6]
    outer = arrays[6] if len(arrays) == 7 else None
    # Exclusive creation protects the PI's source and previous prepared inputs.
    with Path(output).open('xb') as f:
        f.write((MAGIC if outer is None else b'MARSNC02') + struct.pack('<QQQ', len(x), len(cells), len(tags)))
        if outer is not None:
            f.write(struct.pack('<Q', len(outer)))
        for data in (x.astype('<f8'), cells.astype('<i4'), k.astype('<f8'),
                     boundary.astype('<i4'), side.astype('<i4'), tags.astype('<i4')):
            f.write(data.tobytes(order='C'))
        if outer is not None:
            f.write(outer.astype('<i4').tobytes(order='C'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mesh', type=Path)
    p.add_argument('--output', required=True, type=Path)
    for key, value in [('interface-left', 'interface_left'), ('interface-right', 'interface_right'),
                       ('inlet', 'inlet'), ('outlet', 'outlet'), ('wall', 'wall')]:
        p.add_argument('--' + key, default=value)
    p.add_argument('--flow-boundaries', action='store_true', help='retain all face labels for the channel-flow example')
    p.add_argument('--k-left', type=float, default=1.)
    p.add_argument('--k-right', type=float, default=1.)
    a = p.parse_args()
    try:
        arrays = prepare(a.mesh, a.interface_left, a.interface_right, a.inlet, a.outlet, a.wall, a.k_left, a.k_right, a.flow_boundaries)
        write_mesh(a.output, arrays)
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, 'Mesh preparation rejected: ' + str(error) + '\n')
    print('Prepared stationary Tet4 input; original nodes and elements preserved. No solver launched.')


if __name__ == '__main__':
    main()
