#!/usr/bin/env python3
"""Partition a prepared stationary Tet4 channel without merging interface nodes."""
import argparse
import json
from pathlib import Path
import struct
import numpy as np


def read_prepared(path):
    data = Path(path).read_bytes()
    if len(data) < 40 or data[:8] != b'MARSNC02':
        raise ValueError('expected input prepared with --flow-boundaries')
    n, e, t, f = struct.unpack_from('<4Q', data, 8)
    if not 0 < e <= 1000000 or not 0 < n <= 4*e or not 0 < t <= 4*e or not 0 < f <= 4*e:
        raise ValueError('invalid prepared dimensions')
    if len(data) != 40 + 32*n + 24*e + 12*(t+f):
        raise ValueError('prepared length mismatch')
    offset = 40
    arrays = []
    for dtype, shape in [('<f8',(n,3)), ('<i4',(e,4)), ('<f8',(e,)),
                         ('<i4',(n,)), ('<i4',(n,)), ('<i4',(t,3)), ('<i4',(f,3))]:
        count = int(np.prod(shape))
        a = np.frombuffer(data, dtype=dtype, count=count, offset=offset).reshape(shape)
        offset += a.nbytes
        arrays.append(a)
    x, cells, k, boundary, side, tags, outer = arrays
    if not np.isfinite(x).all() or not np.isfinite(k).all() or (k <= 0).any():
        raise ValueError('invalid prepared geometry or diffusivity')
    if cells.min() < 0 or cells.max() >= n or (np.diff(np.sort(cells,axis=1),axis=1) == 0).any():
        raise ValueError('invalid prepared connectivity')
    if (boundary < 0).any() or (boundary > 2).any() or (side < 0).any() or (side > 1).any():
        raise ValueError('invalid node labels')
    for a, limit in ((tags, 1), (outer, 2)):
        if (a < 0).any() or (a[:,0] >= e).any() or (a[:,1] >= 4).any() or (a[:,2] > limit).any():
            raise ValueError('invalid face tags')
    return arrays


def partitions(arrays, ranks):
    x, cells, k, boundary, side, tags, outer = arrays
    n = len(x)
    if not 1 <= ranks <= min(n, 65536):
        raise ValueError('rank count outside supported range')
    # Contiguous global row ranges preserve the existing example's ownership.
    ends = np.arange(1, ranks+1, dtype=np.int64)*n//ranks
    owners = np.searchsorted(ends, np.arange(n), side='right').astype('<i4')
    cell_owners = owners[cells]
    element_owners = cell_owners.min(axis=1)
    incident = [set() for _ in range(ranks)]
    peer = [set() for _ in range(ranks)]
    for e, row in enumerate(cell_owners):
        participants = set(map(int,row))
        for r in participants:
            incident[r].add(e)
            peer[r].update(participants-{r})
    # All corner owners are included: an off-face node can be a gradient donor.
    bounds = [[np.full(3,np.inf),np.full(3,-np.inf)] for _ in range(ranks)]
    for e, opposite, _ in tags:
        vertices = x[np.delete(cells[e], opposite)]
        for r in set(map(int,cell_owners[e])):
            bounds[r][0] = np.minimum(bounds[r][0], vertices.min(axis=0)-1e-11)
            bounds[r][1] = np.maximum(bounds[r][1], vertices.max(axis=0)+1e-11)
    active=[]
    for r in sorted(range(ranks),key=lambda r:bounds[r][0][0]):
        lo,hi=bounds[r]
        if not np.isfinite(lo).all(): continue
        active=[s for s in active if bounds[s][1][0] >= lo[0]]
        for s in active:
            if np.all(bounds[s][1] >= lo) and np.all(hi >= bounds[s][0]):
                peer[r].add(s)
                peer[s].add(r)
        active.append(r)
    for rank in range(ranks):
        elements=np.array(sorted(incident[rank]),dtype='<u8')
        nodes=np.unique(cells[elements].ravel()).astype('<i8')
        local_cells=np.searchsorted(nodes,cells[elements]).astype('<i4')
        local_tags=[]
        for source in (tags,outer):
            mask=np.isin(source[:,0],elements)
            a=source[mask].copy()
            a[:,0]=np.searchsorted(elements,a[:,0])
            local_tags.append(a.astype('<i4'))
        yield (nodes, elements, x[nodes], local_cells, k[elements], boundary[nodes], side[nodes],
               element_owners[elements], *local_tags, np.array(sorted(peer[rank]),dtype='<i4'))


def write_partitions(source, output, ranks):
    arrays=read_prepared(source)
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    counts=[]
    for rank, part in enumerate(partitions(arrays,ranks)):
        nodes,elements,x,cells,k,boundary,side,owners,tags,outer,peer=part
        with (output/f'rank{rank}.mncp').open('xb') as f:
            f.write(b'MARSNP01'+struct.pack('<9Q',ranks,rank,len(arrays[0]),len(nodes),len(elements),
                                            len(tags),len(outer),len(peer),len(arrays[1])))
            for a in part: f.write(a.tobytes(order='C'))
        counts.append({'rank':rank,'local_nodes':len(nodes),'local_cells':len(elements),'peers':len(peer)})
    report={'schema':'mars-nonconformal-partitions-v1','ranks':ranks,'global_nodes':len(arrays[0]),
            'global_cells':len(arrays[1]),'partitions':counts}
    (output/'partition.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ranks',type=int,required=True)
    a=p.parse_args()
    try: report=write_partitions(a.input,a.output,a.ranks)
    except (ValueError,OSError) as e: p.exit(1,f'Partitioning rejected: {e}\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
