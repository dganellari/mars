#!/usr/bin/env python3
"""Prepare native Hex8/Tet4 partitions for the stationary PI channel experiment."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import struct
import numpy as np
from prepare_nonconformal_mesh import read_exodus, names, require

FACES = {4: ((0,1,3),(1,2,3),(0,3,2),(0,2,1)),
         8: ((0,1,5,4),(1,2,6,5),(2,3,7,6),(0,4,7,3),(0,3,2,1),(4,5,6,7))}
LABELS = ('inlet','outlet','wall','interface_hex','interface_tet')

def read(path):
    with read_exodus(path) as f:
        xyz=np.column_stack([f.variables['coord'+c][:] for c in 'xyz']).astype('<f8')
        require(xyz.ndim==2 and xyz.shape[1]==3 and np.isfinite(xyz).all(), 'invalid coordinates')
        cells=[]; counts=[]; blocks=[]
        for b in (1,2):
            v=f.variables['connect'+str(b)]; kind=v.elem_type
            if isinstance(kind,bytes):kind=kind.decode()
            count=8 if kind.strip().upper() in ('HEX','HEX8','HEXAHEDRON') else 4 if kind.strip().upper() in ('TETRA','TETRA4','TET4') else 0
            a=np.asarray(v[:],dtype='<i4')-1
            require(count>0 and a.shape[1]==count and len(a)>0,'only native Hex8 and Tet4 blocks are supported')
            require(a.min()>=0 and a.max()<len(xyz),'invalid node index')
            require(np.all(np.diff(np.sort(a,axis=1),axis=1)>0),'duplicate cell node')
            blocks.append(set(a.ravel()))
            for row in a: cells.append(list(row)+[-1]*(8-count));counts.append(count)
        require(set(counts)=={4,8} and not blocks[0].intersection(blocks[1]),'need disjoint Hex8/Tet4 blocks')
        require(blocks[0].union(blocks[1])==set(range(len(xyz))),'unused node')
        for key,n in [('node_num_map',len(xyz)),('elem_num_map',len(cells))]:
            if key in f.variables:require(np.array_equal(f.variables[key][:],np.arange(1,n+1)),'this driver requires consecutive original Exodus IDs')
        labels=names(f.variables['ss_names']);require(set(labels)==set(LABELS) and len(labels)==5,'expected five PI side sets')
        faces=[]
        for i,label in enumerate(labels,1):
            for e,s in zip(f.variables['elem_ss'+str(i)][:],f.variables['side_ss'+str(i)][:]):
                e,s=int(e)-1,int(s)-1
                require(0<=e<len(cells) and 0<=s<len(FACES[counts[e]]),'invalid side ordinal')
                require(label!='interface_hex' or counts[e]==8,'hex interface is not on a hex')
                require(label!='interface_tet' or counts[e]==4,'tet interface is not on a tet')
                faces.append((e,s,LABELS.index(label)))
    cells=np.asarray(cells,dtype='<i4');counts=np.asarray(counts,dtype='<i4')
    incidence=defaultdict(list)
    for e,c in enumerate(cells):
        for f,corners in enumerate(FACES[counts[e]]):incidence[tuple(sorted(c[list(corners)]))].append((e,f))
    require(all(len(v)<=2 for v in incidence.values()),'nonmanifold face')
    exterior={v[0] for v in incidence.values() if len(v)==1}
    require(len({(e,f) for e,f,k in faces})==len(faces),'duplicate boundary face')
    require({(e,f) for e,f,k in faces}==exterior,'side sets must cover exactly the exterior faces')
    # This is an offline topology preflight. Runtime geometry/search use the GPU.
    parent=list(range(len(cells)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for v in incidence.values():
        if len(v)==2:parent[root(v[0][0])]=root(v[1][0])
    require(len({root(i) for i in range(len(cells))})==2,'each block must be connected')
    return xyz,cells,counts,np.asarray(faces,dtype='<i4')

def prepare(source,output,ranks):
    xyz,cells,counts,faces=read(source); n=len(xyz)
    require(1<=ranks<=n,'invalid MPI size')
    node_owner=np.searchsorted([n*(r+1)//ranks for r in range(ranks)],np.arange(n),side='right')
    corners=[cells[e,list(FACES[counts[e]][f])] for e,f,k in faces]
    boxes=[(xyz[c].min(axis=0),xyz[c].max(axis=0)) for c in corners]
    candidates={}
    for i,(e,f,kind) in enumerate(faces):
        if kind<3:continue
        candidates[i]=[j for j,(_,_,other) in enumerate(faces) if other==7-kind and
                       np.all(boxes[i][0]<=boxes[j][1]+1e-12) and np.all(boxes[j][0]<=boxes[i][1]+1e-12)]
        require(bool(candidates[i]),'interface face has no bounding-box donor')
    partitions=[];peers=[set() for _ in range(ranks)]
    for rank in range(ranks):
        incident={e for e,c in enumerate(cells) if np.any(node_owner[c[:counts[e]]]==rank)}
        needed=set(incident); receiving=[]
        for f,donors in candidates.items():
            if any(node_owner[node]==rank for node in corners[f]):
                receiving.append(f);needed.update(int(faces[d,0]) for d in donors)
        es=sorted(needed); emap={e:i for i,e in enumerate(es)}
        fs=[f for f,(e,_,_) in enumerate(faces) if e in needed];fmap={f:i for i,f in enumerate(fs)}
        ids=sorted({int(k) for e in es for k in cells[e,:counts[e]]});nmap={v:i for i,v in enumerate(ids)}
        owned=set(range(n*rank//ranks,n*(rank+1)//ranks));require(owned.issubset(ids),'missing owned nodes')
        peers[rank].update(int(node_owner[k]) for k in ids if node_owner[k]!=rank)
        cs=np.array([[counts[e],e]+[nmap[int(k)] if k>=0 else -1 for k in cells[e]] for e in es],dtype='<i4')
        # Identity follows original topology, not the order of side-set records.
        ff=np.array([[emap[int(faces[f,0])],int(faces[f,1]),int(faces[f,2]),6*int(faces[f,0])+int(faces[f,1])] for f in fs],dtype='<i4')
        qs=[];ds=[]
        for f in receiving:
            for q,node in enumerate(corners[f]):
                if node_owner[node]!=rank:continue
                qs.append([fmap[f],q,len(ds),len(candidates[f]),nmap[int(node)]])
                ds.extend(fmap[d] for d in candidates[f])
        partitions.append((np.array(ids,dtype='<i8'),xyz[ids],cs,ff,np.asarray(qs,dtype='<i4').reshape(-1,5),np.array(ds,dtype='<i4')))
    for rank in range(ranks):
        for peer in list(peers[rank]):peers[peer].add(rank)
    output.mkdir(parents=True,exist_ok=False)
    for rank,arrays in enumerate(partitions):
        ids,x,c,f,q,d=arrays;ps=np.asarray(sorted(peers[rank]),dtype='<i4')
        with (output/('rank%d.mmx'%rank)).open('xb') as out:
            out.write(b'MARSMX01'+struct.pack('<9Q',ranks,rank,n,len(ids),len(c),len(f),len(q),len(d),len(ps)))
            for a in (ids,x,c,f,q,d,ps):out.write(a.tobytes(order='C'))
    manifest={'source':str(source.resolve()),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'ranks':ranks,'nodes':n,'hexes':int((counts==8).sum()),'tets':int((counts==4).sum()),
              'format':'MARSMX01','geometry':'native elements; no conversion or welding',
              'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.glob('rank*.mmx'))}}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mesh',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--ranks',type=int,required=True)
    a=p.parse_args()
    try:print(json.dumps(prepare(a.mesh,a.output,a.ranks),indent=2))
    except (ValueError,OSError,KeyError) as e:p.exit(1,str(e)+'\n')
