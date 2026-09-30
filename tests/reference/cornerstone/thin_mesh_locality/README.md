# Decomposition locality on one-cell-thick meshes

`thin_mesh_locality.cpp` reproduces, with public cornerstone API only, a locality problem that MARS hit on a
planar channel mesh (one hex cell through z), and shows that the mixed-dimension (MixD) SFC keys on master
(#598, #622) fix it.

## What it measures

Each hex element is one particle, placed at its corner with the smallest SFC key (as MARS does) or at its
centroid. After `Domain::sync`, every mesh node is given to the rank whose SFC range holds the node position.
A finite-element code needs that rank to hold the elements around the node. Printed per run:

- **max owner distance**: largest distance, in cells, from a node to the nearest element assigned to its owner
- **far nodes**: nodes whose owner holds no assigned element within 2 cells
- **missing star**: (node, element) incidences where the owner holds the element neither as assigned element nor
  as halo; each one must be sent to the owner before it can assemble its row
- **halos**: halo particles over all ranks after `sync`

## Build and run

```bash
mpicxx -std=c++20 -O2 -DNDEBUG -I <cornerstone>/include thin_mesh_locality.cpp -o thin_mesh_locality
mpirun -np 16 ./thin_mesh_locality 320 32 1 10 1 0.06 corner   # nx ny nz Lx Ly Lz corner|centroid
```

## Results

Before MixD: master at 9d7e8e30 (parent of #598). After: master at d7fddfb4. The same source builds against
both. CPU domain, bucket size 64, theta 0.5, box padded by 5%.

Thin channel, 10 x 1 x 0.06, one cell through z, corner placement:

| Cells | Ranks | Version | Max owner distance | Far nodes | Missing star | Halos |
|-------|-------|---------|--------------------|-----------|--------------|-------|
| 320 x 32 x 1 | 4 | before | 7 | 395 of 21186 | 1285 | 2373 |
| 320 x 32 x 1 | 4 | after | 0 | 0 | 0 | 802 |
| 320 x 32 x 1 | 16 | before | 60 | 6316 of 21186 | 23410 | 9264 |
| 320 x 32 x 1 | 16 | after | 0 | 0 | 0 | 4305 |
| 1600 x 160 x 1 | 4 | before | 39 | 17437 of 515522 | 75255 | 14677 |
| 1600 x 160 x 1 | 4 | after | 0 | 0 | 0 | 9929 |
| 1600 x 160 x 1 | 16 | before | >= 64 | 183715 of 515522 | 710576 | 81587 |
| 1600 x 160 x 1 | 16 | after | 2 | 0 | 0 | 50150 |

Control, cube 40 x 40 x 40 on 1 x 1 x 1: both versions give identical results (4 ranks: 0 / 0 / 0 / 16898
halos; 16 ranks: 0 / 0 / 2 / 47029 halos).

Why it happens before MixD: the key box is normalized per axis, so the two node layers of the thin mesh sit at
opposite ends of normalized z, and the curve orders them separately. A rank's range then covers the bottom
layer in one region and the top layer in another, and nodes are owned far from their elements. MixD gives short
axes fewer key bits, so the cells stay close to isotropic in physical space.

In MARS on 64 GPUs (a 25298 x 2530 x 1 channel, 1M cells per GPU), the old keys made element-star completion
request 58M of the 64M elements, node halos covered about 70% of each rank, and the run aborted in setup. With
MixD keys the same ownership rule needs no extra elements on these meshes.

## A second observation: planar particle sets

With centroid placement all particles have the same z, the open box fitted in `sync` has zero extent there, and
the halo search finds nothing, without an error, on both versions:

| Cells | Ranks | Version | Halos | Missing star | Box z after sync |
|-------|-------|---------|-------|--------------|------------------|
| 320 x 32 x 1 | 4 | before | 0 | 1406 | [0.03, 0.03] |
| 320 x 32 x 1 | 4 | after | 0 | 384 | [0.03, 0.03] |
| 320 x 32 x 1 | 16 | before | 0 | 4206 | [0.03, 0.03] |
| 320 x 32 x 1 | 16 | after | 0 | 1920 | [0.03, 0.03] |

A floor on the fitted extent, or an error for a degenerate box, would make this visible.
