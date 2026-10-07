# Poiseuille channel mesh

`poiseuille_hex_14k_elem.e` is the mesh used by the
[Poiseuille channel tutorial](../../../docs/poiseuille_tutorial.md) and by the
`marsReleasePoiseuille_np*` release checks. Distributed with MARS with the owner's permission.

Exodus II, 14,751 hexahedra, 30,000 nodes, one element thick in z (a quasi-2D channel):

```
x: 0 .. 10     streamwise, 150 node planes
y: 0 .. 1      wall-normal, 100 node planes, channel height H = 1
z: 0 .. 0.06   one element, 2 node planes
```

The solver finds the boundaries geometrically (inflow at x = xmin, outflow at x = xmax, walls
on the y faces, symmetry on the z faces), so the mesh needs no named side sets.
Reading it requires a MARS build with netCDF.
