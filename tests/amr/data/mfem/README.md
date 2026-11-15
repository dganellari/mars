# Public Escher tetrahedra

`escher.mesh` is the unchanged linear Tet4 fixture from
[MFEM v4.8](https://github.com/mfem/mfem/blob/v4.8/data/escher.mesh), verified against
the upstream file. SHA-256:
`9b57c9d3bf3247fd6646a5b726d199026106cd4e773f5a466623f81a1a0d97c9`.
The accompanying BSD 3-Clause license is retained.

It has 42 tetrahedra, 26 vertices and 48 triangular boundary labels. The test
parser normalizes signed orientation; it does not curve or remesh the geometry.
Host and native tests refine three times, checking positive volume, interior-face
pairing, boundary-label area and affine/constant transfer. Native GPU results are
required separately from the host reference. This is one non-box geometry, not
an arbitrary-mesh guarantee.

## Refined Exodus fixture

`escher-refined.e` is the public fixture used by the v0.2.0 AMR/PDE package
validation. It is the original Escher mesh uniformly red-refined twice, with
orientation normalized and boundary labels retained, then converted to Exodus.
It contains 2688 tetrahedra and 665 nodes. The accompanying MFEM license applies
to this derived fixture as well. SHA-256:
`1d9a26dac3da46f782a945ce528107291459fd2f5878421deb29ff3da585d0d5`.

This file is byte-identical to the input used in the accepted NVIDIA and AMD
release package campaigns. See [the walkthrough](../../../../docs/AMR-Peak-Example.md#reproduce-the-v020-example).
