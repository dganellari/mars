#!/usr/bin/env python3
"""
Convert Exodus mesh to MFEM format for MARS
Usage: python exodus_to_mfem.py input.exo output.mesh
"""

import sys
import numpy as np

try:
    import meshio
except ImportError:
    print("Error: meshio not installed. Install with: pip install meshio")
    sys.exit(1)

def exodus_to_mfem(exo_file, mfem_file):
    """Convert Exodus II mesh to MFEM format"""
    
    print(f"Reading Exodus mesh: {exo_file}")
    mesh = meshio.read(exo_file)
    
    print(f"Mesh info:")
    print(f"  Points: {len(mesh.points)}")
    print(f"  Cell types: {list(mesh.cells_dict.keys())}")
    print(f"  Cell blocks: {len(mesh.cells)}")
    
    # Show field data (parts/blocks)
    if mesh.field_data:
        print(f"  Parts/Blocks:")
        for name, data in mesh.field_data.items():
            print(f"    {name}: dimension={data[0]}, tag={data[1]}")
    
    # Get hexahedral elements (hexahedron in meshio)
    if 'hexahedron' not in mesh.cells_dict:
        print("Error: No hexahedral elements found in mesh")
        sys.exit(1)
    
    hex_cells = mesh.cells_dict['hexahedron']
    print(f"  Hexahedra: {len(hex_cells)}")
    
    # Get boundary quadrilateral faces if present
    boundary_quads = mesh.cells_dict.get('quad', np.array([]))
    print(f"  Boundary quads: {len(boundary_quads)}")
    
    # Write MFEM format
    print(f"Writing MFEM mesh: {mfem_file}")
    with open(mfem_file, 'w') as f:
        # Header
        f.write("MFEM mesh v1.0\n\n")
        
        # Dimension
        f.write("dimension\n")
        f.write("3\n\n")
        
        # Elements (hexahedra)
        f.write("elements\n")
        f.write(f"{len(hex_cells)}\n")
        for i, cell in enumerate(hex_cells):
            # MFEM format: attribute element_type vertex_indices
            # element_type 5 = hexahedron
            f.write(f"1 5 {' '.join(map(str, cell))}\n")
        f.write("\n")
        
        # Boundary elements (quads)
        f.write("boundary\n")
        f.write(f"{len(boundary_quads)}\n")
        for i, cell in enumerate(boundary_quads):
            # MFEM format: attribute element_type vertex_indices
            # element_type 3 = quadrilateral
            # Default boundary attribute = 1
            f.write(f"1 3 {' '.join(map(str, cell))}\n")
        f.write("\n")
        
        # Vertices
        f.write("vertices\n")
        f.write(f"{len(mesh.points)}\n")
        f.write("3\n")  # 3D coordinates
        for point in mesh.points:
            f.write(f"{point[0]:.16e} {point[1]:.16e} {point[2]:.16e}\n")
    
    print(f"Conversion complete!")
    print(f"MFEM mesh written to: {mfem_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python exodus_to_mfem.py input.exo output.mesh")
        sys.exit(1)
    
    exodus_to_mfem(sys.argv[1], sys.argv[2])
