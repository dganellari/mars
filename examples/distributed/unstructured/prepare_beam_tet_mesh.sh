#!/bin/bash
# Download and prepare beam-tet mesh for MARS

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MESH_DIR="${SCRIPT_DIR}/../../../meshes"

echo "=== MFEM beam-tet Mesh Preparation ==="
echo ""

# Create mesh directory if needed
mkdir -p "$MESH_DIR"

# Download beam-tet.mesh from MFEM if not present
if [ ! -f "$MESH_DIR/beam-tet.mesh" ]; then
    echo "Downloading beam-tet.mesh from MFEM repository..."
    curl -L -o "$MESH_DIR/beam-tet.mesh" \
        "https://raw.githubusercontent.com/mfem/mfem/master/data/beam-tet.mesh"
    echo "Downloaded beam-tet.mesh"
else
    echo "beam-tet.mesh already exists"
fi

# Compile converter if needed
if [ ! -f "$SCRIPT_DIR/mfem_to_mars_binary" ]; then
    echo ""
    echo "Compiling MFEM to MARS binary converter..."
    g++ -std=c++17 -O3 \
        -I"${SCRIPT_DIR}/../../../backend/distributed/unstructured/utils" \
        "${SCRIPT_DIR}/mfem_to_mars_binary.cpp" \
        -o "${SCRIPT_DIR}/mfem_to_mars_binary"
    echo "Compiled mfem_to_mars_binary"
fi

# Convert with 3 refinement levels (~24k elements)
echo ""
echo "Converting beam-tet.mesh to MARS binary format with 3 refinement levels..."
mkdir -p "$MESH_DIR/beam-tet-refined"
"${SCRIPT_DIR}/mfem_to_mars_binary" \
    "$MESH_DIR/beam-tet.mesh" \
    "$MESH_DIR/beam-tet-refined" \
    3

echo ""
echo "=== Preparation Complete ==="
echo ""
echo "To run the example:"
echo "  cd ${SCRIPT_DIR}/../../../"
echo "  mpirun -np 1 ./mars_ex_beam_tet --mesh meshes/beam-tet-refined"
echo ""
