#!/bin/bash
# Test the complete workflow

set -e

echo "=== Testing MFEM to MARS Binary Workflow ==="

# 1. Create a simple test mesh
cat > test_mesh.mesh << 'EOF'
MFEM mesh v1.0

dimension
3

elements
2
1 4 0 1 2 3
1 4 1 4 2 3

boundary
4
1 2 0 1 2
1 2 1 4 2
1 2 0 1 3
1 2 1 4 3

vertices
5
3
0 0 0
1 0 0
0 1 0
0 0 1
1 1 0
EOF

echo "Created test_mesh.mesh (2 tets, 5 vertices)"

# 2. Compile converter
echo ""
echo "Compiling converter..."
g++ -std=c++17 -O3 \
    -I../../../backend/distributed/unstructured/utils \
    mfem_to_mars_binary.cpp \
    -o test_converter

echo "Compiled successfully"

# 3. Convert without refinement
echo ""
echo "Testing conversion without refinement..."
mkdir -p test_output_noref
./test_converter test_mesh.mesh test_output_noref 0

# Verify files exist
for file in x.float32 y.float32 z.float32 i0.int32 i1.int32 i2.int32 i3.int32; do
    if [ ! -f "test_output_noref/$file" ]; then
        echo "ERROR: Missing file test_output_noref/$file"
        exit 1
    fi
done
echo "All files created successfully"

# 4. Convert with 1 refinement level
echo ""
echo "Testing conversion with 1 refinement..."
mkdir -p test_output_ref1
./test_converter test_mesh.mesh test_output_ref1 1

echo ""
echo "Checking refined mesh size..."
# Each tet splits into 8, so 2 → 16 elements
ELEM_COUNT=$(wc -c < test_output_ref1/i0.int32)
ELEM_COUNT=$((ELEM_COUNT / 4))  # int32 = 4 bytes
echo "  Elements: $ELEM_COUNT (expected 16)"

if [ "$ELEM_COUNT" -eq 16 ]; then
    echo "  ✓ Refinement worked correctly"
else
    echo "  ✗ Unexpected element count"
    exit 1
fi

# Clean up
echo ""
echo "Cleaning up test files..."
rm -f test_mesh.mesh test_converter
rm -rf test_output_noref test_output_ref1

echo ""
echo "=== All Tests Passed ==="
echo ""
echo "The converter is working correctly."
echo "Now run: ./prepare_beam_tet_mesh.sh"
