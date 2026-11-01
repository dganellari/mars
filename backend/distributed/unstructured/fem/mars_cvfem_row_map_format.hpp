#pragma once

#include <cstdint>

namespace mars {
namespace fem {

// The row map (mars_cvfem_row_map.hpp) as the kernels read it. Block b, 32 elements, adds up in shared
// memory the rows rowDof[blockRowStart[b] .. blockRowStart[b + 1]), in increasing DOF order. Each row takes
// its own length, so that the long rows of unstructured meshes do not widen the short ones: row k starts at
// rowBase[k] with the right-hand side, then its CSR row; blockSize[b] doubles in all, 0 for a block that
// does not fit. elemRow gives each element node's place in the block's rows; elemCols, for each element row,
// the CSR offset in the row that each of the element's 8 columns adds to, a byte each (4 bits in a compact
// map, below).
// lumped: a column outside the pattern went to the diagonal (graph assembly) instead of being dropped.
struct TensorRowsMap
{
    const uint8_t* elemRow;
    const uint8_t* elemCols;
    const int* blockRowStart;
    const uint16_t* blockSize;
    const uint16_t* rowBase;
    const int* rowDof;
    bool lumped;
};

// elemRow: the row is not assembled here. elemCols: the column is dropped (no DOF, or outside the pattern of
// the full assembly), and all 8 bytes for a row not assembled here.
constexpr uint8_t rowMapNone = 0xFF;
// elemRow: the row is added straight to the matrix, because its block's rows do not fit in shared memory or
// the row is too long for elemCols
constexpr uint8_t rowMapDirect = 0xFE;
// elemCols, all 8 bytes: the row is longer than rowMapMaxLen, and the kernels search its CSR row
constexpr uint8_t rowMapSearch = 0xFE;
constexpr int rowMapMaxLen     = 0xFE;
// A compact map (CvfemRowMap::compact) has no row longer than rowMapCompactMaxLen, no direct block and no
// searched row: every graph map so far. Its elemCols take 4 bits per column, one word per element row,
// 32 bytes per element instead of 64, and the kernels for it leave out the direct and search paths.
constexpr unsigned rowMapCompactNone = 0xF;
constexpr int rowMapCompactMaxLen    = 15;

} // namespace fem
} // namespace mars
