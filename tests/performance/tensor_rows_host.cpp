// Included after the kernels by test_cvfem_tensor_rows.py. Hex meshes with scrambled element order, DOF
// numbers and ownership, ghost rows, a node without a DOF and a periodic pair of corners sharing a DOF, and
// unstructured ones: hexes swept from a fan of quads around one node. The rows kernel and the one-thread
// tensor kernel, with and without the row map, must match tensor_faces entry by entry, in the graph and the
// full assembly, including the blocks and rows the map sends straight to the matrix.
using namespace mars::fem;
using Key = uint64_t;

struct Mesh
{
    size_t numElements = 0;
    std::vector<Key> conn[8];
    std::vector<double> x, y, z, gamma, phi, beta, gx, gy, gz, mdot, ax, ay, az;
    std::vector<int> nodeToDof, rowPtr, colInd, diagPtr;
    std::vector<uint8_t> own;
    int numTotal = 0;
};

unsigned next(unsigned& s) { return s = s * 1664525u + 1013904223u; }
template<class T>
void shuffle(std::vector<T>& v, unsigned seed)
{
    for (size_t i = v.size(); i > 1; --i) std::swap(v[i - 1], v[next(seed) % i]);
}

// the graph pattern (hex edges and the diagonal, as buildGraphSparsity) or the full one (every pair of an
// element's nodes, as buildFullSparsity), plus extra entries
void buildPattern(Mesh& m, bool full, const std::vector<std::pair<int, int>>& extra = {})
{
    static const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    std::vector<std::vector<int>> rows(m.numTotal);
    for (int r = 0; r < m.numTotal; ++r) rows[r].push_back(r);
    for (size_t e = 0; e < m.numElements; ++e)
        for (int a = 0; a < 8; ++a)
            for (int b = 0; b < 8; ++b)
            {
                const bool edge = std::any_of(std::begin(edges), std::end(edges), [&](const int* ed) {
                    return (ed[0] == a && ed[1] == b) || (ed[0] == b && ed[1] == a);
                });
                const int da = m.nodeToDof[m.conn[a][e]], db = m.nodeToDof[m.conn[b][e]];
                if (da >= 0 && db >= 0 && (full || edge)) rows[da].push_back(db);
            }
    for (auto [r, c] : extra) rows[r].push_back(c);
    m.rowPtr = {0}, m.colInd.clear(), m.diagPtr.assign(m.numTotal, -1);
    for (int r = 0; r < m.numTotal; ++r)
    {
        std::sort(rows[r].begin(), rows[r].end());
        rows[r].erase(std::unique(rows[r].begin(), rows[r].end()), rows[r].end());
        for (int c : rows[r])
        {
            if (c == r) m.diagPtr[r] = int(m.colInd.size());
            m.colInd.push_back(c);
        }
        m.rowPtr.push_back(int(m.colInd.size()));
    }
}

void fields(Mesh& m, int nodes)
{
    for (int n = 0; n < nodes; ++n)
    {
        m.phi.push_back(std::sin(.37 * n) + .2), m.gamma.push_back(1 + .3 * std::cos(.21 * n));
        m.beta.push_back(.5 + .25 * std::sin(.13 * n)), m.gx.push_back(std::sin(.11 * n + 1));
        m.gy.push_back(std::cos(.07 * n)), m.gz.push_back(std::sin(.05 * n + 2));
    }
    for (size_t f = 0; f < m.numElements * 12; ++f)
    {
        m.mdot.push_back(std::sin(.73 * f + .1)), m.ax.push_back(std::cos(.31 * f));
        m.ay.push_back(std::sin(.17 * f)), m.az.push_back(.5 + .3 * std::cos(.23 * f));
    }
}

// Scrambled DOF numbers, owned rows first, about one node in 5 a ghost; slave takes master's DOF (periodic)
// and noDof gets none (-1: no such node)
void numberDofs(Mesh& m, int nodes, unsigned seed, int master, int slave, int noDof)
{
    std::vector<int> owned, ghost;
    m.own.assign(nodes, 0);
    for (int n = 0; n < nodes; ++n)
    {
        if (n == slave || n == noDof) continue;
        if (n == master || (n * 7919u) % 5 != 0) owned.push_back(n), m.own[n] = 1;
        else ghost.push_back(n);
    }
    shuffle(owned, seed + 1), shuffle(ghost, seed + 2);
    m.nodeToDof.assign(nodes, -1);
    int dof = 0;
    for (int n : owned) m.nodeToDof[n] = dof++;
    for (int n : ghost) m.nodeToDof[n] = dof++;
    if (slave >= 0) m.nodeToDof[slave] = m.nodeToDof[master], m.own[slave] = 1;
    m.numTotal = dof;
}

// nx x ny x nz hexes, elements in scrambled order; corners (0,0,0) and (nx,0,0) share a DOF (periodic)
Mesh structured(int nx, int ny, int nz, unsigned seed)
{
    Mesh m;
    auto node = [&](int i, int j, int k) { return i + (nx + 1) * (j + (ny + 1) * k); };
    std::vector<int> order(nx * ny * nz);
    std::iota(order.begin(), order.end(), 0);
    shuffle(order, seed);
    static const int corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    m.numElements = order.size();
    for (int cell : order)
    {
        const int i = cell % nx, j = (cell / nx) % ny, k = cell / (nx * ny);
        for (int c = 0; c < 8; ++c) m.conn[c].push_back(Key(node(i + corner[c][0], j + corner[c][1], k + corner[c][2])));
    }
    const int nodes = (nx + 1) * (ny + 1) * (nz + 1);
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i)
            {
                m.x.push_back(i + .15 * std::sin(1.3 * j + .7 * k));
                m.y.push_back(j + .1 * std::cos(.9 * i + 1.1 * k));
                m.z.push_back(k + .12 * std::sin(.5 * i + .8 * j));
            }
    fields(m, nodes);

    numberDofs(m, nodes, seed, node(0, 0, 0), node(nx, 0, 0), node(nx / 2, ny / 2, nz / 2));
    return m;
}

// Hexes swept nz layers up from a fan of k quads around the node at the origin (centre, a_i, b_i, a_i+1):
// the swept centre nodes have k + 2 edge neighbours and up to 6k + 8 full-pattern neighbours.
Mesh fan(int k, int nz, unsigned seed)
{
    Mesh m;
    const int perLayer = 1 + 2 * k;
    auto a = [&](int i, int l) { return l * perLayer + 1 + (i % k); };
    auto b = [&](int i, int l) { return l * perLayer + 1 + k + i; };
    std::vector<int> order(k * nz);
    std::iota(order.begin(), order.end(), 0);
    shuffle(order, seed);
    m.numElements = order.size();
    for (int cell : order)
    {
        const int i = cell % k, l = cell / k;
        const int quad[2][4] = {{l * perLayer, a(i, l), b(i, l), a(i + 1, l)},
                                {(l + 1) * perLayer, a(i, l + 1), b(i, l + 1), a(i + 1, l + 1)}};
        for (int c = 0; c < 8; ++c) m.conn[c].push_back(Key(quad[c / 4][c % 4]));
    }
    const int nodes = perLayer * (nz + 1);
    const double pi = 3.14159265358979;
    for (int l = 0; l <= nz; ++l)
    {
        m.x.push_back(.05 * l), m.y.push_back(0), m.z.push_back(l);
        for (int i = 0; i < k; ++i)
            m.x.push_back(std::cos(2 * pi * i / k)), m.y.push_back(std::sin(2 * pi * i / k)), m.z.push_back(l + .1 * std::sin(i));
        for (int i = 0; i < k; ++i)
            m.x.push_back(1.6 * std::cos(2 * pi * (i + .5) / k)), m.y.push_back(1.6 * std::sin(2 * pi * (i + .5) / k)),
                m.z.push_back(l + .1 * std::cos(i));
    }
    fields(m, nodes);
    numberDofs(m, nodes, seed, -1, -1, -1);
    return m;
}

// 32 elements without a shared node: 256 rows in one block, more than the rows kernel takes
Mesh disjoint()
{
    Mesh m;
    m.numElements = 32;
    for (int e = 0; e < 32; ++e)
        for (int c = 0; c < 8; ++c) m.conn[c].push_back(Key(e * 8 + c));
    static const int corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int n = 0; n < 256; ++n)
    {
        const int e = n / 8, c = n % 8;
        m.x.push_back(2 * e + corner[c][0] + .1 * std::sin(n)), m.y.push_back(corner[c][1] + .1 * std::cos(n));
        m.z.push_back(corner[c][2] + .05 * std::sin(.7 * n)), m.nodeToDof.push_back(n);
    }
    m.own.assign(256, 1), m.numTotal = 256;
    fields(m, 256);
    return m;
}

struct System
{
    std::vector<double> values, rhs;
};

template<class Kernel>
System faces(Mesh& m, Kernel kernel, int block)
{
    System s{std::vector<double>(m.colInd.size(), 0), std::vector<double>(m.numTotal, 0)};
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), s.values.data(), m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    hostLaunch((m.numElements + block / 16 - 1) / (block / 16), size_t(block), kernel, m.conn[0].data(), m.conn[1].data(),
               m.conn[2].data(), m.conn[3].data(), m.conn[4].data(), m.conn[5].data(), m.conn[6].data(), m.conn[7].data(),
               m.numElements, m.x.data(), m.y.data(), m.z.data(), m.gamma.data(), m.phi.data(), m.beta.data(), m.gx.data(),
               m.gy.data(), m.gz.data(), m.mdot.data(), m.ax.data(), m.ay.data(), m.az.data(),
               (const int*)m.nodeToDof.data(), (const uint8_t*)m.own.data(), &A, s.rhs.data());
    return s;
}

// the one-thread tensor kernel; elemCols selects the row map scatter
template<class Kernel>
System tensor(Mesh& m, Kernel kernel, const uint8_t* elemCols)
{
    System s{std::vector<double>(m.colInd.size(), 0), std::vector<double>(m.numTotal, 0)};
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), s.values.data(), m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    hostLaunch((m.numElements + 255) / 256, 256, kernel, m.conn[0].data(), m.conn[1].data(), m.conn[2].data(),
               m.conn[3].data(), m.conn[4].data(), m.conn[5].data(), m.conn[6].data(), m.conn[7].data(), m.numElements,
               m.x.data(), m.y.data(), m.z.data(), m.gamma.data(), m.phi.data(), m.beta.data(), m.gx.data(),
               m.gy.data(), m.gz.data(), m.mdot.data(), m.ax.data(), m.ay.data(), m.az.data(),
               (const int*)m.nodeToDof.data(), (const uint8_t*)m.own.data(), &A, s.rhs.data(), elemCols);
    return s;
}

template<class Kernel>
System rows(Mesh& m, Kernel kernel, CvfemRowMap& map)
{
    System s{std::vector<double>(m.colInd.size(), 0), std::vector<double>(m.numTotal, 0)};
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), s.values.data(), m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    hostLaunch(size_t(map.numBlocks), 512, kernel, m.conn[0].data(),
               m.conn[1].data(), m.conn[2].data(), m.conn[3].data(), m.conn[4].data(), m.conn[5].data(), m.conn[6].data(),
               m.conn[7].data(), m.numElements, m.x.data(), m.y.data(), m.z.data(), m.gamma.data(), m.phi.data(),
               m.beta.data(), m.gx.data(), m.gy.data(), m.gz.data(), m.mdot.data(), m.ax.data(), m.ay.data(), m.az.data(),
               (const int*)m.nodeToDof.data(), &A, s.rhs.data(), map.view);
    return s;
}

CvfemRowMap& build(Mesh& m, bool full, CvfemRowMap& map)
{
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), nullptr, m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    buildCvfemRowMap<Key, double>(m.conn[0].data(), m.conn[1].data(), m.conn[2].data(), m.conn[3].data(),
                                  m.conn[4].data(), m.conn[5].data(), m.conn[6].data(), m.conn[7].data(), m.numElements,
                                  m.nodeToDof.data(), m.own.data(), A,
                                  full ? CvfemRowMapPattern::Full : CvfemRowMapPattern::Graph, map);
    return map;
}

bool same(const char* name, const System& a, const System& ref)
{
    double diff = 0, scale = 0;
    size_t nonzero = 0;
    for (size_t k = 0; k < ref.values.size(); ++k)
        diff = std::max(diff, std::abs(a.values[k] - ref.values[k])), scale = std::max(scale, std::abs(ref.values[k])),
        nonzero += ref.values[k] != 0;
    for (size_t k = 0; k < ref.rhs.size(); ++k)
        diff = std::max(diff, std::abs(a.rhs[k] - ref.rhs[k])), scale = std::max(scale, std::abs(ref.rhs[k]));
    const bool ok = nonzero > 0 && std::isfinite(scale) && diff <= 1e-12 * scale;
    std::printf("  %-36s max |diff| %.3e of %.3e (%zu nonzero entries): %s\n", name, diff, scale, nonzero, ok ? "PASS" : "FAIL");
    return ok;
}

// the map paths a check took: blocks sent straight to the matrix for too many rows or too many doubles, rows
// searched, shared-memory rows longer than the 8 doubles of a structured graph row
struct Paths
{
    bool direct = false, tooBig = false, search = false, wide = false, compact = false, bytes = false;
};

// every kernel against tensor_faces on one mesh and pattern; extraPerRow more columns in each row, which no
// element adds to, lengthen the rows
bool check(const char* name, Mesh& m, bool full, Paths& paths, int extraPerRow = 0)
{
    std::vector<std::pair<int, int>> extra;
    for (int r = 0; r < m.numTotal; ++r)
        for (int j = 0; j < extraPerRow; ++j) extra.push_back({r, (r + 1 + 3 * j) % m.numTotal});
    buildPattern(m, full, extra);
    int longest = 0;
    for (int r = 0; r < m.numTotal; ++r) longest = std::max(longest, m.rowPtr[r + 1] - m.rowPtr[r]);
    std::printf("%s, %s pattern (%zu elements, %d DOFs, longest row %d)\n", name, full ? "full" : "graph",
                m.numElements, m.numTotal, longest);
    bool ok = true;
    System ref;
    if (full)
    {
        ref = faces(m, cvfem_hex_assembly_kernel_tensor_faces<Key, double, 256, false, false>, 256);
        ok  = same("tensor", tensor(m, cvfem_hex_assembly_kernel_tensor<Key, double, 256, false>, nullptr), ref) && ok;
    }
    else
    {
        ref = faces(m, cvfem_hex_assembly_kernel_tensor_faces<Key, double, 256, true, false>, 256);
        ok  = same("tensor_faces aggregated",
                   faces(m, cvfem_hex_assembly_kernel_tensor_faces<Key, double, 512, true, true>, 512), ref) && ok;
        ok  = same("tensor", tensor(m, cvfem_hex_assembly_kernel_tensor<Key, double, 256, true>, nullptr), ref) && ok;
    }
    CvfemRowMap map;
    build(m, full, map);
    int widest = 0;
    for (int b = 0; b < map.numBlocks; ++b)
        for (int k = map.view.blockRowStart[b]; k < map.view.blockRowStart[b + 1]; ++k)
        {
            const int end = k + 1 < map.view.blockRowStart[b + 1] ? map.view.rowBase[k + 1] : map.view.blockSize[b];
            widest = std::max(widest, end - int(map.view.rowBase[k]));
        }
    std::printf("  row map: %d blocks (%d direct), %d block rows, widest %d doubles, %d searched rows, %s columns\n",
                map.numBlocks, map.numDirectBlocks, map.numBlockRows, widest, map.numSearchRows,
                map.compact ? "4-bit" : "byte");
    paths.direct |= map.numDirectBlocks > 0, paths.search |= map.numSearchRows > 0, paths.wide |= widest > 8;
    paths.tooBig |= extraPerRow > 0 && map.numDirectBlocks > 0;
    paths.compact |= map.compact, paths.bytes |= !map.compact;
    // the kernels must match the map's format, as the assembler picks them
    const auto mapKernel = full ? (map.compact ? cvfem_hex_assembly_kernel_tensor<Key, double, 256, false, true, true>
                                               : cvfem_hex_assembly_kernel_tensor<Key, double, 256, false, true, false>)
                                : (map.compact ? cvfem_hex_assembly_kernel_tensor<Key, double, 256, true, true, true>
                                               : cvfem_hex_assembly_kernel_tensor<Key, double, 256, true, true, false>);
    const System mapped  = tensor(m, mapKernel, map.view.elemCols);
    ok = same("tensor_map", mapped, ref) && ok;
    const auto rowsKernel = map.compact ? cvfem_hex_assembly_kernel_tensor_rows<Key, double, true>
                                        : cvfem_hex_assembly_kernel_tensor_rows<Key, double, false>;
    ok = same("tensor_rows", rows(m, rowsKernel, map), ref) && ok;
    return ok;
}

int main()
{
    bool ok = true;
    Paths paths;
    for (bool full : {false, true})
    {
        Mesh a = structured(7, 5, 3, 12345u + 7), b = structured(4, 4, 4, 12345u + 4);
        ok = check("7x5x3 hexes", a, full, paths) && ok;
        ok = check("4x4x4 hexes", b, full, paths) && ok;
        // valence 8, as at the irregular nodes of unstructured hex meshes: graph rows of 11 entries
        Mesh f8 = fan(8, 3, 99u);
        ok = check("fan of 8 quads, 3 layers", f8, full, paths) && ok;
        // valence 48: graph rows of 51 entries, full rows of up to 3 * 97
        Mesh f48 = fan(48, 3, 7u);
        ok = check("fan of 48 quads, 3 layers", f48, full, paths) && ok;
        // valence 253, one layer: graph rows of 256 entries, so the lumped graph rows are searched too
        Mesh f253 = fan(253, 1, 5u);
        ok = check("fan of 253 quads, 1 layer", f253, full, paths) && ok;
        Mesh d = disjoint();
        ok = check("32 disjoint hexes", d, full, paths) && ok;
        // about 88 rows per block of 48 doubles or more: more than the shared memory holds
        Mesh w = structured(4, 4, 4, 12345u + 40);
        ok = check("4x4x4 hexes, 40 more columns per row", w, full, paths, 40) && ok;
    }
    std::printf("paths taken: direct blocks %s, blocks too big for shared memory %s, searched rows %s, wide "
                "shared-memory rows %s, 4-bit maps %s, byte maps %s\n",
                paths.direct ? "yes" : "NO", paths.tooBig ? "yes" : "NO", paths.search ? "yes" : "NO",
                paths.wide ? "yes" : "NO", paths.compact ? "yes" : "NO", paths.bytes ? "yes" : "NO");
    ok = ok && paths.direct && paths.tooBig && paths.search && paths.wide && paths.compact && paths.bytes;
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
