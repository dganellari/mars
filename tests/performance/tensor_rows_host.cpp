// Included after the kernels by test_cvfem_tensor_rows.py. Hex meshes with scrambled element order, DOF
// numbers and ownership, ghost rows, a node without a DOF and a periodic pair of corners sharing a DOF: the
// rows kernel and the one-thread tensor kernel, with and without the row map, must match tensor_faces
// (direct and aggregated) entry by entry.
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

// the graph pattern: hex edges and the diagonal, as buildGraphSparsity
void buildPattern(Mesh& m, const std::vector<std::pair<int, int>>& extra)
{
    static const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    std::vector<std::vector<int>> rows(m.numTotal);
    for (int r = 0; r < m.numTotal; ++r) rows[r].push_back(r);
    for (size_t e = 0; e < m.numElements; ++e)
        for (auto& ed : edges)
        {
            const int a = m.nodeToDof[m.conn[ed[0]][e]], b = m.nodeToDof[m.conn[ed[1]][e]];
            if (a >= 0 && b >= 0) rows[a].push_back(b), rows[b].push_back(a);
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

    const int master = node(0, 0, 0), slave = node(nx, 0, 0), noDof = node(nx / 2, ny / 2, nz / 2);
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
    m.nodeToDof[slave] = m.nodeToDof[master], m.own[slave] = 1;
    m.numTotal = dof;
    buildPattern(m, {});
    return m;
}

// 32 elements without a shared node: 256 rows in one block, one more than the rows kernels take
Mesh disjoint()
{
    Mesh m;
    m.numElements = 32;
    for (int e = 0; e < 32; ++e)
        for (int c = 0; c < 8; ++c) m.conn[c].push_back(Key(e * 8 + c));
    for (int n = 0; n < 256; ++n) m.x.push_back(n), m.y.push_back(n % 7), m.z.push_back(n % 3), m.nodeToDof.push_back(n);
    m.own.assign(256, 1), m.numTotal = 256;
    fields(m, 256);
    buildPattern(m, {});
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
System tensor(Mesh& m, Kernel kernel, const unsigned* elemCols)
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

System rows(Mesh& m, CvfemRowMap& map)
{
    System s{std::vector<double>(m.colInd.size(), 0), std::vector<double>(m.numTotal, 0)};
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), s.values.data(), m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    hostLaunch(size_t(map.numBlocks), 512, cvfem_hex_assembly_kernel_tensor_rows<Key, double>, m.conn[0].data(),
               m.conn[1].data(), m.conn[2].data(), m.conn[3].data(), m.conn[4].data(), m.conn[5].data(), m.conn[6].data(),
               m.conn[7].data(), m.numElements, m.x.data(), m.y.data(), m.z.data(), m.gamma.data(), m.phi.data(),
               m.beta.data(), m.gx.data(), m.gy.data(), m.gz.data(), m.mdot.data(), m.ax.data(), m.ay.data(), m.az.data(),
               &A, s.rhs.data(), map.view);
    return s;
}

CvfemRowMap& build(Mesh& m, CvfemRowMap& map)
{
    CSRMatrix<double> A{m.rowPtr.data(), m.colInd.data(), nullptr, m.diagPtr.data(), m.numTotal, int(m.colInd.size())};
    buildCvfemRowMap<Key, double>(m.conn[0].data(), m.conn[1].data(), m.conn[2].data(), m.conn[3].data(),
                                  m.conn[4].data(), m.conn[5].data(), m.conn[6].data(), m.conn[7].data(), m.numElements,
                                  m.nodeToDof.data(), m.own.data(), A, map);
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
    const bool ok = nonzero > 0 && diff <= 1e-12 * scale;
    std::printf("  %-36s max |diff| %.3e of %.3e (%zu nonzero entries): %s\n", name, diff, scale, nonzero, ok ? "PASS" : "FAIL");
    return ok;
}

template<class F>
bool throws(const char* name, F f)
{
    try { f(); }
    catch (const std::runtime_error& e)
    {
        std::printf("  %-36s throws \"%s\": PASS\n", name, e.what());
        return true;
    }
    std::printf("  %-36s does not throw: FAIL\n", name);
    return false;
}

int main()
{
    bool ok = true;
    for (auto [nx, ny, nz] : {std::array<int, 3>{7, 5, 3}, std::array<int, 3>{4, 4, 4}})
    {
        Mesh m = structured(nx, ny, nz, 12345u + nx);
        std::printf("%dx%dx%d hexes (%zu elements, %d DOFs)\n", nx, ny, nz, m.numElements, m.numTotal);
        const System ref = faces(m, cvfem_hex_assembly_kernel_tensor_faces<Key, double, 256, true, false>, 256);
        ok = same("tensor_faces aggregated", faces(m, cvfem_hex_assembly_kernel_tensor_faces<Key, double, 512, true, true>, 512), ref) && ok;
        CvfemRowMap map;
        build(m, map);
        std::printf("  row map: %d blocks, %d block rows\n", map.numBlocks, map.numBlockRows);
        ok = same("tensor", tensor(m, cvfem_hex_assembly_kernel_tensor<Key, double, 256, true>, nullptr), ref) && ok;
        ok = same("tensor_map", tensor(m, cvfem_hex_assembly_kernel_tensor<Key, double, 256, true, true>, map.view.elemCols), ref) && ok;
        ok = same("tensor_rows", rows(m, map), ref) && ok;
    }
    std::printf("limits\n");
    Mesh d = disjoint();
    ok = throws("256 rows in a block", [&] { CvfemRowMap map; build(d, map); }) && ok;
    // the periodic corner's row (owned, 7 entries) gains 3 columns of the far corner: 10 entries
    Mesh m = structured(3, 3, 2, 7u);
    auto dofAt = [&](int i, int j, int k) { return m.nodeToDof[i + 4 * (j + 4 * k)]; };
    const int row = dofAt(0, 0, 0);
    buildPattern(m, {{row, dofAt(3, 3, 2)}, {row, dofAt(2, 3, 2)}, {row, dofAt(3, 2, 2)}});
    ok = throws("a CSR row of more than 8 entries", [&] { CvfemRowMap map; build(m, map); }) && ok;
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
