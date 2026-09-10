// GPT/Codex, 2026-09-09. Production CUDA scatters on public one/two-tet fixtures.
// MPI reductions below deliberately replace the domain halo; this tests kernels
// and owner-row indexing, not ElementDomain communication or a flow solution.
#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"
#include "mars_outlet_boundary_gate.hpp"

namespace {

using mars::outlet_gate::Checks;
using KeyType = unsigned;
using RealType = double;
using Field = std::vector<double>;
using VectorField = std::array<Field, 3>;

void cuda_check(cudaError_t error, const char* operation)
{
    if (error == cudaSuccess) return;
    std::fprintf(stderr, "FAIL %s: %s\n", operation, cudaGetErrorString(error));
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::abort();
}

void complete(const char* operation)
{
    cuda_check(cudaGetLastError(), operation);
    cuda_check(cudaDeviceSynchronize(), operation);
}

template<class T> struct DeviceArray {
    T* data = nullptr;
    std::size_t size = 0;
    DeviceArray() = default;
    DeviceArray(const DeviceArray&) = delete;
    DeviceArray& operator=(const DeviceArray&) = delete;
    ~DeviceArray() { if (data) cuda_check(cudaFree(data), "free fixture"); }
    void resize(std::size_t requested)
    {
        if (requested == size) return;
        if (data) cuda_check(cudaFree(data), "resize fixture");
        data = nullptr;
        size = requested;
        if (size) cuda_check(cudaMalloc(&data, size*sizeof(T)), "allocate fixture");
    }
    void upload(const std::vector<T>& source)
    {
        resize(source.size());
        if (size) cuda_check(cudaMemcpy(data, source.data(), size*sizeof(T), cudaMemcpyHostToDevice), "publish fixture");
    }
    void zero(std::size_t requested)
    {
        resize(requested);
        if (size) cuda_check(cudaMemset(data, 0, size*sizeof(T)), "clear fixture");
    }
    std::vector<T> download() const
    {
        std::vector<T> result(size);
        if (size) cuda_check(cudaMemcpy(result.data(), data, size*sizeof(T), cudaMemcpyDeviceToHost), "read fixture");
        return result;
    }
};

Field sum_ranks(const Field& local)
{
    Field result(local.size());
    MPI_Allreduce(local.data(), result.data(), static_cast<int>(local.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

struct Facet {
    int element, opposite;
    std::array<int, 3> nodes;
    std::array<double, 3> area;
    uint8_t outlet;
};

struct DeviceFacets {
    DeviceArray<int> nodes, elements, opposites;
    DeviceArray<uint8_t> outlets;
    std::array<DeviceArray<double>, 3> area;
    int count = 0;
    void upload(const std::vector<Facet>& facets)
    {
        std::vector<int> node_values, element_values, opposite_values;
        std::vector<uint8_t> outlet_values;
        VectorField areas;
        for (const auto& facet : facets) {
            node_values.insert(node_values.end(), facet.nodes.begin(), facet.nodes.end());
            element_values.push_back(facet.element);
            opposite_values.push_back(facet.opposite);
            outlet_values.push_back(facet.outlet);
            for (int d = 0; d < 3; ++d) areas[d].push_back(facet.area[d]);
        }
        nodes.upload(node_values); elements.upload(element_values); opposites.upload(opposite_values);
        outlets.upload(outlet_values);
        for (int d = 0; d < 3; ++d) area[d].upload(areas[d]);
        count = static_cast<int>(facets.size());
    }
};

struct Fixture {
    int elements, nodes;
    std::array<std::array<int,4>,2> connectivity{{{0,1,2,3},{4,1,3,2}}};
    const double coordinates[5][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
    const double gradients[2][4][3] = {
        {{-1,-1,-1},{1,0,0},{0,1,0},{0,0,1}},
        {{.5,.5,.5},{.5,-.5,-.5},{-.5,-.5,.5},{-.5,.5,-.5}}
    };
    std::vector<Facet> facets;
    Field mass, coefficient, pressure, trace;
    VectorField velocity, reconstructed, scs_area, outlet_area;
    std::vector<uint8_t> fixed;

    explicit Fixture(int num_elements) : elements(num_elements), nodes(num_elements == 1 ? 4 : 5),
        mass(nodes), coefficient(nodes), pressure(nodes), trace(nodes), fixed(nodes)
    {
        for (int d = 0; d < 3; ++d) {
            velocity[d].resize(nodes); reconstructed[d].resize(nodes);
            scs_area[d].resize(6*elements); outlet_area[d].resize(nodes);
        }
        for (int node = 0; node < nodes; ++node) {
            coefficient[node] = .2+.05*node;
            pressure[node] = .3*(node%3)-.1*node;
            trace[node] = .7-.11*node;
            fixed[node] = node == 0 || node == 2 || node == 3;
            for (int d = 0; d < 3; ++d) {
                velocity[d][node] = .13*(node+1)-.07*d;
                reconstructed[d][node] = .21*(d+1)-.09*node;
            }
        }
        for (int element = 0; element < elements; ++element) {
            const double volume = element == 0 ? 1./6 : 1./3;
            for (int local = 0; local < 4; ++local) mass[connectivity[element][local]] += volume/4;
            for (int edge = 0; edge < 6; ++edge) {
                const int left = mars::outlet_gate::edges[edge][0], right = mars::outlet_gate::edges[edge][1];
                for (int d = 0; d < 3; ++d)
                    scs_area[d][6*element+edge] = volume*(gradients[element][right][d]-gradients[element][left][d])/4;
            }
        }
        add_facet(0, 3, true);
        add_facet(0, 1, false);
        if (elements == 2) add_facet(1, 1, true);
    }

    void add_facet(int element, int opposite_local, bool outlet)
    {
        Facet facet{};
        facet.element = element;
        facet.opposite = connectivity[element][opposite_local];
        facet.outlet = outlet;
        int r = 0;
        for (int local = 0; local < 4; ++local)
            if (local != opposite_local) facet.nodes[r++] = connectivity[element][local];
        double edge1[3], edge2[3], inward[3];
        for (int d = 0; d < 3; ++d) {
            edge1[d] = coordinates[facet.nodes[1]][d]-coordinates[facet.nodes[0]][d];
            edge2[d] = coordinates[facet.nodes[2]][d]-coordinates[facet.nodes[0]][d];
            inward[d] = coordinates[facet.opposite][d]-coordinates[facet.nodes[0]][d];
        }
        mars::outlet_gate::cross(edge1, edge2, facet.area.data());
        const double scale = mars::outlet_gate::dot(facet.area.data(), inward) > 0 ? -.5 : .5;
        for (int d = 0; d < 3; ++d) {
            facet.area[d] *= scale;
            if (outlet) for (int node : facet.nodes) outlet_area[d][node] += facet.area[d]/3;
        }
        facets.push_back(facet);
    }

    Field reference_residual(const Field& p, const VectorField& u, bool nodal,
                             double* exterior = nullptr) const
    {
        Field result(nodes);
        for (int element = 0; element < elements; ++element) {
            double compact[3]{};
            for (int local = 0; local < 4; ++local)
                for (int d = 0; d < 3; ++d) compact[d] += p[connectivity[element][local]]*gradients[element][local][d];
            for (int edge = 0; edge < 6; ++edge) {
                const int left = connectivity[element][mars::outlet_gate::edges[edge][0]];
                const int right = connectivity[element][mars::outlet_gate::edges[edge][1]];
                const double diffusion = nodal ? .5*(coefficient[left]+coefficient[right]) : .3;
                double flux = 0;
                for (int d = 0; d < 3; ++d) {
                    const double smooth = .5*((fixed[left] ? 0 : reconstructed[d][left])+
                                              (fixed[right] ? 0 : reconstructed[d][right]));
                    flux += (.5*(u[d][left]+u[d][right])+diffusion*(smooth-compact[d]))*scs_area[d][6*element+edge];
                }
                result[left] += flux; result[right] -= flux;
            }
        }
        double total = 0;
        for (const auto& facet : facets) {
            double stabilization = 0;
            if (facet.outlet) {
                double diffusion = 0;
                for (int node : facet.nodes) diffusion += nodal ? coefficient[node]/3 : .1;
                for (int d = 0; d < 3; ++d) {
                    double compact = 0, face_gradient = 0;
                    for (int local = 0; local < 4; ++local) {
                        const int node = connectivity[facet.element][local];
                        compact += (node == facet.opposite ? p[node] : trace[node])*gradients[facet.element][local][d];
                    }
                    for (int node : facet.nodes) face_gradient += reconstructed[d][node]/3;
                    stabilization += diffusion*(.5*(face_gradient+reconstructed[d][facet.opposite])-compact)*facet.area[d]/3;
                }
            }
            for (int node : facet.nodes) {
                double flux = stabilization;
                for (int d = 0; d < 3; ++d) flux += u[d][node]*facet.area[d]/3;
                result[node] += flux; total += flux;
            }
        }
        if (exterior) *exterior = total;
        return result;
    }

    std::vector<Field> divergence() const
    {
        std::vector<Field> result(nodes, Field(3*nodes));
        for (int element = 0; element < elements; ++element)
            for (int edge = 0; edge < 6; ++edge) {
                const int left = connectivity[element][mars::outlet_gate::edges[edge][0]];
                const int right = connectivity[element][mars::outlet_gate::edges[edge][1]];
                for (int d = 0; d < 3; ++d) {
                    const double value = .5*scs_area[d][6*element+edge];
                    result[left][3*left+d] += value; result[left][3*right+d] += value;
                    result[right][3*left+d] -= value; result[right][3*right+d] -= value;
                }
            }
        for (const auto& facet : facets) if (facet.outlet)
            for (int node : facet.nodes)
                for (int d = 0; d < 3; ++d) result[node][3*node+d] += facet.area[d]/3;
        return result;
    }
};

struct DeviceFixture {
    const Fixture& fixture;
    int rank, ranks, owned_count = 0;
    std::array<DeviceArray<KeyType>, 4> connectivity;
    std::array<DeviceArray<double>, 3> coordinates, area, velocity, reconstructed, outlet_area;
    DeviceArray<double> pressure, trace, coefficient, mass;
    DeviceArray<int> node_to_dof, row_ptr, col_ind;
    DeviceArray<uint8_t> ownership, fixed_node, fixed_dof;
    DeviceFacets local_facets, all_outlets;
    std::vector<int> dof_to_node;

    DeviceFixture(const Fixture& input, int rank_in, int ranks_in, Checks& checks)
        : fixture(input), rank(rank_in), ranks(ranks_in)
    {
        const int n = fixture.nodes;
        std::vector<uint8_t> owner(n), fixed(n);
        for (int node = 0; node < n; ++node) if (node%ranks == rank) {
            owner[node] = 1; dof_to_node.push_back(node); ++owned_count;
        }
        for (int node = 0; node < n; ++node) if (!owner[node]) dof_to_node.push_back(node);
        std::vector<int> map(n), rows(owned_count+1), columns(owned_count*n);
        for (int dof = 0; dof < n; ++dof) {
            map[dof_to_node[dof]] = dof;
            fixed[dof] = fixture.fixed[dof_to_node[dof]];
        }
        for (int row = 0; row <= owned_count; ++row) rows[row] = row*n;
        for (int index = 0; index < owned_count*n; ++index) columns[index] = index%n;
        node_to_dof.upload(map); row_ptr.upload(rows); col_ind.upload(columns);
        ownership.upload(owner); fixed_node.upload(fixture.fixed); fixed_dof.upload(fixed);
        for (int local = 0; local < 4; ++local) {
            std::vector<KeyType> values(fixture.elements);
            for (int element = 0; element < fixture.elements; ++element) values[element] = fixture.connectivity[element][local];
            connectivity[local].upload(values);
        }
        for (int d = 0; d < 3; ++d) {
            Field values(n);
            for (int node = 0; node < n; ++node) values[node] = fixture.coordinates[node][d];
            coordinates[d].upload(values); velocity[d].upload(fixture.velocity[d]);
            reconstructed[d].upload(fixture.reconstructed[d]); outlet_area[d].upload(fixture.outlet_area[d]);
            area[d].zero(6*fixture.elements);
        }
        pressure.upload(fixture.pressure); trace.upload(fixture.trace);
        coefficient.upload(fixture.coefficient); mass.upload(fixture.mass);
        precomputeTetAreaVectorsKernel<KeyType,double><<<1,128>>>(
            connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
            fixture.elements, coordinates[0].data, coordinates[1].data, coordinates[2].data,
            area[0].data, area[1].data, area[2].data);
        complete("production tet SCS geometry");
        for (int d = 0; d < 3; ++d) {
            const auto values = area[d].download();
            for (std::size_t index = 0; index < values.size(); ++index)
                checks.near(values[index], fixture.scs_area[d][index], 1e-13, "CUDA SCS area versus independent simplex formula");
        }
        std::vector<Facet> local, outlets;
        for (const auto& facet : fixture.facets) {
            if (facet.element%ranks == rank) local.push_back(facet);
            if (facet.outlet) outlets.push_back(facet);
        }
        local_facets.upload(local); all_outlets.upload(outlets);
    }

    Field residual(const Field& p, const VectorField& u, bool nodal, double* exterior = nullptr)
    {
        const int n = fixture.nodes;
        pressure.upload(p);
        for (int d = 0; d < 3; ++d) velocity[d].upload(u[d]);
        DeviceArray<double> accumulator;
        accumulator.zero(n);
        const double* diffusion = nodal ? coefficient.data : nullptr;
        for (int element = 0; element < fixture.elements; ++element) if (element%ranks == rank) {
            computeDivergenceVMSTetKernel<KeyType,double><<<1,128>>>(
                connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
                velocity[0].data, velocity[1].data, velocity[2].data, pressure.data,
                reconstructed[0].data, reconstructed[1].data, reconstructed[2].data,
                coordinates[0].data, coordinates[1].data, coordinates[2].data,
                area[0].data, area[1].data, area[2].data, .3, diffusion, true, 1., nullptr,
                fixed_node.data, accumulator.data, element, 1);
            complete("production VMS interior scatter");
        }
        const auto& f = local_facets;
        double local_exterior = 0;
        if (f.count) {
            scatterOutletContinuityKernel<KeyType,double><<<1,128>>>(
                connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
                f.nodes.data, f.elements.data, f.opposites.data, f.outlets.data,
                f.area[0].data, f.area[1].data, f.area[2].data,
                velocity[0].data, velocity[1].data, velocity[2].data, pressure.data,
                reconstructed[0].data, reconstructed[1].data, reconstructed[2].data,
                coordinates[0].data, coordinates[1].data, coordinates[2].data,
                diffusion, .3, trace.data, 0., true, accumulator.data, f.count);
            complete("production outlet/inlet sample scatter");
            if (exterior) {
                DeviceArray<double> partial_in, partial_out;
                partial_in.zero(1); partial_out.zero(1);
                boundaryMassFluxKernel<KeyType,double><<<1,128>>>(
                    connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
                    f.nodes.data, f.elements.data, f.opposites.data, f.outlets.data,
                    f.area[0].data, f.area[1].data, f.area[2].data,
                    velocity[0].data, velocity[1].data, velocity[2].data, pressure.data,
                    reconstructed[0].data, reconstructed[1].data, reconstructed[2].data,
                    coordinates[0].data, coordinates[1].data, coordinates[2].data,
                    diffusion, .3, trace.data, 0., true, partial_in.data, partial_out.data, f.count);
                complete("production boundary reporting");
                local_exterior = partial_in.download()[0]+partial_out.download()[0];
            }
        }
        const Field result = sum_ranks(accumulator.download());
        if (exterior) MPI_Allreduce(&local_exterior, exterior, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return result;
    }

    Field compact_matrix(double rho_over_dt, bool nodal)
    {
        const int n = fixture.nodes;
        DeviceArray<double> matrix;
        matrix.zero(owned_count*n);
        const double* diffusion = nodal ? coefficient.data : nullptr;
        // Row owners need contributions from all incident cells, including halo cells.
        assembleRhieChowSensitivityTetKernel<KeyType,double><<<1,128>>>(
            connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
            coordinates[0].data, coordinates[1].data, coordinates[2].data,
            area[0].data, area[1].data, area[2].data, diffusion, .3, rho_over_dt,
            node_to_dof.data, ownership.data, row_ptr.data, col_ind.data, owned_count,
            matrix.data, 0, fixture.elements);
        complete("production interior compact derivative");
        const auto& f = all_outlets;
        DeviceArray<int> invalid_anchor;
        invalid_anchor.zero(1);
        assemble_outlet_pressure_sensitivity_kernel<KeyType,double><<<1,128>>>(
            connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
            f.nodes.data, f.elements.data, f.opposites.data,
            f.area[0].data, f.area[1].data, f.area[2].data,
            coordinates[0].data, coordinates[1].data, coordinates[2].data,
            diffusion, .3, rho_over_dt, node_to_dof.data, ownership.data,
            row_ptr.data, col_ind.data, owned_count, n, matrix.data, invalid_anchor.data, f.count);
        complete("production outlet compact derivative");
        int local_invalid = invalid_anchor.download()[0], any_invalid = 0;
        MPI_Allreduce(&local_invalid, &any_invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_invalid) {
            std::fprintf(stderr, "FAIL complete fixture rejected by outlet anchor check\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        const auto local_matrix = matrix.download();
        Field global_indices(n*n);
        for (int row = 0; row < owned_count; ++row)
            for (int col = 0; col < n; ++col)
                global_indices[n*dof_to_node[row]+dof_to_node[col]] = local_matrix[row*n+col];
        return sum_ranks(global_indices);
    }

    VectorField increment_gradient(const Field& phi)
    {
        const int n = fixture.nodes;
        pressure.upload(phi);
        std::array<DeviceArray<double>, 3> accumulator, gradient;
        for (int d = 0; d < 3; ++d) { accumulator[d].zero(n); gradient[d].zero(n); }
        for (int element = 0; element < fixture.elements; ++element) if (element%ranks == rank) {
            applyDivTransposePerNodeKernel<KeyType,double,TetTag><<<1,128>>>(
                connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
                nullptr, nullptr, nullptr, nullptr, pressure.data,
                area[0].data, area[1].data, area[2].data,
                accumulator[0].data, accumulator[1].data, accumulator[2].data, element, 1);
            complete("production pressure transpose scatter");
        }
        for (int d = 0; d < 3; ++d) accumulator[d].upload(sum_ranks(accumulator[d].download()));
        normalizeGradientPerNodeKernel<double><<<1,128>>>(
            accumulator[0].data, accumulator[1].data, accumulator[2].data, mass.data,
            node_to_dof.data, ownership.data, gradient[0].data, gradient[1].data, gradient[2].data, n);
        negateThreeOwnedKernel<double><<<1,128>>>(gradient[0].data, gradient[1].data, gradient[2].data,
                                                 node_to_dof.data, ownership.data, n);
        addOutletGradientTermKernel<double><<<1,128>>>(pressure.data, nullptr,
            outlet_area[0].data, outlet_area[1].data, outlet_area[2].data, mass.data,
            node_to_dof.data, ownership.data, fixed_dof.data, owned_count,
            gradient[0].data, gradient[1].data, gradient[2].data, n);
        complete("production normalized boundary increment gradient");
        VectorField result;
        for (int d = 0; d < 3; ++d) result[d] = sum_ranks(gradient[d].download());
        return result;
    }
};

void run_fixture(int num_elements, int rank, int ranks, Checks& checks)
{
    Fixture fixture(num_elements);
    DeviceFixture device(fixture, rank, ranks, checks);
    const int n = fixture.nodes;
    const auto b = fixture.divergence();
    Field direction(n);
    for (int node = 0; node < n; ++node) direction[node] = .17*(node+1)*(node%2 ? -1 : 1);
    const auto gradient = device.increment_gradient(direction);
    for (int node = 0; node < n; ++node) if (!fixture.fixed[node])
        for (int d = 0; d < 3; ++d) {
            double expected = 0;
            for (int row = 0; row < n; ++row) expected -= b[row][3*node+d]*direction[row]/fixture.mass[node];
            checks.near(gradient[d][node], expected, 1e-12, "actual boundary-aware increment gradient");
        }
    for (bool nodal : {false, true}) {
        double exterior = 0, expected_exterior = 0;
        const auto base = device.residual(fixture.pressure, fixture.velocity, nodal, &exterior);
        const auto expected = fixture.reference_residual(fixture.pressure, fixture.velocity, nodal, &expected_exterior);
        double sum = 0;
        for (int node = 0; node < n; ++node) {
            checks.near(base[node], expected[node], 1e-12, "actual interior/boundary scatter per continuity row");
            sum += base[node];
        }
        checks.near(sum, exterior, 1e-12, "actual continuity versus production boundary reporter");
        checks.near(exterior, expected_exterior, 1e-12, "independent exterior sum");
        for (int startup = 0; startup < 2; ++startup) {
            const double dt = .14, rho = 2.7;
            const double h = (startup == 0 ? dt : 2*dt/3)/rho;
            const auto compact = device.compact_matrix(1/h, nodal);
            for (int column = 0; column < n; ++column) {
                auto perturbed = fixture.pressure;
                constexpr double epsilon = 1e-5;
                perturbed[column] += epsilon;
                const auto changed = device.residual(perturbed, fixture.velocity, nodal);
                for (int row = 0; row < n; ++row)
                    checks.near((changed[row]-base[row])/epsilon, h*compact[row*n+column], 1e-9,
                                "actual owner-row compact CSR derivative versus pressure finite difference");
            }
            for (double epsilon : {1e-3,1e-4,1e-5}) {
                auto perturbed_pressure = fixture.pressure;
                auto perturbed_velocity = fixture.velocity;
                for (int node = 0; node < n; ++node) {
                    perturbed_pressure[node] += epsilon*direction[node];
                    if (!fixture.fixed[node]) for (int d = 0; d < 3; ++d)
                        perturbed_velocity[d][node] -= epsilon*h*gradient[d][node];
                }
                const auto changed = device.residual(perturbed_pressure, perturbed_velocity, nodal);
                for (int row = 0; row < n; ++row) {
                    double action = 0;
                    for (int col = 0; col < n; ++col) action += h*compact[row*n+col]*direction[col];
                    for (int node = 0; node < n; ++node) if (!fixture.fixed[node])
                        for (int d = 0; d < 3; ++d) action -= h*b[row][3*node+d]*gradient[d][node];
                    checks.near((changed[row]-base[row])/epsilon, action, 1e-9,
                                "actual CUDA full pressure/velocity correction JVP");
                }
            }
        }
    }
    int local_empty = device.local_facets.count == 0 ? 1 : 0, empty_ranks = 0;
    MPI_Allreduce(&local_empty, &empty_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    checks.require(empty_ranks == std::max(0, ranks-num_elements), "actual MPI empty-facet participation");
}

void check_anchor_failures(int rank, Checks& checks)
{
    std::array<DeviceArray<KeyType>, 4> connectivity;
    for (unsigned k = 0; k < 4; ++k) connectivity[k].upload({k});
    std::array<DeviceArray<double>, 3> coordinates, area;
    for (int d = 0; d < 3; ++d) {
        Field values(4, 0.);
        values[d+1] = 1.;
        coordinates[d].upload(values);
        area[d].upload({.5});
    }
    DeviceArray<int> nodes, element, opposite, map, rows, columns, invalid;
    DeviceArray<uint8_t> ownership;
    DeviceArray<double> matrix;
    nodes.upload({1,2,3}); rows.upload({0,1});
    // Only node 1 owns a row. Its opposite pressure is a valid ghost column (3).
    const std::vector<int> valid_map{3,0,1,2};
    const std::vector<uint8_t> valid_owner{0,1,0,0};
    auto run = [&](std::vector<int> dofs, std::vector<uint8_t> owners, int column,
                   int elem, int opp, bool have_facet, int expected, const char* label) {
        map.upload(dofs); ownership.upload(owners); columns.upload({column});
        element.upload({elem}); opposite.upload({opp});
        invalid.zero(1); matrix.zero(1);
        // Other ranks have no matrix facets but must receive rank 0's failure status.
        const int facets = rank == 0 && have_facet ? 1 : 0;
        if (facets > 0) {
            assemble_outlet_pressure_sensitivity_kernel<KeyType,double><<<1,128>>>(
                connectivity[0].data, connectivity[1].data, connectivity[2].data, connectivity[3].data,
                nodes.data, element.data, opposite.data, area[0].data, area[1].data, area[2].data,
                coordinates[0].data, coordinates[1].data, coordinates[2].data,
                nullptr, 2., 1., map.data, ownership.data, rows.data, columns.data,
                1, 4, matrix.data, invalid.data, facets);
            complete("outlet anchor fault injection");
        }
        int local_invalid = invalid.download()[0], any_invalid = 0;
        MPI_Allreduce(&local_invalid, &any_invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        checks.require(any_invalid == expected, label);
        const double expected_entry = facets && !expected && owners[1] == 1 ? 1. : 0.;
        checks.near(matrix.download()[0], expected_entry, 1e-14, "anchor write or safe rejection");
    };
    run(valid_map, valid_owner, 3, 0, 0, true, 0, "partial face ownership and ghost column accepted");
    run({3,0,-1,-1}, valid_owner, 3, 0, 0, true, 0, "non-owned face rows need no local DOF");
    run(valid_map, valid_owner, 2, 0, 0, true, 1, "missing CSR column rejected collectively");
    run({-1,0,1,2}, valid_owner, 3, 0, 0, true, 1, "missing opposite DOF rejected collectively");
    run({4,0,1,2}, valid_owner, 3, 0, 0, true, 1, "out-of-range opposite DOF rejected collectively");
    run({3,-1,1,2}, valid_owner, 3, 0, 0, true, 1, "missing owned row rejected collectively");
    run({3,1,1,2}, valid_owner, 3, 0, 0, true, 1, "owned node mapped to ghost row rejected");
    run(valid_map, valid_owner, 3, -1, 0, true, 1, "unresolved facet rejected collectively");
    run(valid_map, valid_owner, 3, 0, -1, true, 1, "unresolved opposite rejected collectively");
    run({-1,-1,-1,-1}, {0,0,0,0}, 3, 0, 0, true, 0, "no owned rows is a valid skip");
    run(valid_map, valid_owner, 3, 0, 0, false, 0, "empty facet ranks participate without error");
}

void check_krylov_products(int rank, int ranks, Checks& checks)
{
    constexpr int count = 4, n = 513, stride = 520;
    Field basis(count*stride, 1e90), vector(n), reference(count, 0.);
    for (int i = 0; i < n; ++i)
    {
        vector[i] = std::sin(.03*i);
        for (int j = 0; j < count; ++j)
        {
            basis[j*stride+i] = std::cos(.02*i+j);
            reference[j] += basis[j*stride+i]*vector[i];
        }
    }
    DeviceArray<double> d_basis, d_vector, d_products;
    d_basis.upload(basis); d_vector.upload(vector); d_products.zero(count);
    for (bool empty_peers : {false, true})
    {
        const int local_n = empty_peers && rank != 0 ? 0 : n;
        outlet_krylov_products_kernel<double><<<count,256>>>(
            d_basis.data, d_vector.data, stride, local_n, d_products.data);
        complete("batched Krylov projection");
        const auto result = sum_ranks(d_products.download());
        for (int j = 0; j < count; ++j)
            checks.near(result[j], reference[j]*(empty_peers ? 1 : ranks), 1e-10,
                        "batched Krylov projection, padding and empty-rank reduction");
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (ranks != 1 && ranks != 2 && ranks != 4) {
        if (rank == 0) std::fprintf(stderr, "Use one, two, or four MPI ranks.\n");
        MPI_Finalize();
        return 2;
    }
    // The launcher selects each process's visible GPU, as with the scalar gate.
    cuda_check(cudaFree(nullptr), "initialize visible GPU");
    Checks checks;
    check_anchor_failures(rank, checks);
    check_krylov_products(rank, ranks, checks);
    run_fixture(1, rank, ranks, checks);
    run_fixture(2, rank, ranks, checks);
    int failures = 0;
    MPI_Allreduce(&checks.failures, &failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        std::printf("%s: %d checks per rank, actual CUDA scatters and owner-row derivatives, %d MPI ranks\n",
                    failures ? "FAIL" : "PASS", checks.count, ranks);
        std::printf("Replicated synthetic fixtures use real MPI reductions; ElementDomain halo and flow validation are not covered.\n");
    }
    MPI_Finalize();
    return failures ? 1 : 0;
}
