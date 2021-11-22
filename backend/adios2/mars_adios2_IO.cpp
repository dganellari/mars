
#include "mars_adios2_IO.hpp"

#include "mars_distributed_base_data_management.hpp"
#include "mars_distributed_dof_management.hpp"
#include "mars_distributed_finite_element.hpp"

#include "mars_globals.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_mesh_kokkos.hpp"
#endif

////////////////////////////////////

#include <cmath>
#include <complex>
#include <iostream>

#include "adios2.h"

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

////////////////////////////////////

namespace mars {
    namespace adios2 {

        template <class Mesh>
        class Adios2Helper {
        public:
        };

        template <>
        class Adios2Helper<ParallelMesh3> {
        public:
            static Integer n_nodes(const ParallelMesh3& mesh) { return mesh.n_nodes(); }
            static Integer n_elements(const ParallelMesh3& mesh) { return mesh.n_elements(); }
            static Integer n_nodes_x_element(const ParallelMesh3& mesh) {
                auto elementsmesh = mesh.get_view_elements();
                size_t element_nvertices = static_cast<size_t>(elementsmesh.extent(1));
                return element_nvertices;
            }

            static void write_elements(const ParallelMesh3& mesh,
                                       ::adios2::Variable<uint64_t>::Span& span_connectivity,
                                       ::adios2::Variable<int32_t>::Span& span_element_attribute) {
                size_t elementPosition = 0;
                Integer ne = n_elements(mesh);
                auto elementsmesh = mesh.get_view_elements();

                for (Integer e = 0; e < ne; ++e) {
                    const int meshNNodes = mars::ParallelMesh3::Elem::NNodes;
                    span_element_attribute[e] = static_cast<int32_t>(e);
                    span_connectivity[elementPosition] = meshNNodes;
                    for (int v = 0; v < meshNNodes; ++v) {
                        span_connectivity[elementPosition + v + 1] = elementsmesh(e, v);
                    }
                    elementPosition += meshNNodes + 1;
                }
            }

            static void write_nodes(const ParallelMesh3& mesh, ::adios2::Variable<Real>::Span& span_vertices) {
                auto points = mesh.get_view_points();
                Integer nn = n_nodes(mesh);
                for (Integer v = 0; v < nn; ++v) {
                    const int space_dim = static_cast<int>(mesh.Dim);
                    for (int coord = 0; coord < space_dim; ++coord) {
                        span_vertices[v * space_dim + coord] = points(v, coord);
                    }
                }
            }

            static void write_field(const ParallelMesh3& mesh,
                                    int n_components,
                                    const ViewVectorType<Real>& data,
                                    ::adios2::Variable<Real>::Span& span) {
                Integer nn = n_nodes(mesh);

                for (int v = 0; v < nn; ++v) {
                    for (int d = 0; d < n_components; ++d) {
                        span[v * n_components + d] = data(v * n_components + d);
                    }
                }
            }
        };

        template <class DofHandler>
        class Adios2Helper<FEDofMap<DofHandler>> {
        public:
            using FEDofMap = mars::FEDofMap<DofHandler>;

            static Integer n_nodes(const FEDofMap& fe_dof_map) {
                auto dof = fe_dof_map.get_dof_handler().get_local_dof_enum();
                return dof.get_elem_size();
            }

            static Integer n_elements(const FEDofMap& fe_dof_map) {
                return fe_dof_map.get_dof_handler().get_elem_size();
            }

            static Integer n_nodes_x_element(const FEDofMap& fe_dof_map) {
                return fe_dof_map.get_dof_handler().get_elem_type();
            }

            static void write_elements(const FEDofMap& fe_dof_map,
                                       ::adios2::Variable<uint64_t>::Span& span_connectivity,
                                       ::adios2::Variable<int32_t>::Span& span_element_attribute) {
                auto dof_handler = fe_dof_map.get_dof_handler();
                int block_size = dof_handler.get_block();
                int nne = n_nodes_x_element(dof_handler);
                int ne = n_elements(dof_handler);

                ViewVectorType<int> cells("cells", ne * (nne + 1));
                ViewVectorType<int>::HostMirror cells_host = Kokkos::create_mirror_view(cells);

                dof_handler.elem_iterate([&](const Integer elem_index) {
                    auto offset = elem_index * (nne + 1);
                    cells(offset) = nne;

                    for (int i = 0; i < nne; i++) {
                        const Integer local_dof = fe_dof_map.get_elem_local_dof(elem_index, i * block_size);
                        const Integer base = dof_handler.compute_base(local_dof);
                        assert(i < 8);

                        cells(offset + 1 + i) = base;
                    }
                });

                Kokkos::deep_copy(cells_host, cells);

                for (Integer i = 0; i < ne; ++i) {
                    span_connectivity[i] = i;
                }

                auto len = cells_host.size();
                for (Integer i = 0; i < len; ++i) {
                    span_connectivity[i] = cells_host[i];
                }
            }

            static void write_nodes(const FEDofMap& fe_dof_map, ::adios2::Variable<Real>::Span& span_vertices) {
                auto dof_handler = fe_dof_map.get_dof_handler();
                auto dof = dof_handler.get_local_dof_enum();

                Integer n_nodes = dof.get_elem_size();
                Integer block_size = dof_handler.get_block();
                int dim = DofHandler::Dim;

                ViewVectorType<Real> points("points", n_nodes * dim);
                ViewVectorType<Real>::HostMirror points_host = Kokkos::create_mirror_view(points);

                Kokkos::parallel_for("for", n_nodes, [&](const int i) {
                    const Integer sfc_elem = dof_handler.local_to_sfc(i * block_size);
                    const Integer global_dof = dof_handler.local_to_global(i);

                    Real point[3] = {0., 0., 0.};
                    get_vertex_coordinates_from_sfc<DofHandler::ElemType>(
                        sfc_elem, point, dof.get_XDim(), dof.get_YDim(), dof.get_ZDim());

                    for (int d = 0; d < dim; ++d) {
                        points[i * dim + d] = point[d];
                    }
                });

                Kokkos::deep_copy(points_host, points);

                for (Integer i = 0; i < n_nodes; ++i) {
                    for (int d = 0; d < dim; ++d) {
                        span_vertices[i * dim + d] = points_host[i * dim + d];
                    }
                }
            }

            static void write_field(const FEDofMap& fe_dof_map,
                                    int n_components,
                                    const ViewVectorType<Real>& owned_data,
                                    ::adios2::Variable<Real>::Span& span) {
                auto dof_handler = fe_dof_map.get_dof_handler();
                auto& ctx = dof_handler.get_context();
                int comm_size = ::mars::num_ranks(ctx);

                assert(dof_handler.get_owned_dof_size() == owned_data.extent(0));
                const Integer size = dof_handler.get_owned_dof_size();

                ViewVectorType<Integer> global_to_sfc = dof_handler.get_global_dof_enum().get_view_elements();

                ViewVectorType<Real> local_data("local_data", dof_handler.get_dof_size());

                Kokkos::parallel_for(
                    "set_locally_owned_data", size, MARS_LAMBDA(const Integer i) {
                        const Integer base = dof_handler.compute_base(i);
                        const Integer component = dof_handler.compute_component(i);
                        const Integer sfc = global_to_sfc(base);

                        const Integer local = dof_handler.sfc_to_local(sfc);
                        local_data(dof_handler.compute_block_index(local, component)) = owned_data(i);
                    });

                if (comm_size > 1) {
                    ::mars::gather_ghost_data(dof_handler, local_data);
                }

                typename ViewVectorType<Real>::HostMirror data_host = Kokkos::create_mirror_view(local_data);
                Kokkos::deep_copy(data_host, local_data);
                Integer n = data_host.extent(0);

                // Integer nn = n_nodes(mesh);

                for (int v = 0; v < n; ++v) {
                    span[v] = data_host[v];
                }
            }
        };  // namespace adios2

        template <class Mesh>
        class IO<Mesh>::Impl {
        public:
            using Adios2Helper = mars::adios2::Adios2Helper<Mesh>;

            static constexpr int Dim = Mesh::Dim;

            void header(std::ostream& vtk_schema) const {
                vtk_schema << R"(
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.2" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="NumOfVertices" NumberOfCells="NumOfElements">
      <Points>
        <DataArray Name="vertices" />)";

                vtk_schema << R"(
      </Points>
          <CellData>
        <DataArray Name="attribute" />
    </CellData>
      <Cells>
        <DataArray Name="connectivity" />
        <DataArray Name="types" />
      </Cells>
      <PointData>)";

                if (point_fields.empty()) {
                    vtk_schema << "\n";
                } else {
                    for (const auto& point_datum : point_fields) {
                        vtk_schema << "        <DataArray Name=\"" << point_datum->name << "\"/>\n";
                    }
                }

                vtk_schema << "        <DataArray Name=\"TIME\">\n";
                vtk_schema << "          TIME\n";
                vtk_schema << "        </DataArray>\n";

                vtk_schema << R"(
          </PointData>
          </Piece>
        </UnstructuredGrid>
      </VTKFile>)";
            }

            bool open_output() {
                engine = io.Open(output_path, ::adios2::Mode::Write);

                // FIXME: handle failures
                return true;
            }

            bool write_step() {
                engine.BeginStep();

                // TODO
                for (const auto& point_datum : point_fields) {
                    point_datum->set(mesh, engine, io);
                }

                engine.EndStep();

                // FIXME: handle failures
                return true;
            }

            void close() { engine.Close(); }

            void define_variables() {
                for (const auto& point_datum : point_fields) {
                    point_datum->define(io);
                }
            }

            void define_mesh() {
                uint32_t n_nodes = Adios2Helper::n_nodes(mesh);
                uint64_t n_elements = Adios2Helper::n_elements(mesh);
                uint64_t n_nodes_x_element = Adios2Helper::n_nodes_x_element(mesh);

                ///////////////////////////

                io.DefineVariable<uint32_t>("NumOfElements", {::adios2::LocalValueDim});
                io.DefineVariable<uint32_t>("NumOfVertices", {::adios2::LocalValueDim});
                io.DefineVariable<Real>("vertices", {}, {}, {n_nodes, Dim});
                io.DefineVariable<uint64_t>("connectivity", {}, {}, {n_elements, n_nodes_x_element + 1});
                io.DefineVariable<int32_t>("attribute", {}, {}, {n_elements});

                ///////////////////////////
                io.DefineVariable<uint32_t>("types");

                // VTK type of mesh:
                // 1: POINT, 3:SEGMENT, 5:TRIANGLE, 9:SQUARE, 10:TETRAHEDRON, 12:CUBE, 13:PRISM
                uint32_t vtktype = 12;
                if (Dim == 2) {
                    vtktype = 9;
                }

                ::adios2::Variable<uint32_t> var_types = io.InquireVariable<uint32_t>("types");
                engine.Put(var_types, vtktype);
                ///////////////////////////////////////////////////////////////////////////////

                // More info written in .bp file.
                std::string mesh_type = "VTU";
                std::vector<std::string> viz_tools;
                viz_tools.push_back("Paraview: ADIOS2VTXReader");
                viz_tools.push_back("VTK: vtkADIOS2VTXReader.h");

                io.DefineAttribute<std::string>("format/mars_mesh", mesh_type);
                io.DefineAttribute<std::string>("format/viz_tools", viz_tools.data(), viz_tools.size());

                ///////////////////////////////////////////////////////////////////////////////
                std::stringstream ss;
                header(ss);

                io.DefineAttribute<std::string>("vtk.xml", ss.str(), {}, {});
                ///////////////////////////////////////////////////////////////////////////////

                engine.Put("NumOfElements", static_cast<uint32_t>(n_elements));

                ::adios2::Variable<int32_t> var_element_attribute = io.InquireVariable<int32_t>("attribute");
                ::adios2::Variable<int32_t>::Span span_element_attribute = engine.Put<int32_t>(var_element_attribute);

                ::adios2::Variable<uint64_t> var_connectivity = io.InquireVariable<uint64_t>("connectivity");
                ::adios2::Variable<uint64_t>::Span span_connectivity = engine.Put<uint64_t>(var_connectivity);

                Adios2Helper::write_elements(mesh, span_connectivity, span_element_attribute);

                engine.Put("NumOfVertices", static_cast<uint32_t>(n_nodes));
                ::adios2::Variable<Real> var_vertices = io.InquireVariable<Real>("vertices");
                ::adios2::Variable<Real>::Span span_vertices = engine.Put<Real>(var_vertices);

                Adios2Helper::write_nodes(mesh, span_vertices);
            }

            Impl(Mesh& mesh, MPI_Comm comm) : mesh(mesh), adios(comm), io(adios.DeclareIO("BPFile_SZ")), engine() {}

            class CompareFields {
            public:
                bool operator()(const std::unique_ptr<Field>& l, const std::unique_ptr<Field>& r) const {
                    return l->name < r->name;
                }
            };

            Mesh& mesh;
            ::adios2::ADIOS adios;
            ::adios2::IO io;
            ::adios2::Engine engine;
            std::set<std::unique_ptr<Field>, CompareFields> point_fields;
            std::string output_path{"out.bp"};
        };

        template <class Mesh>
        void IO<Mesh>::RealField::define(::adios2::IO& io) {
            if (this->n_components == 1) {
                io.DefineVariable<Real>(this->name, {}, {}, {this->data.size()});
            } else {
                ::adios2::Dims count =
                    ::adios2::Dims{this->data.size() / this->n_components, size_t(this->n_components)};
                io.DefineVariable<Real>(this->name, {}, {}, count);
            }
        }

        template <class Mesh>
        void IO<Mesh>::RealField::set(const Mesh& mesh, ::adios2::Engine& engine, ::adios2::IO& io) {
            using Adios2Helper = mars::adios2::Adios2Helper<Mesh>;

            ::adios2::Variable<Real> var = io.InquireVariable<Real>(this->name);
            ::adios2::Variable<Real>::Span span = engine.Put<Real>(var);
            Adios2Helper::write_field(mesh, this->n_components, this->data, span);
        }

        template <class Mesh>
        IO<Mesh>::~IO() = default;

        template <class Mesh>
        IO<Mesh>::IO(Mesh& mesh) : impl_(std::make_unique<Impl>(mesh, MPI_COMM_WORLD)) {}

        template <class Mesh>
        void IO<Mesh>::add_field(const std::string& name, const int n_components, const ViewVectorType<Real>& data) {
            impl_->point_fields.insert(std::make_unique<RealField>(name, n_components, data));
        }

        template <class Mesh>
        void IO<Mesh>::set_output_path(const std::string& path) {
            impl_->output_path = path;
        }

        template <class Mesh>
        bool IO<Mesh>::open_output() {
            return impl_->open_output();
        }

        template <class Mesh>
        bool IO<Mesh>::write_step() {
            return impl_->write_step();
        }

        template <class Mesh>
        void IO<Mesh>::close() {
            impl_->close();
        }

        template <class Mesh>
        bool IO<Mesh>::aux_write(const std::string& path, ViewVectorType<Real>& data) {
            // TODO
            set_output_path(path);
            add_field("U", -1, data);
            if (!open_output()) {
                return false;
            }

            write_step();

            close();
            return true;
        }

        template class IO<mars::ParallelMesh3>;

        using DMesh2 = ::mars::DistributedMesh<::mars::ElementType::Quad4>;
        using DMesh3 = ::mars::DistributedMesh<::mars::ElementType::Hex8>;

        template class IO<mars::FEDofMap<mars::DofHandler<DMesh2, 1, 0>>>;
        template class IO<mars::FEDofMap<mars::DofHandler<DMesh3, 1, 0>>>;

    }  // namespace adios2
}  // namespace mars
