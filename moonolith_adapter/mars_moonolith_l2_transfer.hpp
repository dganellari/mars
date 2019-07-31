#ifndef MOONOLITH_MARS_L2_TRANSFER_HPP
#define MOONOLITH_MARS_L2_TRANSFER_HPP

#include "par_moonolith_config.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_duplicate_intersection_avoidance.hpp"
#include "moonolith_one_master_one_slave_algorithm.hpp"
#include "moonolith_single_collection_one_master_one_slave_algorithm.hpp"

#include "moonolith_l2_assembler.hpp"
#include "moonolith_keast_quadrature_rule.hpp"
#include "moonolith_gauss_quadrature_rule.hpp"
#include "mars_moonolith_function_space_adapter.hpp"
#include "mars_moonolith_mesh_adapter.hpp"

#include <vector>
#include <tuple>

namespace moonolith {

    template<class MarsElem>
    class MarsMoonolithConverter {
    public:
        //using Type = ...;
    };

    template<Integer Dim>
    class MarsMoonolithConverter<mars::Simplex<Dim, 1>> {
    public:
        using ElemType = moonolith::Edge2<double, Dim>;
        using QuadratureType = moonolith::Quadrature1<double>;
    };

    template<Integer Dim>
    class MarsMoonolithConverter<mars::Simplex<Dim, 2>> {
    public:
        using ElemType = moonolith::Tri3<double, Dim>;
        using QuadratureType = moonolith::Quadrature2<double>;
    };

    template<Integer Dim>
    class MarsMoonolithConverter<mars::Simplex<Dim, 3>> {
    public:
        using ElemType = moonolith::Tet4<double, Dim>;
        using QuadratureType = moonolith::Quadrature3<double>;
    };

    ////////////////////////////////////////////////////////////
    ///@brief assembles tranfer operator T = Q D^{-1} B
    template<Integer Dim, Integer ManifoldDim>
    class MarsMeshTransfer {
    public:
        using MeshT      = mars::Mesh<Dim, ManifoldDim>;
        // using MeshT      = mars::IMesh<Dim>;
        using SpaceT     = moonolith::FunctionSpace<MeshT>;

        using MarsElemT             = typename MeshT::Elem; 
        using MeshAlgorithmT        = moonolith::OneMasterOneSlaveAlgorithm<Dim, MeshT>;
        using SpaceAlgorithmT       = moonolith::OneMasterOneSlaveAlgorithm<Dim, SpaceT>;
        using SingleSpaceAlgorithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, SpaceT>;

        using MeshElemAdapter         = typename MeshAlgorithmT::Adapter;
        using SpaceElemAdapter        = typename SpaceAlgorithmT::Adapter;
        using SingleSpaceElemAdapter  = typename SingleSpaceAlgorithmT::Adapter;

        using Converter = MarsMoonolithConverter<MarsElemT>;

        using Elem = typename Converter::ElemType;
        using QuadratureT = typename Converter::QuadratureType;
        using L2TransferT = moonolith::L2Transfer<Elem, Elem>;
        
        class MeshAssembler {
        public:

            inline bool operator()(const MeshElemAdapter &master, const MeshElemAdapter &slave)
            {
                auto &e_m = master.elem();
                auto &e_s = slave.elem();

                auto &m_m = master.collection();
                auto &m_s = slave.collection();

                make(m_m, e_m, master_elem);
                make(m_s, e_s, slave_elem);

                if(algo.assemble(master_elem, slave_elem)) {
                    uint n_nodes_master = mars::n_nodes(e_m);
                    uint n_nodes_slave  = mars::n_nodes(e_s);
                    
                    const auto &B_e = algo.coupling_matrix();
                    const auto &D_e = algo.mass_matrix();
                    const auto &Q_e = algo.transformation();



                    for(uint i = 0; i < n_nodes_slave; ++i) {
                        for(uint j = 0; j < n_nodes_master; ++j) {
                            B.add(e_s.nodes[i], e_m.nodes[j], B_e(i, j));
                        }

                        for(uint j = 0; j < n_nodes_slave; ++j) {
                            D.add(e_s.nodes[i], e_s.nodes[j], D_e(i, j));

                            auto Q_val = Q_e(i, j);

                            if(std::abs(Q_val) > 1e-15) {
                                Q.set(e_s.nodes[i], e_s.nodes[j], Q_val);
                            }
                        }
                    }

                    return true;
                } else {
                    return false;
                }
            }

            MeshAssembler(Communicator &comm)
            : B(comm), D(comm), Q(comm)
            {
                Gauss::get(Elem::Order + Elem::Order, q_rule);
                algo.set_quadrature(q_rule);
            }

            QuadratureT q_rule;
            L2TransferT algo;
            Elem master_elem, slave_elem;
            SparseMatrix<double> B, D, Q;
        };

        class SpaceAssembler {
        public:

            inline bool operator()(const SpaceElemAdapter &master, const SpaceElemAdapter &slave)
            {
                auto &e_m = master.elem();
                auto &e_s = slave.elem();

                auto &m_m = master.collection();
                auto &m_s = slave.collection();

                auto &dof_m = m_m.dof_map().dofs(master.element_handle());
                auto &dof_s = m_s.dof_map().dofs(slave.element_handle());

                make(m_m, e_m, master_elem);
                make(m_s, e_s, slave_elem);

                if(algo.assemble(master_elem, slave_elem)) {
                    uint n_nodes_master = mars::n_nodes(e_m);
                    uint n_nodes_slave  = mars::n_nodes(e_s);
                    
                    const auto &B_e = algo.coupling_matrix();
                    const auto &D_e = algo.mass_matrix();
                    const auto &Q_e = algo.transformation();

                    // D_e.describe(std::cout);

                    for(uint i = 0; i < n_nodes_slave; ++i) {
                        for(uint j = 0; j < n_nodes_master; ++j) {
                            B.add(dof_s[i], dof_m[j], B_e(i, j));
                        }

                        for(uint j = 0; j < n_nodes_slave; ++j) {
                            D.add(dof_s[i], dof_s[j], D_e(i, j));

                            auto Q_val = Q_e(i, j);

                            if(std::abs(Q_val) > 1e-15) {
                                Q.set(dof_s[i], dof_s[j], Q_val);
                            }
                        }
                    }

                    return true;
                } else {
                    return false;
                }
            }

            SpaceAssembler(Communicator &comm)
            : B(comm), D(comm), Q(comm)
            {
                //FIXME
                Gauss::get(Elem::Order + Elem::Order, q_rule);
                algo.set_quadrature(q_rule);
            }

            QuadratureT q_rule;
            L2TransferT algo;
            Elem master_elem, slave_elem;

            SparseMatrix<double> B, D, Q;
        };

        MarsMeshTransfer(Communicator &comm) : comm(comm) {}

        bool assemble(
            const MeshT &mesh_master,
            const MeshT &mesh_slave)
        {
            double elapsed = MPI_Wtime();

            MeshAlgorithmT algo(comm, moonolith::make_unique<CollectionManager<MeshT>>());
            
            algo.init_simple(
                mesh_master,
                mesh_slave,
                0.0
            );

            elapsed = MPI_Wtime() - elapsed;
            logger() << "init: " << elapsed  << std::endl;

            ////////////////////////////////////////////////////
            /////////////////// pair-wise method ///////////////////////
            elapsed = MPI_Wtime();
            
            MeshAssembler assembler(comm);
            algo.compute([&](const MeshElemAdapter &master, const MeshElemAdapter &slave) -> bool {
                if(assembler(master, slave)) {
                    return true;
                }

                return false;
            });

            double vol = assembler.algo.intersection_measure();

            comm.all_reduce(&vol, 1, MPISum());

            double sum_mat_B = assembler.B.sum();
            double sum_mat_D = assembler.D.sum();

            elapsed = MPI_Wtime() - elapsed;
            logger() << "time MarsMeshTransfer::assemble: " << elapsed  << std::endl;
            logger() << "vol: " << vol << " sum(B): " << sum_mat_B << " sum(D): " << sum_mat_D << std::endl;
            return vol > 0.0;
        }

        bool assemble(
            const SpaceT &mesh_master,
            const SpaceT &mesh_slave)
        {
            double elapsed = MPI_Wtime();

            SpaceAlgorithmT algo(comm, moonolith::make_unique<CollectionManager<SpaceT>>());
            
            algo.init_simple(
                mesh_master,
                mesh_slave,
                0.0
            );

            elapsed = MPI_Wtime() - elapsed;
            logger() << "init: " << elapsed  << std::endl;

            ////////////////////////////////////////////////////
            /////////////////// pair-wise method ///////////////////////
            elapsed = MPI_Wtime();
            
            SpaceAssembler assembler(comm);
            algo.compute([&](const SpaceElemAdapter &master, const SpaceElemAdapter &slave) -> bool {
                if(assembler(master, slave)) {
                    return true;
                }

                return false;
            });

            double vol = assembler.algo.intersection_measure();

            comm.all_reduce(&vol, 1, MPISum());

            double sum_mat_B = assembler.B.sum();
            double sum_mat_D = assembler.D.sum();

            elapsed = MPI_Wtime() - elapsed;
            logger() << "time MarsMeshTransfer::assemble: " << elapsed  << std::endl;
            logger() << "vol: " << vol << " sum(B): " << sum_mat_B << " sum(D): " << sum_mat_D << std::endl;
            return vol > 0.0;
        }


        bool assemble(
            const std::vector<std::pair<int, int>> &tags,
            const SpaceT &space
            )
        {
            double elapsed = MPI_Wtime();

            SingleSpaceAlgorithmT algo(comm, moonolith::make_unique<CollectionManager<SpaceT>>());
            
            algo.init(
                space,
                tags,
                0.0
            );

            elapsed = MPI_Wtime() - elapsed;
            logger() << "init: " << elapsed  << std::endl;

            ////////////////////////////////////////////////////
            /////////////////// pair-wise method ///////////////////////
            elapsed = MPI_Wtime();
            
            SpaceAssembler assembler(comm);
            algo.compute([&](const SpaceElemAdapter &master, const SpaceElemAdapter &slave) -> bool {
                if(assembler(master, slave)) {
                    return true;
                }

                return false;
            });

            double vol = assembler.algo.intersection_measure();

            comm.all_reduce(&vol, 1, MPISum());

            double sum_mat_B = assembler.B.sum();
            double sum_mat_D = assembler.D.sum();

            elapsed = MPI_Wtime() - elapsed;
            logger() << "time MarsMeshTransfer::assemble: " << elapsed  << std::endl;
            logger() << "vol: " << vol << " sum(B): " << sum_mat_B << " sum(D): " << sum_mat_D << std::endl;
            return vol > 0.0;
        }

        Communicator &comm;
    };

}

#endif //MOONOLITH_MARS_L2_TRANSFER_HPP
