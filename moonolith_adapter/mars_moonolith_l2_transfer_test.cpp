#include "mars_moonolith_l2_transfer_test.hpp"
#include "mars_moonolith_l2_transfer.hpp"
#include "generation/mars_mesh_generation.hpp"
#include "mars_par_mesh.hpp"
#include "mars_moonolith_function_space_adapter.hpp"

#include <cassert>

using namespace moonolith;

namespace mars {

    static void run_mars_l2_transfer_2()
    {
        logger() << "--------------------------------\n";
        logger() << "run_mars_l2_transfer_2" << std::endl;

        moonolith::Communicator comm;

        Integer res = 2;
        Integer n_ref_master = 6;
        Integer n_ref_slave  = 6;

        mars::Mesh2 mesh_master;

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_master);
            mars::Bisection<mars::Mesh2> b(mesh_master);
            b.uniform_refine(res * n_ref_master);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);
        }

        mars::Mesh2 mesh_slave;

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_slave);
            mars::Bisection<mars::Mesh2> b(mesh_slave);
            b.uniform_refine(res * n_ref_slave);
            mesh_slave.clean_up();

            mars::write_mesh("mesh_s", mesh_slave);
        }

        comm.barrier();

        logger() << "n_master " << mesh_master.n_elements() << " n_slave " << mesh_slave.n_elements() << " " << std::endl;

        MarsMeshTransfer<2, 2> transfer(comm);
        if(!transfer.assemble(mesh_master, mesh_slave)) {
            assert(false);
        }

        comm.barrier();
    }


    static const int N_MASTER = 3;
    static const int N_SLAVE  = 4;

    static void run_mars_l2_transfer_3()
    {
        logger() << "--------------------------------\n";
        logger() << "run_mars_l2_transfer_3" << std::endl;

        moonolith::Communicator comm;

       Integer n_master = N_MASTER;
       Integer n_slave  = N_SLAVE;
        
        mars::Mesh3 mesh_master;

        if(comm.rank() == 0 || comm.is_alone()) {
            mars::generation::generate_cube(mesh_master, n_master, n_master, n_master);
            mars::write_mesh("mesh_m", mesh_master);
        }

        mars::Mesh3 mesh_slave;

        if(comm.rank() == 1 || comm.is_alone()) {
            mars::generation::generate_cube(mesh_slave, n_slave, n_slave, n_slave);
            mars::write_mesh("mesh_s", mesh_slave);
        }

        comm.barrier();

        logger() << "n_master " << mesh_master.n_elements() << " n_slave " << mesh_slave.n_elements() << " " << std::endl;

        MarsMeshTransfer<3, 3> transfer(comm);
        if(!transfer.assemble(mesh_master, mesh_slave)) {
            assert(false);
        }

        comm.barrier();
    }

    using ParMesh3 = mars::ParMesh<3, 3>;

    void partition_mesh(mars::Mesh3 &mesh, ParMesh3 &parallel_mesh)
    {
        //Bad partioning
        Integer n_parts = parallel_mesh.comm().size();

        std::vector<Integer> partitioning(mesh.n_elements(), 0);

        Integer n_local_elements = mesh.n_elements()/n_parts;

        for(Integer i = 0; i < n_parts; ++i) {
            Integer begin = i * n_local_elements;
            Integer end   = begin + n_local_elements;

           if(i == n_parts - 1) {
                end = mesh.n_elements();
            }

            for(Integer k = begin; k < end; ++k) {
                partitioning[k] = i;
            }
        }

        parallel_mesh.init(mesh, partitioning);
    }

    void concatenate_mesh(mars::Mesh3 &master_mesh, mars::Mesh3 &slave_mesh, ParMesh3 &parallel_mesh)
    {
        for(Integer i = 0; i < master_mesh.n_elements(); ++i) {
            master_mesh.elem(i).block = 1;
        }

        for(Integer i = 0; i < slave_mesh.n_elements(); ++i) {
            slave_mesh.elem(i).block = 2;
        }

        mars::Mesh3 mesh = master_mesh;
        mesh += slave_mesh;

        partition_mesh(mesh, parallel_mesh);
    }


    void concatenate_mesh(mars::Mesh3 &mesh1, mars::Mesh3 &mesh2, mars::Mesh3 &mesh3, ParMesh3 &parallel_mesh)
    {
        for(Integer i = 0; i < mesh1.n_elements(); ++i) {
            mesh1.elem(i).block = 1;
        }

        for(Integer i = 0; i < mesh2.n_elements(); ++i) {
            mesh2.elem(i).block = 2;
        }

        for(Integer i = 0; i < mesh3.n_elements(); ++i) {
            mesh3.elem(i).block = 3;
        }

        mars::Mesh3 mesh = mesh1;
        mesh += mesh2;
        mesh += mesh3;

        partition_mesh(mesh, parallel_mesh);
    }

    void run_mars_l2_transfer_function_space_test()
    {
        logger() << "--------------------------------\n";
        logger() << "run_mars_l2_transfer_function_space_test" << std::endl;

        moonolith::Communicator comm;

        Integer n_master = N_MASTER;
        Integer n_slave  = N_SLAVE;
        
        ParMesh3 mesh_master(comm.get());

        mars::Mesh3 mesh;
        mars::generation::generate_cube(mesh, n_master, n_master, n_master);
        partition_mesh(mesh, mesh_master);
        
        mars::write_mesh(
            "mesh_m" + std::to_string(comm.rank()),
            mesh_master.get_serial_mesh()
        );

        ParMesh3 mesh_slave(comm.get());

        mesh.clear();
        mars::generation::generate_cube(mesh, n_slave, n_slave, n_slave);
        partition_mesh(mesh, mesh_slave);
        
        mars::write_mesh(
            "mesh_s" + std::to_string(comm.rank()),
            mesh_slave.get_serial_mesh()
        );

        ////////////////////////////////////////////////////////////////

        auto space_master = make_function_space(mesh_master);
        auto space_slave  = make_function_space(mesh_slave);

        comm.barrier();

        logger() << "n_master " << mesh_master.n_local_elements() << " n_slave " << mesh_slave.n_local_elements() << " " << std::endl;

        MarsMeshTransfer<3, 3> transfer(comm);
        if(!transfer.assemble(*space_master, *space_slave)) {
            assert(false);
        }

        comm.barrier();
    }


    void run_mars_l2_transfer_function_space_test_with_concatenate_mesh()
    {
        logger() << "--------------------------------\n";
        logger() << "run_mars_l2_transfer_function_space_test_with_concatenate_mesh" << std::endl;

        moonolith::Communicator comm;

        // Integer n1 = N_MASTER;
        // Integer n2 = N_SLAVE;
        // Integer n3 = n1 + n2;


        Integer n1 = 2;
        Integer n2 = 2;
        Integer n3 = 2;
        
        ParMesh3 mesh(comm.get());

        mars::Mesh3 mesh1;
        mars::generation::generate_cube(mesh1, n1, n1, n1);
        
        mars::Mesh3 mesh2(comm.get());
        mars::generation::generate_cube(mesh2, n2, n2, n2);

        mars::Mesh3 mesh3(comm.get());
        mars::generation::generate_cube(mesh3, n3, n3, n3);
        
        concatenate_mesh(mesh1, mesh2, mesh3, mesh);

        ////////////////////////////////////////////////////////////////

        auto space = make_function_space(mesh);

        comm.barrier();

        logger() << "n_elements " << mesh.n_local_elements() << std::endl;

        MarsMeshTransfer<3, 3> transfer(comm);
        if(!transfer.assemble({{1, 2}, {1, 3}}, *space)) {
            assert(false);
        }

        comm.barrier();
    }

    void run_mars_l2_transfer_test()
    {
        run_mars_l2_transfer_2();
        run_mars_l2_transfer_3();
        run_mars_l2_transfer_function_space_test();
        run_mars_l2_transfer_function_space_test_with_concatenate_mesh();
    }
}
