#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include "mars.hpp"
#include "mars_quad4.hpp"

void mesh_test(int &argc, char **&argv, const int level) {
    using namespace mars;
    try {
        mars::proc_allocation resources;
        /*
        // try to detect how many threads can be run on this system
        resources.num_threads = marsenv::thread_concurrency();

        // override thread count if the user set MARS_NUM_THREADS
        if (auto nt = marsenv::get_env_num_threads())
        {
            resources.num_threads = nt;
        } */

#ifdef WITH_MPI
        // initialize MPI
        marsenv::mpi_guard guard(argc, argv, false);

        // assign a unique gpu to this rank if available
        /*  resources.gpu_id = marsenv::find_private_gpu(MPI_COMM_WORLD); */

        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();

        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef WITH_KOKKOS_KERNELS

        Kokkos::Timer timer;
        // create the quad mesh distributed through the mpi procs.
        DistributedQuad4Mesh mesh(context);
        generate_distributed_cube(mesh, level, level, 0);
        /* mesh.set_periodic(); //set the domain to be periodic */
        /*
                const Integer xDim = mesh.get_XDim();
                const Integer yDim = mesh.get_YDim();
                const Integer zDim = mesh.get_ZDim();
         */
        constexpr Integer Dim = DistributedQuad4Mesh::Dim;

        using DOFHandler = DofHandler<DistributedQuad4Mesh, 1>;
        using Elem = typename DistributedQuad4Mesh::Elem;
        // the type of the mesh elements. In this case quad4 (Type=4)
        constexpr Integer Type = Elem::ElemType;

        DOFHandler dof_handler(&mesh);
        dof_handler.enumerate_dofs();
        // create the dm object
        // DMQ2 dm(dof_handler);
        // enumerate the dofs locally and globally. The ghost dofs structures
        // are now created and ready to use for the gather and scatter ops.
        //
        // print locally owned dof numbering
        /* print_dof<Type>(dm.get_global_dof_enum(), proc_num); */

        // print local dof numbering
        /* print_dof<Type>(dm.get_local_dof_enum(), proc_num); */
        /* print_dofs<Type>(dm, proc_num); */

        // print the global dofs for each element's local dof
        // print_elem_global_dof(dm, proc_num);

        // auto fe = build_fe_dof_map(dof_handler);

        // PVTUMeshWriter<DMQ2, Type> w;
        // std::string path = "io_test.vtu";
        // w.write_vtu(path, dm, fe);

        //     double time = timer.seconds();
        //     std::cout << "Setup took: " << time << " seconds." << std::endl;

        //     // use as more readable tuple index to identify the data
        //     constexpr int u = 0;
        //     constexpr int v = 1;
        //     // local dof enum size
        //     const Integer dof_size = dm.get_dof_size();

        //     // initialize the values by iterating through local dofs
        //     Kokkos::parallel_for(
        //         "initdatavalues", dof_size, MARS_LAMBDA(const Integer i) {
        //           dm.get_dof_data<INPUT>(i) = 1.0;
        //           // dm.get_dof_data<INPUT>(i) = i;
        //           dm.get_dof_data<OUTPUT>(i) = 0.0;
        //         });

        //     // specify the tuple indices of the tuplelements that are needed to gather.
        //     // if no index specified it gathers all views of the tuple. All data.
        //     dm.gather_ghost_data<INPUT>();

        //     if (proc_num == 0) {
        //       std::cout << "form_operator..." << std::flush;
        //     }

        //     form_operator<INPUT, OUTPUT>(dm, proc_num);

        //     if (proc_num == 0) {
        //       std::cout << "DONE" << std::endl;
        //     }

        //     // iterate through the local dofs and print the local number and the data
        //     /* dm.dof_iterate(
        //         MARS_LAMBDA(const Integer i) {
        //             printf("lid: %li, u: %lf, v: %lf, rank: %i\n", i,
        //                    dm.get_dof_data<u>(i), dm.get_dof_data<v>(i), proc_num);
        //         }); */

        //     scatter_add_ghost_data<OUTPUT>(dm, context);

        //     std::cout << "[" << proc_num
        //               << "] ndofs : " << dm.get_local_dof_enum().get_elem_size()
        //               << std::endl;

        //     // dm.dof_iterate(MARS_LAMBDA(const Integer i) {
        //     //   printf("ggid: %li, INPUT: %lf, OUTPUT: %lf, rank: %i\n", i,
        //     //          dm.get_dof_data<INPUT>(i), dm.get_dof_data<OUTPUT>(i),
        //     //          proc_num);
        //     // });

        //     time = timer.seconds();
        //     std::cout << "Matrix free took: " << time << " seconds." << std::endl;

#endif
    } catch (std::exception &e) {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
    }
}
