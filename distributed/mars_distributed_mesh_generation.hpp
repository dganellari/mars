#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "KokkosKernels_Sorting.hpp"
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"
namespace mars {

    template <Integer Dim, Integer ManifoldDim, Integer Type>
    using DMesh = Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;

    template <typename V>
    auto scan_count(const V &GpNp, const std::vector<int> &counts, const Integer size) {
        using namespace Kokkos;
        auto GpNp_host = create_mirror_view(GpNp);

        for (int i = 0; i < size; ++i) {
            // counts[0] is always chunk size. Kind of scan with the first and the last elem of the sfc
            // acc sum scan for the count
            GpNp_host(2 * (i + 1) + 1) = GpNp_host(2 * i + 1) + counts[i];
        }

        deep_copy(GpNp, GpNp_host);

        return GpNp_host;

        /* parallel_for(
           "print_elem_gp:", size+1, KOKKOS_LAMBDA(const int i) {
           printf(" elch: (%li-%li) - %i - %i\n", GpNp(2*i),  GpNp(2*i+1), i, proc_num);
           }); */
    }

    template <typename V, typename M>
    void build_gp_np(const V &first_sfc,
                     const V &GpNp,
                     const M &GpNp_host,
                     const Integer last_sfc) {
        using namespace Kokkos;
        Unsigned size = first_sfc.extent(0);
        auto first_sfc_mirror = create_mirror_view(first_sfc);
        deep_copy(first_sfc_mirror, first_sfc);

        for (Unsigned i = 0; i < size; ++i) {
            GpNp_host(2 * i) = first_sfc_mirror(i);
        }
        /* insert the last element of the sfc adding 1 to it to make the last element not part of the linearization.
        In this way the binary search works properly. */
        GpNp_host(2 * size) = last_sfc;
        deep_copy(GpNp, GpNp_host);
    }

    template <typename V>
    MARS_INLINE_FUNCTION bool inside_proc_local(const Integer local, const V &gp, const Integer rank) {
        return (local >= gp(2 * rank + 1) && local < gp(2 * (rank + 1) + 1));
    }

    template <typename V>
    MARS_INLINE_FUNCTION Integer
    get_local_index(const bool is_sfc, const Integer local, const V &gp, const Integer rank) {
        Integer index = -1;

        if (is_sfc && inside_proc_local(local, gp, rank)) {
            index = local - gp(2 * rank + 1);
        }
        return index;
    }

    template <typename V>
    MARS_INLINE_FUNCTION Integer
    get_first_sfc_rank(const bool is_sfc, const Integer local, const V &gp, const Integer size) {
        Integer index = -1;

        for (Integer i = 0; i < size; ++i) {
            if (is_sfc && local == gp(2 * i + 1)) {
                index = i;
            }
        }

        return index;
    }

    template <typename V>
    inline void build_first_sfc(const V &sfc_to_local,
                                const ViewVectorType<bool> all_elements,
                                const V &first_sfc,
                                const V &gp_np,
                                const Integer size) {
        Kokkos::parallel_for(
            sfc_to_local.extent(0), KOKKOS_LAMBDA(const Integer i) {
                const Integer index = get_first_sfc_rank(all_elements(i), sfc_to_local(i), gp_np, size);
                if (index >= 0) {
                    first_sfc(index) = i;
                }
            });
    }

    template <typename VW, typename V>
    inline void build_first_sfc(const VW &elem_sfc, const V &first_sfc, const V &gp_np) {
        auto size = first_sfc.extent(0);
        Kokkos::parallel_for(
            size, KOKKOS_LAMBDA(const Integer i) {
                const Integer index = gp_np(2 * i + 1);
                first_sfc(i) = elem_sfc(index);
            });
    }

    template <typename V>
    inline void compact_into_local(const V &sfc_to_local,
                                   const ViewVectorType<bool> all_elements,
                                   const V &local,
                                   const V &gp_np,
                                   const Integer rank) {
        using namespace Kokkos;

        auto all_range = sfc_to_local.extent(0);
        exclusive_bool_scan(0, all_range, sfc_to_local, all_elements);

        parallel_for(
            all_range, KOKKOS_LAMBDA(const Integer i) {
                const Integer index = get_local_index(all_elements(i), sfc_to_local(i), gp_np, rank);
                if (index >= 0) {
                    local(index) = i;
                }
            });
    }

    template <typename V, typename M>
    ViewVectorType<Integer> compact_into_local(const V &elem_sfc, M &mesh) {
        using namespace Kokkos;

        ViewVectorType<Integer> local("local_partition_sfc", mesh.get_chunk_size());
        auto rank = mesh.get_proc();

        auto gp = mesh.get_view_gp();
        const Integer start = gp(2 * rank + 1);
        const Integer end = gp(2 * (rank + 1) + 1);
        parallel_for(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const Integer i) {
                const Integer index = i - start;
                local(index) = elem_sfc(i);
            });
        return local;
    }

    struct GenerateSFC {
        ViewVectorType<bool> predicate;
        Integer max_oc;
        Integer xDim;
        Integer yDim;
        Integer zDim;

        GenerateSFC(ViewVectorType<bool> el, Integer m, Integer xdm, Integer ydm)
            : predicate(el), max_oc(m), xDim(xdm), yDim(ydm) {}
        GenerateSFC(ViewVectorType<bool> el, Integer m, Integer xdm, Integer ydm, Integer zdm)
            : predicate(el), max_oc(m), xDim(xdm), yDim(ydm), zDim(zdm) {}
        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const {
            const Integer enc_oc = encode_morton_2D(i, j);
            // set to true only those elements from the vector that are generated.
            assert(enc_oc < max_oc);
            /* if (enc_oc >= max_oc) {
                const Integer chunk_size = xDim * yDim;
                printf("You have reached the mesh generation limit size. Can not generate %li elements mesh\n",
                       chunk_size);
                exit(1);
            } */
            predicate(enc_oc) = 1;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int k, int j, int i) const {
            const Integer enc_oc = encode_morton_3D(i, j, k);
            // set to true only those elements from the vector that are generated.
            assert(enc_oc < max_oc);
            /* if (enc_oc >= max_oc) {
                const Integer chunk_size = xDim * yDim * zDim;
                printf("You have reached the mesh generation limit size. Can not generate %li elements mesh\n",
                       chunk_size);
                exit(1);
            } */
            predicate(enc_oc) = 1;
        }
    };

    template <Integer Type, typename M>
    auto generate_local_sfc(M &mesh, const std::vector<int> &counts) {
        auto rank_size = num_ranks(mesh.get_context());

        ViewVectorType<Integer> GpNp = ViewVectorType<Integer>("global_static_partition", 2 * (rank_size + 1));
        auto GpNp_host = scan_count(GpNp, counts, rank_size);
        mesh.set_view_gp(GpNp);

        auto elem_sfc = generate_elements_sfc(mesh);

        auto local = compact_into_local(elem_sfc, mesh);
        mesh.set_view_sfc(local);

        ViewVectorType<Integer> first_sfc("first_sfc_per_rank", rank_size);
        build_first_sfc(elem_sfc, first_sfc, mesh.get_view_gp());

        auto last_sfc = compute_all_range<Type>(mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim());
        build_gp_np(first_sfc, mesh.get_view_gp(), GpNp_host, last_sfc - 1);
    }

    using unsigned_l = unsigned long;

    template <Integer Dim, Integer ManifoldDim, Integer Type>
    ViewVectorType<unsigned_l> generate_elements_sfc(DMesh<Dim, ManifoldDim, Type> &mesh) {
        using namespace Kokkos;
        const Integer xDim = mesh.get_XDim();
        const Integer yDim = mesh.get_YDim();
        const Integer zDim = mesh.get_ZDim();


        Integer number_of_elements = get_number_of_elements(mesh);
        ViewVectorType<unsigned_l> element_sfc("element_sfc", number_of_elements);

        switch (Type) {
            case ElementType::Quad4: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim == 0);

                Kokkos::parallel_for(
                    MDRangePolicy<Rank<2>>({0, 0}, {xDim, yDim}), MARS_LAMBDA(const int i, const int j) {
                        const Integer oc = encode_morton_2D(i, j);
                        const Integer index = xDim * j + i;
                        element_sfc(index) = oc;
                    });
                break;
            }
            case ElementType::Hex8: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);

                Kokkos::parallel_for(
                    MDRangePolicy<Rank<3>>({0, 0, 0}, {xDim, yDim, zDim}),
                    MARS_LAMBDA(const int i, const int j, const int k) {
                        const Integer oc = encode_morton_3D(i, j, k);
                        const Integer index = xDim * (j + k * yDim) + i;
                        element_sfc(index) = oc;
                    });
                break;
            }
        }

        ViewVectorType<unsigned_l> aux_elem_sfc("aux_elem_sfc", number_of_elements);
        /* Kokkos::Impl::sort(element_sfc, 0, element_sfc.extent(0)); */
        KokkosKernels::Impl::SerialRadixSort<Integer, unsigned_l>(
            element_sfc.data(), aux_elem_sfc.data(), number_of_elements);

        return element_sfc;
    }

    template <Integer Type>
    auto generate_sfc(const SFC<Type> &morton) {
        using namespace Kokkos;

        const Integer xDim = morton.get_XDim();
        const Integer yDim = morton.get_YDim();
        const Integer zDim = morton.get_ZDim();

        ViewVectorType<bool> all_elements("all_sfc_predicate", morton.get_all_range());

        switch (Type) {
            case ElementType::Quad4: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim == 0);

                const Integer max_oc = encode_morton_2D(xDim, yDim);
                parallel_for(MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}),
                             GenerateSFC(all_elements, max_oc, xDim, yDim));
                break;
            }
            case ElementType::Hex8: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);

                const Integer max_oc = encode_morton_3D(xDim, yDim, zDim);
                parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                             GenerateSFC(all_elements, max_oc, xDim, yDim, zDim));
                break;
            }
        }

        return all_elements;
    }

    template <Integer Dim, Integer ManifoldDim, Integer Type>
    Integer get_number_of_elements(DMesh<Dim, ManifoldDim, Type> &mesh) {
        Integer number_of_elements = 0;

        switch (Type) {
            case ElementType::Quad4: {
                number_of_elements = mesh.get_XDim() * mesh.get_YDim();
                break;
            }
            case ElementType::Hex8: {
                number_of_elements = mesh.get_XDim() * mesh.get_YDim() * mesh.get_ZDim();
                break;
            }
            default: {
                errx(1,
                     "Unknown elemnt type for the mesh generation. Implemented only for Quad4 and Hex8 element types.");
            }
        }
        return number_of_elements;
    }

    template <Integer Dim, Integer ManifoldDim, Integer Type>
    void partition_mesh(DMesh<Dim, ManifoldDim, Type> &mesh) {
        using namespace Kokkos;

        Timer timer;

        const context &context = mesh.get_context();

        int proc_num = rank(context);
        mesh.set_proc(proc_num);
        // std::cout << "rank -:    " << proc_num << std::endl;
        int size = num_ranks(context);
        // std::cout << "size - :    " << size << std::endl;

        Integer number_of_elements = get_number_of_elements(mesh);

        // Integer chunk_size = (Integer)ceil((double)number_of_elements / size);
        // Integer chunk_size = number_of_elements / size + (number_of_elements % size != 0);
        // Integer last_chunk_size = chunk_size - (chunk_size * size - number_of_elements);
        Integer chunk_size = number_of_elements / size;
        Integer last_chunk_size = chunk_size - (chunk_size * size - number_of_elements);

        printf("chunk_size: %li, number_of_elements: %li, rank: %i\n", chunk_size, number_of_elements, proc_num);
        assert(chunk_size > 0);

        if (chunk_size <= 0) {
            errx(1, " Invalid number of mpi processes. Defined more mpi processes than mesh elements to be generated!");
        }

        bool root = mars::rank(context) == 0;
        if (root) {
            std::cout << "Mesh Number of Elements: " << number_of_elements << std::endl;
            std::cout << "Mesh Chunk Size: " << chunk_size << std::endl;
            std::cout << "Mesh Last Chunk Size: " << last_chunk_size << std::endl;
        }

        std::vector<int> counts(size);
        for (int i = 0; i < size; ++i) {
            if (i == size - 1)
                counts[i] = last_chunk_size;
            else
                counts[i] = chunk_size;
        }
        // set the chunk size to the remainder for the last mpi processes.
        if (proc_num == size - 1) {
            chunk_size = last_chunk_size;
        }
        mesh.set_chunk_size(chunk_size);

        ViewVectorType<Integer> GpNp = ViewVectorType<Integer>("global_static_partition", 2 * (size + 1));
        auto GpNp_host = scan_count(GpNp, counts, size);
        mesh.set_view_gp(GpNp);

        // generate the SFC linearization
        SFC<Type> morton(mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim());
        auto all_sfc_elements_predicate = generate_sfc(morton);
        morton.reserve_elements(mesh.get_chunk_size());

        // compacting sfc predicate and inserting the morton code corresponding to a true value in the predicate
        // leaves the sfc elements array sorted.
        /* ViewVectorType<Integer> local("local_partition_sfc", mesh.get_chunk_size()); */
        ViewVectorType<Integer> sfc_to_local("sfc_to_local_mesh_generation", morton.get_all_range());
        compact_into_local(
            sfc_to_local, all_sfc_elements_predicate, morton.get_view_elements(), mesh.get_view_gp(), proc_num);

        ViewVectorType<Integer> first_sfc("first_sfc_per_rank", size);
        build_first_sfc(sfc_to_local, all_sfc_elements_predicate, first_sfc, mesh.get_view_gp(), size);

        auto all_range = morton.get_all_range();
        build_gp_np(first_sfc, mesh.get_view_gp(), GpNp_host, all_range - 1);

        mesh.set_view_sfc(morton.get_view_elements());
        morton.generate_sfc_to_local_map();
        mesh.set_sfc_to_local_map(morton.get_sfc_to_local_map());

        double time_gen = timer.seconds();
        std::cout << "SFC Generation and Partition took: " << time_gen << " seconds. Process: " << proc_num
                  << std::endl;

        /* //Classical method using the radix sort. 10x slower than the current one due to the sort.
        generate_local_sfc<Type>(mesh, counts); */
    }

    // the points and elements can be generated on the fly from the sfc code in case meshless true.
    template <Integer Dim, Integer ManifoldDim, Integer Type, bool Meshless = true>
    void generate_distributed_cube(DMesh<Dim, ManifoldDim, Type> &mesh,
                                   const Integer xDim,
                                   const Integer yDim,
                                   const Integer zDim) {
        using namespace Kokkos;

        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        const context &context = mesh.get_context();

        int proc_num = rank(context);

        mesh.set_XDim(xDim);
        mesh.set_YDim(yDim);
        mesh.set_ZDim(zDim);

        // partition the mesh and then generate the points and elements.
        partition_mesh(mesh);

        mesh.template create_ghost_layer<Type>();

        if (!Meshless) {
            Kokkos::Timer timer_gen;

            // the mesh construct depends on template parameters.
            auto gen_pts = mesh.template generate_points<Type>();

            auto gen_elm = mesh.template generate_elements<Type>();

            double time_gen = timer_gen.seconds();
            std::cout << "Distributed Mesh Generation took: " << time_gen << " seconds. Process: " << proc_num
                      << std::endl;

            if (!gen_pts || !gen_elm) std::cerr << "Not implemented for other dimensions yet" << std::endl;
        }
    }
}  // namespace mars

#endif
#endif
#endif  // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
