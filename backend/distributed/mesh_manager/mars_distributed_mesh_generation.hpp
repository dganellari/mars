#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#include "mars_err.hpp"
#ifdef MARS_ENABLE_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef MARS_ENABLE_KOKKOS
#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "KokkosKernels_Sorting.hpp"
#endif
#include "mars_distributed_mesh_kokkos.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
    using DMesh =
        Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>, KeyType>;

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
    void build_gp_np(const V &first_sfc, const V &GpNp, const M &GpNp_host, const Integer last_sfc) {
        using namespace Kokkos;
        auto size = first_sfc.extent(0);
        auto first_sfc_mirror = create_mirror_view(first_sfc);
        deep_copy(first_sfc_mirror, first_sfc);

        for (int i = 0; i < size; ++i) {
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
    inline void build_first_sfc(const VW &elem_sfc, const V &gp_np) {
        auto size = gp_np.extent(0) / 2;
        Kokkos::parallel_for(
            size, KOKKOS_LAMBDA(const Integer i) {
                auto index = gp_np(2 * i + 1);
                // The last index is the total size of the sfc vector
                int inc = 0;
                if (i == size - 1) {
                    --index;
                    ++inc;
                }
                gp_np(2 * i) = elem_sfc(index) + inc;
            });
    }

    template <typename V, typename L>
    inline void compact_into_local(const V &sfc_to_local,
                                   const ViewVectorType<bool> all_elements,
                                   const L &local,
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

    template <typename V, typename S, typename H>
    void compact_into_local(const V &elem_sfc, const S &local, const H &gp, const Integer rank) {
        using namespace Kokkos;

        const Integer start = gp(2 * rank + 1);
        const Integer end = gp(2 * (rank + 1) + 1);
        parallel_for(
            RangePolicy<>(start, end), KOKKOS_LAMBDA(const Integer i) {
                const Integer index = i - start;
                local(index) = elem_sfc(i);
            });
    }

    template<class SfcKeyType>
    struct GenerateSFC {
        using KeyType = typename SfcKeyType::ValueType;

        ViewVectorType<bool> predicate;
        KeyType max_oc;
        Integer xDim;
        Integer yDim;
        Integer zDim;

        GenerateSFC(ViewVectorType<bool> el, KeyType m, Integer xdm, Integer ydm)
            : predicate(el), max_oc(m), xDim(xdm), yDim(ydm) {}
        GenerateSFC(ViewVectorType<bool> el, KeyType m, Integer xdm, Integer ydm, Integer zdm)
            : predicate(el), max_oc(m), xDim(xdm), yDim(ydm), zDim(zdm) {}
        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const {
            const auto enc_oc = encode_sfc_2D<SfcKeyType>(i, j);
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
            const auto enc_oc = encode_sfc_3D<SfcKeyType>(i, j, k);
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

    template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
    ViewVectorType<typename KeyType::ValueType> generate_elements_sfc(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh) {
        using namespace Kokkos;
        const Integer xDim = mesh.get_XDim();
        const Integer yDim = mesh.get_YDim();
        const Integer zDim = mesh.get_ZDim();

        Integer number_of_elements = get_number_of_elements(mesh);
        ViewVectorType<typename KeyType::ValueType> element_sfc("element_sfc", number_of_elements);

        switch (Type) {
            case ElementType::Quad4: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim == 0);

                Kokkos::parallel_for(
                    MDRangePolicy<Rank<2>>({0, 0}, {xDim, yDim}), MARS_LAMBDA(const int i, const int j) {
                        const auto oc = encode_sfc_2D<KeyType>(i, j);
                        const auto index = xDim * j + i;
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
                        const auto oc = encode_sfc_3D<KeyType>(i, j, k);
                        const auto index = xDim * (j + k * yDim) + i;
                        element_sfc(index) = oc;
                    });
                break;
            }
        }

        return element_sfc;
    }

    template <Integer Type, class KeyType>
    auto generate_sfc(const SFC<Type, KeyType> &morton) {
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

                const auto max_oc = encode_sfc_2D<KeyType>(xDim, yDim);
                parallel_for(MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}),
                             GenerateSFC<KeyType>(all_elements, max_oc, xDim, yDim));
                break;
            }
            case ElementType::Hex8: {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);

                const auto max_oc = encode_sfc_3D<KeyType>(xDim, yDim, zDim);
                parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                             GenerateSFC<KeyType>(all_elements, max_oc, xDim, yDim, zDim));
                break;
            }
        }

        return all_elements;
    }

    template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
    Integer get_number_of_elements(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh) {
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
                errorx(
                    1,
                    "Unknown elemnt type for the mesh generation. Implemented only for Quad4 and Hex8 element types.");
            }
        }
        return number_of_elements;
    }

//write the generate_sfc function assuming you have a square with quad4 elements
//which does not need sorting and the sfc code can be generated directly from the
//element indices. The sfc code is generated using the hilbert code and


//write the print_mesh_info function to print the mesh information
void print_mesh_info(Integer num_elements, Integer chunk_size, Integer last_chunk_size) {
    std::cout << "Mesh Number of Elements: " << num_elements << std::endl;
    std::cout << "Mesh Chunk Size: " << chunk_size << std::endl;
    std::cout << "Mesh Last Chunk Size: " << last_chunk_size << std::endl;
}

//write the calculate_counts function to calculate the counts for each mpi process
std::vector<int> calculate_counts(int size, int chunk_size, int last_chunk_size) {
    std::vector<int> counts(size);
    for (int i = 0; i < size; ++i) {
        if (i == size - 1)
            counts[i] = last_chunk_size;
        else
            counts[i] = chunk_size;
    }
    return counts;
}

#ifdef MARS_ENABLE_CUDA
// This function processes the mesh using CUDA. It's only compiled if MARS_ENABLE_CUDA is defined.
template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType, class H>
void process_with_cuda(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh, SFC<Type, KeyType> &sfc_generator, const H &global_partition_host, int rank) {
    // Generate SFC elements for each mesh element
    auto elem_sfc = generate_elements_sfc(mesh);
    // Sort the SFC elements using radix sort
    auto sorted_elem_sfc = buffer_cub_radix_sort(elem_sfc);
    // Build the first SFC using the sorted elements
    build_first_sfc(sorted_elem_sfc, mesh.get_view_gp());
    // Copy the global partition to the mesh's view of the global partition
    deep_copy(global_partition_host, mesh.get_view_gp());
    // Compact the sorted elements into the local partition
    compact_into_local(sorted_elem_sfc, sfc_generator.get_view_elements(), global_partition_host, rank);
}
#else
// This function processes the mesh without using CUDA. It's only compiled if MARS_ENABLE_CUDA is not defined.
template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType, class H>
void process_without_cuda(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh, SFC<Type, KeyType> &sfc_generator, const H &global_partition_host, int size, int rank) {
    // Create a vector to hold the first SFC per rank
    ViewVectorType<Integer> first_sfc("first_sfc_per_rank", size);
    // Generate all SFC elements
    auto all_sfc_elements_predicate = generate_sfc<Type, KeyType>(sfc_generator);
    // Create a vector to map SFC to local elements
    ViewVectorType<Integer> sfc_to_local("sfc_to_local_mesh_generation", sfc_generator.get_all_range());
    // Compact the SFC to local elements
    compact_into_local(
        sfc_to_local, all_sfc_elements_predicate, sfc_generator.get_view_elements(), mesh.get_view_gp(), rank);

    // Build the first SFC using the local elements
    build_first_sfc(sfc_to_local, all_sfc_elements_predicate, first_sfc, mesh.get_view_gp(), size);
    // Get the total range of the SFC
    auto all_range = sfc_generator.get_all_range();
    // Build the global partition using the first SFC and the total range
    build_gp_np(first_sfc, mesh.get_view_gp(), global_partition_host, all_range - 1);
}
#endif

template <typename IntegerType, int D>
IntegerType generateHilbertCurve(int L) {
    // The number of points on the curve is 2^(DL)
    return pow(2, D*L);
}

template <Integer Type, typename Oc, typename T>
MARS_INLINE_FUNCTION bool is_inside_subdomain(const Oc &coords,
                                            const T domain_max_x,
                                            const T domain_max_y,
                                            const T domain_max_z) {
    Integer domain_min_x = 0, domain_min_y = 0, domain_min_z = 0;
    // Check if the coordinates are inside the domain
    if constexpr (Type == ElementType::Quad4) {
        return coords.x >= domain_min_x && coords.x < domain_max_x && coords.y >= domain_min_y &&
               coords.y < domain_max_y;
    } else if constexpr (Type == ElementType::Hex8) {
        return coords.x >= domain_min_x && coords.x < domain_max_x && coords.y >= domain_min_y &&
               coords.y < domain_max_y && coords.z >= domain_min_z && coords.z < domain_max_z;
    }
}

template <Integer Dim, Integer ManifoldDim, Integer Type, typename KeyType>
void process_segment(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh,
                    const typename KeyType::ValueType start,
                    const typename KeyType::ValueType end) {
    using ValueType = typename KeyType::ValueType;
    static_assert(std::is_same<KeyType, HilbertKey<ValueType>>::value, "Key type must be HilbertKey");

    const auto rank_elements_size = end - start;
    auto rank_elements_predicate = ViewVectorType<bool>("rank_elements", rank_elements_size);
    auto rank_elements_scan = ViewVectorType<ValueType>("rank_elements_scan", rank_elements_size + 1);

    const Integer xDim = mesh.get_XDim();
    const Integer yDim = mesh.get_YDim();
    const Integer zDim = mesh.get_ZDim();
    // Each rank works on its segment of the Hilbert curve
    Kokkos::parallel_for(
        "workOnSegment", Kokkos::RangePolicy<>(start, end), MARS_LAMBDA(const Integer i) {
            auto coords = get_octant_from_sfc<Type, KeyType>(i);
            if (is_inside_subdomain<Type>(coords, xDim, yDim, zDim)) {
                rank_elements_predicate(i - start) = true;
            }
        });
    // Wait for all parallel operations to complete
    Kokkos::fence();

    // Perform an inclusive exclusive scan
    incl_excl_scan(0, rank_elements_size, rank_elements_predicate, rank_elements_scan);

    auto index_subview = subview(rank_elements_scan, rank_elements_size);
    auto index_host = create_mirror_view(index_subview);

    // Deep copy device view to host view.
    deep_copy(index_host, index_subview);
    auto local_size = index_host();
    auto rank_elements = ViewVectorType<ValueType>("rank_elements", local_size);

    printf("Hilbert chunk_size: %i, Full Hilbert size: %li\n", local_size, rank_elements_size);

    Kokkos::parallel_for(
        "compactElements", Kokkos::RangePolicy<>(start, end), KOKKOS_LAMBDA(const ValueType i) {
            auto index = i - start;
            if (rank_elements_predicate(index)) {
                auto local_index = rank_elements_scan(index);
                rank_elements(local_index) = i;
            }
        });

    mesh.set_view_sfc(rank_elements);
    mesh.set_chunk_size(local_size);
    mesh.generate_sfc_to_local_map();
}

// to aovid predicate, scan and compact and so reduce memory footprint when full hilbert curve is processed
template <Integer Dim, Integer ManifoldDim, Integer Type, typename KeyType>
void process_full_segment(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh,
                    const typename KeyType::ValueType start,
                    const typename KeyType::ValueType end) {
    using ValueType = typename KeyType::ValueType;
    static_assert(std::is_same<KeyType, HilbertKey<ValueType>>::value, "Key type must be HilbertKey");

    const auto rank_elements_size = end - start;
    auto rank_elements = ViewVectorType<ValueType>("rank_elements", rank_elements_size);

    printf("Hilbert chunk_size: %li\n", rank_elements_size);

    Kokkos::parallel_for(
        "compactElements", Kokkos::RangePolicy<>(start, end), KOKKOS_LAMBDA(const ValueType i) {
            auto index = i - start;
            rank_elements(index) = i;
        });

    mesh.set_view_sfc(rank_elements);
    mesh.set_chunk_size(rank_elements_size);
    mesh.generate_sfc_to_local_map();
}

template <typename V>
void build_first_hilbert_sfc(const V &first_sfc, int size, Integer portion_size) {
    // Each rank calculates its first SFC in parallel
    Kokkos::parallel_for(
        "update_first_sfc", Kokkos::RangePolicy<>(0, size + 1), KOKKOS_LAMBDA(const Integer rank) {
            // Compute the start of the SFC for this rank
            Integer start = rank * portion_size;
            // Update the start in the view
            first_sfc(2 * rank) = start;
        });

    // Wait for all parallel operations to complete
    Kokkos::fence();
}

template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
bool is_full_hilbert_curve(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh) {
    // Check if all dimensions are equal
    if (mesh.get_XDim() != mesh.get_YDim() || mesh.get_XDim() != mesh.get_YDim()) {
        return false;
    }
    // Calculate L based on one of the dimensions
    int L = log2(mesh.get_XDim());
    // Check if the dimension is a power of 2
    if (pow(2, L) != mesh.get_XDim()) {
        return false;
    }
    // Get the number of elements in the mesh
    Integer num_elements = get_number_of_elements(mesh);
    // Generate the Hilbert curve range for this process
    auto num_points = generateHilbertCurve<Integer, Dim>(L);
    // Check if the number of elements is equal to the number of points
    return num_elements == num_points;
}

template <Integer Dim, Integer ManifoldDim, Integer Type, class IntegerType>
void partition_mesh(DMesh<Dim, ManifoldDim, Type, HilbertKey<IntegerType>> &mesh) {
    Integer num_elements = get_number_of_elements(mesh);
    printf("partition_mesh for hilbert curve with %li elements\n", num_elements);
    const context &context = mesh.get_context();
    int rank = mars::rank(context);
    mesh.set_proc(rank);
    int size = num_ranks(context);
    int max_dim = std::max({mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim()});
    int L = ceil(log2(max_dim));
    // Generate the Hilbert curve range for this process
    auto num_points = generateHilbertCurve<IntegerType, Dim>(L);

    // Calculate the portion of the mesh for this process
    IntegerType portion_size = num_points / size;
    IntegerType start = rank * portion_size;
    IntegerType end = (rank == size - 1) ? num_points : start + portion_size;

    // Check if the Hilbert curve is full. If it is, we can process the segment without predicate, scan and compact
    /* if (is_full_hilbert_curve(mesh)) {
        process_full_segment(mesh, start, end);
    } else { */
        process_segment(mesh, start, end);
    // }

    ViewVectorType<Integer> global_partition = ViewVectorType<Integer>("global_static_partition", 2 * (size + 1));
    mesh.set_view_gp(global_partition);

    build_first_hilbert_sfc(mesh.get_view_gp(), size, portion_size);
}

template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
void partition_mesh(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh) {
    using namespace Kokkos;
    Timer timer;

    const context &context = mesh.get_context();
    int rank = mars::rank(context);
    mesh.set_proc(rank);
    int size = num_ranks(context);
    Integer num_elements = get_number_of_elements(mesh);
    Integer chunk_size = num_elements / size;
    Integer last_chunk_size = chunk_size - (chunk_size * size - num_elements);

    printf("chunk_size: %li, number_of_elements: %li, rank: %i\n", chunk_size, num_elements, rank);

    if (chunk_size <= 0) {
        errorx(1, " Invalid number of mpi processes. Defined more mpi processes than mesh elements to be generated!");
    }

    bool is_root = mars::rank(context) == 0;
    if (is_root) {
        print_mesh_info(num_elements, chunk_size, last_chunk_size);
    }

    std::vector<int> counts = calculate_counts(size, chunk_size, last_chunk_size);
    if (rank == size - 1) {
        chunk_size = last_chunk_size;
    }
    mesh.set_chunk_size(chunk_size);

    ViewVectorType<Integer> global_partition = ViewVectorType<Integer>("global_static_partition", 2 * (size + 1));
    auto global_partition_host = scan_count(global_partition, counts, size);
    mesh.set_view_gp(global_partition);

    SFC<Type, KeyType> sfc_generator(mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim());
    sfc_generator.reserve_elements(mesh.get_chunk_size());

    Timer time;
#ifdef MARS_ENABLE_CUDA
    process_with_cuda(mesh, sfc_generator, global_partition_host, rank);
#else
    process_without_cuda(mesh, sfc_generator, global_partition_host, size, rank);
#endif
    double time_radix = time.seconds();
    std::cout << "Radix sort method took: " << time_radix<< " seconds." << std::endl;

    mesh.set_view_sfc(sfc_generator.get_view_elements());

    sfc_generator.generate_sfc_to_local_map();
    mesh.set_sfc_to_local_map(sfc_generator.get_sfc_to_local_map());

    double time_gen = timer.seconds();
    std::cout << "SFC Generation and Partition took: " << time_gen << " seconds. Process: " << rank << std::endl;
}

    // write the same partition_mesh function for the distributed mesh with the difference that the sfc elements are generated independtly using
//the hilbert curve taking advantage of the property of the hilbert curve that the elements are generated in a sorted order and there are no jumps
//in the sfc code which allows independent generation of the sfc code for each mpi process. What the current implementation does is each process to generate the
//sfc code for the whole mesh and then each process takes it cut. This is not efficient and it is not scalable. The sfc code should be generated independently
//for each mpi process and then the partitioning should be done.
// the points and elements can be generated on the fly from the sfc code in case meshless true.

    template <class KeyType, Integer Dim, Integer ManifoldDim, Integer Type, bool Meshless = true>
    void generate_distributed_cube(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh,
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

        printf("xDim %i, yDim %i, zDim %i\n", xDim, yDim, zDim);
        // partition the mesh and then generate the points and elements.
        partition_mesh(mesh);

        mesh.template create_ghost_layer<Type>();

        if (!Meshless) {
            Kokkos::Timer timer_gen;

            /* the mesh construct depends on template parameters. */
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
