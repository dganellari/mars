#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_

#ifdef WITH_MPI
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"
namespace mars
{

template <Integer Dim, Integer ManifoldDim, Integer Type>
using DMesh = Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;

template <typename V>
void scan_count(const V &GpNp, const std::vector<int> &counts, const Integer size) {
    auto GpNp_host = create_mirror_view(GpNp);

    for (int i = 0; i < size; ++i) {
        // counts[0] is always chunk size. Kind of scan with the first and the last elem of the sfc
        // acc sum scan for the count
        GpNp_host(2 * (i + 1) + 1) = GpNp_host(2 * i + 1) + counts[i];
    }
    // insert the last element of the sfc adding 1 to it (to make the last element not part of the linearization) for
    // the binary search to work properly
    deep_copy(GpNp, GpNp_host);

    /* parallel_for(
       "print_elem_gp:", size+1, KOKKOS_LAMBDA(const int i) {
       printf(" elch: (%li-%li) - %i - %i\n", GpNp(2*i),  GpNp(2*i+1), i, proc_num);
       }); */
}

/* template <Integer Dim, Integer ManifoldDim, Integer Type>
void broadcast_gp_np(const context &context,
                     const std::vector<int> &counts,
                     ViewVectorType<Integer> elems,
                     Integer n__anchor_nodes) {
    using namespace Kokkos;
    int proc_num = rank(context);
    int size = num_ranks(context);

    auto GpNp_host = create_mirror_view(GpNp);

    auto local = subview(elems, 0);
    auto local_host = create_mirror_view(local);
    deep_copy(local_host, local);

    [>mpi_all_gather();<]
    for (int i = 0; i < size; ++i) {
        // counts[0] is always chunk size. Kind of scan with the first and the last elem of the sfc
        GpNp_host(2 * i) = elem_view_host(i * counts[0]);
    }
    // insert the last element of the sfc adding 1 to it (to make the last element not part of the linearization) for
    // the binary search to work properly
    GpNp_host(2 * size) = elem_view_host(n__anchor_nodes - 1) + 1;
    deep_copy(GpNp, GpNp_host);

    [>parallel_for(
       "print_elem_gp:", size+1, KOKKOS_LAMBDA(const int i) {
       printf(" elch: (%li-%li) - %i - %i\n", GpNp(2*i),  GpNp(2*i+1), i, proc_num);
       });<]
} */

struct GenerateSFC {
    ViewVectorType<bool> predicate;
    Integer xDim;
    Integer yDim;
    Integer zDim;

    GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm) : predicate(el), xDim(xdm), yDim(ydm) {}
    GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm, Integer zdm)
        : predicate(el), xDim(xdm), yDim(ydm), zDim(zdm) {}
    KOKKOS_INLINE_FUNCTION
    void operator()(int j, int i) const {
        // set to true only those elements from the vector that are generated.
        // in this way the array is already sorted and you just compact it using scan which is much faster in
        // parallel.
        assert(encode_morton_2D(i, j) < encode_morton_2D(xDim, yDim));
        if (encode_morton_2D(i, j) >= encode_morton_2D(xDim, yDim)) {
            const Integer chunk_size = xDim * yDim;
            printf("You have reached the mesh genration limit size. Can not generate mesh %li elements\n", chunk_size);
            exit(1);
        }

        predicate(encode_morton_2D(i, j)) = 1;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int k, int j, int i) const {
        // set to true only those elements from the vector that are generated.
        // in this way the array is already sorted and you just compact it using scan which is much faster in
        // parallel.
        assert(encode_morton_3D(i, j, k) < encode_morton_3D(xDim, yDim, zDim));
        if (encode_morton_3D(i, j, k) >= encode_morton_3D(xDim, yDim, zDim)) {
            const Integer chunk_size = xDim * yDim * zDim;
            printf("You have reached the mesh genration limit size. Can not generate mesh %li elements\n", chunk_size);
            exit(1);
        }

        predicate(encode_morton_3D(i, j, k)) = 1;
    }
};

template <typename V>
MARS_INLINE_FUNCTION Integer
get_first_sfc_rank(const bool is_sfc, const Integer sfc_to_local, const V &gp, const Integer rank) {
    const Integer index = -1;
    bool inside_rank_segment = (sfc_to_local >= gp(2 * rank + 1) && sfc_to_local < gp(2 * (rank + 1) + 1));

    if (is_sfc && inside_rank_segment) {
        index = sfc_to_local - gp(2 * rank + 1);
    }

    return index;
}

template <typename V>
MARS_INLINE_FUNCTION Integer
get_local_index(const bool is_sfc, const Integer sfc_to_local, const V &gp, const Integer size) {
    const Integer index = -1;

    for(Integer i = 0; i< size; ++i) {
        if (is_fsc && sfc_to_local == gp(2 * i + 1)) {
            index = i;
        }
    }

    return index;
}

template <Integer Type, typename V>
inline void build_gp_np(const SFC<Type> &morton,
                        const ViewVectorType<bool> all_elements,
                        const V &first_sfc,
                        const V &gp_np,
                        const Integer size) {
    using namespace Kokkos;

    ViewVectorType<Integer> sfc_to_local = morton.get_view_sfc_to_local();
    parallel_for(
        morton.get_all_range(), KOKKOS_LAMBDA(const Integer i) {
            const Integer index = get_first_sfc_rank(all_elements(i), sfc_to_local(i), gp_np, size);
            if (index >= 0) {
                first_sfc(index) = i;
            }
        });
}

template <Integer Type, typename V>
inline void compact_into_local(const SFC<Type> &morton,
                               const ViewVectorType<bool> all_elements,
                               const V &local,
                               const V &gp_np,
                               const Integer rank) {
    using namespace Kokkos;

    exclusive_bool_scan(0, morton.get_all_range(), morton.get_view_sfc_to_local(), all_elements);

    // otherwise kokkos lambda will not work with CUDA
    ViewVectorType<Integer> sfc_to_local = morton.get_view_sfc_to_local();
    parallel_for(
        morton.get_all_range(), KOKKOS_LAMBDA(const Integer i) {
            const Integer index = get_local_index(all_elements(i), sfc_to_local(i), gp_np, rank);
            if (index >= 0) {
                local(index) = i;
            }
        });
}

template <Integer Type, typename V>
inline void generate_sfc(const SFC<Type> &morton, const V &local, const V &GpNp, const context &context) {
    using namespace Kokkos;

    int rank = rank(context);
    int size = num_ranks(context);

    const Integer xDim = morton.get_XDim();
    const Integer yDim = morton.get_YDim();
    const Integer zDim = morton.get_ZDim();

    ViewVectorType<bool> all_elements("predicate", morton.get_all_range());
    Integer n__anchor_nodes;

    switch (Type) {
        case ElementType::Quad4: {
            assert(xDim != 0);
            assert(yDim != 0);
            assert(zDim == 0);

            n__anchor_nodes = xDim * yDim;

            parallel_for(MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}), GenerateSFC(all_elements, xDim, yDim));
            break;
        }
        case ElementType::Hex8: {
            assert(xDim != 0);
            assert(yDim != 0);
            assert(zDim != 0);

            n__anchor_nodes = xDim * yDim * zDim;
            parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                         GenerateSFC(all_elements, xDim, yDim, zDim));
            break;
        }
        default: {
            return false;
        }
    }
    // compacting the 1 and 0 array and inserting the "true" index of the all elements
    // which is the correct morton code leaving the sfc elements array sorted.
    compact_into_local(morton, all_elements, local, GpNp, rank);
    /* assert(n__anchor_nodes == get_elem_size()); */

    V first_sfc("first_sfc_per_rank", size);
    build_gp_np(morton, all_elements, first_sfc, GpNp, size);
    auto first_sfc_mirror = create_mirror_view(first_sfc);
    deep_copy(first_sfc_mirror, first_sfc);

    //TODO: Fill the gpnp host from first_sfc mirror and deep copy it on the device. The current data of the gpnp host
    //mirror should be reused.

}

template <Integer Dim, Integer ManifoldDim, Integer Type>
void partition_mesh(const context &context, DMesh<Dim, ManifoldDim, Type> &mesh)
{
    using namespace Kokkos;

    Kokkos::Timer timer;

    int proc_num = rank(context);
    // std::cout << "rank -:    " << proc_num << std::endl;

    int size = num_ranks(context);
    // std::cout << "size - :    " << size << std::endl;

    Integer n__anchor_nodes = 0;

    switch (Type)
    {
    case ElementType::Quad4:
    {
        n__anchor_nodes = mesh.get_XDim() * mesh.get_YDim();
        break;
    }
    case ElementType::Hex8:
    {
        n__anchor_nodes = mesh.get_XDim() * mesh.get_YDim() * mesh.get_ZDim();
        break;
    }
    default:
    {
        std::cout << "Not yet implemented for other element types!" << std::endl;
        return;
    }
    }

    //Integer chunk_size = (Integer)ceil((double)n__anchor_nodes / size);
    //Integer chunk_size = n__anchor_nodes / size + (n__anchor_nodes % size != 0);
    //Integer last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);
    Integer chunk_size = n__anchor_nodes / size;
    Integer last_chunk_size = chunk_size - (chunk_size * size - n__anchor_nodes);

    printf("chunk_size: %li, n__anchor_nodes: %li, rank: %i\n", chunk_size, n__anchor_nodes, proc_num);
    assert(chunk_size > 0);

    if (chunk_size <= 0)
    {
        errx(1, " Invalid number of mpi processes. Defined more mpi processes than mesh elements to be generated!");
    }

    bool root = mars::rank(context) == 0;
    if (root)
    {
        std::cout << "chunk_size - :    " << chunk_size << std::endl;
        std::cout << "n__anchor_nodes size:: - :    " << n__anchor_nodes << std::endl;
        std::cout << "last_chunk_size size:: - :    " << last_chunk_size << std::endl;
    }

    std::vector<int> counts(size);

    for (int i = 0; i < size; ++i)
    {
        if (i == size - 1)
            counts[i] = last_chunk_size;
        else
            counts[i] = chunk_size;
    }


    ViewVectorType<Integer> GpNp = ViewVectorType<Integer>("global_static_partition", 2 * (size + 1));

    scan_count(GpNp, counts, size);

    //set the chunk size to the remainder for the last mpi processes.
    if (proc_num == size - 1)
    {
        chunk_size = last_chunk_size;
    }


    SFC<Type> morton(mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim());

    ViewVectorType<Integer> local("local_partition_sfc", chunk_size);
    generate_sfc(morton, local, GpNp, proc_num);

    /* context->distributed->scatterv_gids(morton.get_view_elements(), local, counts); */

    /* std::cout << "MPI Scatter ended!" << std::endl; */

    /* std::cout << "Broadcasting the sfc_to_local..." << std::endl; */

    /* context->distributed->broadcast(morton.get_view_sfc_to_local()); //broadcast to all processors. */
    /* std::cout << "Broadcasting the sfc_to_local ended!" << std::endl; */

    /*  parallel_for(
        "print_elem_chunk",chunk_size, KOKKOS_LAMBDA(const int i) {
        printf(" elch: %u-%i\n", local(i), proc_num);
        }); */

    mesh.set_view_gp(GpNp);
    mesh.set_view_sfc(local);
    mesh.set_view_sfc_to_local(morton.get_view_sfc_to_local());
    mesh.set_chunk_size(chunk_size);
    mesh.set_proc(proc_num);

    /* broadcast_gp_np(context, mesh, counts, local, n__anchor_nodes); */

    double time_gen = timer.seconds();
    std::cout << "Paritioning took: " << time_gen << " seconds. Process: " << proc_num << std::endl;
}

//the points and elements can be generated on the fly from the sfc code in case meshless true.
template <Integer Dim, Integer ManifoldDim, Integer Type, bool Meshless = true>
bool generate_distributed_cube(const context &context, DMesh<Dim, ManifoldDim, Type> &mesh,
                               const Integer xDim, const Integer yDim, const Integer zDim)
{
    using namespace Kokkos;

    assert(Dim <= 3);
    assert(ManifoldDim <= Dim);

    int proc_num = rank(context);

    mesh.set_XDim(xDim);
    mesh.set_YDim(yDim);
    mesh.set_ZDim(zDim);

    //partition the mesh and then generate the points and elements.
    partition_mesh(context, mesh);

    bool gen_pts = true, gen_elm = true;

    if (!Meshless)
    {
        Kokkos::Timer timer_gen;

        //the mesh construct depends on template parameters.
        gen_pts = mesh.template generate_points<Type>();

        gen_elm = mesh.template generate_elements<Type>();

        double time_gen = timer_gen.seconds();
        std::cout << "Distributed Generation kokkos took: " << time_gen << " seconds. Process: " << proc_num << std::endl;

        if (!gen_pts || !gen_elm)
            std::cerr << "Not implemented for other dimensions yet" << std::endl;
    }

    /* mesh.template build_boundary_element_sets<Type>(); */

    return (gen_pts && gen_elm);
}

} // namespace mars

#endif
#endif

#endif // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
