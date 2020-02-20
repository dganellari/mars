#ifndef MARS_DIST_MESH_KOKKOS_HPP
#define MARS_DIST_MESH_KOKKOS_HPP

#include "mars_point.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "mars_distributed_non_simplex_kokkos.hpp"
#include "mars_distributed_simplex_kokkos.hpp"
#include "mars_imesh_kokkos.hpp"

#ifdef WITH_MPI
#include "mars_sfc_generation.hpp"
#endif
#include "mars_utils_kokkos.hpp"

namespace mars
{

template <Integer Dim_, Integer ManifoldDim_, class Simplex_>
class Mesh<Dim_, ManifoldDim_, DistributedImplementation, Simplex_> : public ParallelIMesh<Dim_>
{
public:
    static constexpr Integer Dim = Dim_;
    static constexpr Integer ManifoldDim = ManifoldDim_;

    using Elem = Simplex_;
    using Point = mars::Point<Real, Dim>;
    using Comb = Combinations<ManifoldDim + 1, 2, DistributedImplementation>;

    MARS_INLINE_FUNCTION Mesh()
        : ParallelIMesh<Dim_>(), elements_size_(0), points_size_(0)
    //, combinations(nullptr)
    {
    }

    void reserve(const std::size_t n_elements, const std::size_t n_points) override
    {
        elements_size_ = n_elements;
        points_size_ = n_points;

        elements_ = ViewMatrixType<Integer>("elems", n_elements, Elem::ElemType);
        active_ = ViewVectorType<bool>("active_", n_elements);
        points_ = ViewMatrixType<Real>("pts", n_points, Dim);
    }

    void reserve_points(const std::size_t n_points)
    {
        points_size_ = n_points;
        points_ = ViewMatrixType<Real>("pts", n_points, Dim);
    }

    void reserve_elements(const std::size_t n_elements)
    {
        elements_size_ = n_elements;
        elements_ = ViewMatrixType<Integer>("elems", n_elements,
                                            Elem::ElemType);
        active_ = ViewVectorType<bool>("active_", n_elements);
    }

    /*inline Elem &elem(const Integer id) override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return simplices_(id);
	}

	inline const Elem &elem(const Integer id) const override
	{
		assert(id >= 0);
		assert(id < n_elements());
		return simplices_(id);
	}*/

    MARS_INLINE_FUNCTION Elem elem(const Integer id) // override
    {
        assert(id >= 0);
        assert(id < n_elements());
        Elem e = Elem(SubView<Integer, Elem::ElemType>(&elements_, id), combinations);
        e.id = id;
        return e;
    }

    MARS_INLINE_FUNCTION const Elem elem(const Integer id) const //override
    {
        assert(id >= 0);
        assert(id < n_elements());
        Elem e = Elem(SubView<Integer, Elem::ElemType>(&elements_, id), combinations);
        e.id = id;
        return e;
    }

    MARS_INLINE_FUNCTION bool is_active(const Integer id) const override
    {
        assert(id >= 0);
        assert(id < n_elements());
        return active_(id);
    }

    MARS_INLINE_FUNCTION bool is_valid(const Integer id) const
    {
        return id >= 0 && id < n_elements();
    }

    MARS_INLINE_FUNCTION bool is_node_valid(const Integer id) const
    {
        return id >= 0 && id < n_nodes();
    }

    MARS_INLINE_FUNCTION
    const ViewMatrixTextureC<Integer, Comb::value, 2> &combs() const
    {
        return combinations;
    }

    MARS_INLINE_FUNCTION
    void set_combs(const ViewMatrixTextureC<Integer, Comb::value, 2> &c)
    {
        combinations = c;
    }

    MARS_INLINE_FUNCTION
    void set_active(const Integer id, const bool val = true)
    {
        assert(id >= 0);
        assert(id < n_elements());
        active_(id) = val;
    }

    MARS_INLINE_FUNCTION
    void add_point(const Point &point, const Integer index)
    {
        for (Integer i = 0; i < Dim; ++i)
        {
            points_(index, i) = point[i];
        }
    }

    MARS_INLINE_FUNCTION
    void add_point(const TempArray<Real, Dim> &point, const Integer index)
    {
        for (Integer i = 0; i < Dim; ++i)
        {
            points_(index, i) = point[i];
        }
    }

    MARS_INLINE_FUNCTION
    Point point(const Integer i) override
    {
        assert(i >= 0);
        assert(i < points_size_);
        return Point(points_, i);
    }

    MARS_INLINE_FUNCTION
    const Point point(const Integer i) const override
    {
        assert(i >= 0);
        assert(i < points_size_);
        return Point(points_, i);
    }

    const ViewMatrixType<Real> &points() const //override
    {
        return points_;
    }

    ViewMatrixType<Real> get_view_points() const //override
    {
        return points_;
    }

    ViewMatrixType<Integer> get_view_elements() const //override
    {
        return elements_;
    }

    MARS_INLINE_FUNCTION
    ViewVectorType<bool> get_view_active() const //override
    {
        return active_;
    }

    MARS_INLINE_FUNCTION
    const ViewVectorType<unsigned int> &get_view_sfc() const
    {
        return local_sfc_;
    }

    MARS_INLINE_FUNCTION
    void set_view_sfc(const ViewVectorType<unsigned int> &local)
    {
        local_sfc_ = local;
    }

    void resize_points(const Integer size)
    {
        points_size_ += size;
        resize(points_, points_size_, Dim);
    }

    void resize_elements(const Integer size)
    {
        elements_size_ += size;
        resize(elements_, elements_size_, Elem::ElemType);
        resize(active_, elements_size_);
    }

    void resize_active(const Integer size)
    {
        resize(active_, elements_size_ + size);
    }

    MARS_INLINE_FUNCTION
    Integer n_nodes() const override
    {
        return points_size_;
    }

    MARS_INLINE_FUNCTION
    void set_n_nodes(Integer size)
    {
        points_size_ = size;
    }

    MARS_INLINE_FUNCTION
    Integer n_elements() const override
    {
        return elements_size_;
    }

    MARS_INLINE_FUNCTION
    void set_n_elements(Integer size)
    {
        elements_size_ = size;
    }

    MARS_INLINE_FUNCTION
    void set_chunk_size(unsigned int size)
    {
        chunk_size_ = size;
    }

    MARS_INLINE_FUNCTION
    unsigned int get_chunk_size()
    {
        return chunk_size_;
    }

    MARS_INLINE_FUNCTION
    void set_global_to_local_map(const UnorderedMap<unsigned int,unsigned int>& gl_map)
    {
        global_to_local_map_ = gl_map;
    }

    inline Integer n_active_elements(const Integer N)
    {
        using namespace Kokkos;

        Timer timer;

        ViewVectorType<bool> active_ = get_view_active();

        double result = 0;
        parallel_reduce(
            "Active_all", N, KOKKOS_LAMBDA(const int &i, double &lsum) {
                lsum += active_(i);
            },
            result);

        double time = timer.seconds();
        std::cout << "Active_all Reduction took: " << time << " seconds." << std::endl;

        printf("Result: %i %lf\n", N, result);

        return result;
    }

    inline Integer n_active_elements(const ViewVectorType<Integer> elements)
    {
        using namespace Kokkos;

        Timer timer;

        ViewVectorType<bool> active_ = get_view_active();

        Integer N = elements.extent(0);

        double result = 0;
        parallel_reduce(
            "Active_elem", N, KOKKOS_LAMBDA(const int &i, double &lsum) {
                lsum += active_(elements(i));
            },
            result);

        double time = timer.seconds();
        std::cout << "Active_elements Reduction took: " << time << " seconds." << std::endl;

        printf("Result: %i %lf\n", N, result);

        return result;
    }

    MARS_INLINE_FUNCTION Integer type() const override
    {
        return ManifoldDim;
    }

    //add point functor
    struct AddNonSimplexPoint
    {
        ViewMatrixType<Real> points;
        ViewVectorType<unsigned int> encoded_points;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm) : points(pts), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm, Integer zdm) : points(pts), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        AddNonSimplexPoint(ViewMatrixType<Real> pts, ViewVectorType<unsigned int> epts, Integer xdm, Integer ydm) : points(pts), encoded_points(epts), xDim(xdm), yDim(ydm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(unsigned int index) const
        {
            unsigned int gl_index = encoded_points(index);

            points(index, 0) = static_cast<Real>(decode_morton_2DX(gl_index)) / static_cast<Real>(xDim);
            points(index, 1) = static_cast<Real>(decode_morton_2DY(gl_index)) / static_cast<Real>(yDim);
        }

        /* KOKKOS_INLINE_FUNCTION
        void operator()(int z, int y, int x) const
        {

            int index = (xDim + 1) * (yDim + 1) * z + (xDim + 1) * y + x;

            points(index, 0) = static_cast<Real>(x) / static_cast<Real>(xDim);
            points(index, 1) = static_cast<Real>(y) / static_cast<Real>(yDim);
            points(index, 2) = static_cast<Real>(z) / static_cast<Real>(zDim);
        } */
    };

    template <Integer Type>
    struct BuildGlobalPointIndex
    {

        ViewVectorType<bool> predicate;
        ViewVectorType<unsigned int> global;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        BuildGlobalPointIndex(ViewVectorType<bool> el, Integer xdm, Integer ydm) : predicate(el), xDim(xdm), yDim(ydm)
        {
        }
        BuildGlobalPointIndex(ViewVectorType<bool> el, Integer xdm, Integer ydm, Integer zdm) : predicate(el), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        BuildGlobalPointIndex(ViewVectorType<bool> el, ViewVectorType<unsigned int> gl,
                              Integer xdm, Integer ydm) : predicate(el), global(gl), xDim(xdm), yDim(ydm)
        {
        }
        KOKKOS_INLINE_FUNCTION
        void operator()(int i) const
        {
            switch (Type)
            {
            case ElementType::Quad4:
            {
                // set to true only those elements from the vector that are generated.
                // in this way the array is already sorted and you just compact it using scan which is much faster in parallel.

                unsigned int gl_index = global(i);
                assert(gl_index < encode_morton_2D(xDim + 1, yDim + 1));

                // gl_index defines one element by its corner node. Then we add all the other nodes for that element.
                const unsigned int x = decode_morton_2DX(gl_index);
                const unsigned int y = decode_morton_2DY(gl_index);
                predicate(gl_index) = 1;
                predicate(encode_morton_2D(x + 1, y)) = 1;
                predicate(encode_morton_2D(x, y + 1)) = 1;
                predicate(encode_morton_2D(x + 1, y + 1)) = 1;
                break;

            }
            case ElementType::Hex8:
            {
                // set to true only those elements from the vector that are generated.
                // in this way the array is already sorted and you just compact it using scan which is much faster in parallel.

                unsigned int gl_index = global(i);
                assert(gl_index < encode_morton_3D(xDim + 1, yDim + 1, zDim + 1));

                // gl_index defines one element by its corner node. Then we add all the other nodes for that element.
                const unsigned int x = decode_morton_3DX(gl_index);
                const unsigned int y = decode_morton_3DY(gl_index);
                const unsigned int z = decode_morton_3DZ(gl_index);

                predicate(gl_index) = 1;
                predicate(encode_morton_3D(x + 1, y, z)) = 1;
                predicate(encode_morton_3D(x, y + 1, z)) = 1;
                predicate(encode_morton_3D(x + 1, y + 1, z)) = 1;

                predicate(encode_morton_3D(x , y , z + 1)) = 1;
                predicate(encode_morton_3D(x + 1, y, z + 1)) = 1;
                predicate(encode_morton_3D(x, y + 1, z + 1)) = 1;
                predicate(encode_morton_3D(x + 1, y + 1, z + 1)) = 1;
                break;
            }
            }
        }
    };

    template <Integer Type>
    inline ViewVectorType<bool> build_global_elements(const unsigned int allrange, 
                const int xDim, const int yDim, const int zDim)
    {
        using namespace Kokkos;

        Timer timer;

        ViewVectorType<bool> all_elements("predicate", allrange);

        parallel_for("local_sfc_range", get_chunk_size(),
                     BuildGlobalPointIndex<Type>(all_elements, local_sfc_, xDim, yDim));

        return all_elements;                     
    }

    template <Integer Type>
    inline unsigned int compact_elements(ViewVectorType<unsigned int> &ltg, const unsigned int allrange, 
                const int xDim, const int yDim, const int zDim)
    {
        using namespace Kokkos;

        Timer timer;

        const ViewVectorType<bool>& all_elements = build_global_elements<Type>(allrange, xDim, yDim, zDim);
        ViewVectorType<unsigned int> scan_indices("scan_indices", allrange + 1);

        incl_excl_scan(0, allrange, all_elements, scan_indices);

        auto index_subview = subview(scan_indices, allrange);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        //std::cout << "Hyper count result: " << h_ic(0)<< std::endl;

        ltg = ViewVectorType<unsigned int>("local_to_global", h_ic());

        parallel_for(
            allrange, KOKKOS_LAMBDA(const unsigned int i) {
                if (all_elements(i) == 1)
                {
                    unsigned int k = scan_indices(i);
                    ltg(k) = i;
                }
            });

        return h_ic();
    }

    template <Integer Type>
    inline bool generate_points(const int xDim, const int yDim, const int zDim)
    {
        using namespace Kokkos;

        switch (ManifoldDim_)
        {

        case 2:
        {
            switch (Type)
            {

            case ElementType::Quad4:
            {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim == 0);

                ViewVectorType<unsigned int> local_to_global;

                const unsigned int allrange = encode_morton_2D(xDim + 1, yDim + 1); //TODO : check if enough. Test with xdim != ydim.

                const unsigned int nr_points = compact_elements<Type>(local_to_global, allrange, xDim, yDim, zDim);
                printf("nr_p: %u\n", nr_points);

                reserve_points(nr_points);

                parallel_for("non_simplex_point", nr_points,
                             AddNonSimplexPoint(points_, local_to_global, xDim, yDim));

                //build global_to_local map for use in generate elements.
                UnorderedMap<unsigned int, unsigned int> global_to_local_map(nr_points);

                parallel_for(
                    "local_global_map_for", nr_points, KOKKOS_LAMBDA(const int i) {
                        const int offset = xDim + 1;

                        const unsigned int gl_index = local_to_global(i);

                        const unsigned int x = decode_morton_2DX(gl_index);
                        const unsigned int y = decode_morton_2DY(gl_index);

                        global_to_local_map.insert(x + offset * y, i);
                    });

                set_global_to_local_map(global_to_local_map);

                return true;
            }
            default:
            {
                return false;
            }
            }
        }
        case 3:
        {
            switch (Type)
            {
            case ElementType::Hex8:
            {
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);

               /*  ViewVectorType<unsigned int> local_to_global;

                unsigned int nr_points = compact_elements<Type>(local_to_global, xDim, yDim, zDim);
                printf("nr_p: %u\n", nr_points); */

                /*   const int n_nodes = (xDim + 1) * (yDim + 1) * (zDim + 1);

                reserve_points(n_nodes);

                parallel_for(
                    MDRangePolicy<Rank<3>>({0, 0, 0},
                                           {zDim + 1, yDim + 1, xDim + 1}),
                    AddNonSimplexPoint(points_, xDim, yDim, zDim)); */
                return true;
            }
            default:
            {
                return false;
            }
            }
        }
        default:
        {
            return false;
        }
        }
    }

    //add elem functor
    struct AddNonSimplexElem
    {

        using UMap = UnorderedMap<unsigned int, unsigned int>;
        ViewMatrixType<Integer> elem;
        ViewVectorType<bool> active;
        ViewVectorType<unsigned int> global;
        UMap map;
        
        Integer xDim;
        Integer yDim;
        Integer zDim;

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
                          Integer xdm) : elem(el), active(ac), xDim(xdm)
        {
        }

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
                          Integer xdm, Integer ydm) : elem(el), active(ac), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<unsigned int> gl, UMap mp,
                          ViewVectorType<bool> ac, Integer xdm, Integer ydm) : 
                                    elem(el), global(gl), map(mp), active(ac), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
                          Integer xdm, Integer ydm, Integer zdm) : elem(el), active(ac), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int index) const
        {
            const int offset = xDim + 1;

            const unsigned int gl_index = global(index);

            const unsigned int i = decode_morton_2DX(gl_index);
            const unsigned int j = decode_morton_2DY(gl_index);

            elem(index, 0) = map.value_at(map.find(i + offset * j));
            elem(index, 1) = map.value_at(map.find((i + 1) + offset * j));
            elem(index, 2) = map.value_at(map.find((i + 1) + offset * (j + 1)));
            elem(index, 3) = map.value_at(map.find(i + offset * (j + 1)));

            active(index) = true;
        }

        /*   KOKKOS_INLINE_FUNCTION
        void operator()(int k, int j, int i) const
        {
            int index = k * xDim * yDim + j * xDim + i;

            elem(index, 0) = elem_index(i, j, k, xDim, yDim);
            elem(index, 1) = elem_index(i + 1, j, k, xDim, yDim);
            elem(index, 2) = elem_index(i + 1, j + 1, k, xDim, yDim);
            elem(index, 3) = elem_index(i, j + 1, k, xDim, yDim);
            elem(index, 4) = elem_index(i, j, k + 1, xDim, yDim);
            elem(index, 5) = elem_index(i + 1, j, k + 1, xDim, yDim);
            elem(index, 6) = elem_index(i + 1, j + 1, k + 1, xDim, yDim);
            elem(index, 7) = elem_index(i, j + 1, k + 1, xDim, yDim);

            active(index) = true;
        } */
    };

    inline bool generate_elements(const int xDim, const int yDim,
                                  const int zDim, Integer type)
    {

        using namespace Kokkos;

        switch (ManifoldDim_)
        {

        case 2:
        {

            switch (type)
            {

            case ElementType::Quad4:
            {
                reserve_elements(get_chunk_size());

                parallel_for("generate_elements", get_chunk_size(),
                             AddNonSimplexElem(elements_, local_sfc_, global_to_local_map_, active_, xDim, yDim));
                return true;
            }
            default:
            {
                return false;
            }
            }
        }
        case 3:
        {

            switch (type)
            {

            case ElementType::Hex8:
            {
                /*    const int n_elements = xDim * yDim * zDim;
                reserve_elements(n_elements);

                parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                             AddNonSimplexElem(elements_, active_, xDim, yDim, zDim)); */
                return true;
            }
            default:
            {
                return false;
            }
            }
        }
        default:
        {
            return false;
        }
        }
    }

private:
    ViewMatrixTextureC<Integer, Comb::value, 2> combinations;

    ViewMatrixType<Integer> elements_;
    ViewMatrixType<Real> points_;
    ViewVectorType<bool> active_;
    Integer elements_size_;
    Integer points_size_;

    ViewVectorType<unsigned int> local_sfc_;
    unsigned int chunk_size_;

    UnorderedMap<unsigned int, unsigned int> global_to_local_map_;
};

using DistributedMesh1 = mars::Mesh<1, 1, DistributedImplementation, Simplex<1, 1, DistributedImplementation>>;
using DistributedMesh2 = mars::Mesh<2, 2, DistributedImplementation, Simplex<2, 2, DistributedImplementation>>;
using DistributedMesh3 = mars::Mesh<3, 3, DistributedImplementation, Simplex<3, 3, DistributedImplementation>>;
using DistributedMesh4 = mars::Mesh<4, 4, DistributedImplementation, Simplex<4, 4, DistributedImplementation>>;

using DistributedQuad4Mesh = mars::Mesh<2, 2, DistributedImplementation, Quad4DElem>;
using DistributedHex8Mesh = mars::Mesh<3, 3, DistributedImplementation, Hex8DElem>;

template <Integer Type>
using DistributedNSMesh2 = mars::Mesh<2, 2, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;
} // namespace mars
#endif //MARS_MESH_HPP
