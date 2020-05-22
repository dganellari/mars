#ifndef MARS_DIST_MESH_KOKKOS_HPP
#define MARS_DIST_MESH_KOKKOS_HPP

#include "mars_base.hpp"
#include "mars_globals.hpp"
#include "mars_point.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "mars_distributed_non_simplex_kokkos.hpp"
#include "mars_distributed_octant.hpp"
#include "mars_distributed_simplex_kokkos.hpp"
#include "mars_imesh_kokkos.hpp"
#include "mars_distributed_utils.hpp"

#ifdef WITH_MPI
#include "mars_sfc_generation.hpp"
#endif
#include "mars_utils_kokkos.hpp"
#define ASSERT(left, operator, right)                                                                                                                                                                  \
    {                                                                                                                                                                                                  \
        if (!((left) operator(right)))                                                                                                                                                                 \
        {                                                                                                                                                                                              \
            std::cerr << "ASSERT FAILED: " << #left << #operator<< #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" <<(left) << "; " << #right << "=" << (right) << std::endl; \
        }                                                                                                                                                                                              \
    }

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
    const ViewVectorType<Integer> &get_view_sfc() const
    {
        return local_sfc_;
    }

    MARS_INLINE_FUNCTION
    void set_view_sfc(const ViewVectorType<Integer> &local)
    {
        local_sfc_ = local;
    }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> &get_view_gp() const
    {
        return gp_np;
    }

    MARS_INLINE_FUNCTION
    void set_view_gp(const ViewVectorType<Integer> &gp)
    {
        gp_np = gp;
    }


    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> &get_view_boundary() const
    {
        return boundary_;
    }

    MARS_INLINE_FUNCTION
    void set_view_boundary(const ViewVectorType<Integer> &b)
    {
        boundary_ = b;
    }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> &get_view_boundary_sfc_index() const
    {
        return boundary_lsfc_index_;
    }

    MARS_INLINE_FUNCTION
    void set_view_boundary_sfc_index(const ViewVectorType<Integer> &b)
    {
        boundary_lsfc_index_ = b;
    }

    MARS_INLINE_FUNCTION
    const ViewVectorType<Integer> &get_view_scan_boundary() const
    {
        return scan_boundary_;
    }

    MARS_INLINE_FUNCTION
    void set_view_scan_boundary(const ViewVectorType<Integer> &b)
    {
        scan_boundary_ = b;
    }

    MARS_INLINE_FUNCTION
    void set_XDim(const Integer xd)
    {
        xDim =xd;
    }

    MARS_INLINE_FUNCTION
    Integer get_XDim() const
    {
        return xDim;
    }

    MARS_INLINE_FUNCTION
    void set_YDim(const Integer yd)
    {
        yDim =yd;
    }

    MARS_INLINE_FUNCTION
    Integer get_YDim() const
    {
        return yDim;
    }

    MARS_INLINE_FUNCTION
    void set_ZDim(const Integer zd)
    {
        zDim =zd;
    }

    MARS_INLINE_FUNCTION
    Integer get_ZDim() const
    {
        return zDim;
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
    void set_chunk_size(Integer size)
    {
        chunk_size_ = size;
    }

    MARS_INLINE_FUNCTION
    Integer get_chunk_size()
    {
        return chunk_size_;
    }

    MARS_INLINE_FUNCTION
    void set_proc(Integer p)
    {
        proc = p;
    }

    MARS_INLINE_FUNCTION
    Integer get_proc()
    {
        return proc;
    }

    MARS_INLINE_FUNCTION
    void set_global_to_local_map(const UnorderedMap<Integer, Integer> &gl_map)
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
    template <Integer Type>
    struct AddNonSimplexPoint
    {
        ViewMatrixType<Real> points;
        ViewVectorType<Integer> encoded_points;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        AddNonSimplexPoint(ViewMatrixType<Real> pts, ViewVectorType<Integer> epts, Integer xdm, Integer ydm) : points(pts), encoded_points(epts), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexPoint(ViewMatrixType<Real> pts, ViewVectorType<Integer> epts, Integer xdm, Integer ydm, Integer zdm) : points(pts), encoded_points(epts), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(Integer index) const
        {
            Integer gl_index = encoded_points(index);

            switch (Type)
            {
            case ElementType::Quad4:
            {
                points(index, 0) = static_cast<Real>(decode_morton_2DX(gl_index)) / static_cast<Real>(xDim);
                points(index, 1) = static_cast<Real>(decode_morton_2DY(gl_index)) / static_cast<Real>(yDim);
                break;
            }
            case ElementType::Hex8:
            {
                points(index, 0) = static_cast<Real>(decode_morton_3DX(gl_index)) / static_cast<Real>(xDim);
                points(index, 1) = static_cast<Real>(decode_morton_3DY(gl_index)) / static_cast<Real>(yDim);
                points(index, 2) = static_cast<Real>(decode_morton_3DZ(gl_index)) / static_cast<Real>(zDim);
                break;
            }
            }
        }
    };

    template <Integer Type>
    struct BuildGlobalPointIndex
    {

        ViewVectorType<bool> predicate;
        ViewVectorType<Integer> global;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        BuildGlobalPointIndex(ViewVectorType<bool> el, ViewVectorType<Integer> gl,
                              Integer xdm, Integer ydm) : predicate(el), global(gl), xDim(xdm), yDim(ydm)
        {
        }

        BuildGlobalPointIndex(ViewVectorType<bool> el, ViewVectorType<Integer> gl,
                              Integer xdm, Integer ydm, Integer zdm) : predicate(el), global(gl), xDim(xdm), yDim(ydm), zDim(zdm)
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

                Integer gl_index = global(i);
                assert(gl_index < encode_morton_2D(xDim + 1, yDim + 1));

                // gl_index defines one element by its corner node. Then we add all the other nodes for that element.
                const Integer x = decode_morton_2DX(gl_index);
                const Integer y = decode_morton_2DY(gl_index);
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

                Integer gl_index = global(i);

                assert(gl_index < encode_morton_3D(xDim + 1, yDim + 1, zDim + 1));

                // gl_index defines one element by its corner node. Then we add all the other nodes for that element.
                const Integer x = decode_morton_3DX(gl_index);
                const Integer y = decode_morton_3DY(gl_index);
                const Integer z = decode_morton_3DZ(gl_index);

                predicate(gl_index) = 1;
                predicate(encode_morton_3D(x + 1, y, z)) = 1;
                predicate(encode_morton_3D(x, y + 1, z)) = 1;
                predicate(encode_morton_3D(x + 1, y + 1, z)) = 1;

                predicate(encode_morton_3D(x, y, z + 1)) = 1;
                predicate(encode_morton_3D(x + 1, y, z + 1)) = 1;
                predicate(encode_morton_3D(x, y + 1, z + 1)) = 1;
                predicate(encode_morton_3D(x + 1, y + 1, z + 1)) = 1;
                break;
            }
            }
        }
    };

    template <Integer Type>
    inline ViewVectorType<bool> build_global_elements(const Integer allrange,
                                                      const int xDim, const int yDim, const int zDim)
    {
        using namespace Kokkos;

        Timer timer;

        ViewVectorType<bool> all_elements("predicate", allrange);

        parallel_for("local_sfc_range", get_chunk_size(),
                     BuildGlobalPointIndex<Type>(all_elements, local_sfc_, xDim, yDim, zDim));

        return all_elements;
    }

    template <Integer Type>
    inline Integer compact_elements(ViewVectorType<Integer> &ltg, const Integer allrange,
                                         const int xDim, const int yDim, const int zDim)
    {
        using namespace Kokkos;

        Timer timer;

        const ViewVectorType<bool> &all_elements = build_global_elements<Type>(allrange, xDim, yDim, zDim);
        ViewVectorType<Integer> scan_indices("scan_indices", allrange + 1);

        incl_excl_scan(0, allrange, all_elements, scan_indices);

        auto index_subview = subview(scan_indices, allrange);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        //std::cout << "Hyper count result: " << h_ic(0)<< std::endl;

        ltg = ViewVectorType<Integer>("local_to_global", h_ic());

        parallel_for(
            allrange, KOKKOS_LAMBDA(const Integer i) {
                if (all_elements(i) == 1)
                {
                    Integer k = scan_indices(i);
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

                ViewVectorType<Integer> local_to_global;

                const Integer allrange = encode_morton_2D(xDim + 1, yDim + 1); //TODO : check if enough. Test with xdim != ydim.

                const Integer nr_points = compact_elements<Type>(local_to_global, allrange, xDim, yDim, zDim);
                printf("nr_p: %u\n", nr_points);

                reserve_points(nr_points);

                parallel_for("non_simplex_point", nr_points,
                             AddNonSimplexPoint<Type>(points_, local_to_global, xDim, yDim));

                //build global_to_local map for use in generate elements.
                UnorderedMap<Integer, Integer> global_to_local_map(nr_points);

                parallel_for(
                    "local_global_map_for", nr_points, KOKKOS_LAMBDA(const int i) {
                        const int offset = xDim + 1;

                        const Integer gl_index = local_to_global(i);

                        const Integer x = decode_morton_2DX(gl_index);
                        const Integer y = decode_morton_2DY(gl_index);

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
            break;
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

                ViewVectorType<Integer> local_to_global;
                const Integer allrange = encode_morton_3D(xDim + 1, yDim + 1, zDim + 1); //TODO : check if enough. Test with xdim != ydim.

                const Integer nr_points = compact_elements<Type>(local_to_global, allrange, xDim, yDim, zDim);
                printf("nr_p 3D: %u\n", nr_points);

                reserve_points(nr_points);

                parallel_for("non_simplex_point", nr_points,
                             AddNonSimplexPoint<Type>(points_, local_to_global, xDim, yDim, zDim));

                //build global_to_local map for use in generate elements.
                UnorderedMap<Integer, Integer> global_to_local_map(nr_points);

                parallel_for(
                    "local_global_map_for", nr_points, KOKKOS_LAMBDA(const int i) {
                        const Integer gl_index = local_to_global(i);

                        const Integer x = decode_morton_3DX(gl_index);
                        const Integer y = decode_morton_3DY(gl_index);
                        const Integer z = decode_morton_3DZ(gl_index);

                        global_to_local_map.insert(elem_index(x, y, z, xDim, yDim), i);
                    });

                set_global_to_local_map(global_to_local_map);
                return true;
            }
            default:
            {
                return false;
            }
            }
            break;
        }
        default:
        {
            return false;
        }
        }
    }

    //add elem functor
    template <Integer Type>
    struct AddNonSimplexElem
    {

        using UMap = UnorderedMap<Integer, Integer>;
        ViewMatrixType<Integer> elem;
        ViewVectorType<bool> active;
        ViewVectorType<Integer> global;
        UMap map;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<Integer> gl, UMap mp,
                          ViewVectorType<bool> ac, Integer xdm, Integer ydm) : elem(el), global(gl), map(mp), active(ac), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<Integer> gl, UMap mp,
                          ViewVectorType<bool> ac, Integer xdm, Integer ydm, Integer zdm) : elem(el), global(gl), map(mp), active(ac), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int index) const
        {

            switch (Type)
            {
            case ElementType::Quad4:
            {
                const int offset = xDim + 1;

                const Integer gl_index = global(index);

                const Integer i = decode_morton_2DX(gl_index);
                const Integer j = decode_morton_2DY(gl_index);

                elem(index, 0) = map.value_at(map.find(i + offset * j));
                elem(index, 1) = map.value_at(map.find((i + 1) + offset * j));
                elem(index, 2) = map.value_at(map.find((i + 1) + offset * (j + 1)));
                elem(index, 3) = map.value_at(map.find(i + offset * (j + 1)));

                active(index) = true;
                break;
            }
            case ElementType::Hex8:
            {
                const Integer gl_index = global(index);

                const Integer i = decode_morton_3DX(gl_index);
                const Integer j = decode_morton_3DY(gl_index);
                const Integer k = decode_morton_3DZ(gl_index);

                elem(index, 0) = map.value_at(map.find(elem_index(i, j, k, xDim, yDim)));
                elem(index, 1) = map.value_at(map.find(elem_index(i + 1, j, k, xDim, yDim)));
                elem(index, 2) = map.value_at(map.find(elem_index(i + 1, j + 1, k, xDim, yDim)));
                elem(index, 3) = map.value_at(map.find(elem_index(i, j + 1, k, xDim, yDim)));
                elem(index, 4) = map.value_at(map.find(elem_index(i, j, k + 1, xDim, yDim)));
                elem(index, 5) = map.value_at(map.find(elem_index(i + 1, j, k + 1, xDim, yDim)));
                elem(index, 6) = map.value_at(map.find(elem_index(i + 1, j + 1, k + 1, xDim, yDim)));
                elem(index, 7) = map.value_at(map.find(elem_index(i, j + 1, k + 1, xDim, yDim)));

                active(index) = true;
                break;
            }
            }
        }
    };

    template <Integer Type>
    inline bool generate_elements(const int xDim, const int yDim,
                                  const int zDim)
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
                reserve_elements(get_chunk_size());

                parallel_for("generate_elements", get_chunk_size(),
                             AddNonSimplexElem<Type>(elements_, local_sfc_, global_to_local_map_, active_, xDim, yDim));
                return true;
            }
            default:
            {
                return false;
            }
            }
            break;
        }
        case 3:
        {

            switch (Type)
            {

            case ElementType::Hex8:
            {
                reserve_elements(get_chunk_size());

                parallel_for("generate_elements", get_chunk_size(),
                             AddNonSimplexElem<Type>(elements_, local_sfc_, global_to_local_map_, active_,
                                                     xDim, yDim, zDim));
                return true;
            }
            default:
            {
                return false;
            }
            }
            break;
        }
        default:
        {
            return false;
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION static Octant face_nbh(const Octant &ref_octant, const int face,
                                                const Integer xDim, const Integer yDim, const Integer zDim)
    {

        switch (Type)
        {
        case ElementType::Quad4:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x = ref_octant.x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
            const Integer y = ref_octant.y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);

            Octant o(x, y);
            if (o.x < 0 || o.y < 0 || o.x >= xDim || o.y >= yDim)
                o.set_invalid();

            return o;
        }
        case ElementType::Hex8:
        {

            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x = ref_octant.x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
            const Integer y = ref_octant.y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);
            const Integer z = ref_octant.z + ((face == 4) ? -1 : (face == 5) ? 1 : 0);

            Octant o(x, y, z);
            if (o.x < 0 || o.y < 0 || o.z < 0 || o.x >= xDim || o.y >= yDim || o.z >= zDim)
                o.set_invalid();
            return o;
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION static Octant corner_nbh(const Octant &ref_octant, const int corner,
                                                  const Integer xDim, const Integer yDim, const Integer zDim)
    {

        switch (Type)
        {
        case ElementType::Quad4:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x = ref_octant.x + 2 * (corner & 1) - 1;
            const Integer y = ref_octant.y + (corner & 2) - 1;

            Octant o(x, y);
            if (o.x < 0 || o.y < 0 || o.x >= xDim || o.y >= yDim)
                o.set_invalid();

            return o;
        }
        case ElementType::Hex8:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x = ref_octant.x + 2 * (corner & 1) - 1;
            ;
            const Integer y = ref_octant.y + (corner & 2) - 1;
            const Integer z = ref_octant.z + (corner & 4) / 2 - 1;

            Octant o(x, y, z);
            if (o.x < 0 || o.y < 0 || o.z < 0 || o.x >= xDim || o.y >= yDim || o.z >= zDim)
                o.set_invalid();

            return o;
        }
        }
    }

    template <Integer Type>
    struct CountGhostNeighbors
    {
        ViewVectorType<Integer> global;
        ViewVectorType<Integer> count;
        ViewVectorType<Integer> gp;

        Integer proc;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        CountGhostNeighbors(ViewVectorType<Integer> gl, ViewVectorType<Integer> ct, ViewVectorType<Integer> g,
                            Integer p, Integer xdm, Integer ydm, Integer zdm) : global(gl), count(ct), gp(g), proc(p), xDim(xdm),
                                                                                yDim(ydm), zDim(zdm)
        {
        }

        MARS_INLINE_FUNCTION
        void increment(const Octant &o) const
        {
            Integer enc_oc = get_sfc_from_octant<Type>(o);

            assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
            Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

            //printf("Proc: %li, %li\n", proc, owner_proc);
            //the case when the neihgbor is a ghost element.
            if (proc != owner_proc)
            {
                Kokkos::atomic_increment(&count(owner_proc));
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(int index) const
        {
            const Integer gl_index = global(index);
            Octant ref_octant = get_octant_from_sfc<Type>(gl_index);

            const int offset = xDim + 1;

            for (int face = 0; face < 2 * ManifoldDim; ++face)
            {
                Octant o = face_nbh<Type>(ref_octant, face, xDim, yDim, zDim);
                if (o.is_valid())
                {

                    printf("face Nbh of %li (%li) is : %li with--- x and y: %li - %li\n", gl_index,
                            elem_index(ref_octant.x, ref_octant.y, ref_octant.z, xDim, yDim),elem_index(o.x, o.y, o.z, xDim, yDim), o.x, o.y);
                    increment(o);
                }
            }

            for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner)
            {
                Octant o = corner_nbh<Type>(ref_octant, corner, xDim, yDim, zDim);
                if (o.is_valid())
                {
                    printf("Corner Nbh of %li (%li) is : %li with--- x and y: %li - %li\n", gl_index,
                            elem_index(ref_octant.x, ref_octant.y, ref_octant.z, xDim, yDim),elem_index(o.x, o.y, o.z, xDim, yDim), o.x, o.y);
                    increment(o);
                }
            }
        }
    };

    template <Integer Type>
    struct BuildBoundarySets
    {
        ViewVectorType<Integer> global;
        ViewVectorType<Integer> gp;

        ViewVectorType<Integer> set;
        ViewVectorType<Integer> scan;
        //ViewVectorType<Integer> index;

        Integer proc;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        BuildBoundarySets(ViewVectorType<Integer> gl, ViewVectorType<Integer> g, ViewVectorType<Integer> st,
                          ViewVectorType<Integer> sc,                                                               //ViewVectorType<Integer> in,
                          Integer p, Integer xdm, Integer ydm, Integer zdm) : global(gl), gp(g), set(st), scan(sc), //index(in),
                                                                              proc(p), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        MARS_INLINE_FUNCTION
        void scatter(const Octant &o, const Integer ref) const
        {
            Integer enc_oc = get_sfc_from_octant<Type>(o);

            assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
            Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

            //the case when the neihgbor is a ghost element.
            if (proc != owner_proc)
            {
                // printf("Proc: %li, %li , %li\n", proc, owner_proc, i);
                //get the starting index before incrementing it.
                Integer i = Kokkos::atomic_fetch_add(&scan(owner_proc), 1);
                set(i) = ref;
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(int index) const
        {
            const Integer gl_index = global(index);
            Octant ref_octant = get_octant_from_sfc<Type>(gl_index);

            const int offset = xDim + 1;
            for (int face = 0; face < 2 * ManifoldDim; ++face)
            {
                Octant o = face_nbh<Type>(ref_octant, face, xDim, yDim, zDim);
                if (o.is_valid())
                {
                    /* printf("face Nbh of %u is %li: x and y: %li - %li\n",
                           ref_octant.x + offset * ref_octant.y,
                           o.x + offset * o.y, o.x, o.y); */
                    scatter(o, gl_index);
                }
            }

            for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner)
            {
                Octant o = corner_nbh<Type>(ref_octant, corner, xDim, yDim, zDim);
                if (o.is_valid())
                {
                    /* printf("Corner Nbh of %u is %li: x and y: %li - %li\n",
                           ref_octant.x + offset * ref_octant.y,
                           o.x + offset * o.y, o.x, o.y); */
                    scatter(o, gl_index);
                }
            }
        }
    };

    //another possible way of building the ghost layer directly without using the boundary layers and then send with mpi. This would be an mpi less method howerver it requires extra compute time for sorting"
    template <Integer Type>
    inline void build_ghost_element_sets(const int xDim, const int yDim,
                                         const int zDim)
    {
        using namespace Kokkos;

        const Integer rank_size = gp_np.extent(0) / 2 - 1;

        ViewVectorType<Integer> count("count_per_proc", rank_size);

        parallel_for("count_ghost_nbhs", get_chunk_size(),
                     CountGhostNeighbors<Type>(local_sfc_, count, gp_np,
                                               proc, xDim, yDim, zDim));

        parallel_for(
            "print acc", rank_size, KOKKOS_LAMBDA(const int i) {
                printf(" count : %i-%li\n", i, count(i));
            });

        ViewVectorType<Integer> scan("scan", rank_size + 1);
        incl_excl_scan(0, rank_size, count, scan);

        auto index_subview = subview(scan, rank_size);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);

        //the set containing all the ghost elements for all processes. Scan helps to identify part of the set array to which process it belongs.
        ViewVectorType<Integer> set("build_set", h_ic());

        parallel_for("build_set_kernel", get_chunk_size(),
                     BuildBoundarySets<Type>(local_sfc_, gp_np, set, scan,
                                             proc, xDim, yDim, zDim));
        parallel_for(
            "print set", h_ic(), KOKKOS_LAMBDA(const int i) {
                printf(" set : %i - %li\n", i, set(i));
            });
    }

    template <Integer Type>
    struct IdentifyBoundaryPerRank
    {
        ViewVectorType<Integer> global;
        ViewMatrixType<bool> predicate;
        ViewVectorType<Integer> gp;

        Integer proc;

        Integer xDim;
        Integer yDim;
        Integer zDim;

        IdentifyBoundaryPerRank(ViewVectorType<Integer> gl, ViewMatrixType<bool> pr, ViewVectorType<Integer> g,
                                Integer p, Integer xdm, Integer ydm, Integer zdm) : global(gl), predicate(pr), gp(g), proc(p), xDim(xdm),
                                                                                    yDim(ydm), zDim(zdm)
        {
        }

        MARS_INLINE_FUNCTION
        void setPredicate(const int index, const Octant &o) const
        {
            Integer enc_oc = get_sfc_from_octant<Type>(o);

            assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
            Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

            //printf("Proc: %li, %li\n", proc, owner_proc);
            //the case when the neihgbor is a ghost element.
            if (proc != owner_proc)
            {
                predicate(owner_proc, index) = 1;
            }
        }

        MARS_INLINE_FUNCTION
        void operator()(int index) const
        {
            const Integer gl_index = global(index);
            Octant ref_octant = get_octant_from_sfc<Type>(gl_index);

            const int offset = xDim + 1;

            for (int face = 0; face < 2 * ManifoldDim; ++face)
            {
                Octant o = face_nbh<Type>(ref_octant, face, xDim, yDim, zDim);
                if (o.is_valid())
                {
                    /*  printf("face Nbh of %li (%li) is : %li with--- x and y: %li - %li\n", gl_index,
                           ref_octant.x + offset * ref_octant.y, o.x + offset * o.y, o.x, o.y); */
                    setPredicate(index, o);
                }
            }

            for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner)
            {
                Octant o = corner_nbh<Type>(ref_octant, corner, xDim, yDim, zDim);
                if (o.is_valid())
                {
                    /* printf("corner Nbh of %li (%li) is : %li with--- x and y: %li - %li\n", gl_index,
                           ref_octant.x + offset * ref_octant.y, o.x + offset * o.y, o.x, o.y); */
                    setPredicate(index, o);
                }
            }
        }
    };

    inline void compact_boundary_elements(const ViewVectorType<Integer> scan_indices,
                                          const ViewMatrixType<bool> predicate, const ViewMatrixType<Integer> predicate_scan,
                                          const Integer rank_size)
    {
        using namespace Kokkos;

        Timer timer;

        /* the only way to work using the lambda instead of the functor on c++11 */
        ViewVectorType<Integer> boundary = boundary_;
        ViewVectorType<Integer> local_sfc = local_sfc_;
        ViewVectorType<Integer> boundary_lsfc_index = boundary_lsfc_index_;

        parallel_for(
            MDRangePolicy<Rank<2>>({0, 0}, {rank_size, chunk_size_}),
            KOKKOS_LAMBDA(const Integer i, const Integer j) {
                if (predicate(i, j) == 1)
                {
                    Integer index = scan_indices(i) + predicate_scan(i, j);
                    boundary(index) = local_sfc(j);
                    boundary_lsfc_index(index) = j;
                }
            });
    }

    template <Integer Type>
    inline void build_boundary_element_sets(const int xDim, const int yDim,
                                            const int zDim)
    {
        using namespace Kokkos;

        const Integer rank_size = gp_np.extent(0) / 2 - 1;

        ViewMatrixType<bool> rank_boundary("count_per_proc", rank_size, chunk_size_);
        scan_boundary_ = ViewVectorType<Integer>("scan_boundary_", rank_size + 1);

        parallel_for("IdentifyBoundaryPerRank", get_chunk_size(),
                     IdentifyBoundaryPerRank<Type>(local_sfc_, rank_boundary, gp_np,
                                                   proc, xDim, yDim, zDim));

        /* perform a scan for each row with the sum at the end for each rank */
        ViewMatrixType<Integer> rank_scan("rank_scan", rank_size, chunk_size_ + 1);
        for (int i = 0; i < rank_size; ++i)
        {
            if (i != proc)
            {
                auto row_predicate = subview(rank_boundary, i, ALL);
                auto row_scan = subview(rank_scan, i, ALL);
                incl_excl_scan(0, chunk_size_, row_predicate, row_scan);
/* 
                parallel_for(
                    "print scan", chunk_size_, KOKKOS_LAMBDA(const int i) {
                        printf(" boundary -inside: %i-%i", i, row_predicate(i));
                    });

                printf("\n"); */

            }
        }

        //perform a scan on the last column to get the total sum.
        column_scan(rank_size, chunk_size_, rank_scan, scan_boundary_);

        auto index_subview = subview(scan_boundary_, rank_size);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        std::cout << "boundary_ count result: " << h_ic() << std::endl;

        boundary_ = ViewVectorType<Integer>("boundary_", h_ic());
        boundary_lsfc_index_ = ViewVectorType<Integer>("boundary_lsfc_index_", h_ic());

        /*   parallel_for(
            "print scan", rank_size + 1, KOKKOS_LAMBDA(const int i) {
                printf(" scan boundary: %i-%li\n", i, scan_boundary_(i));
            }); */

        /* We use this strategy so that the compacted elements from the local_sfc
        would still be sorted and unique. */
        compact_boundary_elements(scan_boundary_, rank_boundary, rank_scan, rank_size);

        /* parallel_for(
            "print set", h_ic(), KOKKOS_LAMBDA(const Integer i) {

                const Integer rank = find_owner_processor(scan_boundary_, i, 1, proc);

                printf(" boundary_ : %i - %li (%li) - proc: %li - rank: %li\n", i, boundary_(i),
                            get_octant_from_sfc<Type>(boundary_(i)).template get_global_index<Type>(xDim, yDim),
                            rank , proc);
            }); */
    }



private:
    ViewMatrixTextureC<Integer, Comb::value, 2> combinations;

    ViewMatrixType<Integer> elements_;
    ViewMatrixType<Real> points_;
    ViewVectorType<bool> active_;
    Integer elements_size_;
    Integer points_size_;

    ViewVectorType<Integer> local_sfc_;
    ViewVectorType<Integer> gp_np; // parallel partition info shared among all processes.
    Integer xDim, yDim, zDim;
    Integer chunk_size_;
    Integer proc;

    UnorderedMap<Integer, Integer> global_to_local_map_;

    ViewVectorType<Integer> boundary_;
    ViewVectorType<Integer> boundary_lsfc_index_;
    ViewVectorType<Integer> scan_boundary_;
};

using DistributedMesh1 = mars::Mesh<1, 1, DistributedImplementation, Simplex<1, 1, DistributedImplementation>>;
using DistributedMesh2 = mars::Mesh<2, 2, DistributedImplementation, Simplex<2, 2, DistributedImplementation>>;
using DistributedMesh3 = mars::Mesh<3, 3, DistributedImplementation, Simplex<3, 3, DistributedImplementation>>;
using DistributedMesh4 = mars::Mesh<4, 4, DistributedImplementation, Simplex<4, 4, DistributedImplementation>>;

using DistributedQuad4Mesh = mars::Mesh<2, 2, DistributedImplementation, Quad4DElem>;
using DistributedHex8Mesh = mars::Mesh<3, 3, DistributedImplementation, Hex8DElem>;

template <Integer Type>
using DistributedNSMesh2 = mars::Mesh<2, 2, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;

template <Integer Dim, Integer ManifoldDim, Integer Type>
using DistributedNSMesh = mars::Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>>;
} // namespace mars
#endif //MARS_MESH_HPP
