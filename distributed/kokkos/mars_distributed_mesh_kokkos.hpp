#ifndef MARS_DIST_MESH_KOKKOS_HPP
#define MARS_DIST_MESH_KOKKOS_HPP

#include "mars_point.hpp"
#include "mars_visualization.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "mars_distributed_non_simplex_kokkos.hpp"
#include "mars_distributed_simplex_kokkos.hpp"
#include "mars_imesh_kokkos.hpp"

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
    ViewVectorType<unsigned int> &get_view_sfc() const
    {
        return local_sfc_;
    }

    MARS_INLINE_FUNCTION
    void set_view_sfc(ViewVectorType<unsigned int> &local) const
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
        Integer xDim;
        Integer yDim;
        Integer zDim;

        AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm) : points(pts), xDim(xdm), yDim(ydm)
        {
        }

        AddNonSimplexPoint(ViewMatrixType<Real> pts, Integer xdm, Integer ydm, Integer zdm) : points(pts), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const
        {

            int index = j * (xDim + 1) + i;

            points(index, 0) = static_cast<Real>(i) / static_cast<Real>(xDim);
            points(index, 1) = static_cast<Real>(j) / static_cast<Real>(yDim);
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int z, int y, int x) const
        {

            int index = (xDim + 1) * (yDim + 1) * z + (xDim + 1) * y + x;

            points(index, 0) = static_cast<Real>(x) / static_cast<Real>(xDim);
            points(index, 1) = static_cast<Real>(y) / static_cast<Real>(yDim);
            points(index, 2) = static_cast<Real>(z) / static_cast<Real>(zDim);
        }
    };

    inline bool generate_points(const int xDim, const int yDim, const int zDim,
                                Integer type)
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
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim == 0);

                const int n_nodes = local_sfc_.extent(0);
                reserve_points(n_nodes);

                parallel_for(
                    MDRangePolicy<Rank<2>>({0, 0}, {yDim + 1, xDim + 1}),
                    AddNonSimplexPoint(points_, xDim, yDim));
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
                assert(xDim != 0);
                assert(yDim != 0);
                assert(zDim != 0);

                const int n_nodes = (xDim + 1) * (yDim + 1) * (zDim + 1);

                reserve_points(n_nodes);

                parallel_for(
                    MDRangePolicy<Rank<3>>({0, 0, 0},
                                           {zDim + 1, yDim + 1, xDim + 1}),
                    AddNonSimplexPoint(points_, xDim, yDim, zDim));
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

        ViewMatrixType<Integer> elem;
        ViewVectorType<bool> active;

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

        AddNonSimplexElem(ViewMatrixType<Integer> el, ViewVectorType<bool> ac,
                          Integer xdm, Integer ydm, Integer zdm) : elem(el), active(ac), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const
        {
            const int offset = xDim + 1;
            int index = j * xDim + i;

            elem(index, 0) = i + offset * j;
            elem(index, 1) = (i + 1) + offset * j;
            elem(index, 2) = (i + 1) + offset * (j + 1);
            elem(index, 3) = i + offset * (j + 1);

            active(index) = true;
        }

        KOKKOS_INLINE_FUNCTION
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
        }
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
                const int n_elements = xDim * yDim;
                reserve_elements(n_elements);

                parallel_for(MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}),
                             AddNonSimplexElem(elements_, active_, xDim, yDim));
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
                const int n_elements = xDim * yDim * zDim;
                reserve_elements(n_elements);

                parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                             AddNonSimplexElem(elements_, active_, xDim, yDim, zDim));
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
