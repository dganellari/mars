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
#include <type_traits>
#include <vector>

#ifdef MARS_ENABLE_MPI
#include "mars_env.hpp"
#include "mars_distributed_non_simplex_kokkos.hpp"
#include "mars_distributed_octant_utils.hpp"
#include "mars_distributed_simplex_kokkos.hpp"
#include "mars_distributed_utils.hpp"
#include "mars_imesh_kokkos.hpp"
#include "mars_sfc_generation.hpp"
#endif
#include "mars_utils_kokkos.hpp"
#include <Kokkos_Sort.hpp>
#define ASSERT(left, operator, right)                                                                                  \
    {                                                                                                                  \
        if (!((left) operator(right))) {                                                                               \
            std::cerr << "ASSERT FAILED: " << #left                                                                    \
                      << #operator<< #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" <<(left) \
                      << "; " << #right << "=" << (right) << std::endl;                                                \
        }                                                                                                              \
    }

namespace mars {

    template <Integer Dim_, Integer ManifoldDim_, class Simplex_, class SfcKey>
    class StructuredMeshView {
    public:

        static constexpr Integer Dim = Dim_;
        static constexpr Integer ManifoldDim = ManifoldDim_;

        using Elem = typename DeviceElement<Simplex_>::type;
        using Point = mars::Point<Real, Dim>;
        using Comb = Combinations<ManifoldDim + 1, 2, DistributedImplementation>;

        using KeyType = typename SfcKey::ValueType;
        using SfcKeyType = SfcKey;

        MARS_INLINE_FUNCTION Elem elem(const Integer id) {
            return static_cast<const StructuredMeshView &>(*this).elem(id);
        }

        MARS_INLINE_FUNCTION const Elem elem(const Integer id) const {
            assert(id >= 0);
            assert(id < n_elements());
            const SubView<Integer, Elem::ElemType> nodes(elements_, id);
            if constexpr (std::is_constructible_v<Elem, decltype(nodes), decltype(combinations)>) {
                Elem element(nodes, combinations);
                element.id = id;
                return element;
            } else {
                Elem element(nodes);
                element.set_elem_id(id);
                return element;
            }
        }

        MARS_INLINE_FUNCTION bool is_active(const Integer id) const {
            assert(id >= 0);
            assert(id < n_elements());
            return active_(id);
        }

        MARS_INLINE_FUNCTION bool is_valid(const Integer id) const { return id >= 0 && id < n_elements(); }

        MARS_INLINE_FUNCTION bool is_node_valid(const Integer id) const { return id >= 0 && id < n_nodes(); }

        MARS_INLINE_FUNCTION
        const ViewMatrixTextureC<Integer, Comb::value, 2> &combs() const { return combinations; }

        MARS_INLINE_FUNCTION
        void set_combs(const ViewMatrixTextureC<Integer, Comb::value, 2> &c) { combinations = c; }

        MARS_INLINE_FUNCTION
        void set_active(const Integer id, const bool val = true) {
            assert(id >= 0);
            assert(id < n_elements());
            active_(id) = val;
        }

        MARS_INLINE_FUNCTION
        void add_point(const Point &point, const Integer index) {
            for (Integer i = 0; i < Dim; ++i) {
                points_(index, i) = point[i];
            }
        }

        MARS_INLINE_FUNCTION
        void add_point(const TempArray<Real, Dim> &point, const Integer index) {
            for (Integer i = 0; i < Dim; ++i) {
                points_(index, i) = point[i];
            }
        }

        MARS_INLINE_FUNCTION
        Point point(const Integer i) {
            assert(i >= 0);
            assert(i < points_size_);
            return Point(points_, i);
        }

        MARS_INLINE_FUNCTION
        const Point point(const Integer i) const {
            assert(i >= 0);
            assert(i < points_size_);
            return Point(points_, i);
        }

        const ViewMatrixType<Real> &points() const  //
        {
            return points_;
        }

        ViewMatrixType<Real> get_view_points() const  //
        {
            return points_;
        }

        ViewMatrixType<Integer> get_view_elements() const  //
        {
            return elements_;
        }

        MARS_INLINE_FUNCTION
        ViewVectorType<bool> get_view_active() const  //
        {
            return active_;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<KeyType> &get_view_sfc() const { return local_sfc_; }

        MARS_INLINE_FUNCTION
        void set_view_sfc(const ViewVectorType<KeyType> &local) { local_sfc_ = local; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<KeyType> &get_view_gp() const { return gp_np; }

        MARS_INLINE_FUNCTION
        void set_view_gp(const ViewVectorType<KeyType> &gp) { gp_np = gp; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<KeyType> &get_view_boundary() const { return boundary_; }

        MARS_INLINE_FUNCTION
        KeyType get_boundary_sfc(const KeyType sfc_index) const { return boundary_(sfc_index); }

        MARS_INLINE_FUNCTION
        void set_view_boundary(const ViewVectorType<KeyType> &b) { boundary_ = b; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<KeyType> &get_view_boundary_sfc_index() const { return boundary_lsfc_index_; }

        MARS_INLINE_FUNCTION
        void set_view_boundary_sfc_index(const ViewVectorType<KeyType> &b) { boundary_lsfc_index_ = b; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_boundary() const { return scan_boundary_; }

        MARS_INLINE_FUNCTION
        void set_view_scan_boundary(const ViewVectorType<Integer> &b) { scan_boundary_ = b; }

        MARS_INLINE_FUNCTION
        void set_XDim(const Integer xd) { xDim = xd; }

        MARS_INLINE_FUNCTION
        Integer get_XDim() const { return xDim; }

        MARS_INLINE_FUNCTION
        void set_YDim(const Integer yd) { yDim = yd; }

        MARS_INLINE_FUNCTION
        Integer get_YDim() const { return yDim; }

        MARS_INLINE_FUNCTION
        void set_ZDim(const Integer zd) { zDim = zd; }

        MARS_INLINE_FUNCTION
        Integer get_ZDim() const { return zDim; }

        MARS_INLINE_FUNCTION
        Integer n_nodes() const { return points_size_; }

        MARS_INLINE_FUNCTION
        void set_n_nodes(Integer size) { points_size_ = size; }

        MARS_INLINE_FUNCTION
        Integer n_elements() const { return elements_size_; }

        MARS_INLINE_FUNCTION
        void set_n_elements(Integer size) { elements_size_ = size; }

        MARS_INLINE_FUNCTION
        void set_chunk_size(Integer size) { chunk_size_ = size; }

        MARS_INLINE_FUNCTION
        Integer get_ghost_size() const { return get_view_ghost().extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_boundary_size() const { return get_view_boundary().extent(0); }

        MARS_INLINE_FUNCTION
        Integer get_chunk_size() const { return chunk_size_; }

        MARS_INLINE_FUNCTION
        void set_proc(Integer p) { proc = p; }

        MARS_INLINE_FUNCTION
        Integer get_proc() const { return proc; }

        MARS_INLINE_FUNCTION Integer type() const { return ManifoldDim; }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) {
            const Integer size = get_chunk_size();
            Kokkos::parallel_for("init_initial_cond", size, f);
        }

        template <typename H>
        struct FaceIterate {
            using simplex_type = typename StructuredMeshView::Elem;
            FaceIterate(StructuredMeshView m,
                        H f,
                        ViewVectorType<KeyType> gl,
                        ViewVectorType<Integer> sg,
                        Integer p,
                        Integer x,
                        Integer y,
                        Integer z)
                : mesh(m), func(f), ghost_layer(gl), scan_ghost(sg), proc(p), xDim(x), yDim(y), zDim(z) {}

            template <Integer dir>
            MARS_INLINE_FUNCTION void iterate(const Integer i) const {
                // side  0 means origin side and 1 destination side.
                for (int side = 0; side < 2; ++side) {
                    Integer face_nr;

                    if (side == 0)
                        face_nr = 2 * dir + 1;
                    else
                        face_nr = 2 * dir;

                    Octant nbh_oc = mesh.get_octant_face_nbh(i, face_nr);

                    bool ghost = false;
                    Integer index;

                    if (nbh_oc.is_valid()) {
                        auto enc_oc = get_sfc_from_octant<simplex_type::ElemType, SfcKeyType>(nbh_oc);

                        Integer owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, proc);
                        assert(owner_proc >= 0);

                        if (proc != owner_proc) {
                            ghost = true;

                            const int start_index = scan_ghost(owner_proc);
                            const int last_index = scan_ghost(owner_proc + 1) - 1;

                            index = binary_search(ghost_layer.data(), start_index, last_index, enc_oc);
                            assert(index >= 0);
                        } else {
                            // using the sfc (global) to local mapping of the mesh.
                            index = mesh.get_index_of_sfc_elem(enc_oc);
                            assert(index >= 0);
                        }
                    }

                    bool boundary = nbh_oc.shares_boundary();

                    Face<simplex_type::ElemType, dir> face;

                    if (side == 1 && mesh.is_periodic() && boundary && !ghost) {
                        nbh_oc.set_invalid();
                        face.invalidate();
                    }

                    if (face.is_valid() && ((side == 0 && nbh_oc.is_valid()) || ghost || boundary)) {
                        int origin_side = side;

                        if (boundary && !mesh.is_periodic()) origin_side = 0;

                        face.get_side(origin_side).set_elem_id(i);
                        face.get_side(origin_side).set_boundary(boundary);

                        face.get_side(origin_side).set_origin();

                        if (!boundary || mesh.is_periodic()) {
                            int otherside = origin_side ^ 1;

                            face.get_side(otherside).set_elem_id(index);
                            face.get_side(otherside).set_ghost(ghost);
                        }

                        if (boundary && mesh.is_periodic()) face.swap_sides();

                        func(face);
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                iterate<0>(i);
                iterate<1>(i);
                // TODO: 3D part
            }

            StructuredMeshView mesh;
            H func;

            ViewVectorType<KeyType> ghost_layer;
            ViewVectorType<Integer> scan_ghost;

            Integer proc;
            Integer xDim;
            Integer yDim;
            Integer zDim;
        };

        template <typename H>
        void face_iterate(H f) {
            Integer xDim = get_XDim();
            Integer yDim = get_YDim();
            Integer zDim = get_ZDim();

            Kokkos::parallel_for(
                "face_iterate",
                get_chunk_size(),
                FaceIterate<H>(*this, f, get_view_ghost(), get_view_scan_ghost(), get_proc(), xDim, yDim, zDim));
        }

        void print_sfc() {
            using namespace Kokkos;

            Integer xDim = get_XDim();
            Integer yDim = get_YDim();
            Integer zDim = get_ZDim();

            auto sfcv = get_view_sfc();
            auto proc_num = get_proc();
            parallel_for(
                "print set", get_chunk_size(), KOKKOS_LAMBDA(const Integer i) {
                    const Integer sfc = sfcv(i);
                    double point[3];
                    get_vertex_coordinates_from_sfc<Elem::ElemType, SfcKeyType>(sfc, point, xDim, yDim, zDim);

                    Octant o = get_octant_from_sfc<Elem::ElemType, SfcKeyType>(sfc);
                    printf("mesh sfc : %li - %li - %li - (%lf, %lf, %lf) -rank: %i\n",
                           i,
                           sfc,
                           elem_index(o.x, o.y, o.z, xDim, yDim),
                           point[0],
                           point[1],
                           point[2],
                           proc_num);
                });
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc_index(const Integer enc_oc) {
            return binary_search(get_view_sfc().data(), 0, get_chunk_size() - 1, enc_oc);
        }

        MARS_INLINE_FUNCTION
        Integer get_index_of_sfc_elem(const Integer enc_oc) const {

            auto index = INVALID_INDEX;
            const auto it = get_sfc_to_local_map().find(enc_oc);
            if (get_sfc_to_local_map().valid_at(it)) {
                index = get_sfc_to_local_map().value_at(it);
            }
            assert(get_sfc_to_local_map().valid_at(it));
            return index;
        }

        MARS_INLINE_FUNCTION
        void set_periodic() { periodic = true; }

        MARS_INLINE_FUNCTION
        bool is_periodic() const { return periodic; }

        MARS_INLINE_FUNCTION
        KeyType get_ghost_sfc(const Integer index) const { return ghost_(index); }

        MARS_INLINE_FUNCTION
        Octant get_ghost_octant(const Integer index) const {
            const KeyType sfc_code = ghost_(index);
            return get_octant_from_sfc<Elem::ElemType, SfcKeyType>(sfc_code);
        }

        MARS_INLINE_FUNCTION
        KeyType get_sfc(const Integer sfc_index) const { return local_sfc_(sfc_index); }

        MARS_INLINE_FUNCTION
        Octant get_octant(const Integer sfc_index) const {
            const KeyType sfc_code = local_sfc_(sfc_index);
            return get_octant_from_sfc<Elem::ElemType, SfcKeyType>(sfc_code);
        }

        MARS_INLINE_FUNCTION
        Octant octant_from_sfc(const KeyType sfc) const { return get_octant_from_sfc<Elem::ElemType, SfcKeyType>(sfc); }

        MARS_INLINE_FUNCTION
        Octant get_octant_face_nbh(const Octant &oc, const Integer face_nr) const {
            return oc.face_nbh<Elem::ElemType>(face_nr, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_face_nbh(const Integer sfc_index, const Integer face_nr) const {
            Octant oc = get_octant(sfc_index);

            return oc.face_nbh<Elem::ElemType>(face_nr, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_edge_start(const Octant &oc, const Integer edge) const {
            return oc.edge_start<Elem::ElemType>(edge, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_edge_start(const Integer sfc_index, const Integer edge) const {
            Octant oc = get_octant(sfc_index);

            return oc.edge_start<Elem::ElemType>(edge, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_edge_nbh(const Octant &oc, const Integer edge) const {
            return oc.edge_nbh<Elem::ElemType>(edge, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_edge_nbh(const Integer sfc_index, const Integer edge) const {
            Octant oc = get_octant(sfc_index);

            return oc.edge_nbh<Elem::ElemType>(edge, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_corner_nbh(const Octant &oc, const Integer corner_nr) {
            return oc.corner_nbh<Elem::ElemType>(corner_nr, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        Octant get_octant_corner_nbh(const Integer sfc_index, const Integer corner_nr) {
            Octant oc = get_octant(sfc_index);

            return oc.corner_nbh<Elem::ElemType>(corner_nr, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        template <typename F>
        MARS_INLINE_FUNCTION void get_one_ring_edge_nbhs(const Octant &oc, const Integer direction, F f) const {
            oc.one_ring_edge_nbhs<Elem::ElemType, F>(f, direction, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        template <typename F>
        MARS_INLINE_FUNCTION void get_one_ring_corner_nbhs(const Octant &oc, F f) const {
            oc.one_ring_corner_nbhs<Elem::ElemType, F>(f, get_XDim(), get_YDim(), get_ZDim(), is_periodic());
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_ghost() const { return scan_ghost_; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<KeyType> &get_view_ghost() const { return ghost_; }

        MARS_INLINE_FUNCTION
        bool is_owned(const Integer sfc) const {
            const auto it = get_sfc_to_local_map().find(sfc);
            return get_sfc_to_local_map().valid_at(it);
        }

        MARS_INLINE_FUNCTION
        bool is_owned_index(const Integer sfc_index) const {
            auto sfc = get_sfc(sfc_index);
            return is_owned(sfc);
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<KeyType, Integer> &get_sfc_to_local_map() const { return sfc_to_local_map_; }

        MARS_INLINE_FUNCTION
        void set_sfc_to_local_map(const UnorderedMap<KeyType, Integer> &map) { sfc_to_local_map_ = map; }
        MARS_INLINE_FUNCTION StructuredMeshView device_view() const { return *this; }

    protected:

        // This block represents the SFC data and the data structures needed to manage it in a distributed manner.
        ViewVectorType<KeyType> local_sfc_;

        UnorderedMap<KeyType, Integer> sfc_to_local_map_;
        ViewVectorType<KeyType> gp_np;  // parallel partition info shared among all processes.
        Integer xDim = 0, yDim = 0, zDim = 0;
        Integer chunk_size_ = 0;
        Integer proc = 0;

        // Periodic mesh feature supported.
        bool periodic = false;

        ViewMatrixTextureC<Integer, Comb::value, 2> combinations;

        ViewMatrixType<Integer> elements_;
        ViewMatrixType<Real> points_;
        ViewVectorType<bool> active_;
        Integer elements_size_ = 0;
        Integer points_size_ = 0;

        // Boundary and ghost layer data
        ViewVectorType<KeyType> boundary_;             // sfc code for the ghost layer
        ViewVectorType<KeyType> boundary_lsfc_index_;  // view index of the previous
        ViewVectorType<Integer> scan_boundary_;
        // ghost data
        ViewVectorType<KeyType> ghost_;
        ViewVectorType<Integer> scan_ghost_;
    };

    template <Integer Dim_, Integer ManifoldDim_, class Simplex_, class SfcKey>
    class Mesh<Dim_, ManifoldDim_, DistributedImplementation, Simplex_, SfcKey>
        : public StructuredMeshView<Dim_, ManifoldDim_, Simplex_, SfcKey>, public ParallelIMesh<Dim_> {
    public:
        using DeviceView = StructuredMeshView<Dim_, ManifoldDim_, Simplex_, SfcKey>;

        static constexpr Integer Dim = Dim_;
        static constexpr Integer ManifoldDim = ManifoldDim_;

        using Elem = Simplex_;
        using Point = mars::Point<Real, Dim>;
        using Comb = Combinations<ManifoldDim + 1, 2, DistributedImplementation>;

        using KeyType = typename SfcKey::ValueType;
        using SfcKeyType = SfcKey;
        using DeviceView::add_point;
        using DeviceView::combs;
        using DeviceView::elem;
        using DeviceView::elem_iterate;
        using DeviceView::face_iterate;
        using DeviceView::get_XDim;
        using DeviceView::get_YDim;
        using DeviceView::get_ZDim;
        using DeviceView::get_boundary_sfc;
        using DeviceView::get_boundary_size;
        using DeviceView::get_chunk_size;
        using DeviceView::get_ghost_octant;
        using DeviceView::get_ghost_sfc;
        using DeviceView::get_ghost_size;
        using DeviceView::get_index_of_sfc_elem;
        using DeviceView::get_octant;
        using DeviceView::get_octant_corner_nbh;
        using DeviceView::get_octant_edge_nbh;
        using DeviceView::get_octant_edge_start;
        using DeviceView::get_octant_face_nbh;
        using DeviceView::get_one_ring_corner_nbhs;
        using DeviceView::get_one_ring_edge_nbhs;
        using DeviceView::get_proc;
        using DeviceView::get_sfc;
        using DeviceView::get_sfc_index;
        using DeviceView::get_sfc_to_local_map;
        using DeviceView::get_view_active;
        using DeviceView::get_view_boundary;
        using DeviceView::get_view_boundary_sfc_index;
        using DeviceView::get_view_elements;
        using DeviceView::get_view_ghost;
        using DeviceView::get_view_gp;
        using DeviceView::get_view_points;
        using DeviceView::get_view_scan_boundary;
        using DeviceView::get_view_scan_ghost;
        using DeviceView::get_view_sfc;
        using DeviceView::is_node_valid;
        using DeviceView::is_owned;
        using DeviceView::is_owned_index;
        using DeviceView::is_periodic;
        using DeviceView::is_valid;
        using DeviceView::octant_from_sfc;
        using DeviceView::points;
        using DeviceView::print_sfc;
        using DeviceView::set_XDim;
        using DeviceView::set_YDim;
        using DeviceView::set_ZDim;
        using DeviceView::set_active;
        using DeviceView::set_chunk_size;
        using DeviceView::set_combs;
        MARS_INLINE_FUNCTION
        void set_global_to_local_map(const UnorderedMap<Integer, Integer> &gl_map) { global_to_local_map_ = gl_map; }
        using DeviceView::set_n_elements;
        using DeviceView::set_n_nodes;
        using DeviceView::set_periodic;
        using DeviceView::set_proc;
        using DeviceView::set_sfc_to_local_map;
        using DeviceView::set_view_boundary;
        using DeviceView::set_view_boundary_sfc_index;
        using DeviceView::set_view_gp;
        using DeviceView::set_view_scan_boundary;
        using DeviceView::set_view_sfc;
        explicit Mesh(const context& c) : ctx(c) {}
        Mesh(context&&) = delete;
        Mesh(const context&&) = delete;
        DeviceView device_view() const { return static_cast<const DeviceView&>(*this); }
        MARS_INLINE_FUNCTION Integer n_nodes() const override { return DeviceView::n_nodes(); }
        MARS_INLINE_FUNCTION Integer n_elements() const override { return DeviceView::n_elements(); }
        MARS_INLINE_FUNCTION Integer type() const override { return DeviceView::type(); }
        MARS_INLINE_FUNCTION bool is_active(Integer i) const override { return DeviceView::is_active(i); }
        MARS_INLINE_FUNCTION Point point(Integer i) override { return DeviceView::point(i); }
        MARS_INLINE_FUNCTION const Point point(Integer i) const override { return DeviceView::point(i); }
        void set_partition_local_counts(std::size_t generated, std::size_t received) {
            partition_generated_ = generated;
            partition_received_ = received;
        }
        std::size_t partition_generated_keys() const { return partition_generated_; }
        std::size_t partition_received_keys() const { return partition_received_; }

        void reserve(const std::size_t n_elements, const std::size_t n_points) override {
            elements_size_ = n_elements;
            points_size_ = n_points;

            elements_ = ViewMatrixType<Integer>("elems", n_elements, Elem::ElemType);
            active_ = ViewVectorType<bool>("active_", n_elements);
            points_ = ViewMatrixType<Real>("pts", n_points, Dim);
        }

        void reserve_points(const std::size_t n_points) {
            points_size_ = n_points;
            points_ = ViewMatrixType<Real>("pts", n_points, Dim);
        }

        void reserve_elements(const std::size_t n_elements) {
            elements_size_ = n_elements;
            elements_ = ViewMatrixType<Integer>("elems", n_elements, Elem::ElemType);
            active_ = ViewVectorType<bool>("active_", n_elements);
        }

        void resize_points(const Integer size) {
            points_size_ += size;
            resize(points_, points_size_, Dim);
        }

        void resize_elements(const Integer size) {
            elements_size_ += size;
            resize(elements_, elements_size_, Elem::ElemType);
            resize(active_, elements_size_);
        }

        void resize_active(const Integer size) { resize(active_, elements_size_ + size); }

        inline Integer n_active_elements(const Integer N) {
            using namespace Kokkos;
            ViewVectorType<bool> active_ = get_view_active();

            double result = 0;
            parallel_reduce(
                "Active_all", N, KOKKOS_LAMBDA(const int &i, double &lsum) { lsum += active_(i); }, result);

            return result;
        }

        inline Integer n_active_elements(const ViewVectorType<Integer> elements) {
            using namespace Kokkos;

            ViewVectorType<bool> active_ = get_view_active();

            Integer N = elements.extent(0);

            double result = 0;
            parallel_reduce(
                "Active_elem", N, KOKKOS_LAMBDA(const int &i, double &lsum) { lsum += active_(elements(i)); }, result);

            return result;
        }

        // add point functor
        template <Integer Type>
        struct AddNonSimplexPoint {
            ViewMatrixType<Real> points;
            ViewVectorType<KeyType> encoded_points;

            Integer xDim;
            Integer yDim;
            Integer zDim;

            AddNonSimplexPoint(ViewMatrixType<Real> pts, ViewVectorType<KeyType> epts, Integer xdm, Integer ydm)
                : points(pts), encoded_points(epts), xDim(xdm), yDim(ydm) {}

            AddNonSimplexPoint(ViewMatrixType<Real> pts,
                               ViewVectorType<KeyType> epts,
                               Integer xdm,
                               Integer ydm,
                               Integer zdm)
                : points(pts), encoded_points(epts), xDim(xdm), yDim(ydm), zDim(zdm) {}

            KOKKOS_INLINE_FUNCTION
            void operator()(Integer index) const {
                KeyType gl_index = encoded_points(index);

                switch (Type) {
                    case ElementType::Quad4: {
                        auto octant = decode_sfc_2D<SfcKeyType>(gl_index);

                        points(index, 0) = static_cast<Real>(octant.x) / static_cast<Real>(xDim);
                        points(index, 1) = static_cast<Real>(octant.y) / static_cast<Real>(yDim);
                        break;
                    }
                    case ElementType::Hex8: {
                        auto octant = decode_sfc_3D<SfcKeyType>(gl_index);
                        points(index, 0) = static_cast<Real>(octant.x) / static_cast<Real>(xDim);
                        points(index, 1) = static_cast<Real>(octant.y) / static_cast<Real>(yDim);
                        points(index, 2) = static_cast<Real>(octant.z) / static_cast<Real>(zDim);
                        break;
                    }
                }
            }
        };

        template <Integer Type>
        ViewVectorType<KeyType> build_local_point_keys() const {
            static_assert(Type == ElementType::Quad4 || Type == ElementType::Hex8);
            constexpr Integer corners = Type == ElementType::Quad4 ? 4 : 8;
            using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
            const auto cell_keys = get_view_sfc();
            const Integer count = corners * get_chunk_size();
            ViewVectorType<KeyType> corner_keys("local_corner_keys", count);

            // A Hilbert corner is not a key-space bound. Enumerate only local cell corners,
            // then form their union without concurrent writes to a shared predicate.
            Kokkos::parallel_for("local_corner_keys", Policy(0, count), KOKKOS_LAMBDA(const Integer i) {
                const Integer corner = i % corners;
                if constexpr (Type == ElementType::Quad4) {
                    const auto cell = decode_sfc_2D<SfcKeyType>(cell_keys(i / corners));
                    corner_keys(i) = encode_sfc_2D<SfcKeyType>(cell.x + (corner & 1),
                                                              cell.y + ((corner >> 1) & 1));
                } else {
                    const auto cell = decode_sfc_3D<SfcKeyType>(cell_keys(i / corners));
                    corner_keys(i) = encode_sfc_3D<SfcKeyType>(cell.x + (corner & 1),
                                                              cell.y + ((corner >> 1) & 1),
                                                              cell.z + ((corner >> 2) & 1));
                }
            });
            Kokkos::sort(corner_keys);
            Integer unique_count = 0;
            Kokkos::parallel_reduce("count_local_points", Policy(0, count),
                KOKKOS_LAMBDA(const Integer i, Integer& total) {
                    if (i == 0 || corner_keys(i) != corner_keys(i - 1)) ++total;
                }, unique_count);
            ViewVectorType<KeyType> point_keys("local_point_keys", unique_count);
            Kokkos::parallel_scan("compact_local_points", Policy(0, count),
                KOKKOS_LAMBDA(const Integer i, Integer& offset, const bool final) {
                    if (i == 0 || corner_keys(i) != corner_keys(i - 1)) {
                        if (final) point_keys(offset) = corner_keys(i);
                        ++offset;
                    }
                });
            return point_keys;
        }

        template <Integer Type>
        inline bool generate_points() {
            using namespace Kokkos;

            switch (ManifoldDim_) {
                case 2: {
                    switch (Type) {
                        case ElementType::Quad4: {
                            assert(xDim != 0);
                            assert(yDim != 0);
                            assert(zDim == 0);

                            const auto local_to_global = build_local_point_keys<Type>();
                            const Integer nr_points = local_to_global.extent(0);

                            reserve_points(nr_points);

                            parallel_for("non_simplex_point",
                                         nr_points,
                                         AddNonSimplexPoint<Type>(points_, local_to_global, xDim, yDim));

                            // build global_to_local map for use in generate elements.
                            UnorderedMap<Integer, Integer> global_to_local_map(nr_points);

                            // otherwise lambda does not see the xDim
                            const Integer xD = xDim;

                            parallel_for(
                                "local_global_map_for", nr_points, KOKKOS_LAMBDA(const int i) {
                                    const int offset = xD + 1;

                                    const auto gl_index = local_to_global(i);

                                    const auto octant = decode_sfc_2D<SfcKeyType>(gl_index);

                                    global_to_local_map.insert(octant.x + offset * octant.y, i);
                                });

                            set_global_to_local_map(global_to_local_map);

                            return true;
                        }
                        default: {
                            return false;
                        }
                    }
                }
                case 3: {
                    switch (Type) {
                        case ElementType::Hex8: {
                            assert(xDim != 0);
                            assert(yDim != 0);
                            assert(zDim != 0);

                            const auto local_to_global = build_local_point_keys<Type>();
                            const Integer nr_points = local_to_global.extent(0);
                            reserve_points(nr_points);
                            parallel_for("non_simplex_point",
                                         nr_points,
                                         AddNonSimplexPoint<Type>(points_, local_to_global, xDim, yDim, zDim));
                            // build global_to_local map for use in generate elements.
                            UnorderedMap<Integer, Integer> global_to_local_map(nr_points);

                            // otherwise lambda does not see the xDim or ydim
                            const Integer xD = xDim;
                            const Integer yD = yDim;

                            parallel_for(
                                "local_global_map_for", nr_points, KOKKOS_LAMBDA(const int i) {
                                    const auto gl_index = local_to_global(i);

                                    const auto octant = decode_sfc_3D<SfcKeyType>(gl_index);

                                    global_to_local_map.insert(elem_index(octant.x, octant.y, octant.z, xD, yD), i);
                                });

                            set_global_to_local_map(global_to_local_map);
                            return true;
                        }
                        default: {
                            return false;
                        }
                    }
                }
                default: {
                    return false;
                }
            }
        }

        // add elem functor
        template <Integer Type>
        struct AddNonSimplexElem {
            using UMap = UnorderedMap<Integer, Integer>;
            ViewMatrixType<Integer> elem;
            ViewVectorType<bool> active;
            ViewVectorType<KeyType> global;
            UMap map;

            Integer xDim;
            Integer yDim;
            Integer zDim;

            AddNonSimplexElem(ViewMatrixType<Integer> el,
                              ViewVectorType<KeyType> gl,
                              UMap mp,
                              ViewVectorType<bool> ac,
                              Integer xdm,
                              Integer ydm)
                : elem(el), global(gl), map(mp), active(ac), xDim(xdm), yDim(ydm) {}

            AddNonSimplexElem(ViewMatrixType<Integer> el,
                              ViewVectorType<KeyType> gl,
                              UMap mp,
                              ViewVectorType<bool> ac,
                              Integer xdm,
                              Integer ydm,
                              Integer zdm)
                : elem(el), global(gl), map(mp), active(ac), xDim(xdm), yDim(ydm), zDim(zdm) {}

            KOKKOS_INLINE_FUNCTION
            void operator()(int index) const {
                switch (Type) {
                    case ElementType::Quad4: {
                        const int offset = xDim + 1;

                        const auto gl_index = global(index);
                        const auto octant = decode_sfc_2D<SfcKeyType>(gl_index);

                        elem(index, 0) = map.value_at(map.find(octant.x + offset * octant.y));
                        elem(index, 1) = map.value_at(map.find((octant.x + 1) + offset * octant.y));
                        elem(index, 2) = map.value_at(map.find((octant.x + 1) + offset * (octant.y + 1)));
                        elem(index, 3) = map.value_at(map.find(octant.x + offset * (octant.y + 1)));

                        active(index) = true;
                        break;
                    }
                    case ElementType::Hex8: {
                        const auto gl_index = global(index);

                        const auto octant = decode_sfc_3D<SfcKeyType>(gl_index);

                        elem(index, 0) = map.value_at(map.find(elem_index(octant.x, octant.y, octant.z, xDim, yDim)));
                        elem(index, 1) =
                            map.value_at(map.find(elem_index(octant.x + 1, octant.y, octant.z, xDim, yDim)));
                        elem(index, 2) =
                            map.value_at(map.find(elem_index(octant.x + 1, octant.y + 1, octant.z, xDim, yDim)));
                        elem(index, 3) =
                            map.value_at(map.find(elem_index(octant.x, octant.y + 1, octant.z, xDim, yDim)));
                        elem(index, 4) =
                            map.value_at(map.find(elem_index(octant.x, octant.y, octant.z + 1, xDim, yDim)));
                        elem(index, 5) =
                            map.value_at(map.find(elem_index(octant.x + 1, octant.y, octant.z + 1, xDim, yDim)));
                        elem(index, 6) =
                            map.value_at(map.find(elem_index(octant.x + 1, octant.y + 1, octant.z + 1, xDim, yDim)));
                        elem(index, 7) =
                            map.value_at(map.find(elem_index(octant.x, octant.y + 1, octant.z + 1, xDim, yDim)));

                        active(index) = true;
                        break;
                    }
                }
            }
        };

        template <Integer Type>
        inline bool generate_elements() {
            using namespace Kokkos;

            switch (ManifoldDim_) {
                case 2: {
                    switch (Type) {
                        case ElementType::Quad4: {
                            reserve_elements(get_chunk_size());

                            parallel_for("generate_elements",
                                         get_chunk_size(),
                                         AddNonSimplexElem<Type>(
                                             elements_, local_sfc_, global_to_local_map_, active_, xDim, yDim));
                            return true;
                        }
                        default: {
                            return false;
                        }
                    }
                }
                case 3: {
                    switch (Type) {
                        case ElementType::Hex8: {
                            reserve_elements(get_chunk_size());

                            parallel_for("generate_elements",
                                         get_chunk_size(),
                                         AddNonSimplexElem<Type>(
                                             elements_, local_sfc_, global_to_local_map_, active_, xDim, yDim, zDim));
                            return true;
                        }
                        default: {
                            return false;
                        }
                    }
                }
                default: {
                    return false;
                }
            }
        }

        template <Integer Type>
        struct CountGhostNeighbors {
            ViewVectorType<Integer> global;
            ViewVectorType<KeyType> count;
            ViewVectorType<Integer> gp;

            Integer proc;

            Integer xDim;
            Integer yDim;
            Integer zDim;

            bool periodic;

            CountGhostNeighbors(ViewVectorType<KeyType> gl,
                                ViewVectorType<Integer> ct,
                                ViewVectorType<Integer> g,
                                Integer p,
                                Integer xdm,
                                Integer ydm,
                                Integer zdm,
                                bool period)
                : global(gl), count(ct), gp(g), proc(p), xDim(xdm), yDim(ydm), zDim(zdm), periodic(period) {}

            MARS_INLINE_FUNCTION
            void increment(const Octant &o) const {
                auto enc_oc = get_sfc_from_octant<Type, SfcKeyType>(o);

                assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
                Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

                // the case when the neihgbor is a ghost element.
                if (proc != owner_proc) {
                    Kokkos::atomic_increment(&count(owner_proc));
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(int index) const {
                const auto gl_index = global(index);
                Octant ref_octant = get_octant_from_sfc<Type, SfcKeyType>(gl_index);

                const int offset = xDim + 1;

                for (int face = 0; face < 2 * ManifoldDim; ++face) {
                    Octant o = ref_octant.face_nbh<Type>(face, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        increment(o);
                    }
                }

                for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner) {
                    Octant o = ref_octant.corner_nbh<Type>(corner, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        increment(o);
                    }
                }
            }
        };

        template <Integer Type>
        struct BuildBoundarySets {
            ViewVectorType<KeyType> global;
            ViewVectorType<Integer> gp;

            ViewVectorType<Integer> set;
            ViewVectorType<Integer> scan;
            // ViewVectorType<Integer> index;

            Integer proc;

            Integer xDim;
            Integer yDim;
            Integer zDim;

            bool periodic;

            BuildBoundarySets(ViewVectorType<KeyType> gl,
                              ViewVectorType<Integer> g,
                              ViewVectorType<Integer> st,
                              ViewVectorType<Integer> sc,  // ViewVectorType<Integer> in,
                              Integer p,
                              Integer xdm,
                              Integer ydm,
                              Integer zdm,
                              bool period)
                : global(gl),
                  gp(g),
                  set(st),
                  scan(sc),  // index(in),
                  proc(p),
                  xDim(xdm),
                  yDim(ydm),
                  zDim(zdm),
                  periodic(period) {}

            MARS_INLINE_FUNCTION
            void scatter(const Octant &o, const Integer ref) const {
                auto enc_oc = get_sfc_from_octant<Type, SfcKeyType>(o);

                assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
                Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

                // the case when the neihgbor is a ghost element.
                if (proc != owner_proc) {
                    // get the starting index before incrementing it.
                    Integer i = Kokkos::atomic_fetch_add(&scan(owner_proc), 1);
                    set(i) = ref;
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(int index) const {
                const auto gl_index = global(index);
                Octant ref_octant = get_octant_from_sfc<Type, SfcKeyType>(gl_index);

                const int offset = xDim + 1;
                for (int face = 0; face < 2 * ManifoldDim; ++face) {
                    Octant o = ref_octant.face_nbh<Type>(face, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        scatter(o, gl_index);
                    }
                }

                for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner) {
                    Octant o = ref_octant.corner_nbh<Type>(corner, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        scatter(o, gl_index);
                    }
                }
            }
        };

        // another possible way of building the ghost layer directly without using the boundary layers and then send
        // with mpi. This would be an mpi less method howerver it requires extra compute time for sorting. The other
        // method is currently used.
        template <Integer Type>
        inline void build_ghost_element_sets() {
            using namespace Kokkos;

            const Integer rank_size = gp_np.extent(0) / 2 - 1;

            ViewVectorType<Integer> count("count_per_proc", rank_size);

            parallel_for("count_ghost_nbhs",
                         get_chunk_size(),
                         CountGhostNeighbors<Type>(get_view_sfc(), count, gp_np, proc, xDim, yDim, zDim));

            ViewVectorType<Integer> scan("scan", rank_size + 1);
            incl_excl_scan(0, rank_size, count, scan);

            auto index_subview = subview(scan, rank_size);
            auto h_ic = create_mirror_view(index_subview);

            // Deep copy device view to host view.
            deep_copy(h_ic, index_subview);

            // the set containing all the ghost elements for all processes. Scan helps to identify part of the set array
            // to which process it belongs.
            ViewVectorType<Integer> set("build_set", h_ic());

            parallel_for("build_set_kernel",
                         get_chunk_size(),
                         BuildBoundarySets<Type>(get_view_sfc(), gp_np, set, scan, proc, xDim, yDim, zDim));
        }

        template <Integer Type, bool OP = 0>
        struct CountOrInsertBoundaryPerRank {
            ViewVectorType<KeyType> global;
            ViewVectorType<KeyType> data;
            ViewVectorType<Integer> scan;
            ViewVectorType<KeyType> gp;

            Integer proc;

            Integer xDim;
            Integer yDim;
            Integer zDim;

            bool periodic;

            CountOrInsertBoundaryPerRank(ViewVectorType<KeyType> gl,
                                    ViewVectorType<KeyType> d,
                                    ViewVectorType<KeyType> g,
                                    Integer p,
                                    Integer xdm,
                                    Integer ydm,
                                    Integer zdm,
                                    bool period)
                : global(gl), data(d), gp(g), proc(p), xDim(xdm), yDim(ydm), zDim(zdm), periodic(period) {}

            CountOrInsertBoundaryPerRank(ViewVectorType<KeyType> gl,
                                         ViewVectorType<KeyType> d,
                                         ViewVectorType<Integer> s,
                                         ViewVectorType<KeyType> g,
                                         Integer p,
                                         Integer xdm,
                                         Integer ydm,
                                         Integer zdm,
                                         bool period)
                : global(gl), data(d), scan(s), gp(g), proc(p), xDim(xdm), yDim(ydm), zDim(zdm), periodic(period) {}

            MARS_INLINE_FUNCTION
            void increment_count(const int index, const Octant &o) const {
                auto enc_oc = get_sfc_from_octant<Type, SfcKeyType>(o);

                assert(find_owner_processor(gp, enc_oc, 2, proc) >= 0);
                Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);

                // the case when the neihgbor is a ghost element.
                if (owner_proc >= 0 && proc != owner_proc) {
                    if constexpr (OP) {
                        // get the starting index before incrementing it.
                        Integer i = Kokkos::atomic_fetch_add(&scan(owner_proc), 1);
                        data(i) = index;
                    } else {
                        Kokkos::atomic_increment(&data(owner_proc));
                    }
                }
            }

            MARS_INLINE_FUNCTION void predicate_face(const Octant ref_octant,
                                                     const Integer index,
                                                     const Integer xDim,
                                                     const Integer yDim,
                                                     const Integer zDim,
                                                     bool periodic) const {
                for (int face = 0; face < 2 * ManifoldDim; ++face) {
                    Octant o = ref_octant.face_nbh<Type>(face, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        increment_count(index, o);
                    }
                }
            }

            MARS_INLINE_FUNCTION void predicate_corner(const Octant ref_octant,
                                                       const Integer index,
                                                       const Integer xDim,
                                                       const Integer yDim,
                                                       const Integer zDim,
                                                       bool periodic) const {
                for (int corner = 0; corner < power_of_2(ManifoldDim); ++corner) {
                    Octant o = ref_octant.corner_nbh<Type>(corner, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        increment_count(index, o);
                    }
                }
            }

            template <Integer T = Type>
            MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void> predicate_edge(const Octant ref_octant,
                                                                                               const Integer index,
                                                                                               const Integer xDim,
                                                                                               const Integer yDim,
                                                                                               const Integer zDim,
                                                                                               bool periodic) const {
                for (int edge = 0; edge < 4 * ManifoldDim; ++edge) {
                    Octant o = ref_octant.edge_nbh<Type>(edge, xDim, yDim, zDim, periodic);
                    if (o.is_valid()) {
                        increment_count(index, o);
                    }
                }
            }

            template <Integer T = Type>
            MARS_INLINE_FUNCTION std::enable_if_t<T != ElementType::Hex8, void> predicate_edge(const Octant ref_octant,
                                                                                               const Integer index,
                                                                                               const Integer xDim,
                                                                                               const Integer yDim,
                                                                                               const Integer zDim,
                                                                                               bool periodic) const {}

            MARS_INLINE_FUNCTION
            void operator()(int index) const {
                const auto gl_index = global(index);
                Octant ref_octant = get_octant_from_sfc<Type, SfcKeyType>(gl_index);

                predicate_face(ref_octant, index, xDim, yDim, zDim, periodic);
                predicate_corner(ref_octant, index, xDim, yDim, zDim, periodic);
                predicate_edge(ref_octant, index, xDim, yDim, zDim, periodic);
            }
        };

        struct BoundaryCell {
            Integer peer;
            KeyType local;
        };

        struct BoundaryCellLess {
            MARS_INLINE_FUNCTION bool operator()(const BoundaryCell& a, const BoundaryCell& b) const {
                return a.peer < b.peer || (a.peer == b.peer && a.local < b.local);
            }
        };

        template <Integer Type, Integer Device = 0>
        void build_boundary_element_sets() {
            const Integer rank_size = gp_np.extent(0) / 2 - 1;
            ViewVectorType<KeyType> counts("boundary_candidates_per_peer", rank_size);
            Kokkos::parallel_for("count_boundary_candidates", get_chunk_size(),
                CountOrInsertBoundaryPerRank<Type>(get_view_sfc(), counts, gp_np, proc, xDim, yDim, zDim, periodic));

            scan_boundary_ = ViewVectorType<Integer>("scan_boundary", rank_size + 1);
            incl_excl_scan(0, rank_size, counts, scan_boundary_);
            const Integer candidates = get_last_value(scan_boundary_);
            ViewVectorType<KeyType> local_indices("boundary_candidate_indices", candidates);
            ViewVectorType<Integer> insertion_offsets("boundary_insertion_offsets", rank_size + 1);
            Kokkos::deep_copy(insertion_offsets, scan_boundary_);
            Kokkos::parallel_for("insert_boundary_candidates", get_chunk_size(),
                CountOrInsertBoundaryPerRank<Type, true>(get_view_sfc(), local_indices, insertion_offsets,
                                                        gp_np, proc, xDim, yDim, zDim, periodic));

            // Sort explicit pairs: packing peer IDs into cell keys can overflow or collide.
            ViewVectorType<BoundaryCell> pairs("boundary_peer_cells", candidates);
            const auto scan = scan_boundary_;
            const Integer my_rank = proc;
            Kokkos::parallel_for("boundary_pairs", candidates, KOKKOS_LAMBDA(const Integer i) {
                pairs(i) = {find_owner_processor(scan, i, 1, my_rank), local_indices(i)};
            });
            Kokkos::sort(pairs, BoundaryCellLess{});
            ViewVectorType<bool> unique("unique_boundary_cells", candidates);
            Kokkos::deep_copy(counts, KeyType(0));
            Integer unique_count = 0;
            Kokkos::parallel_reduce("count_unique_boundary_cells", candidates,
                KOKKOS_LAMBDA(const Integer i, Integer& total) {
                    const auto cell = pairs(i);
                    const bool first = i == 0 || cell.peer != pairs(i - 1).peer || cell.local != pairs(i - 1).local;
                    unique(i) = first;
                    if (first) {
                        ++total;
                        Kokkos::atomic_increment(&counts(cell.peer));
                    }
                }, unique_count);
            ViewVectorType<Integer> unique_offsets("boundary_unique_offsets", candidates);
            exclusive_bool_scan(0, candidates, unique_offsets, unique);
            boundary_ = ViewVectorType<KeyType>("boundary", unique_count);
            boundary_lsfc_index_ = ViewVectorType<KeyType>("boundary_local_indices", unique_count);
            const auto boundary = boundary_;
            const auto boundary_indices = boundary_lsfc_index_;
            const auto local_sfc = local_sfc_;
            Kokkos::parallel_for("compact_unique_boundary_cells", candidates, KOKKOS_LAMBDA(const Integer i) {
                if (unique(i)) {
                    const Integer out = unique_offsets(i);
                    boundary_indices(out) = pairs(i).local;
                    boundary(out) = local_sfc(pairs(i).local);
                }
            });
            incl_excl_scan(0, rank_size, counts, scan_boundary_);
        }

        void exchange_ghost_counts(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            int proc_num = rank(context);
            int size = num_ranks(context);
            // copy the boundary elements to the host
            scan_send_mirror = create_mirror_view(get_view_scan_boundary());
            Kokkos::deep_copy(scan_send_mirror, get_view_scan_boundary());

            std::vector<Integer> send_count(size, 0);
            std::vector<Integer> receive_count(size, 0);

            for (int i = 0; i < size; ++i) {
                Integer count = scan_send_mirror(i + 1) - scan_send_mirror(i);
                if (count > 0) {
                    send_count[i] = count;
                    receive_count[i] = count;
                }
            }

            context->distributed->i_send_recv_vec(send_count, receive_count);

            // create the scan recv mirror view from the receive count
            reserve_scan_ghost(size + 1);

            scan_recv_mirror = create_mirror_view(get_view_scan_ghost());
            make_scan_index_mirror(scan_recv_mirror, receive_count);
            Kokkos::deep_copy(get_view_scan_ghost(), scan_recv_mirror);
        }

        void exchange_ghost_layer(const context &context) {
            using namespace Kokkos;

            Kokkos::Timer timer;

            int size = num_ranks(context);

            Integer ghost_size = scan_recv_mirror(size);

            reserve_ghost(ghost_size);

            context->distributed->i_send_recv_view(
                get_view_ghost(), scan_recv_mirror.data(), get_view_boundary(), scan_send_mirror.data());

            std::cout << "MPI send receive for the mesh ghost layer done." << std::endl;
        }

        template <Integer Type>
        void create_ghost_layer() {
            const context &context = get_context();

            build_boundary_element_sets<Type>();

            Kokkos::fence();

            exchange_ghost_counts(context);
            exchange_ghost_layer(context);

            std::cout << "Finished building the ghost layer (boundary element set). Rank: " << get_proc() << std::endl;
        }

        void print_ghost_layer() {
            using namespace Kokkos;

            const context &context = get_context();
            int proc_num = rank(context);
            int rank_size = num_ranks(context);

            Integer ghost_size = scan_recv_mirror(rank_size);

            auto sv = get_view_ghost();
            auto scv = get_view_scan_ghost();
            parallel_for(
                "print set", ghost_size, KOKKOS_LAMBDA(const Integer i) {
                    const Integer rank = find_owner_processor(scv, i, 1, proc_num);

                    printf(" ghost_sfc: %i - %li - proc: %li - rank: %li\n", i, sv(i), rank, proc_num);
                });
        }

        void reserve_ghost(const KeyType n_elements) { ghost_ = ViewVectorType<KeyType>("ghost_", n_elements); }

        void reserve_scan_ghost(const Integer n_elements) {
            scan_ghost_ = ViewVectorType<Integer>("scan_ghost_", n_elements);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_recv_mirror() const { return scan_recv_mirror; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const { return scan_send_mirror; }

        const context &get_context() const { return ctx; }

        [[nodiscard]] bool generate_sfc_to_local_map() {
            auto map = UnorderedMap<KeyType, Integer>(get_chunk_size());
            set_sfc_to_local_map(map);
            // copies needed because of cuda lambda functions. Issue: Copied as *this.{} host pointer.
            auto element_view = get_view_sfc();
            auto sfc_map = get_sfc_to_local_map();
            Kokkos::parallel_for(
                "generate_sfc_to_local_map", get_chunk_size(), MARS_LAMBDA(const Integer i) {
                    sfc_map.insert(element_view(i), i);
                });
            return !sfc_map.failed_insert() && sfc_map.size() == get_chunk_size();
        }
    private:
        using DeviceView::elements_size_;
        using DeviceView::points_size_;
        using DeviceView::boundary_;
        using DeviceView::boundary_lsfc_index_;
        using DeviceView::scan_boundary_;
        using DeviceView::ghost_;
        using DeviceView::scan_ghost_;
        using DeviceView::local_sfc_;
        using DeviceView::sfc_to_local_map_;
        using DeviceView::gp_np;
        using DeviceView::xDim;
        using DeviceView::yDim;
        using DeviceView::zDim;
        using DeviceView::chunk_size_;
        using DeviceView::proc;
        using DeviceView::periodic;
        using DeviceView::combinations;
        using DeviceView::elements_;
        using DeviceView::points_;
        using DeviceView::active_;

        // careful: not a device pointer!
        const context &ctx;
        /* global to local map for the mesh elem indices. */
        UnorderedMap<Integer, Integer> global_to_local_map_;
        // mirror view on the mesh scan boundary view used for the mpi send
        ViewVectorType<Integer>::HostMirror scan_send_mirror;
        // mirror view on the scan_ghost view for the mpi receive
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;
        std::size_t partition_generated_ = 0, partition_received_ = 0;
    };

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedMesh1 =
        mars::Mesh<1, 1, DistributedImplementation, Simplex<1, 1, DistributedImplementation>, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedMesh2 =
        mars::Mesh<2, 2, DistributedImplementation, Simplex<2, 2, DistributedImplementation>, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedMesh3 =
        mars::Mesh<3, 3, DistributedImplementation, Simplex<3, 3, DistributedImplementation>, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedMesh4 =
        mars::Mesh<4, 4, DistributedImplementation, Simplex<4, 4, DistributedImplementation>, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedQuad4Mesh = mars::Mesh<2, 2, DistributedImplementation, Quad4DElem, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>>
    using DistributedHex8Mesh = mars::Mesh<3, 3, DistributedImplementation, Hex8DElem, KeyType>;

    template <typename KeyType = MortonKey<Unsigned>,
              Integer Type = ElementType::Hex8,
              Integer Dim = (Type == ElementType::Quad4) ? 2 : 3,
              Integer ManifoldDim = Dim>
    using DistributedMesh =
        mars::Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>, KeyType>;
}  // namespace mars
#endif  // MARS_MESH_HPP
