#ifndef GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_
#define GENERATION_MARS_DISTRIBUTED_DofHandler_HPP_

#include <cstdlib>
#include <type_traits>
#include "mars_base.hpp"
#include "mars_distributed_mesh_generation.hpp"
#include "mars_utils_kokkos.hpp"
#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_dof.hpp"

namespace mars {

    template <Integer degree, Integer Label, Integer Type>
    struct NumDofs {};

    template <Integer degree, Integer Label>
    struct NumDofs<degree, Label, ElementType::Quad4> {
        static constexpr Integer Dim = 2;

        static constexpr Integer num_corners = power_of_2(Dim);
        static constexpr Integer num_faces = 2 * Dim;
        static constexpr Integer num_edges = 0;

        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = 0;
        static constexpr Integer total_elem_dofs = (degree + 1) * (degree + 1);

        static constexpr Integer elem_dofs() noexcept {
            auto result = 0;

            if (Label & DofLabel::lCorner) {
                result += corner_dofs * num_corners;
            }
            if (Label & DofLabel::lVolume) {
                result += volume_dofs;
            }
            if (Label & DofLabel::lFace) {
                result += face_dofs * num_faces;
            }
            return result;
        }
    };

    template <Integer degree, Integer Label>
    struct NumDofs<degree, Label, ElementType::Hex8> {
        static constexpr Integer Dim = 3;

        static constexpr Integer num_corners = power_of_2(Dim);
        static constexpr Integer num_faces = 2 * Dim;
        static constexpr Integer num_edges = 4 * Dim;

        static constexpr Integer corner_dofs = 1;
        static constexpr Integer face_dofs = (degree - 1) * (degree - 1);
        static constexpr Integer volume_dofs = (degree - 1) * (degree - 1) * (degree - 1);
        static constexpr Integer edge_dofs = (degree - 1);
        static constexpr Integer total_elem_dofs = (degree + 1) * (degree + 1) * (degree + 1);

        static constexpr Integer elem_dofs() noexcept {
            auto result = 0;

            if (Label & DofLabel::lCorner) {
                result += corner_dofs * num_corners;
            }
            if (Label & DofLabel::lEdge) {
                result += edge_dofs * num_edges;
            }
            if (Label & DofLabel::lVolume) {
                result += volume_dofs;
            }
            if (Label & DofLabel::lFace) {
                result += face_dofs * num_faces;
            }
            return result;
        }
    };

    class IDofHandler {
    public:
        virtual MARS_INLINE_FUNCTION ~IDofHandler() {}
        // virtual MARS_FUNCTION const Integer get_ZMax() const;
    };

    template <class Mesh_, Integer degree, Integer Block_ = 1>
    class DofDeviceView {
    public:

        static_assert(degree >= 1, "Structured DOF degree must be positive.");
        static_assert(Block_ >= 0, "Use a positive fixed block or zero for a runtime block.");

        using Mesh = typename Mesh_::DeviceView;
        using KeyType = typename Mesh::KeyType;
        using SfcKeyType = typename Mesh::SfcKeyType;

        using simplex_type = typename Mesh::Elem;

        static constexpr Integer ElemType = simplex_type::ElemType;

        static constexpr Integer dofLabel = DofLabel::lAll;

        static constexpr Integer Block = Block_;

        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        using NDofs = NumDofs<degree, dofLabel, ElemType>;

        static constexpr Integer corner_dofs = NDofs::corner_dofs;
        static constexpr Integer num_corners = NDofs::num_corners;

        static constexpr Integer face_dofs = NDofs::face_dofs;
        static constexpr Integer num_faces = NDofs::num_faces;

        static constexpr Integer volume_dofs = NDofs::volume_dofs;

        static constexpr Integer edge_dofs = NDofs::edge_dofs;
        static constexpr Integer num_edges = NDofs::num_edges;

        static constexpr Integer elem_dofs = NDofs::elem_dofs();

        MARS_INLINE_FUNCTION
        constexpr Integer get_elem_type() const { return simplex_type::ElemType; }

        // computes the component (block) from the block local
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type get_component(
            const Integer local) const {
            return 0;
        }

        // B==0 goes for runtime block input. Useful for interfacing with UtopiaFE.
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type get_component(
            const Integer local) const {
            return local % get_block<B>();
        }

        // computes the base (scalar) local id from the block local
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type compute_base(const Integer local) const {
            return local;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type compute_base(const Integer local) const {
            return local / get_block<B>();
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer compute_block_index(const Integer base_local, const Integer component) const {
            return base_local * get_block<B>() + component;
        }

        template <typename H>
        void owned_iterate(H f) const {
            owned_dof_iterate(f);
        }

        template <typename H, Integer B = Block_>
        void owned_dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", get_block<B>() * global_dof_enum.get_elem_size(), f);
        }

        template <typename H, Integer B = Block_>
        void dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", get_block<B>() * local_dof_enum.get_elem_size(), f);
        }

        template <typename H>
        void base_dof_iterate(H f) const {
            Kokkos::parallel_for("init_initial_cond", local_dof_enum.get_elem_size(), f);
        }

        template <typename H>
        MARS_INLINE_FUNCTION void elem_iterate(H f) const {
            get_mesh().elem_iterate(f);
        }

        MARS_INLINE_FUNCTION Integer get_elem_size() const {
            return get_mesh().get_chunk_size();
        }

        template <typename H>
        MARS_INLINE_FUNCTION void face_iterate(H f) const {
            get_mesh().face_iterate(f);
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc_face_nbh(const Octant &oc, const Integer face_nr) const {
            Octant o = oc.sfc_face_nbh<simplex_type::ElemType>(face_nr);
            return get_sfc_from_octant(o);
        }

        template <typename MT, typename ML>
        static MARS_INLINE_FUNCTION bool exists_in_map(const UnorderedMap<MT, ML> map, const MT value) {
            const auto it = map.find(value);
            return map.valid_at(it);
        }

        MARS_INLINE_FUNCTION
        bool is_local(const KeyType sfc) const {
            auto map = get_local_dof_enum().get_sfc_to_local_map();
            return exists_in_map(map, sfc);
        }

        MARS_INLINE_FUNCTION
        bool is_owned_dof_sfc(const KeyType sfc) const {
            auto map = get_global_dof_enum().get_sfc_to_local_map();
            return exists_in_map(map, sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION bool is_owned(const Integer local) const {
            const auto sfc = local_to_sfc<B>(local);
            return is_owned_dof_sfc(sfc);
        }

        MARS_INLINE_FUNCTION Integer get_eval_value_in_local_map(const KeyType sfc) const {
            auto map = get_local_dof_enum().get_sfc_to_local_map();
            const auto it = map.find(sfc);
            if (map.valid_at(it)) {
                return map.value_at(it);
            }
            return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_eval_value_in_global_map(const KeyType sfc) const {
            auto map = get_global_dof_enum().get_sfc_to_local_map();
            const auto it = map.find(sfc);
            if (map.valid_at(it)) {
                return map.value_at(it);
            }
            return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Integer>::type local_to_owned_dof(
            const Integer local) const {
            const Integer comp_local = get_component<B>(local);
            const auto sfc = local_to_sfc<B>(local);
            const Integer base_owned = get_eval_value_in_global_map(sfc);
            return compute_block_index<B>(base_owned, comp_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Integer>::type local_to_owned_dof(
            const Integer local) const {
            const auto sfc = local_to_sfc<B>(local);
            return get_eval_value_in_global_map(sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_owned_index(const Integer local) const {
            return local_to_owned_dof<B>(local);
        }

        MARS_INLINE_FUNCTION
        Dof sfc_to_global_dof(const KeyType sfc) const {
            Dof dof;
            const auto sfc_lid = get_eval_value_in_global_map(sfc);
            if (sfc_lid != INVALID_INDEX) {
                const Integer proc = get_mesh().get_proc();
                const auto gid = sfc_lid + global_dof_offset(proc);
                dof.set_gid(gid);
                dof.set_proc(proc);
            } else {
                const auto it = ghost_local_to_global_map.find(sfc);
                if (!ghost_local_to_global_map.valid_at(it)) {
                    dof.set_invalid();
                } else {
                    dof = ghost_local_to_global_map.value_at(it);
                }
            }
            return dof;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 1, Dof>::type local_to_global_dof(const Integer local) const {
            const auto local_sfc = local_to_sfc<B>(local);
            return sfc_to_global_dof(local_sfc);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B != 1, Dof>::type local_to_global_dof(const Integer local) const {
            // use the local to sfc view (elements) to get the sfc of the local numbering.
            const Integer component = get_component<B>(local);

            const auto local_sfc = local_to_sfc<B>(local);
            auto dof = sfc_to_global_dof(local_sfc);
            dof.set_gid(compute_block_index<B>(dof.get_gid(), component));
            return dof;
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_global(const KeyType sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_global(const Integer local) const {
            Dof dof = local_to_global_dof<B>(local);
            if (dof.is_valid())
                return dof.get_gid();
            else
                return INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_global_proc(const KeyType sfc) const {
            Dof dof = sfc_to_global_dof(sfc);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_global_proc(const Integer local) const {
            Dof dof = local_to_global_dof<B>(local);
            if (dof.is_valid())
                return dof.get_proc();
            else
                return INVALID_INDEX;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_local_label(const Integer local) const {
            const Integer base_local = compute_base<B>(local);
            return get_local_dof_enum().get_view_element_labels()(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION KeyType local_to_sfc(const Integer local) const {
            auto base = compute_base<B>(local);
            return get_local_dof_enum().get_view_elements()(base);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION KeyType owned_to_sfc(const Integer owned) const {
            auto base = compute_base<B>(owned);
            return get_global_dof_enum().get_view_elements()(base);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer local_to_owned(const Integer local) const {
            return local_to_owned_dof<B>(local);
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_owned(const KeyType sfc) const { return get_eval_value_in_global_map(sfc); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_owned(const KeyType sfc, const Integer component) const {
            return compute_block_index<Block>(get_eval_value_in_global_map(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const KeyType sfc) const { return get_eval_value_in_local_map(sfc); }

        MARS_INLINE_FUNCTION
        Integer sfc_to_local(const KeyType sfc, const Integer component) const {
            return compute_block_index<Block>(get_eval_value_in_local_map(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer get_local_index(const KeyType sfc) const { return sfc_to_local(sfc); }

        MARS_INLINE_FUNCTION
        Integer get_local_index(const KeyType sfc, const Integer component) const {
            return sfc_to_local(sfc, component);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_sfc(const KeyType sfc, double *point) const {
            get_vertex_coordinates_from_sfc<Type, SfcKeyType>(sfc,
                                                              point,
                                                              get_local_dof_enum().get_XDim(),
                                                              get_local_dof_enum().get_YDim(),
                                                              get_local_dof_enum().get_ZDim());
        }

        template <Integer Type, Integer B = Block_>
        MARS_INLINE_FUNCTION void get_dof_coordinates_from_local(const Integer local, double *point) const {
            const auto sfc = local_to_sfc<B>(local);
            return get_dof_coordinates_from_sfc<Type>(sfc, point);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION void get_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_coordinates_from_local<ElemType, B>(local, point);
        }

        MARS_INLINE_FUNCTION void get_base_local_dof_coordinates(const Integer local, double *point) const {
            get_dof_coordinates_from_local<ElemType, 1>(local, point);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION bool boundary_sfc(const KeyType sfc, const Integer FaceNr = -1) const {
            const Integer xdim = get_local_dof_enum().get_XDim();
            const Integer ydim = get_local_dof_enum().get_YDim();
            const Integer zdim = get_local_dof_enum().get_ZDim();

            return is_boundary_sfc<Type, SfcKeyType>(sfc, xdim, ydim, zdim, FaceNr);
        }

        template <Integer Type, Integer B = Block_>
        MARS_INLINE_FUNCTION bool is_boundary(const Integer local, const Integer FaceNr = -1) const {
            const auto sfc = local_to_sfc<B>(local);
            return boundary_sfc<Type>(sfc, FaceNr);
        }

        MARS_INLINE_FUNCTION bool is_boundary_dof(const Integer local, const Integer FaceNr = -1) const {
            return is_boundary<ElemType>(local, FaceNr);
        }

        MARS_INLINE_FUNCTION
        static bool is_separated_block(const Integer block, const Integer component) {
            return (block == -1 || block == component);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_dof_size() const {
            return get_block<B>() * local_dof_enum.get_elem_size();
        }

        MARS_INLINE_FUNCTION Integer get_base_dof_size() const { return local_dof_enum.get_elem_size(); }

        // no need to be vector valued since it will give the same result.
        MARS_INLINE_FUNCTION
        Integer get_local_dof(const Integer i) const {
            assert(i < get_dof_size());
            return i;
        }

        // get the local dof of the owned index.
        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_owned_dof(const Integer owned_dof) const {
            auto base = compute_base<B>(owned_dof);
            auto component = get_component<B>(owned_dof);
            const auto sfc = get_global_dof_enum().get_view_elements()(base);
            return compute_block_index<B>(sfc_to_local(sfc), component);
        }

        MARS_INLINE_FUNCTION
        Integer get_base_owned_dof_size() const { return global_dof_enum.get_elem_size(); }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_owned_dof_size() const {
            return get_block<B>() * global_dof_enum.get_elem_size();
        }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType, SfcKeyType> &get_local_dof_enum() const { return local_dof_enum; }

        MARS_INLINE_FUNCTION
        const SFC<simplex_type::ElemType, SfcKeyType> &get_global_dof_enum() const { return global_dof_enum; }

        MARS_INLINE_FUNCTION
        Integer get_global_dof_offset(const Integer proc) const { return global_dof_offset(proc); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> get_global_dof_offset() const { return global_dof_offset; }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Dof> &get_ghost_lg_map() const { return ghost_local_to_global_map; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_dof(const Integer i) const {
            const auto base = compute_base<B>(i);
            const auto component = get_component<B>(i);
            const auto local = sfc_to_local(get_boundary_dofs()(base));
            assert(INVALID_INDEX != local);
            return compute_block_index<B>(local, component);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_local_index(const Integer i) const {
            return get_boundary_dof<B>(i);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_boundary_dof_size() const {
            return get_block<B>() * get_boundary_dofs().extent(0);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_dof(const Integer i) const {
            const auto base = compute_base<B>(i);
            const auto component = get_component<B>(i);
            const auto local = sfc_to_local(get_ghost_dofs()(base));
            assert(INVALID_INDEX != local);
            return compute_block_index<B>(local, component);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_local_index(const Integer i) const {
            return get_ghost_dof<B>(i);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_ghost_dof_size() const {
            return get_block<B>() * get_ghost_dofs().extent(0);
        }

        template <typename F>
        void ghost_iterate(F f) const {
            Kokkos::parallel_for("ghost_dof_iter", get_ghost_dof_size(), f);
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<KeyType, Integer> get_local_dof_map() const {
            return get_local_dof_enum().get_sfc_to_local_map();
        }

        MARS_INLINE_FUNCTION
        Mesh get_mesh() const { return mesh; }

        MARS_INLINE_FUNCTION
        Integer get_proc() const { return get_mesh().get_proc(); }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_recv() const { return scan_recv; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_scan_send() const { return scan_send; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_orientation(const Integer local_dof) const {
            const Integer base_local = compute_base<B>(local_dof);
            return get_local_dof_enum().get_orientation(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_label(const Integer local_dof) const {
            const Integer base_local = compute_base<B>(local_dof);
            return get_local_dof_enum().get_label(base_local);
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Integer get_owned_label(const Integer owned_dof) const {
            const Integer base_owned = compute_base<B>(owned_dof);
            return get_global_dof_enum().get_label(base_owned);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_boundary_dofs() const { return boundary_dofs_sfc; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_ghost_dofs() const { return ghost_dofs_sfc; }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION Octant get_octant_from_local(const Integer local) const {
            const auto sfc = local_to_sfc<B>(local);
            return get_octant_from_sfc(sfc);
        }

        MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const KeyType sfc) const {
            return mars::get_octant_from_sfc<ElemType, SfcKeyType>(sfc);
        }

        MARS_INLINE_FUNCTION KeyType get_sfc_from_octant(const Octant &o) const {
            return mars::get_sfc_from_octant<ElemType, SfcKeyType>(o);
        }

        MARS_INLINE_FUNCTION Integer get_local_sfc_from_octant(const Octant &o) const {
            const auto sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o) const {
            const auto sfc = get_sfc_from_octant(o);
            return is_local(sfc) ? sfc_to_local(sfc) : INVALID_INDEX;
        }

        MARS_INLINE_FUNCTION Integer get_local_from_octant(const Octant &o, const Integer component) const {
            auto local = get_local_from_octant(o);
            return compute_block_index<Block>(local, component);
        }

        MARS_INLINE_FUNCTION
        Integer get_XMax() const { return Degree * get_mesh().get_XDim(); }

        MARS_INLINE_FUNCTION
        Integer get_YMax() const { return Degree * get_mesh().get_YDim(); }

        MARS_INLINE_FUNCTION
        Integer get_ZMax() const { return Degree * get_mesh().get_ZDim(); }

        MARS_INLINE_FUNCTION
        void set_block(const Integer b) {
            if (b <= 0) Kokkos::abort("Structured DOF block size must be positive.");
            block = b;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION typename std::enable_if<B == 0, Integer>::type get_block() const {
            if (block <= 0) Kokkos::abort("Structured DOF block size must be positive.");
            return block;
        }

        template <Integer B = Block_>
        MARS_INLINE_FUNCTION constexpr typename std::enable_if<(B > 0), Integer>::type get_block() const {
            return B;
        }

        template <class VIEW, class V>
        MARS_INLINE_FUNCTION void vector_apply_constraints(const Integer row, VIEW v, const V value) const {
            auto diag_row = local_to_owned_index(row);
            v(diag_row, 0) = value;
        }

        template <class VIEW>
        MARS_INLINE_FUNCTION void apply_zero_constraints(const Integer row, VIEW v) const {
            vector_apply_constraints(row, v, 0);
        }
        explicit DofDeviceView(Mesh m) : mesh(m) {}
        MARS_INLINE_FUNCTION DofDeviceView device_view() const { return *this; }

    protected:

        // manages host and device mesh objects.

        Mesh mesh;

        // local dof enumeration containing the local to sfc and sfc to local views.
        SFC<simplex_type::ElemType, SfcKeyType> local_dof_enum;
        // locally owned dof enumeration containing the global to sfc and sfc to global views.
        SFC<simplex_type::ElemType, SfcKeyType> global_dof_enum;
        // global offset used to calc the global numbering of the dofs.
        ViewVectorType<Integer> global_dof_offset;

        // the local to global mapp for the locally owned is mapped from the sfc_to global view
        // for the ghost dofs is used the following kokkos unordered map.
        UnorderedMap<Integer, Dof> ghost_local_to_global_map;

        // dofs sfc for the boundary dofs.
        ViewVectorType<Integer> boundary_dofs_sfc;
        // mirror view on the dof scan boundary view used as offset for the mpi send dofs.
        ViewVectorType<Integer> scan_send;

        // ghost dofs received from other procs's boundary sfcs.
        ViewVectorType<Integer> ghost_dofs_sfc;
        // mirror view on the scan_ghost view used as an offset to receive ghost dofs.
        ViewVectorType<Integer> scan_recv;

        Integer block = 1;
    };

    template <class Mesh_, Integer degree, Integer Block_ = 1>
    class DofHandler : public DofDeviceView<Mesh_, degree, Block_>, public IDofHandler {
    public:
        using DeviceView = DofDeviceView<Mesh_, degree, Block_>;

        static_assert(degree >= 1, "Structured DOF degree must be positive.");
        static_assert(Block_ >= 0, "Use a positive fixed block or zero for a runtime block.");

        using Mesh = typename Mesh_::DeviceView;
        using KeyType = typename Mesh::KeyType;
        using SfcKeyType = typename Mesh::SfcKeyType;

        using simplex_type = typename Mesh::Elem;

        static constexpr Integer ElemType = simplex_type::ElemType;

        static constexpr Integer dofLabel = DofLabel::lAll;

        static constexpr Integer Block = Block_;

        static constexpr Integer Degree = degree;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;

        using NDofs = NumDofs<degree, dofLabel, ElemType>;

        static constexpr Integer corner_dofs = NDofs::corner_dofs;
        static constexpr Integer num_corners = NDofs::num_corners;

        static constexpr Integer face_dofs = NDofs::face_dofs;
        static constexpr Integer num_faces = NDofs::num_faces;

        static constexpr Integer volume_dofs = NDofs::volume_dofs;

        static constexpr Integer edge_dofs = NDofs::edge_dofs;
        static constexpr Integer num_edges = NDofs::num_edges;

        static constexpr Integer elem_dofs = NDofs::elem_dofs();
        using DeviceView::apply_zero_constraints;
        using DeviceView::base_dof_iterate;
        using DeviceView::boundary_sfc;
        using DeviceView::compute_base;
        using DeviceView::compute_block_index;
        using DeviceView::dof_iterate;
        using DeviceView::elem_iterate;
        using DeviceView::exists_in_map;
        using DeviceView::face_iterate;
        using DeviceView::get_XMax;
        using DeviceView::get_YMax;
        using DeviceView::get_ZMax;
        using DeviceView::get_base_dof_size;
        using DeviceView::get_base_local_dof_coordinates;
        using DeviceView::get_base_owned_dof_size;
        using DeviceView::get_block;
        using DeviceView::get_boundary_dof;
        using DeviceView::get_boundary_dof_size;
        using DeviceView::get_boundary_dofs;
        using DeviceView::get_boundary_local_index;
        using DeviceView::get_component;
        using DeviceView::get_dof_coordinates_from_local;
        using DeviceView::get_dof_coordinates_from_sfc;
        using DeviceView::get_dof_size;
        using DeviceView::get_elem_size;
        using DeviceView::get_elem_type;
        using DeviceView::get_eval_value_in_global_map;
        using DeviceView::get_eval_value_in_local_map;
        using DeviceView::get_ghost_dof;
        using DeviceView::get_ghost_dof_size;
        using DeviceView::get_ghost_dofs;
        using DeviceView::get_ghost_lg_map;
        using DeviceView::get_ghost_local_index;
        using DeviceView::get_global_dof_enum;
        using DeviceView::get_global_dof_offset;
        using DeviceView::get_label;
        using DeviceView::get_local_dof;
        using DeviceView::get_local_dof_coordinates;
        using DeviceView::get_local_dof_enum;
        using DeviceView::get_local_dof_map;
        using DeviceView::get_local_from_octant;
        using DeviceView::get_local_index;
        using DeviceView::get_local_label;
        using DeviceView::get_local_sfc_from_octant;
        using DeviceView::get_mesh;
        using DeviceView::get_octant_from_local;
        using DeviceView::get_octant_from_sfc;
        using DeviceView::get_orientation;
        using DeviceView::get_owned_dof;
        using DeviceView::get_owned_dof_size;
        using DeviceView::get_owned_label;
        using DeviceView::get_proc;
        using DeviceView::get_sfc_face_nbh;
        using DeviceView::get_sfc_from_octant;
        using DeviceView::get_view_scan_recv;
        using DeviceView::get_view_scan_send;
        using DeviceView::ghost_iterate;
        using DeviceView::is_boundary;
        using DeviceView::is_boundary_dof;
        using DeviceView::is_local;
        using DeviceView::is_owned;
        using DeviceView::is_owned_dof_sfc;
        using DeviceView::is_separated_block;
        using DeviceView::local_to_global;
        using DeviceView::local_to_global_dof;
        using DeviceView::local_to_global_proc;
        using DeviceView::local_to_owned;
        using DeviceView::local_to_owned_dof;
        using DeviceView::local_to_owned_index;
        using DeviceView::local_to_sfc;
        using DeviceView::owned_dof_iterate;
        using DeviceView::owned_iterate;
        using DeviceView::owned_to_sfc;
        using DeviceView::set_block;
        using DeviceView::sfc_to_global;
        using DeviceView::sfc_to_global_dof;
        using DeviceView::sfc_to_global_proc;
        using DeviceView::sfc_to_local;
        using DeviceView::sfc_to_owned;
        using DeviceView::vector_apply_constraints;
        explicit DofHandler(const Mesh_& m) : DeviceView(m.device_view()), ctx(m.get_context()) {}
        DeviceView device_view() const { return static_cast<const DeviceView&>(*this); }
        const context& get_context() const { return ctx; }

        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer> enum_corner(const Octant &oc,
                                                                                               const int i,
                                                                                               const int j) const {
            Octant o;
            // get the next corner using the elem sfc
            o.x = oc.x + i;
            o.y = oc.y + j;
            // convert the octant value into the new nodal sfc system
            o.x *= degree;
            o.y *= degree;
            auto sfc = get_sfc_from_octant(o);

            return sfc_to_local(sfc);
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer> enum_corner(const Octant &oc,
                                                                                              const int i,
                                                                                              const int j,
                                                                                              const int k) const {
            Octant o;
            // get the next corner using the elem sfc
            o.x = oc.x + i;
            o.y = oc.y + j;
            o.z = oc.z + k;
            // convert the octant value into the new nodal sfc system
            o.x *= degree;
            o.y *= degree;
            o.z *= degree;
            auto sfc = get_sfc_from_octant(o);

            return sfc_to_local(sfc);
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Octant> enum_face_corner(
            const Octant &oc,
            const int dir) {
            Octant face_cornerA;

            if (part == 0) {
                face_cornerA.x = oc.x + dir;
                face_cornerA.y = oc.y;
            }
            if (part == 1) {
                face_cornerA.x = oc.x + (dir ^ 1);
                face_cornerA.y = oc.y + 1;
            }

            return face_cornerA;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Octant> enum_face_corner(
            const Octant &oc,
            const int dir) {
            Octant face_cornerA;
            build_starting_face_corner<Type, part>(face_cornerA, dir, oc);
            return face_cornerA;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer>
        sfc_face_node(const Octant &face_cornerA, const int j, const int dir) {
            int sign = (part == 1) ? -1 : 1;
            // dir is 1 for y direction and 0 for x direction
            Octant o;
            if (dir == 1) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + sign * (j + 1);
            }
            if (dir == 0) {
                o.x = degree * face_cornerA.x + sign * (j + 1);
                o.y = degree * face_cornerA.y;
            }
            auto sfc = mars::get_sfc_from_octant<Type, SfcKeyType>(o);

            return sfc;
        }

        template <Integer part, Integer Type>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        sfc_face_node(const Octant &face_cornerA, const int j, const int dir) {
            return build_face_node<Type, part>(mars::get_sfc_from_octant<Type, SfcKeyType>, face_cornerA, j);
        }

        template <Integer part, Integer Type>
        MARS_INLINE_FUNCTION Integer enum_face_node(const Octant &face_cornerA, const int j, const int dir) const {
            return sfc_to_local(sfc_face_node<part, Type>(face_cornerA, j, dir));
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void>
        build_starting_face_corner(Octant &face_cornerA, const int side, const Octant &oc) {
            // find the starting corner "face_cornerA" of the face and use the direction
            const int val = side ^ 1;
            if (dir == 0) {
                face_cornerA.x = oc.x + val;
                face_cornerA.y = oc.y;
                face_cornerA.z = oc.z;
            }
            if (dir == 1) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y + val;
                face_cornerA.z = oc.z;
            }
            if (dir == 2) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y;
                face_cornerA.z = oc.z + val;
            }
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, void>
        build_starting_face_corner(Octant &face_cornerA, const int side, const Octant &oc) {
            // find the starting corner "face_cornerA" of the face and use the direction
            const int val = side ^ 1;
            if (dir == 0) {
                face_cornerA.x = oc.x + val;
                face_cornerA.y = oc.y;
            }
            if (dir == 1) {
                face_cornerA.x = oc.x;
                face_cornerA.y = oc.y + val;
            }
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer
        process_face_corner(Octant &face_cornerA, const Mesh mesh, const int side, const Octant &oc) {
            Integer face_nr;
            if (side == 0)
                face_nr = 2 * dir + 1;
            else
                face_nr = 2 * dir;

            Octant nbh_oc = mesh.get_octant_face_nbh(oc, face_nr);
            Integer owner_proc = -1;

            // find owner proc method returns an invalid result with invalid octant.
            if (nbh_oc.is_valid()) {
                Integer enc_oc = mars::get_sfc_from_octant<Type, SfcKeyType>(nbh_oc);
                owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, mesh.get_proc());
            }

            build_starting_face_corner<Type, dir>(face_cornerA, side, oc);

            return owner_proc;
        }

        template <Integer Type, Integer dir, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Integer>
        build_face_node(F f, const Octant &face_cornerA, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
            }
            return f(o);
        }

        template <Integer Type, Integer dir, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        build_face_node(F f, const Octant &face_cornerA, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z + j + 1;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z + j + 1;
            }
            if (dir == 2) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z;
            }
            return f(o);
        }

        template <Integer Type, Integer dir>
        static MARS_INLINE_FUNCTION Integer process_face_node(const Octant &face_cornerA, const int j) {
            return build_face_node<Type, dir>(mars::get_sfc_from_octant<Type, SfcKeyType>, face_cornerA, j);
        }

        template <Integer Type, typename F>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Integer>
        process_edge_node(F f, const Octant &face_cornerA, const Integer dir, const int j) {
            Octant o;
            if (dir == 0) {
                o.x = degree * face_cornerA.x + j + 1;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z;
            }
            if (dir == 1) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y + j + 1;
                o.z = degree * face_cornerA.z;
            }
            if (dir == 2) {
                o.x = degree * face_cornerA.x;
                o.y = degree * face_cornerA.y;
                o.z = degree * face_cornerA.z + j + 1;
            }
            return f(o);
        }

        // TODO: Use this instead of the following functor when switching to c++17 and Cuda 11.

        struct MaxOneRingProc {
            Integer &max_proc;
            ViewVectorType<KeyType> gp;
            Integer proc;

            MARS_INLINE_FUNCTION
            MaxOneRingProc(ViewVectorType<KeyType> g, Integer p, Integer &mx) : gp(g), proc(p), max_proc(mx) {}

            MARS_INLINE_FUNCTION
            void operator()(const Octant &nbh_oc) const {
                // find owner proc method returns an invalid result with invalid octant.
                if (nbh_oc.is_valid()) {
                    Integer enc_oc = mars::get_sfc_from_octant<ElemType, SfcKeyType>(nbh_oc);
                    Integer owner_proc = find_owner_processor(gp, enc_oc, 2, proc);
                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }
        };

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T == ElementType::Hex8, void> edge_iterate(const KeyType sfc,
                                                                                                const Mesh mesh,
                                                                                                const Integer i,
                                                                                                F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            for (int edge = 0; edge < 4 * ManifoldDim; ++edge) {
                // get the starting node of the edge to be used to compute the dofs contained in that edge.
                const Integer direction = oc.get_edge_direction(edge);
                const Octant start = mesh.get_octant_edge_start(oc, edge);

                Integer max_proc = -1;
                mesh.get_one_ring_edge_nbhs(
                    start, direction, MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));

                for (int j = 0; j < edge_dofs; j++) {
                    Integer dof_sfc =
                        process_edge_node<ElemType>(mars::get_sfc_from_octant<ElemType, SfcKeyType>, start, direction, j);
                    f(mesh, i, dof_sfc, max_proc, direction);
                }
            }
        }

        template <typename F, Integer T = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<T != ElementType::Hex8, void> edge_iterate(const KeyType sfc,
                                                                                                const Mesh mesh,
                                                                                                const Integer i,
                                                                                                F f) {}

        template <typename CornerF, typename FaceF, typename VolumeF, typename EdgeF>
        static MARS_INLINE_FUNCTION void elem_dof_iterate(const KeyType sfc,
                                                          const Mesh mesh,
                                                          const Integer i,
                                                          CornerF cf,
                                                          FaceF ff,
                                                          VolumeF vf,
                                                          EdgeF ef) {
            if (corner_dofs > 0) {
                corner_iterate(sfc, mesh, i, cf);
            }
            if (face_dofs > 0) {
                // templetized with 1 value means filtered on the direction (No dir is forwarded to the f call).
                face_dir_iterate<1>(sfc, mesh, i, ff);
            }
            if (volume_dofs > 0) {
                volume_iterate(sfc, mesh, i, vf);
            }
            if (edge_dofs > 0) {
                edge_iterate(sfc, mesh, i, IterateFilter<EdgeF, 1>(ef));
            }
        }
        struct BoundaryDof {
            Integer peer;
            Integer owned;
        };
        struct BoundaryDofLess {
            MARS_INLINE_FUNCTION bool operator()(const BoundaryDof& a, const BoundaryDof& b) const {
                return a.peer < b.peer || (a.peer == b.peer && a.owned < b.owned);
            }
        };
        struct InsertBoundaryDof {
            ViewVectorType<BoundaryDof> pairs;
            UnorderedMap<KeyType, Integer> owned;
            Integer peer;
            Integer offset;
            Integer& used;
            MARS_INLINE_FUNCTION void operator()(const Mesh&, Integer, KeyType key, Integer = 0) const {
                const auto slot = owned.find(key);
                if (owned.valid_at(slot)) pairs(offset + used++) = {peer, owned.value_at(slot)};
            }
        };

        void build_boundary_dof_sets(ViewVectorType<Integer>& scan_boundary,
                                     ViewVectorType<Integer>& boundary_lsfc_index,
                                     ViewVectorType<Integer>& sender_ranks_scan,
                                     Integer nbh_rank_size) {
            using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
            const auto mesh_view = get_mesh();
            const Integer cells = mesh_view.get_boundary_size();
            if (cells > std::numeric_limits<Integer>::max() / elem_dofs)
                mpi::abort_invalid_count(get_context()->distributed->mpi_comm(), "boundary DOF scratch overflow");
            const Integer capacity = cells * elem_dofs;
            ViewVectorType<BoundaryDof> pairs("boundary_dof_candidates", capacity);
            const auto owned = get_global_dof_enum().get_sfc_to_local_map();
            Kokkos::parallel_for("boundary_dof_candidates", Policy(0, cells), KOKKOS_LAMBDA(Integer i) {
                const Integer peer = find_owner_processor(mesh_view.get_view_scan_boundary(), i, 1,
                                                          mesh_view.get_proc());
                Integer used = 0;
                const InsertBoundaryDof insert{pairs, owned, sender_ranks_scan(peer), i * elem_dofs, used};
                elem_dof_iterate(mesh_view.get_boundary_sfc(i), mesh_view, i, insert, insert, insert, insert);
                for (; used < elem_dofs; ++used) pairs(i * elem_dofs + used) = {nbh_rank_size, 0};
            });
            Kokkos::sort(pairs, BoundaryDofLess{});
            const auto distinct = KOKKOS_LAMBDA(Integer i) {
                return pairs(i).peer < nbh_rank_size && (i == 0 || pairs(i).peer != pairs(i - 1).peer ||
                                                         pairs(i).owned != pairs(i - 1).owned);
            };
            Integer unique_count = 0;
            Kokkos::parallel_reduce("count_boundary_dofs", Policy(0, capacity),
                KOKKOS_LAMBDA(Integer i, Integer& n) { if (distinct(i)) ++n; }, unique_count);
            boundary_dofs_sfc = ViewVectorType<Integer>("boundary_sfc_dofs", unique_count);
            boundary_lsfc_index = ViewVectorType<Integer>("boundary_owned_dofs", unique_count);
            ViewVectorType<Integer> peer_counts("boundary_dof_peer_counts", nbh_rank_size);
            auto boundary_keys = boundary_dofs_sfc;
            auto owned_keys = get_global_dof_enum().get_view_elements();
            Kokkos::parallel_scan("compact_boundary_dofs", Policy(0, capacity),
                KOKKOS_LAMBDA(Integer i, Integer& n, bool final) {
                    if (distinct(i)) {
                        if (final) {
                            boundary_lsfc_index(n) = pairs(i).owned;
                            boundary_keys(n) = owned_keys(pairs(i).owned);
                            Kokkos::atomic_increment(&peer_counts(pairs(i).peer));
                        }
                        ++n;
                    }
                });
            scan_boundary = ViewVectorType<Integer>("boundary_dof_offsets", nbh_rank_size + 1);
            incl_excl_scan(0, nbh_rank_size, peer_counts, scan_boundary);
        }

        template <Integer Type>
        static MARS_INLINE_FUNCTION Integer process_corner(const Mesh mesh, const Octant &o) {
            Integer xDim = mesh.get_XDim();
            Integer yDim = mesh.get_YDim();
            Integer zDim = mesh.get_ZDim();

            Integer max_proc = -1;

            Octant one_ring[simplex_type::ElemType];
            o.one_ring_corner_nbhs<Type>(one_ring, xDim, yDim, zDim, mesh.is_periodic());

            for (int k = 0; k < Type; k++) {
                if (one_ring[k].is_valid()) {
                    auto enc_oc = mars::get_sfc_from_octant<Type, SfcKeyType>(one_ring[k]);

                    Integer owner_proc = find_owner_processor(mesh.get_view_gp(), enc_oc, 2, mesh.get_proc());

                    if (owner_proc > max_proc) {
                        max_proc = owner_proc;
                    }
                }
            }

            return max_proc;
        }

        template <typename F>
        static MARS_INLINE_FUNCTION void corner_iterate(const KeyType sfc, const Mesh mesh, const Integer index, F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            corner_octant_transform(oc, mesh, index, f);
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        corner_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    Octant o;
                    // get the next corner using the elem sfc
                    o.x = oc.x + i;
                    o.y = oc.y + j;

                    Integer max_proc = -1;
                    mesh.get_one_ring_corner_nbhs(o, MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));

                    // convert the octant value into the new nodal sfc system
                    o.x *= degree;
                    o.y *= degree;
                    auto dof_sfc = mars::get_sfc_from_octant<ElemType, SfcKeyType>(o);

                    f(mesh, index, dof_sfc, max_proc);
                }
            }
        }

        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        corner_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f) {
            for (int i = 0; i < corner_dofs + 1; i++) {
                for (int j = 0; j < corner_dofs + 1; j++) {
                    for (int k = 0; k < corner_dofs + 1; k++) {
                        Octant o;
                        // get the next corner using the elem sfc
                        o.x = oc.x + i;
                        o.y = oc.y + j;
                        o.z = oc.z + k;

                        Integer max_proc = -1;
                        mesh.get_one_ring_corner_nbhs(o,
                                                       MaxOneRingProc(mesh.get_view_gp(), mesh.get_proc(), max_proc));

                        // convert the octant value into the new nodal sfc system
                        o.x *= degree;
                        o.y *= degree;
                        o.z *= degree;
                        auto dof_sfc = mars::get_sfc_from_octant<ElemType, SfcKeyType>(o);

                        f(mesh, index, dof_sfc, max_proc);
                    }
                }
            }
        }

       // sfc | octant_from_sfc | volume_octant_transform. Composable functions.
        template <typename F, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION void volume_iterate(const KeyType sfc, const Mesh mesh, const Integer index, F f) {
            Octant oc = mesh.octant_from_sfc(sfc);
            volume_octant_transform(oc, mesh, index, f, mars::get_sfc_from_octant<ElemType, SfcKeyType>);
        }

        template <typename F, typename G, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Quad4, void>
        volume_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f, G g) {
            Octant o;
            // go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++) {
                for (int j = 0; j < degree - 1; j++) {
                    o.x = degree * oc.x + i + 1;
                    o.y = degree * oc.y + j + 1;
                    f(mesh, index, g(o));
                }
            }
        }

        template <typename F, typename G, Integer ET = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<ET == ElementType::Hex8, void>
        volume_octant_transform(const Octant &oc, const Mesh mesh, const Integer index, F f, G g) {
            Octant o;
            // go through all the inside dofs for the current element
            for (int i = 0; i < degree - 1; i++) {
                for (int j = 0; j < degree - 1; j++) {
                    for (int k = 0; k < degree - 1; k++) {
                        o.x = degree * oc.x + i + 1;
                        o.y = degree * oc.y + j + 1;
                        o.z = degree * oc.z + k + 1;
                        f(mesh, index, g(o));
                    }
                }
            }
        }

        // Iterate filter redirecting the call with filtered number of parameters to different functors.
        template <typename F, Integer Params = 0>
        struct IterateFilter {
            F f;

            MARS_INLINE_FUNCTION
            IterateFilter(F fun) : f(fun) {}

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const KeyType dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                f(mesh, i, dof_sfc, owner_proc, dir);
            }
        };

        // Could have specialized on a bool for dir or not dir but left it open for further filters.
        // For example one might need a face_iterate with only mesh, i, dofs_sfc parameteres.
        template <typename F>
        struct IterateFilter<F, 1> {
            F f;

            MARS_INLINE_FUNCTION
            IterateFilter(F fun) : f(fun) {}

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const KeyType dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                f(mesh, i, dof_sfc, owner_proc);
            }
        };

        template <Integer Params, typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, void>
        face_dir_iterate(const KeyType sfc, const Mesh mesh, const Integer i, F f) {
            auto filter = IterateFilter<F, Params>(f);
            // dir 0 = x, 1 = y
            face_iterate<0>(sfc, mesh, i, filter);
            face_iterate<1>(sfc, mesh, i, filter);
        }

        template <Integer Params, typename F, Integer Type = ElemType>
        static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void>
        face_dir_iterate(const KeyType sfc, const Mesh mesh, const Integer i, F f) {
            auto filter = IterateFilter<F, Params>(f);
            // dir 0 = x, 1 = y
            // dir 0 = x, 1 = y, 2 = z
            face_iterate<0>(sfc, mesh, i, filter);
            face_iterate<1>(sfc, mesh, i, filter);
            face_iterate<2>(sfc, mesh, i, filter);
        }

        // careful: do not use this function for face iteration when expensive operations are involved. Useful mainly
        // for predicate builder. instead use the locally_owned_face_dof array to iterate over face dof sfcs.
        template <Integer dir, typename F>
        static MARS_INLINE_FUNCTION void face_iterate(const KeyType sfc, const Mesh mesh, const Integer index, F f) {
            // side  0 means origin side and 1 destination side.
            Octant oc = mesh.octant_from_sfc(sfc);

            for (int side = 0; side < 2; ++side) {
                Octant face_cornerA;
                Integer owner_proc = process_face_corner<simplex_type::ElemType, dir>(face_cornerA, mesh, side, oc);

                for (int j = 0; j < face_dofs; j++) {
                    auto dof_sfc = process_face_node<simplex_type::ElemType, dir>(face_cornerA, j);
                    f(mesh, index, dof_sfc, owner_proc, dir);
                }
            }
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == true, Octant>::type get_octant_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return mesh.get_ghost_octant(index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == false, Octant>::type get_octant_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return mesh.get_octant(index);
        }

        static MARS_INLINE_FUNCTION Integer get_ghost_sfc(const Mesh mesh, const Integer index) {
            return mesh.get_ghost_sfc(index);
        }
        static MARS_INLINE_FUNCTION KeyType get_local_sfc(const Mesh mesh, const Integer index) {
            return mesh.get_sfc(index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == true, KeyType>::type get_sfc_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return get_ghost_sfc(mesh, index);
        }

        template <bool G>
        static MARS_INLINE_FUNCTION typename std::enable_if<G == false, KeyType>::type get_sfc_ghost_or_local(
            const Mesh mesh,
            const Integer index) {
            return get_local_sfc(mesh, index);
        }

        // generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <typename T, bool Ghost, Integer Label>
        struct LocalOwnedDofsFill{
            ViewVectorType<KeyType> owned_dofs;
            Integer &owned_index;
            ViewVectorType<KeyType> dofs;
            Integer &index;
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> local_thread_count;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalOwnedDofsFill(ViewVectorType<KeyType> od,
                               Integer &ow_idx,
                               ViewVectorType<KeyType> d,
                               Integer &idx,
                               ViewVectorType<T> gc,
                               ViewVectorType<T> lc,
                               Integer p)
                : owned_dofs(od),
                  owned_index(ow_idx),
                  dofs(d),
                  index(idx),
                  owned_thread_count(gc),
                  local_thread_count(lc),
                  proc(p) {}

            MARS_INLINE_FUNCTION
            void lg_insert(const Mesh mesh,
                          const Integer i,
                          const KeyType dof_sfc,
                          const Integer owner_proc,
                          std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                auto owned_size = mesh.get_chunk_size();
                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    auto scan_idx = local_thread_count(owned_size + i - 1);
                    dofs(scan_idx + index++) = dof_sfc;
                }
            }

            MARS_INLINE_FUNCTION
            void lg_insert(const Mesh mesh,
                              const Integer i,
                              const KeyType dof_sfc,
                              const Integer owner_proc,
                              std::false_type) const {

                if (proc >= owner_proc) {
                    //since owned_thread_count is an inclusive scan the address starts at 0
                    auto owned_scan_idx = i == 0 ? 0 : owned_thread_count(i - 1);
                    owned_dofs(owned_scan_idx + owned_index++) = dof_sfc;
                }
                // since owned_thread_count is an inclusive scan the address starts at 0
                auto scan_idx = i == 0 ? 0 : local_thread_count(i - 1);
                dofs(scan_idx + index++) = dof_sfc;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const KeyType dof_sfc,
                            const Integer owner_proc) const {
                lg_insert(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <typename T, bool Ghost>
        struct LocalOwnedDofsFill<T, Ghost, DofLabel::lVolume> {
            ViewVectorType<KeyType> owned_dofs;
            Integer &owned_index;
            ViewVectorType<KeyType> dofs;
            Integer &index;
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> local_thread_count;

            MARS_INLINE_FUNCTION
            LocalOwnedDofsFill(ViewVectorType<KeyType> od,
                               Integer &ow_idx,
                               ViewVectorType<KeyType> d,
                               Integer &idx,
                               ViewVectorType<T> gc,
                               ViewVectorType<T> lc)
                : owned_dofs(od),
                  owned_index(ow_idx),
                  dofs(d),
                  index(idx),
                  owned_thread_count(gc),
                  local_thread_count(lc) {}

            MARS_INLINE_FUNCTION
            void volume_insert(const Mesh mesh,
                                  const Integer i,
                                  const KeyType dof_sfc,
                                  std::true_type) const {
                auto owned_size = mesh.get_chunk_size();
                auto scan_idx = local_thread_count(owned_size + i - 1);
                dofs(scan_idx + index++) = dof_sfc;
            }

            MARS_INLINE_FUNCTION
            void volume_insert(const Mesh mesh,
                                  const Integer i,
                                  const KeyType dof_sfc,
                                  std::false_type) const {
                //since owned_thread_count is an inclusive scan the address starts at 0
                auto owned_scan_idx = i == 0 ? 0 : owned_thread_count(i - 1);
                owned_dofs(owned_scan_idx + owned_index++) = dof_sfc;
                // since owned_thread_count is an inclusive scan the address starts at 0
                auto scan_idx = i == 0 ? 0 : local_thread_count(i - 1);
                dofs(scan_idx + index++) = dof_sfc;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const KeyType dof_sfc) const {
                volume_insert(mesh, i, dof_sfc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost, typename T>
        struct InsertLocalGlobalDofs {
            MARS_INLINE_FUNCTION
            InsertLocalGlobalDofs(Mesh m,
                                  ViewVectorType<KeyType> od,
                                  ViewVectorType<KeyType> d,
                                  ViewVectorType<T> gc,
                                  ViewVectorType<T> lc)
                : mesh(m), owned_dofs(od), dofs(d), owned_thread_count(gc), thread_count(lc) {}

            Mesh mesh;
            ViewVectorType<KeyType> owned_dofs;
            ViewVectorType<KeyType> dofs;
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> thread_count;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const KeyType sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();

                Integer owned_index = 0;
                Integer index = 0;

                auto cp = LocalOwnedDofsFill<T, Ghost, DofLabel::lCorner>(
                    owned_dofs, owned_index, dofs, index, owned_thread_count, thread_count, proc);

                auto fp = LocalOwnedDofsFill<T, Ghost, DofLabel::lFace>(
                    owned_dofs, owned_index, dofs, index, owned_thread_count, thread_count, proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalOwnedDofsFill<T, Ghost, DofLabel::lEdge>(
                    owned_dofs, owned_index, dofs, index, owned_thread_count, thread_count, proc);

                auto vp = LocalOwnedDofsFill<T, Ghost, DofLabel::lVolume>(
                    owned_dofs, owned_index, dofs, index, owned_thread_count, thread_count);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        // generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <typename T, bool Ghost, Integer Label>
        struct LocalOwnedThreadCount {
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> local_thread_count;
            ViewVectorType<bool> nbh_proc_predicate_recv;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalOwnedThreadCount(ViewVectorType<T> gc,
                                  ViewVectorType<T> lc,
                                  ViewVectorType<bool> nppr,
                                  Integer p)
                : owned_thread_count(gc),
                  local_thread_count(lc),
                  nbh_proc_predicate_recv(nppr),
                  proc(p) {}

            MARS_INLINE_FUNCTION
            void lg_count(const Mesh mesh,
                          const Integer i,
                          const Integer owner_proc,
                          std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                auto owned_size = mesh.get_chunk_size();
                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    Kokkos::atomic_exchange(&nbh_proc_predicate_recv(elem_sfc_proc), true);
                    local_thread_count(owned_size + i)++;
                }
            }

            MARS_INLINE_FUNCTION
            void lg_count(const Mesh mesh,
                              const Integer i,
                              const Integer owner_proc,
                              std::false_type) const {

                if (proc >= owner_proc) {
                    owned_thread_count(i)++;
                }
                local_thread_count(i)++;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const KeyType dof_sfc, const Integer owner_proc) const {
                lg_count(mesh, i, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <typename T, bool Ghost>
        struct LocalOwnedThreadCount<T, Ghost, DofLabel::lVolume> {
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> local_thread_count;

            MARS_INLINE_FUNCTION
            LocalOwnedThreadCount(ViewVectorType<T> gc, ViewVectorType<T> lc)
                : owned_thread_count(gc), local_thread_count(lc) {}

            MARS_INLINE_FUNCTION
            void volume_predicate(const Mesh mesh, const Integer i, std::true_type) const {
                auto owned_size = mesh.get_chunk_size();
                local_thread_count(owned_size + i)++;
            }

            MARS_INLINE_FUNCTION
            void volume_predicate(const Mesh mesh, const Integer i, std::false_type) const {
                owned_thread_count(i)++;
                local_thread_count(i)++;
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const KeyType dofs_sfc) const {
                volume_predicate(mesh, i, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost, typename T>
        struct CountLocalGlobalDofs {
            MARS_INLINE_FUNCTION
            CountLocalGlobalDofs(Mesh m,
                                 ViewVectorType<T> gc,
                                 ViewVectorType<T> lc,
                                 ViewVectorType<bool> npbr)
                : mesh(m),
                  owned_thread_count(gc),
                  thread_count(lc),
                  nbh_proc_predicate_recv(npbr) {}

            Mesh mesh;
            ViewVectorType<T> owned_thread_count;
            ViewVectorType<T> thread_count;
            ViewVectorType<bool> nbh_proc_predicate_recv;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const KeyType sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();
                auto cp = LocalOwnedThreadCount<T, Ghost, DofLabel::lCorner>(
                    owned_thread_count, thread_count, nbh_proc_predicate_recv, proc);

                auto fp = LocalOwnedThreadCount<T, Ghost, DofLabel::lFace>(
                    owned_thread_count, thread_count, nbh_proc_predicate_recv, proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalOwnedThreadCount<T, Ghost, DofLabel::lEdge>(
                    owned_thread_count, thread_count, nbh_proc_predicate_recv, proc);

                auto vp = LocalOwnedThreadCount<T, Ghost, DofLabel::lVolume>(owned_thread_count, thread_count);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        void print_dofs(const int rank = 0) {
            auto handler = device_view();

            printf("owned and local size: %li, %li\n",
                   get_global_dof_enum().get_elem_size(),
                   get_local_dof_enum().get_elem_size());

            SFC<ElemType, SfcKeyType> gdof = get_global_dof_enum();
            Kokkos::parallel_for(
                "for", gdof.get_elem_size(), MARS_LAMBDA(const int i) {
                    const auto sfc_elem = handler.owned_to_sfc(i);
                    const Dof global_dof = handler.sfc_to_global_dof(sfc_elem);

                    Octant o = handler.get_octant_from_sfc(gdof.get_sfc(i));

                    printf("i: %d global sfc: %ld gdof: %ld octant: [%ld, %ld, %ld] -  rank: %i\n",
                           i,
                           gdof.get_sfc(i),
                           sfc_elem,
                           o.x,
                           o.y,
                           o.z,
                           global_dof.get_proc());
                });

        }

        template <typename KeyType>
        void process_full_segment(const typename KeyType::ValueType start, const typename KeyType::ValueType end) {
            using ValueType = typename KeyType::ValueType;
            static_assert(std::is_same<KeyType, HilbertKey<ValueType>>::value, "Key type must be HilbertKey");

            const auto rank_elements_size = end - start;
            auto rank_elements = ViewVectorType<ValueType>("rank_elements", rank_elements_size);

            printf("Hilbert DOF chunk_size:%li,  start: %li end: %li\n", rank_elements_size, start, end);

            Kokkos::parallel_for(
                "compactSegment",
                Kokkos::RangePolicy<Kokkos::IndexType<ValueType>>(start, end),
                KOKKOS_LAMBDA(const ValueType i) {
                    auto index = i - start;
                    rank_elements(index) = i;
                });

            mesh.set_view_sfc(rank_elements);
            mesh.set_chunk_size(rank_elements_size);
        }

        void build_local_global_dofs(const context &context, ViewVectorType<bool> &nbhr) {
            const int rank_size = num_ranks(context);

            Integer xDim = get_mesh().get_XDim() * degree;
            Integer yDim = get_mesh().get_YDim() * degree;
            Integer zDim = get_mesh().get_ZDim() * degree;

            global_dof_enum = SFC<simplex_type::ElemType, SfcKeyType>(xDim, yDim, zDim);
            local_dof_enum = SFC<simplex_type::ElemType, SfcKeyType>(xDim, yDim, zDim);

            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_ghost_size();

            // Both curves enumerate the same local/ghost DOF union. Sorting the
            // candidates avoids allocating predicates over the global SFC range.
            ViewVectorType<size_t> owned_thread_count("thread_count", size);
            ViewVectorType<size_t> local_thread_count("thread_count", size + ghost_size);
            //count the number of local and global dofs
            // Iterate through the local sfc and enumerate
            Kokkos::parallel_for(
                "lg_predicate",
                size,
                CountLocalGlobalDofs<false, size_t>(get_mesh(), owned_thread_count, local_thread_count, nbhr));
            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for(
                "lg_predicate_from_ghost",
                ghost_size,
                CountLocalGlobalDofs<true, size_t>(get_mesh(), owned_thread_count, local_thread_count, nbhr));

            inclusive_scan(0, owned_thread_count.size(), owned_thread_count);
            inclusive_scan(0, local_thread_count.size(), local_thread_count);
            //get the last value of the thread count
            const Integer owned_size = get_last_value(owned_thread_count);
            const Integer local_size = get_last_value(local_thread_count);

            //reserve the space for the dofs
            ViewVectorType<KeyType> owned_dofs("owned_dofs", owned_size);
            ViewVectorType<KeyType> dofs("dofs", local_size);

            // Insert the local and global dofs
            // Iterate through the local sfc and enumerate
            Kokkos::parallel_for(
                "lg_predicate",
                size,
                InsertLocalGlobalDofs<false, size_t>(get_mesh(), owned_dofs, dofs, owned_thread_count, local_thread_count));
            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for(
                "lg_predicate_from_ghost",
                ghost_size,
                InsertLocalGlobalDofs<true, size_t>(get_mesh(), owned_dofs, dofs, owned_thread_count, local_thread_count));

            kokkos_sort(owned_dofs);
            kokkos_sort(dofs);
            auto owned_unique_dofs = kokkos_unique(owned_dofs);
            auto unique_dofs = kokkos_unique(dofs);

            global_dof_enum.set_elements(owned_unique_dofs);
            local_dof_enum.set_elements(unique_dofs);

            const bool owned_map_valid = global_dof_enum.generate_sfc_to_local_map();
            const bool local_map_valid = local_dof_enum.generate_sfc_to_local_map();
            const int maps_valid = owned_map_valid && local_map_valid;
            if (!context->distributed->min(maps_valid))
                throw std::runtime_error("Structured DOF SFC map construction failed");

            ViewVectorType<Integer> global_dof_size_per_rank("global_dof_size_per_rank", rank_size);
            global_dof_offset = ViewVectorType<Integer>("global_dof_offset", rank_size + 1);

            context->distributed->gather_all_view(global_dof_enum.get_elem_size(), global_dof_size_per_rank);
            incl_excl_scan(0, rank_size, global_dof_size_per_rank, global_dof_offset);
        }

        void exchange_ghost_counts(const context &context,
                                   ViewVectorType<Integer> &scan_boundary,
                                   ViewVectorType<Integer> &sender_ranks,
                                   ViewVectorType<Integer> &recver_ranks,
                                   std::vector<Integer> &send_count,
                                   std::vector<Integer> &receive_count)

        {
            using namespace Kokkos;

            Kokkos::Timer timer;

            int proc_num = rank(context);
            int size = num_ranks(context);

            auto scan_snd = create_mirror_view(scan_boundary);
            Kokkos::deep_copy(scan_snd, scan_boundary);

            auto s_rank_mirror = create_mirror_view(sender_ranks);
            Kokkos::deep_copy(s_rank_mirror, sender_ranks);
            const Integer nbh_rank_size = sender_ranks.extent(0);

            for (int i = 0; i < nbh_rank_size; ++i) {
                Integer count = scan_snd(i + 1) - scan_snd(i);
                if (count > 0) {
                    send_count[s_rank_mirror(i)] = count;

                }
            }

            auto r_rank_mirror = create_mirror_view(recver_ranks);
            Kokkos::deep_copy(r_rank_mirror, recver_ranks);
            const Integer nbh_rank_recv_size = recver_ranks.extent(0);

            for (int i = 0; i < nbh_rank_recv_size; ++i) {
                receive_count[r_rank_mirror(i)] = 1;
            }

            context->distributed->i_send_recv_vec(send_count, receive_count);

        }

        void exchange_ghost_dofs(const context &context,
                                 ViewVectorType<Integer> boundary_lsfc_index,
                                 ViewVectorType<Integer> &ghost_dofs_index) {
            int rank_size = num_ranks(context);
            Integer ghost_size = scan_recv_mirror(rank_size);
            std::cout << "DofHandler:MPI send receive for the ghost dofs layer started..." << std::endl;

            ghost_dofs_sfc = ViewVectorType<Integer>("ghost_dofs", ghost_size);
            ghost_dofs_index = ViewVectorType<Integer>("ghost_dofs_index", ghost_size);
            // do it again to have all the process range to make it fit the i_send_recv_view

            context->distributed->i_send_recv_view(
                ghost_dofs_sfc, scan_recv_mirror.data(), boundary_dofs_sfc, scan_send_mirror.data());

            context->distributed->i_send_recv_view(
                ghost_dofs_index, scan_recv_mirror.data(), boundary_lsfc_index, scan_send_mirror.data());
            std::cout << "DofHandler:MPI send receive for the ghost dofs layer done." << std::endl;
        }

        template <class UMAP>
        static MARS_INLINE_FUNCTION Integer get_value_in_map(const UMAP map, const KeyType sfc) {
            const auto it = map.find(sfc);
            if (map.valid_at(it)) {
                return map.value_at(it);
            } else {
                return INVALID_INDEX;
            }
        }

        template <bool Ghost>
        struct FaceOrientDof {
            ViewVectorType<Integer> dir;
            UnorderedMap<KeyType, Integer> sfc_to_local;
            Integer proc;

            MARS_INLINE_FUNCTION
            FaceOrientDof(ViewVectorType<Integer> d, UnorderedMap<KeyType, Integer> l, Integer p)
                : dir(d), sfc_to_local(l), proc(p) {}

            MARS_INLINE_FUNCTION
            Integer get_local_from_sfc(const KeyType sfc) const { return get_value_in_map(sfc_to_local, sfc); }

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh mesh,
                                 const Integer i,
                                 const KeyType sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    Integer index = get_local_from_sfc(sfc);
                    if (index > INVALID_INDEX) Kokkos::atomic_exchange(&dir(index), d);
                }
            }

            MARS_INLINE_FUNCTION
            void face_orient_dof(const Mesh mesh,
                                 const Integer i,
                                 const KeyType sfc,
                                 const Integer owner_proc,
                                 const Integer d,
                                 std::false_type) const {
                Integer index = get_local_from_sfc(sfc);
                if (index > INVALID_INDEX) Kokkos::atomic_exchange(&dir(index), d);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh,
                            const Integer i,
                            const KeyType dof_sfc,
                            const Integer owner_proc,
                            const Integer dir) const {
                face_orient_dof(mesh, i, dof_sfc, owner_proc, dir, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct OrientDofs {
            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const auto sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();

                if (face_dofs > 0) {
                    FaceOrientDof<Ghost> fp = FaceOrientDof<Ghost>(face_dir, stl_map, proc);
                    // templatized with 0 values means no filter
                    face_dir_iterate<0>(sfc, mesh, i, fp);
                }
            }

            MARS_INLINE_FUNCTION
            OrientDofs(Mesh m, ViewVectorType<Integer> sd, UnorderedMap<KeyType, Integer> sl)
                : mesh(m), face_dir(sd), stl_map(sl) {}

            Mesh mesh;
            ViewVectorType<Integer> face_dir;
            UnorderedMap<KeyType, Integer> stl_map;
        };

        void build_local_orientation() {
            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_view_ghost().extent(0);

            const Integer local_size = get_local_dof_enum().get_elem_size();

            local_dof_enum.init_element_orientations(local_size);

            Kokkos::parallel_for("orient_dofs",
                                 size,
                                 OrientDofs<false>(get_mesh(),
                                                   get_local_dof_enum().get_view_element_orientations(),
                                                   get_local_dof_enum().get_sfc_to_local_map()));

            Kokkos::parallel_for("orient_ghost_dofs",
                                 ghost_size,
                                 OrientDofs<true>(get_mesh(),
                                                  get_local_dof_enum().get_view_element_orientations(),
                                                  get_local_dof_enum().get_sfc_to_local_map()));
        }

        // generic local predicate functor to be used for edge, face and corners. The volume is specialized.
        template <bool Ghost, Integer Label>
        struct LocalGlobalLabel {
            SFC<ElemType, SfcKeyType> local_enum;
            SFC<ElemType, SfcKeyType> global_enum;
            Integer proc;

            MARS_INLINE_FUNCTION
            LocalGlobalLabel(SFC<ElemType, SfcKeyType> le, SFC<ElemType, SfcKeyType> ge, Integer p)
                : local_enum(le), global_enum(ge), proc(p) {}

            MARS_INLINE_FUNCTION
            void set_label_from_sfc(const Mesh mesh,
                                    const Integer i,
                                    const KeyType sfc,
                                    const Integer owner_proc,
                                    std::true_type) const {
                Integer elem_sfc_proc =
                    find_owner_processor(mesh.get_view_gp(), mesh.get_ghost_sfc(i), 2, mesh.get_proc());

                // if the ghost elem owns the dof then he is able to send it.
                if (elem_sfc_proc >= owner_proc) {
                    auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                    if (index > INVALID_INDEX) local_enum.set_label(index, Label);
                }
            }

            MARS_INLINE_FUNCTION
            void set_label_from_sfc(const Mesh mesh,
                                    const Integer i,
                                    const KeyType sfc,
                                    const Integer owner_proc,
                                    std::false_type) const {

                if (proc >= owner_proc) {
                    auto gindex = get_value_in_map(global_enum.get_sfc_to_local_map(), sfc);
                    if (gindex > INVALID_INDEX)
                        global_enum.set_label(gindex, Label);
                    else
                        printf("Setting global label: %li, %li, %li, %li, %li\n", Label, gindex, i, sfc, owner_proc);
                }
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                if (index > INVALID_INDEX)
                    local_enum.set_label(index, Label);
                else
                    printf("Setting global label: %li, %li, %li, %li, %li\n", Label, index, i, sfc, owner_proc);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const KeyType dof_sfc, const Integer owner_proc) const {
                set_label_from_sfc(mesh, i, dof_sfc, owner_proc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct LocalGlobalLabel<Ghost, DofLabel::lVolume> {
            SFC<ElemType, SfcKeyType> local_enum;
            SFC<ElemType, SfcKeyType> global_enum;

            MARS_INLINE_FUNCTION
            LocalGlobalLabel(SFC<ElemType, SfcKeyType> le, SFC<ElemType, SfcKeyType> ge)
                : local_enum(le), global_enum(ge) {}

            MARS_INLINE_FUNCTION
            void volume_label(const KeyType sfc, std::true_type) const {
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                if (index > INVALID_INDEX) local_enum.set_label(index, DofLabel::lVolume);
            }

            MARS_INLINE_FUNCTION
            void volume_label(const KeyType sfc, std::false_type) const {
                auto gindex = get_value_in_map(global_enum.get_sfc_to_local_map(), sfc);
                if (gindex > INVALID_INDEX) global_enum.set_label(gindex, DofLabel::lVolume);
                auto index = get_value_in_map(local_enum.get_sfc_to_local_map(), sfc);
                if (index > INVALID_INDEX) local_enum.set_label(index, DofLabel::lVolume);
            }

            MARS_INLINE_FUNCTION
            void operator()(const Mesh mesh, const Integer i, const Integer enc_oc) const {
                volume_label(enc_oc, std::integral_constant<bool, Ghost>{});
            }
        };

        template <bool Ghost>
        struct BuildLocalGlobalLabel {
            MARS_INLINE_FUNCTION
            BuildLocalGlobalLabel(Mesh m, SFC<ElemType, SfcKeyType> le, SFC<ElemType, SfcKeyType> ge)
                : mesh(m), local_enum(le), global_enum(ge) {}

            Mesh mesh;
            SFC<ElemType, SfcKeyType> local_enum;
            SFC<ElemType, SfcKeyType> global_enum;

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {
                const auto sfc = get_sfc_ghost_or_local<Ghost>(mesh, i);
                const Integer proc = mesh.get_proc();
                auto cp = LocalGlobalLabel<Ghost, DofLabel::lCorner>(local_enum, global_enum, proc);

                auto fp = LocalGlobalLabel<Ghost, DofLabel::lFace>(local_enum, global_enum, proc);

                // fp logic is the same as the edge one. Just different label.
                auto ep = LocalGlobalLabel<Ghost, DofLabel::lEdge>(local_enum, global_enum, proc);

                auto vp = LocalGlobalLabel<Ghost, DofLabel::lVolume>(local_enum, global_enum);

                elem_dof_iterate(sfc, mesh, i, cp, fp, vp, ep);
            }
        };

        void label_local_global_dofs() {
            const Integer size = get_mesh().get_chunk_size();
            const Integer ghost_size = get_mesh().get_view_ghost().extent(0);

            local_dof_enum.reserve_element_labels(local_dof_enum.get_elem_size());
            global_dof_enum.reserve_element_labels(global_dof_enum.get_elem_size());

            Kokkos::parallel_for("label_local",
                                 size,
                                 BuildLocalGlobalLabel<false>(get_mesh(), get_local_dof_enum(), get_global_dof_enum()));

            // Iterate through ghost sfc and enumerate
            Kokkos::parallel_for("label_ghost",
                                 ghost_size,
                                 BuildLocalGlobalLabel<true>(get_mesh(), get_local_dof_enum(), get_global_dof_enum()));
        }

        virtual void enumerate_dofs() {
            printf("Enumerating the Degrees of Freedom on the Mesh with Degree: %li\n", degree);

            const auto &context = get_context();

            const auto coordinate_limit = structured_coordinate_limit<ElemType, SfcKeyType>();
            const auto maximum_dimension = (coordinate_limit - 1) / degree;
            const auto& current_mesh = get_mesh();
            const bool dimensions_valid = current_mesh.get_chunk_size() > 0 &&
                current_mesh.get_XDim() > 0 && current_mesh.get_YDim() > 0 &&
                static_cast<std::uint64_t>(current_mesh.get_XDim()) <= maximum_dimension &&
                static_cast<std::uint64_t>(current_mesh.get_YDim()) <= maximum_dimension &&
                (ElemType == ElementType::Quad4 || (current_mesh.get_ZDim() > 0 &&
                 static_cast<std::uint64_t>(current_mesh.get_ZDim()) <= maximum_dimension));
            if (!context->distributed->min(static_cast<int>(dimensions_valid))) {
                throw std::invalid_argument("Generate a mesh whose degree-scaled coordinates fit the DOF SFC key.");
            }

            const Integer rank_size = num_ranks(context);
            const int proc_num = rank(context);

            ViewVectorType<bool> nbh_proc_predicate_send("send_to", rank_size);
            ViewVectorType<bool> nbh_proc_predicate_recv("receive_from", rank_size);
            ViewVectorType<Integer> proc_scan_send("nbh_scan_send", rank_size + 1);
            ViewVectorType<Integer> proc_scan_recv("nbh_scan_recv", rank_size + 1);

            build_local_global_dofs(context, nbh_proc_predicate_recv);

            label_local_global_dofs();
            build_local_orientation();

            // A neighbor may receive our DOFs without owning any that we need.
            // Use mesh peers for send candidates; the boundary sets give the actual counts.
            const auto mesh_send_scan = get_mesh().get_view_scan_boundary();
            const auto mesh_recv_scan = get_mesh().get_view_scan_ghost();
            Kokkos::parallel_for("dof_peer_candidates", rank_size, KOKKOS_LAMBDA(const Integer peer) {
                nbh_proc_predicate_send(peer) = mesh_send_scan(peer + 1) > mesh_send_scan(peer);
                if constexpr (volume_dofs > 0) {
                    if (mesh_recv_scan(peer + 1) > mesh_recv_scan(peer)) nbh_proc_predicate_recv(peer) = true;
                }
            });

            incl_excl_scan(0, rank_size, nbh_proc_predicate_send, proc_scan_send);
            incl_excl_scan(0, rank_size, nbh_proc_predicate_recv, proc_scan_recv);

            auto ps = get_last_value(proc_scan_send);
            auto pr = get_last_value(proc_scan_recv);

            ViewVectorType<Integer> sender_ranks("sender_ranks", ps);
            ViewVectorType<Integer> recver_ranks("receiver_ranks", pr);
            compact_scan(nbh_proc_predicate_send, proc_scan_send, sender_ranks);
            compact_scan(nbh_proc_predicate_recv, proc_scan_recv, recver_ranks);

            ViewVectorType<Integer> scan_boundary, boundary_dofs_index;
            build_boundary_dof_sets(scan_boundary, boundary_dofs_index, proc_scan_send, ps);

            Kokkos::fence();

            std::vector<Integer> s_count(rank_size, 0);
            std::vector<Integer> r_count(rank_size, 0);
            exchange_ghost_counts(context, scan_boundary, sender_ranks, recver_ranks, s_count, r_count);

            scan_recv = ViewVectorType<Integer>("scan_recv_", rank_size + 1);
            host_scan_to_device(scan_recv, scan_recv_mirror, r_count);
            scan_send = ViewVectorType<Integer>("scan_send_", rank_size + 1);
            host_scan_to_device(scan_send, scan_send_mirror, s_count);

            ViewVectorType<Integer> ghost_dofs_index;
            exchange_ghost_dofs(context, boundary_dofs_index, ghost_dofs_index);

            build_ghost_local_global_map(context, ghost_dofs_index);
        }

        struct BuildLocalToGlobalGhostMap {
            Mesh mesh;
            UnorderedMap<Integer, Dof> glgm;
            ViewVectorType<Integer> global_dof_offset;
            ViewVectorType<Integer> scan_recv_proc;

            ViewVectorType<Integer> ghost_dofs;
            ViewVectorType<Integer> ghost_dofs_index;

            MARS_INLINE_FUNCTION
            BuildLocalToGlobalGhostMap(Mesh m,
                                       UnorderedMap<Integer, Dof> g,
                                       ViewVectorType<Integer> gdo,
                                       ViewVectorType<Integer> srm,
                                       ViewVectorType<Integer> gd,
                                       ViewVectorType<Integer> gdi)
                : mesh(m),
                  glgm(g),
                  global_dof_offset(gdo),
                  scan_recv_proc(srm),
                  ghost_dofs(gd),
                  ghost_dofs_index(gdi) {}

            MARS_INLINE_FUNCTION
            void operator()(const Integer i) const {

                const auto ghost_sfc = ghost_dofs(i);
                const auto ghost_sfc_lid = ghost_dofs_index(i);

                const Integer owner_proc = find_owner_processor(scan_recv_proc, i, 1, mesh.get_proc());
                const auto gid = ghost_sfc_lid + global_dof_offset(owner_proc);

                // build the ghost dof object and insert into the map
                const auto result = glgm.insert(ghost_sfc, Dof(gid, owner_proc));
                assert(result.success());
            }
        };

        void build_ghost_local_global_map(const context &context, const ViewVectorType<Integer> &ghost_dofs_index) {
            int rank_size = num_ranks(context);
            Integer size = scan_recv_mirror(rank_size);

            ghost_local_to_global_map = UnorderedMap<Integer, Dof>(size);

            Kokkos::parallel_for("BuildLocalGlobalmap",
                                 size,
                                 BuildLocalToGlobalGhostMap(get_mesh(),
                                                            ghost_local_to_global_map,
                                                            global_dof_offset,
                                                            scan_recv,
                                                            ghost_dofs_sfc,
                                                            ghost_dofs_index));

            assert(size == ghost_local_to_global_map.size());
        }

        template <typename F>
        void boundary_owned_dof_iterate(F f, const std::string side = "all") {
            const Integer side_value = map_side_to_value<ElemType>(side);

            auto handler = device_view();
            owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const auto sfc = handler.owned_to_sfc(i);
                if (handler.template boundary_sfc<ElemType>(sfc, side_value)) {
                    f(i, sfc);
                }
            });
        }

        template <typename F>
        void boundary_dof_iterate(F f, const Integer bl) {
            boundary_dof_iterate(f, "all", bl);
        }

        template <typename F>
        void boundary_dof_iterate(F f, const std::string side = "all", const Integer bl = -1) {
            const Integer side_value = map_side_to_value<ElemType>(side);

            auto handler = device_view();
            owned_dof_iterate(MARS_LAMBDA(const Integer i) {
                const auto sfc = handler.owned_to_sfc(i);
                auto comp = handler.get_component(i);
                if (is_separated_block(bl, comp) && handler.template boundary_sfc<ElemType>(sfc, side_value)) {
                    const Integer local = handler.sfc_to_local(sfc);
                    f(handler.compute_block_index(local, comp));
                }
            });
        }

        Integer get_global_base_dof_size() const {
            const Integer rank_size = num_ranks(get_context());

            auto ss = subview(global_dof_offset, rank_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            return h_ss();
        }

        template <Integer B = Block_>
        Integer get_global_dof_size() const {
            return this->template get_block<B>() * get_global_base_dof_size();
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_recv_mirror() const { return scan_recv_mirror; }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer>::HostMirror &get_view_scan_send_mirror() const { return scan_send_mirror; }
    private:
        const context& ctx;
        using DeviceView::ghost_dofs_sfc;
        using DeviceView::scan_recv;
        using DeviceView::block;
        using DeviceView::mesh;
        using DeviceView::local_dof_enum;
        using DeviceView::global_dof_enum;
        using DeviceView::global_dof_offset;
        using DeviceView::ghost_local_to_global_map;
        using DeviceView::boundary_dofs_sfc;
        using DeviceView::scan_send;

        ViewVectorType<Integer>::HostMirror scan_send_mirror;
        ViewVectorType<Integer>::HostMirror scan_recv_mirror;};

}  // namespace mars

#endif
#endif

#endif
