#ifndef GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_fd_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <class Mesh, bool Complete = true, typename... T>
    class StagDM : public FDDM<Mesh, 2, T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        static constexpr Integer Degree = 2;
        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer VWidth = 1;
        static constexpr Integer FWidth = 2;
        static constexpr bool VComplete = false;

        using SuperFDDM = FDDM<Mesh, Degree, T...>;
        using SuperDM = DM<Mesh, Degree, T...>;

        using VolumeStencil = Stencil<SuperFDDM::Dim, Degree, VWidth, VComplete>;
        using FaceStencil = Stencil<SuperFDDM::Dim, Degree, FWidth, Complete>;

        MARS_INLINE_FUNCTION
        StagDM(Mesh *mesh, const context &c) : FDDM<Mesh, Degree, T...>(mesh, c) {}

        virtual void enumerate_dofs(const context &context) override {
            SuperFDDM::enumerate_dofs(context);
            build_stencils();
        }

        template <typename F>
        void face_dof_iterate(F f) {
            Kokkos::parallel_for("face_dof_iter", SuperFDDM::get_locally_owned_face_dofs().extent(0), f);
        }

        template <typename F>
        void volume_dof_iterate(F f) {
            Kokkos::parallel_for("volume_dof_iter", SuperFDDM::get_locally_owned_volume_dofs().extent(0), f);
        }


        /*
         * Face numbering on the stencil => ordering in the stencil stencil[1,0,3,2]
                ----3----
                |       |
                0   x   1
                |       |
                ----2---- */
        // building the stencil is the responsibility of the specialized DM.
        void build_volume_stencil() {
            volume_stencil.reserve_stencil(SuperFDDM::get_volume_dof_size());

            auto vstencil = volume_stencil;
            auto dm = *this;
            volume_dof_iterate(MARS_LAMBDA(const Integer i) {
                const Integer localid = dm.get_volume_dof(i);
                SuperFDDM::fill_stencil(dm, vstencil, localid, i);
            });
        }

        void build_face_stencil() {
            face_stencil.reserve_stencil(SuperFDDM::get_face_dof_size());

            auto fstencil = face_stencil;
            auto dm = *this;
            face_dof_iterate(MARS_LAMBDA(const Integer i) {
                const Integer localid = dm.get_face_dof(i);
                SuperFDDM::fill_stencil(dm, fstencil, localid, i);
            });
        }

        void build_stencils() {
            build_volume_stencil();
            build_face_stencil();
        }

        MARS_INLINE_FUNCTION
        const VolumeStencil get_volume_stencil() const {
                return volume_stencil; }

        MARS_INLINE_FUNCTION
        const FaceStencil get_face_stencil() const {
                return face_stencil; }

    private:
        // the pressure stencil for the continuity equation
        VolumeStencil volume_stencil;
        // the stokes stencil for the stokes equation
        FaceStencil face_stencil;
    };

}  // namespace mars

#endif
#endif

#endif
