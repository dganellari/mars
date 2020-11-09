#ifndef GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_fd_data_management.hpp"

namespace mars {

    template <Integer Dim>
    using StagStencil = Stencil<Dim, 2>;

    template <class Mesh, typename FaceStencil = StagStencil<Mesh::Dim>, typename... T>
    class StagDM : public FDDM<Mesh, 2, T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer Degree = 2;
        static constexpr Integer VWidth = 1;
        /* static constexpr Integer FWidth = 2; */

        using SuperFDDM = FDDM<Mesh, Degree, T...>;
        using SuperDM = DM<Mesh, Degree, T...>;

        using VolumeStencil = Stencil<SuperFDDM::Dim, VWidth>;

        MARS_INLINE_FUNCTION
        StagDM(Mesh *mesh, const context &c) : SuperFDDM(mesh, c) {}

        virtual void enumerate_dofs(const context &context) override {
            SuperFDDM::enumerate_dofs(context);
            /* build_stencils(); */
        }

        void build_pressure_stencil() { volume_stencil = mars::build_volume_stencil<VolumeStencil>(*this); }

        //TODO: remove the orientation dependency on the corner extension then add the template here.
        /* template<typename Orient = false> */
        void build_stokes_stencil() { face_stencil = mars::build_face_stencil<FaceStencil>(*this); }

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
