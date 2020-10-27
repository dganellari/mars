#ifndef GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SFDDM_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_fd_data_management.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <class Mesh, Integer Complete = true, typename... T>
    class StagDM : public FDDM<Mesh, 2, T...> {
    public:
        /* using UD = UserData<Mesh, double>; */
        using UD = UserData<Mesh>;
        using simplex_type = typename Mesh::Elem;

        using user_tuple = ViewsTuple<T...>;
        using tuple = std::tuple<T...>;

        template <Integer idx>
        using UserDataType = typename std::tuple_element<idx, tuple>::type;
        constexpr Integer Degree = 2;

        using SuperFDDM= FDDM<Mesh, Degree, T...>;

        MARS_INLINE_FUNCTION
        FDDM(Mesh *mesh, const context &c) : FDDM<Mesh, Degree, T...>(mesh, c) {}

        virtual void enumerate_dofs(const context &context) override {
            SuperFDDM::enumerate_dofs(context);
            build_stencils();
        }

        void build_volume_stencil(){

        }

        void build_face_stencil(){

        }

    private:

        //the pressure stencil for the continuity equation
        Stencil<Dim, Degree, 1, false> volume_stencil;
        //the stokes stencil for the stokes equation
        Stencil<Dim, Degree, 2, Complete> face_stencil;
    };

}  // namespace mars

#endif
#endif

#endif
