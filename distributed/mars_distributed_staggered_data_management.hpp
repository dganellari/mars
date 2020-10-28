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
        static constexpr Integer Degree = 2;

        using SuperFDDM = FDDM<Mesh, Degree, T...>;
        using SuperDM = DM<Mesh, Degree, T...>;

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
            ViewVectorType<Integer> lovd = SuperFDDM::get_locally_owned_volume_dofs();
            volume_stencil.reserve_stencil(lovd.extent(0));

            volume_dof_iterate(MARS_LAMBDA(const Integer i) {
                const Integer localid = lovd(i);
                Dof d = SuperDM::local_to_global_dof(localid);
                printf("localid: i: %li, local: %li, global: %li, proc: %li\n", i, localid, d.get_gid(), d.get_proc());

                const Integer sfc = SuperDM::local_to_sfc(localid);
                Octant oc = get_octant_from_sfc<simplex_type::ElemType>(sfc);

                volume_stencil.stencil(i, 0) = localid;

                Integer face_nr;
                for (int dir = 0; dir < 2; ++dir) {
                    for (int side = 0; side < 2; ++side) {
                        if (side == 0)
                            face_nr = 2 * dir + 1;
                        else
                            face_nr = 2 * dir;

                        const Integer nbh_sfc = SuperDM::get_sfc_face_nbh(oc, face_nr);
                        Integer nbh_id = SuperDM::is_local(nbh_sfc) ? SuperDM::sfc_to_local(nbh_sfc) : -1;
                        Integer index = 2 * dir + side;
                        volume_stencil.stencil(i, index) = nbh_id;
                    }
                }
            });
        }

        void build_face_stencil() {}

        void build_stencils() { build_volume_stencil(); }

    private:
        // the pressure stencil for the continuity equation
        Stencil<SuperFDDM::Dim, Degree, 1, false> volume_stencil;
        // the stokes stencil for the stokes equation
        Stencil<SuperFDDM::Dim, Degree, 2, Complete> face_stencil;
    };

}  // namespace mars

#endif
#endif

#endif
