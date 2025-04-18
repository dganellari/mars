#include "mars_base.hpp"
#include "mars_env.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>

#if defined (MARS_ENABLE_AMR_BACKEND) && defined (MARS_ENABLE_KOKKOS_KERNELS)
#include <KokkosBlas1_nrm1.hpp>
#include <KokkosBlas1_nrminf.hpp>

#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#include "mars_mesh_kokkos.hpp"
#include "mars_model_test.hpp"
#include "mars_spacetime_ex.hpp"
#include "mars_umesh_laplace.hpp"
#include "mars_umesh_st_heat_equation.hpp"

#include "mars_poisson.hpp"
#include "mars_poisson_operator.hpp"

namespace mars {

    // using MeshQuad4 = mars::ParallelQuad4Mesh;
    // using DMQuad4 = mars::DM<MeshQuad4, 1, Real, Real, Real>;

    // template <>
    // class FEValues<MeshQuad4> : public FEDMValues<DMQuad4> {};

    // template <>
    // class UMeshLaplace<MeshQuad4> : public PoissonOperator<INPUT, OUTPUT, RHSD, DMQuad4> {};

    // template <>
    // class Interpolate<MeshQuad4> : public DMInterpolate<DMQuad4> {};

    class ST1Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_st_exact(p); }
    };

    class ST1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_st_spacetime(p); }
    };

    template <class Mesh>
    class ST1BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex1_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST2Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_st_exact(p); }
    };

    class ST2RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_st_spacetime(p); }
    };

    template <class Mesh>
    class ST2BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex2_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST3Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_exact(p); }
    };

    class ST3RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_st_spacetime(p); }
    };

    template <class Mesh>
    class ST3BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex3_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST3D1Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex4_3D_st_exact(p); }
    };

    class ST3D1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex4_3D_st_spacetime(p); }
    };

    template <class Mesh>
    class ST3D1BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex4_3D_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////
    class ST3D2Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex5_3D_st_exact(p); }
    };

    class ST3D2RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex5_3D_st_spacetime(p); }
    };

    template <class Mesh>
    class ST3D2BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex4_3D_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

    class ST4D1Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex7_4D_st_exact(p); }
    };

    class ST4D1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex7_4D_st_spacetime(p); }
    };

    template <class Mesh>
    class ST4D1BC {
    public:
        /* BC --> zero dirichlet + natural neumann on upper bound */
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex4_3D_st_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) {
            bool ret = false;
            for (int d = 0; d < Dim; ++d) {
                if (p[d] <= 1e-14) {
                    ret = true;
                    break;
                }

                if (d < Dim - 1 && p[d] >= 1 - 1e-14) {
                    ret = true;
                    break;
                }
            }

            return ret;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace mars

int main(int argc, char *argv[]) {
    using namespace mars;

    Env env(argc, argv);

    {

    using ModelTest_t = ModelTest<ParallelQuad4Mesh,
                  UMeshSTHeatEquation<ParallelQuad4Mesh>,
                  ST3BC<ParallelQuad4Mesh>,
                  ST3RHS,
                  ST3Analitcal>;

#ifdef MARS_ENABLE_CXXOPTS
        using namespace cxxopts;
        Options options("./st_example", "Run M.A.R.S. based applications.");

        options.add_options()("d,debug", "Enable debugging")                                     //
            ("l,level", "Number of levels", value<Integer>()->default_value("1"))                //
            ("x,nx", "Number of elements in x direction", value<Integer>()->default_value("6"))  //
            ("y,ny", "Number of elements in y direction", value<Integer>()->default_value("6"))  //
            ("z,nz", "Number of elements in z direction", value<Integer>()->default_value("6"))  //
            ("t,nt", "Number of elements in t direction", value<Integer>()->default_value("6"))  //
            ("a,adaptive", "Adaptivity", value<bool>()->default_value("false"))                  //
            ("o,output", "Enable output", value<bool>()->default_value("true"))                  //
            ("r,refine_level",
             "Number of refinements",
             value<Integer>()->default_value("1"))                                  //
            ("v,verbose", "Verbose output", value<bool>()->default_value("false"))  //
            ("h,help", "Print usage");

        auto args = options.parse(argc, argv);

        // ModelTest<ParallelMesh2, UMeshSTHeatEquation<ParallelMesh2>, ST1BC<ParallelMesh2>, ST1RHS,
        // ST1Analitcal>().run(
        //     args);

        // ModelTest<ParallelMesh2, UMeshLaplace<ParallelMesh2>, Example2Dirichlet, Example2RHS,
        // Example2Analitcal>().run(
        //     args);

        // ModelTest<ParallelQuad4Mesh,
        //           UMeshLaplace<ParallelQuad4Mesh>,
        //           Example2Dirichlet,
        //           Example2RHS,
        //           Example2Analitcal>()
        //     .run(args);

        // ModelTest<ParallelMesh2, UMeshLaplace<ParallelMesh2>, Example3Dirichlet, Example3RHS,
        // Example3Analitcal>().run(
        //     args);

        // ModelTest<ParallelMesh2, UMeshSTHeatEquation<ParallelMesh2>, ST2BC<ParallelMesh2>, ST2RHS,
        // ST2Analitcal>().run(
        //     args);

        // ModelTest<ParallelMesh2, UMeshSTHeatEquation<ParallelMesh2>, ST3BC<ParallelMesh2>, ST3RHS,
        // ST3Analitcal>().run(
        //     args);

        // ModelTest<ParallelMesh3,
        //           UMeshLaplace<ParallelMesh3>,
        //           ZeroDirchletOnUnitCube<ParallelMesh3>,
        //           One<ParallelMesh3>,
        //           One<ParallelMesh3>>()
        //     .run(args);

        // // 3D example 1
        // ModelTest<ParallelMesh3, UMeshSTHeatEquation<ParallelMesh3>, ST3D1BC<ParallelMesh3>, ST3D1RHS,
        // ST3D1Analitcal>()
        //     .run(args);

        // 3D example 2
        // ModelTest<ParallelMesh3, UMeshSTHeatEquation<ParallelMesh3>, ST3D2BC<ParallelMesh3>, ST3D2RHS, ST3D2Analitcal>()
        //     .run(args);

        // 4D example 1
        // ModelTest<ParallelMesh4, UMeshSTHeatEquation<ParallelMesh4>, ST4D1BC<ParallelMesh4>, ST4D1RHS,
        // ST4D1Analitcal>()
        //     .run(args);

        // ModelTest<ParallelMesh4,
        //           UMeshLaplace<ParallelMesh4>,
        //           ZeroDirchletOnUnitCube<ParallelMesh4>,
        //           One<ZeroDirchletOnUnitCube<ParallelMesh4>>,
        //           One<ZeroDirchletOnUnitCube<ParallelMesh4>>>()
        //     .run(args);

        ModelTest_t().run(args);
#else
        ModelTest_t().run();
#endif
    }

    return env.exit_code();
}
#endif
