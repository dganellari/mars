#ifndef mars_boundary_conditions_hpp
#define mars_boundary_conditions_hpp

#include "mars_base.hpp"
#include "mars_globals.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "mars_laplace_ex.hpp"

namespace mars {

    template <class Mesh>
    class BoundaryConditions {
    public:
        static const int Dim = Mesh::Dim;

        BoundaryConditions(Mesh &mesh) : mesh_(mesh) {}

        template <typename F>
        void apply(ViewVectorType<Real> &v, F fun) {
            ViewMatrixType<Real> points = mesh_.get_view_points();

            Kokkos::parallel_for(
                "BoundaryConditions::apply", mesh_.n_nodes(), MARS_LAMBDA(const Integer i) {
                    Real p[Dim];

                    for (int d = 0; d < Dim; ++d) {
                        p[d] = points(i, d);
                    }

                    fun(p, v(i));
                });
        }

    private:
        Mesh &mesh_;
    };

    ////////// Extras ////////////////////////////////////////

    template <int Dim>
    MARS_INLINE_FUNCTION bool is_boundary_of_unit_cube(const Real *p) {
        bool ret = false;
        for (int d = 0; d < Dim; ++d) {
            if (p[d] <= 1e-14 || p[d] >= 1 - 1e-14) {
                ret = true;
                break;
            }
        }

        return ret;
    }

    template <class Mesh>
    class ZeroDirchletOnUnitCube {
    public:
        static const int Dim = Mesh::Dim;

        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = 0.0;
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<Dim>(p); }
    };

    class Example1Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex1_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example1RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex1_laplacian(p); }
    };

    class Example2Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex2_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example2Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_exact(p); }
    };

    class Example2RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex2_laplacian(p); }
    };

    ////////////////////////////////

    class Example3Dirichlet {
    public:
        MARS_INLINE_FUNCTION void operator()(const Real *p, Real &val) const {
            if (is_boundary(p)) {
                val = ex3_exact(p);
            }
        }

        MARS_INLINE_FUNCTION static bool is_boundary(const Real *p) { return is_boundary_of_unit_cube<2>(p); }
    };

    class Example3Analitcal {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_exact(p); }
    };

    class Example3RHS {
    public:
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const { return ex3_laplacian(p); }
    };

    template <class Mesh>
    class One {
    public:
        static const int Dim = Mesh::Dim;
        MARS_INLINE_FUNCTION Real operator()(const Real *) const { return 1.0; }
    };

    template <class Mesh>
    class Norm2Squared {
    public:
        static const int Dim = Mesh::Dim;
        MARS_INLINE_FUNCTION Real operator()(const Real *p) const {
            Real ret = 0.0;

            for (int i = 0; i < Dim; ++i) {
                const Real x = p[i];
                ret += x * x;
            }

            return ret;
        }
    };

}  // namespace mars

#endif
#endif  // mars_boundary_conditions_hpp