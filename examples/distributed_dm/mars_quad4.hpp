#ifndef MARS_QUAD4_HPP
#define MARS_QUAD4_HPP

#include "mars_globals.hpp"

#include "mars_distributed_utils.hpp"
#include "mars_fe_simplex.hpp"

namespace mars {

    template <typename T_, int PhysicalDim_ = 2>
    class FEQuad4 final {
    public:
        static const int Order = 1;
        static const int Dim = 2;
        static const int PhysicalDim = PhysicalDim_;
        static const int NNodes = 4;

        using T = T_;

        using Vector = T[Dim];
        using CoVector = T[PhysicalDim];

        using Point = T[Dim];
        using CoPoint = T[PhysicalDim];

        class Fun final {
        public:
            MARS_INLINE_FUNCTION static auto f0(const Point &p) -> T { return (1.0 - p[0]) * (1.0 - p[1]); }
            MARS_INLINE_FUNCTION static auto f1(const Point &p) -> T { return p[0] * (1.0 - p[1]); }
            MARS_INLINE_FUNCTION static auto f2(const Point &p) -> T { return p[0] * p[1]; }
            MARS_INLINE_FUNCTION static auto f3(const Point &p) -> T { return (1.0 - p[0]) * p[1]; }

            MARS_INLINE_FUNCTION static Real f(const int i, Vector &p) {
                switch (i) {
                    case 0: {
                        return (1.0 - p[0]) * (1.0 - p[1]);
                    }
                    case 1: {
                        return p[0] * (1.0 - p[1]);
                    }
                    case 2: {
                        return p[0] * p[1];
                    }
                    case 3: {
                        return (1.0 - p[0]) * p[1];
                    }
                }
            }
        };

        class Grad final {
        public:
            // f:= (1.0 - p[0]) * (1.0 - p[1])
            MARS_INLINE_FUNCTION static void f0(const T *p, Vector &g) {
                g[0] = p[1] - 1.0;
                g[1] = p[0] - 1.0;
            }

            // f:= p[0] * (1.0 - p[1])
            MARS_INLINE_FUNCTION static void f1(const T *p, Vector &g) {
                g[0] = 1.0 - p[1];
                g[1] = -p[0];
            }

            // f := p[0] * p[1]
            MARS_INLINE_FUNCTION static void f2(const T *p, Vector &g) {
                g[0] = p[1];
                g[1] = p[0];
            }

            // f := (1.0 - p[0]) * p[1]
            MARS_INLINE_FUNCTION static void f3(const T *p, Vector &g) {
                g[0] = -p[1];
                g[1] = 1.0 - p[0];
            }

            MARS_INLINE_FUNCTION static void affine_f(const int i, const T *J_inv, const T *p, Vector &g) {
                T gi[2];

                switch (i) {
                    case 0: {
                        f0(p, gi);
                        break;
                    }
                    case 1: {
                        f1(p, gi);
                        break;
                    }
                    case 2: {
                        f2(p, gi);
                        break;
                    }
                    case 3: {
                        f3(p, gi);
                        break;
                    }
                    default: {
                        break;
                    }
                }

                Algebra<2>::m_t_v_mult(J_inv, gi, g);
            }

            MARS_INLINE_FUNCTION static void ref(const T *q, const T *u, Vector &g) {
                T gi[2];

                f0(q, gi);
                g[0] = u[0] * gi[0];
                g[1] = u[0] * gi[1];

                f1(q, gi);
                g[0] += u[1] * gi[0];
                g[1] += u[1] * gi[1];

                f2(q, gi);
                g[0] += u[2] * gi[0];
                g[1] += u[2] * gi[1];

                f3(q, gi);
                g[0] += u[3] * gi[0];
                g[1] += u[3] * gi[1];
            }
        };

        static constexpr int n_qp = 6;
        static T_ q_weights[6];
        static T_ q_points[6][PhysicalDim_];

        class Quadrature final {
        public:
            static Quadrature make() {
                Quadrature q;
                q.init();
                return q;
            }

            MARS_INLINE_FUNCTION static constexpr int n_points() { return n_qp; }
            MARS_INLINE_FUNCTION static constexpr int dim() { return PhysicalDim_; }

            ViewMatrixTextureC<T_, n_qp, PhysicalDim_> q_p;
            ViewVectorTextureC<T_, n_qp> q_w;

            Quadrature() : q_p("q_p"), q_w("q_w") {}

            void init() {
                fill_view_matrix(q_p, q_points);
                fill_view_vector(q_w, q_weights);
            }
        };

    };  // namespace mars

    // a singleton pattern might be a better solution

    template <typename T_, int PhysicalDim_>
    T_ FEQuad4<T_, PhysicalDim_>::q_points[6][PhysicalDim_] = {
        {0.5, 0.5},
        {0.98304589153964795245728880523899, 0.5},
        {0.72780186391809642112479237299488, 0.074042673347699754349082179816666},
        {0.72780186391809642112479237299488, 0.92595732665230024565091782018333},
        {0.13418502421343273531598225407969, 0.18454360551162298687829339850317},
        {0.13418502421343273531598225407969, 0.81545639448837701312170660149683}};

    template <typename T_, int PhysicalDim_>
    T_ FEQuad4<T_, PhysicalDim_>::q_weights[6] = {0.28571428571428571428571428571428,
                                                  0.10989010989010989010989010989011,
                                                  0.14151805175188302631601261486295,
                                                  0.14151805175188302631601261486295,
                                                  0.16067975044591917148618518733485,
                                                  0.16067975044591917148618518733485};
    /*
    template<typename T_, int PhysicalDim_>
    ViewMatrixTextureC<T_, 6, PhysicalDim_> FEQuad4<T_, PhysicalDim_>::Quadrature::q_p;

    template<typename T_, int PhysicalDim_>
    ViewVectorTextureC<T_, 6> FEQuad4<T_, PhysicalDim_>::Quadrature::q_w; */

}  // namespace mars

#endif  // MARS_QUAD4_HPP
