#ifndef MARS_QUAD4_HPP
#define MARS_QUAD4_HPP

#include "mars_globals.hpp"

namespace mars {

template <typename T_, int PhysicalDim_ = 2> class FEQuad4 final {
public:
  static const int Order = 1;
  static const int Dim = 2;
  static const int PhysicalDim = PhysicalDim_;
  static const int NNodes = 4;

  static const int n_qp = 6;
  T_ q_points[n_qp][2] = {
      {0.5, 0.5},
      {0.98304589153964795245728880523899, 0.5},
      {0.72780186391809642112479237299488, 0.074042673347699754349082179816666},
      {0.72780186391809642112479237299488, 0.92595732665230024565091782018333},
      {0.13418502421343273531598225407969, 0.18454360551162298687829339850317},
      {0.13418502421343273531598225407969, 0.81545639448837701312170660149683}};

  T_ q_weights[n_qp] = {
      0.28571428571428571428571428571428, 0.10989010989010989010989010989011,
      0.14151805175188302631601261486295, 0.14151805175188302631601261486295,
      0.16067975044591917148618518733485, 0.16067975044591917148618518733485};


  static ViewMatrixTextureC<T_, 6, 2> q_p;
  fill_view_matrix(q_p, q_points);

  static ViewVectorTextureC<T_, 6> q_w;
  fill_view_vector(q_w, q_weights);

  using T = T_;

  using Vector = T[Dim];
  using CoVector = T[PhysicalDim];

  using Point = T[Dim];
  using CoPoint = T[PhysicalDim];

  class Fun final {
  public:
    MARS_INLINE_FUNCTION static auto f0(const Point &p) -> T {
      return (1.0 - p[0]) * (1.0 - p[1]);
    }
    MARS_INLINE_FUNCTION static auto f1(const Point &p) -> T {
      return p[0] * (1.0 - p[1]);
    }
    MARS_INLINE_FUNCTION static auto f2(const Point &p) -> T {
      return p[0] * p[1];
    }
    MARS_INLINE_FUNCTION static auto f3(const Point &p) -> T {
      return (1.0 - p[0]) * p[1];
    }
  };

  class Grad final {
  public:
    // f:= (1.0 - p[0]) * (1.0 - p[1])
    MARS_INLINE_FUNCTION static void f0(const Point &p, Vector &g) {
      g[0] = p[1] - 1.0;
      g[1] = p[0] - 1.0;
    }

    // f:= p[0] * (1.0 - p[1])
    MARS_INLINE_FUNCTION static void f1(const Point &p, Vector &g) {
      g[0] = 1.0 - p[1];
      g[1] = -p[0];
    }

    // f := p[0] * p[1]
    MARS_INLINE_FUNCTION static void f2(const Point &p, Vector &g) {
      g[0] = p[1];
      g[1] = p[0];
    }

    // f := (1.0 - p[0]) * p[1]
    MARS_INLINE_FUNCTION static void f3(const Point &p, Vector &g) {
      g[0] = -p[1];
      g[1] = 1.0 - p[0];
    }

    MARS_INLINE_FUNCTION static void affine_f(const int i, const T *J_inv,
                                              const Point &p, Vector &g) {
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

      g[0] = J_inv[0] * gi[0] + J_inv[2] * gi[1];
      g[1] = J_inv[1] * gi[0] + J_inv[3] * gi[1];
    }

    MARS_INLINE_FUNCTION static void ref(const Point &q, const T *u,
                                         Vector &g) {
      T gi[2];

      f0(q, gi);
      g[0] += u[0] * gi[0];
      g[1] += u[0] * gi[1];

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
}; // namespace mars

} // namespace mars

#endif // MARS_QUAD4_HPP