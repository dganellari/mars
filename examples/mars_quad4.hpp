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
}; // namespace mars

} // namespace mars

#endif // MARS_QUAD4_HPP
