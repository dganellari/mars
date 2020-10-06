#ifndef MARS_EDGE_HPP
#define MARS_EDGE_HPP

#include "mars_base.hpp"
#include "mars_fwd.hpp"
#include "mars_stream.hpp"

#include <algorithm>
#include <array>
#include <initializer_list>
#include <vector>

namespace mars {

template <Integer N, class Implementation_> class Side {
public:
  static_assert(N > 0, "N cannot be zero");

  std::array<Integer, N> nodes;

  virtual ~Side() {}

  Integer &operator[](const Integer i) {
    assert(i >= 0);
    assert(i < N);
    return nodes[i];
  }

  const Integer &operator[](const Integer i) const {
    assert(i >= 0);
    assert(i < N);
    return nodes[i];
  }

  inline bool has_node(const Integer v) const {
    for (auto n : nodes) {
      if (v == n)
        return true;
    }

    return false;
  }

  Side() { std::fill(nodes.begin(), nodes.end(), INVALID_INDEX); }

  Side(const std::array<Integer, N> &in) : nodes(in) { fix_ordering(); }

  inline bool is_valid() const {
    for (auto n : nodes) {
      if (n == INVALID_INDEX)
        return false;
    }

    for (Integer i = 1; i < N; ++i) {
      if (nodes[i - 1] >= nodes[i])
        return false;
    }

    return true;
  }

  void fix_ordering() { std::sort(std::begin(nodes), std::end(nodes)); }

  Side(const std::vector<Integer> &in) {
    assert(N == in.size());

    std::copy(std::begin(in), std::end(in), std::begin(nodes));
    std::sort(std::begin(nodes), std::end(nodes));
  }

  Side(std::initializer_list<Integer> in) {
    assert(N == in.size());

    std::copy(std::begin(in), std::end(in), std::begin(nodes));
    std::sort(std::begin(nodes), std::end(nodes));
  }

  inline bool operator==(const Side &other) const {
    for (Integer i = 0; i < N; ++i) {
      if (nodes[i] != other.nodes[i])
        return false;
    }

    return true;
  }

  inline bool operator!=(const Side &other) const {
    return !((*this) == other);
  }

  inline bool operator<(const Side &other) const {
    for (Integer i = 0; i < N - 1; ++i) {
      if (nodes[i] < other.nodes[i]) {
        return true;
      }

      if (nodes[i] > other.nodes[i]) {
        return false;
      }
    }

    return nodes[N - 1] < other.nodes[N - 1];
  }

  void describe(std::ostream &os) const {
    os << "(";
    for (Integer i = 0; i < N - 1; ++i) {
      os << nodes[i] << ",";
    }

    os << nodes[N - 1] << ")";
  }
};

class Edge : public Side<2> {
public:
  Edge() : Side<2>() {}
  Edge(const Integer a_node, const Integer another_node)
      : Side({a_node, another_node}) {}
};

template <Integer N> void write(const Side<N> &side, std::ostream &os) {
  write(&side.nodes[0], side.nodes.size(), os);
}

template <Integer N> void read(Side<N> &side, std::istream &is) {
  read(&side.nodes[0], side.nodes.size(), is);
}

inline void write(const Edge &edge, std::ostream &os) {
  write(static_cast<const Side<2> &>(edge), os);
}

inline void read(Edge &edge, std::istream &is) {
  read(static_cast<Side<2> &>(edge), is);
}

} // namespace mars

#endif // MARS_EDGE_HPP
