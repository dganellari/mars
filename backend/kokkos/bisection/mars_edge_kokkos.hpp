#ifndef MARS_EDGE_KOKKOS_HPP
#define MARS_EDGE_KOKKOS_HPP

#include "mars_base.hpp"
// #include "mars_stream.hpp"

#include <algorithm>
#include <array>
#include <initializer_list>
#include <vector>

#include "mars_interfaces.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template <Integer N>
    class Side<N, KokkosImplementation> {
    public:
        static_assert(N > 0, "N cannot be zero");

        TempArray<Integer, N> nodes;

        MARS_INLINE_FUNCTION
        virtual ~Side() {}

        MARS_INLINE_FUNCTION
        Integer &operator[](const Integer i) {
            assert(i >= 0);
            assert(i < N);
            return nodes[i];
        }

        MARS_INLINE_FUNCTION
        const Integer &operator[](const Integer i) const {
            assert(i >= 0);
            assert(i < N);
            return nodes[i];
        }

        MARS_INLINE_FUNCTION
        bool has_node(const Integer v) const {
            for (Integer i = 0; i < N; ++i) {
                if (v == nodes[i]) return true;
            }

            return false;
        }

        MARS_INLINE_FUNCTION Side() { nodes.set(INVALID_INDEX); }

        MARS_INLINE_FUNCTION
        Side(const TempArray<Integer, N> &in) {
            nodes = in;
            fix_ordering();
        }

        MARS_INLINE_FUNCTION
        Side(const Integer a_node, const Integer another_node) {
            static_assert(N == 2, "This constructor can only be used for a 2-side");

            nodes[0] = a_node;
            nodes[1] = another_node;
            fix_ordering();
        }

        MARS_INLINE_FUNCTION
        bool is_valid() const {
            for (Integer j = 0; j < N; ++j) {
                if (nodes[j] == INVALID_INDEX) return false;
            }

            for (Integer i = 1; i < N; ++i) {
                if (nodes[i - 1] >= nodes[i]) return false;
            }

            return true;
        }

        MARS_INLINE_FUNCTION
        void fix_ordering() { quick_sort(nodes, 0, N - 1); }

        /*Side(const std::vector<Integer> &in)
        {
                assert(N == in.size());

                std::copy(std::begin(in), std::end(in), std::begin(nodes));
                std::sort(std::begin(nodes), std::end(nodes));
        }*/

        /*Side(std::initializer_list<Integer> in)
        {
                assert(N == in.size());

                std::copy(std::begin(in), std::end(in), std::begin(nodes));
                std::sort(std::begin(nodes), std::end(nodes));
        }*/

        MARS_INLINE_FUNCTION
        bool operator==(const Side<N, KokkosImplementation> &other) const {
            for (Integer i = 0; i < N; ++i) {
                if (nodes[i] != other.nodes[i]) return false;
            }

            return true;
        }

        MARS_INLINE_FUNCTION
        bool operator!=(const Side<N, KokkosImplementation> &other) const { return !((*this) == other); }

        MARS_INLINE_FUNCTION
        bool operator<(const Side<N, KokkosImplementation> &other) const {
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

        MARS_INLINE_FUNCTION
        void describe() const {
            printf("(");

            for (Integer i = 0; i < N - 1; ++i) {
                printf("%i,", nodes[i]);
            }

            printf("%i)", nodes[N - 1]);
        }
    };

    /*class ParallelEdge : public Side<2,KokkosImplementation> {
    public:
            ParallelEdge() : Side<2, KokkosImplementation>() {}
            ParallelEdge(const Integer a_node, const Integer another_node)
            {
                    nodes[0] = a_node;
                    nodes[1] = another_node;
                    fix_ordering();
            }
    };*/
    /*template<Integer N>
    void write(
        const Side<N> &side,
        std::ostream &os)
    {
        write(&side.nodes[0], side.nodes.size(), os);
    }

    template<Integer N>
    void read(
        Side<N> &side,
        std::istream &is)
    {
        read(&side.nodes[0], side.nodes.size(), is);
    }

    inline void write(
        const Edge &edge,
        std::ostream &os)
    {
        write(static_cast<const Side<2> &>(edge), os);
    }

    inline void read(
        Edge &edge,
        std::istream &is)
    {
        read(static_cast<Side<2> &>(edge), is);
    }*/

}  // namespace mars

#endif  // MARS_EDGE_HPP
