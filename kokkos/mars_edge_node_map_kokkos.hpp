#ifndef MARS_EDGE_NODE_MAP_KOKKOS_HPP
#define MARS_EDGE_NODE_MAP_KOKKOS_HPP

#include "mars_edge_kokkos.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <vector>

namespace mars {

    template <Integer N>
    class DeviceEdgeNodeMap {
        using Edge = Side<N, KokkosImplementation>;

    public:
        MARS_INLINE_FUNCTION
        static bool update(const Edge& edge,
                           const UnorderedMap<Edge, Integer>& mapping_,
                           const ViewObject<Integer> node_start_index,
                           Integer& midpoint) {
            assert(edge.is_valid());

            auto result = mapping_.insert(edge);  // to make it atomic and insert only for one thread.

            if (result.success()) {
                midpoint = Kokkos::atomic_fetch_add(&node_start_index(0), 1);

                // Kokkos::volatile_store(&mapping_.value_at(result.index()), &midpoint);

                mapping_.value_at(result.index()) = midpoint;

                // Kokkos::memory_fence();

                return true;
            } else if (result.existing()) {
                // wait until the other thread writes the midpoint to map and the value is visible.
                while (midpoint == 0) {
                    // volatile needed  otherwise since value_at is const it assumes that it should not read it again.
                    midpoint = Kokkos::volatile_load(&mapping_.value_at(result.index()));
                }
                return false;
            } else {
                printf("Edge Node Map: Exceeded UnorderedMap capacity\n");
                midpoint = INVALID_INDEX;
                return false;
            }
        }

        MARS_INLINE_FUNCTION
        static Integer get(const Integer an_edge_node,
                           const Integer another_edge_node,
                           const UnorderedMap<Edge, Integer>& mapping_) {
            TempArray<Integer, 2> nodes;
            nodes[0] = an_edge_node;
            nodes[1] = another_edge_node;

            return get(Side<N, KokkosImplementation>(nodes), mapping_);
        }

        MARS_INLINE_FUNCTION
        static Integer get(const Edge& edge, const UnorderedMap<Edge, Integer>& mapping_) {
            auto it = mapping_.find(edge);
            if (!mapping_.valid_at(it)) {
                return INVALID_INDEX;
            }

            return mapping_.value_at(it);
        }

        void describe(UnorderedMap<Edge, Integer>& mapping_) const {
            parallel_for(
                mapping_.capacity(), KOKKOS_LAMBDA(uint32_t i) {
                    if (mapping_.valid_at(i)) {
                        auto key = mapping_.key_at(i);
                        auto value = mapping_.value_at(i);
                        printf("(%li, %li) -> %li\n", key.nodes[0], key.nodes[1], value);
                    }
                });
        }

        void reserve_map(Integer capacity) { mapping_ = UnorderedMap<Edge, Integer>(capacity); }

        void rehash_map(Integer capacity) { mapping_.rehash(capacity); }

        void clear() { mapping_.clear(); }

        Integer size() const { return mapping_.size(); }

        inline bool empty() const { return (mapping_.size() == 0); }

        MARS_INLINE_FUNCTION
        UnorderedMap<Edge, Integer>& mapping() { return mapping_; }

        UnorderedMap<Edge, Integer> mapping_;
    };
}  // namespace mars

#endif  // MARS_EDGE_NODE_MAP_HPP
