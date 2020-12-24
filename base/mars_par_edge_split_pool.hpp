#ifndef MARS_PAR_EDGE_SPLIT_POOL_HPP
#define MARS_PAR_EDGE_SPLIT_POOL_HPP

#include <iomanip>
#include <sstream>
#include <string>

#include "mars_communicator.hpp"
#include "mars_edge.hpp"
#include "mars_longest_edge.hpp"

namespace mars {
template <class Mesh, class EdgeSelect>
class Bisection;

template <Integer Dim, Integer ManifoldDim>
class ParMesh;

class ParEdgeSplitPool {
 public:
  using OutputStream = std::ostringstream;
  using InputStream = std::istringstream;
  using BufferObject = std::string;

  template <class ParMesh, class EdgeSelect>
  void build_edge_interface(
      ParMesh &p_mesh,
      Bisection<typename ParMesh::Mesh, EdgeSelect> &bisection) {
    const auto &mesh = p_mesh.get_serial_mesh();
    std::vector<std::set<Edge> > edges(comm_.size());

    std::vector<Integer> sharing;
    for (Integer i = 0; i < mesh.n_elements(); ++i) {
      if (!mesh.is_active(i)) continue;

      const auto &e = mesh.elem(i);

      for (Integer i = 0; i < n_edges(e); ++i) {
        Edge edge;
        e.edge(i, edge.nodes[0], edge.nodes[1]);
        edge.fix_ordering();

        assert(edge.is_valid());

        Edge global_edge(p_mesh.node_map().global(edge.nodes[0]),
                         p_mesh.node_map().global(edge.nodes[1]));

        auto &inter = edge_interface_[global_edge];
        inter.clear();

        if (!inter.empty()) continue;

        p_mesh.node_map().intersect_partitions(std::begin(edge.nodes),
                                               std::end(edge.nodes), sharing);

        if (sharing.size() == 1) {
          assert(sharing[0] == comm_.rank());
          inter.push_back(comm_.rank());
          continue;
        }

        for (auto s : sharing) {
          if (s != comm_.rank()) {
            edges[s].insert(global_edge);
          } else {
            // This ensures that this edge is considered only once
            inter.push_back(comm_.rank());
          }
        }
      }
    }

    // serialize
    for (Integer r = 0; r < comm_.size(); ++r) {
      if (r == comm_.rank()) continue;
      output_[r].clear();
      write(edges[r], output_[r]);
    }

    for (Integer r = 0; r < comm_.size(); ++r) {
      if (r == comm_.rank()) continue;

      send_buff_[r] = output_[r].str();
      assert(!send_buff_[r].empty());

      comm_.i_send(&send_buff_[r][0], send_buff_[r].size(), r, r);
    }

    Integer n_connections = comm_.size() - 1;

    for (Integer i = 0; i < n_connections; ++i) {
      Integer rank, size;
      while (!comm_.i_probe_any<byte>(&rank, &size)) {
      }
      recv_buff_[rank].resize(size);
      comm_.i_recv(&recv_buff_[rank][0], size, rank, comm_.rank());
    }

    for (Integer i = 0; i < n_connections; ++i) {
      Integer rank, index;
      while (!comm_.test_recv_any(&rank, &index)) {
      }

      Integer n_edges = 0;
      InputStream is(recv_buff_[rank]);
      read(n_edges, is);

      for (Integer k = 0; k < n_edges; ++k) {
        Edge global_edge;
        read(global_edge, is);

        if (edges[rank].find(global_edge) != edges[rank].end()) {
          is_edge_interfaced_[rank] = true;

          auto &inter = edge_interface_[global_edge];
          inter.push_back(rank);
        }
      }
    }

    comm_.wait_all();
  }

  void synchronize() {}

  template <Integer Dim, Integer ManifoldDim>
  void update_midpoint_parts(const EdgeNodeMap &edge_node_map,
                             ParMesh<Dim, ManifoldDim> &mesh) {}

  template <Integer Dim, Integer ManifoldDim>
  void update(const ParMesh<Dim, ManifoldDim> &mesh,
              const EdgeElementMap &edge_elem_map,
              const EdgeNodeMap &edge_node_map,
              const std::vector<Edge> &bisected_edges) {}

  template <Integer Dim, Integer ManifoldDim>
  void write_to_mesh(const EdgeNodeMap &edge_node_map,
                     ParMesh<Dim, ManifoldDim> &mesh) {}

  template <Integer Dim, Integer ManifoldDim>
  void read_from_mesh(const EdgeNodeMap &edge_node_map,
                      const ParMesh<Dim, ManifoldDim> &mesh) {}

  template <Integer Dim, Integer ManifoldDim>
  void collect_splits_to_local_edges(const ParMesh<Dim, ManifoldDim> &mesh,
                                     std::vector<Edge> &local_splits) {
    local_splits.clear();
    local_splits.reserve(splits_.size());

    for (auto gs_it = splits_.begin(); gs_it != splits_.end(); ++gs_it) {
      const auto &gs = *gs_it;
      Edge le(gs.edge);
      mesh.make_local(le.nodes);
      le.fix_ordering();
      local_splits.push_back(le);
    }
  }

  bool empty() const {
    return comm_.reduce_and(splits_.empty() && to_communicate_.empty());
  }

  void describe(std::ostream &os) const {
    serial_apply(comm_, [&os, this]() {
      os << "---------------------------\n";
      os << comm_ << " -> ";

      os << "(";
      for (Integer i = 0; i < comm_.size(); ++i) {
        if (is_edge_interfaced_[i]) {
          os << " " << i;
        }
      }

      os << " )\n";

      for (const auto &ei : edge_interface_) {
        ei.first.describe(os);
        os << " ->";
        for (auto r : ei.second) {
          os << " " << r;
        }
        os << "\n";
      }

      os << std::flush;
    });
  }

  /////////////////////////////////////////////////////

  bool add_split(const EdgeSplit &edge_split) {
    assert(!edge_split.partitions.empty());

    auto &es = get_split(edge_split.edge);
    if (es.midpoint != INVALID_INDEX) {
      // does not need to comunicate
      // is global already
      return true;
    }

    // put the split in the global pool
    const EdgeSplit &gs = add_global_split(comm_.rank(), edge_split);

    if (gs.only_on_partition(comm_.rank())) {
      return false;
    }

    if (!es.is_valid()) {
      to_communicate_.push_back(gs);
    }

    return false;
  }

  const EdgeSplit &add_global_split(const Integer originator,
                                    const EdgeSplit &edge_split) {
    assert(!edge_split.partitions.empty());
    assert(edge_to_split_map_.size() == splits_.size());

    auto e_temp = edge_split.edge;
    e_temp.fix_ordering();

    auto it = edge_to_split_map_.find(e_temp);

    if (it == edge_to_split_map_.end()) {
      Integer split_id = splits_.size();
      edge_to_split_map_[e_temp] = split_id;
      splits_.push_back(edge_split);
      auto &es = splits_.back();
      es.owner = originator;
      assert(edge_to_split_map_.size() == splits_.size());

      return es;
    }

    EdgeSplit &g_split = get_split(it->second);
    if (g_split.midpoint != INVALID_INDEX) {
      // ownership and global id already determined
      return g_split;
    }

    // update ids
    if (edge_split.midpoint != INVALID_INDEX) {
      assert(edge_split.owner != INVALID_INDEX);
      g_split = edge_split;
      return g_split;
    }

    // determine owner
    if (g_split.owner == INVALID_INDEX) {
      g_split.owner = originator;
    } else {
      g_split.owner = std::min(originator, g_split.owner);
    }

    assert(edge_to_split_map_.size() == splits_.size());
    return g_split;
  }

  inline const EdgeSplit &get_split(const Edge &e) const {
    static const EdgeSplit null_;

    Edge e_temp = e;
    e_temp.fix_ordering();

    auto it = edge_to_split_map_.find(e_temp);
    if (it == edge_to_split_map_.end()) {
      return null_;
    }

    return get_split(it->second);
  }

  inline const EdgeSplit &get_split(const Integer split_id) const {
    assert(split_id >= 0);
    assert(split_id < splits_.size());
    assert(edge_to_split_map_.size() == splits_.size());

    return splits_[split_id];
  }

  inline EdgeSplit &get_split(const Edge &e) {
    static EdgeSplit null_;

    Edge e_temp = e;
    e_temp.fix_ordering();

    auto it = edge_to_split_map_.find(e_temp);
    if (it == edge_to_split_map_.end()) {
      return null_;
    }

    return get_split(it->second);
  }

  inline EdgeSplit &get_split(const Integer split_id) {
    assert(split_id >= 0);
    assert(split_id < splits_.size());
    assert(edge_to_split_map_.size() == splits_.size());

    return splits_[split_id];
  }

  inline Integer owner(const Edge &e) const {
    auto it = edge_to_split_map_.find(e);
    if (it == edge_to_split_map_.end()) {
      return INVALID_INDEX;
    }

    return get_split(it->second).owner;
  }

  inline void resolve_midpoint_id(const Edge &e, const Integer local_m_id,
                                  const Integer m_id, Map &node_map) {
    auto &s = get_split(e);
    assert(s.is_valid());
    assert(s.midpoint == INVALID_INDEX);
    s.midpoint = m_id;

    if (!s.only_on_partition(comm_.rank())) {
      to_communicate_.push_back(s);
    }

    update_edge_interface(e, m_id);

    std::vector<Integer> parts;
    edge_interface(e, parts);
    assert(!parts.empty());

    node_map.set_partitions(local_m_id, std::begin(parts), std::end(parts));
  }

  void edge_interface(const Edge &global_edge,
                      std::vector<Integer> &partitions) const {
    auto it = edge_interface_.find(global_edge);
    if (it == edge_interface_.end()) {
      partitions.clear();
      return;
    }

    partitions = it->second;
  }

  void update_edge_interface(const Edge &global_edge,
                             const Integer global_midpoint) {
    Edge e1(global_edge.nodes[0], global_midpoint);
    Edge e2(global_edge.nodes[1], global_midpoint);

    auto it = edge_interface_.find(global_edge);
    const auto parent_interface = it->second;

    assert(it != edge_interface_.end());
    assert(!it->second.empty());
    // assert(edge_interface_[e1].empty());
    // assert(edge_interface_[e2].empty());

    auto &ei1 = edge_interface_[e1];
    if (ei1.empty()) {
      ei1 = parent_interface;
    }

    auto &ei2 = edge_interface_[e2];
    if (ei2.empty()) {
      ei2 = parent_interface;
    }
  }

  inline Integer midpoint(const Edge &e) const {
    auto it = edge_to_split_map_.find(e);
    if (it == edge_to_split_map_.end()) {
      return INVALID_INDEX;
    }

    return get_split(it->second).midpoint;
  }

  bool check_consistent() const {
    assert(edge_to_split_map_.size() == splits_.size());

    if (splits_.empty()) return true;

    std::vector<bool> visited(splits_.size(), false);

    for (auto etoid : edge_to_split_map_) {
      auto &s = get_split(etoid.second);
      assert(s.edge == etoid.first);

      if (s.edge != etoid.first) return false;
      visited[etoid.second] = true;
    }

    for (Integer i = 0; i < visited.size(); ++i) {
      if (!visited[i]) {
        std::cerr << "split(" << i << ") not visited: ";
        get_split(i).describe(std::cerr);
        std::cerr << std::endl;
        return false;
      }
    }

    return true;
  }

  template <class ParMesh>
  void set_midpoint_ids_to(const EdgeNodeMap &enm, ParMesh &part) {
    assert(check_consistent());
    if (splits_.empty()) return;

    Integer remove_index = 0;
    for (auto gs_it = splits_.begin(); gs_it != splits_.end();) {
      const auto &gs = *gs_it;
      auto map_it = edge_to_split_map_.find(gs.edge);
      map_it->second -= remove_index;

      Edge e = part.local_edge(gs.edge);

      if (!e.is_valid()) {
        std::cout << "invalid local edge: ";
        gs.describe(std::cout);
        std::cout << "\n";

        ++gs_it;
        continue;
      }

      Integer local_mp = enm.get(e);

      if (local_mp == INVALID_INDEX) {
        ++gs_it;
        continue;
      }

      assert(local_mp != INVALID_INDEX);

      part.assign_node_owner(local_mp, gs.owner);

      if (gs.midpoint == INVALID_INDEX) {
        ++gs_it;
        continue;
      }

      part.assign_node_local_to_global(local_mp, gs.midpoint);
      update_edge_interface(gs.edge, gs.midpoint);

      std::vector<Integer> parts;
      edge_interface(gs.edge, parts);
      assert(!parts.empty());
      part.node_map().set_partitions(local_mp, std::begin(parts),
                                     std::end(parts));

      if (map_it != edge_to_split_map_.end()) {
        edge_to_split_map_.erase(map_it);
      } else {
        assert(false);
      }

      splits_.erase(gs_it);
      remove_index++;
    }

    assert(edge_to_split_map_.size() == splits_.size());
  }

  template <class ParMesh>
  void collect_splits_to_apply(const ParMesh &part,
                               std::vector<EdgeSplit> &local_splits) {
    assert(edge_to_split_map_.size() == splits_.size());

    local_splits.clear();
    local_splits.reserve(splits_.size());

    for (auto gs_it = splits_.begin(); gs_it != splits_.end();) {
      const auto &gs = *gs_it;
      local_splits.push_back(gs);
    }
  }

  inline std::vector<EdgeSplit>::iterator begin() { return splits_.begin(); }

  inline std::vector<EdgeSplit>::iterator end() { return splits_.end(); }

  inline std::vector<EdgeSplit>::const_iterator begin() const {
    return splits_.begin();
  }

  inline std::vector<EdgeSplit>::const_iterator end() const {
    return splits_.end();
  }

  /////////////////////////////////////////////////////

  ParEdgeSplitPool(const Communicator &comm)
      : comm_(comm),
        output_(comm.size()),
        send_buff_(comm.size()),
        recv_buff_(comm.size()),
        is_edge_interfaced_(comm.size(), false) {}

  /////////////////////////////////////////////////////

 private:
  Communicator comm_;
  std::map<Edge, std::vector<Integer> > edge_interface_;

  std::vector<OutputStream> output_;
  std::vector<BufferObject> send_buff_;
  std::vector<BufferObject> recv_buff_;
  std::vector<bool> is_edge_interfaced_;

  std::vector<EdgeSplit> splits_;
  std::map<Edge, Integer> edge_to_split_map_;
  std::vector<EdgeSplit> to_communicate_;
};
}  // namespace mars

#endif  // MARS_PAR_EDGE_SPLIT_POOL_HPP
