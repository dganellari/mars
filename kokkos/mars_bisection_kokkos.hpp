#ifndef MARS_BISECTION_KOKKOS_HPP
#define MARS_BISECTION_KOKKOS_HPP

#include <iostream>
#include <typeinfo>
// #include "mars_dof_map.hpp"
#include "mars_edge_element_map_kokkos.hpp"
#include "mars_edge_node_map_kokkos.hpp"
#include "mars_edge_select_kokkos.hpp"
#include "mars_longest_edge_kokkos.hpp"
#include "mars_mark_kokkos.hpp"
// #include "mars_tracker.hpp"

namespace mars {
    template <class Mesh_, class EdgeSelect_ = LongestEdgeSelect<Mesh_, KokkosImplementation>>
    class ParallelBisection {
    public:
        using Mesh = Mesh_;
        using Elem = typename Mesh::Elem;
        using SideElem = typename Mesh::SideElem;
        // struct Params { static constexpr Integer ManifoldDim=Mesh::ManifoldDim; using Impl=KokkosImplementation; };
        using ParallelEdgeElementMap = ParallelSubManifoldElementMap<2, Mesh::ManifoldDim>;
        using ParallelEdgeNodeMap = DeviceEdgeNodeMap<2>;
        using Edge = Side<2, KokkosImplementation>;
        using Comb = Combinations<Mesh::ManifoldDim + 1, 2, KokkosImplementation>;
        using CombView = ViewMatrixTextureC<Integer, Comb::value, 2>;

        using ElementVector = typename ParallelEdgeElementMap::ElementVector;

        static constexpr Integer Dim = Mesh_::Dim;
        static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

        virtual ~ParallelBisection() {
            if (verbose)
                std::cout << "Calling the ParallelBisection destructor to deallocate the composition singleton!"
                          << std::endl;
            // Not possible to call the destructor of the view since it will be called again
            // Instead an empty view with data()=nullptr is assigned to it to avoid
            // the deallocation of the static view after the kokkos initialize.
            // Comb::instance().combs.~View();
            Comb::instance(false).combs = decltype(Comb::instance().combs)();
        }

        ParallelBisection(Mesh* mesh) : mesh(nullptr), host_mesh(mesh), verbose(false), fail_if_not_refine(false) {}

        void set_fail_if_not_refine(const bool val) { fail_if_not_refine = val; }

        bool get_fail_if_not_refine() { return fail_if_not_refine; }

        Mesh* get_mesh() const { return mesh; }

        Mesh* get_host_mesh() const { return host_mesh; }

        void set_verbose(const bool val) { verbose = val; }

        /*

        inline Integer side_num(
                const Integer element_id,
                const SideElem &side) const
        {
                auto nodes = side.nodes;
                std::sort(nodes.begin(), nodes.end());

                const auto &e = mesh.elem(element_id);

                SideElem e_side;

                for(Integer i = 0; i < n_sides(e); ++i) {
                        e.side(i, e_side);
                        std::sort(std::begin(e_side.nodes), std::end(e_side.nodes));

                        bool same_side = true;
                        for(Integer k = 0; k < nodes.size(); ++k) {
                                assert(nodes[k] != INVALID_INDEX);
                                assert(e_side.nodes[k] != INVALID_INDEX);

                                if(nodes[k] != e_side.nodes[k]) {
                                        same_side = false;
                                        break;
                                }
                        }

                        if(same_side) {
                                return i;
                        }
                }

                return INVALID_INDEX;
        }


        void bisect_side_tags(
                const Integer element_id,
                const Edge &edge,
                const Integer midpoint_id)
        {
                const auto &e = mesh.elem(element_id);

                for(auto c : e.children) {
                        std::fill(mesh.elem(c).side_tags.begin(),
                                mesh.elem(c).side_tags.end(),
                                INVALID_INDEX);
                }

                SideElem side;
                SideElem child_side;

                for(Integer i = 0; i < n_sides(e); ++i) {
                        e.side(i, side);
                        const Integer tag = e.side_tags[i];
                        if(tag == INVALID_INDEX) continue;

                        if(has_edge(side, edge[0], edge[1])) {
                                //tag of split sides
                                child_side.nodes[0] = midpoint_id;

                                for(Integer j = 0; j < 2; ++j) {
                                        const Integer vj = edge[j];
                                        child_side.nodes[1] = vj;

                                        Integer local_ind = 2;
                                        for(auto s : side.nodes) {
                                                if(s == edge[0] || s == edge[1]) {
                                                        continue;
                                                }
                                                child_side.nodes[local_ind++] = s;
                                        }

                                        bool found_side = false;
                                        for(auto c : e.children) {
                                                auto sn = side_num(c, child_side);

                                                if(INVALID_INDEX != sn) {
                                                        assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX ||
                                                                mesh.elem(c).side_tags[sn] == tag);

                                                        mesh.elem(c).side_tags[sn] = tag;
                                                        found_side = true;
                                                }
                                        }

                                        assert(found_side);
                                }
                        } else {
                                bool found_side = false;

                                for(auto c : e.children) {
                                        auto sn = side_num(c, side);

                                        if(INVALID_INDEX != sn) {
                                                assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX ||
                                                        mesh.elem(c).side_tags[sn] == tag);

                                                mesh.elem(c).side_tags[sn] = tag;
                                                found_side = true;
                                        }
                                }

                                assert(found_side);
                        }
                }
        }

        const EdgeNodeMap &edge_node_map() const
        {
                return edge_node_map_;
        }

        const EdgeElementMap &edge_element_map() const
        {
                return edge_element_map_;
        }

        EdgeElementMap &edge_element_map()
        {
                return edge_element_map_;
        }


        const Mesh &get_mesh() const
        {
                return mesh;
        }

        Mesh &get_mesh()
        {
                return mesh;
        }

        const std::vector<Edge> &bisected_edges() const
        {
                return bisected_edges_;
        }

        void clear_bisected_edges()
        {
                bisected_edges_.clear();
        }

        void refine_edges(const std::vector<Edge> &edges)
        {
                if(flags.empty()) {
                        flags.resize(mesh.n_elements(), NONE);
                        level.resize(mesh.n_elements(), 0);
                        edge_element_map_.update(mesh);
                        mesh.update_dual_graph();
                }

                for(auto e : edges) {
                        refine_edge(e);
                }

                mesh.update_dual_graph();
                mesh.tags() = level;
        }


        void clear()
        {
                flags.clear();
                level.clear();
                side_flags.clear();
                edge_node_map_.clear();
                edge_element_map_.clear();
        }

        void tracking_begin()
        {
                tracker_.begin_iterate();
        }

        void tracking_end()
        {
                tracker_.end_iterate();
        }

        void undo()
        {
                tracker_.undo_last_iterate(mesh);
        }*/

        // does ont perform a deep copy of the view containted in the mesh. Just the mesh object.
        void copy_mesh_to_device() {
            Mesh* tmp = (Mesh*)Kokkos::kokkos_malloc(sizeof(Mesh));
            Mesh mCopy = *host_mesh;
            Mesh* oldDeviceMesh = mesh;
            Kokkos::parallel_for(
                "CreateMeshObject", 1, KOKKOS_LAMBDA(const int&) {
                    new ((Mesh*)tmp)
                        Mesh(mCopy);  // two local copies for m and tmp since this->m, this->mesh host pointers
                    if (oldDeviceMesh) oldDeviceMesh->~Mesh();
                    // it only works on a copy since **m is still a host pointer and fails on the device.
                });

            mesh = tmp;  // make the mesh pointer a device one so that this init func is not neccessary anymore
        }

        Integer euler_graph_formula(Mesh* mesh) {
            /*Mesh3 sMesh3;
            read_mesh("../data/write/euler.MFEM", sMesh3);*/

            typename Mesh::SerialMesh sMesh;
            convert_parallel_mesh_to_serial(sMesh, *mesh);

            Kokkos::Timer timer;

            double time = timer.seconds();
            std::cout << "serial kokkos took: " << time << " seconds." << std::endl;

            // map_.describe(std::cout);

            sMesh.update_dual_graph();

            // nr_faces per eleme counted twicee * nr elements
            Integer bound = Combinations<ManifoldDim + 1, ManifoldDim>::value * sMesh.n_elements();

            // Euler's formula for graphs
            Integer faces_nr = (bound + sMesh.n_boundary_sides()) / 2;  // boundary sides counted only once
            // Integer edges_nr = sMesh3.n_nodes() + faces_nr - 2;

            Integer edges_nr = 0;

            if (ManifoldDim == 2)
                edges_nr = sMesh.n_nodes() + sMesh.n_elements() + 1 - 2;  // +1 because also the outer face is counted
            else
                edges_nr = sMesh.n_nodes() + faces_nr - 1 - sMesh.n_elements();

            // edges_nr = 2 * mesh->n_elements();

            return edges_nr;
        }

        MARS_INLINE_FUNCTION
        static void insert_longest_edge(const Mesh_* mesh,
                                        const UnorderedMap<Edge, ElementVector>& mapping,
                                        const Integer element_id) {
            EdgeSelect_ es;
            const Integer edge_num = es.stable_select(*mesh, element_id);
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id,
            // Bisection<Mesh>::edge_element_map());
            Edge edge;
            mesh->elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
            edge.fix_ordering();

            auto result = mapping.insert(edge);

            if (result.failed()) printf("Exceeded UnorderedMap capacity: Edge Element Map\n");
        }

        struct BuildLEPMap {
            Mesh_* mesh;
            UnorderedMap<Edge, ElementVector> mapping;
            ViewVectorType<Integer> active_elems;

            BuildLEPMap(UnorderedMap<Edge, ElementVector> mp, Mesh_* ms, ViewVectorType<Integer> ae)
                : mapping(mp), mesh(ms), active_elems(ae) {}

            MARS_INLINE_FUNCTION
            void operator()(int i) const {
                const Integer element_id = active_elems(i);
                insert_longest_edge(mesh, mapping, element_id);
            }
        };

        void reserve_tree(Integer capacity) { tree_ = UnorderedMap<Integer, ElementVector>(capacity); }

        void reserve_leppOccupied(Integer size) { lepp_occupied = ViewVectorType<bool>("lepp_occupied_fill", size); }

        void precompute_lepp_incidents(Mesh_* mesh, const ViewVectorType<Integer> active_elems) {
            Kokkos::parallel_for(active_elems.extent(0), BuildLEPMap(edge_element_map_.mapping_, mesh, active_elems));
        }

        MARS_INLINE_FUNCTION
        static Integer is_leaf(const Integer element_id,
                               Mesh_* mesh,
                               const UnorderedMap<Integer, ElementVector>& tree) {
            Integer index = tree.find(element_id);

            if (tree.valid_at(index)) {
                for (auto i = 0; i < tree.value_at(index).index; ++i) {
                    if (mesh->is_active(tree.value_at(index)[i])) {
                        return -1;
                    }
                }
            }

            return index;
        }

        struct CountLepp {
            Mesh_* mesh;
            ViewVectorType<Integer> elements;
            UnorderedMap<Edge, ElementVector> mapping;
            UnorderedMap<Integer, ElementVector> tree;
            ViewVectorType<bool> lepp_occupied;

            ViewVectorType<Integer> index_count;
            ViewVectorType<Integer> pt_count;

            CombView combs;

            CountLepp(UnorderedMap<Edge, ElementVector> mp,
                      Mesh_* ms,
                      ViewVectorType<Integer> elems,
                      UnorderedMap<Integer, ElementVector> tr,
                      ViewVectorType<bool> lo,
                      ViewVectorType<Integer> ic,
                      ViewVectorType<Integer> pc,
                      CombView cmbs)
                : mapping(mp),
                  mesh(ms),
                  elements(elems),
                  tree(tr),
                  lepp_occupied(lo),
                  index_count(ic),
                  pt_count(pc),
                  combs(cmbs) {}

            CountLepp(UnorderedMap<Edge, ElementVector> mp,
                      Mesh_* ms,
                      ViewVectorType<Integer> elems,
                      UnorderedMap<Integer, ElementVector> tr,
                      ViewVectorType<bool> lo,
                      CombView cmbs)
                : mapping(mp), mesh(ms), elements(elems), tree(tr), lepp_occupied(lo), combs(cmbs) {}

            MARS_INLINE_FUNCTION
            bool is_terminal(const Integer map_index, const Edge& edge, const ElementVector& incidents) const {
                bool terminal =
                    true;  // in case the elements share the longest edge or there is only one incident (itself)
                           // - meaning the longest edge is on the boundary.

                for (auto i = 0; i < incidents.index; ++i) {
                    if (mesh->is_active(incidents[i])) {
                        Edge new_edge;

                        EdgeSelect_ es;
                        const Integer edge_num = es.stable_select(*mesh, incidents[i]);
                        // const Integer edge_num = es.select(*mesh, edge, incidents[i]);

                        mesh->elem(incidents[i]).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
                        new_edge.fix_ordering();

                        //	if(edge != tree.edges[i]){
                        if (edge != new_edge) {
                            terminal = false;
                            // insert child in tree
                            auto result = tree.insert(incidents[i]);
                            if (result.failed()) printf("Exceeded UnorderedMap: tree capacity\n");
                            tree.value_at(map_index).insert(incidents[i]);
                        }
                    }
                }

                return terminal;
            }

            MARS_INLINE_FUNCTION
            void compute_lepp(const Integer element_id,
                              const Integer map_index,
                              Integer& count,
                              Integer& pt_count) const {
                Edge edge;

                EdgeSelect_ es;
                const Integer edge_num = es.stable_select(*mesh, element_id);
                // const Integer edge_num = es.select(*mesh, edge, incidents[i]);
                mesh->elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
                edge.fix_ordering();

                Integer index = mapping.find(edge);

                if (is_terminal(map_index, mapping.key_at(index), mapping.value_at(index))) {
                    ++pt_count;

                    for (auto i = 0; i < mapping.value_at(index).index; ++i) {
                        if (mesh->is_active(mapping.value_at(index)[i])) {
                            ++count;
                        }
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void depth_first(const Integer node, Integer& count, Integer& pt_count) const {
                // Avoids lepp path collisions. If the value is alrady set to 1 the threads returns.
                if (!Kokkos::atomic_compare_exchange_strong(&lepp_occupied(node), false, true)) {
                    return;
                }
                Integer index = is_leaf(node, mesh, tree);

                if (index >= 0) {
                    compute_lepp(node, index, count, pt_count);
                } else {
                    index = tree.find(node);
                }

                for (auto i = 0; i < tree.value_at(index).index; ++i) {
                    const auto& child = tree.value_at(index)[i];
                    if (mesh->is_active(child)) {
                        depth_first(child, count, pt_count);
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const int i) const {
                Integer element_id = elements(i);
                tree.insert(element_id);
                // Integer terminal_elem_count=0;
                Integer terminal_incident_count = 0;
                Integer terminal_pt_count = 0;

                depth_first(element_id, terminal_incident_count, terminal_pt_count);

                index_count(i + 1) = terminal_incident_count;
                pt_count(i + 1) = terminal_pt_count;
                //+1 for leaving the first cell 0 and performing an inclusive scan on the rest
                // to have both exclusive and inclusive and the total on the last cell.
            }
        };

        struct ScatterElem : CountLepp {
            UnorderedMap<Integer, Edge> lepp_incidents_map;
            UnorderedMap<Edge, Integer> lepp_edge_map;

            ScatterElem(UnorderedMap<Edge, ElementVector> mp,
                        Mesh_* ms,
                        ViewVectorType<Integer> elems,
                        UnorderedMap<Integer, ElementVector> tr,
                        ViewVectorType<bool> lo,
                        CombView cmbs,
                        UnorderedMap<Integer, Edge> lim,
                        UnorderedMap<Edge, Integer> lem)
                : CountLepp(mp, ms, elems, tr, lo, cmbs),  // lepp_elem_index(lei),
                  lepp_incidents_map(lim),
                  lepp_edge_map(lem) {}

            MARS_INLINE_FUNCTION
            bool is_terminal(const Edge& edge, const ElementVector& incidents) const {
                bool terminal =
                    true;  // in case the elements share the longest edge or there is only one incident (itself)
                           // - meaning the longest edge is on the boundary.

                for (auto i = 0; i < incidents.index; ++i) {
                    if (this->mesh->is_active(incidents[i])) {
                        Edge new_edge;

                        EdgeSelect_ es;
                        const Integer edge_num = es.stable_select(*this->mesh, incidents[i]);
                        // const Integer edge_num = es.select(*mesh, edge, incidents[i]);

                        this->mesh->elem(incidents[i]).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
                        new_edge.fix_ordering();

                        //	if(edge != tree.edges[i]){
                        if (edge != new_edge) {
                            terminal = false;
                        }
                    }
                }

                return terminal;
            }

            MARS_INLINE_FUNCTION
            void compute_lepp(const Integer element_id, const Integer map_index) const {
                Edge edge;

                EdgeSelect_ es;
                const Integer edge_num = es.stable_select(*this->mesh, element_id);
                // const Integer edge_num = es.select(*mesh, edge, incidents[i]);

                this->mesh->elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
                edge.fix_ordering();

                Integer index = this->mapping.find(edge);

                if (is_terminal(this->mapping.key_at(index), this->mapping.value_at(index))) {
                    auto res = lepp_edge_map.insert(this->mapping.key_at(index), index);
                    assert(!res.failed());

                    /*if (res.failed())
                            printf("Exceeded UnorderedMap: lepp_edge_map capacity\n");*/

                    for (auto i = 0; i < this->mapping.value_at(index).index; ++i) {
                        auto& element = this->mapping.value_at(index)[i];
                        if (this->mesh->is_active(element)) {
                            auto result = lepp_incidents_map.insert(element, this->mapping.key_at(index));
                            assert(!result.failed());

                            /*if (result.failed())
                                    printf("Exceeded UnorderedMap: lepp_incidents_map capacity\n");*/
                        }
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void depth_first(const Integer node) const {
                // Avoids lepp path collisions. If the value is alrady set to 1 the threads returns.
                if (!Kokkos::atomic_compare_exchange_strong(&this->lepp_occupied(node), true, false)) {
                    return;
                }

                Integer index = is_leaf(node, this->mesh, this->tree);

                if (index >= 0) {
                    compute_lepp(node, index);
                } else {
                    index = this->tree.find(node);
                }

                for (auto i = 0; i < this->tree.value_at(index).index; ++i) {
                    const auto& child = this->tree.value_at(index)[i];
                    if (this->mesh->is_active(child)) {
                        depth_first(child);
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void operator()(const int i) const {
                Integer element_id = this->elements(i);
                // this->tree.insert(element_id);
                depth_first(element_id);
            }
        };

        struct BisectEdge {
            Mesh_* mesh;
            bool verbose;
            using UEMap = UnorderedMap<Edge, ElementVector>;
            using UNMap = UnorderedMap<Edge, Integer>;

            ViewMatrixType<Integer> lepp_point_index;
            UEMap edge_element_map;
            UNMap edge_node_map;

            Integer node_start_index;

            BisectEdge(Mesh_* ms, ViewMatrixType<Integer> lpi, UEMap mp, UNMap nmp, bool v, Integer nsi)
                : mesh(ms),
                  lepp_point_index(lpi),
                  edge_element_map(mp),
                  edge_node_map(nmp),
                  verbose(v),
                  node_start_index(nsi) {}

            MARS_INLINE_FUNCTION
            void operator()(const int i) const {
                const Integer v1 = edge_element_map.key_at(lepp_point_index(i, 0))[0];
                const Integer v2 = edge_element_map.key_at(lepp_point_index(i, 0))[1];

                Integer midpoint = i + node_start_index;

                mesh->add_point((mesh->point(v1) + mesh->point(v2)) / 2., midpoint);

                auto index = edge_node_map.find(edge_element_map.key_at(lepp_point_index(i, 0)));
                if (edge_node_map.valid_at(index)) {
                    edge_node_map.value_at(index) = midpoint;
                }
            }
        };

        struct BisectElem {
            Mesh_* mesh;
            ViewMatrixType<Integer> lepp_incident_index;
            static constexpr Integer value = Factorial<Mesh_::ManifoldDim + 1>::value /
                                             (Factorial<2>::value * Factorial<Mesh_::ManifoldDim + 1 - 2>::value);
            ViewMatrixTextureC<Integer, value, 2> combs;

            using UEMap = UnorderedMap<Edge, ElementVector>;
            using UNMap = UnorderedMap<Edge, Integer>;

            UEMap edge_element_map;
            UNMap edge_node_map;

            Integer elem_start_index;
            Integer lepp_size;

            BisectElem(Mesh_* ms,
                       ViewMatrixType<Integer> lii,
                       ViewMatrixTextureC<Integer, value, 2> cmb,
                       UEMap mp,
                       UNMap nmp,
                       Integer esi,
                       Integer csi)
                : mesh(ms),
                  lepp_incident_index(lii),
                  combs(cmb),
                  edge_element_map(mp),
                  edge_node_map(nmp),
                  elem_start_index(esi),
                  lepp_size(csi) {}

            MARS_INLINE_FUNCTION
            void other_nodes(const SubView<Integer, ManifoldDim + 1>& nodes,
                             const Integer v1,
                             const Integer v2,
                             TempArray<Integer, ManifoldDim - 1>& opposite_nodes) const {
                Integer i = 0;
                for (Integer j = 0; j < ManifoldDim + 1; ++j) {
                    Integer n = nodes[j];
                    if (n != v1 && n != v2) {
                        opposite_nodes[i++] = n;
                    }
                }
            }

            MARS_INLINE_FUNCTION
            void edge_index_in_elem(const SubView<Integer, ManifoldDim + 1>& nodes,
                                    const Integer v1,
                                    const Integer v2,
                                    TempArray<Integer, 2>& opposite_nodes) const {
                Integer i = 0;
                for (Integer j = 0; j < ManifoldDim + 1; ++j) {
                    Integer n = nodes[j];
                    if (n == v1 || n == v2) {
                        opposite_nodes[i++] = j;
                    }
                }
            }

            /*MARS_INLINE_FUNCTION
            void update_elem(Integer elem_new_id) const
            {

                    flags.push_back(NONE);
                    level.push_back(
                            (mesh.elem(id).parent_id == INVALID_INDEX) ?
                                    0 : (level[mesh.elem(id).parent_id] + 1));
            }*/

            MARS_INLINE_FUNCTION
            void add_new_elem(const Integer elem_new_id,
                              Elem& old_el,
                              const Integer nodes_new_id,
                              const Integer index) const {
                Elem new_el = mesh->elem(elem_new_id);
                new_el.set_parent_id(old_el.id);

                // in order to keep the counterclockwise ordering of the element nodes.
                for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                    if (i == index)
                        new_el.nodes[i] = nodes_new_id;
                    else
                        new_el.nodes[i] = old_el.nodes[i];
                }

                mesh->set_active(elem_new_id);
                ParallelEdgeElementMap::update_elem(edge_element_map, new_el, combs);

                new_el.block = old_el.block;
            }

            /*MARS_INLINE_FUNCTION
            Integer add_node(const Integer v1,	const Integer v2) const
            {
                    TempArray<Integer, 2> nodes;
                    nodes[0] = v1;
                    nodes[1] = v2;

                    Edge edge(nodes);

                    Integer midpoint=0;

                    if (ParallelEdgeNodeMap::update(edge, edge_node_map, node_start_index, midpoint))
                    {
                            assert(midpoint!=INVALID_INDEX);

                            mesh->add_point((mesh->point(v1) + mesh->point(v2)) / 2.,
                                            midpoint);
                    }

                    return midpoint;
            }*/

            MARS_INLINE_FUNCTION
            Integer get_new_node(const Integer v1, const Integer v2) const {
                TempArray<Integer, 2> nodes;
                nodes[0] = v1;
                nodes[1] = v2;

                Edge edge(nodes);

                Integer midpoint = INVALID_INDEX;

                Integer index = edge_node_map.find(edge);
                if (edge_node_map.valid_at(index)) {
                    midpoint = edge_node_map.value_at(index);
                }

                return midpoint;
            }

            MARS_INLINE_FUNCTION
            void operator()(const int i) const {
                Integer element_id = lepp_incident_index(i, 0);
                Integer v1 = lepp_incident_index(i, 1);
                Integer v2 = lepp_incident_index(i, 2);

                /*if (verbose)
                 {
                 printf("bisect(%li , %li) for %li\n", v1, v2, element_id);
                 }*/

                Integer nodes_new_id = get_new_node(v1, v2);

                Integer elem_new_id = elem_start_index + 2 * i;
                Elem old_el = mesh->elem(element_id);

                TempArray<Integer, 2> edge_indices;
                edge_index_in_elem(old_el.nodes, v1, v2, edge_indices);

                add_new_elem(elem_new_id, old_el, nodes_new_id, edge_indices[0]);
                add_new_elem(++elem_new_id, old_el, nodes_new_id, edge_indices[1]);
            }
        };

        struct HandleParentElem {
            Mesh_* mesh;
            ViewMatrixType<Integer> lepp_incident_index;

            Integer elem_start_index;
            Integer child_start_index;

            HandleParentElem(Mesh_* ms, ViewMatrixType<Integer> lii, Integer esi, Integer csi)
                : mesh(ms), lepp_incident_index(lii), elem_start_index(esi), child_start_index(csi) {}

            MARS_INLINE_FUNCTION
            void operator()(const int i) const {
                Integer element_id = lepp_incident_index(i, 0);
                mesh->set_active(element_id, false);

                Integer elem_new_id = elem_start_index + 2 * i;
                Integer childrens_id = child_start_index + i;

                Elem old_el = mesh->elem(element_id, childrens_id);

                old_el.children(0) = elem_new_id;
                old_el.children(1) = ++elem_new_id;
                const auto result = mesh->get_elem_children_map().insert(element_id, childrens_id);
                assert(result.success());
            }
        };

        inline void free_mesh() {
            /*const Mesh* tmp = ref.mesh;

             parallel_for("DestroyMeshObject",1, KOKKOS_LAMBDA (const int&) {
             mesh->~Mesh();
             });

             fence();

             kokkos_free(mesh);
             */
        }

        /*void compact_map_to_view(const UnorderedMap<Integer, Edge>& lepp_incidents_map,
                        ViewMatrixType<Integer> lepp_incident_index)
        {
                using namespace Kokkos;

                Timer timer;

                ViewObject<Integer> global_index("global_index");

                parallel_for(lepp_incidents_map.capacity(), KOKKOS_LAMBDA(uint32_t i)
                {
                  if( lepp_incidents_map.valid_at(i) )
                  {
                        Integer k = Kokkos::atomic_fetch_add(&global_index(0), 1);

                        lepp_incident_index(k,0) = lepp_incidents_map.key_at(i);
                        lepp_incident_index(k,1) = lepp_incidents_map.value_at(i).nodes(0);
                        lepp_incident_index(k,2) = lepp_incidents_map.value_at(i).nodes(1);

                  }
                });

                double time = timer.seconds();
                if (verbose) std::cout << "compact_map_to_view took: " << time << " seconds." << std::endl;

        }*/

        template <typename T>
        void compact_map_to_view(const UnorderedMap<T, Edge>& lepp_incidents_map,
                                 ViewMatrixType<T> lepp_incident_index,
                                 const bool value = true) {
            using namespace Kokkos;

            Timer timer;

            const Integer capacity = lepp_incidents_map.capacity();

            ViewVectorType<Integer> indices = ViewVectorType<Integer>("indices", capacity);

            parallel_for(
                capacity, KOKKOS_LAMBDA(uint32_t i) { indices[i] = lepp_incidents_map.valid_at(i); });

            exclusive_scan(0, capacity, indices);

            parallel_for(
                capacity, KOKKOS_LAMBDA(uint32_t i) {
                    if (lepp_incidents_map.valid_at(i)) {
                        Integer k = indices(i);
                        lepp_incident_index(k, 0) = lepp_incidents_map.key_at(i);

                        if (value) {
                            lepp_incident_index(k, 1) = lepp_incidents_map.value_at(i).nodes(0);
                            lepp_incident_index(k, 2) = lepp_incidents_map.value_at(i).nodes(1);
                        }
                    }
                });

            double time = timer.seconds();
            if (verbose) std::cout << "compact_map_to_view took: " << time << " seconds." << std::endl;
        }

        template <typename U>
        void compact_map_to_view(const UnorderedMap<Edge, U>& lepp_incidents_map,
                                 ViewMatrixType<U> lepp_incident_index,
                                 const bool key = false) {
            using namespace Kokkos;

            Timer timer;

            const Integer capacity = lepp_incidents_map.capacity();

            ViewVectorType<Integer> indices = ViewVectorType<Integer>("indices", capacity);

            parallel_for(
                capacity, KOKKOS_LAMBDA(uint32_t i) { indices[i] = lepp_incidents_map.valid_at(i); });

            exclusive_scan(0, capacity, indices);
            parallel_for(
                capacity, KOKKOS_LAMBDA(uint32_t i) {
                    if (lepp_incidents_map.valid_at(i)) {
                        Integer k = indices(i);
                        lepp_incident_index(k, 0) = lepp_incidents_map.value_at(i);

                        if (key) {
                            lepp_incident_index(k, 1) = lepp_incidents_map.key_at(i).nodes(0);
                            lepp_incident_index(k, 2) = lepp_incidents_map.key_at(i).nodes(1);
                        }
                    }
                });

            double time = timer.seconds();
            if (verbose) std::cout << "compact_map_to_view took: " << time << " seconds." << std::endl;
        }

        inline std::array<Integer, 2> count_lepp(const ViewVectorType<Integer> elements) {
            using namespace Kokkos;

            std::array<Integer, 2> res;

            Integer nr_elements = elements.extent(0);

            /*+1 for the incident total sum for both inclusive and exclusive sum scan at the same result*/
            ViewVectorType<Integer> index_count_ = ViewVectorType<Integer>("index_count", nr_elements + 1);

            ViewVectorType<Integer> pt_count_ = ViewVectorType<Integer>(
                "pt_count",
                nr_elements + 1);  // TODO: try to use the results of the index_count avoiding the uniqueness of the map
                                   // since the atomic lepp should avoid it.

            /*ViewVectorType<bool> lepp_occupied = ViewVectorType<bool>(
                                                                    "lepp_occupied_count", host_mesh->n_elements());*/
            lepp_occupied = ViewVectorType<bool>("lepp_occupied", host_mesh->n_elements());

            Timer timer1;

            parallel_for(nr_elements,
                         CountLepp(edge_element_map_.mapping_,
                                   mesh,
                                   elements,
                                   tree_,
                                   lepp_occupied,
                                   index_count_,
                                   pt_count_,
                                   Comb::instance().combs));

            double time1 = timer1.seconds();
            if (verbose) std::cout << "Count took: " << time1 << " seconds." << std::endl;

            Timer timer2;

            complex_inclusive_scan(1, nr_elements + 1, index_count_, pt_count_);

            double time2 = timer2.seconds();
            if (verbose) std::cout << "Scan took: " << time2 << " seconds." << std::endl;

            Timer timer3;
            auto index_subview = subview(index_count_, nr_elements);
            auto h_iac = create_mirror_view(index_subview);

            // Deep copy device view to host view.
            deep_copy(h_iac, index_subview);

            auto pt_subview = subview(pt_count_, nr_elements);
            auto h_pac = create_mirror_view(pt_subview);

            // Deep copy device view to host view.
            deep_copy(h_pac, pt_subview);

            res[0] = h_iac();
            res[1] = h_pac();

            double time3 = timer3.seconds();
            if (verbose) std::cout << "Deep copy subview took: " << time3 << " seconds." << std::endl;

            /*printf("lepp_incidents_count: %li\n", h_iac(0));
            printf("lepp_node_count: %li\n", h_pac(0));*/

            return res;
        }

        void fill_lepp(const ViewVectorType<Integer> elements,
                       UnorderedMap<Integer, Edge>& lepp_incidents_map,
                       UnorderedMap<Edge, Integer>& lepp_edge_map) {
            using namespace Kokkos;
            Timer timer;

            Integer nr_elements = elements.extent(0);
            // reserve_tree(nr_elements);

            /*ViewVectorType<bool> lepp_occupied = ViewVectorType<bool>("lepp_occupied_fill",
                            host_mesh->n_elements());*/

            parallel_for(nr_elements,
                         ScatterElem(edge_element_map_.mapping_,
                                     mesh,
                                     elements,
                                     tree_,
                                     lepp_occupied,
                                     Comb::instance().combs,
                                     lepp_incidents_map,
                                     lepp_edge_map));

            double time = timer.seconds();
            if (verbose) std::cout << "Scatter/Fill took: " << time << " seconds." << std::endl;
        }

        void resize_mesh_and_update_map(const Integer lip_size, const Integer points_count) {
            using namespace Kokkos;

            Timer timer1;

            host_mesh->resize_children(lip_size);
            host_mesh->resize_elements(2 * lip_size);
            host_mesh->resize_points(points_count);

            copy_mesh_to_device();  // only the object. No deep copy of the member views.

            double time1 = timer1.seconds();
            if (verbose) std::cout << "Resize mesh took: " << time1 << " seconds." << std::endl;
        }

        inline void copy_to_device(ViewObject<Integer> node_start_index) {
            using namespace Kokkos;

            Timer timer;

            auto h_aci = create_mirror_view(node_start_index);
            h_aci(0) = host_mesh->n_nodes();
            deep_copy(node_start_index, h_aci);

            double time = timer.seconds();
            if (verbose) std::cout << "copy_to_device took: " << time << " seconds." << std::endl;
        }

        inline void update_maps(const Integer it_count) {
            using namespace Kokkos;

            Timer timer;

            if (it_count == 0) {
                ViewVectorType<Integer> active_elems = mark_active(host_mesh->get_view_active());

                const Integer nr_active_elements = active_elems.extent(0);

                reserve_tree(nr_active_elements);
                // const Integer nr_active_elements = euler_graph_formula(host_mesh);
                // currently an approximation in the absence of the parallel dual graph need for the euler formula.
                edge_element_map_.reserve_map(4 * nr_active_elements);
                edge_element_map_.update(mesh, active_elems);
            }

            double time = timer.seconds();
            if (verbose) std::cout << "precompute_lepp_incidents took: " << time << " seconds." << std::endl;
        }

        inline void bisect_elements(const Integer lip_size,
                                    const Integer pt_size,
                                    UnorderedMap<Integer, Edge>& lepp_incidents_map,
                                    UnorderedMap<Edge, Integer>& lepp_edge_node_map,
                                    const Integer node_start_index,
                                    const Integer elem_start_index,
                                    const Integer child_start_index) {
            using namespace Kokkos;

            ViewMatrixType<Integer> lepp_incident_index = ViewMatrixType<Integer>("lepp_incidents_index", lip_size, 3);

            ViewMatrixType<Integer> lepp_point_index = ViewMatrixType<Integer>("lepp_point_index", pt_size, 1);

            compact_map_to_view(lepp_incidents_map, lepp_incident_index);
            compact_map_to_view(lepp_edge_node_map, lepp_point_index);

            Timer timer;

            parallel_for(lip_size, HandleParentElem(mesh, lepp_incident_index, elem_start_index, child_start_index));

            double time = timer.seconds();
            if (verbose) std::cout << "HandleParentElem took: " << time << " seconds." << std::endl;

            Timer timer1;

            parallel_for(
                pt_size,
                BisectEdge(
                    mesh, lepp_point_index, edge_element_map_.mapping_, lepp_edge_node_map, verbose, node_start_index));

            double time1 = timer1.seconds();
            if (verbose) std::cout << "BisectEdge took: " << time1 << " seconds." << std::endl;

            Timer timer2;

            parallel_for(lip_size,
                         BisectElem(mesh,
                                    lepp_incident_index,
                                    Comb::instance().combs,
                                    edge_element_map_.mapping_,
                                    lepp_edge_node_map,
                                    elem_start_index,
                                    lip_size));

            double time2 = timer2.seconds();
            if (verbose) std::cout << "BisectElem took: " << time2 << " seconds." << std::endl << std::endl;

            /* std::cout << "------------------------------------------------------------" << std::endl; */
        }

        inline void refine_elements(ViewVectorType<Integer>& elements) {
            using namespace Kokkos;

            Timer timer_refine;

            Integer it_count = 0;
            while (compact(mesh, elements)) {
                update_maps(it_count);

                ++it_count;

                auto h_ac = count_lepp(elements);

                UnorderedMap<Integer, Edge> lepp_incidents_map = UnorderedMap<Integer, Edge>(h_ac[0]);

                UnorderedMap<Edge, Integer> lepp_edge_node_map = UnorderedMap<Edge, Integer>(h_ac[1]);

                fill_lepp(elements, lepp_incidents_map, lepp_edge_node_map);

                const Integer lip_size = lepp_incidents_map.size();
                const Integer pt_size = lepp_edge_node_map.size();

                if (verbose) {
                    std::cout << "Lip_size: " << lip_size << std::endl;
                    std::cout << "Point_size: " << pt_size << std::endl;
                }

                // keep track of the start index before the resize of the mesh.
                const Integer node_start_index = host_mesh->n_nodes();
                const Integer elem_start_index = host_mesh->n_elements();
                const Integer child_start_index = host_mesh->n_childrens();

                host_mesh->reserve_elem_children_map(elem_start_index);
                resize_mesh_and_update_map(lip_size, pt_size);

                bisect_elements(lip_size,
                                pt_size,
                                lepp_incidents_map,
                                lepp_edge_node_map,
                                node_start_index,
                                elem_start_index,
                                child_start_index);
            }

            // free_mesh();

            double time = timer_refine.seconds();
            std::cout << "Refine Mesh took: " << time << " seconds. In " << it_count << " iterations." << std::endl;
        }

        void refine(ViewVectorType<Integer>& elements) {
            if (edge_element_map_.empty()) {
                host_mesh->set_combs(Comb::instance().combs);
                copy_mesh_to_device();
            }

            refine_elements(elements);
        }

        void uniform_refine(const Integer n_levels) {
            Kokkos::Timer timer_refine;

            if (edge_element_map_.empty()) {
                host_mesh->set_combs(Comb::instance().combs);
                copy_mesh_to_device();
            }

            for (Integer i = 0; i < n_levels; ++i) {
                ViewVectorType<Integer> elements = mark_active(host_mesh->get_view_active());

                std::cout << "\nn_marked(" << (i + 1) << "/" << n_levels << ") : " << elements.extent(0) << std::endl;

                refine_elements(elements);
            }

            double time = timer_refine.seconds();
            std::cout << "parallel Uniform LEPP took: " << time << " seconds." << std::endl;
        }

    private:
        Mesh* mesh;
        Mesh* host_mesh;
        ParallelEdgeElementMap edge_element_map_;
        ParallelEdgeNodeMap edge_node_map_;

        UnorderedMap<Integer, ElementVector> tree_;
        ViewVectorType<bool> lepp_occupied;
        // ViewVectorType<uint32_t> tree_;

        /*std::vector<Integer> level;
        std::vector<std::array<Integer, ManifoldDim+1> > side_flags;


        //tracking the refinement
        std::vector<Edge> bisected_edges_;

        std::vector<Edge>    incomplete_edges_;
        std::vector<Integer> incomplete_elements_;
        Tracker tracker_;*/
        bool verbose;
        bool fail_if_not_refine;
    };
}  // namespace mars

#endif  // MARS_BISECTION_HPP
