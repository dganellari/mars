#ifndef MARS_HPP
#define MARS_HPP

#include "mars_base.hpp"

#include "mars_bisection.hpp"
#include "mars_mesh.hpp"

#include "mars_benchmark.hpp"
#include "mars_dual_graph.hpp"
#include "mars_edge_node_map.hpp"
#include "mars_edge_split.hpp"
#include "mars_edge_split_pool.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_longest_edge.hpp"
#include "mars_matrix.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_newest_vertex.hpp"
#include "mars_node_rank.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_quality.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_red_green_refinement.hpp"
#include "mars_static_math.hpp"
#include "mars_test.hpp"
#include "mars_utils.hpp"
#include "mars_vector.hpp"
#include "mars_visualization.hpp"
#include "mars_vtk_writer.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_longest_edge_kokkos.hpp"
#include "mars_mesh_kokkos.hpp"
#include "mars_utils_kokkos.hpp"
#endif  // WITH_KOKKOS

#ifdef WITH_MPI
#include "mars_mpi_guard.hpp"
#ifdef WITH_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_mesh_generation.hpp"
#include "mars_distributed_sparsity_pattern.hpp"
#include "mars_distributed_staggered_data_management.hpp"
#include "mars_distributed_user_data.hpp"
#endif  // WITH_KOKKOS
#endif  // WITH_MPI
#endif  // MARS_HPP
