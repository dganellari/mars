/*
 * mars_edge_element_map_kokkos.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: gandanie
 */


#include "mars_simplex.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_mesh.hpp"
#include "mars_bisection.hpp"
#include "mars_vtk_writer.hpp"
#include "mars_quality.hpp"
#include "mars_utils.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_benchmark.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "generation/mars_lepp_benchmark_kokkos.hpp"

#include "mars_test.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_longest_edge.hpp"
#include "generation/mars_memory.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"

template<>
mars::Combinations<4, 2, mars::KokkosImplementation>
mars::ParallelSubManifoldElementMap<2l, 3l>::combinations = mars::Combinations<4, 2, mars::KokkosImplementation>();
