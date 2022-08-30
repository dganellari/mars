#ifndef MARS_HPP
#define MARS_HPP

#include "mars_base.hpp"

#ifdef MARS_ENABLE_KOKKOS
#include "mars_utils_kokkos.hpp"
#ifdef MARS_ENABLE_AMR_BACKEND
#include "mars_mesh_kokkos.hpp"
#include "mars_lepp_benchmark_kokkos.hpp"
#endif  // MARS_ENABLE_AMR
#endif  // MARS_ENABLE_KOKKOS

#ifdef MARS_ENABLE_MPI
#include "mars_mpi_guard.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_mesh_generation.hpp"
#include "mars_distributed_sparsity_pattern.hpp"
#include "mars_distributed_staggered_data_management.hpp"
#include "mars_distributed_user_data.hpp"
#endif  // MARS_ENABLE_KOKKOS
#endif  // MARS_ENABLE_MPI
#endif  // MARS_HPP
