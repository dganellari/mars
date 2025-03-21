option(MARS_ENABLE_MPI "Build mars with MPI support" ON)

option(MARS_ENABLE_KOKKOS "Use -DMARS_ENABLE_KOKKOS=ON for enabling mesh transfer functions." ON)
option(MARS_ENABLE_KOKKOS_KERNELS "Use -DMARS_ENABLE_KOKKOS_KERNELS=ON for enabling kokkos_kernels" ON)
option(MARS_ENABLE_CUDA "Build mars with cuda support" OFF)
option(MARS_ENABLE_CUDAUVM "Use as default memory space CUDAUVM. This has only an effect if MARS_ENABLE_CUDA is used too." OFF)
option(MARS_ENABLE_RDMA "Build mars with RDMA support" ON)

option(MARS_ENABLE_HIP "Build mars with cuda support" OFF)

option(MARS_ENABLE_ADIOS2 "Uses ADIOS2 IO" OFF)

# maintainace
option(MARS_ENABLE_BENCHMARK "Enable benchmarks" OFF)
option(MARS_ENABLE_TESTS "Enable tests" ON)

# FIXME Extra backends 
option(MARS_ENABLE_MOONOLITH "Use -DMARS_ENABLE_MOONOLITH=ON for enabling mesh transfer functions." OFF)
option(MARS_ENABLE_VTK "Uses VTK IO" OFF)

# FIXME  and turn on!
option(MARS_ENABLE_SERIAL_BACKEND "Use -DMARS_ENABLE_SERIAL_BACKEND=ON for enabling the serial backend" OFF)
option(MARS_ENABLE_DISTRIBUTED_BACKEND "Use -DMARS_ENABLE_DISTRIBUTED_BACKEND=ON for enabling the distributed backend" ON)
option(MARS_ENABLE_AMR_BACKEND "Use -DMARS_ENABLE_AMR_BACKEND=ON for enabling the AMR backend" OFF)

option(MARS_ENABLE_EXAMPLES "Enable examples" OFF)
option(MARS_ENABLE_CXXOPTS "Enable cxxopts" ON)

# New options for unstructured and octree modules
option(MARS_ENABLE_UNSTRUCTURED "Enable unstructured module" ON)
option(MARS_ENABLE_OCTREE "Enable octree module" ON)