option(MARS_ENABLE_VTK "Uses VTK IO" OFF)
option(MARS_ENABLE_ADIOS2 "Uses ADIOS2 IO" OFF)
option(MARS_ENABLE_BENCHMARK "Enable benchmarks" ON)
option(MARS_ENABLE_CXXOPTS "Enable cxxopts" OFF)
option(MARS_ENABLE_TESTING "Enable tests" ON)
option(MARS_ENABLE_RDMA "Build mars with RDMA support" ON)
option(MARS_ENABLE_MPI "Build mars with MPI support" ON)
option(MARS_ENABLE_CUDA "Build mars with cuda support" OFF)
option(MARS_ENABLE_CUDAUVM "Use as default memory space CUDAUVM. This has only an effect if MARS_ENABLE_CUDA is used too." OFF)

option(MARS_ENABLE_MOONOLITH "Use -DMARS_ENABLE_MOONOLITH=ON for enabling mesh transfer functions." OFF)
option(MARS_ENABLE_KOKKOS "Use -DMARS_ENABLE_KOKKOS=ON for enabling mesh transfer functions." ON)
option(MARS_ENABLE_KOKKOS_KERNELS "Use -DMARS_ENABLE_KOKKOS_KERNELS=ON for enabling kokkos_kernels" ON)

# Fixme and turn on!
option(MARS_ENABLE_SERIAL_BACKEND "Use -DMARS_ENABLE_SERIAL_BACKEND=ON for enabling serial backend" OFF)