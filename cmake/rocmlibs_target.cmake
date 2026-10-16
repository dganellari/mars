find_package(hip REQUIRED)
# Thrust, CUB, the solvers' SpMV and BLAS on AMD (backend/distributed/unstructured/utils/mars_gpu_*.hpp)
find_package(rocthrust REQUIRED)
find_package(hipcub REQUIRED)
find_package(hipsparse REQUIRED)
find_package(hipblas REQUIRED)

# Only include paths and host libraries, as in cornerstone: the rocthrust/hipcub targets carry
# hip::device, which puts --offload-arch on every C++ compile and link line, and GCC rejects it.
if (NOT TARGET mars::rocmlibs)
  add_library(mars::rocmlibs INTERFACE IMPORTED)
  foreach(lib roc::rocthrust hip::hipcub roc::rocprim roc::hipsparse roc::hipblas)
    get_property(inc TARGET ${lib} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(mars::rocmlibs INTERFACE ${inc})
  endforeach()
  target_link_libraries(mars::rocmlibs INTERFACE hip::host $<TARGET_FILE:roc::hipsparse> $<TARGET_FILE:roc::hipblas>)
endif()
