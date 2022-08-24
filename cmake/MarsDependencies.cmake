# MarsDependencies.cmake

# ##############################################################################
if(MARS_ENABLE_MPI)
  find_package(MPIExtended REQUIRED)
  # find_package(MPI REQUIRED)

  if(MPI_FOUND)
    set(MARS_ENABLE_MPI ON)

    if(MPI_C_INCLUDE_PATH)
      set(MARS_DEP_INCLUDES "${MARS_DEP_INCLUDES};${MPI_C_INCLUDE_PATH}")
    endif()

    if(MPI_CXX_INCLUDE_PATH)
      set(MARS_DEP_INCLUDES "${MARS_DEP_INCLUDES};${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(MPI_LIBRARIES)
      set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${MPI_LIBRARIES}")
    endif()

    if(MPI_C_LIBRARIES)
      set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${MPI_C_LIBRARIES}")
    endif()

    if(MPI_CXX_LIBRARIES)
      set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${MPI_CXX_LIBRARIES}")
    endif()
  else()
    message(
      FATAL_ERROR
        "We should never end up here, because find_package above is REQUIRED")
  endif()
endif()

# MPI

# ##############################################################################

# Kokkos
cmake_dependent_option(MARS_ENABLE_CUDA "Enable CUDA backend for Kokkos"
                       OFF "MARS_ENABLE_KOKKOS" OFF)
cmake_dependent_option(MARS_ENABLE_HIP "Enable HIP backend for Kokkos"
                       OFF "MARS_ENABLE_KOKKOS" OFF)
cmake_dependent_option(
  MARS_ENABLE_CUDAUVM "Kokkos use as default memory space CUDAUVM" OFF
  "MARS_ENABLE_KOKKOS;MARS_ENABLE_CUDA" OFF)
cmake_dependent_option(MARS_ENABLE_KOKKOS_KERNELS "Enable Kokkos Kernels" ON
                       "MARS_ENABLE_KOKKOS" OFF)

if(MARS_ENABLE_KOKKOS)
  message(STATUS "Setup Kokkos")
  list(APPEND CMAKE_MESSAGE_INDENT "${MARS_CMAKE_INDENT}")

  if(NOT TRILINOS_DIR)
    message(DEBUG "Setting TRILINOS_DIR to $ENV{TRILINOS_DIR}")
    set(TRILINOS_DIR
        $ENV{TRILINOS_DIR}
        CACHE PATH "Directory where Kokkos is installed")
  endif()

  if(NOT KOKKOS_DIR)
    message(DEBUG "Setting KOKKOS_DIR to $ENV{KOKKOS_DIR}")
    set(KOKKOS_DIR
        $ENV{KOKKOS_DIR}
        CACHE PATH "Directory where Kokkos is installed")
  endif()

  if(WIN32)
    find_package(Kokkos HINTS C:/projects/installations/kokkos/lib/cmake/Kokkos
                 ${Kokkos_DIR} $ENV{KOKKOS_DIR} REQUIRED)
  else()
    find_package(
      Kokkos
      HINTS
      ${KOKKOS_DIR}
      ${KOKKOS_DIR}/lib/CMake/Kokkos
      ${KOKKOS_DIR}/lib64/CMake/Kokkos
      ${TRILINOS_DIR}
      ${TRILINOS_DIR}/lib/cmake/Kokkos
      ${TRILINOS_DIR}/lib64/cmake/Kokkos
      REQUIRED)
  endif()


  message(VERBOSE "Found Kokkos")
  # set(MARS_ENABLE_KOKKOS ON)
  # check what was found
  message(STATUS "Kokkos_CXX_FLAGS: ${Kokkos_CXX_FLAGS}")
  message(STATUS "Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
  message(STATUS "Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
  message(STATUS "Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
  message(STATUS "Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
  message(STATUS "Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")


  # if the cxx compiler was not set at command line set it to the kokkos compiler
  if(Kokkos_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_SET_EXTERNALLY)
    message(STATUS "[Status] Setting CMAKE_CXX_COMPILER=${Kokkos_CXX_COMPILER}")
    set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER} CACHE FILEPATH "CXX nvcc compiler" FORCE)
    set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${Kokkos_C_COMPILER})
  endif()

  # _KK_TARGET is set as a local variable do not use outside this file
  set(_KK_TARGET "Kokkos::kokkos")

  if(Kokkos_ENABLE_OPENMP)
    set(_openmp "-fopenmp")
    # we need to be sure that all targets link against opemp
    add_link_options(${_openmp})
  endif()

  if(NOT TARGET ${_KK_TARGET})
    message(DEBUG "Kokkos target is not defined")
    add_library(${_KK_TARGET} INTERFACE IMPORTED)
    set_property(
      TARGET ${_KK_TARGET}
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Kokkos_INCLUDE_DIRS}
               ${Kokkos_TPL_INCLUDE_DIRS})
    set_property(
      TARGET ${_KK_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES
                                    ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES})
    set_property(TARGET ${_KK_TARGET} PROPERTY INTERFACE_LINK_DIRECTORIES
                                               ${Kokkos_LIBRARY_DIRS})
    set_property(TARGET ${_KK_TARGET} PROPERTY INTERFACE_COMPILE_OPTIONS
                                               ${_openmp})
  else()
    message(DEBUG "Kokkos target is defined")
  endif()

  # Check what the (imported) target does
  get_target_property(Kokkos_INTERFACE_COMPILE_OPTIONS ${_KK_TARGET}
                      INTERFACE_COMPILE_OPTIONS)
  message(
    DEBUG
    "Kokkos_INTERFACE_COMPILE_OPTIONS: ${Kokkos_INTERFACE_COMPILE_OPTIONS}")
  get_target_property(Kokkos_INTERFACE_LINK_LIBRARIES ${_KK_TARGET}
                      INTERFACE_LINK_LIBRARIES)
  message(DEBUG
          "Kokkos_INTERFACE_LINK_LIBRARIES: ${Kokkos_INTERFACE_LINK_LIBRARIES}")
  get_target_property(Kokkos_INTERFACE_INCLUDE_DIRECTORIES ${_KK_TARGET}
                      INTERFACE_INCLUDE_DIRECTORIES)
  message(
    DEBUG
    "Kokkos_INTERFACE_INCLUDE_DIRECTORIES: ${Kokkos_INTERFACE_INCLUDE_DIRECTORIES}"
  )

  # perhaps later we can attach this to the target
  # add_compile_definitions("MARS_ENABLE_KOKKOS")

  if(MARS_ENABLE_CUDA)
    if(NOT DEFINED Kokkos_ENABLE_CUDA OR NOT ${Kokkos_ENABLE_CUDA})
      message(
        FATAL_ERROR
        "Enable Kokkos CUDA or unset MARS_ENABLE_CUDA to continue with OpenMP!")
    endif()
    message(VERBOSE "Kokkos CUDA Enabled = ${Kokkos_ENABLE_CUDA}")
    # target_compile_definitions(${_KK_TARGET} INTERFACE
    # MARS_ENABLE_CUDA)
    # add_compile_definitions("MARS_ENABLE_CUDA")
    kokkos_check(OPTIONS CUDA_LAMBDA)

    # get cuda flags from the wrapper alternatively we can strip
    # Kokkos_INTERFACE_COMPILE_OPTIONS when defined
    execute_process(
      COMMAND ${Kokkos_CXX_COMPILER} --show
      OUTPUT_VARIABLE _wrapper_command
      ERROR_QUIET)
    string(REGEX REPLACE [[\n\v\c\c]] "" _wrapper_flags ${_wrapper_command})
    string(STRIP "${_wrapper_flags}" _wrapper_flags)
    message(DEBUG "_wrapper_flags ${_wrapper_flags}")

    # this could be done per target if we need to compile other parts of QuICC
    # with different CUDA settings
    set(CMAKE_CUDA_FLAGS "${_wrapper_flags} ${_openmp}")

    elseif(MARS_ENABLE_HIP)
        if(NOT DEFINED Kokkos_ENABLE_HIP OR NOT ${Kokkos_ENABLE_HIP})
      message(
        FATAL_ERROR
        "Enable Kokkos HIP or unset MARS_ENABLE_HIP to continue with OpenMP!")
        endif()
    message(VERBOSE "Kokkos HIP Enabled = ${Kokkos_ENABLE_HIP}")

  else()
    string(FIND "${CMAKE_CXX_FLAGS}" "${_openmp}" _pos)
    if(_pos EQUAL -1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_openmp}")
    endif()

  endif()

  get_target_property(Kokkos_INTERFACE_LINK_LIBRARIES ${_KK_TARGET}
                      INTERFACE_LINK_LIBRARIES)
  message(DEBUG
          "Kokkos_INTERFACE_LINK_LIBRARIES: ${Kokkos_INTERFACE_LINK_LIBRARIES}")

  add_library(mars::Kokkos INTERFACE IMPORTED)
  set_target_properties(
    mars::Kokkos
    PROPERTIES INTERFACE_LINK_LIBRARIES "${Kokkos_INTERFACE_LINK_LIBRARIES}"
               INTERFACE_INCLUDE_DIRECTORIES
               "${Kokkos_INTERFACE_INCLUDE_DIRECTORIES}"
               INTERFACE_COMPILE_OPTIONS "${_openmp}")

  target_compile_definitions(mars::Kokkos INTERFACE COMPILE_FOR_KOKKOS)

  set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};Kokkos::kokkos")

  # done with setting up Kokkos target
  unset(_KK_TARGET)

  # list(POP_BACK CMAKE_MESSAGE_INDENT)

  #
  # to be separated
  #

  message(STATUS "Setup Kokkos Kernels")
  list(APPEND CMAKE_MESSAGE_INDENT "${MARS_CMAKE_INDENT}")

  # Kokkos Kernels
  if(MARS_ENABLE_KOKKOS_KERNELS)

    if(WIN32)
      find_package(
        KokkosKernels HINTS
        C:/projects/installations/kokkos-kernels/lib/cmake/KokkosKernels
        ${KokkosKernels_DIR} $ENV{KokkosKernels_DIR} REQUIRED)
    else()
      find_package(
        KokkosKernels
        HINTS
        ${KOKKOS_DIR}
        ${KOKKOS_DIR}/lib/CMake/KokkosKernels
        ${KOKKOS_DIR}/lib64/cmake/KokkosKernels
        ${TRILINOS_DIR}
        ${TRILINOS_DIR}/lib/cmake/KokkosKernels
        ${TRILINOS_DIR}/lib64/cmake/KokkosKernels
        REQUIRED)
    endif()

    message(VERBOSE "Found Kokkos Kernels")
    # set(MARS_ENABLE_KOKKOS_KERNELS ON)
    message(VERBOSE "KokkosKernels_LIBRARIES = ${KokkosKernels_LIBRARIES}")
    message(VERBOSE
            "KokkosKernels_LIBRARY_DIRS = ${KokkosKernels_LIBRARY_DIRS}")

    # _KKK_TARGET is set as a local variable do not use outside this file
    set(_KKK_TARGET "Kokkos::kokkoskernels")
    if(NOT TARGET ${_KKK_TARGET})
      message(DEBUG "Kokkos kernel target is not defined")
      add_library(${_KKK_TARGET} INTERFACE IMPORTED)
      set_property(
        TARGET ${_KKK_TARGET}
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${KokkosKernels_INCLUDE_DIRS}
                 ${KokkosKernels_TPL_INCLUDE_DIRS})
      set_property(
        TARGET ${_KKK_TARGET}
        PROPERTY INTERFACE_LINK_LIBRARIES ${KokkosKernels_LIBRARIES}
                 ${KokkosKernels_TPL_LIBRARIES})
      set_property(TARGET ${_KKK_TARGET} PROPERTY INTERFACE_LINK_DIRECTORIES
                                                  ${KokkosKernels_LIBRARY_DIRS})

      set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${KokkosKernels_LIBRARIES}")
    else()
        set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};Kokkos::kokkoskernels")
    endif()



    # done with setting up Kokkos Kernels target
    unset(_KKK_TARGET)
  endif()

  # list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

# ##############################################################################
# Moonolith

if(MARS_ENABLE_MOONOLITH)
  add_subdirectory(moonolith_adapter)
endif()

if(MARS_ENABLE_CXXOPTS)
  include(cxxopts/cxxopts.cmake)
  # set(MARS_ENABLE_CXXOPTS ON)
endif()

if(MARS_ENABLE_ADIOS2)
  find_package(ADIOS2 REQUIRED)
  if(ADIOS2_FOUND)
    message(STATUS "Adios2 found.")
    # set(MARS_ENABLE_ADIOS2 ON)
    add_subdirectory(backend/adios2)

    set(ADIOS2_INCLUDE_DIRS ${ADIOS2_DIR}/../../../include)
    set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${ADIOS2_LIBRARIES}")
    set(MARS_DEP_INCLUDES "${MARS_DEP_INCLUDES};${ADIOS2_INCLUDE_DIRS}")
  else()
    message(STATUS "Adios2 not found.")
  endif()
endif()

