# MarsDependencies.cmake

# ##############################################################################
if(MARS_ENABLE_MPI AND MARS_ENABLE_KOKKOS)
  find_package(MPI REQUIRED COMPONENTS CXX)
  # Keep compiler flags and install-time dependency lookup attached to one target.
  list(APPEND MARS_DEP_LIBRARIES MPI::MPI_CXX)
elseif(MARS_ENABLE_MPI)
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

if(MARS_ENABLE_KOKKOS)
  find_package(Kokkos CONFIG REQUIRED
      HINTS "${KOKKOS_DIR}" "$ENV{KOKKOS_DIR}" "${TRILINOS_DIR}" "$ENV{TRILINOS_DIR}"
      PATH_SUFFIXES lib/cmake/Kokkos lib64/cmake/Kokkos)
  include(${CMAKE_CURRENT_LIST_DIR}/MarsKokkos.cmake)
  mars_check_kokkos_backend("${MARS_GPU_LANGUAGE}")
  # Kokkos exports the flags for its compiler language and enabled host backends.
  list(APPEND MARS_DEP_LIBRARIES Kokkos::kokkos)

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
if(MARS_ENABLE_HYPRE)
  message(STATUS "Setup Hypre")
  list(APPEND CMAKE_MESSAGE_INDENT "${MARS_CMAKE_INDENT}")

  # Try to find Hypre using pkg-config first
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(HYPRE_PC QUIET hypre)
  endif()

  if(HYPRE_PC_FOUND)
    # Found via pkg-config
    set(HYPRE_FOUND TRUE)
    set(HYPRE_INCLUDE_DIRS ${HYPRE_PC_INCLUDE_DIRS})
    set(HYPRE_LIBRARIES ${HYPRE_PC_LIBRARIES})
    set(HYPRE_LIBRARY_DIRS ${HYPRE_PC_LIBRARY_DIRS})
  else()
    # Try manual finding
    find_path(HYPRE_INCLUDE_DIR
      NAMES HYPRE.h
      PATHS $ENV{HYPRE_DIR}/include ${HYPRE_DIR}/include
      PATH_SUFFIXES hypre
    )

    find_library(HYPRE_LIBRARY
      NAMES HYPRE hypre
      PATHS $ENV{HYPRE_DIR} ${HYPRE_DIR}
      PATH_SUFFIXES lib lib64
    )

    if(HYPRE_INCLUDE_DIR AND HYPRE_LIBRARY)
      set(HYPRE_FOUND TRUE)
      set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
      set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
    endif()
  endif()

  if(NOT HYPRE_FOUND)
    message(FATAL_ERROR "Hypre not found. Please install Hypre and set HYPRE_DIR or ensure pkg-config can find it.")
  endif()

  message(VERBOSE "Found Hypre")
  message(STATUS "HYPRE_INCLUDE_DIRS = ${HYPRE_INCLUDE_DIRS}")
  message(STATUS "HYPRE_LIBRARIES = ${HYPRE_LIBRARIES}")

  set(MARS_DEP_INCLUDES "${MARS_DEP_INCLUDES};${HYPRE_INCLUDE_DIRS}")
  set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};${HYPRE_LIBRARIES}")

  list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

# ##############################################################################
