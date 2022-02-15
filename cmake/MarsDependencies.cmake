# MarsDependencies.cmake

# ##############################################################################

find_package(MPIExtended REQUIRED)

if(MPI_FOUND)
    set(WITH_MPI ON)

    if(MPI_C_INCLUDE_PATH)
        set(MARS_DEP_INCLUDES
            "${MARS_DEP_INCLUDES};${MPI_C_INCLUDE_PATH}")
    endif()

    if(MPI_CXX_INCLUDE_PATH)
        set(MARS_DEP_INCLUDES
            "${MARS_DEP_INCLUDES};${MPI_CXX_INCLUDE_PATH}")
    endif()

    if(MPI_LIBRARIES)
        set(MARS_DEP_LIBRARIES
            "${MARS_DEP_LIBRARIES};${MPI_LIBRARIES}")
    endif()

    if(MPI_C_LIBRARIES)
        set(MARS_DEP_LIBRARIES
            "${MARS_DEP_LIBRARIES};${MPI_C_LIBRARIES}")
    endif()

    if(MPI_CXX_LIBRARIES)
        set(MARS_DEP_LIBRARIES
            "${MARS_DEP_LIBRARIES};${MPI_CXX_LIBRARIES}")
    endif()
else()
    message(
        FATAL_ERROR
            "We should never end up here, because find_package above is REQUIRED"
    )
endif()

#MPI

# ##############################################################################

#Kokkos

if(NOT TRILINOS_DIR)
        message(STATUS "Setting TRILINOS_DIR to $ENV{TRILINOS_DIR}")
        set(TRILINOS_DIR
            $ENV{TRILINOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    if(NOT KOKKOS_DIR)
        message(STATUS "Setting KOKKOS_DIR to $ENV{KOKKOS_DIR}")
        set(KOKKOS_DIR
            $ENV{KOKKOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    # FIND_PACKAGE(Trilinos PATHS ${TRILINOS_DIR}/lib/cmake/Trilinos QUIET)

    find_package(
        Kokkos
        HINTS
        ${KOKKOS_DIR}
        ${KOKKOS_DIR}/lib/CMake/Kokkos
        ${KOKKOS_DIR}/lib64/CMake/Kokkos
        ${TRILINOS_DIR}
        ${TRILINOS_DIR}/lib/cmake/Kokkos
        ${TRILINOS_DIR}/lib64/cmake/Kokkos
        QUIET)

    if(Kokkos_FOUND)
        if (TARGET Kokkos::kokkos)
            message(STATUS "Kokkos::kokkos is a target, get include directories from the target")
            get_target_property(Kokkos_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
            get_target_property(Kokkos_LIBRARIES Kokkos::kokkos INTERFACE_LINK_LIBRARIES)
            get_target_property(Kokkos_LIBRARY_DIRS Kokkos::kokkos INTERFACE_LINK_DIRECTORIES)
        endif()
        message("\nFound Kokkos!  Here are the details: ")
        message(" Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
        message(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
        message(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
        message(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
        message(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")

        if(MARS_USE_CUDA)
            message(" Kokkos CUDA Enabled = ${Kokkos_ENABLE_CUDA}")
        endif()

        # message("   Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}") message("
        # Kokkos_C_COMPILER = ${Kokkos_C_COMPILER}")

        if(Kokkos_CXX_COMPILER)
            set(CMAKE_C_COMPILER ${Kokkos_C_COMPILER})
            set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
            # message( "Setting CMAKE_CXX_COMPILER to
            # Kokkos_CXX_COMPILER=${Kokkos_CXX_COMPILER}" )
        endif()

        set(WITH_KOKKOS ON)

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

        find_package(OpenMP)

        if(OPENMP_FOUND)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_LINK_FLAGS
                "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()

    else()
        message(WARNING "Could not find Kokkos!")
        if(MARS_USE_CUDA)
            message(FATAL_ERROR "Could not use CUDA without Kokkos!")
        endif()
endif()

if(Kokkos_FOUND)

    message(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
    message(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
    message(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
    message(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")

    if(MARS_USE_CUDA AND (NOT DEFINED Kokkos_ENABLE_CUDA OR NOT ${Kokkos_ENABLE_CUDA}))
            message(
                FATAL_ERROR
                    "Enable Kokkos Cuda or unset MARS_USE_CUDA to continue with OpenMP!"
            )
    endif()

    if(Kokkos_ENABLE_CUDA)
        kokkos_check(OPTIONS CUDA_LAMBDA)
    endif()

    target_include_directories(mars SYSTEM PUBLIC ${Kokkos_TPL_INCLUDE_DIRS} ${Kokkos_INCLUDE_DIRS})

    message(
        STATUS
            "Kokkos_INCLUDE_DIRS=${Kokkos_INCLUDE_DIRS}, ${Kokkos_TPL_INCLUDE_DIRS}"
    )

    if(Kokkos_CXX_COMPILER)
        target_link_libraries(mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES})
    else()
        target_link_libraries(mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES}
                              -L${Kokkos_LIBRARY_DIRS})
    endif()

    if(KokkosKernels_FOUND)
        if (TARGET Kokkos::kokkoskernels)
            target_link_libraries(mars Kokkos::kokkoskernels)
        else()
            target_include_directories(mars PUBLIC ${KokkosKernels_TPL_INCLUDE_DIRS}
                                                           ${KokkosKernels_INCLUDE_DIRS})

            if(Kokkos_CXX_COMPILER)
                target_link_libraries(mars ${KokkosKernels_LIBRARIES}
                                      ${KokkosKernels_TPL_LIBRARIES})
            else()
                target_link_libraries(
                    mars ${KokkosKernels_LIBRARIES} ${KokkosKernels_TPL_LIBRARIES}
                    -L${KokkosKernels_LIBRARY_DIRS})
            endif()
        endif()
    endif()

endif()
