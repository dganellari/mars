# First with try with clang compile
find_library(
    MPI_C_CLANG_LIBRARY
    NAMES mpi mpich
    PATHS /opt/local/lib/openmpi-mp/ /opt/local/lib/mpich-mp/
          /opt/local/lib/mpich-clang/ /opt/local/lib
    DOC "The MPI_CXX_CLANG_LIBRARY library to link against")

find_library(
    MPI_CXX_CLANG_LIBRARY
    NAMES mpi_cxx mpicxx-mpich-clang mpichcxx mpicxx
    PATHS /opt/local/lib/openmpi-mp/ /opt/local/lib/mpich-mp/
          /opt/local/lib/mpich-clang/ /opt/local/lib
    DOC "The MPI_CXX_CLANG_LIBRARY library to link against")

if(MPI_CXX_CLANG_LIBRARY)

    set(MPI_CXX_LIBRARIES ${MPI_CXX_CLANG_LIBRARY})

    get_filename_component(MPI_LIB_DIR ${MPI_CXX_CLANG_LIBRARY} PATH)

    find_path(
        MPI_CLANG_HEADERS mpi.h
        HINTS ${MPI_LIB_DIR}/../../include
              ${MPI_LIB_DIR}/../include
              /opt/local/include/openmpi-mp/
              /opt/local/include/mpich-mp
              ${MPI_LIB_DIR}/../../include/openmpi-mp/
              ${MPI_LIB_DIR}/../../include/mpich-clang
              ${MPI_LIB_DIR}/../include/mpich-clang
              /opt/local/include/mpich-clang
        DOC "The MPI_CLANG_HEADERS path")

    if(MPI_CLANG_HEADERS)
        find_file(
            MPI_CXX_COMPILER
            NAMES mpicxx mpicxx-openmpi-mp mpicxx-mpich-clang
            HINTS ${MPI_CLANG_HEADERS}/../bin ${MPI_CLANG_HEADERS}/../../bin
                  ${MPI_LIB_DIR}/../bin ${MPI_LIB_DIR}/../../bin /opt/local/bin/
            DOC "the MPI_COMPILER_PATH dir")

        find_file(
            MPI_C_COMPILER
            NAMES mpicc mpicc-openmpi-mp mpicc-mpich-clang
            HINTS ${MPI_CLANG_HEADERS}/../bin ${MPI_CLANG_HEADERS}/../../bin
                  ${MPI_LIB_DIR}/../bin ${MPI_LIB_DIR}/../../bin /opt/local/bin/
            DOC "the MPI_COMPILER_PATH dir")

        if(MPI_CXX_COMPILER AND MPI_C_COMPILER)
            # set variables
            set(MPI_FOUND TRUE)
            set(MPI_C_LIBRARIES ${MPI_C_CLANG_LIBRARY})
            set(MPI_CXX_LIBRARIES ${MPI_CXX_CLANG_LIBRARY})

            set(MPI_C_INCLUDE_PATH ${MPI_CLANG_HEADERS})
            set(MPI_CXX_INCLUDE_PATH ${MPI_CLANG_HEADERS})

            message(STATUS "MPI: ${MPI_C_LIBRARIES}; ${MPI_CXX_LIBRARIES}")

            include_directories(${MPI_C_INCLUDE_PATH})
            link_libraries(${MPI_C_LIBRARIES})
            link_libraries(${MPI_CXX_LIBRARIES})

        endif()
    endif()
endif()

# MESSAGE(STATUS "${MPI_CXX_CLANG_LIBRARY} ${MPI_CLANG_HEADERS}
# ${MPI_CXX_COMPILER}")

if(NOT MPI_FOUND)
    find_package(MPI)
endif()
