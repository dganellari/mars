set(MARS_PETSC_ROOT "$ENV{PETSC_DIR}" CACHE PATH "PETSc installation or source build root")
set(MARS_PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "PETSc architecture directory for an uninstalled build")
find_path(MARS_PETSC_INCLUDE_DIR petscksp.h HINTS "${MARS_PETSC_ROOT}/include")
find_path(MARS_PETSC_CONFIG_DIR petscconf.h HINTS
    "${MARS_PETSC_ROOT}/${MARS_PETSC_ARCH}/include" "${MARS_PETSC_ROOT}/include")
find_library(MARS_PETSC_LIBRARY NAMES petsc HINTS
    "${MARS_PETSC_ROOT}/${MARS_PETSC_ARCH}/lib" "${MARS_PETSC_ROOT}/lib")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MarsPETSc REQUIRED_VARS
    MARS_PETSC_INCLUDE_DIR MARS_PETSC_CONFIG_DIR MARS_PETSC_LIBRARY)
if(NOT MarsPETSc_FOUND)
    message(FATAL_ERROR "PETSc was not found; set MARS_PETSC_ROOT and optionally MARS_PETSC_ARCH")
endif()
if(MARS_PETSC_LIBRARY MATCHES "\\.a$")
    message(FATAL_ERROR "SIMPLE requires shared PETSc; configure PETSc with --with-shared-libraries=1")
endif()
if(NOT TARGET mars_petsc)
    add_library(mars_petsc UNKNOWN IMPORTED GLOBAL)
    set_target_properties(mars_petsc PROPERTIES IMPORTED_LOCATION "${MARS_PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MARS_PETSC_INCLUDE_DIR};${MARS_PETSC_CONFIG_DIR}")
endif()
