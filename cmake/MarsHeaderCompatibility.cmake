# Keep source-relative includes working while each header has one implementation.
function(mars_compatibility_headers target)
    set(relocations
        "kokkos|backend/kokkos/common"
        "backend/amr|backend/kokkos/bisection"
        "backend/distributed/mesh|backend/distributed/structured/mesh"
        "backend/distributed/mesh_manager|backend/distributed/structured/mesh_manager"
        "backend/distributed/sfc|backend/distributed/structured/sfc")
    set(compat "${CMAKE_CURRENT_BINARY_DIR}/compat/include")
    foreach(header IN LISTS MARS_HEADERS)
        file(RELATIVE_PATH canonical "${CMAKE_CURRENT_SOURCE_DIR}" "${header}")
        foreach(relocation IN LISTS relocations)
            string(REPLACE "|" ";" paths "${relocation}")
            list(GET paths 0 previous)
            list(GET paths 1 current)
            if(canonical MATCHES "^${current}/(.+)$")
                set(old_header "${previous}/${CMAKE_MATCH_1}")
                get_filename_component(old_dir "${old_header}" DIRECTORY)
                file(MAKE_DIRECTORY "${compat}/${old_dir}")
                file(WRITE "${compat}/${old_header}"
                    "#pragma once\n#include \"${canonical}\"\n")
                install(FILES "${compat}/${old_header}" DESTINATION "include/${old_dir}")
            endif()
        endforeach()
    endforeach()
    target_include_directories(${target} PUBLIC "$<BUILD_INTERFACE:${compat}>")
endfunction()
