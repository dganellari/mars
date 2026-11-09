function(mars_check_kokkos_backend gpu_language)
    if(NOT TARGET Kokkos::kokkos)
        message(FATAL_ERROR "MARS requires a Kokkos package exporting Kokkos::kokkos.")
    endif()
    foreach(backend CUDA HIP)
        if(Kokkos_ENABLE_${backend} AND NOT gpu_language STREQUAL backend)
            message(FATAL_ERROR
                "The Kokkos package enables ${backend}, but MARS selects '${gpu_language}'. "
                "Use -DMARS_ENABLE_${backend}=ON with this package, or select a matching Kokkos_DIR.")
        endif()
    endforeach()
    if(gpu_language)
        if(NOT Kokkos_ENABLE_${gpu_language})
            message(FATAL_ERROR
                "MARS selects ${gpu_language}, but the Kokkos package does not enable it. "
                "Select a Kokkos_DIR built with -DKokkos_ENABLE_${gpu_language}=ON.")
        endif()
        if(NOT Kokkos_COMPILE_LANGUAGE STREQUAL gpu_language)
            message(FATAL_ERROR
                "MARS GPU sources use CMake's native ${gpu_language} language, but Kokkos exports "
                "'${Kokkos_COMPILE_LANGUAGE}'. Rebuild Kokkos with "
                "-DKokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON and -DKokkos_ENABLE_${gpu_language}=ON; "
                "do not change CMAKE_CXX_COMPILER after configuring MARS.")
        endif()
    endif()
endfunction()
