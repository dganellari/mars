# MarsClangFormat.cmake - Adds convenient formatting targets for C/C++ files
# Two functions are provided:
#   mars_format_add_sources - Format individual source files
#   mars_format_add_target - Format all sources in specified targets

# Find clang-format once and cache the result
if(NOT CLANG_FORMAT_EXECUTABLE)
    # Allow user-specified override
    if(CLANG_FORMAT_EXE)
        set(CLANG_FORMAT_EXECUTABLE ${CLANG_FORMAT_EXE} CACHE FILEPATH "Path to clang-format executable")
    else()
        # Try to find clang-format with a preference for newer versions
        find_program(CLANG_FORMAT_EXECUTABLE
            NAMES
                clang-format
                clang-format-16 clang-format-15 clang-format-14 clang-format-13
                clang-format-12 clang-format-11 clang-format-10 clang-format-9
                clang-format-mp-16 clang-format-mp-15 clang-format-mp-14
                clang-format-mp-13
            DOC "Path to clang-format executable"
            HINTS
                ${CLANG_FORMAT_ROOT_DIR}/bin
                /usr/bin
                /usr/local/bin
                /opt/local/bin
        )
    endif()
    # Provide information about the found executable
    if(CLANG_FORMAT_EXECUTABLE)
        execute_process(
            COMMAND ${CLANG_FORMAT_EXECUTABLE} --version
            OUTPUT_VARIABLE CLANG_FORMAT_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        message(STATUS "Found clang-format: ${CLANG_FORMAT_EXECUTABLE} (${CLANG_FORMAT_VERSION})")
    else()
        message(STATUS "clang-format not found. Code formatting targets will be disabled.")
    endif()
endif()

# Create main format target if it doesn't exist
if(NOT TARGET format AND CLANG_FORMAT_EXECUTABLE)
    add_custom_target(format)
    message(STATUS "Created main 'format' target")
endif()

# Set up formatting for source files
function(mars_format_add_sources)
    if(NOT CLANG_FORMAT_EXECUTABLE)
        return()
    endif()

    # Process each source argument
    set(format_sources "")
    foreach(file ${ARGN})
        # Get absolute path
        get_filename_component(file_abs ${file} ABSOLUTE)
        # Only add C/C++ source files
        get_filename_component(file_ext ${file} EXT)
        if(file_ext MATCHES "\\.(c|cpp|cxx|cc|h|hpp|hxx|hh|cu|cuh)$")
            list(APPEND format_sources ${file_abs})
        endif()
    endforeach()
    # Only add target if we have sources to format
    if(format_sources)
        # Create unique target name based on directory
        string(MD5 target_hash "${format_sources}")
        set(format_target_name "mars_format_${target_hash}")
        # Create the format target
        add_custom_target(
            ${format_target_name}
            COMMAND ${CLANG_FORMAT_EXECUTABLE}
                    -style=file
                    -fallback-style=none
                    -i
                    ${format_sources}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMENT "Formatting ${format_sources}"
            VERBATIM
        )
        # Add as dependency to the main format target
        add_dependencies(format ${format_target_name})
    endif()
endfunction()

# Set up formatting for target sources
function(mars_format_add_target)
    if(NOT CLANG_FORMAT_EXECUTABLE)
        return()
    endif()

    # Process each target and collect its sources
    set(all_sources "")
    foreach(target_name ${ARGN})
        if(TARGET ${target_name})
            # Get target sources
            get_target_property(target_sources ${target_name} SOURCES)
            if(target_sources)
                # Get target source directory to resolve relative paths
                get_target_property(target_dir ${target_name} SOURCE_DIR)
                # Add each source with full path
                foreach(source ${target_sources})
                    if(NOT IS_ABSOLUTE ${source})
                        set(source "${target_dir}/${source}")
                    endif()
                    list(APPEND all_sources ${source})
                endforeach()
            endif()
        else()
            message(WARNING "Target ${target_name} does not exist, skipping formatting")
        endif()
    endforeach()
    # Format the collected sources
    if(all_sources)
        mars_format_add_sources(${all_sources})
    endif()
endfunction()

# Add compatibility functions for backward compatibility with older code
function(target_format_setup)
    mars_format_add_target(${ARGN})
endfunction()

function(format_setup)
    mars_format_add_sources(${ARGN})
endfunction()