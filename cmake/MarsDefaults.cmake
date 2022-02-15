if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR AND NOT MSVC_IDE)
    message(FATAL_ERROR "In-source builds are not allowed.")
endif()

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(LINUX)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(WIN32)
    add_definitions(/bigobj)
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG}   -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -MP -DWIN32_LEAN_AND_MEAN -DNOMINMAX")
endif()

# ##############################################################################
# CMake policies
if(CMAKE_VERSION VERSION_GREATER "3.13.0")
    cmake_policy(SET CMP0079 NEW)
endif()

set(CMAKE_MACOSX_RPATH 1)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ##############################################################################

# if(NOT MOONOLITH_LAUNCH_EXE) set(MOONOLITH_LAUNCH_EXE "") endif()