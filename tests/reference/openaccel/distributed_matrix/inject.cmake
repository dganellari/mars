# Pass as -DCMAKE_PROJECT_mars_INCLUDE=<this file>. CMake includes it right after project(mars);
# the deferred includes add the distributed matrix and distributed SIMPLE gate targets once
# MARS's top-level CMakeLists has defined `mars`. Deferred calls may not create subdirectories
# (hence include) and expand their arguments at call time (hence EVAL, which bakes the paths
# in now). Remove with -UCMAKE_PROJECT_mars_INCLUDE.
foreach(_gates "${CMAKE_CURRENT_LIST_DIR}/gates.cmake" "${CMAKE_CURRENT_LIST_DIR}/../distributed_simple/gates.cmake")
    get_filename_component(_gates "${_gates}" ABSOLUTE)
    cmake_language(EVAL CODE "cmake_language(DEFER DIRECTORY [[${CMAKE_SOURCE_DIR}]] CALL include [[${_gates}]])")
endforeach()
