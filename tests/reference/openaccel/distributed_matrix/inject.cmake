# Pass as -DCMAKE_PROJECT_mars_INCLUDE=<this file>. CMake includes it right after project(mars);
# the deferred include adds the gate targets once MARS's top-level CMakeLists has defined `mars`.
# Deferred calls may not create subdirectories (hence include), and expand their arguments at
# call time (hence EVAL, which bakes in this directory now). Remove with -UCMAKE_PROJECT_mars_INCLUDE.
cmake_language(EVAL CODE
    "cmake_language(DEFER DIRECTORY [[${CMAKE_SOURCE_DIR}]] CALL include [[${CMAKE_CURRENT_LIST_DIR}/gates.cmake]])")
