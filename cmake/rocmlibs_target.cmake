find_package(hip REQUIRED)

if (NOT TARGET mars::rocmlibs)
  add_library(mars::rocmlibs INTERFACE IMPORTED)
  target_link_libraries(mars::rocmlibs INTERFACE hip::device hip::host)
endif()
