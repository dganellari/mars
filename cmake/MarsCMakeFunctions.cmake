# MarsCMakeFunctions.cmake

# ##############################################################################

macro(mars_add_library libraryRootDir subDirs)
  set(LOCAL_HEADERS "")
  set(LOCAL_SOURCES "")
  find_project_files(${libraryRootDir} "${subDirs}" LOCAL_HEADERS LOCAL_SOURCES)

  target_sources(mars PRIVATE ${LOCAL_SOURCES})

  install(FILES ${LOCAL_HEADERS} DESTINATION include)
  foreach(MODULE ${subDirs})
    target_include_directories(
      mars BEFORE PUBLIC $<BUILD_INTERFACE:${libraryRootDir}/${MODULE}>)
  endforeach(MODULE)
endmacro()

# ##############################################################################

function(find_project_files rootPath dirPaths headers sources)
  set(verbose TRUE)

  set(theaders ${${headers}})
  set(tsources ${${sources}})

  foreach(INCLUDE_PATH ${dirPaths})

    file(GLOB TEMP_HPPSRC CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.cpp")
    file(GLOB TEMP_SRC CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.c")
    file(GLOB TEMP_HPPDR CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.hpp")
    file(GLOB TEMP_HDR CONFIGURE_DEPENDS "${rootPath}/${INCLUDE_PATH}/*.h")

    source_group(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR};
                                       ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})

    set(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
    set(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
  endforeach(INCLUDE_PATH)

  set(${headers}
      ${theaders}
      PARENT_SCOPE)
  set(${sources}
      ${tsources}
      PARENT_SCOPE)
endfunction()

# ##############################################################################

function(scan_directories in_root_dir in_dirs_to_be_scanned out_includes
         out_headers out_sources)

  # APPEND the results
  set(local_includes ${${out_includes}})
  set(local_headers ${${out_headers}})
  set(local_sources ${${out_sources}})

  find_project_files(${in_root_dir} "${in_dirs_to_be_scanned}" local_headers
                     local_sources)

  foreach(dir ${in_dirs_to_be_scanned})
    list(APPEND local_includes ${in_root_dir}/${dir})
  endforeach(dir)

  set(${out_includes}
      ${local_includes}
      PARENT_SCOPE)
  set(${out_headers}
      ${local_headers}
      PARENT_SCOPE)
  set(${out_sources}
      ${local_sources}
      PARENT_SCOPE)

endfunction()
