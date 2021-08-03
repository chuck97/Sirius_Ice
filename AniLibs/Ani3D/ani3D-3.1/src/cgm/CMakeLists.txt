project(CGMA)
cmake_minimum_required(VERSION 3.1)
SET(PACKAGE_NAME "CGM")
SET(PACKAGE_VERSION "15.1")

list(INSERT CMAKE_MODULE_PATH 0
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Choose static or shared libraries:
option(BUILD_SHARED_LIBS "Build CGMA with shared libraries." OFF)
option(USE_MPI "Use MPI" OFF)

option(CGM_USE_CUBIT "Use CUBIT" OFF)
if (CGM_USE_CUBIT)
  if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
  endif ()

  find_package(CUBIT REQUIRED)
endif ()

function (cgm_add_library name)
  add_library("${name}" ${ARGN})
  set_target_properties("${name}"
    PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${CGMA_BINARY_DIR}/lib"
      LIBRARY_OUTPUT_DIRECTORY "${CGMA_BINARY_DIR}/lib"
      RUNTIME_OUTPUT_DIRECTORY "${CGMA_BINARY_DIR}/bin")
  install(
    TARGETS   "${name}"
    RUNTIME   DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}
    LIBRARY   DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
    ARCHIVE   DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
    COMPONENT Runtime)
  set_property(GLOBAL APPEND
    PROPERTY cgma_libs
    "${name}")
  set_property(GLOBAL APPEND
    PROPERTY cgma_export_targets
    "${name}")
  get_target_property(target_type
    "${name}" TYPE)
  if (BUILD_SHARED_LIBS AND target_type STREQUAL "STATIC_LIBRARY")
    set_target_properties("${name}"
      PROPERTIES
        POSITION_INDEPENDENT_CODE TRUE)
  endif ()
endfunction ()

function (cgm_source_interface var)
  set(all)

  foreach (src IN LISTS ${ARGN})
    if (IS_ABSOLUTE "${src}")
      list(APPEND all
        "${src}")
    else ()
      list(APPEND all
        "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
    endif ()
  endforeach ()

  set(${var} "${all}" PARENT_SCOPE)
endfunction ()

function (cgm_install_headers)
  install(
    FILES       ${ARGN}
    DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}/CGM
    COMPONENT   Development)
endfunction ()

function (cgm_target_link_libraries name)
  target_link_libraries("${name}" ${ARGN})
  set_property(GLOBAL APPEND
    PROPERTY cgma_deplibs
    ${ARGN})
endfunction ()

add_subdirectory(util)
add_subdirectory(geom)

get_property(CGMA_LIBS GLOBAL
  PROPERTY cgma_libs)
get_property(CGMA_DEPLIBS GLOBAL
  PROPERTY cgma_deplibs)

add_subdirectory(init)
add_subdirectory(itaps)

set(export_sym "")
if (WIN32)
  set(export_sym "__declspec(dllexport) ")
endif ()
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp" "${export_sym}int cgm_dummy() { return 0; }")
add_library(cgm "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp")
target_link_libraries(cgm
  PUBLIC
    cgma_init
    cubit_geom
    cubit_geom_compat
    cubit_util)
install(TARGETS cgm
  RUNTIME   DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}
  LIBRARY   DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE   DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
  COMPONENT Runtime)

get_property(cgma_exports GLOBAL
  PROPERTY cgma_export_targets)

if (WIN32)
  set(cmakedir cmake)
else ()
  set(cmakedir lib/cmake/CGM)
endif ()

# Install top-level headers:
install(
  FILES
    ${CGMA_SOURCE_DIR}/cgm_version.h
  DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}/CGM
  COMPONENT Development
)

configure_file(
  "${CGMA_SOURCE_DIR}/CGMConfig.cmake.in"
  "${CGMA_BINARY_DIR}/CGMConfig.cmake"
  @ONLY)
install(
  FILES       "${CGMA_BINARY_DIR}/CGMConfig.cmake"
  DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
  COMPONENT   Development)
