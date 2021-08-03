project(cubit_facetbool)

set(facetbool_srcs
  FBClassify.cpp
  FBDataUtil.cpp
  FBImprint.cpp
  FBIntersect.cpp
  FBPolyhedron.cpp
  FBRetriangulate.cpp
  FBTiler.cpp
  IntegerHash.cpp
  KdTree.cpp)

set(facetbool_headers
  FBDefines.hpp
  FBStructs.hpp)

foreach (src IN LISTS facetbool_srcs)
  string(REGEX REPLACE ".cpp" ".hpp" header "${src}")
  list(APPEND facetbool_headers
    "${header}")
endforeach ()

cgm_source_interface(all_facetbool_srcs
  facetbool_srcs
  facetbool_headers)

add_library(cubit_facetboolstub INTERFACE)
target_sources(cubit_facetboolstub
  INTERFACE
    ${all_facetbool_srcs})
target_include_directories(cubit_facetboolstub
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(cubit_facetboolstub
  INTERFACE
    cubit_util)
cgm_install_headers(${facetbool_headers})

set_property(GLOBAL APPEND
  PROPERTY cgm_libs
  cubit_facetboolstub)
