project(cubit_facet)

set(facet_srcs
  CubitOctreeCell.cpp
  CubitOctree.cpp
  CubitOctreeGenerator.cpp
  CubitOctreeGeneratorVolumes.cpp
  CubitOctreeNode.cpp
  FacetAttrib.cpp
  FacetAttribSet.cpp
  FacetBody.cpp
  FacetCoEdge.cpp
  FacetCurve.cpp
  FacetLoop.cpp
  FacetLump.cpp
  FacetModifyEngine.cpp
  FacetParamTool.cpp
  FacetPoint.cpp
  FacetProjectTool.cpp
  FacetQueryEngine.cpp
  FacetShell.cpp
  FacetSurface.cpp
  FacetboolInterface.cpp
  GridSearchTree.cpp
  OctreeFacetPointData.cpp
  OctreeIntersectionData.cpp
  PST_Data.cpp
  TDOctreeRefEdge.cpp
  TDOctreeRefFace.cpp)

set(facet_headers
  CubitOctreeCell.hpp
  CubitOctreeConstants.hpp
  CubitOctreeGenerator.hpp
  CubitOctreeGeneratorVolumes.hpp
  CubitOctree.hpp
  CubitOctreeNode.hpp
  FacetAttrib.hpp
  FacetAttribSet.hpp
  FacetBody.hpp
  FacetCoEdge.hpp
  FacetCurve.hpp
  FacetGeometryCreator.hpp
  FacetLoop.hpp
  FacetLump.hpp
  FacetModifyEngine.hpp
  FacetParamTool.hpp
  FacetPoint.hpp
  FacetProjectTool.hpp
  FacetQueryEngine.hpp
  FacetShell.hpp
  FacetSurface.hpp
  FacetboolInterface.hpp
  GridSearchTree.hpp
  GridSearchTreeNode.hpp
  OctreeFacetPointData.hpp
  OctreeIntersectionData.hpp
  PST_Data.hpp
  TDOctreeRefEdge.hpp
  TDOctreeRefFace.hpp)

cgm_source_interface(all_facet_srcs
  facet_srcs
  facet_headers)

add_library(cubit_facet INTERFACE)
target_sources(cubit_facet
  INTERFACE
    ${all_facet_srcs})
target_include_directories(cubit_facet
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(cubit_facet
  INTERFACE
    cubit_facetboolstub
    cubit_geom
    cholla)
cgm_install_headers(${facet_headers})

set_property(GLOBAL APPEND
  PROPERTY cgm_libs
  cubit_facet)
