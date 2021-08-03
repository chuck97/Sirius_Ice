project(cubit_virtual)

set(virtual_srcs
  CACompositeVG.cpp
  CAPartitionVG.cpp
  CAVirtualVG.cpp
  CollapseAngleTool.cpp
  CompSurfFacets.cpp
  CompositeAttrib.cpp
  CompositeBody.cpp
  CompositeCoEdge.cpp
  CompositeCoSurf.cpp
  CompositeCurve.cpp
  CompositeEngine.cpp
  CompositeGeom.cpp
  CompositeLoop.cpp
  CompositeLump.cpp
  CompositePoint.cpp
  CompositeShell.cpp
  CompositeSurface.cpp
  CompositeTool.cpp
  Faceter.cpp
  FaceterFacetData.cpp
  FaceterPointData.cpp
  HiddenEntitySet.cpp
  ImprintBoundaryTool.cpp
  ImprintLineSegment.cpp
  ImprintPointData.cpp
  PartPTCurve.cpp
  PartSurfFacetTool.cpp
  PartitionBody.cpp
  PartitionCoEdge.cpp
  PartitionCoSurf.cpp
  PartitionCurve.cpp
  PartitionEngine.cpp
  PartitionEntity.cpp
  PartitionLoop.cpp
  PartitionLump.cpp
  PartitionLumpImprint.cpp
  PartitionPoint.cpp
  PartitionShell.cpp
  PartitionSurface.cpp
  PartitionTool.cpp
  SegmentedCurve.cpp
  SubCurve.cpp
  SubEntitySet.cpp
  SubSurface.cpp
  TDVGFacetOwner.cpp
  TDVGFacetSplit.cpp
  VirtualImprintTool.cpp
  VirtualQueryEngine.cpp)

if (ALPHA_SPLIT_VIRTUAL)
  list(APPEND virtual_srcs
    SplitSurfaceVirtual.cpp)
endif ()

set(extra_virtual_srcs
  AllocMemManagersVirtual.cpp)

set(virtual_headers
  CACompositeVG.hpp
  CAPartitionVG.hpp
  CAVirtualVG.hpp
  CollapseAngleTool.hpp
  CompSurfFacets.hpp
  CompositeAttrib.hpp
  CompositeBody.hpp
  CompositeCoEdge.hpp
  CompositeCoSurf.hpp
  CompositeCurve.hpp
  CompositeEngine.hpp
  CompositeGeom.hpp
  CompositeLoop.hpp
  CompositeLump.hpp
  CompositePoint.hpp
  CompositeShell.hpp
  CompositeSurface.hpp
  CompositeTool.hpp
  Faceter.hpp
  FaceterFacetData.hpp
  FaceterPointData.hpp
  HiddenEntitySet.hpp
  ImprintBoundaryTool.hpp
  ImprintLineSegment.hpp
  ImprintMatchData.hpp
  ImprintPointData.hpp
  PartPTCurve.hpp
  PartSurfFacetTool.hpp
  PartitionBody.hpp
  PartitionCoEdge.hpp
  PartitionCoSurf.hpp
  PartitionCurve.hpp
  PartitionEngine.hpp
  PartitionEntity.hpp
  PartitionLoop.hpp
  PartitionLump.hpp
  PartitionLumpImprint.hpp
  PartitionPoint.hpp
  PartitionShell.hpp
  PartitionSurface.hpp
  PartitionTool.hpp
  SegmentedCurve.hpp
  SplitSurfaceVirtual.hpp
  SubCurve.hpp
  SubEntitySet.hpp
  SubSurface.hpp
  TDVGFacetOwner.hpp
  TDVGFacetSplit.hpp
  VGArray.cpp
  VGArray.hpp
  VGDefines.h
  VGLoopTool.cpp
  VGLoopTool.hpp
  VirtualImprintTool.hpp
  VirtualQueryEngine.hpp)

cgm_source_interface(all_virtual_srcs
  virtual_srcs
  extra_virtual_srcs)
  #virtual_headers) # FIXME: template recompilation

add_library(cubit_virtual INTERFACE)
target_sources(cubit_virtual
  INTERFACE
    ${all_virtual_srcs})
target_include_directories(cubit_virtual
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(cubit_virtual
  INTERFACE
    cubit_facet
    cubit_geom
    cubit_util
    cholla)
cgm_install_headers(${virtual_headers})

set_property(GLOBAL APPEND
  PROPERTY cgm_libs
  cubit_virtual)
