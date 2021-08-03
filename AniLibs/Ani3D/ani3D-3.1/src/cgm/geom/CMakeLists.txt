project(cubit_geom)

option(CGM_PARALLEL "Build CGM with parallel support" OFF)
option(CGM_OCC "Build CGM with Open Cascade support" ON)

set(geom_srcs
  AnalyticGeometryTool.cpp
  AutoMidsurfaceTool.cpp
  BasicTopologyEntity.cpp
  Body.cpp
  BodySM.cpp
  BoundingBoxTool.cpp
  BridgeManager.cpp
  CAActuateSet.cpp
  CADeferredAttrib.cpp
  CAEntityColor.cpp
  CAEntityId.cpp
  CAEntityName.cpp
  CAEntitySense.cpp
  CAEntityTol.cpp
  CAGroup.cpp
  CAMergePartner.cpp
  CAMergeStatus.cpp
  CASourceFeature.cpp
  CAUniqueId.cpp
  CGMApp.cpp
  CGMEngineDynamicLoader.cpp
  CGMHistory.cpp
  Chain.cpp
  CoEdge.cpp
  CoEdgeSM.cpp
  CoFace.cpp
  CoVertex.cpp
  CoVolume.cpp
  CollectionEntity.cpp
  CubitAttrib.cpp
  CubitAttribManager.cpp
  CubitAttribUser.cpp
  CubitCompat.cpp
  CubitPolygon.cpp
  CubitSimpleAttrib.cpp
  Curve.cpp
  CurveOverlapFacet.cpp
  CurveSM.cpp
  CylinderEvaluator.cpp
  DAG.cpp
  DagDrawingTool.cpp
  GSaveOpen.cpp
  GeomDataObserver.cpp
  GeomMeasureTool.cpp
  GeometryEntity.cpp
  GeometryEvent.cpp
  GeometryFeatureEngine.cpp
  GeometryFeatureTool.cpp
  GeometryHealerEngine.cpp
  GeometryHealerTool.cpp
  GeometryModifyEngine.cpp
  GeometryModifyTool.cpp
  GeometryQueryEngine.cpp
  GeometryQueryTool.cpp
  GeometryUtil.cpp
  GfxPreview.cpp
  GroupingEntity.cpp
  LocalToleranceTool.cpp
  Loop.cpp
  LoopSM.cpp
  Lump.cpp
  LumpSM.cpp
  MedialTool2D.cpp
  MedialTool3D.cpp
  MergeTool.cpp
  MergeToolAssistant.cpp
  MidPlaneTool.cpp
  ModelQueryEngine.cpp
  OffsetSplitTool.cpp
  OldUnmergeCode.cpp
  Point.cpp
  PointSM.cpp
  RefCollection.cpp
  RefEdge.cpp
  RefEntity.cpp
  RefEntityFactory.cpp
  RefEntityName.cpp
  RefFace.cpp
  RefGroup.cpp
  RefVertex.cpp
  RefVolume.cpp
  SenseEntity.cpp
  Shell.cpp
  ShellSM.cpp
  SphereEvaluator.cpp
  SplitSurfaceTool.cpp
  SurfParamTool.cpp
  Surface.cpp
  SurfaceOverlapFacet.cpp
  SurfaceOverlapTool.cpp
  SurfaceSM.cpp
  TBOwner.cpp
  TBOwnerSet.cpp
  TDCAGE.cpp
  TDSourceFeature.cpp
  TDSplitSurface.cpp
  TDSurfaceOverlap.cpp
  TDUniqueId.cpp
  TopologyBridge.cpp
  TopologyEntity.cpp)

if (CAT_BUILD)
  list(APPEND geom_srcs
    CAProWeld.cpp
    TDProWeld.cpp)
endif ()

set(extra_geom_srcs
  AllocMemManagersGeom.cpp)

# Create a list of header-only files:
set(geom_headers
  AnalyticGeometryTool.hpp
  AutoMidsurfaceTool.hpp
  BasicTopologyEntity.hpp
  Body.hpp
  BodySM.hpp
  BoundingBoxTool.hpp
  BridgeManager.hpp
  CAActuateSet.hpp
  CADeferredAttrib.hpp
  CADefines.hpp
  CAEntityColor.hpp
  CAEntityId.hpp
  CAEntityName.hpp
  CAEntityTol.hpp
  CAEntitySense.hpp
  CAGroup.hpp
  CAMergePartner.hpp
  CAMergeStatus.hpp
  CASourceFeature.hpp
  CAUniqueId.hpp
  CGMApp.hpp
  CGMEngineDynamicLoader.hpp
  CGMHistory.hpp
  Chain.hpp
  CoEdge.hpp
  CoEdgeSM.hpp
  CoFace.hpp
  CoVertex.hpp
  CoVolume.hpp
  CollectionEntity.hpp
  CubitAttrib.hpp
  CubitAttribManager.hpp
  CubitAttribUser.hpp
  CubitCompat.hpp
  CubitCompat.h
  CubitEvaluator.hpp
  CubitGeomConfigure.h
  CubitPolygon.hpp
  CubitSimpleAttrib.hpp
  Curve.hpp
  CurveOverlapFacet.hpp
  CurveSM.hpp
  CylinderEvaluator.hpp
  DAG.hpp
  DagDrawingTool.hpp
  DagType.hpp
  GSaveOpen.hpp
  GeomDataObserver.hpp
  GeomMeasureTool.hpp
  GeomPoint.hpp
  GeomSeg.hpp
  GeometryEntity.hpp
  GeometryFeatureEngine.hpp
  GeometryFeatureTool.hpp
  GeometryHealerEngine.hpp
  GeometryHealerTool.hpp
  GeometryModifyEngine.hpp
  GeometryModifyTool.hpp
  GeometryQueryEngine.hpp
  GeometryQueryTool.hpp
  GeometryUtil.hpp
  GfxPreview.hpp
  GroupingEntity.hpp
  IntermediateGeomEngine.hpp
  LocalToleranceTool.hpp
  Loop.hpp
  LoopSM.hpp
  Lump.hpp
  LumpSM.hpp
  MedialTool2D.hpp
  MedialTool3D.hpp
  MergeTool.hpp
  MergeToolAssistant.hpp
  MidPlaneTool.hpp
  ModelQueryEngine.hpp
  OffsetSplitTool.hpp
  OldUnmergeCode.hpp
  Point.hpp
  PointSM.hpp
  RefCollection.hpp
  RefEdge.hpp
  RefEntity.hpp
  RefEntityFactory.hpp
  RefEntityName.hpp
  RefFace.hpp
  RefGroup.hpp
  RefVertex.hpp
  RefVolume.hpp
  SenseEntity.hpp
  Shell.hpp
  ShellSM.hpp
  SphereEvaluator.hpp
  SplitSurfaceTool.hpp
  SurfParamTool.hpp
  Surface.hpp
  SurfaceOverlapFacet.hpp
  SurfaceOverlapTool.hpp
  SurfaceSM.hpp
  TBOwner.hpp
  TBOwnerSet.hpp
  TDCAGE.hpp
  TDCompare.hpp
  TDSourceFeature.hpp
  TDSplitSurface.hpp
  TDSurfaceOverlap.hpp
  TDUniqueId.hpp
  TopologyBridge.hpp
  TopologyEntity.hpp
  # Generated files:
  "${CMAKE_CURRENT_BINARY_DIR}/CGMConfigure.h"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitGeomConfigure.h")

set(CUBIT_GEOM_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CubitGeomConfigure.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitGeomConfigure.h"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CGMConfigure.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CGMConfigure.h"
  @ONLY)

# Now algorithmically generate header names
# from implementation filenames:
foreach (src IN LISTS geom_srcs)
  string(REGEX REPLACE ".cpp" ".hpp" header "${src}")
  list(APPEND geom_headers
    "${header}")
endforeach ()

if (CUBIT_FOUND)
  # Make an interface target which causes users to link to CUBIT's cubit_geom
  # library instead.
  set(cubit_geom_target CUBIT::Geom PARENT_SCOPE)
  set(cubit_geom_target CUBIT::Geom)
  # INTERFACE targets cannot deal with this variable existing from a previous
  # non-CUBIT build.
  unset(cubit_geom_LIB_DEPENDS CACHE)
  add_library(cubit_geom INTERFACE)
  set(cubit_geom_include_type INTERFACE)
  target_link_libraries(cubit_geom
    INTERFACE
      cubit_util
      CUBIT::Geom)
  install(
    TARGETS   cubit_geom
    EXPORT    cgm
    COMPONENT Runtime)
  set_property(GLOBAL APPEND
    PROPERTY cgma_export_targets
      cubit_geom)
else ()
  set(cubit_geom_target cubit_geom PARENT_SCOPE)
  set(cubit_geom_target cubit_geom)
  cgm_add_library(cubit_geom 
    ${geom_srcs}
    ${extra_geom_srcs})
  set(cubit_geom_include_type PUBLIC)
  if(CUBIT_LIBRARY_PROPERTIES)
    set_target_properties(cubit_geom
      PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
  endif()

  target_link_libraries(cubit_geom
    PUBLIC
      cubit_util)
endif ()

cgm_install_headers(${geom_headers})

target_include_directories(cubit_geom
  ${cubit_geom_include_type}
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")

# Make a library with the compat functions.
cgm_add_library(cubit_geom_compat STATIC
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h"
  CubitCompat.hpp
  CubitCompat.cpp)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CubitCompat.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h"
  @ONLY)

target_link_libraries(cubit_geom_compat
  PUBLIC
    cubit_geom)

cgm_install_headers(
  "${CMAKE_CURRENT_BINARY_DIR}/CubitCompat.h"
  CubitCompat.hpp)

add_subdirectory(Cholla)
add_subdirectory(facetbool)
add_subdirectory(facet)
add_subdirectory(virtual)

if (CGM_OCC)
  add_subdirectory(OCC)
endif ()

if (CGM_PARALLEL)
  add_subdirectory(parallel)
endif ()
