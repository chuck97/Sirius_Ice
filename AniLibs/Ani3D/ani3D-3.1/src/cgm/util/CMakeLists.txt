project(cubit_util)

# User option to build with concurrent execution support.
option(BUILD_WITH_CONCURRENT_SUPPORT "Build Cubit with concurrent execution support" OFF)

set(cgma_util_srcs
  AllocMemManagersList.cpp
  AppUtil.cpp
  ArrayBasedContainer.cpp
  CCubitFile.cpp
  CpuTimer.cpp
  Cubit2DPoint.cpp
  CubitBox.cpp
  CubitColor.cpp
  CubitCoordinateSystem.cpp
  CubitDirIterator.cpp
  CubitDynamicLoader.cpp
  CubitEntity.cpp
  CubitEvent.cpp
  CubitEventDispatcher.cpp
  CubitFile.cpp
  CubitFileFEModel.cpp
  CubitFileIOWrapper.cpp
  CubitFileMetaData.cpp
  CubitFileSimModel.cpp
  CubitFileUtil.cpp
  CubitInstrumentation.cpp
  CubitMatrix.cpp
  CubitMessage.cpp
  CubitMessageHandler.cpp
  CubitObservable.cpp
  CubitObserver.cpp
  CubitOperationEvent.cpp
  CubitPlane.cpp
  CubitProcess.cpp
  CubitSparseMatrix.cpp
  CubitString.cpp
  CubitTransformMatrix.cpp
  CubitUndo.cpp
  CubitUtil.cpp
  CubitVector.cpp
  DLList.cpp
  GMem.cpp
  GetLongOpt.cpp
  GfxDebug.cpp
  GlobalCommandFeedback.cpp
  IntersectionTool.cpp
  MemoryBlock.cpp
  MemoryManager.cpp
  ParamCubitPlane.cpp
  PlanarParamTool.cpp
  RandomMersenne.cpp
  SDLList.cpp
  SettingHandler.cpp
  SettingHolder.cpp
  StubProgressTool.cpp
  TDUPtr.cpp
  TextProgressTool.cpp
  ToolData.cpp
  ToolDataUser.cpp
  TtyProgressTool.cpp)

set(cgma_util_hdrs
  AbstractTree.hpp
  AppUtil.hpp
  ArrayBasedContainer.hpp
  CCubitFile.hpp
  CastTo.hpp
  CommandFeedback.hpp
  CpuTimer.hpp
  Cubit2DPoint.hpp
  CubitBox.hpp
  CubitBoxStruct.h
  CubitColorConstants.hpp
  CubitCoordEvent.hpp
  CubitCoordinateSystem.hpp
  CubitDefines.h
  CubitDirIterator.hpp
  CubitDynamicLoader.hpp
  CubitEntity.hpp
  CubitEvent.hpp
  CubitEventDispatcher.hpp
  CubitFile.hpp
  CubitFileFEModel.hpp
  CubitFileIOWrapper.hpp
  CubitFileMetaData.hpp
  CubitFileSimModel.hpp
  CubitFileUtil.hpp
  CubitInputFile.hpp
  CubitInstrumentation.hpp
  CubitLoops.hpp
  CubitMatrix.hpp
  CubitMessage.hpp
  CubitObservable.hpp
  CubitObserver.hpp
  CubitOperationEvent.hpp
  CubitPlane.hpp
  CubitPlaneStruct.h
  CubitProcess.hpp
  CubitSparseMatrix.hpp
  CubitString.hpp
  CubitTransformMatrix.hpp
  CubitUndo.hpp
  CubitUtil.hpp
  CubitVector.hpp
  CubitVectorStruct.h
  DLIList.hpp
  DLList.hpp
  DrawingToolDefines.h
  DynamicDLIIterator.hpp
  DynamicTreeIterator.hpp
  ElementType.h
  FacetShapeDefs.hpp
  GMem.hpp
  GeometryDefines.h
  GetLongOpt.hpp
  GfxDebug.hpp
  GlobalCommandFeedback.hpp
  IGUIObservers.hpp
  IndexedDouble.hpp
  IntersectionTool.hpp
  InvalidEntity.hpp
  KDDTree.cpp
  KDDTree.hpp
  KDDTreeNode.cpp
  KDDTreeNode.hpp
  LocalStart.h
  ManagedPtrVector.hpp
  MemoryBlock.hpp
  MemoryManager.hpp
  OctTree.cpp
  OctTree.hpp
  OctTreeCell.cpp
  OctTreeCell.hpp
  ParamCubitPlane.hpp
  ParamTool.hpp
  PlanarParamTool.hpp
  PriorityQueue.cpp
  PriorityQueue.hpp
  ProgressTool.hpp
  RStarTree.cpp
  RStarTree.hpp
  RStarTreeNode.cpp
  RStarTreeNode.hpp
  RTree.cpp
  RTree.hpp
  RTreeNode.cpp
  RTreeNode.hpp
  RandomMersenne.hpp
  SDLCAMergePartnerList.hpp
  SDLList.hpp
  SettingHandler.hpp
  SettingHolder.hpp
  StubProgressTool.hpp
  TDCellIndex.hpp
  TDUPtr.hpp
  TDVector.hpp
  TextProgressTool.hpp
  ToolData.hpp
  ToolDataUser.hpp
  TtyProgressTool.hpp
  CubitMessageHandler.hpp
  CubitColor.hpp
  # Generated files:
  "${CMAKE_CURRENT_BINARY_DIR}/CubitUtilConfigure.h")

set(CUBIT_UTIL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CubitUtilConfigure.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CubitUtilConfigure.h"
  @ONLY)

if (CUBIT_FOUND)
  # Make an interface target which causes users to link to CUBIT's cubit_util
  # library instead.
  set(cubit_util_target CUBIT::Util PARENT_SCOPE)
  # INTERFACE targets cannot deal with this variable existing from a previous
  # non-CUBIT build.
  unset(cubit_util_LIB_DEPENDS CACHE)
  add_library(cubit_util INTERFACE)
  set(cubit_util_include_type INTERFACE)
  target_link_libraries(cubit_util
    INTERFACE
      CUBIT::Util)
  target_compile_definitions(cubit_util
    INTERFACE
      -DHAVE_ACIS)
  install(
    TARGETS   cubit_util
    EXPORT    cgm
    COMPONENT Runtime)
  set_property(GLOBAL APPEND
    PROPERTY cgma_export_targets
      cubit_util)
else ()
  if (BUILD_WITH_CONCURRENT_SUPPORT)
    list(APPEND cgma_util_srcs
      CubitConcurrentApi.cpp
      CubitQtConcurrentApi.cpp)
    list(APPEND cgma_util_hdrs
      CubitConcurrentApi.h
      CubitQtConcurrentApi.h)
    find_package(Qt4 REQUIRED QtCore)
    include(${QT_USE_FILE})
  endif ()

  set(cubit_util_target cubit_util PARENT_SCOPE)
  cgm_add_library(cubit_util
    ${cgma_util_srcs})
    #${cgma_util_hdrs}) # FIXME: template recompilation
  set(cubit_util_include_type PUBLIC)
  target_link_libraries(cubit_util
    PRIVATE
      ${CMAKE_DL_LIBS})

  if (CUBIT_UTIL_NAME)
    set_target_properties(cubit_util
      PROPERTIES OUTPUT_NAME "${CUBIT_UTIL_NAME}")
  endif ()

  if (WIN32)
    cgm_target_link_libraries(cubit_util shlwapi)
  endif ()

  if (CUBIT_LIBRARY_PROPERTIES)
    set_target_properties(cubit_util
      PROPERTIES ${CUBIT_LIBRARY_PROPERTIES})
  endif ()

  if (BUILD_WITH_CONCURRENT_SUPPORT)
    target_link_libraries(cubit_util
      PRIVATE
        ${QT_LIBRARIES})
  endif ()
endif ()

target_include_directories(cubit_util
  ${cubit_util_include_type}
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
cgm_install_headers(${cgma_util_hdrs})
