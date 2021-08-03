# The name of the library to build
noinst_LTLIBRARIES = 
if BUILD_CGM
  noinst_LTLIBRARIES += libcubit_util.la
endif

# The directory where headers will be installed
libcubit_util_la_includedir = $(includedir)
AM_CPPFLAGS += -I$(top_builddir)/util

# The non-template sources
# If any of these do not need to be installed, move them
# to the _SOURCES list above.
libcubit_util_la_SOURCES = \
  AllocMemManagersList.cpp \
  AppUtil.cpp \
  ArrayBasedContainer.cpp \
  CCubitFile.cpp \
  CommandFeedback.hpp \
  CpuTimer.cpp \
  Cubit2DPoint.cpp \
  CubitBox.cpp \
  CubitColor.cpp \
  CubitCoordinateSystem.cpp \
  CubitDirIterator.cpp \
  CubitDynamicLoader.cpp \
  CubitEntity.cpp \
  CubitEvent.cpp \
  CubitEventDispatcher.cpp \
  CubitFile.cpp \
  CubitFileFEModel.cpp \
  CubitFileIOWrapper.cpp \
  CubitFileMetaData.cpp \
  CubitFileSimModel.cpp \
  CubitFileUtil.cpp \
  CubitMatrix.cpp \
  CubitMessage.cpp \
  CubitMessageHandler.cpp \
  CubitObservable.cpp \
  CubitObserver.cpp \
  CubitOperationEvent.cpp \
  CubitPlane.cpp \
  CubitSparseMatrix.cpp \
  CubitString.cpp \
  CubitTransformMatrix.cpp \
  CubitUndo.cpp \
  CubitUtil.cpp \
  CubitVector.cpp \
  DLList.cpp \
  GMem.cpp \
  GetLongOpt.cpp \
  GfxDebug.cpp \
  IntersectionTool.cpp \
  MemoryBlock.cpp \
  MemoryManager.cpp \
  ParamCubitPlane.cpp \
  PlanarParamTool.cpp \
  RandomMersenne.cpp \
  SDLList.cpp \
  SettingHandler.cpp \
  SettingHolder.cpp \
  StubProgressTool.cpp \
  TDUPtr.cpp \
  TextProgressTool.cpp \
  ToolData.cpp \
  ToolDataUser.cpp \
  TtyProgressTool.cpp 

#if WITH_CUBIT
  libcubit_util_la_SOURCES += GlobalCommandFeedback.cpp \
      CubitProcess.cpp \
      CubitInstrumentation.cpp
#endif 

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_util_la_include_HEADERS = \
  AbstractTree.hpp \
  AppUtil.hpp \
  ArrayBasedContainer.hpp \
  CCubitFile.hpp \
  CastTo.hpp \
  CpuTimer.hpp \
  Cubit2DPoint.hpp \
  CubitBox.hpp \
  CubitBoxStruct.h \
  CubitColorConstants.hpp \
  CubitCoordEvent.hpp \
  CubitCoordinateSystem.hpp \
  CubitDefines.h \
  CubitDirIterator.hpp \
  CubitDynamicLoader.hpp \
  CubitEntity.hpp \
  CubitEvent.hpp \
  CubitEventDispatcher.hpp \
  CubitFile.hpp \
  CubitFileFEModel.hpp \
  CubitFileIOWrapper.hpp \
  CubitFileMetaData.hpp \
  CubitFileSimModel.hpp \
  CubitLoops.hpp \
  CubitUndo.hpp \
  CubitFileUtil.hpp \
  CubitInputFile.hpp \
  CubitMatrix.hpp \
  CubitMessage.hpp \
  CubitObservable.hpp \
  CubitObserver.hpp \
  CubitOperationEvent.hpp \
  CubitPlane.hpp \
  CubitPlaneStruct.h \
  CubitSparseMatrix.hpp \
  CubitString.hpp \
  CubitTransformMatrix.hpp \
  CubitUtil.hpp \
  CubitUtilConfigure.h \
  CubitVector.hpp \
  CubitVectorStruct.h \
  DLIList.hpp \
  DLList.hpp \
  DrawingToolDefines.h \
  DynamicDLIIterator.hpp \
  DynamicTreeIterator.hpp \
  ElementType.h \
  FacetShapeDefs.hpp \
  GMem.hpp \
  GeometryDefines.h \
  GetLongOpt.hpp \
  GfxDebug.hpp \
  IGUIObservers.hpp \
  IndexedDouble.hpp \
  IntersectionTool.hpp \
  InvalidEntity.hpp \
  KDDTree.hpp \
  KDDTreeNode.hpp \
  LocalStart.h \
  ManagedPtrVector.hpp \
  MemoryBlock.hpp \
  MemoryManager.hpp \
  OctTree.hpp \
  OctTreeCell.hpp \
  ParamCubitPlane.hpp \
  ParamTool.hpp \
  PlanarParamTool.hpp \
  PriorityQueue.hpp \
  ProgressTool.hpp \
  RStarTree.hpp \
  RStarTreeNode.hpp \
  RTree.hpp \
  RTreeNode.hpp \
  RandomMersenne.hpp \
  SDLCAMergePartnerList.hpp \
  SDLList.hpp \
  SettingHandler.hpp \
  SettingHolder.hpp \
  StubProgressTool.hpp \
  TDCellIndex.hpp \
  TDUPtr.hpp \
  TDVector.hpp \
  TextProgressTool.hpp \
  ToolData.hpp \
  ToolDataUser.hpp \
  TtyProgressTool.hpp \
  CubitUtilConfigure.h \
  CubitMessageHandler.hpp \
  CubitColor.hpp

#if WITH_CUBIT
  libcubit_util_la_include_HEADERS += GlobalCommandFeedback.hpp \
           CubitProcess.hpp \
           CubitInstrumentation.hpp
#endif 
 
# If template defs are included, then the template definitions
# need to be installed with the headers.  Otherwise they need
# to be compiled.
if INCLUDE_TEMPLATE_DEFS
  libcubit_util_la_include_HEADERS += \
      KDDTree.cpp \
      KDDTreeNode.cpp \
      OctTreeCell.cpp \
      OctTree.cpp \
      PriorityQueue.cpp \
      RTree.cpp \
      RTreeNode.cpp \
      RStarTree.cpp \
      RStarTreeNode.cpp
else
  libcubit_util_la_SOURCES += \
      KDDTree.cpp \
      KDDTreeNode.cpp \
      OctTreeCell.cpp \
      OctTree.cpp \
      PriorityQueue.cpp \
      RTree.cpp \
      RTreeNode.cpp \
      RStarTree.cpp \
      RStarTreeNode.cpp
endif

libcubit_util_la_LIBADD = ${CGM_EXT_LDFLAGS} -ldl 
