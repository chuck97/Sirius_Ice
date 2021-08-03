project(cubit_occ)

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(OpenCascade REQUIRED)

SET(occ_srcs
  OCCAttribSet.cpp
  OCCAttribSet.hpp
  OCCBody.cpp
  OCCBody.hpp
  OCCCoEdge.cpp
  OCCCoEdge.hpp
  OCCCurve.cpp
  OCCCurve.hpp
  OCCDrawTool.cpp
  OCCDrawTool.hpp
  OCCGeometryCreator.hpp
  OCCLoop.cpp
  OCCLoop.hpp
  OCCLump.cpp
  OCCLump.hpp
  OCCModifyEngine.cpp
  OCCModifyEngine.hpp
  OCCPoint.cpp
  OCCPoint.hpp
  OCCQueryEngine.cpp
  OCCQueryEngine.hpp
  OCCShapeAttributeSet.cpp
  OCCShapeAttributeSet.hpp
  OCCShell.cpp
  OCCShell.hpp
  OCCSurface.cpp
  OCCSurface.hpp)

set(occ_headers
  OCCAttribSet.hpp
  OCCBody.hpp
  OCCCoEdge.hpp
  OCCCurve.hpp
  OCCDrawTool.hpp
  OCCGeometryCreator.hpp
  OCCLoop.hpp
  OCCLump.hpp
  OCCModifyEngine.hpp
  OCCPoint.hpp
  OCCQueryEngine.hpp
  OCCShapeAttributeSet.hpp
  OCCShell.hpp
  OCCSurface.hpp)

set(OCC_DEFINES HAVE_OCC LIST)

if(OPENCASCADE_TKSTL_LIBRARY)
  list(APPEND OCC_DEFINES HAVE_OCC_STL)
endif()
if(OPENCASCADE_TKIGES_LIBRARY)
  list(APPEND OCC_DEFINES HAVE_OCC_IGES)
endif()
if(OPENCASCADE_TKSTEP_LIBRARY AND OPENCASCADE_TKSTEP209_LIBRARY AND OPENCASCADE_TKSTEPAttr_LIBRARY AND OPENCASCADE_TKSTEPBase_LIBRARY)
  list(APPEND OCC_DEFINES HAVE_OCC_STEP)
endif()

# Needed definitions for OCC versions < 6.8
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("iostream" HAVE_IOSTREAM)
CHECK_INCLUDE_FILE_CXX("iostream.h" HAVE_IOSTREAM_H)
CHECK_INCLUDE_FILE_CXX("iomanip" HAVE_IOMANIP)
CHECK_INCLUDE_FILE_CXX("iomanip.h" HAVE_IOMANIP_H)
CHECK_INCLUDE_FILE_CXX("fstream" HAVE_FSTREAM)
CHECK_INCLUDE_FILE_CXX("fstream.h" HAVE_FSTREAM_H)
CHECK_INCLUDE_FILE_CXX("limits.h" HAVE_LIMITS_H)
if (HAVE_IOSTREAM)
  list(APPEND OCC_DEFINES HAVE_IOSTREAM)
endif ()
if (HAVE_IOSTREAM_H)
  list(APPEND OCC_DEFINES HAVE_IOSTREAM_H)
endif ()
if (HAVE_IOMANIP)
  list(APPEND OCC_DEFINES HAVE_IOMANIP)
endif ()
if (HAVE_IOMANIP_H)
  list(APPEND OCC_DEFINES HAVE_IOMANIP_H)
endif ()
if (HAVE_FSTREAM)
  list(APPEND OCC_DEFINES HAVE_FSTREAM)
endif ()
if (HAVE_FSTREAM_H)
  list(APPEND OCC_DEFINES HAVE_FSTREAM_H)
endif ()
if (HAVE_LIMITS_H)
  list(APPEND OCC_DEFINES HAVE_LIMITS_H)
endif ()

if (WIN32)
  list(APPEND OCC_DEFINES WNT HAVE_OCC_DEF OCC_INC_FLAG)
endif ()

set(CGMA_OCC_DEFINES ${OCC_DEFINES} CACHE LIST "")
mark_as_advanced(CGMA_OCC_DEFINES)

# build cubit occ library from sources
cgm_source_interface(all_occ_srcs
  occ_srcs
  occ_headers)

add_library(cubit_occ INTERFACE)
target_sources(cubit_occ
  INTERFACE
    ${all_occ_srcs})
target_link_libraries(cubit_occ
  INTERFACE
    cubit_geom_compat
    cubit_geom
    cubit_util
    ${OPENCASCADE_LIBRARIES})
target_include_directories(cubit_occ
  INTERFACE
    ${OPENCASCADE_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:install>)
target_compile_definitions(cubit_occ
  INTERFACE
    ${CGMA_OCC_DEFINES})
cgm_install_headers(${occ_headers})

set_property(GLOBAL APPEND
  PROPERTY cgm_libs
  cubit_occ)
