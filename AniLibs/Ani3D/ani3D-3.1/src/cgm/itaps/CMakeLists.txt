option(ITAPS_USE_BABEL "Should ITAPS use Babel to generate C and FORTRAN clients?" "OFF")

if (ITAPS_USE_BABEL) # TODO: cmakeify
  add_subdirectory(SIDL)
endif ()

set(cgma_igeom_srcs
  CATag.cpp
  iGeomError.cc
  iGeom_CGMA.cc)

set(cgma_igeom_hdrs
  iBase.h
  iBase_f.h
  iGeom.h
  iGeom_FCDefs.h
  iGeom_f.h
  iGeom_protos.h
  # Generated files:
  "${CMAKE_CURRENT_BINARY_DIR}/CgmaIGeomConfigure.h")

set(CGMA_IGEOM_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CgmaIGeomConfigure.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CgmaIGeomConfigure.h"
  @ONLY)

cgm_add_library(iGeom
  ${cgma_igeom_srcs}
  ${cgma_igeom_hdrs})
target_include_directories(iGeom
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(iGeom
  PUBLIC
    cubit_geom
    cubit_util
    ${CGMA_LIBS}
    cgma_init)
cgm_install_headers(${cgma_igeom_hdrs})
