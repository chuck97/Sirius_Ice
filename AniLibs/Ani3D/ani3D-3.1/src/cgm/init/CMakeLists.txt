project(cgma_init)

set(extra_libs)

if (CUBIT_FOUND)
  add_definitions(-DCUBIT_CGM)
  list(APPEND extra_libs
    CUBIT::CUBIT)
endif ()

include_directories(${OCE_INCLUDE_DIRS})

get_property(CGM_LIBS GLOBAL
  PROPERTY cgm_libs)

set(cgma_init_srcs
  InitCGMA.cpp)
set(cgma_init_hdrs
  InitCGMA.hpp
  # Generated files:
  "${cgma_init_BINARY_DIR}/CgmaInitConfigure.h")

cgm_add_library(cgma_init 
  ${cgma_init_srcs}
  ${cgma_init_hdrs})

set(CGMA_INIT_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/CgmaInitConfigure.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CgmaInitConfigure.h"
  @ONLY)
target_include_directories(cgma_init
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(cgma_init
  PUBLIC
    cubit_util
    cubit_geom
  PRIVATE
    ${CGM_LIBS}
    ${extra_libs})
cgm_install_headers(${cgma_init_hdrs})
