project(cubit_parallel)

if (USE_MPI)
  find_package(MPI REQUIRED)
endif ()

set(cubit_parallel_srcs
  CABodies.cpp
  CGMFileOptions.cpp
  CGMProcConfig.cpp
  TDParallel.cpp)

set(cubit_parallel_headers
  CABodies.hpp
  CGMFileOptions.hpp
  CGMParallelConventions.h
  CGMProcConfig.hpp
  CGMmpi.h
  CGMmpi_config.h
  TDParallel.hpp)

if (USE_MPI)
  list(APPEND cubit_parallel_srcs
    CGMParallelComm.cpp
    CGMReadParallel.cpp)
  list(APPEND cubit_parallel_headers
    CGMParallelComm.hpp
    CGMReadParallel.hpp)
endif ()

cgm_source_interface(all_cubit_parallel_srcs
  cubit_parallel_srcs
  cubit_parallel_headers)

add_library(cubit_parallel INTERFACE)
target_sources(cubit_parallel
  INTERFACE
    ${all_cubit_parallel_srcs})
target_include_directories(cubit_parallel
  INTERFACE
    "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(cubit_parallel
  INTERFACE
    cubit_geom
    cubit_util
    cubit_geom_compat)
if (USE_MPI)
  target_compile_definitions(cubit_parallel
    INTERFACE
      USE_MPI)
  target_include_directories(cubit_parallel
    INTERFACE
      ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(cubit_parallel
    INTERFACE
      ${MPI_CXX_LIBRARIES})
endif ()
cgm_install_headers(${cubit_parallel_headers})

set_property(GLOBAL APPEND
  PROPERTY cgm_libs
  cubit_parallel)
