cmake_minimum_required(VERSION 3.17.0)

project(GRID_TO_PMF)
enable_language (CXX C)
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_CXX_STANDARD 17)

set(INMOST_SRC "...")
set(INMOST_BUILD "...")

include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})

list(APPEND SOURCE_FILES main.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES})

target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
target_link_libraries(${EXEC_NAME} -ldl)
