cmake_minimum_required(VERSION 3.17.0)

project(OUT_TO_TXT)
enable_language (CXX)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)

set(CMAKE_CXX_STANDARD 17)

list(APPEND SOURCE_FILES main.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )

