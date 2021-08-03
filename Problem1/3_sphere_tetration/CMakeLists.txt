cmake_minimum_required(VERSION 3.17.0)

project(SPHERE_TETRATION)
enable_language (CXX C)
set(CMAKE_C_COMPILER gcc)

set(ANI_DIR "...")

set(ANI_LIB ${ANI_DIR}/lib)
set(ANI_INC ${ANI_DIR}/include)

include_directories(${ANI_INC})

list(APPEND SOURCES main.c)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCES})

target_link_libraries(${EXEC_NAME} ${ANI_LIB}/libaft3D-3.1.a )
target_link_libraries(${EXEC_NAME} ${ANI_LIB}/libfrtscg-3.1.a)
target_link_libraries(${EXEC_NAME} ${ANI_LIB}/libfrtprm-3.1.a)
target_link_libraries(${EXEC_NAME} -lm)
target_link_libraries(${EXEC_NAME} -lstdc++)

