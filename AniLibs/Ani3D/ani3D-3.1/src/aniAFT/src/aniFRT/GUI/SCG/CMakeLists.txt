# Find GLUT and OpenGL libraries and headers
cmake_minimum_required(VERSION 3.1)
find_package(GLUT)
find_package(OpenGL)

if(GLUT_FOUND)
# Make sure the compiler can find include files
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})

# Make sure the linker can find the aniAFT library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 gui_glut${SUFFICS_EXE})
add_executable(${TEST1}  "gui_glut.cpp" "gl2ps.c")
target_link_libraries(${TEST1} frtscg-${VERSION} frtprm-${VERSION} aft3D-${VERSION} GLUT::GLUT m GL) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})
endif(GLUT_FOUND)
