# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniPRJ library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniPRJ${SUFFICS_EXE})
add_executable(${TEST1} "main.f")
target_link_libraries(${TEST1} view3D-${VERSION} ilu-${VERSION} prj3D-${VERSION}
                               fem3D-${VERSION} mba3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})






