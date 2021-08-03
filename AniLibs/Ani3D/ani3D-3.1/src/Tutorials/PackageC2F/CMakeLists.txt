# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniC2F library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniC2F_analytic${SUFFICS_EXE})
add_executable(${TEST1} "main_analytic.c")
target_link_libraries(${TEST1} c2f3D-${VERSION} mba3D-${VERSION} lmr3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniC2F_nodal${SUFFICS_EXE})
add_executable(${TEST2} "main_nodal.c")
target_link_libraries(${TEST2} c2f3D-${VERSION} mba3D-${VERSION} lmr3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})






