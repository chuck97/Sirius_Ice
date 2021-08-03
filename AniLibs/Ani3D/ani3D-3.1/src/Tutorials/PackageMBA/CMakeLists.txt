# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniMBA library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniMBA_analytic${SUFFICS_EXE})
add_executable(${TEST1} "main_analytic.f")
target_link_libraries(${TEST1} mba3D-${VERSION} view3D-${VERSION} ${LAPACK_LIBS} m)
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniMBA_nodal${SUFFICS_EXE})
add_executable(${TEST2} "main_nodal.f")
target_link_libraries(${TEST2} view3D-${VERSION} mba3D-${VERSION} ${LAPACK_LIBS} m)
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test 
set(TEST3 aniMBA_fixshape${SUFFICS_EXE})
add_executable(${TEST3} "main_fixshape.f")
target_link_libraries(${TEST3} mba3D-${VERSION} view3D-${VERSION} ${LAPACK_LIBS} m) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add fourth test
set(TEST4 aniMBA_readAFTmesh${SUFFICS_EXE})
add_executable(${TEST4} "main_read_aft_mesh.f")
target_link_libraries(${TEST4} mba3D-${VERSION} view3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST4} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add fifth test
set(TEST5 aniMBA_refine${SUFFICS_EXE})
add_executable(${TEST5} "main_refine.f")
target_link_libraries(${TEST5} mba3D-${VERSION} view3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST5} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add sixth test
set(TEST6 aniMBA_split_unite${SUFFICS_EXE})
add_executable(${TEST6} "main_split_unite.f")
target_link_libraries(${TEST6} mba3D-${VERSION} view3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST6} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

