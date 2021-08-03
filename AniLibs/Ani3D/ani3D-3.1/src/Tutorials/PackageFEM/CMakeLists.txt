# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniFEM library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniFEM_BC${SUFFICS_EXE})
add_executable(${TEST1} "main_bc.f")
target_link_libraries(${TEST1} fem3D-${VERSION} mba3D-${VERSION}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniFEM_tetrahedron${SUFFICS_EXE})
add_executable(${TEST2} "main_tetrahedron.f")
target_link_libraries(${TEST2} fem3D-${VERSION} mba3D-${VERSION}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniFEM_template${SUFFICS_EXE})
add_executable(${TEST3} "main_template.f")
target_link_libraries(${TEST3} fem3D-${VERSION} mba3D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})





