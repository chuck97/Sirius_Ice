# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the libraries once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 multi_StokesNavier${SUFFICS_EXE})
add_executable(${TEST1} "main.f" "forlibfem.f" "forlibinb.f")
target_link_libraries(${TEST1} lmr3D-${VERSION} mba3D-${VERSION} 
                              view3D-${VERSION} inb-${VERSION} 
                               fem3D-${VERSION} lu-${UMFPACK_VER} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})




