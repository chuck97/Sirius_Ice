# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniAFT library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniAFT_front${SUFFICS_EXE})
add_executable(${TEST1}  "main_aft.c")
target_link_libraries(${TEST1} aft3D-${VERSION} m) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniAFT_param${SUFFICS_EXE})
add_executable(${TEST2}  "main_prm.c")
target_link_libraries(${TEST2} aft3D-${VERSION} frtprm-${VERSION} m) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniAFT_mdf${SUFFICS_EXE})
add_executable(${TEST3}  "main_mdf.c")
target_link_libraries(${TEST3} aft3D-${VERSION} frtmdf-${VERSION} m) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add fourth test
set(TEST4 aniAFT_scg${SUFFICS_EXE})
add_executable(${TEST4}  "main_scg.c")
target_link_libraries(${TEST4} aft3D-${VERSION} frtscg-${VERSION} frtprm-${VERSION} m) 
install(TARGETS ${TEST4} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add fifth test
if (ENABLE_OPENCASCADE)
  set(TEST5 aniAFT_cad${SUFFICS_EXE})
  add_executable(${TEST5}  "main_cad.c")
  target_link_libraries(${TEST5} aft3D-${VERSION} frtcad-${VERSION} frtprm-${VERSION} m) 
  install(TARGETS ${TEST5} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})
endif()


