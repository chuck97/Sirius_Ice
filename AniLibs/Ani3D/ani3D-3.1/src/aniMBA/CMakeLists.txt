#
#  Ani3D 3-X
#  Metric Based Adaptation
#

set(CORE "mba_analytic.f"  "mba_nodal.f"  "mba_fixshape.f" 
         "ani2.f"  "forlibmba.f"  
         "makM.f"  "makQ.f"  "tangled.f"  "check_mesh.f" 
         "list.f"  "lintrp3D.f"  "dsort.f" 
         "nlnfnc.f"  "minim.f"  "missing.f"
         "E2F.f"  "F2E.f"  "swapR.f"  "insrtP.f"  "collapse_edge.f"  
         "move_point.f"  "delete_point.f"  "split_tet.f" 
         "auxSP.f"  "auxSR.f"  "auxSF.f"  "auxSE.f"  "checks.f" 
         "tetrahedron.f"  "triangle.f"
         "update.f"  "utils.f"  "utils_add.f"  "output.f"
         "refine.f"  "smoothing.f"  
         "ZZ.f"  "datacopy.f"
         "error.f"  "time.c"  "time_new.f"  "isnan.c"  "random.c")

set(IO  "loadM.f"  "saveM.f")


add_library(mba3D-${VERSION} STATIC ${CORE} ${IO} ${DEBUG})

install(TARGETS mba3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "status.fd" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





