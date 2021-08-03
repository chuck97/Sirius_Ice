#
#  Ani3D 3-X
#  Finite Element Method
#

add_library(fem3D-${VERSION} STATIC
            "fem3Dtet.f" "fem3Dface.f" "fem3Dtri.f" "fem3Derr.f" "fem3Dsub.f"
            "opIDEN.f" "opGRAD.f" "opDIV.f" "opCURL.f" "opDUDX.f" 
            "template.f" "bc.f" "algebra.f" "assemble.f" 
            "forlibfem.f" "isys.f" "dof.f" "error.f")

target_link_libraries(fem3D-${VERSION} ${LAPACK_LIBS})

install(TARGETS fem3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "assemble.fd"  "fem3Dtet.fd"  "fem3Dtet.h"  "fem3Dtri.fd" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





