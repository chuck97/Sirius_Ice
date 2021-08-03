#
#  Ani3D 3-X
#  Refinening/Coarsening Bisection
#

add_library(rcb3D-${VERSION} STATIC
           "aux.f"  "refine.f"  "coarse.f")

install(TARGETS rcb3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




