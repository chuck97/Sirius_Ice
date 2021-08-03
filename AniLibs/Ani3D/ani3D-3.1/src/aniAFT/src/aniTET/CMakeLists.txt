#
#  Ani3D 3-X
#  Advance Front Technique
#

add_library(aft3D-${VERSION} STATIC
            "aft3d.c"  "aft2d.c"    "refine.c" 
            "check.c"  "delaunay.c"  "libaft.c"  "ugly.c"
            "cdt.c"  "qual.c"  "helper.c"  "det.c")

install(TARGETS aft3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




