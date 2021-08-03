#
#  Ani3D 3-X
#  Advance Front Technique - support of SCG
#

add_library(frtscg-${VERSION} STATIC
            "libfrtscg.cpp"  "libprim.c"  "tritri1.c"  "mesh.cpp"  "raytri.c")

install(TARGETS frtscg-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




