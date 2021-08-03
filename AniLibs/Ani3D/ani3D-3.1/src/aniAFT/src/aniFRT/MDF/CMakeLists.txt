#
#  Ani3D 3-X
#  Advance Front Technique - support of MDF
#

add_library(frtmdf-${VERSION} STATIC
            "support.c"  "surfrf.cpp"  "3Dto2D1.cpp"  
            "polyreg.cpp"  "libfrtmdf.c")

install(TARGETS frtmdf-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




