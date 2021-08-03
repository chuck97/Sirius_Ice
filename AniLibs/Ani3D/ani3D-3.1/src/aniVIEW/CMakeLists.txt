#
#  Ani3D 3-X
#  VIzualization 
#

add_library(view3D-${VERSION} STATIC "drawGMV.f")

install(TARGETS view3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




