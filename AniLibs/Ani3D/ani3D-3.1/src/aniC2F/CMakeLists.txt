#
#  Ani3D 3-X
#  Interface to C-code
#

add_library(c2f3D-${VERSION} STATIC "ani3D.c")

install(TARGETS c2f3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
install(FILES "ani3D.h" DESTINATION ${INCLUDE_OUTPUT_DIRECTORY})





