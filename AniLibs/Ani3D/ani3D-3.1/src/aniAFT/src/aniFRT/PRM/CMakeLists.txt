#
#  Ani3D 3-X
#  Advance Front Technique - support of PRM
#

add_library(frtprm-${VERSION} STATIC
	    "error3.c"  "libfrtprm.c"  "memory3.c"  "refine3.c"
	    "region3.c"  "struct3.c"  "support.c"  "tree32.c"
	    "tree3.c"  "tria32.c"  "user3.c"
	    )

install(TARGETS frtprm-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




