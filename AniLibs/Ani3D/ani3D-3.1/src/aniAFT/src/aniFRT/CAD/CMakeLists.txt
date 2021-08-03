#
#  Ani3D 3-X
#  Advance Front Technique - support of CAD 
#

include_directories(${cgma_init_SOURCE_DIR})
include_directories(${cgma_init_BINARY_DIR})
include_directories(${cubit_util_SOURCE_DIR})
include_directories(${cubit_util_BINARY_DIR})
include_directories(${cubit_geom_SOURCE_DIR})
include_directories(${cubit_geom_BINARY_DIR})
include_directories(${cubit_occ_SOURCE_DIR})
include_directories(${OPENCASCADE_INC_DIR})

add_library(frtcad-${VERSION} STATIC
            "libfrtcad.c"  "region3c.c"  "supportc.c"  "user3c.c"  "cgmwrap.cpp")
target_link_libraries(frtcad-${VERSION} cgm)
target_compile_definitions(frtcad-${VERSION} PUBLIC ${CGMA_OCC_DEFINES})

install(TARGETS frtcad-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




