#
#  Ani3D 3-X
#  Projection from one mesh to another mesh
#

add_library(prj3D-${VERSION} STATIC
            "comp_geom.f" "metamesh_fast.f" "metamesh_slow.f" "rhs.f" "tet2tet.f")

install(TARGETS prj3D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




