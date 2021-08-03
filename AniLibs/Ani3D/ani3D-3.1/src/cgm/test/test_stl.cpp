#undef NDEBUG

#include <vector>

#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "CubitCompat.hpp"
#include "FacetQueryEngine.hpp"

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S
#ifndef SRCDIR
#  define SRCDIR "."
#endif

int main( int argc, char** argv ) {

  std::vector<CubitString> args;
  CGMApp::instance()->startup(args);

  const char* filename = STRINGIFY(SRCDIR) "/unit_cube.stl";

  CubitStatus status = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != status) return 1;

  DLIList<CubitQuadFacet*> quad_facet_list;
  DLIList<CubitFacet*> tri_facet_list;
  DLIList<Surface*> surface_list;
  FacetQueryEngine* fqe = FacetQueryEngine::instance();
  status = fqe->import_facets(filename, CUBIT_TRUE, 135, 1.e-6, 4, CUBIT_TRUE,
    CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE, quad_facet_list, tri_facet_list,
    surface_list, STL_FILE);
  if (CUBIT_SUCCESS != status) return 1;
 
  return 0;
}
