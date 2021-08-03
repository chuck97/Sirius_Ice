#include "iGeom.h"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "Body.hpp"
#include "CubitVector.hpp"
#include "ModelQueryEngine.hpp"
#include "GeometryQueryTool.hpp"
#include <iostream>
#define CHECK( STR ) if (err != iBase_SUCCESS) return print_error( STR, err, geom, __FILE__, __LINE__ )

#ifdef HAVE_ACIS
#  define ENGINE "ACIS"
#  define FORMAT "ACIS_SAT"
#  define FILE_NAME "brick_2.sat"
#elif defined (HAVE_OCC)
#  define ENGINE "OCC"
#  define FORMAT "OCC"
#  define FILE_NAME "brick_2.stp"
#  define FILE_NAME1 "ilc_13body.stp"
#  define FILE_NAME2 "ilc_1body.stp"
#  define FILE_NAME3 "ilc_problem_surf8.stp"
#else
#  error "Which engine to test?"
#endif

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S

int getFirstVolume(iGeom_Instance geom, iBase_EntityHandle & volume);
int checkEdgeOrientation(iGeom_Instance geom, iBase_EntityHandle & volume);

static int print_error( const char* desc,
                         int err,
                         iGeom_Instance geom,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';
  
  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;
  
  return 1; // must always return false or CHECK macro will break
}

int main(int argc, char *argv[])
{
  // initialize the Mesh
  int i, err;
  iGeom_Instance geom;
  std::string engine_opt = ";engine=";
  engine_opt += ENGINE;
  iGeom_newGeom(engine_opt.c_str(), &geom, &err, engine_opt.length());

  // read in the geometry file, which is the geometry of a brick
  std::string input_file;
  input_file = STRINGIFY(SRCDIR)"/";
  input_file += FILE_NAME;
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry" );

  iBase_EntityHandle volume;
  err = getFirstVolume(geom, volume);
  if (err)
    return err;

  // get the bounding box of the brick and calculate its approximate center
  double minX, minY, minZ, maxX, maxY, maxZ;
  iGeom_getEntBoundBox(geom, volume, &minX, &minY, &minZ,
      &maxX, &maxY, &maxZ, &err);
  CHECK("Failed to get bounding box.");
  double cntrX = (minX + maxX) / 2.0;
  double cntrY = (minY + maxY) / 2.0;
  double cntrZ = (minZ + maxZ) / 2.0;

  // get brick faces
  iBase_EntityHandle* all_faces = NULL;
  int af_alloc = 0;
  int af_size = 0;
  iGeom_getEntAdj(geom, volume, iBase_FACE, &all_faces, &af_alloc,
                &af_size, &err);
  CHECK("Failed to get faces.");
  
#ifdef HAVE_ACIS
  // This is a check that is specific to ACIS.
  // Check that all face senses are FORWARD with respect to the parent brick
  for (i = 0; i < af_size; i++) { // for all faces
    // get face sense compared by volume
    int face_sense;
    iGeom_getEntNrmlSense(geom, all_faces[i], volume, &face_sense, &err);
    CHECK("Failed to get face sense.");
    
    // check if face sense is FORWARD respect to parent brick
    if (face_sense != 1) {
      std::cerr << "Error: face sense is not FORWARD." << std::endl;
      return 1;
    }
  }
#endif

  // This is a check that is always applicable for a brick.
  // Check that the face sense relative to the brick is FORWARD if the normal
  // points outward and REVERSED if the normal points inward
  for (i = 0; i < af_size; i++) { // for all faces
    // get face sense compared by volume
    int face_sense;
    iGeom_getEntNrmlSense(geom, all_faces[i], volume, &face_sense, &err);
    CHECK("Failed to get face sense.");
    double clsstX, clsstY, clsstZ, nrmlX, nrmlY, nrmlZ;
    iGeom_getEntNrmlPlXYZ(geom, all_faces[i], cntrX, cntrY, cntrZ,
        &clsstX, &clsstY, &clsstZ, &nrmlX, &nrmlY, &nrmlZ, &err);
    CHECK("Failed to get closest point and normal direction.");

    double dotPrdct = nrmlX * (cntrX - clsstX) + nrmlY * (cntrY - clsstY)
        + nrmlZ * (cntrZ - clsstZ);
    int normalDir = 1;
    if (dotPrdct > 0) // normal is same dir as vector from closest to center
      normalDir = -1;

    // check that face sense is forward if the normal points outward from
    // the brick and reversed if the normal points inward
    if (face_sense != normalDir) {
      std::cerr << "Error: face sense does not match direction of normal."
          << std::endl;
      return 1;
    }
  }

  err = checkEdgeOrientation(geom, volume);
  if (err)
    return err;

#if defined (HAVE_OCC)
  input_file = STRINGIFY(SRCDIR)"/";
  input_file  += FILE_NAME1;
  iGeom_deleteAll(geom, &err);
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry" );
  err = getFirstVolume(geom, volume);
  if (err)
    return err;
  err = checkEdgeOrientation(geom, volume);
  if (err)
    return err;


  input_file = STRINGIFY(SRCDIR)"/";
  input_file += FILE_NAME2;
  iGeom_deleteAll(geom, &err);
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry" );
  err = getFirstVolume(geom, volume);
  if (err)
    return err;
  err = checkEdgeOrientation(geom, volume);
  if (err)
    return err;

  input_file = STRINGIFY(SRCDIR)"/";
  input_file += FILE_NAME3;
  iGeom_deleteAll(geom, &err);
  iGeom_load(geom, input_file.c_str(), 0, &err, input_file.length(), 0);
  CHECK( "ERROR : can not load a geometry" );
  err = getFirstVolume(geom, volume);
  if (err)
    return err;
  err = checkEdgeOrientation(geom, volume);
  if (err)
    return err;
#endif

  std::cout << "All tests are passed." << std::endl;

  return 0;
}

int getFirstVolume(iGeom_Instance geom, iBase_EntityHandle & volume)
{
  int err;

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet(geom, &root_set, &err);
  CHECK("Failed to get root set.");
  
  // get the (first) volume
  iBase_EntityHandle* vols = NULL;
  int v_alloc = 0;
  int v_size = 0;
  iGeom_getEntities(geom, root_set, iBase_REGION, &vols,
                    &v_alloc, &v_size, &err);
  CHECK("Failed to get volumes.");

  volume = vols[0];

  return 0;
}

int checkEdgeOrientation(iGeom_Instance geom, iBase_EntityHandle & volume)
{

  int i, j, err;

  // get all edges
  iBase_EntityHandle* edges = NULL;
  int e_alloc = 0;
  int e_size = 0;
  iGeom_getEntAdj(geom, volume, iBase_EDGE, &edges, &e_alloc,
                  &e_size, &err);
  CHECK("Failed to get edges.");

  // Check that the edge senses to the two parent faces, each multiplied by
  // its parent face's sense to the volume, should be opposite
  for (i = 0; i < e_size; i++) { // for all edges
    // get parent faces 
    iBase_EntityHandle* faces = NULL;
    int f_alloc = 0;
    int f_size = 0;
    iGeom_getEntAdj(geom, edges[i], iBase_FACE, &faces, &f_alloc,
                    &f_size, &err);
    CHECK("Failed to get edges.");

    // check if # of parent faces of the edges on the volume is 2
    if (f_size != 2) {
      std::cerr << "Error: # of parent faces of solid edges should be 2."
          << std::endl;
      return 1;
    }

    // compute the edge to volume senses by multiplying the edge to face sense
    // by the parent face sense with respect to the volume.
    int face_sense;
    int sense[2];
    for (j = 0; j < 2; j++) {
      iGeom_getEgFcSense(geom, edges[i], faces[j], &sense[j], &err);
      CHECK("Failed to get edge to face sense.");
      iGeom_getEntNrmlSense(geom, faces[j], volume, &face_sense, &err);
      CHECK("Failed to get face to volume sense.");
      sense[j] *= face_sense;
    }
    
    // check if the edge to volume senses are opposite
    if (sense[0]*sense[1] != -1) {
      std::cerr << "Error: Edge senses to 2 parent faces should be opposite."
          << std::endl;
      return 1;
    }
  }

  std::cout << "Verified opposite edge to face senses for " << e_size
      << " edges." << std::endl;

  return 0;
}
