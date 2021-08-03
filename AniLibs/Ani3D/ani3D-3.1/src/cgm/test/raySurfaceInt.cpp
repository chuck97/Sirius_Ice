/**
 * \file raySurfaceInt.cpp
 *
 * \brief a test of projecting a point onto a surface along a direction
 */

#include <iostream>

#include "Body.hpp"
#include "CubitCompat.hpp"
#include "DLIList.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "InitCGMA.hpp"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "Surface.hpp"

#ifdef HAVE_OCC
#define FILE_EXT "stp" 
#define GEO_TYPE "STEP"
#else
#define FILE_EXT  "sat"
#define GEO_TYPE "ACIS_SAT"
#endif

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S

bool testDirectedProjection();
bool testPieceOfTorusFromFile();

int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != status)
  {
    return 1;
  }

  if (!testDirectedProjection())
  {
    return 2;
  }

  std::cout << std::endl;
  if (!testPieceOfTorusFromFile())
  {
    return 3;
  }

  return 0;
}

bool testDirectedProjection()
{
  GeometryQueryTool* gqt = GeometryQueryTool::instance();
  GeometryModifyTool* gmt = GeometryModifyTool::instance();

  // create a sphere
  Body* sphereBody = gmt->sphere(5.0);
  std::cout << "Created sphere of radius 5 centered at (0, 0, 0)."
      << std::endl;

  // extract the surface of the sphere
  DLIList<RefEntity*> sphereChildren;
  sphereBody->get_child_ref_entities(sphereChildren);
  RefVolume* sphereVol = dynamic_cast<RefVolume*>(sphereChildren.get());
  if (!sphereVol)
  {
    std::cout << "The child of the sphere body is not a RefVolume."
        << std::endl;
    return false;
  }
  sphereVol->get_child_ref_entities(sphereChildren);
  RefFace* sphereFace = dynamic_cast<RefFace*>(sphereChildren.get());
  if (!sphereFace)
  {
    std::cout << "The child of the sphere volume is not a RefFace."
        << std::endl;
    return false;
  }
  Surface* sphere = sphereFace->get_surface_ptr();

  // test firing ray directly down at the north pole of the sphere
  CubitVector raySource01(0, 0, 10);
  CubitVector rayDir01(0, 0, -2); // NOTE: not a unit vector
  CubitVector rayHit01;
  CubitStatus rayResult01 =
      sphere->closest_point_along_vector(raySource01, rayDir01, rayHit01);
  if (rayResult01 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (0, 0, 10) in direction (0, 0, -1) "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 0, 10) in direction (0, 0, -1) hits at ("
      << rayHit01.x() << ", " << rayHit01.y() << ", "
      << rayHit01.z() << ")." << std::endl;
  if ((fabs(rayHit01.x() - 0) + fabs(rayHit01.y() - 0) +
      fabs(rayHit01.z() - 5)) > 1e-12)
  {
    std::cout << "FAIL: The ray should hit at (0, 0, 5)." << std::endl;
    return false;
  }

  // test firing ray from above the north pole at an angle that should
  // hit the sphere where x = 2.5 and y = 0.
  CubitVector raySource02(0, 0, 10);
  CubitVector rayDir02(0.4034491106775367, 0, -0.9150020847481741);
  CubitVector rayHit02;
  CubitStatus rayResult02 =
      sphere->closest_point_along_vector(raySource02, rayDir02, rayHit02);
  if (rayResult02 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (0, 0, 10)\n  in direction "
        << "(0.4034491106775367, 0, -0.9150020847481741)\n  "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 0, 10)\n  in direction "
      << "(0.4034491106775367, 0, -0.9150020847481741)\n  hits at ("
      << rayHit02.x() << ", " << rayHit02.y() << ", "
      << rayHit02.z() << ")." << std::endl;
  if ((fabs(rayHit02.x() - 2.5) + fabs(rayHit02.y() - 0) +
      fabs(rayHit02.z() - 4.330127018922193)) > 1e-12)
  {
    std::cout << "FAIL: The ray should hit at approximately "
        << "(2.5, 0, 4.330127018922193)." << std::endl;
    return false;
  }

  // test firing ray directly down from the center of the sphere
  CubitVector raySource03(0, 0, 0);
  CubitVector rayDir03(0, 0, -1);
  CubitVector rayHit03;
  CubitStatus rayResult03 =
      sphere->closest_point_along_vector(raySource03, rayDir03, rayHit03);
  if (rayResult03 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (0, 0, 0) in direction (0, 0, -1) "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 0, 0) in direction (0, 0, -1) hits at ("
      << rayHit03.x() << ", " << rayHit03.y() << ", "
      << rayHit03.z() << ")." << std::endl;
  if ((fabs(rayHit03.x() - 0) + fabs(rayHit03.y() - 0) +
      fabs(rayHit03.z() + 5)) > 1e-12)
  {
    std::cout << "FAIL: The ray should hit at (0, 0, -5)." << std::endl;
    return false;
  }
  std::cout << std::endl;

  // delete the sphere body
  gqt->delete_Body(sphereBody);


  // create a torus
  Body* torusBody = gmt->torus(3.0, 1.0);
  std::cout << "Created torus of major radius 3 and minor radius 1 "
      << "centered at (0, 0, 0)\n  with normal along the z-axis."
      << std::endl;

  // extract the surface of the torus
  DLIList<RefEntity*> torusChildren;
  torusBody->get_child_ref_entities(torusChildren);
  RefVolume* torusVol = dynamic_cast<RefVolume*>(torusChildren.get());
  if (!torusVol)
  {
    std::cout << "The child of the torus body is not a RefVolume."
        << std::endl;
    return false;
  }
  torusVol->get_child_ref_entities(torusChildren);
  RefFace* torusFace = dynamic_cast<RefFace*>(torusChildren.get());
  if (!torusFace)
  {
    std::cout << "The child of the torus volume is not a RefFace."
        << std::endl;
    return false;
  }
  Surface* torus = torusFace->get_surface_ptr();

  // test firing ray directly down from (0, 3.5, 1)
  CubitVector raySource04(0, 3.5, 1);
  CubitVector rayDir04(0, 0, -1);
  CubitVector rayHit04;
  CubitStatus rayResult04 =
      torus->closest_point_along_vector(raySource04, rayDir04, rayHit04);
  if (rayResult04 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (0, 3.5, 1) in direction (0, 0, -1) "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 3.5, 1) in direction (0, 0, -1)\n  hits at ("
      << rayHit04.x() << ", " << rayHit04.y() << ", "
      << rayHit04.z() << ")." << std::endl;
  if ((fabs(rayHit04.x() - 0) + fabs(rayHit04.y() - 3.5) +
      fabs(rayHit04.z() - 0.8660254037844386)) > 1e-12)
  {
    std::cout << "FAIL: The ray should hit at (0, 3.5, 0.8660254037844386)."
        << std::endl;
    return false;
  }

  // test firing ray directly down from (0, 4.5, 1)
  CubitVector raySource05(0, 4.5, 1);
  CubitVector rayDir05(0, 0, -1);
  CubitVector rayHit05;
  CubitStatus rayResult05 =
      torus->closest_point_along_vector(raySource05, rayDir05, rayHit05);
  if (rayResult05 != CUBIT_FAILURE)
  {
    std::cout << "FAIL: Fire from (0, 4.5, 1) in direction (0, 0, -1) "
        << "did not produce CUBIT_FAILURE." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 4.5, 1) in direction (0, 0, -1) "
      << "did not intersect torus,\n  as expected." << std::endl;

  // test firing ray towards (-1, 1, 0) from (0, 0, 0)
  CubitVector raySource06(0, 0, 0);
  CubitVector rayDir06(-1, 1, 0);
  CubitVector rayHit06;
  CubitStatus rayResult06 =
      torus->closest_point_along_vector(raySource06, rayDir06, rayHit06);
  if (rayResult06 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (0, 0, 0) in direction (-1, 1, 0) "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (0, 0, 0) in direction (-1, 1, 0)\n  hits at ("
      << rayHit06.x() << ", " << rayHit06.y() << ", "
      << rayHit06.z() << ")." << std::endl;
  if ((fabs(rayHit06.x() + 1.414213562373095) +
      fabs(rayHit06.y() - 1.414213562373095) + fabs(rayHit06.z() - 0)) > 1e-12)
  {
    std::cout << "FAIL: The ray should hit\n  "
        << "at (-1.414213562373095, 1.414213562373095, 0)."
        << std::endl;
    return false;
  }

  // delete the torus body
  gqt->delete_Body(torusBody);

  return true;
}

bool testPieceOfTorusFromFile()
{
  GeometryQueryTool* gqt = GeometryQueryTool::instance();

  std::string fileName(STRINGIFY(SRCDIR) "/pieceOfTorus01.");
  fileName = fileName + FILE_EXT;
  CubitStatus readGeoResult =
      CubitCompat_import_solid_model(fileName.c_str(), GEO_TYPE);
  if (readGeoResult != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Attempt to read geometry from " << fileName
        << " did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  
  // extract the torus sheet body
  DLIList<Body*> bodies;
  gqt->bodies(bodies);
  Body* pieceOfTorusBody = bodies.pop();
  std::cout << "Extracted piece of torus body from file" << std::endl;

  // extract the surface of the torus
  DLIList<RefEntity*> torusPieceChildren;
  pieceOfTorusBody->get_child_ref_entities(torusPieceChildren);
  RefVolume* torusPieceVol =
      dynamic_cast<RefVolume*>(torusPieceChildren.get());
  if (!torusPieceVol)
  {
    std::cout << "The child of the piece of torus body is not a RefVolume."
        << std::endl;
    return false;
  }
  torusPieceVol->get_child_ref_entities(torusPieceChildren);
  RefFace* torusPieceFace = dynamic_cast<RefFace*>(torusPieceChildren.get());
  if (!torusPieceFace)
  {
    std::cout << "The child of the piece of torus volume is not a RefFace."
        << std::endl;
    return false;
  }
  Surface* torusPiece = torusPieceFace->get_surface_ptr();

  // test firing ray in the direction (0.707107, 0, -0.707107) from
  // (3.17919, 0.1, 1.07919)
  CubitVector raySource01(3.17919, 0.1, 1.07919);
  CubitVector rayDir01(0.707107, 0, -0.707107);
  CubitVector rayHit01;
  CubitStatus rayResult01 =
      torusPiece->closest_point_along_vector(raySource01, rayDir01, rayHit01);
  if (rayResult01 != CUBIT_SUCCESS)
  {
    std::cout << "FAIL: Fire from (3.17919, 0.1, 1.07919) "
        << "in direction (0.707107, 0, -0.707107) "
        << "did not produce CUBIT_SUCCESS." << std::endl;
    return false;
  }
  std::cout << "Fire from (3.17919, 0.1, 1.07919) "
      << "in direction (0.707107, 0, -0.707107)\n  hits at ("
      << rayHit01.x() << ", " << rayHit01.y() << ", "
      << rayHit01.z() << ")." << std::endl;
  if ((fabs(rayHit01.x() - 3.30723) +
      fabs(rayHit01.y() - 0.1) + fabs(rayHit01.z() - 0.951148)) > 1e-4)
  {
    std::cout << "FAIL: The ray should hit\n  "
        << "at (3.30723, 0.1, 0.951148)."
        << std::endl;
    return false;
  }

  return true;
}
