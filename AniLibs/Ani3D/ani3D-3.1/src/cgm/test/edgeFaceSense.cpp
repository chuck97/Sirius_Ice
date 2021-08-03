/**
 * \file edgeFaceSense.cpp
 *
 * \brief a test of whether the edge to face sense is correct for
 *    a cylinder created by the geometry engine
 *
 * This test uses the geometry engine to create a cylinder and checks
 * to see whether the coedge sense entities have the correct senses.  It
 * is a regression test, since this was a problem for cylinders created
 * with the OCC geometry engine.
 */

#include <cmath>
#include <iostream>

#include "BasicTopologyEntity.hpp"
#include "Body.hpp"
#include "CoEdge.hpp"
#include "CubitBox.hpp"
#include "GeometryModifyTool.hpp"
#include "GroupingEntity.hpp"
#include "InitCGMA.hpp"
#include "RefEdge.hpp"
#include "SenseEntity.hpp"

bool testSense();

int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != status)
  {
    return 1;
  }

  if (!testSense())
  {
    return 2;
  }

  return 0;
}

bool testSense()
{
  GeometryModifyTool *modToolPtr = GeometryModifyTool::instance();

  double height = 10.0;
  double radius = 2.5;
  double tol = 4e-5;
  Body *cylBodyPtr = modToolPtr->cylinder(height, radius, radius, radius);
  // assumes that the body contains one covolume corresponding to one volume
  // and that the volume has one shell
  GroupingEntity *cylShellPtr = cylBodyPtr->get_first_sense_entity_ptr()->
      get_basic_topology_entity_ptr()->get_first_grouping_entity_ptr();

  BasicTopologyEntity *vertSurfPtr = NULL;
  double topZVal = 0.0;
  double bottomZVal = 0.0;
  SenseEntity* cylCoFacePtr = cylShellPtr->get_first_sense_entity_ptr();
  while (cylCoFacePtr != NULL)
  {
    CubitBox boundingBox =
        cylCoFacePtr->get_basic_topology_entity_ptr()->bounding_box();
    if (fabs(boundingBox.z_range() - height) < tol)
    {
      topZVal = boundingBox.max_z();
      bottomZVal = boundingBox.min_z();
      vertSurfPtr = cylCoFacePtr->get_basic_topology_entity_ptr();
      break;
    }
  }

  if (vertSurfPtr == NULL)
  {
    std::cout << "Failed to find vertical face." << std::endl;
    return false;
  }

  CoEdge *topCoEdgePtr = NULL;
  CoEdge *bottomCoEdgePtr = NULL;
  GroupingEntity *vsLoopPtr = vertSurfPtr->get_first_grouping_entity_ptr();
  while (vsLoopPtr != NULL &&
      (topCoEdgePtr == NULL || bottomCoEdgePtr == NULL))
  {
    SenseEntity *coEdgePtr = vsLoopPtr->get_first_sense_entity_ptr();
    while (coEdgePtr != NULL &&
        (topCoEdgePtr == NULL || bottomCoEdgePtr == NULL))
    {
      CubitBox edgeBoundBox =
          coEdgePtr->get_basic_topology_entity_ptr()->bounding_box();
      if (edgeBoundBox.z_range() < tol)
      {
        if (fabs(edgeBoundBox.max_z() - topZVal) < tol)
        {
          topCoEdgePtr = dynamic_cast<CoEdge*>(coEdgePtr);
        }
        else if (fabs(edgeBoundBox.max_z() - bottomZVal) < tol)
        {
          bottomCoEdgePtr = dynamic_cast<CoEdge*>(coEdgePtr);
        }
      }
      coEdgePtr = coEdgePtr->next();
    }
    vsLoopPtr = vsLoopPtr->next();
  }

  if (topCoEdgePtr == NULL || bottomCoEdgePtr == NULL)
  {
    std::cout << "Failed to find top and/or bottom coedge." << std::endl;
    return false;
  }

  // determine the orientation of the top edge from parametrization
  bool topRHCC = false; // top right-handed counterclockwise
  RefEdge* topEdgePtr = topCoEdgePtr->get_ref_edge_ptr();
  double startParam = topEdgePtr->start_param();
  double endParam = topEdgePtr->end_param();
  double startDeltaParam = (endParam - startParam) / 64;
  CubitVector startVec, startDeltaVec;
  topEdgePtr->position_from_u(startParam, startVec);
  topEdgePtr->position_from_u(startDeltaParam, startDeltaVec);
  CubitVector normal = startVec * startDeltaVec; // cross product
  if (normal.z() > 0)
  {
    std::cout << "Top edge is parametrized/oriented counterclockwise "
        << std::endl
        << "    around top face looking down from above." << std::endl;
    topRHCC = true;
  }
  else
  {
    std::cout << "Top edge is parametrized/oriented clockwise "
        << std::endl
        << "    around top face looking down from above." << std::endl;
  }

  // determine the orientation of the bottom edge from parametrization
  bool bottomLHCC = false; // bottom left-handed counterclockwise
  RefEdge* bottomEdgePtr = bottomCoEdgePtr->get_ref_edge_ptr();
  startParam = bottomEdgePtr->start_param();
  endParam = bottomEdgePtr->end_param();
  startDeltaParam = (endParam - startParam) / 64;
  bottomEdgePtr->position_from_u(startParam, startVec);
  bottomEdgePtr->position_from_u(startDeltaParam, startDeltaVec);
  normal = startVec * startDeltaVec; // cross product
  if (normal.z() > 0)
  {
    std::cout << "Bottom edge is parametrized/oriented clockwise "
        << std::endl
        << "    around bottom face looking up from below." << std::endl;
  }
  else
  {
    std::cout << "Bottom edge is parametrized/oriented counterclockwise "
        << std::endl
        << "    around bottom face looking up from below." << std::endl;
    bottomLHCC = true;
  }

  // report what is expected and what the result is
  bool expectSameSense = (topRHCC == bottomLHCC);
  std::cout << "The coedges of the vertical cylinder surface should have"
      << std::endl << "    " << (expectSameSense ? "the same sense." :
      "different senses.") << std::endl;
  if (topCoEdgePtr->get_sense() == bottomCoEdgePtr->get_sense())
  {
    std::cout << "The senses are the same." << std::endl;
  }
  else
  {
    std::cout << "The senses are different." << std::endl;
  }

  // fail the test if necessary
  if (expectSameSense &&
      topCoEdgePtr->get_sense() != bottomCoEdgePtr->get_sense())
  {
    std::cout << "FAILED TEST." << std::endl;
    return false;
  }
  else if (!expectSameSense &&
      topCoEdgePtr->get_sense() == bottomCoEdgePtr->get_sense())
  {
    std::cout << "FAILED TEST." << std::endl;
    return false;
  }

  // pass the test
  std::cout << "PASSED TEST." << std::endl;
  return true;
}
