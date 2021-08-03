/**
 * \file split_circle.cpp
 *
 * \brief split_circle, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitMessage.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CastTo.hpp"
#include "InitCGMA.cpp"
#include "CubitCompat.hpp"
#include "Point.hpp"
#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus split_circle();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma("OCC");
  if (CUBIT_SUCCESS != status) return 1;

  status = split_circle();
  if (status == CUBIT_FAILURE) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val != 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  else
    ret_val = 0;

  return ret_val;
  
}

CubitStatus split_circle()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();
  CubitVector center(0,0,0);
  CubitVector plane (0,0,1);
  RefEdge* re = gmti->create_arc(center, 3, 0, 360, plane);
  double length = re->measure();
  std::cout << "Cicle's parameter = " << length << "\n";
/*
  CubitVector p1(0,3,0), p2(3, 0, 0), p3(0, -3, 0);
 
  RefVertex* v1 = gmti->make_RefVertex(p1);
  RefVertex* v2 = gmti->make_RefVertex(p2);
  RefVertex* v3 = gmti->make_RefVertex(p3);
  RefEdge* re = gmti->create_arc_three(v1, v2, v3, true);
*/ 
  //Curve* curve = re->get_curve_ptr() ; 
  CubitVector point(0, 3, 0);
  DLIList<CubitVector> sp_pts;
  sp_pts.append(point);
  DLIList<RefEdge*> new_ref_edges;
  gmti->split_free_curve(re, sp_pts, new_ref_edges);

  DLIList<RefEntity*>  free_entities;
  gti->get_free_ref_entities(free_entities);
  int count = 0;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* ref_edge = dynamic_cast<RefEdge*> (free_entities.get_and_step());
    if (ref_edge)
    {
      double d = ref_edge->measure();
      std::cout << "edge length = " << d << "\n";
      count++;
    }
  }
  assert( count == 2);
  std::cout << "number of free edges: " << count << " <---should be 2 free edges."<<"\n";

  return CUBIT_SUCCESS;  
}
