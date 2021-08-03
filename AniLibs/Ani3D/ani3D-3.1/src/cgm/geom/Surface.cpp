//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Surface.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Many
//
// Creation Date : 08/02/96
//
// Owner         : Timothy J. Tautges
//-------------------------------------------------------------------------


#include "Surface.hpp"
#include "RefFace.hpp"
#include "GeometryQueryEngine.hpp"

//-------------------------------------------------------------------------
// Purpose       : The default constructor. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------

Surface::Surface() 
{
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------

Surface::~Surface() 
{}

CubitStatus Surface::closest_points(DLIList<CubitVector *> &location_list,
                                    DLIList<CubitVector *> *closest_location_list,
                                    DLIList<CubitVector *> *unit_normal_list,
                                    DLIList<CubitVector *> *curvature1_list,
                                    DLIList<CubitVector *> *curvature2_list)
{
  CubitVector *curvature1, *curvature2;
  CubitVector *unit_normal;
  CubitVector *closest_location;
  CubitVector *location;
  CubitStatus stat;
  location_list.reset();
  if (closest_location_list) closest_location_list->reset();
  if (unit_normal_list) unit_normal_list->reset();
  if (curvature1_list) curvature1_list->reset();
  if (curvature2_list) curvature2_list->reset();
  for (int i=0; i<location_list.size(); i++)
  {
    location = location_list.get_and_step();
    if (closest_location_list == NULL)
      closest_location = NULL;
    else
      closest_location = closest_location_list->get_and_step();
    if (unit_normal_list == NULL)
      unit_normal = NULL;
    else
      unit_normal = unit_normal_list->get_and_step();
    if (curvature1_list == NULL)
      curvature1 = NULL;
    else
      curvature1 = curvature1_list->get_and_step();
    if (curvature2_list == NULL)
      curvature2 = NULL;
    else
      curvature2 = curvature2_list->get_and_step();
    stat = closest_point( *location, closest_location, unit_normal, curvature1, curvature2 );
    if (stat != CUBIT_SUCCESS)
      return stat;
  }
  return CUBIT_SUCCESS;
}

void Surface::closest_points_trimmed( std::vector<CubitVector> &from_points_list, 
                                      std::vector<CubitVector> &points_on_surface_list)
{
  CubitVector from_point;
  CubitVector point_on_surface;  
  for (unsigned int i=0; i<from_points_list.size(); i++)
  {
    from_point = from_points_list[i];    
    closest_point_trimmed( from_point, point_on_surface );
    points_on_surface_list.push_back( point_on_surface );
  }
}

void Surface::are_positions_on( DLIList<CubitVector *> &test_position_list,
                              DLIList<CubitBoolean *> &is_on_list )
{
  CubitVector *test_position;
  CubitBoolean *is_on;
  test_position_list.reset();
  is_on_list.reset();
  for (int i=0; i<test_position_list.size(); i++)
  {
    test_position = test_position_list.get_and_step();
    is_on = is_on_list.get_and_step();
    *is_on = is_position_on( *test_position );
  }
}


CubitStatus Surface::closest_point_along_vector(CubitVector& from_point, 
                                         CubitVector& along_vector,
                                         CubitVector& point_on_surface)
{
  return CUBIT_FAILURE;
}

CubitStatus Surface::evaluate( double u, double v, CubitVector& pos, CubitVector deriv1[2], CubitVector deriv2[3])
{
  PRINT_ERROR("evaluate is not implemented.");
  return CUBIT_FAILURE;
}
