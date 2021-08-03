//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Surface.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef SURFACE_HPP
#define SURFACE_HPP


#include <assert.h>


#include "CubitDefines.h"
#include "GeometryEntity.hpp"


class CubitVector;
template <class X> class DLIList;
class ShellSM;

class CUBIT_GEOM_EXPORT Surface : public GeometryEntity
{
   public :

      Surface() ;
      //- The default constructor

      virtual ~Surface() ;
      //- The destructor

      typedef Curve ChildType;
  
      virtual CubitSense get_shell_sense( ShellSM* shell_ptr ) const = 0;
  
      virtual void closest_point_trimmed(CubitVector from_point, 
                                         CubitVector& point_on_surface) = 0;
      //R void
      //I CubitVector
      //I- point from which to find closest point on trimmed surface
      //O CubitVector
      //O- point on trimmed surface closest to passed-in point
      //- This function finds the closest point on a TRIMMED surface to the
      //- passed-in point.

      virtual void closest_points_trimmed(std::vector<CubitVector> &from_points_list,
                                          std::vector<CubitVector> &points_on_surface_list);

      virtual CubitStatus closest_point_along_vector(CubitVector& from_point,
                                         CubitVector& along_vector,
                                         CubitVector& point_on_surface) = 0;

      virtual CubitStatus get_point_normal( CubitVector& point,
                                            CubitVector& normal ) = 0;
      //- Only valid for planar surfaces
      //- Finds the underlying plane's origin and normal (unit) vector
      //- Returns CubitFailure if not a plane.

      virtual void param_dir(CubitVector &/*unit_dir_in_world_space*/,
        CubitVector &/*pos_on_surf*/, double &/*du*/, double &/*dv*/){}

      virtual CubitStatus closest_point_uv_guess(  
          CubitVector const& location,
          double& u_guess, double& v_guess,
          CubitVector* closest_location = NULL,
          CubitVector* unit_normal = NULL ) = 0;
        //R CubitStatus
        //I location
        //I- Position to evaluate from.
        //I u, v
        //I- As input, hint as to location of result.
        //I- As output, u and v of the result.
        //I closest_location
        //I- If not null, the closest point to 'location' on the surface.
        //I unit_normal
        //I- If not null, set to normal at result position.
        //- Find closest point on surface, passing a parameter pair
        //- near the result point as a hint to the surface evaluator
        //- and passing back the paramter pair of the resuling 
        //- position on the surface.


      virtual CubitStatus closest_point(  
          CubitVector const& location, 
          CubitVector* closest_location = NULL,
          CubitVector* unit_normal = NULL,
          CubitVector* curvature1 = NULL,
          CubitVector* curvature2 = NULL) = 0;
      //R  CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //I location
      //I- The point to which the closest point on the surface is desired.
      //O closest_location
      //O- The point on the Surface, closest to the 
      //O- input location (which might not be on the Surface).  This is
      //O- input as a reference so that the function can modify its
      //O- contents.
/*BWC
It says that this is passed in as a reference but it isn't.  My guess is
that it should be a reference and not a pointer.  I will leave it for
now in case I am missing something, but, maybe it should be changed.
BWC*/
      //I refvolume_ptr
      //I- The first underlying geometric entity is used to compute the
      //I- normal.
      //O unit_normal
      //O- The normal (represented as a unit vector) at the closest_location.
      //O- If this pointer is NULL, the normal is not returned.
      //O curvature1
      //O- The first principal curvature of the Surface at closest_location.
      //O- If this pointer is NULL, this curvature is not returned.
      //O curvature2
      //O- The second principal curvature of the Surface at closest_location.
      //O- If this pointer is NULL, this curvature is not returned.
      //- This function computes the point on the Surface closest to the input 
      //- location -- i.e., closest_location. 
      //-
      //- The first ACIS FACE in the list
      //- is queried.
      //-
      //- If the input pointer values of unit_normal, curvature1 and curvature2
      //- are non-NULL, the normal and principal curvatures, too, are
      //- returned.  These are computed at closest_location, not at the
      //- input location.
      //-
      //- NOTE:
      //- It is assumed that if the calling code needs the normal or the 
      //- principal curvatures, it will *allocate* space for the CubitVectors
      //- before sending in the pointers.

      virtual CubitStatus closest_points(DLIList<CubitVector *> &location_list,
                                   DLIList<CubitVector *> *closest_location_list = NULL,
                                   DLIList<CubitVector *> *unit_normal_list = NULL,
                                   DLIList<CubitVector *> *curvature1_list = NULL,
                                   DLIList<CubitVector *> *curvature2_list = NULL);

      virtual CubitStatus principal_curvatures(
          CubitVector const& location, 
          double& curvature_1,
          double& curvature_2,
          CubitVector* closest_location = NULL ) = 0;
      //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //I location
      //I- The point at which the curvatures are being requested -- it is also
      //I- the point to which the closest point on the surface is returned.
      //I refvolume_ptr
      //O closest_location
      //O- The point on the surface, closest to the input location (this
      //O- might not be on the surface).  This is input as a reference 
      //O- so that the function can modify its contents.
      //O curvature_1/2
      //O- Returned principal curvature magnitudes.
      //- This functions computes the point on the surface that is closest
      //- to the input location and then calculates the magnitudes of the 
      //- principal curvatures at this (possibly, new) point on the surface. 


      virtual CubitStatus evaluate( double u, double v,
                            CubitVector *position,
                            CubitVector *normal,
                            CubitVector *curvature1,
                            CubitVector *curvature2 )=0;

      // evaluate parameters on a surface
      // returns the position on the surface
      // returns the first order derivatives: du, dv
      // returns the second order derivatives: duu, duv, dvv
      virtual CubitStatus evaluate( double u, double v, CubitVector& pos, CubitVector deriv1[2], CubitVector deriv2[3]);

     //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
           //I u
      //I- The u coordinate value (local parametric space).
      //I v
      //I- The v coordinate value (local parametric space).
      //O- The coordinates of a point in global
      //- (world) space that correspond to the input {u,v} point in
      //- local parametric space. This is an optional output that is 
      //- only calculated if the supplied CubitVector is not NULL;
      //O-The normal to the surface at the point with the given parametric position.
      //- This is an optional output that is only calculated if the 
      //- supplied CubitVector is not NULL;
      //O- The principle curvature directions at the point with the given parametric position.
      //- These are optional outputs that are only calculated if the 
      //- supplied CubitVectors curvature1 and curvature2 are both not NULL;


      virtual CubitVector position_from_u_v (double u, double v) = 0;
      //R CubitVector
      //R- The returned position in global space
      //I u
      //I- The u coordinate value (local parametric space).
      //I v
      //I- The v coordinate value (local parametric space).
      //- This function returns the coordinates of a point in global
      //- (world) space that correspond to the input {u,v} point in
      //- local parametric space.

      virtual CubitStatus u_v_from_position (
                                     CubitVector const& location,
                                     double& u, 
                                     double& v,
                                     CubitVector* closest_location = NULL ) = 0;
      //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //I location
      //I- The input point in global space
      //O closest_point
      //O- The point on the Surface closest to the input location
      //O u, v
      //O- The returned u, v coordinate values (in local parametric space)
      //O- of the closest_point
      //- This function returns the {u, v} coordinates of the point 
      //- on the Surface closest to the input point (specified in global
      //- space). The closest_location is also returned.

      virtual CubitBoolean is_periodic() = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //- This function determines whether the underlying geometry of the
      //- Surface is periodic or not.  Returns CUBIT_TRUE if it is and 
      //- CUBIT_FALSE if it is not.

      virtual CubitBoolean is_periodic_in_U( double& period ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //O period
      //O- The value of the period in the U direction.
      //- Determines whether the surface object is 
      //- periodic in the U direction or not.  If it is, it
      //- returns CUBIT_TRUE and the value of the period. Otherwise,
      //- it returns CUBIT_FALSE and a value of 0.0 or the period.

      virtual CubitBoolean is_periodic_in_V( double& period ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //O period
      //O- The value of the period in the V direction.
      //- Determines whether the surface object is 
      //- periodic in the V direction or not.  If it is, it
      //- returns CUBIT_TRUE and the value of the period. Otherwise,
      //- it returns CUBIT_FALSE and a value of 0.0 or the period.

      virtual CubitBoolean is_singular_in_U( double u_param ) = 0;
      virtual CubitBoolean is_singular_in_V( double v_param ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //I double u parameter value.
      //- Determines if the surface is singular in a given direction
      //- at a given parameter value.

      virtual CubitBoolean is_closed_in_U(){return CUBIT_FALSE;}
      virtual CubitBoolean is_closed_in_V(){return CUBIT_FALSE;}
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //- Determines if the surface is closed, smoothly or not in the
      //- given parameter direction.
      //- A periodic surface is always closed but a closed surface is
      //- is not always periodic.
      //- For modelars that do allow closed surfaces they will need
      //- to implement this functionality.  For instance, Pro does
      //- not allow such surfaces so it will always return false.
      //- ACIS will so SurfaceACIS implements its own function.
  
      virtual CubitStatus uv_derivitives( double u_param,
                                          double v_param,
                                          CubitVector &du,
                                          CubitVector &dv ) = 0;
      //R CubitStatus
      //R- CUBIT_SUCCESS/CUBIT_FAILURE
      //O- du, dv
      //- Determines the u and v derivitives from the given parameter
      //- values.

      virtual CubitBoolean is_parametric() = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //- This function determines whether the underlying geometry of the
      //- Surface is parametrically defined or not.  Returns CUBIT_TRUE if 
      //- it is and CUBIT_FALSE if it is not.

 

      virtual CubitBoolean get_param_range_U( double& lower_bound,
                                              double& upper_bound ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //O lower_bound
      //O- The lower bound of the parametric range in the U direction.
      //O- This is set to 0.0 if the surface is not parametric.
      //O upper_bound
      //O- The upper bound of the parametric range in the U direction.
      //O- This is set to 0.0 if the surface is not parametric.
      //- Returns the lower and upper parametric bounds of the 
      //- surface in U, if it is parametric.  Otherwise, it returns
      //- CUBIT_FALSE and zeroes for the upper and lower parametric
      //- bounds.

      virtual CubitBoolean get_param_range_V( double& lower_bound,
                                              double& upper_bound ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //O lower_bound
      //O- The lower bound of the parametric range in the V direction.
      //O- This is set to 0.0 if the surface is not parametric.
      //O upper_bound
      //O- The upper bound of the parametric range in the V direction.
      //O- This is set to 0.0 if the surface is not parametric.
      //- Returns the lower and upper parametric bounds of the 
      //- surface in V, if it is parametric.  Otherwise, it returns
      //- CUBIT_FALSE and zeroes for the upper and lower parametric
      //- bounds.

      virtual CubitBoolean is_position_on( CubitVector &test_position ) = 0;
      //R CubitBoolean
      //R- CUBIT_TRUE/CUBIT_FALSE
      //I CubitVector
      //I- position, point where we want to test, whether or not it
      //- is on the surface.
      virtual void are_positions_on( DLIList<CubitVector *> &test_position_list,
                                     DLIList<CubitBoolean *> &is_on_list );

      virtual CubitPointContainment point_containment( const CubitVector &point ) = 0;
      virtual CubitPointContainment point_containment( double u, double v ) = 0;
      //R CubitPointContainment - is the point outside, inside or on the boundary?
      //R- CUBIT_PNT_OUTSIDE, CUBIT_PNT_INSIDE, CUBIT_PNT_BOUNDARY, 
      //   CUBIT_PNT_UNKNOWN
      //I CubitVector
      //I- position to check, known to be on the Surface
      //I double
      //I- u coordinate, if known (significantly faster, if this is known - however
      //                           if not known let the function figure it out)
      //I double
      //I- v coordinate, if known (significantly faster, if this is known - however
      //                           if not known let the function figure it out)

      virtual CubitSense get_geometry_sense() = 0;
      //- Returns the relative sense of the Surface with respect
      //- to the geometry underneath.  This is geometry engine dependent.
      //- Currently this is used for the tet mesher...

      virtual CubitSense get_uv_sense()
      {
        return CUBIT_UNKNOWN;
      }

      virtual CubitStatus get_projected_distance_on_surface( CubitVector *pos1,
                                                             CubitVector *pos2, 
                                                             double &distance ) = 0;
      

      virtual GeometryType geometry_type()
      {return UNDEFINED_SURFACE_TYPE;};
      //R GeometryType (enum)
      //R- The enumerated type of the geometric representation

      virtual GeometryType is_cylindrical()
      {return UNDEFINED_SURFACE_TYPE;};
      //R GeometryType (enum) 
      //R- CYLINDRICAL_SURFACE_TYPE if true UNDEFINED_SURFACE_TYPE if false

      // Now handled at RefFace level.  -- j.k. Oct, 2003
      //virtual void reverse_sense() = 0;
      //- Switch the sense of this Surface wrt the RefFace that owns it.

      virtual CubitStatus get_nurb_params( bool &rational,
                                           int &degree_u,
                                           int &degree_v,
                                           int &num_cntrl_pts_u,
                                           int &num_cntrl_pts_v,
                                           DLIList<CubitVector> &cntrl_pts,
                                           DLIList<double> &weights,
                                           DLIList<double> &u_knots,
                                           DLIList<double> &v_knots ) const = 0;
      //- Only valid for nurbs surfaces
      //O rational
      //O-   True if the nurb is rational
      //O degree_u
      //O-   The degree of the nurb in the u direction
      //O degree_v
      //O-   The degree of the nurb in the v direction
      //O num_cntrl_pts_u
      //O-   Number of control points in the u direction
      //O num_cntrl_pts_v
      //O-   Number of control points in the v direction
      //O cntrl_pts
      //O-   The control points stored as
      //O-           cntrl_pts[0                ] = pt[u=0][v=0]
      //O-           cntrl_pts[1                ] = pt[u=1][v=0]
      //O-               ...
      //O-           cntrl_pts[num_cntrl_pts_u-1] = pt[u=?][v=0]
      //O-           cntrl_pts[num_cntrl_pts_u  ] = pt[u=0][v=1]
      //O-               ...
      //O weights
      //O-   If rational, weights for each control point, stored in the same
      //O-   order as the control points.  No weights are returned if
      //O-   rational == false
      //O u_knots
      //O-   knot vector in the u direction
      //O v_knots
      //O-   knot vector in the v direction

      virtual CubitStatus get_sphere_params( CubitVector &center, 
                                             double &radius ) const = 0;
      //- Only valid for spherical surfaces
      //O center
      //O- The center of the sphere
      //O radius
      //O- The radius of the sphere

      virtual CubitStatus get_cone_params( CubitVector &center,
                                           CubitVector &normal,
                                           CubitVector &major_axis,
                                           double &radius_ratio,
                                           double &sine_angle,
                                           double &cos_angle ) const = 0;
      //- Only valid for conical surfaces.  Cylinders are a special case of conicals.
      //O center
      //O- 
      //O normal
      //O- 
      //O major_axis
      //O- 
      //O radius_ratio
      //O- 
      //O sine_angle
      //O- 
      //O cos_angle
      //O- 

      virtual CubitStatus get_torus_params( CubitVector &center,
                                            CubitVector &normal,
                                            double &major_radius,
                                            double &minor_radius ) const = 0;
      //- Only valid for torus surfaces.
      //O center
      //O- 
      //O normal
      //O- 
      //O major_radius
      //O- 
      //O minor_radius
      //O- 

   protected: 

   private:
  
} ;


#endif

