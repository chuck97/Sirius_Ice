//-------------------------------------------------------------------------
// Filename      : FacetSurface.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 06/06/00
//
// Owner         : David R. White
//-------------------------------------------------------------------------

#ifndef SURFACE_Facet_HPP
#define SURFACE_Facet_HPP

#include "CubitDefines.h"
#include "Surface.hpp"
#include "FacetAttribSet.hpp"
#include <map>

class TopologyEntity;
class RefVolume;
class RefFace;
class RefVolume;
class FacetEvalTool;
class FacetShell;
class FacetAttrib;

class FacetBody;
class FacetLump;
class FacetLoop;
class FacetCoEdge;
class FacetCurve;
class FacetPoint;
class CubitFacetEdge;
class CubitFacet;
class CubitPoint;
class CubitTransformMatrix;
class CubitEvaluator;
class CubitEvaluatorData;
class SphereEvaluatorData;
class CylinderEvaluatorData;

class FacetSurface : public Surface
{

private:
    //CubitSense sense_;
    //- The sense of the RefFace that owns this Surface with respect
    //- to the positive sense of the first FACE in FACEPtrList_.
    //- When a Surface is first constructed, this value is arbitrarily
    //- set to CUBIT_FORWARD.
    //- In the case of Surfaces, the normal is used in determining
    //- the relative sense value.
    //- MJP NOTE:
    //- Not only does the RefFace have a sense wrt its Surface, but each
    //- Facet FACE has a sense wrt its underlying "surface" object.

  //sjowen FacetEvalTool *facetEvalTool;
  FacetEvalTool *facetEvalTool;
    //For topology traversals we need to have connections to the
    //entitities connected to this surface.  From these we can get the rest.
  DLIList<LoopSM*> myLoops;
  DLIList<ShellSM*> myShells;

  FacetAttribSet attribSet;
    //List of FacetAttrib*'s instead of CubitSimpleAttribs 

  // the sense of the surface with respect to the shells in the myShells list
  CubitSense myShellSense;

  CubitEvaluator *myEvaluator;

protected: 

public :
  
  FacetSurface(FacetEvalTool *facet_eval_tool_ptr,
               DLIList<ShellSM*> &shellsms,DLIList<LoopSM*> &loopsms );
    //I- facet_eval_tool pointer
    //I- A pointer to the set of facets that define this surface.
 
  FacetSurface(FacetEvalTool *facet_eval_tool_ptr,
               CubitSense sense,
               CubitSense shell_sense0,
               CubitBoolean use_facets,
               DLIList<LoopSM*> &loopsms );


  FacetSurface( const CylinderEvaluatorData *cylinder_data,
                FacetEvalTool *facet_tool,
                DLIList<ShellSM*> &shellsms,
                DLIList<LoopSM*> &loopsms );
  //-  Constructor used to create a faceted surface representing a cylinder.
  //I-  eval_data - radius, base, etc. of cylinder.
  //I-  facet_eval_tool_ptr - evaluator to evaluate directly on facets.
  //I-  shellsms - the shells in this facet model
  //I-  loopsms - the loops in this facet model.

  FacetSurface( const SphereEvaluatorData *eval_data,
                FacetEvalTool *facet_eval_tool_ptr,
                DLIList<ShellSM*> &shellsms,
                DLIList<LoopSM*> &loopsms );
  //-  Constructor used to create a faceted surface representing a sphere.
  //I-  eval_data - radius, center, etc. of sphere.
  //I-  facet_eval_tool_ptr - evaluator to evaluate directly on facets.
  //I-  shellsms - the shells in this facet model
  //I-  loopsms - the loops in this facet model.

  virtual ~FacetSurface() ;
    //- The destructor
  
  void add_shell(ShellSM *shell_ptr)
    {myShells.append(shell_ptr);}
    
  CubitStatus remove_shell(FacetShell* shell_ptr);
  
  void disconnect_all_loops();
  
  inline bool has_parent_shell() { return myShells.size() > 0; }
      
    //CubitSense get_relative_surface_sense();
    //- Return the relative surface sense. (see below)

  virtual void append_simple_attribute_virt(const CubitSimpleAttrib&);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(const CubitSimpleAttrib&);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
  
  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  // Added by CAT
  virtual CubitStatus get_point_normal( CubitVector& ,
                                        CubitVector& );
    //- Only valid for planar surfaces
    //- Finds the underlying plane's origin and normal vector
    //- Returns CubitFailure if not a plane.  The origin and normal 
    //- are returned directly from the underlying format for a plane.
  
  virtual CubitStatus closest_point_uv_guess(  
      CubitVector const& location,
      double &u, double &v,
      CubitVector* closest_location = NULL,
      CubitVector* unit_normal = NULL );


  virtual void closest_point_trimmed(CubitVector from_point,
                                     CubitVector& point_on_surface);
    //R void
    //I CubitVector
    //I- point from which to find closest point on trimmed surface
    //O CubitVector
    //O- point on trimmed surface closest to passed-in point
    //- This function finds the closest point on a TRIMMED surface to the
    //- passed-in point.

  virtual CubitStatus closest_point_along_vector(CubitVector& from_point, 
                                         CubitVector& along_vector,
                                         CubitVector& point_on_surface);

  
  virtual CubitStatus closest_point(  
    CubitVector const& location, 
    CubitVector* closest_location = NULL,
    CubitVector* unit_normal_ptr = NULL,
    CubitVector* curvature1_ptr = NULL,
    CubitVector* curvature2_ptr = NULL);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point to which the closest point on the surface is desired.
    //O closest_location
    //O- The point on the Surface, closest to the 
    //O- input location (which might not be on the Surface).  This is
    //O- input as a reference so that the function can modify its
    //O- contents.
    //O unit_normal_ptr
    //O- The normal (represented as a unit vector) at the closest_location.
    //O- If this pointer is NULL, the normal is not returned.
    //O curvature1_ptr
    //O- The first principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //O curvature2_ptr
    //O- The second principal curvature of the surface at closest_location.
    //O- If this pointer is NULL, this curvature is not returned.
    //- This function computes the point on the surface closest to the input 
    //- location -- i.e., closest_location. 
    //- The first Facet FACE in the list
    //- is queried.
    //-
    //- If the input pointer values of unit_normal, curvature1 and
    //- curvature2
    //- are non-NULL, the normal and principal curvatures, too, are
    //- returned.  These are computed at closest_location, not at the
    //- input location.
    //-
    //- NOTE:
    //- It is assumed that if the calling code needs the normal or the 
    //- principal curvatures, it will *allocate* space for the CubitVectors
    //- before sending in the pointers.
  
  virtual CubitStatus principal_curvatures(
    CubitVector const& location, 
    double& curvature_1,
    double& curvature_2,
    CubitVector* closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The point at which the curvatures are being requested -- it is also
    //I- the point to which the closest point on the surface is returned.
    //I- curvatures.
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
    CubitVector *curvature2 );



  virtual CubitVector position_from_u_v (double u, double v);
    //R CubitVector
    //R- Returned position vector.
    //I u, v
    //I- Input point in {u.v} space
    //- This function returns the coordinates in world space of a point
    //- in the parameter space of this Surface object.
  
  virtual CubitStatus u_v_from_position (CubitVector const& location,
                                         double& u, 
                                         double& v,
                                         CubitVector*
                                         closest_location = NULL );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The input point in global space
    //O closest_point
    //O- The point on the Surface closest to the input location
    //O u, v
    //O- The returned u, v coordinate values (in local parametric space)
    //O- of the closest_point
    //I refvolume_ptr
    //- This function returns the {u, v} coordinates of the point 
    //- on the Surface closest to the input point (specified in global
    //- space). The closest_location is also returned.
  
  virtual CubitBoolean is_periodic();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- FacetSurface is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
    virtual CubitBoolean is_periodic_in_U( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the Facet surface object associated
    //- with one of the FACEs of this FacetSurface object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_periodic_in_V( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the Facet surface object associated
    //- with one of the FACEs of this FacetSurface object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_singular_in_U( double u_param );
  virtual CubitBoolean is_singular_in_V( double v_param );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I double u/v parameter value.
    //- Determines if the surface is singular in a given direction
    //- at a given parameter value.
  
  virtual CubitBoolean is_closed_in_U();  
  virtual CubitBoolean is_closed_in_V();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Determines if the surface is closed, smoothly or not in the
    //- given parameter direction.
    //- A periodic surface is always closed but a closed surface is
    //- is not always periodic.
  
  virtual CubitStatus uv_derivitives( double u_param,
                                      double v_param,
                                      CubitVector &du,
                                      CubitVector &dv );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O- du, dv
    //- Determines the u and v derivitives from the given parameter
    //- values.
  
  virtual CubitBoolean is_parametric();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
 			//- This method returns CUBIT_TRUE if a parametric representation
			//- is available for the surface
  
  virtual CubitBoolean get_param_range_U( double& lower_bound,
                                          double& upper_bound );
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
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean get_param_range_V( double& lower_bound,
                                          double& upper_bound );
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
    //- MJP NOTE: 
    //- The first Facet FACE in the list is queried.  It is assumed
    //- that all the FACEs have the same underlying surface.
  
  virtual CubitBoolean is_position_on( CubitVector &test_position );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I CubitVector
    //I- position, point where we want to test, whether or not it
    //- is on the surface.

  virtual CubitPointContainment point_containment( const CubitVector &point );
  virtual CubitPointContainment point_containment( double u, double v );
//  virtual CubitPointContainment point_containment( const CubitVector &point, 
//                                                   double u, double v );
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
    // NOTE: POINT MUST LIE ON THE SURFACE FOR THIS FUNCTION TO WORK PROPERLY.
  
  GeometryType geometry_type();
    //R GeometryType (enum)
    //R- The enumerated type of the geometric representation
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point

  
  virtual CubitStatus get_projected_distance_on_surface( CubitVector *pos1,
                                                         CubitVector *pos2, 
                                                         double &distance );

  void update_measurement();
    //Make sure we don't retain an out-dated measurement.
  
  virtual CubitSense get_geometry_sense();
    //- Return the relative surface sense. (see below)
  
  virtual void reverse_sense();
    //- Switch the sense of this Surface wrt the RefFace that owns it:
    //- For Facet, this means switch the sense of the RefFace that
    //- owns this Surface with respect to the positive sense of the
    //- first FACE in FACEPtrList_.
  
  CubitStatus save_attribs( FILE* file_ptr );
    // Write FactAttribs out to file

  CubitStatus restore_attribs( FILE* file_ptr, unsigned int endian );
    // Read FactAttribs from file
  
  void get_bodies  ( DLIList<FacetBody   *>& bodies   );
  void get_lumps   ( DLIList<FacetLump   *>& lumps    );
  void get_shells  ( DLIList<FacetShell  *>& shells   );
  void get_loops   ( DLIList<FacetLoop   *>& loops    );
  void get_coedges ( DLIList<FacetCoEdge *>& coedges  );
  void get_curves  ( DLIList<FacetCurve  *>& curves   );

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  CubitStatus get_my_facets(DLIList<CubitFacet*>& facet_list,
                            DLIList<CubitPoint*>& point_list);
    //- Gets the list of facets describing this surface.
  void tris(DLIList<CubitFacet*> &facet_list);
  void get_my_points(DLIList<CubitPoint*>& point_list);
    //- Gets the list of points describing this surface.
  void get_my_facetedges(DLIList<CubitFacetEdge*>& edge_list);
    //- Gets the list of points describing this surface.
  FacetEvalTool *get_eval_tool()
    { return facetEvalTool; }
  const FacetEvalTool *get_eval_tool() const
    { return facetEvalTool; }
    //- return the facet evaluation tool for this surface
    
  CubitSense get_shell_sense( ShellSM *facet_shell ) const;
    // return the sense with respect to the given shell
  
  void get_shell_sense( CubitSense &sense0 ); 
    // return senses 
    
  void set_shell_sense( FacetShell *facet_shell, 
                        CubitSense thesense );
    // set the sense of the surface with respect to the shell

  CubitStatus copy_facets(DLIList<CubitFacet*>&copy_facet_list,
                          DLIList<CubitPoint*>&copy_point_list,
                          std::map<CubitPoint*, CubitPoint*>&old_to_new_cubit_pts);
    // create a copy of the points and facets

  int interp_order();
  double min_dot();

  CubitBoolean is_flat(); 
  CubitBoolean is_spherical();
  CubitBoolean is_conical();

  virtual CubitStatus get_nurb_params( bool &rational,
                                       int &degree_u,
                                       int &degree_v,
                                       int &num_cntrl_pts_u,
                                       int &num_cntrl_pts_v,
                                       DLIList<CubitVector> &cntrl_pts,
                                       DLIList<double> &weights,
                                       DLIList<double> &u_knots,
                                       DLIList<double> &v_knots ) const;
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
                                         double &radius ) const;
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
                                       double &cos_angle ) const;
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
                                        double &minor_radius ) const;
    //- Only valid for torus surfaces.
    //O center
    //O- 
    //O normal
    //O- 
    //O major_radius
    //O- 
    //O minor_radius
    //O- 

  const CubitEvaluatorData *evaluator_data();
  void add_transformation( CubitTransformMatrix &tfmat );

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

