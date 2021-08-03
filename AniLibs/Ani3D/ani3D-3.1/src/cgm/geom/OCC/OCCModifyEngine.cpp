//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Filename      : OCCModifyEngine.cpp
//
// Purpose       : ModifyEngine for OCC geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Author        : Jane Hu
//
// Creation Date : 1/08
//
//-------------------------------------------------------------------------
#include "gp_Pnt.hxx"
#include "gp_Ax2.hxx"
#include "gp_Dir.hxx"
#include "gp_Hypr.hxx"
#include "gp_Parab.hxx"
#include "gp_Elips.hxx"
#include "gp_Pln.hxx"
#include "gp_Circ.hxx"
#include "gp_Cylinder.hxx"
#include "gp_Cone.hxx"
#include "gp_Sphere.hxx"
#include "gp_Torus.hxx"
#include "BRepOffsetAPI_MakeOffset.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_MakeDraft.hxx"
#include "BRepBuilderAPI_TransitionMode.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepPrimAPI_MakeHalfSpace.hxx"
#include "BRepBuilderAPI_MakePolygon.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepBndLib.hxx"
#include "IntersectionTool.hpp"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Compound.hxx"
#include "TopAbs_Orientation.hxx"
#include "TopOpeBRep_Point2d.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "TColgp_HArray1OfPnt.hxx"
#include "TColStd_HArray1OfBoolean.hxx"
#include "TColgp_Array1OfVec.hxx"
#include "TColStd_Array1OfReal.hxx"
#include "TColStd_Array1OfInteger.hxx"
#include "GC_MakeArcOfCircle.hxx"
#include "GC_MakeCircle.hxx"
#include "Geom_Circle.hxx"
#include "Geom_SurfaceOfLinearExtrusion.hxx"
#include "Geom_RectangularTrimmedSurface.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Handle_Geom_RectangularTrimmedSurface.hxx"
#include "GC_MakeArcOfHyperbola.hxx"
#include "GC_MakeArcOfParabola.hxx"
#include "GC_MakeArcOfEllipse.hxx"
#include "GC_MakeSegment.hxx"
#include "GC_MakeTrimmedCone.hxx"
#include "GC_MakeTrimmedCylinder.hxx"
#include "gce_MakeElips.hxx"
#include "BRepFilletAPI_MakeFillet.hxx"
#include "BRepFilletAPI_MakeChamfer.hxx"
#include "BRepAdaptor_CompCurve.hxx"
#include "GeomAPI_Interpolate.hxx"
#include "BRepFilletAPI_MakeFillet2d.hxx"
#include "ChFi2d_ConstructionError.hxx"
#include "Geom_BezierCurve.hxx"
#include "BndLib_AddSurface.hxx"
#include "Handle_Geom_Plane.hxx"
#include "Handle_Geom_OffsetCurve.hxx"
#include "Geom_OffsetCurve.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "Extrema_ExtPC.hxx"
#include "BRepPrimAPI_MakePrism.hxx"
#include "BRepSweep_Revol.hxx"
#include "BRepPrimAPI_MakeCone.hxx"
#include "BRepOffsetAPI_ThruSections.hxx"
#include "BRepLib_FuseEdges.hxx"
#include "BRepOffsetAPI_MakePipe.hxx"
#include "BRepPrimAPI_MakeTorus.hxx"
#include "BRepPrimAPI_MakeCylinder.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GC_MakeEllipse.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "ShapeExtend_Status.hxx"
#include "BRepOffsetAPI_MakeThickSolid.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepBuilderAPI_Copy.hxx"
#include "LocOpe_SplitShape.hxx"
#include "BRep_Tool.hxx"
#include "BRep_Builder.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "TopoDS.hxx"
#include "ShapeFix.hxx"
#include "TopologyBridge.hpp"
#include "ProgressTool.hpp"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepAlgoAPI_Common.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRepPrimAPI_MakeBox.hxx"
#include "BRepPrimAPI_MakeWedge.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "Handle_Geom_TrimmedCurve.hxx"
#include "Handle_ShapeBuild_ReShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "Handle_Geom_RectangularTrimmedSurface.hxx"
#include "BndLib_Add3dCurve.hxx"
#include "TopOpeBRep_EdgesIntersector.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"
#ifndef OCC_VERSION_MINOR
#include "Standard_Version.hxx"
#endif

#include "OCCDrawTool.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "TopOpeBRep_ShapeIntersector.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "CubitUtil.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"
#include "OCCAttribSet.hpp"
#include "CubitFileIOWrapper.hpp"
#include "Body.hpp"
#include "GfxDebug.hpp"
#include "RefFace.hpp"
#include "CpuTimer.hpp"
#include "AppUtil.hpp"
#include "SphereEvaluator.hpp"
#include "CylinderEvaluator.hpp"
#include "GfxPreview.hpp"
#include <vector>
#include "CGMEngineDynamicLoader.hpp"

CGM_ENGINE_EXPORT_CREATE_GME(OpenCascade)
{
  return OCCModifyEngine::instance();
}

OCCModifyEngine* OCCModifyEngine::instance_ = 0;
typedef std::map<OCCSurface*, std::pair<CubitVector, int> >::value_type valType;
typedef std::map<OCCCurve*, std::pair<CubitVector, int> >::value_type valType2;
#define DEBUG
double TOL = 0.0;
//===============================================================================
// Function   : OCCModifyEngine
// Member Type: PUBLIC
// Description: constructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::OCCModifyEngine()
{
//  assert( !instance_ );

    // add this modify engine to geometrymodifytool
  GeometryModifyTool::instance()->add_gme(this);
  TOL = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
}


//===============================================================================
// Function   : ~OCCModifyEngine
// Member Type: PUBLIC
// Description: destructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::~OCCModifyEngine() 
{
        instance_ = 0;
}

//===============================================================================
// Function   : make_Point
// Member Type: PUBLIC
// Description: make a geometric entity point
// Author     : Jane Hu 
// Date       : 10/07
//===============================================================================
TBPoint* OCCModifyEngine::make_Point( CubitVector const& point) const
{
  gp_Pnt pt = gp_Pnt( point.x(), point.y(), point.z());
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);

  // Create a new PointACIS object
  return OCCQueryEngine::instance()->populate_topology_bridge( theVertex, true );
}

//===============================================================================
// Function   : make_Curve
//              This function creates a curve given an existing curve, copy. 
// Member Type: PUBLIC
// Description: make a curve
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve(Curve * curve_ptr, std::map<TopologyBridge*, TopologyBridge*> * /*old_tb_to_new_tb*/) const
{
  OCCCurve* occ_curve = CAST_TO(curve_ptr, OCCCurve);
  if (!occ_curve)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
 
  TopoDS_Edge *theEdge = occ_curve->get_TopoDS_Edge();  
 
  BRepBuilderAPI_Copy api_copy(*theEdge);

  TopoDS_Shape newShape = api_copy.ModifiedShape(*theEdge);
 
  TopoDS_Edge newEdge = TopoDS::Edge(newShape);

  OCCQueryEngine::instance()->copy_attributes(*theEdge, newEdge);

  return OCCQueryEngine::instance()->populate_topology_bridge(newEdge, true);
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve by projecting a straight line defined by 
//              point1_ptr, and point2_ptr onto face_ptr, third_point
//              is used for curves that could be periodic to dertermine
//              the correct direction.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( TBPoint const* point1_ptr,
                             TBPoint const* point2_ptr,
                             Surface* face_ptr,
                             const CubitVector * third_point) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  GeometryType type = STRAIGHT_CURVE_TYPE;
  CubitBoolean closed = CUBIT_FALSE;
  DLIList<CubitVector*> mid_points;
  Curve* curve = NULL;
  if (point1_ptr != point2_ptr)
    curve = make_Curve(type, point1_ptr, point2_ptr, mid_points);
  else //could be a closed shape
  {
    if(third_point != NULL && face_ptr != NULL) 
    {
       closed = CUBIT_TRUE;
       TBPoint * Pnt = make_Point(*third_point);
       curve = make_Curve(type, point1_ptr, Pnt, mid_points);
    }
    else
    {
       PRINT_ERROR("Cannot create an OCC curve from the given duplicated points.\n");
       return (Curve *)NULL;
    }
  }

  Curve* new_curve = NULL;
  if(face_ptr == NULL)
    return curve;
 
  DLIList<TBPoint*> points;
  
  new_curve = 
    CAST_TO(curve, OCCCurve)->
           project_curve(face_ptr, points, closed, third_point);
  OCCQueryEngine::instance()->delete_solid_model_entities( curve );
  return new_curve;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a  spline curve by using the points and tangents.
//              size of point list and tangents must be the same.
//              values in the tangent list may be null.
// Author     : Jane Hu
// Date       : 01/11
//===============================================================================
Curve* OCCModifyEngine::make_Curve( DLIList<CubitVector*>& point_list,
                             DLIList<CubitVector*>& point_tangents) const
{
    if (point_list.size() != point_tangents.size())
    {
        PRINT_ERROR("    point list and tangent list must have same size\n");
        return (Curve *)NULL;
    }

    if(point_list.size() < 2)
    {
        PRINT_ERROR("    Can't create a curve with less than 2 points\n");
        return (Curve *)NULL;
    }

    int size = point_list.size();
    Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, size);
    TColgp_Array1OfVec tangents(1, size);
    Handle(TColStd_HArray1OfBoolean) tangentFlags = 
                                     new TColStd_HArray1OfBoolean(1,size);
    gp_Pnt pt, pt1;
    CubitVector *pt_vec, *tangent_vec;
    gp_Vec tangent;
    for (int i = 1 ; i <= size; i++)
    {
        pt_vec = point_list.get_and_step();
        pt.SetCoord(pt_vec->x(), pt_vec->y(), pt_vec->z()); 
        if(i == 1)
          pt1 = pt; 
        points->SetValue(i, pt);

        tangent_vec = point_tangents.get_and_step();
        if (!tangent_vec)
        {
          tangents.SetValue(i,tangent);
          tangentFlags->SetValue(i, CUBIT_FALSE);
        }
        else 
        {
          tangent.SetCoord(tangent_vec->x(), tangent_vec->y(),
                           tangent_vec->z());
          tangents.SetValue(i,tangent);
          tangentFlags->SetValue(i, CUBIT_TRUE);
        }
    }
    GeomAPI_Interpolate interpolater(points, CUBIT_FALSE, TOL);
    interpolater.Load(tangents, tangentFlags);
    interpolater.Perform() ;
    Handle(Geom_BSplineCurve) pcurve;
    if(interpolater.IsDone())
      pcurve = interpolater.Curve();
    
    else
    {
      PRINT_ERROR("Can't create a curve using provided points and tangents.\n");
      return (Curve *)NULL;
    }
   
    TopoDS_Edge topo_edge = BRepBuilderAPI_MakeEdge(pcurve, pt1, pt);
    return OCCQueryEngine::instance()->populate_topology_bridge(topo_edge, true);
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a  spline curve by using the points on surface.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             TBPoint const* point1_ptr,
                             TBPoint const* point2_ptr,
                             DLIList<CubitVector*>& vector_list,
                             Surface* face_ptr) const
{
  assert(point1_ptr != NULL && point2_ptr != NULL);
  
  if (curve_type != SPLINE_CURVE_TYPE
      && curve_type != STRAIGHT_CURVE_TYPE)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve_type.\n"
                 "Candidates are SPLINE_CURVE_TYPE and STRAIGHT_CURVE_TYPE.\n");
     return (Curve *)NULL;
  }

  if (curve_type == STRAIGHT_CURVE_TYPE)
    return make_Curve(curve_type, point1_ptr, point2_ptr, NULL);

  OCCPoint* occ_point1 = CAST_TO(const_cast<TBPoint*>(point1_ptr), OCCPoint);
  OCCPoint* occ_point2 = CAST_TO(const_cast<TBPoint*>(point2_ptr), OCCPoint);

  if (occ_point1 == NULL || occ_point2 == NULL)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given points.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
    
  //project all points on the surface if possible
  OCCSurface* occ_face = NULL;
  if (face_ptr != NULL)
     occ_face = CAST_TO(face_ptr, OCCSurface);
 
  gp_Pnt pt;
  int size = vector_list.size();
  Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, size);
  CubitVector vector;
  CubitVector closest_location;
  for(int i = 1; i <= size; i++)
  {
     vector = *vector_list.get_and_step();
     pt.SetCoord(vector.x(), vector.y(), vector.z());

     if (occ_face != NULL)
     {
       occ_face->closest_point(vector, &closest_location);
       pt.SetCoord(closest_location.x(), closest_location.y(), closest_location.z()) ;  	 
     }

     points->SetValue(i, pt);
  }    
     
  //make curve according to the curve type.
  if(curve_type == SPLINE_CURVE_TYPE)
  {
    if (size < 3)
    {
      PRINT_ERROR(" Must have at least 3 points to make a spline. \n");
      return (Curve*) NULL;
    }
  
    GeomAPI_Interpolate spline(points, CUBIT_FALSE, TOL);
    spline.Perform();
    if(spline.IsDone())
    {
      Handle_Geom_BSplineCurve curve = spline.Curve();
      TopoDS_Vertex * vt1 = occ_point1->get_TopoDS_Vertex();
      TopoDS_Vertex * vt2 = occ_point2->get_TopoDS_Vertex(); 
      TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve, *vt1, *vt2);
      return OCCQueryEngine::instance()->populate_topology_bridge(new_edge, true); 
    }
  }

  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// For STRAIGHT_CURVE_TYPE:
//    intermediate_point_ptr  is not used
//
// For PARABOLA_CURVE_TYPE
//    intermediate_point_ptr is the tip of the parabola
//
// For HYPERBOLA_CURVE_TYPE
//    intermediate_point_ptr is the center of its two foci
//
// For ELLIPSE_CURVE_TYPE
//    intermediate_point_ptr is the center of the ellipse
//    the two points are vertices, one gives the major radius, 
//    the other point gives the minor radius.
//
// For ARC_CURVE_TYPE
//    arc passes three points
//
// Author     : Jane Hu 
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             TBPoint const* point1_ptr,
                             TBPoint const* point2_ptr,
                             CubitVector const* intermediate_point_ptr) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  DLIList<CubitVector*> mid_points;
  if (intermediate_point_ptr)
  {
    CubitVector mid_point = *intermediate_point_ptr;
    mid_points.append(&mid_point);
  }

  CubitVector v1(point1_ptr->coordinates());
  CubitVector v2(point2_ptr->coordinates());

  gp_Pnt pt1(v1.x(),v1.y(), v1.z());
  gp_Pnt pt2(v2.x(),v2.y(), v2.z());

  CubitVector v3;
  gp_Pnt pt3;

  Handle(Geom_TrimmedCurve) curve_ptr;
  if(intermediate_point_ptr != NULL)
  {
    v3 = *intermediate_point_ptr;
    pt3.SetCoord(v3.x(),v3.y(), v3.z());
  }

  if (curve_type == STRAIGHT_CURVE_TYPE)
  {
     //make sure the two points are not coincident
     if(v1.about_equal(v2))
     {
        PRINT_ERROR("Can't create a line from two identical points.\n");
        return (Curve *)NULL;
     }
     curve_ptr = GC_MakeSegment(pt1,pt2);
  }

  else if (curve_type == ARC_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     curve_ptr = GC_MakeArcOfCircle(pt1, pt3, pt2);
  }

  else if (curve_type == ELLIPSE_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     
     gp_Pnt center(v3.x(), v3.y(), v3.z());

     gp_Elips ellipse;
     gce_MakeElips ellipse1(pt1	, pt2	, center);
     if(ellipse1.IsDone())
       ellipse = ellipse1.Value();
     else if(!ellipse1.IsDone() && ellipse1.Status() == gce_InvertAxis)
     {
        gce_MakeElips ellipse2(pt2, pt1, center);
        if(ellipse2.IsDone())
          ellipse = ellipse2.Value();
        else
        {
          PRINT_ERROR("Can't create an ellipse from give 3 points.\n");
          return (Curve *)NULL;
        }      
     } 
     else
     {
        PRINT_ERROR("Can't create an ellipse from give 3 points.\n");
        return (Curve *)NULL;
     }
     curve_ptr = GC_MakeArcOfEllipse(ellipse, pt1, pt2, CUBIT_TRUE);
  }

  else if(curve_type == PARABOLA_CURVE_TYPE || 
          curve_type == HYPERBOLA_CURVE_TYPE)
  {
    assert(intermediate_point_ptr != NULL);

    //find the directrix and focus of the parabola
    //or the axis, major radius and minor radius of the hyperbola
    CubitVector width_vec = v2 - v1;
    if(width_vec.length_squared() < TOL * TOL)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "2 end points are the same.\n");
       return (Curve *)NULL;
    }

    CubitVector midpoint_vec = (v1 + v2)/2.0;
    CubitVector height_vec = midpoint_vec - v3;
    gp_Pnt center(v3.x(), v3.y(), v3.z());
 
    if (height_vec.length_squared() < TOL * TOL)
    { 
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    CubitVector x = height_vec;
    x.normalize();
    gp_Dir x_dir(x.x(), x.y(), x.z());
 
    CubitVector N = x * (v2 - v1);  
    if (N.length_squared() < TOL * TOL)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    N.normalize();
    gp_Dir N_dir(N.x(), N.y(), N.z());

    gp_Ax2 axis(center, N_dir, x_dir);  

    if(curve_type == HYPERBOLA_CURVE_TYPE)
    { 
       //    (focus2) (v3) . (v2)
       //          .   .   . (midpoint = focus1)
       //                  . (v1)
       CubitVector focus2 = 2 * v3 - midpoint_vec;

       //according to the definition of hyperbola,
       //2 * a = length(v2 - focus2)-length(v2 - focus1)

       double major = (v2 - focus2).length()/2.0 - (v2 - midpoint_vec).length()/2.0;

       // if a = 1/2 length major axis, b = 1/2 length minor axis and
       // c = distance center to focus, then a*a + b*b = c*c

       double c_squared = (midpoint_vec - v3).length_squared();
       double minor = sqrt(c_squared  - major*major );
       gp_Hypr hypt(axis, major, minor);
       curve_ptr =
             GC_MakeArcOfHyperbola(hypt, pt1, pt2, CUBIT_TRUE);
    }

    else
    {
       // Find the focus of this parabola.
       // Since for a parabola with its peak at the origin, y = (1/(4*a))*x^2,
       // and since we have restricted this parabola to be symmetric (per the 
       // FastQ method, see the FastQ file getwt.f), we can use the following 
       // relationship to
       // determine "a", the distance the focus lies from the peak on the line
       // formed by the peak and the midpoint of the start and end points`
       double a = width_vec.length_squared()/(16. * height_vec.length()); 
       gp_Parab parab(axis, a);
       curve_ptr =
		GC_MakeArcOfParabola(parab, pt1, pt2, CUBIT_TRUE);
    } 
  }

  else
  {
      PRINT_ERROR("In OCCModifyEngine::make_Curve\n"
                  "       Invalid curve type.\n");
      return (Curve *)NULL;
  }

  OCCPoint* occ_pt1 = CAST_TO(const_cast<TBPoint*>(point1_ptr),OCCPoint);
  OCCPoint* occ_pt2 = CAST_TO(const_cast<TBPoint*>(point2_ptr),OCCPoint);
  TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
  TopoDS_Vertex * vt2 = occ_pt2->get_TopoDS_Vertex();
  TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt2);
  return OCCQueryEngine::instance()->populate_topology_bridge(new_edge, true);
}

Surface* OCCModifyEngine::make_Surface( Surface *surface_ptr,
    std::map< TopologyBridge*, TopologyBridge* > * /*old_tb_to_new_tb*/) const
{
  // Set extended_from argument to CUBIT_FALSE by default
  return make_Surface(surface_ptr, CUBIT_FALSE);
}

//===============================================================================
// Function   : make_Surface
//              This function creates a surface given an existing surface, copy.
// Member Type: PRIVATE
// Description: make a surface, OCC allows to create a stand along surface,
//              however the CGM has a design of making all surfaces in a (sheet)
//              body. This will add complexity on all free surface related 
//              calculation and modification, and adding potential bugs too. 
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( Surface * surface_ptr,
                                 CubitBoolean extended_from) const
{
  OCCSurface* occ_surface = CAST_TO(surface_ptr, OCCSurface);
  if (!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  //Start of the codes
  double UMax, VMax, UMin, VMin;
  occ_surface->get_param_range_U(UMin, UMax);
  occ_surface->get_param_range_V(VMin, VMax);

  TopoDS_Face *theFace = occ_surface->get_TopoDS_Face();

  if( !theFace)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  TopoDS_Face newFace;
  BRepAdaptor_Surface asurface(*theFace);

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  double const height = 2*(bounding_box.diagonal()).length();
  CubitBox box = occ_surface->bounding_box();
  double ratio = height/(box.diagonal().length());

  double middleU = (UMin + UMax)/2.0;
  double middleV = (VMin + VMax)/2.0;
  double U1 = middleU - (UMax-UMin)/2.0 * ratio;
  double U2 = middleU + (UMax-UMin)/2.0 * ratio;
  double V1 = middleV - (VMax - VMin)/2.0 * ratio;
  double V2 = middleV + (VMax - VMin)/2.0 * ratio;

  if (extended_from == CUBIT_TRUE)
  {
     // We need to get the type of surface.
     GeometryType type = occ_surface->geometry_type();
     if (type  == PLANE_SURFACE_TYPE)
     {
        gp_Pln plane = asurface.Plane();
        newFace = BRepBuilderAPI_MakeFace(plane, U1, U2, V1, V2);
     }
     else if(type == CONE_SURFACE_TYPE)
     {
       //make an infinite cone.
       //Given this lets create another face that is extended from it.
       if(asurface.GetType() == GeomAbs_Cone)
       {
         gp_Cone cone = asurface.Cone();
#if OCC_VERSION_MINOR > 5 
        newFace = BRepBuilderAPI_MakeFace(cone, U1, U2, V1, V2);
#else
         gp_Pnt Apex = cone.Apex();
         double semi_angle = cone.SemiAngle();
         gp_Pnt p2;
         double radius2;
         gp_XYZ xyz;
         if (semi_angle > 0)
           xyz = Apex.XYZ() + cone.Position().Direction().XYZ()*height;
         else
           xyz = Apex.XYZ() - cone.Position().Direction().XYZ()*height;

         p2.SetXYZ(xyz);
         radius2 = height * tan(fabs(semi_angle));
         Handle(Geom_RectangularTrimmedSurface) trimmed_cone;
         trimmed_cone = GC_MakeTrimmedCone(Apex, p2, 0, radius2); 
  #if OCC_VERSION_MAINTENANCE < 2 
         newFace = BRepBuilderAPI_MakeFace(trimmed_cone);
  #else
         newFace = BRepBuilderAPI_MakeFace(trimmed_cone, TOL);
  #endif
#endif
       }
       else
       {
         gp_Cylinder cylinder = asurface.Cylinder();
#if OCC_VERSION_MINOR > 5
         newFace = BRepBuilderAPI_MakeFace(cylinder, U1, U2, V1, V2);
#else
         double radius = cylinder.Radius();
         gp_Ax1 axis = cylinder.Axis(); 
         Handle(Geom_RectangularTrimmedSurface) trimmed_cyl;
         trimmed_cyl = GC_MakeTrimmedCylinder(axis, radius, height);
  #if OCC_VERSION_MAINTENANCE < 2 
         newFace = BRepBuilderAPI_MakeFace(trimmed_cyl);
  #else
         newFace = BRepBuilderAPI_MakeFace(trimmed_cyl, TOL);
  #endif
#endif
       } 
     }
     else if(type == SPHERE_SURFACE_TYPE)
     {
       //make a whole sphere.
       gp_Sphere sphere = asurface.Sphere();
       newFace = BRepBuilderAPI_MakeFace(sphere);
     }
     else if(type == TORUS_SURFACE_TYPE)
     {
       //make a whole torus
       gp_Torus torus = asurface.Torus();
       newFace = BRepBuilderAPI_MakeFace(torus);
     }
     else if(type == SPLINE_SURFACE_TYPE ) 
     {
       //extend the surfaces using the equation if possible.
       Handle(Geom_BezierSurface) bezier = asurface.Bezier();
       Handle(Geom_Surface) p_surf = bezier;
#if OCC_VERSION_MINOR > 5
       newFace = BRepBuilderAPI_MakeFace(p_surf, U1, U2, V1, V2, TOL);
#else
  #if  OCC_VERSION_MAINTENANCE < 2                           
       newFace = BRepBuilderAPI_MakeFace(p_surf, U1, U2, V1, V2);
  #else
       newFace = BRepBuilderAPI_MakeFace(p_surf, U1, U2, V1, V2, TOL);
  #endif
#endif
     }
  }
 
  else
  {
    BRepBuilderAPI_Copy api_copy(*theFace);
    TopoDS_Shape newShape = api_copy.ModifiedShape(*theFace);
    newFace = TopoDS::Face(newShape);
  }
  
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               newFace, CUBIT_TRUE);

  return surface;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface of type surface_type, given the list of curves.
//              check edges option is done in GeometryModifyTool level, so 
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( GeometryType surface_type,
                                 DLIList<Curve*>& curve_list,
                                 Surface * old_surface_ptr,
                                 bool check_edges) const
{
  //Create TopoDS_Edge list to make a surface.
  DLIList<DLIList<TopoDS_Edge*>*> topo_edges_loops;
  curve_list.reset() ;
    
  //check no intersections of the TopoDS_Edge's.
  //need to check that no intersection in the middle of the curves, not at
  //vertices or out of boundary.

  int count = 0; //intersection point should be the same as curve_list size.
  for ( int i = 0 ; i < curve_list.size()-1 ; i++ )
  {
     for(int j = i+1; j < curve_list.size(); j ++)
     {
        DLIList<CubitVector> intscts;
 	CubitBoolean bounded = CUBIT_TRUE;//dummy arg.
	CubitBoolean closest = CUBIT_TRUE;//dummy arg.
        CubitStatus yes_int = 
              OCCQueryEngine::instance()->get_intersections(curve_list[i],
				curve_list[j], intscts, bounded, closest);
        if(yes_int)
        {
           //check intscts point should be vertex or outside boundary.
 	   if (intscts.size() > 2 )  
	   {
	     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
             return (Surface *)NULL;
           }
           else
           {
             for(int k = 0; k < intscts.size(); k++)
             {
               CubitVector &v = intscts.get_and_step();
	       CubitPointContainment is_on = CAST_TO(curve_list[i],OCCCurve)->
					point_containment(v);
               if (is_on == CUBIT_PNT_BOUNDARY)
               {
	 	 is_on = CAST_TO(curve_list[j],OCCCurve)->
				point_containment(v);
		 if (is_on == CUBIT_PNT_BOUNDARY)
                   count++;
               }
               else if(is_on == CUBIT_PNT_INSIDE)
               {
                 PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
                 return (Surface *)NULL;
               }
	     }
	   }
        }
     }
  }
 
  if (count > curve_list.size()) 
  {
      PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                "       Cannot make Surface with intersecting curves.\n");
      return (Surface *)NULL;
  }

  CubitStatus stat = sort_curves(curve_list, topo_edges_loops); 
  if( stat == CUBIT_FAILURE ) //case of one disconnected curve , open wires
  {
     //loft curves.
     BRepOffsetAPI_ThruSections loft(CUBIT_FALSE);
     CubitStatus stat = do_loft(loft, topo_edges_loops);
     for (int i = 0; i < topo_edges_loops.size(); i++)
     {
       DLIList<TopoDS_Edge*>* topo_edges = topo_edges_loops.get_and_step();
       for(int j = 0; j < topo_edges->size(); j++)
         topo_edges->pop();
       delete topo_edges;
       topo_edges = NULL;
     }
     if(!stat)
       return (Surface*) NULL;

     TopoDS_Shape shape = loft.Shape();
     TopoDS_Shell shell = TopoDS::Shell(shape);
     TopExp_Explorer Ex;
     int num_surfaces = 0;
     TopoDS_Face topo_face ;
     for (Ex.Init(shell, TopAbs_FACE); Ex.More(); Ex.Next())
     {
       topo_face = TopoDS::Face(Ex.Current());
       num_surfaces ++;
     }

     if(num_surfaces != 1)
     {
       PRINT_ERROR("In OCCModifyEngine::skin_surface\n"
                 "   Cannot create a skin surface for given curves.\n");
       return (Surface*) NULL;
     }
   
     Surface* surf = OCCQueryEngine::instance()->populate_topology_bridge(topo_face, CUBIT_TRUE);
     if (surf == NULL)
     {
       PRINT_ERROR("In OCCModifyEngine::skin_surfaces\n"
                   "   Cannot create a skin surface for given curves.\n");
       return (Surface*) NULL;
     }
 
     return surf;
  }

  // Use the topo_edges to make a topo_face
  TopoDS_Face* topo_face;
  topo_face = make_TopoDS_Face(surface_type,topo_edges_loops, old_surface_ptr);
 
  for (int i = 0; i < topo_edges_loops.size(); i++)
  {
    DLIList<TopoDS_Edge*>* topo_edges = topo_edges_loops.get_and_step();
    for(int j = 0; j < topo_edges->size(); j++)
      topo_edges->pop();
    delete topo_edges;
    topo_edges = NULL;
  }
  
  if(!topo_face)
  {
     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface object.\n");
     return (Surface *)NULL;
  }

  // make the topology bridges for the face
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               *topo_face, CUBIT_TRUE); 
  topo_face->Nullify();
  delete topo_face;
  topo_face = NULL;

  //Created new surface uses existing OCCPoints on the curves, but created
  //new curves, so remove the curves from the curvelist on those points.
  if(surface)
  {
    for(int i = 0; i <  curve_list.size(); i++)
    {
      OCCCurve* test_c = CAST_TO(curve_list.get_and_step(), OCCCurve);
      DLIList<OCCPoint*> points;
      CAST_TO(surface, OCCSurface)->get_points(points);
      for(int j = 0; j <  points.size(); j ++)
        points.get_and_step()->remove_curve(test_c);      
    }
  }
  return surface ;
}

CubitStatus OCCModifyEngine::do_loft(BRepOffsetAPI_ThruSections& loft,
                                     DLIList<DLIList<TopoDS_Edge*>*> loops) const
{
   TopoDS_Edge  new_edge;
   for(int i = 0; i < loops.size(); i++)
   {
     BRepBuilderAPI_MakeWire aWire;
     DLIList<TopoDS_Edge*> edges = *(loops.get_and_step());

     for(int j = 0; j <  edges.size(); j++)
     {
       TopoDS_Edge* topo_edge = edges.get_and_step();

       BRepBuilderAPI_Copy api_copy(*topo_edge);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_edge);
       new_edge = TopoDS::Edge(newShape);
       aWire.Add(new_edge);
     }
     loft.AddWire(aWire.Wire());
   }
   loft.Build();
   if(!loft.IsDone())
   {
     PRINT_ERROR("Curves can't be loft into a surface.\n");
     return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
} 

//===============================================================================
// Function   : sort_curves
// Member Type: PROTECTED
// Description: sort the curves so they are in order and make closed loop 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::sort_curves(DLIList<Curve*> curve_list,
                        DLIList<DLIList<TopoDS_Edge*>*>& topo_edges_loops)const
{
  topo_edges_loops.clean_out();
  CubitStatus stat = CUBIT_SUCCESS;
  std::vector< DLIList<TopoDS_Edge*>* > topo_edges(curve_list.size());
  int size_in = curve_list.size();
  for(int i = 0; i < size_in; i++)
    topo_edges[i] = new DLIList<TopoDS_Edge*>;

  curve_list.reset() ;
  Curve const* curve_ptr = NULL ;
  OCCCurve* occ_curve = NULL;
  TopoDS_Edge* topo_edge = NULL;

  OCCPoint* start = NULL;
  OCCPoint* end = NULL;
  DLIList<OCCPoint*> point_list;
  CubitBoolean new_end = CUBIT_FALSE;
  int size = curve_list.size();

  int count = 0;
  for ( int i = 0 ; i < size ; i++ )
  {
     if (i == 0)
       new_end = CUBIT_TRUE;
     for(int j = 0; j < curve_list.size(); j ++)
     {
        curve_ptr = curve_list.get() ; 
        occ_curve = CAST_TO(const_cast<Curve*>(curve_ptr), OCCCurve);

        if(occ_curve ==  NULL)
        {
           PRINT_ERROR("In OCCModifyEngine::sort_curves\n"
                       "       Got a NULL pointer to OCCCurve\n") ;
           return CUBIT_FAILURE;
        }

        point_list.clean_out();
        occ_curve->get_points(point_list);
        //assert(point_list.size()==2);

        if (i == 0)
        {
          start = point_list.get();
          end = point_list.pop();  
          break;
        }

        if(end->is_equal(*(point_list.get()), TOL) ||
           end->is_equal(*(point_list.step_and_get()),TOL)) 
        {
           end = point_list.step_and_get();
           new_end = CUBIT_TRUE;
           break;
        }
   
        else if(start->is_equal(*(point_list.get()), TOL) ||
           start->is_equal(*(point_list.step_and_get()),TOL))
        {
           start = end;
           end = point_list.step_and_get(); 
           new_end = CUBIT_TRUE;
           break;
        }
        curve_list.step();
     }

     if (new_end)//found next curve 
     {
        topo_edge = occ_curve->get_TopoDS_Edge();
        topo_edges[count]->append(topo_edge);
        curve_list.remove();
        if(start->is_equal( *end, TOL))  //formed a closed loop
        {
          i = -1;
          size = curve_list.size() ;
          topo_edges_loops.append(topo_edges[count]);
          count++;
        }
        else
          new_end = CUBIT_FALSE;
     }
     else
     {
        stat = CUBIT_FAILURE; 
        i = -1;
        size = curve_list.size();
        topo_edges_loops.append(topo_edges[count]);
        count++;
     }
  }

  if( new_end == CUBIT_FALSE ) //case of one disconnected curve
  {
    topo_edges_loops.append(topo_edges[count]); 
    stat = CUBIT_FAILURE;
  }

  for(int i = 0; i < size_in; i++)
  {
     if(topo_edges[i]->size() == 0)
     {
       delete topo_edges[i];
       topo_edges[i] = NULL;
     }
  }
  return stat;
} 
//===============================================================================
// Function   : make_TopoDS_Face
// Member Type: PROTECTED
// Description: make a opoDS_Face of type surface_type, given the list of 
//              TopoDS_Edge. the TopoDS_Edge's should be in order in loops.
//              check edges option is done in GeometryModifyTool level, so
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
TopoDS_Face* OCCModifyEngine::make_TopoDS_Face(GeometryType surface_type,
			      DLIList<DLIList<TopoDS_Edge*>*> topo_edges_list,
			      Surface * old_surface_ptr)const
{
  TopoDS_Face* topo_face = NULL;
  // Make sure a supported type of surface is being requested.
  if ( surface_type != PLANE_SURFACE_TYPE  &&
       surface_type != BEST_FIT_SURFACE_TYPE)
  {
      PRINT_WARNING("In OCCGeometryEngine::make_TopoDS_Face\n"
                    "       At this time, cannot make a TopoDS_Face that isn't"
                    " planar or best fit.\n");
      return topo_face;
  }
 
  // Set the TopoDS_Face pointer, if requested.
  TopoDS_Face *fit_Face = NULL;
  Handle_Geom_Surface S;
  if ( old_surface_ptr != NULL )
  {
      OCCSurface *surf = CAST_TO(old_surface_ptr, OCCSurface );
      fit_Face = surf->get_TopoDS_Face();
      S = BRep_Tool::Surface(*fit_Face);
  }
 
  // Make a wire from the topo_edges.
  // Coincident TopoDS_Vertex will be deleted by OCC.
  if(topo_edges_list.size() == 0)
      return topo_face;

  DLIList<TopoDS_Wire*> wires;
  GProp_GProps myProps;
  double max_area  = 0.0;
  TopoDS_Wire* out_Wire = NULL;
  TopoDS_Wire test_Wire;
  DLIList<TopoDS_Edge*>* topo_edges; 
  //check and make sure the outer loop is in the first
  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    topo_edges = topo_edges_list.get_and_step();
    BRepBuilderAPI_MakeWire aWire(*(topo_edges->get()));
    for(int j = 1; j < topo_edges->size(); j++)
      aWire.Add(*(topo_edges->step_and_get()));

    test_Wire = aWire.Wire();
    wires.append(&test_Wire);
   
    if (topo_edges_list.size() == 1)
      break;

    BRepBuilderAPI_MakeFace made_face(test_Wire);

    if (!made_face.IsDone())
    {
       PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                   "   Cannot find the best fit surface for given curves.\n");
       return topo_face;
    }
    TopoDS_Face test_face = made_face.Face();
    BRepGProp::SurfaceProperties(test_face, myProps); 
    double area = myProps.Mass();
    out_Wire = max_area > area ? out_Wire : &test_Wire;
    max_area = max_area > area ? max_area : area;
  } 

  if (out_Wire)
  {
    wires.remove(out_Wire);
    wires.insert_first(out_Wire);
  }

  //create the TopoDS_Face
  CubitBoolean error = CUBIT_FALSE;

  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    TopoDS_Wire *the_wire = wires.get_and_step();
    if (i == 0)
    {
      if( old_surface_ptr != NULL )
      {
        BRepBuilderAPI_MakeFace made_face(S, *the_wire);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }
        topo_face = new TopoDS_Face(made_face.Face());
      }
      else
      {
        CubitBoolean is_planar = (surface_type == PLANE_SURFACE_TYPE ?
				  CUBIT_TRUE : CUBIT_FALSE); 
        BRepBuilderAPI_MakeFace made_face(*the_wire, is_planar);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }

        topo_face = new TopoDS_Face(made_face.Face());
      }
    }
    else
    {
      BRepBuilderAPI_MakeFace made_face(*topo_face, *the_wire);
      if (!made_face.IsDone())
      {
        error = CUBIT_TRUE;
        break;
      }
      delete topo_face;
      topo_face = new TopoDS_Face(made_face.Face());
    }
  } 

  if(error)
  {
    PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                 "   Cannot find the best fit surface for given curves.\n");
    return (TopoDS_Face*) NULL;
  }

  return topo_face;
}
//===============================================================================
// Function   : make_Lump
// Member Type: PUBLIC
// Description: make a lump of one shell
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Lump* OCCModifyEngine::make_Lump( DLIList<Surface*>& surface_list ) const
{
  if (surface_list.size() < 2) 
    return (Lump*) NULL;

  //all surfaces should be stand along surface bodies or shell bodies' surface
  DLIList<BodySM*> body_list;
  for(int i = 0; i < surface_list.size(); i++)
  {
    OCCSurface* occ_surface = CAST_TO(surface_list.get_and_step(), OCCSurface);
    if (occ_surface == NULL)
    {
       PRINT_ERROR("Cannot create an OCC lump from the given surfaces.\n"
                 "Possible incompatible geometry engines.\n");
       return (Lump *)NULL;
    }
    OCCBody* occ_body = occ_surface->my_body();
    if(!occ_body)
    {
      OCCShell* occ_shell = occ_surface->my_shell();
      if(occ_shell)
        occ_body = occ_shell->my_body();
    }
    if(!occ_body)
    {
      DLIList<OCCBody*> original_bodies;
      occ_surface->get_bodies(original_bodies);
      if(original_bodies.size() > 1)
      {
        PRINT_ERROR( "Cannot make lump in non-mainfold solids. \n");
        return (Lump*) NULL;
      }
      else if(original_bodies.size() == 0)
      {
        PRINT_ERROR( "Interal error: Can't find associated solid. \n");
        return (Lump*) NULL;
      }
      occ_body = original_bodies.get();
      assert(occ_body != NULL);
    }
    DLIList<Lump*> lumps;
    DLIList<OCCShell*> shells;
    DLIList<OCCSurface*> surfaces;
    surfaces = occ_body->my_sheet_surfaces();
    shells = occ_body->shells();
    lumps = occ_body->lumps();
    if(lumps.size() > 0 || shells.size() + surfaces.size() > 1)
    {
      PRINT_ERROR("Cannot create an OCC lump from the given surfaces.\n"
               "The surfaces are not free.\n");
      return (Lump *)NULL;
    }
    body_list.append_unique(occ_body);
  }

  TopoDS_Shape aShape;
  CubitStatus stat = stitch_surfs(body_list, aShape);
  if(!stat)
  {
    PRINT_ERROR("The surfaces are not all connected, can't make a lump. \n");
    return (Lump*)NULL;
  }

  TopExp_Explorer Ex, Ex2;
  TopoDS_Shell aShell ;
  for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); Ex.More()&& stat; Ex.Next())
    aShell = TopoDS::Shell(Ex.Current());
 
  //check to make sure the aShell is closed.
  int num_edges = 0;
  int pairs = 0;
  //sometimes there's duplicate TopoDS_Edges in the shell.
  DLIList<TopoDS_Edge*> edge_list;
  for (Ex.Init(aShell, TopAbs_EDGE); Ex.More()&& stat; Ex.Next())
  {
    TopoDS_Edge edge1 = TopoDS::Edge(Ex.Current());
    TopoDS_Edge* new_edge = new TopoDS_Edge(edge1);
    edge_list.append(new_edge);
  }

  int size = edge_list.size();
  for (int i = 0; i < size; i++)
  {
    TopoDS_Edge edge1 = *edge_list[i];
    int same = 0;
    for (int j = i+1; j < edge_list.size(); j++)
    {
      TopoDS_Edge edge2 = *edge_list[j];
      if(edge1.IsEqual(edge2))
      {
           same ++;
           edge_list.remove(&edge1);
           i--;
           size--;
           break;
      }
    }
    if(same > 0)
      continue;

    else
      num_edges++;
  
    for (int j = 0; j < size; j++)  
    {
      TopoDS_Edge edge2 = *edge_list[j];    
      if (!edge1.IsEqual(edge2)&& edge1.IsSame(edge2))
      {
        pairs++;
        break;
      }
    }
  }

  for (int k = 0; k < edge_list.size(); k++)
  {
    TopoDS_Edge* edge = edge_list.get_and_step();
    edge->Nullify();
    delete edge;
    edge = NULL;
  }

  if (num_edges == pairs )
    aShell.Closed(CUBIT_TRUE);

  else
    PRINT_ERROR("Surfaces must make a water-tight shape to make a lump.\n");
  
  if(aShell.Closed())
  {
    BRepBuilderAPI_MakeSolid aMakeSolid(aShell);
    if (!aMakeSolid.IsDone())
    {
       PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                   "OCC internal error.\n");
       return (Lump *)NULL;
    }

    TopoDS_Solid aSolid = aMakeSolid.Solid();

    return
      OCCQueryEngine::instance()->populate_topology_bridge(aSolid, CUBIT_TRUE); 
  }

  return (Lump*) NULL;
}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM from a surface
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( Surface *surface ) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  if(!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC body from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  OCCBody* occ_body = occ_surface->my_body();
  TopoDS_Face* face = occ_surface->get_TopoDS_Face();
  TopoDS_Face newFace;
  if(!occ_body)
  {
    DLIList<OCCBody*> original_bodies;
    occ_surface->get_bodies(original_bodies);
    if(original_bodies.size() > 0)
      occ_body = original_bodies.get();
  }
  if(occ_body)
  {
     //copy the surface to make a sheet body.
     BRepBuilderAPI_Copy api_copy(*face);
     TopoDS_Shape newShape = api_copy.ModifiedShape(*face);
     newFace = TopoDS::Face(newShape);
     face = new TopoDS_Face(newFace);
  }

  surface = OCCQueryEngine::instance()->populate_topology_bridge(*face, CUBIT_TRUE);
   
  return CAST_TO(surface, OCCSurface)->my_body();
}



//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM given a list of Lumps.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( DLIList<Lump*>& lump_list ) const
{
  if (lump_list.size() == 0)
    return (BodySM*) NULL;
/*
  //Create a compsolid shape, copy all BodySM's solids to create new compbody 
  DLIList<BodySM*> bodysm_list;
  TopoDS_Compound CS;
  BRep_Builder B;
  B.MakeCompound(CS);

  //Add every shape to the CompSolid
  for(int i = 0; i < lump_list.size(); i++)
  {
     Lump* lump = lump_list.get_and_step();
     OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     if(!occ_lump)
     {
        PRINT_ERROR("Cannot create an OCC BodySM from the given lumps.\n"
                    "Possible incompatible geometry engines.\n");
        return (BodySM *)NULL;
     }
     TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();
     BRepBuilderAPI_Copy api_copy(*solid);
     TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
     B.Add(CS, newShape);
  }
*/
 
  //check if the lumps are already in some bodies, in this case, have to use 
  //unite operation to create compound.
  for (int i = 0 ; i < lump_list.size(); i++)
  {
    OCCLump* lump = CAST_TO(lump_list.get_and_step(), OCCLump);
    if(lump == NULL)
    {
      PRINT_ERROR("Incompatible engines.\n");
      return (BodySM*) NULL;
    }
    BodySM* body = lump->get_body();
    if(body != NULL)
      OCCQueryEngine::instance()->delete_body(body, CUBIT_FALSE);
  }
  TopoDS_Compound* Co;
  DLIList<OCCShell*> shells;
  DLIList<OCCSurface*> surfaces;
  Co = OCCBody::make_Compound(lump_list, shells, surfaces); 
  assert (Co != NULL);
  BodySM* bodysm = OCCQueryEngine::instance()->populate_topology_bridge(*Co);
  Co->Nullify();
  delete Co;
  Co = NULL;
  return bodysm;

}


//===============================================================================
// Function   : sphere
// Member Type: PUBLIC
// Description: build an OCC sphere
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::sphere(double radius) const
{
  if (radius <= 0)
    return (BodySM*) NULL; 

  TopoDS_Solid S = BRepPrimAPI_MakeSphere(radius);
  
  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S, 
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->get_body();
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: build an OCC brick 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( double wid, double dep, double hi ) const
{
  if (wid <= 0 || dep <=0 || hi <= 0)
    return (BodySM*)NULL;
  
  TopoDS_Solid S = BRepPrimAPI_MakeBox(wid, dep, hi);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);

  if (lump == NULL)
    return (BodySM*)NULL;

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  if(body)
    CAST_TO(body,OCCBody)->move(-wid/2.0, -dep/2.0, -hi/2.0);
  return body;
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: create an OCC brick given center axes and extension
//              extension is equvlent to (wid/2, dep/2, hi/2)
//              center should be given the coordinates in global system
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( const CubitVector& center, 
                                const CubitVector axes[3],
                                const CubitVector &extension) const
{
  gp_Pnt left_point(0,0,0);
  gp_Dir main_dir(axes[2].x(), axes[2].y(), axes[2].z());
  gp_Dir x_dir(axes[0].x(), axes[0].y(), axes[0].z());
  gp_Ax2 Axis(left_point, main_dir, x_dir);
  TopoDS_Solid S = BRepPrimAPI_MakeBox( Axis, extension.x()*2, extension.y()*2,
					extension.z()*2);

  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  if(body)
  {
    CubitVector center_point;
    double volume;
    OCCBody* occ_body = CAST_TO(body,OCCBody);
    occ_body->mass_properties(center_point, volume); 
    CubitVector move_vec = center - center_point;
    occ_body->move(move_vec.x(), move_vec.y(), move_vec.z());
    return body;
  }
  return (BodySM*)NULL;
}

//===============================================================================
// Function   : prism
// Member Type: PUBLIC
// Description: create an OCC prism 
// Author     : Jane Hu  
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::prism( double height, int sides, double major,
                               double minor) const
{
  if(major <= 0. || minor <= 0. || (major - minor) <=  -TOL)
  {
    PRINT_ERROR("Major and minor radii must be greater than zero.\n");
    return (BodySM*)NULL;
  }

  if (sides == 4)
    return brick(2 * major, 2 * minor, height); 

  TopoDS_Wire wire ;
    make_base_for_prim_pyramid(major, minor, height, sides, wire);

  TopoDS_Face base = BRepBuilderAPI_MakeFace(wire, Standard_True);
  gp_Dir main_dir(0.0, 0.0, 1.0);
  gp_Vec norm(main_dir);
  norm *= height;
  BRepSweep_Prism swept(base, norm);
  TopoDS_Shape new_shape = swept.Shape();
  DLIList<TopologyBridge*> tbs;
  tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
  assert(tbs.size() == 1);

  BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
  return bodysm;
}

void OCCModifyEngine::make_base_for_prim_pyramid(double major,
                                                 double minor,
                                                 double height,
                                                 int sides,
                                                 TopoDS_Wire& wire)const
{
  //One of the polygon side will be perpendicular to positive x-axis.
  double y = major * sin(CUBIT_PI/sides);
  double x = sqrt(major * major - y * y);
  gp_Pnt start(x, y, -height/2.0);

  DLIList<gp_Pnt> point_list;
  double theta = 2.0/sides*CUBIT_PI;
  for(int n =1 ; n < sides; n++)
  {
    double angle = theta * (n + 0.5);
    gp_Pnt v(major * cos(angle), major * sin(angle), -height/2.0);
    point_list.append(v);
  }

  TopoDS_Edge new_edge;
  BRepBuilderAPI_MakePolygon poly_maker;
  gp_Dir main_dir(0.0, 0.0, 1.0);
  point_list.append(start);

  if (fabs(major - minor) < TOL)
    for (int i = 0; i <sides; i++)
      poly_maker.Add(point_list[i]);
  else
  {
    for (int i = 0; i <sides; i++)
    {
      x = point_list[i].X();
      if((y=point_list[i].Y()) > 0.0)
        y = sqrt((1-x*x/major/major)*minor*minor);
      else
        y = -sqrt((1-x*x/major/major)*minor*minor);
      point_list[i].SetY(y);
      poly_maker.Add(point_list[i]);
    }
  }
  poly_maker.Close();

  wire = poly_maker.Wire();
}

//===============================================================================
// Function   : pyramid
// Member Type: PUBLIC
// Description: create an OCC pyramid 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::pyramid( double height, int sides, double major,
                                 double minor, double top) const
{
  TopoDS_Solid S;

  //build the top and bottom shapes.
  TopoDS_Wire wire_bottom ; 
  make_base_for_prim_pyramid(major, minor, height, sides, wire_bottom);
  double top_minor = top * minor / major;
  TopoDS_Wire wire_top ; 
  BRepOffsetAPI_ThruSections builder(CUBIT_TRUE, CUBIT_TRUE);
  builder.AddWire(wire_bottom);
  if(top > TOL)
  {
    make_base_for_prim_pyramid(top, top_minor, -height, sides, wire_top); 
    builder.AddWire(wire_top);
  }
  else
  {
    gp_Pnt pt = gp_Pnt( 0.0, 0.0, height/2.0);
    TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
    builder.AddVertex(theVertex);
  }
  builder.Build() ;
  S = TopoDS::Solid(builder.Shape());
 
  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  return body;
  
}

//===============================================================================
// Function   : cylinder
// Member Type: PUBLIC
// Description: create an OCC cylinder, its base shape can be ellipse with
//		r1, r2 while r3 = 0; or it can be a cone with r1, 
//		r3 while r2 = 0. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::cylinder( double hi, double r1, double r2, double r3 ) const
{
  TopoDS_Solid S;
  if(r2 != r1)//elliptical based cylinder
  {
    gp_Pnt center(0.0, 0.0, 0.0);
    gp_Dir main_dir(0.0, 0.0, 1.0);
    gp_Dir x_dir(1.0, 0.0, 0.0);
    gp_Ax2 Axis(center, main_dir, x_dir); 
    if(r1 < r2)
    {  
      double temp_r = r1;
      r1 = r2;
      r2 = temp_r;
    }
    Handle(Geom_Curve) curve_ptr = GC_MakeEllipse(Axis, r1, r2); 
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
    BRepBuilderAPI_MakeWire aWire(new_edge);

    TopoDS_Wire test_Wire = aWire.Wire();
    
    BRepOffsetAPI_ThruSections builder(CUBIT_TRUE, CUBIT_TRUE);
    builder.AddWire(test_Wire);
    if (r3 == 0)
    {
      gp_Pnt pt = gp_Pnt( 0.0, 0.0, hi);
      TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
      builder.AddVertex(theVertex);
    }
    else
    {
      gp_Pnt center2(0.0, 0.0,hi);
      gp_Ax2 Axis2(center2, main_dir, x_dir);
      Handle(Geom_Curve) curve_ptr = GC_MakeEllipse(Axis2, r3, r3*r2/r1);
      TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
      BRepBuilderAPI_MakeWire aWire(new_edge);
      TopoDS_Wire test_Wire = aWire.Wire();
      builder.AddWire(test_Wire);
    }
    builder.Build() ;
    S = TopoDS::Solid(builder.Shape());
  }

  else // cone
  {
    if(r1 == r3) //cylinder
      S = BRepPrimAPI_MakeCylinder(r1, hi);
    else
      S = BRepPrimAPI_MakeCone(r1, r3, hi);
  }

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::cylinder\n"
                "   Cannot create a cylinder for given radii.\n");
    return (BodySM*)NULL;
  }

  BodySM* body = CAST_TO(lump, OCCLump)->get_body();
  if(body)
    CAST_TO(body, OCCBody)->move(0, 0, -hi/2.0);
  return body;
}

//===============================================================================
// Function   : torus
// Member Type: PUBLIC
// Description: create an OCC torus
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::torus( double r1, double r2 ) const
{
  if (r1 <= 0 || r2 <= 0)
    return (BodySM*) NULL;
 
  TopoDS_Solid S = BRepPrimAPI_MakeTorus(r1, r2);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::torus\n"
                "   Cannot create a torus for given radii.\n");
    return (BodySM*)NULL;
  }

  return CAST_TO(lump, OCCLump)->get_body();
}

//===============================================================================
// Function   : planar_sheet
// Member Type: PUBLIC
// Description: create an OCC planar_sheet with four vectors. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::planar_sheet ( const CubitVector& p1,
                                       const CubitVector& p2,
                                       const CubitVector& p3,
                                       const CubitVector& p4) const
{
  TBPoint* point1 = make_Point(p1);
  TBPoint* point2 = make_Point(p2);
  TBPoint* point3 = make_Point(p3);
  TBPoint* point4 = make_Point(p4);
  Curve * curve1 = make_Curve( point1, point2);
  if (curve1 == NULL)
	return (BodySM*) NULL;
  Curve * curve2 = make_Curve( point2, point3); 
  if (curve2 == NULL)
        return (BodySM*) NULL;
  Curve * curve3 = make_Curve( point3, point4);
  if (curve3 == NULL)
        return (BodySM*) NULL;
  Curve * curve4 = make_Curve( point4, point1);
  if (curve4 == NULL)
        return (BodySM*) NULL;
  DLIList<Curve*> curves;
  curves.append(curve1);
  curves.append(curve2);
  curves.append(curve3);
  curves.append(curve4);
  Surface* surface = make_Surface(PLANE_SURFACE_TYPE, curves);
  if (surface == NULL)
	return (BodySM*) NULL;

  return CAST_TO(surface,OCCSurface)->my_body();
}

//===============================================================================
// Function   : copy_body
// Member Type: PUBLIC
// Description: copy an OCC-based body
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::copy_body ( BodySM* bodyPtr, std::map<TopologyBridge*, TopologyBridge*> * /*old_tb_to_new_tb*/ ) const
{
  BodySM* new_body = NULL;
  OCCBody* occ_body = CAST_TO(bodyPtr, OCCBody);
  if (!occ_body)
  {
     PRINT_ERROR("Cannot create an OCC bodySM from the given bodySM.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  DLIList<CubitSimpleAttrib> list;
  TopoDS_Compound *theCS = occ_body->get_TopoDS_Shape();
  
  if (theCS == NULL) //single lump or shell or surface body
  {
    DLIList<OCCShell*> shells = occ_body->shells();
    assert(shells.size() < 2);
    for(int i = 0 ; i < shells.size(); i++)
    {
      TopoDS_Shell* shell = shells.get_and_step()->get_TopoDS_Shell();
      BRepBuilderAPI_Copy api_copy(*shell);
      TopoDS_Shape newShape = api_copy.ModifiedShape(*shell);
      TopoDS_Shell newShell = TopoDS::Shell(newShape);
      new_body = OCCQueryEngine::instance()->populate_topology_bridge(newShell, CUBIT_TRUE)->my_body();
      copy_body_attributes((TopoDS_Shape)(*shell), api_copy);
    }
 
    DLIList<OCCSurface*> surfaces = occ_body->my_sheet_surfaces();
    assert(surfaces.size() < 2);
    
    for(int i = 0 ; !new_body && i < surfaces.size(); i++)
    {
      OCCSurface *occ_surface = CAST_TO(surfaces.get_and_step(), OCCSurface);
      TopoDS_Face *theFace = occ_surface->get_TopoDS_Face();
      BRepBuilderAPI_Copy api_copy(*theFace);
      TopoDS_Shape newShape = api_copy.ModifiedShape(*theFace);
      TopoDS_Face newFace = TopoDS::Face(newShape);
      Surface* surface = OCCQueryEngine::instance()->populate_topology_bridge(
                           newFace, CUBIT_TRUE );
      OCCSurface* occ_surf = CAST_TO(surface, OCCSurface);
      new_body = occ_surf->my_body();
      copy_body_attributes((TopoDS_Shape)(*theFace) , api_copy);
    }

    //single lump body
    if (!new_body)
    {
      Lump *lump = occ_body->lumps().get();
      TopoDS_Solid solid = *(CAST_TO(lump, OCCLump)->get_TopoDS_Solid());
      BRepBuilderAPI_Copy api_copy(solid);
      TopoDS_Shape newShape = api_copy.ModifiedShape(solid);
      TopoDS_Solid newSolid = TopoDS::Solid(newShape);
      lump = OCCQueryEngine::instance()->populate_topology_bridge(newSolid,
                                                        CUBIT_TRUE);

      new_body = CAST_TO(lump, OCCLump)->get_body();
      copy_body_attributes((TopoDS_Shape)solid, api_copy);
    }
  }

  if(!new_body && theCS && !theCS->IsNull() &&
     OCCQueryEngine::instance()->OCCMap->IsBound(*theCS))
  {
    BRepBuilderAPI_Copy api_copy(*theCS);

    TopoDS_Shape newShape = api_copy.ModifiedShape(*theCS);

    TopoDS_Compound newCS = TopoDS::Compound(newShape);

    new_body = OCCQueryEngine::instance()->populate_topology_bridge(newCS);
    copy_body_attributes((TopoDS_Shape)(*theCS), api_copy);
    OCCAttribSet::get_attributes(newCS, list);
    for(int kk = 0; kk < list.size(); kk++)
      new_body->append_simple_attribute_virt(list.get_and_step());
  }
  
  return new_body;
}

BodySM* OCCModifyEngine::create_body( VolumeFacets& volume, 
                           std::map<FacetShapes*, GeometryEntity*>& entity_map,
                           const FacetPointSet& points, 
                           int interp_order) const
{
  return (BodySM*) NULL;
}

CubitStatus OCCModifyEngine::copy_body_attributes(TopoDS_Shape orig_shape,
                                           BRepBuilderAPI_Copy& api_copy)const 
{
  DLIList<CubitSimpleAttrib> list;
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(orig_shape, TopAbs_SOLID, M);
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Solid solid = TopoDS::Solid(M(ii));
    TopoDS_Solid new_solid = TopoDS::Solid(api_copy.ModifiedShape(solid));
    if(!new_solid.IsNull())
      OCCQueryEngine::instance()->copy_attributes(solid, new_solid);
    list.clean_out();
    OCCAttribSet::get_attributes(new_solid,list);
    OCCLump *lump = NULL;
    int k = OCCQueryEngine::instance()->OCCMap->Find(new_solid);
    lump = (OCCLump*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;         
    for(int kk = 0; kk < list.size(); kk++)
      lump->append_simple_attribute_virt(list.get_and_step());
    
    if (list.size() == 0)
    {
      k = OCCQueryEngine::instance()->OCCMap->Find(solid);
      OCCLump *orig_lump = (OCCLump*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;
      OCCBody* body = (OCCBody*)orig_lump->get_body();
      if(body)
      {
        body->get_simple_attribute(list);
        for (int kk = 0; kk < list.size(); kk++) 
          lump->get_body()->append_simple_attribute_virt(list.get_and_step()); 
      }
    }
  }

  M.Clear();
  TopExp::MapShapes(orig_shape, TopAbs_FACE, M);
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Face face = TopoDS::Face(M(ii));
    TopoDS_Face new_face = TopoDS::Face(api_copy.ModifiedShape(face));
    if(!new_face.IsNull())
      OCCQueryEngine::instance()->copy_attributes(face, new_face);
    list.clean_out();
    OCCAttribSet::get_attributes(new_face, list);
    OCCSurface* surf = NULL;
    int k = OCCQueryEngine::instance()->OCCMap->Find(new_face);
    surf = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    for(int kk = 0; kk < list.size(); kk++)
      surf->append_simple_attribute_virt(list.get_and_step());
  
    if(list.size() == 0)
    {
      k = OCCQueryEngine::instance()->OCCMap->Find(orig_shape);
      OCCBody* body = NULL;
      if(orig_shape.ShapeType() == TopAbs_FACE)
      {
        OCCSurface *orig_surf = (OCCSurface*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;
        body = orig_surf->my_body();
      }
      else if(orig_shape.ShapeType() == TopAbs_SHELL)
      {
        OCCShell* orig_shell = (OCCShell*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;
        body = orig_shell->my_body();
      }
      //Solid and Compound case has been considered in the above cases.
      if(body)
      {
        body->get_simple_attribute(list);
        for (int kk = 0; kk < list.size(); kk++) 
          surf->my_body()->append_simple_attribute_virt(list.get_and_step());
      }
    }
  }

  M.Clear();
  TopExp::MapShapes(orig_shape, TopAbs_EDGE, M);
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Edge edge = TopoDS::Edge(M(ii));
    TopoDS_Edge new_edge = TopoDS::Edge(api_copy.ModifiedShape(edge));
    if(!new_edge.IsNull())
      OCCQueryEngine::instance()->copy_attributes(edge, new_edge);
    list.clean_out();
    OCCAttribSet::get_attributes(new_edge, list);
    if(list.size() > 0)
    {
      int k = OCCQueryEngine::instance()->OCCMap->Find(new_edge);
      OCCCurve* curve = (OCCCurve*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;
      for(int kk = 0; kk < list.size(); kk++)
        curve->append_simple_attribute_virt(list.get_and_step());
    }
  }

  M.Clear();
  TopExp::MapShapes(orig_shape, TopAbs_VERTEX, M);
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Vertex vertex = TopoDS::Vertex(M(ii));
    TopoDS_Vertex new_vertex = TopoDS::Vertex(api_copy.ModifiedShape(vertex));
    if(!new_vertex.IsNull())
      OCCQueryEngine::instance()->copy_attributes(vertex, new_vertex);
    list.clean_out();
    OCCAttribSet::get_attributes(new_vertex, list);
    if(list.size() > 0)
    {
      int k = OCCQueryEngine::instance()->OCCMap->Find(new_vertex);
      OCCPoint* point = (OCCPoint*) (OCCQueryEngine::instance()->OccToCGM->find(k))->second;
      for(int kk = 0; kk < list.size(); kk++)
        point->append_simple_attribute_virt(list.get_and_step());
    }
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : stitch
// Member Type: PUBLIC
// Description: stitch all surfs and try to make a shell body.
//              tighten_gaps and tolerance are not used here.
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::stitch(
                      DLIList<BodySM*>& surf_bodies,
                      DLIList<BodySM*> &new_bodies,
                      bool tighten_gaps,
                      double tolerance) const
{
  if (surf_bodies.size()==0)
    return CUBIT_SUCCESS;

  if (surf_bodies.size()==1)
  {
    new_bodies = surf_bodies;
    return CUBIT_SUCCESS;
  }

  TopoDS_Shape fuse;
  CubitStatus stat =  stitch_surfs(surf_bodies, fuse);
    
  TopExp_Explorer Ex;
  DLIList<TopologyBridge*> tbs;
  for (Ex.Init(fuse, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
  {
    TopoDS_Shape shape = Ex.Current();
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(shape);
  }

  if (stat)
  {
    BodySM* body = CAST_TO(tbs.get(), BodySM);
    new_bodies.append(body);
  }

  else
  {
    for(int i= 0; i<tbs.size(); i++)
      new_bodies.append(CAST_TO(tbs.get_and_step(), OCCBody));
  }

  PRINT_WARNING("Occ engine doesn't consider tighten_gaps and tolerance in stitch operation. \n");
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : stitch_surfs
// Member Type: PUBLIC
// Description: stitch all surfs and try to make a TopoDS_Shell .
//              called by stitch into surface body
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::stitch_surfs(
                      DLIList<BodySM*>& surf_bodies,
                      TopoDS_Shape& fuse) const
{
  if (surf_bodies.size() < 2)
    return CUBIT_SUCCESS;

  DLIList<TopoDS_Shape*> faces_to_stitch;
  DLIList<OCCSurface*> original_surfaces;
  for (int i = 0; i < surf_bodies.size(); i++)
  {
     BodySM * tool_body = surf_bodies.get_and_step();
     OCCBody* occ_body = CAST_TO(tool_body, OCCBody);
     DLIList<OCCSurface*> surfaces = occ_body->my_sheet_surfaces();
     DLIList<OCCShell*> shells = occ_body->shells();
     DLIList<Lump*> lumps = occ_body->lumps();
     if (surfaces.size()+shells.size() != 1 || lumps.size() > 0)
     {
       PRINT_ERROR("Can't stitch non-sheet bodySM's. \n");
       return CUBIT_FAILURE;
     }

     original_surfaces += surfaces;
     delete occ_body;
     OCCQueryEngine::instance()->BodyList->remove(occ_body);
     if (surfaces.size() == 1)
     {
       OCCSurface* surface = surfaces.get();
       delete surface->my_shell();
       delete surface->my_lump();
       surface->set_shell(NULL);
       surface->set_lump(NULL);
       surface->set_body(NULL);

       TopoDS_Face* topods_face = surface->get_TopoDS_Face();
       if (topods_face != NULL)
         faces_to_stitch.append(topods_face);
     }
     else
     {
       OCCShell* shell = shells.get();
       delete shell->my_lump();
       shell->set_body(NULL);
       shell->set_lump(NULL);

       TopoDS_Shell* topods_shell = shell->get_TopoDS_Shell();
       if(topods_shell)
          faces_to_stitch.append(topods_shell);
     }
  }

  faces_to_stitch.reset();

  BRepBuilderAPI_Sewing sew;

  for( int i = faces_to_stitch.size()-1; i >= 0; i --)
  {
      TopoDS_Shape* face = faces_to_stitch[i];
      sew.Add(*face);
  }

  sew.Perform();
  TopoDS_Shape sewn_shape=sew.SewedShape();

  fuse=sewn_shape;

  TopExp_Explorer Ex;
  int count_shell = 0;
  for (Ex.Init(fuse, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
    count_shell++;

  if ( count_shell != 1)
  {
     PRINT_ERROR("Can't stitch all surfaces into one BodySM's. \n");
     return CUBIT_FAILURE;
  }

  for( int i = original_surfaces.size()-1; i >= 0; i --)
    OCCQueryEngine::instance()->
      delete_solid_model_entities( (Surface*)original_surfaces.pop());

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : subtract
// Member Type: PUBLIC
// Description: subtract boolean operation on OCC-based bodies
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus     OCCModifyEngine::subtract(DLIList<BodySM*> &tool_body_list,
                                          DLIList<BodySM*> &from_bodies,
                                          DLIList<BodySM*> &new_bodies,
                                          bool imprint,
                                          bool keep_old) const
{
  // copy the bodies in case keep_old is true
  DLIList<TopoDS_Shape*> tool_bodies_copy;

  //for subtract function, tool-body has to be solid, 
  //otherwise it's just imprint
  DLIList<CubitBox*>* tool_boxes = new DLIList<CubitBox*>();
  DLIList<CubitBoolean> is_tool_volume;
  //keep the tool_body untouched
  CubitStatus stat = 
    get_shape_list(tool_body_list, tool_bodies_copy, is_tool_volume, CUBIT_TRUE, tool_boxes);

  if(!stat)
    return stat;

  stat = do_subtract(from_bodies, tool_bodies_copy, is_tool_volume,
                     tool_boxes, new_bodies, keep_old, imprint) ;

  //ok, we're done with all cuts, delete unnecessaries.
  CubitBoolean delete_tool_boxes = CUBIT_FALSE;
  if(tool_boxes->size() > 0)
    delete_tool_boxes = CUBIT_TRUE;
  while (tool_boxes->size())
    delete tool_boxes->pop();
  if(delete_tool_boxes)
    delete tool_boxes;
  while (tool_bodies_copy.size())
    delete tool_bodies_copy.pop();
  if(!keep_old) //delete tool_bodies
    OCCQueryEngine::instance()->delete_solid_model_entities(tool_body_list);
  return stat;
}

CubitStatus OCCModifyEngine::do_subtract(DLIList<BodySM*> &from_bodies,
                                      DLIList<TopoDS_Shape*> &tool_bodies_copy,
                                      DLIList<CubitBoolean> &is_tool_volume,
                                      DLIList<CubitBox*>* tool_boxes,
                                      DLIList<BodySM*> &new_bodies,
                                      bool keep_old,
                                      bool imprint) const
{
  DLIList<TopoDS_Shape*> from_bodies_copy;
  DLIList<CubitBoolean> is_volume;
  //get the from_bodies underling shapes
  CubitStatus stat = get_shape_list(from_bodies, from_bodies_copy, is_volume, keep_old);
  if(!stat)
  {
    for (int i = 0; i < tool_bodies_copy.size(); i++)
    {
       TopoDS_Shape* shape = tool_bodies_copy.get_and_step();
       delete shape;
       shape = NULL;
    }
    tool_bodies_copy.clean_out();
    return CUBIT_FAILURE;
  } 

  //check that tool_bodies are all solid, shell and surface body can't be used
  //for subtracting solids.
  if(is_tool_volume.is_in_list(CUBIT_FALSE) && !is_volume.is_in_list(CUBIT_FALSE))
  {
     PRINT_WARNING("Surfaces or Shells can't be used to cut a solid.\n");
     while (tool_bodies_copy.size())
       delete tool_bodies_copy.pop();
     while (from_bodies_copy.size())
       delete from_bodies_copy.pop(); 
     return CUBIT_FAILURE;
  }

  int fraction_remaining = 100;

  // subtract the tool body from each body in the list

  TopoDS_Shape*  from_shape = from_bodies_copy.get();
  DLIList<TopologyBridge*> tbs;
  for (int i = 0; i < from_bodies_copy.size(); i++)
  {
    CubitBoolean from_volume = is_volume.get_and_step();
    BodySM* from_body = from_bodies.get();
    CubitBox box1 = CAST_TO(from_body, OCCBody)->get_bounding_box();
    int count = 0;  //count for not preforming cut
    for(int j = 0; j < tool_bodies_copy.size(); j ++)
    {
      CubitBoolean tool_volume = is_tool_volume.get_and_step();
      if(tool_volume == CUBIT_FALSE && from_volume == CUBIT_TRUE)
      {
        PRINT_WARNING("Surfaces or Shells can't be used to cut a solid.\n");
        continue;
      }
      if (AppUtil::instance()->interrupt())
      {
         PRINT_ERROR("Subtraction interrupted.  Aborting...\n");
         while (tool_bodies_copy.size())
            delete tool_bodies_copy.pop();
         while (from_bodies_copy.size())
            delete from_bodies_copy.pop();
         return CUBIT_FAILURE;
      }
      CubitBox tool_box = *tool_boxes->get_and_step();  
      if(!tool_box.overlap(TOL,box1))
      {
        count++;
        continue;
      } 
      TopoDS_Shape* tool_shape = tool_bodies_copy.get_and_step();

      //bodies overlap, proceed with the subtract
      BRepAlgoAPI_Cut cutter(*from_shape, *tool_shape);
      TopoDS_Shape cut_shape = cutter.Shape(); 

      //compare to see if the from_shape has gotten cut.
      CubitBoolean has_changed = CUBIT_FALSE;
      double after_mass = 0.0;
      GProp_GProps myProps;
      if(is_volume[i])
        BRepGProp::VolumeProperties(cut_shape, myProps);

      else
        BRepGProp::SurfaceProperties(cut_shape, myProps);
      after_mass = myProps.Mass();
      if(after_mass > TOL)
        check_operation(cut_shape, from_shape, is_volume[i], has_changed,
                    &cutter, keep_old);
      else
      {
/*
        if(!keep_old)
          OCCQueryEngine::instance()->delete_solid_model_entities(from_body);
        from_shape->Nullify();
*/
        from_shape = NULL;
      }

      int stat;
      if(!has_changed && from_shape && !from_shape->IsNull())
      {
        //Add imprint code here 
        DLIList<TopoDS_Face*> face_list;
        if(imprint)
        {
          stat = imprint_toposhapes(from_shape, tool_shape, face_list);
          if(stat)
          {
            PRINT_ERROR("Can't do imprint operation on the body. \n");
            count++;
          }
          continue;
        }
      }
    }

    //ok, we're done with all cuts, construct new Body'
    if (count < tool_bodies_copy.size() && from_shape && !from_shape->IsNull())
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*from_shape);
    else if (from_shape && !from_shape->IsNull())
    {
      PRINT_INFO("The %d body did not change because cutting tools are not interscting with it.\n", i+1);
      from_bodies.change_to(NULL);
    }
    from_bodies.step();
    from_shape = from_bodies_copy.step_and_get();

    // done with this j iteration; write out count, if necessary
    if (from_bodies.size() * tool_bodies_copy.size() > 1)
    {
       int frac_done = (100 * (i+1)) / (from_bodies.size());
       if ((100 - frac_done) < fraction_remaining)
       {
          fraction_remaining = 100 - frac_done;
          PRINT_INFO("%d%% remaining.\n ", fraction_remaining+1);
       }
    }
  }

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }    

  //ok, we're done with all cuts, delete unnecessaries. 
  if(keep_old)
  {
    int size  = from_bodies_copy.size();
    for (int i = 0; i < size; i++)
    {
      TopoDS_Shape* shape = from_bodies_copy.pop();
      if(shape != NULL )
      {
        shape->Nullify();
        //delete shape; - Must not be deleted, causes random failure on OSX
        shape = NULL;
      }
    }
  } 
  return CUBIT_SUCCESS; 
}

//===============================================================================
// Function   : imprint_toposhapes
// Member Type: PROTECTED
// Description: imprint boolean operation on OCC-based bodies.
//              from_shape must be TopoDS_Face or above, tool_shape can be
//              TopoDS_Edge or above. 
//              on_faces works only when tool_shape is an Edge, indicates that
//              those edges only imprint on the on_faces.
//              success: return 0
//              needs intersect-subtract-union boolean: return 1
//              fail:     return 2
//              needs face intersect for periodic geometry: return 3 
//              the cutting surface's topoDS_Face in returned in the
//              the last of on_faces list.
// Author     : Jane HU
// Date       : 03/08
//===============================================================================
int OCCModifyEngine::imprint_toposhapes(TopoDS_Shape*& from_shape, 
                                        TopoDS_Shape* tool_shape,
                                        DLIList<TopoDS_Face*>&on_faces)const
{
  int count = 0;   //number of imprinting
 
  //indicate if there's more faces to be imprinted
  CubitBoolean more_face = CUBIT_TRUE; 
  CubitBoolean face_created = CUBIT_FALSE;

  //list of face on from_shape that has been imprinted
  DLIList<TopoDS_Face*> from_faces; 
  while( more_face)
    {
      TopoDS_Face from_face,tool_face;
      TopoDS_Edge* common_edge = NULL;
      DLIList<TopoDS_Shape*> tool_faces_edges;
      TopTools_ListOfShape list_of_edges;
      BRepFeat_SplitShape splitor(*from_shape);
      CubitBoolean qualified = CUBIT_TRUE;
      if (tool_shape->ShapeType() == TopAbs_EDGE)
	{
	  if(count == 1)
	    break;
      
	  DLIList<TopoDS_Face*> faces;
	  //need to delete TopoDS_Face* in faces
	  common_edge = find_imprinting_edge(*from_shape, TopoDS::Edge(*tool_shape),faces);
	  if (common_edge)
	    {
	      if (on_faces.size() > 0)
		qualified = CUBIT_FALSE;
	      for(int j = 0; j < faces.size(); j++)
		{
		  from_face = *faces.get();
		  for (int i = 0; i < on_faces.size(); i++)
		    {
		      if (from_face.IsSame(*(on_faces.get_and_step())))
			{
			  qualified = CUBIT_TRUE; 
			  break;
			}
		    }
		  faces.get()->Nullify();
		  delete faces.get();
		  faces.step();
		}
	      if (qualified && (from_faces.size() == 0 || (from_faces.size() && !from_face.IsSame(*from_faces.get()))) )
		list_of_edges.Append(*common_edge);
	      else
		from_face.Nullify();
	      common_edge->Nullify();
	      delete common_edge;
              common_edge = NULL;
	    }
	}
      else 
	{
	  TopOpeBRep_ShapeIntersector intersector;
	  intersector.InitIntersection(*from_shape, *tool_shape);

	  //find the intersecting edges and faces.
	  int max_edge = 0;
                     
	  for(; intersector.MoreIntersection(); )
	    {
	      TopoDS_Shape face1;
	      face1 = intersector.ChangeFacesIntersector().Face(1);
	      CubitBoolean has_imprinted = CUBIT_FALSE;
	      for (int j = 0; j < from_faces.size(); j++)
		{
		  TopoDS_Face* topo_face = from_faces.get_and_step();
		  if(face1.IsSame(*topo_face))
		    {
		      has_imprinted = CUBIT_TRUE;
		      break;
		    }
		}

	      if (has_imprinted == CUBIT_TRUE)
		{
		  intersector.NextIntersection();
		  continue;
		}
	      TopoDS_Shape edge_face;

	      edge_face = intersector.ChangeFacesIntersector().Face(2);
	      BRepAlgoAPI_Section section(face1, edge_face);

	      //intersection edges between face1 and edge_face
	      TopTools_ListOfShape temp_list_of_edges;
#if OCC_VERSION_MINOR > 5
              TopoDS_Shape s_shape = section.Shape();
              TopExp_Explorer Ex;
              for (Ex.Init(s_shape, TopAbs_EDGE, TopAbs_WIRE); Ex.More(); Ex.Next())
                temp_list_of_edges.Append(TopoDS::Edge( Ex.Current())); 
#else
              temp_list_of_edges.Assign(section.SectionEdges());
#endif
	      int num_edges = temp_list_of_edges.Extent();
  
	      CubitBoolean is_same = face1.IsSame(from_face);
	      CubitBoolean is_same_tool = CUBIT_FALSE;
	      for (int j = 0; j < tool_faces_edges.size(); j++)
		{
		  TopoDS_Shape* topo_face_edge = tool_faces_edges.get_and_step();
		  if(edge_face.IsSame(*topo_face_edge))
		    {
		      is_same_tool = CUBIT_TRUE;
		      break;
		    }
		}
	      if (max_edge < num_edges )
		{
		  list_of_edges.Assign(temp_list_of_edges);  
		  max_edge =  num_edges ;
		  from_face = TopoDS::Face(face1);
		  TopoDS_Shape* topo_shape = new TopoDS_Shape(edge_face);
		  DLIList<TopoDS_Shape*> shape_list;
		  for(int iii = 0; iii < tool_faces_edges.size(); iii++)
		    {
		      int size = shape_list.size();
		      shape_list.append_unique(tool_faces_edges.get_and_step());
		      if (size < shape_list.size())
			{
			  shape_list.last();
			  shape_list.get()->Nullify();
			  delete shape_list.get();
			}
		    }
		  tool_faces_edges.clean_out();
		  for(int i = 0 ; i < num_edges; i++)
		    //later has to use it num_edges times 
		    tool_faces_edges.append(topo_shape);
		}
	      else if(num_edges == max_edge && is_same && !is_same_tool) 
		//multi tool faces cut the same face
		{
		  TopTools_ListIteratorOfListOfShape Itor, temp_Itor;
		  temp_Itor.Initialize(temp_list_of_edges);
		  for(; temp_Itor.More(); temp_Itor.Next())
		    {
		      TopoDS_Edge temp_edge = TopoDS::Edge(temp_Itor.Value());
		      Itor.Initialize(list_of_edges);
		      CubitBoolean same_edge = CUBIT_FALSE;
              
		      GProp_GProps myProps1;
		      BRepGProp::LinearProperties(temp_edge, myProps1);
		      gp_Pnt center1 = myProps1.CentreOfMass();
		      for(; Itor.More(); Itor.Next())
			{
			  TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
			  GProp_GProps myProps2;
			  BRepGProp::LinearProperties(edge, myProps2);
			  gp_Pnt center2 = myProps2.CentreOfMass();
			  if(center1.IsEqual(center2, TOL))
			    {
			      same_edge = CUBIT_TRUE;
			      break;
			    }
			}
		      if(!same_edge)
			{
			  list_of_edges.Append(temp_edge);
			  TopoDS_Shape* topo_shape = new TopoDS_Shape(edge_face);
			  tool_faces_edges.append(topo_shape);
			}
		    }//end 'for'
		}//end  'else if'
              intersector.NextIntersection();
	    } //end 'for'
	}//end 'else'
      if (from_face.IsNull())
	{
	  more_face = CUBIT_FALSE;
	  DLIList<TopoDS_Shape*> shape_list;
	  for(int iii = 0; iii < tool_faces_edges.size(); iii++)
	    {
	      int size = shape_list.size();
	      shape_list.append_unique(tool_faces_edges.get_and_step());
	      if (size < shape_list.size())
		{
		  shape_list.last();
		  shape_list.get()->Nullify();
		  delete shape_list.get();
		}
	    }
	  tool_faces_edges.clean_out();

	  for (int iii=0; iii < from_faces.size(); iii++)
	    {
	      TopoDS_Face* topo_face = from_faces.get_and_step();
	      topo_face->Nullify();
	      delete topo_face;
              topo_face = NULL;
	    }
	  continue;
	}
  
      TopTools_ListIteratorOfListOfShape Itor;

      //list_of_edges is the intersection edges on tool_faces_edges 
      Itor.Initialize(list_of_edges);
      int total_edges = list_of_edges.Extent();
      DLIList<Curve*> curve_list;
      DLIList<OCCCurve*> occ_curves;
      CubitBoolean topo_changed = CUBIT_FALSE;
      tool_faces_edges.reset();

      //check to see if the intersection edge is:
      //1. on from_face: if not, skip it
      //2. overlap with from_edges : if not, add the edge for splitting face
      //3. if overlap, is it the same edge:if not add it for splitting edge
      // if yes, skip it too

      Surface* face = NULL;
      if (OCCQueryEngine::instance()->OCCMap->IsBound(from_face))
	{
	  int i = OCCQueryEngine::instance()->OCCMap->Find(from_face);
	  face = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(i))->second;
	}
      if (face == NULL)
      {
        face_created = CUBIT_TRUE;
        face = OCCQueryEngine::instance()->populate_topology_bridge(from_face);
      }
      OCCSurface* occ_face = CAST_TO(face, OCCSurface);

      GeometryType type = occ_face->geometry_type();
      CubitBoolean periodic = occ_face->is_periodic();

      DLIList<Curve*> common_curves;
      for(; Itor.More(); Itor.Next())
      {
        TopoDS_Edge edge = TopoDS::Edge(Itor.Value());

        //copy the edge for imprinting.
        BRepBuilderAPI_Copy api_copy(edge);
        TopoDS_Shape newShape = api_copy.ModifiedShape(edge);
        edge = TopoDS::Edge(newShape);
        Curve* curve = NULL;
        if(!OCCQueryEngine::instance()->OCCMap->IsBound(edge))
        {
          curve = OCCQueryEngine::instance()->populate_topology_bridge(edge, true);
          DLIList<OCCPoint*> points;
          OCCCurve* occ_c = CAST_TO(curve, OCCCurve);
          occ_c->get_points(points);
          for(int i = 0 ; i <  points.size(); i++)
            points.get_and_step()->remove_curve(occ_c);
        }
        else
          curve = OCCQueryEngine::instance()->populate_topology_bridge(edge);
        if(curve)
          common_curves.append(curve);
      }
      DLIList<DLIList<TopoDS_Edge*>*> temp_edge_lists;
      if (common_curves.size() >= 1)
         sort_curves(common_curves, temp_edge_lists);

      if ( common_curves.size() >= 2 && 
           (type == CONE_SURFACE_TYPE || type == SPHERE_SURFACE_TYPE ||
            type == TORUS_SURFACE_TYPE ||type == UNDEFINED_SURFACE_TYPE ||
            type == SPLINE_SURFACE_TYPE))
      {
        //if the two shapes has common volume, do boolean operations
        BRepAlgoAPI_Common intersector(*from_shape, *tool_shape);
        TopTools_ListOfShape shapes;
        shapes.Assign(intersector.Modified(*tool_shape));
     
        TopoDS_Shape common_shape;
        if (shapes.IsEmpty())
          common_shape = intersector.Shape();
        else 
          common_shape = shapes.First();

        if(!common_shape.IsNull()) 
        { 
          TopTools_IndexedMapOfShape M;
          TopExp::MapShapes(common_shape, TopAbs_SOLID, M);
          if(M.Extent() > 0)
          {
            GProp_GProps myProps;
            BRepGProp::VolumeProperties(common_shape, myProps);
            double after_mass = myProps.Mass();
            BRepGProp::VolumeProperties(*from_shape, myProps);
            double orig_mass = myProps.Mass();
            if(fabs(-after_mass + orig_mass) > TOL && after_mass > TOL)
              return 1; 
          }
          //have to use boolean operation, see 
          //http://www.opencascade.org/org/forum/thread_20672/ for more info
        }
      }

      list_of_edges.Clear(); 
      if (common_curves.size() >= 1)
	{
	  DLIList<TopoDS_Edge*>* edge_list;
	  int size = temp_edge_lists.size();
	  for(int i = 0; i < size; i++)
	    {
	      edge_list = temp_edge_lists.get_and_step();
	      //make sure the copied edges are sharing vertices.
	      BRepBuilderAPI_MakeWire myWire;
	      edge_list->reset();
	      for(int j = 0; j < edge_list->size(); j++)
		{
		  TopoDS_Edge e = *(edge_list->get_and_step());
                  //Don't include zero length edge.
                  GProp_GProps myProps;
                  BRepGProp::LinearProperties(e, myProps);
                  double d = myProps.Mass();
                  if(d > TOL) 
		    myWire.Add(e);
		}
	      TopoDS_Wire wire = myWire.Wire();
	      BRepTools_WireExplorer Ex(wire); 
	      for(; Ex.More(); Ex.Next())
		list_of_edges.Append(Ex.Current());
              edge_list->clean_out();
              delete edge_list;
	    }
	}
      for(Itor.Initialize(list_of_edges); Itor.More(); Itor.Next())
	{
	  TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
          CubitBoolean added = CUBIT_FALSE;
	  //check to see if the intersection edge is on from_face
	  TopExp_Explorer Ex;
	  CubitBoolean skipped = CUBIT_FALSE;
	  GProp_GProps myProps1;
	  BRepGProp::LinearProperties(edge, myProps1);
	  double d1 = myProps1.Mass();
          Curve* curve = NULL;
          //edge should all be bounded, it comes from common_curves
          if (OCCQueryEngine::instance()->OCCMap->IsBound(edge))
          {
            int i = OCCQueryEngine::instance()->OCCMap->Find(edge);
            curve = (Curve*)
                    (OCCQueryEngine::instance()->OccToCGM->find(i))->second;
          }
          else 
          {
            curve = OCCQueryEngine::instance()->populate_topology_bridge(edge, true);
            DLIList<OCCPoint*> points;
            OCCCurve* occ_c = CAST_TO(curve, OCCCurve);
            occ_c->get_points(points);
            for(int i = 0 ; i <  points.size(); i++)
              points.get_and_step()->remove_curve(occ_c);
          }
	  gp_Pnt pt = myProps1.CentreOfMass();
	  CubitVector p = curve->center_point();

	  CubitVector point_on_surf;
	  occ_face->closest_point_trimmed(p, point_on_surf);

	  if(p.distance_between(point_on_surf) > TOL) //edge not on from_face
	    {
	      skipped = CUBIT_TRUE;
	      total_edges--;
	    }

	  else 
	    {
	      for (Ex.Init(from_face, TopAbs_EDGE); Ex.More(); Ex.Next())
		{
		  //check if the edge is on from_face edges, add such edge on existing
		  //edge to split it.
		  TopoDS_Edge from_edge = TopoDS::Edge(Ex.Current());
         
		  GProp_GProps myProps2;
		  BRepGProp::LinearProperties(from_edge, myProps2);
		  double d2 = myProps2.Mass();
                  gp_Pnt pt2 = myProps2.CentreOfMass();
                  CubitBoolean found_ = CUBIT_FALSE;
                  Curve* from_curve;
		  if (OCCQueryEngine::instance()->OCCMap->IsBound(from_edge))
		  {
                    int i = OCCQueryEngine::instance()->OCCMap->Find(from_edge);
		    from_curve = (OCCCurve*)
                        (OCCQueryEngine::instance()->OccToCGM->find(i))->second;
                    found_ = CUBIT_TRUE;
		  }
		  else
		    from_curve = OCCQueryEngine::instance()->populate_topology_bridge(from_edge, true);
                  if(!from_curve)
                    continue;

		  OCCCurve* occ_curve = CAST_TO(from_curve, OCCCurve);

                  if(pt.IsEqual(pt2, TOL) && fabs(d2-d1)< TOL) //overlap
                  {
                    skipped = CUBIT_TRUE;
                    total_edges--;
                    break;
                  }

		  CubitPointContainment pc = CUBIT_PNT_OFF;
                  pc = occ_curve->point_containment(p);
                  if((d2 - d1) > TOL && pc == CUBIT_PNT_ON) 
                  {
                    added = CUBIT_TRUE;
                    splitor.Add(edge, from_edge);
                    total_edges--;
                    break;
                  }

                  if(!found_ && from_curve)
                    OCCQueryEngine::instance()->
                      delete_solid_model_entities(from_curve);
		} 
	      if(list_of_edges.Extent() == 1 && !skipped) 
		{
		  added = CUBIT_TRUE;
		  curve_list.append(curve); 
		}
	    } 
	  if(added)
	    {
	      topo_changed = CUBIT_TRUE;
	      continue;
	    }
	  if (!skipped)
	    {
	      //check if edge's inside from_face by checking bounding boxes  
	      BRepAdaptor_Curve acurve(edge);
	      BRepAdaptor_Surface asurface( from_face);
	      Bnd_Box aBox_edge, aBox_face;
	      BndLib_Add3dCurve::Add(acurve, Precision::Approximation(), aBox_edge);
	      BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox_face);
	      double min[3], max[3];
	      aBox_edge.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
	      CubitBox aBox_e(min, max);
	      aBox_face.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
	      CubitBox aBox_f(min, max);
              //hexlat3 has tolerance issue where aBox_e.x_min is within 
              //tolerance and greater than aBox_f.x_min, causing no edge
              //imprint of the faces. change to add consideration of tolerance
              // lose the bounding box critiral a little more.
              int num_satisfied = 0;
              
              if (aBox_f.min_x() >= 0 && aBox_e.min_x() >= aBox_f.min_x() *0.9)
                  num_satisfied ++;
              if(aBox_f.min_x() < 0 && aBox_e.min_x() >= aBox_f.min_x() *1.1)
		  num_satisfied ++;
     
              if (aBox_f.min_y() >= 0 && aBox_e.min_y() >= aBox_f.min_y() *0.9)
                  num_satisfied ++;
              if(aBox_f.min_y() < 0 && aBox_e.min_y() >= aBox_f.min_y() *1.1)
                  num_satisfied ++;

              if(aBox_f.min_z() >= 0 && aBox_e.min_z() >= aBox_f.min_z() *0.9)
                  num_satisfied ++;
              if (aBox_f.min_z() < 0 && aBox_e.min_z() >= aBox_f.min_z() *1.1)
                  num_satisfied ++;

              if (aBox_f.max_x() > 0  && aBox_e.max_x() <= aBox_f.max_x() *1.1)
                num_satisfied ++;
              if(aBox_f.max_x() <= 0 && aBox_e.max_x() <= aBox_f.max_x() *0.9)
                num_satisfied ++;
              if( aBox_f.max_y() > 0 && aBox_e.max_y() <= aBox_f.max_y() *1.1) 
                num_satisfied ++;
              if( aBox_f.max_y() <= 0 && aBox_e.max_y() <= aBox_f.max_y() *0.9)
                num_satisfied ++;
              if( aBox_f.max_z() >= 0 && aBox_e.max_z() <= aBox_f.max_z() *1.1)
                num_satisfied ++;  
              if( aBox_f.max_z() < 0 && aBox_e.max_z() <= aBox_f.max_z() *0.9)
                num_satisfied ++;

              if(num_satisfied == 6) 
		{
		  curve_list.append_unique(curve);
		}
              else
              {
                curve_list.remove(curve);
                OCCQueryEngine::instance()->delete_solid_model_entities( curve );
                total_edges--;
              }
	    }
	}

      if(face_created)
      {
        OCCQueryEngine::instance()->delete_solid_model_entities(occ_face);
        face_created = CUBIT_FALSE;
      }
      DLIList<DLIList<TopoDS_Edge*>*> edge_lists;
      if ( total_edges >= 1) 
	{      
	  CubitStatus stat = CUBIT_SUCCESS;
          if(curve_list.size() > 0)
	    stat = sort_curves(curve_list, edge_lists); 
          else
	    {
	      TopoDS_Face* topo_face = new TopoDS_Face(from_face);
	      from_faces.append(topo_face);
	      continue;
	    }
	  DLIList<TopoDS_Edge*>* edge_list;
          int size = edge_lists.size();
        
          //if size > 1 , outer wire first. 
          DLIList<CubitBox*> bs;
          edge_lists.reset();
          CubitBox box;
          for(int i = 0; i < size && size > 1; i++)
          {
            edge_list = edge_lists.get_and_step();
            for(int j = 0; j < edge_list->size(); j++)
            {
              TopoDS_Edge e = *(edge_list->get_and_step());
              CubitBoolean bound = CUBIT_FALSE;
              Curve* curve = NULL;
              if(OCCQueryEngine::instance()->OCCMap->IsBound(e))
              {
                int k = OCCQueryEngine::instance()->OCCMap->Find(e);
                bound = CUBIT_TRUE;
                curve = (Curve*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
              }
              else
                curve = OCCQueryEngine::instance()->
                         populate_topology_bridge(e, true);
              if(j == 0)
                box = curve->bounding_box();
              else
                box |= curve->bounding_box();
              if(!bound)
                OCCQueryEngine::instance()->delete_solid_model_entities(curve);
            }
            CubitBox* pbox = new CubitBox(box);
            bs.append(pbox);
          }
          bs.append((CubitBox*) NULL);

          for(int i = 0; i < size-1; i++)
          {
            for(int j = i + 1; j < size; j ++)
            {
              edge_lists.reset();
              if(*bs[j] <= *bs[i])
              {
                edge_lists.step(j);
                edge_lists.get()->clean_out();
                delete edge_lists.get();
                edge_lists.change_to((DLIList<TopoDS_Edge*>*)NULL);
                continue;
              }
              if(*bs[i] <= *bs[j])
              {
                edge_lists.step(i);
                edge_lists.get()->clean_out();
                delete edge_lists.get();
                edge_lists.change_to((DLIList<TopoDS_Edge*>*)NULL);
                break;
              }
            }
          }
          for(int i = 0; i < size && size > 1; i++)
            delete bs.pop();

          edge_lists.remove_all_with_value((DLIList<TopoDS_Edge*>*)NULL);

          size = edge_lists.size();
          for (int iii = 0; iii < size; iii++)
	    {
	      edge_list = edge_lists.pop();

	      //check if the edges make a split of the surface, error out if it's
	      //just a scar on the surface
	      int count_intersection = 0;
	      if (stat == CUBIT_FAILURE) //Open wire
		{
		  count_intersection = check_intersection(edge_list, from_face);
            
		  if (count_intersection == 1 )
                  {
		    PRINT_WARNING("Cant make a scar on existing face without splitting it. \n");
                    edge_list->clean_out();
                    delete edge_list;
                  }
		} 
	      if (stat || count_intersection == 2)
		{
		  BRepBuilderAPI_MakeWire myWire;
		  edge_list->reset(); 
                  DLIList<Curve*> wire_curves;
		  for(int i = 0; i < edge_list->size(); i++)
		    {
		      TopoDS_Edge e = *(edge_list->get_and_step());
                      //remove curve from its vertice's curve list
                      if(OCCQueryEngine::instance()->OCCMap->IsBound(e))
                      {
                        int j = OCCQueryEngine::instance()->OCCMap->Find(e);
                        Curve* curve = (Curve*)(OCCQueryEngine::instance()->OccToCGM->find(j))->second;
                        OCCCurve* curve_to_remove = (OCCCurve*) curve;
                        wire_curves.append(curve);
                        DLIList<OCCPoint*> points;
                        curve_to_remove->get_points(points);
                        for (int k = 0; k < points.size(); k++)
                          points.get_and_step()->remove_curve(curve_to_remove);
                      }
                      else if(periodic && count_intersection != 2)
                      {
                        Curve* curve = OCCQueryEngine::instance()->
                               populate_topology_bridge(e,true);
                        wire_curves.append(curve);
                      }
		      myWire.Add(e); 
		    }
                  edge_list->clean_out();
                  delete edge_list;
                  if(!periodic || count_intersection == 2 ||
                     type == PLANE_SURFACE_TYPE)
                  {
		    splitor.Add(myWire.Wire(),from_face);
		    topo_changed = CUBIT_TRUE; 
                     
                  }
                  else 
                  {
                    //use the myWire to create a surface and webcut the 
                    //periodic body.
                    Surface *wire_surf = make_Surface(BEST_FIT_SURFACE_TYPE,
                       wire_curves);
                    if(wire_surf == NULL)
                      wire_surf = make_Surface(PLANE_SURFACE_TYPE, wire_curves);
                    if(wire_surf)
                    {
                      OCCSurface* occ_wire_s = CAST_TO(wire_surf, OCCSurface);
                      TopoDS_Face *topo_face = occ_wire_s->get_TopoDS_Face();
                      on_faces.append(topo_face);
                      return 3; 
                    }
                  } 
		}
	    }
	} 
      if(topo_changed)
	{
	  splitor.Build();
	  topo_changed = CUBIT_FALSE;
	  if(splitor.IsDone())
	    {
	      //take care of on_faces list first:after operation, the on_faces
	      // will have at least one face changed, update the pointer.
              int size = on_faces.size();
	      if (size > 0)
		{
		  for(int k = 0; k < size; k++)
		    {
		      TopoDS_Face* compare_face = on_faces.get();
                      TopTools_ListOfShape shapes;
                      shapes.Assign(splitor.Modified(*compare_face));
                      if(shapes.Extent() > 0 && 
                         !compare_face->IsSame(shapes.First()))
                      {
		        on_faces.change_to(NULL);
		        while(shapes.Extent() > 0)
		        {
		          TopoDS_Face* face = 
				new TopoDS_Face(TopoDS::Face(shapes.First())); 
			  shapes.RemoveFirst();
			  on_faces.append(face);
			}
                      }
                      on_faces.step();
		    }
                    on_faces.remove_all_with_value(NULL);
		}

	      TopoDS_Shape new_from_shape = splitor.Shape();
	      if(from_shape->ShapeType() == TopAbs_COMPOUND)
		{
		  TopoDS_Compound old_csolid = TopoDS::Compound(*from_shape);
		  OCCBody::update_OCC_entity(old_csolid, new_from_shape, &splitor);
                  if(!old_csolid.IsEqual(new_from_shape))
                  {
                     from_shape->Nullify();
                     *from_shape = new_from_shape;
                  }                     
		}

	      else if(from_shape->ShapeType() == TopAbs_SOLID)
		{
		  TopoDS_Solid old_solid = TopoDS::Solid(*from_shape);
		  OCCLump::update_OCC_entity(old_solid, new_from_shape, &splitor);
                  if(!old_solid.IsEqual(new_from_shape))
                  { 
                     from_shape->Nullify();
                     *from_shape = new_from_shape;
                  }
		}
	      else if(from_shape->ShapeType() == TopAbs_SHELL)
		{
		  TopoDS_Shell old_shell = TopoDS::Shell(*from_shape);
		  OCCShell::update_OCC_entity(old_shell,new_from_shape, &splitor);
                  if(!old_shell.IsEqual(new_from_shape))
                  {
                     from_shape->Nullify();
                     *from_shape = new_from_shape;
                  }
		}
	      else if(from_shape->ShapeType() == TopAbs_FACE)
		{
		  TopoDS_Face old_face = TopoDS::Face(*from_shape);
		  OCCSurface::update_OCC_entity(old_face,new_from_shape, &splitor);
                  if(!old_face.IsEqual(new_from_shape))
                  {
                     from_shape->Nullify();
                     *from_shape = new_from_shape;
                  }
		}

	      TopTools_ListOfShape shapes;
	      for(int i = 0; i < from_faces.size(); i++)
		{
		  TopoDS_Face* topo_face = from_faces.get();
		  shapes.Assign(splitor.Modified(*topo_face));
		  topo_face = new TopoDS_Face(TopoDS::Face(shapes.First()));
		  from_faces.get()->Nullify();
		  delete from_faces.get();
		  from_faces.change_to(topo_face);
		  from_faces.step();
		} 
	      count++;
	    }
	}
      else
	{
	  TopoDS_Face* topo_face = new TopoDS_Face(from_face);
	  from_faces.append(topo_face);
	} 
    }
  
  TopExp_Explorer Ex;
  int num_face = 0;
  for (Ex.Init(*from_shape, TopAbs_FACE); Ex.More(); Ex.Next())
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      num_face++;
    }
  
#ifdef DEBUG  
  PRINT_INFO("Total %d cuts performed, with from_shape having %d faces.\n", count, num_face); 
#endif

  if (count > 0)
    return 0;
  return 2;
}

//===============================================================================
// Function   : find_imprinting_edge
// Member Type: PROTECTED
// Description: imprint boolean operation on OCC-based bodies.
//              from_shape must be TopoDS_Face or above, tool_shape must be
//              TopoDS_Edge.
// Author     : Jane HU
// Date       : 05/08
//===============================================================================
TopoDS_Edge* OCCModifyEngine::find_imprinting_edge(TopoDS_Shape& from_shape,
                                        TopoDS_Edge& tool_shape,
                                        DLIList<TopoDS_Face*>& from_faces)const
{
  TopoDS_Edge* edge = NULL;
  //list of face on from_shape that has been imprinted
  from_faces.clean_out();
  BRepAdaptor_Curve acurve(tool_shape);
  Bnd_Box aBox_edge, aBox_face;
  BndLib_Add3dCurve::Add(acurve, Precision::Approximation(), aBox_edge);

  TopExp_Explorer Ex;
  for (Ex.Init(from_shape, TopAbs_FACE); Ex.More(); Ex.Next())
  {
    TopoDS_Face face = TopoDS::Face(Ex.Current());
    BRepAdaptor_Surface asurface( face);
    aBox_face.SetVoid();
    BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox_face);
    if (aBox_face.IsOut(aBox_edge))
      continue;

    BRepAlgoAPI_Common intersector(face, tool_shape);
    TopTools_ListOfShape shapes;
    shapes.Assign(intersector.Generated(tool_shape));
    if(shapes.IsEmpty())
      shapes.Assign(intersector.Modified(tool_shape));
    if (shapes.IsEmpty())
      continue;
    if ( shapes.Extent() > 1)
    {
      PRINT_ERROR("Edge has multiple intersection with the shape, make it simpler. \n");
      continue;
    }
    if (shapes.First().ShapeType() != TopAbs_EDGE  && 
        shapes.First().ShapeType() != TopAbs_COMPOUND)
      continue;

    TopoDS_Shape result = shapes.First();
    TopoDS_Edge common_edge;
    if(edge == NULL && result.ShapeType() == TopAbs_EDGE)
      common_edge = TopoDS::Edge(result);
    
    else if(edge == NULL && 
            result.ShapeType() == TopAbs_COMPOUND)
    {
      TopExp_Explorer Ex(result, TopAbs_EDGE);
      if(Ex.More())
        common_edge = TopoDS::Edge(Ex.Current());
    }  
    if (common_edge.IsNull())
      continue;
    BRepBuilderAPI_Copy api_copy(common_edge);
    TopoDS_Shape newShape = api_copy.ModifiedShape(common_edge);
    edge = new TopoDS_Edge(TopoDS::Edge(newShape));

    from_faces.append(new TopoDS_Face(face));
  }
  return edge;
}

int OCCModifyEngine::check_intersection(DLIList<TopoDS_Edge*>*& edge_list,
 				        TopoDS_Face from_face)const
{
  int  count_intersection = 0;

  gp_Pnt newP[2] , p_test ;
  for(int kk = 0; kk < edge_list->size(); kk++)
  {
    TopoDS_Edge* edge = edge_list->get_and_step();
    BRepAdaptor_Curve acurve(*edge);
    double lower_bound = acurve.FirstParameter();
    double upper_bound = acurve.LastParameter();
    TopExp_Explorer Ex;
    for (Ex.Init(from_face, TopAbs_EDGE); Ex.More(); Ex.Next())
    {
      TopoDS_Edge from_edge = TopoDS::Edge(Ex.Current());
      BRepAdaptor_Curve acurve2(from_edge);
      double lower_bound2 = acurve2.FirstParameter();
      double upper_bound2 = acurve2.LastParameter();
      BRepExtrema_DistShapeShape distShapeShape(*edge, from_edge);

      if (distShapeShape.IsDone() && distShapeShape.Value() < TOL)
      {
        //double check that the point is on the edges.
        double newVal;
        for(int j =1; j <= distShapeShape.NbSolution(); j++)
        {
          if(count_intersection == 2)
            break;

          p_test = distShapeShape.PointOnShape1(j);
          Extrema_ExtPC ext(p_test, acurve, Precision::Approximation());
          // At this time, there must be a intersection point at least.
          if (ext.IsDone() && (ext.NbExt() > 0)) {
            for ( int i = 1 ; i <= ext.NbExt() ; i++ ) {
              newVal = ext.Point(i).Parameter();
              if ((newVal-lower_bound) >= -TOL &&
                  (upper_bound - newVal) >= -TOL)
              {
                Extrema_ExtPC ext(p_test, acurve2, Precision::Approximation());
                if (ext.IsDone() && (ext.NbExt() > 0)) {
                  for ( int k = 1 ; k <= ext.NbExt() ; k++ ) {
                    newVal = ext.Point(i).Parameter();
                    if ((newVal-lower_bound2) >= -TOL &&
                        (upper_bound2 - newVal) >= -TOL)
                    {
                      newP[count_intersection] = p_test;
                      count_intersection ++;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }

      if (count_intersection == 2)
      {
         //make sure the two intersect point are not the same one
         if (newP[0].IsEqual(newP[1], TOL))
           count_intersection--;
      }
      if (count_intersection == 2)
        break;
    } //for loop
  }
  return count_intersection;
}
//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on OCC-based bodies
// Author     : Jane Hu
// Date       : 04/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(BodySM* BodyPtr1, BodySM* BodyPtr2,
                                         BodySM*& newBody1, BodySM*& newBody2,
                                         bool  keep_old) const
{
  newBody1 = NULL;
  newBody2 = NULL;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_volume;
  
  DLIList<BodySM*> bodysm_list;
  bodysm_list.append(BodyPtr1);
  bodysm_list.append(BodyPtr2);
  
  CubitStatus stat = get_shape_list(bodysm_list,shape_list,is_volume,keep_old);

  if(!stat)
    return stat;

  TopoDS_Shape* shape1 = shape_list.get();
  TopoDS_Shape* shape2 = shape_list.step_and_get();
  DLIList<TopologyBridge*> tbs;
  DLIList<TopoDS_Face*> face_list;
  int result = imprint_toposhapes(shape1, shape2, face_list);
  if(result == 0)
  {
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape1); 
    newBody1 = CAST_TO(tbs.get(),BodySM);
  }

  else if(result == 1)
  {
    //for cylinder overlapping cases, doing boolean operation is necessary for 
    //now. 10/25/11
    BodySM* tool_copy = copy_body(BodyPtr2);
/*
    DLIList<BodySM*> tools;
    tools.append(tool_copy);
    DLIList<BodySM*> from_bodies ;
    BodySM* from_copy = copy_body(BodyPtr1);
    from_bodies.append(from_copy);
    DLIList<BodySM*> new_bodies;
    stat = subtract(tools, from_bodies, new_bodies, false, false);
    if (stat)
    {
      tool_copy = copy_body(BodyPtr2);
      from_bodies.clean_out();
      from_bodies.append(BodyPtr1);
      intersect(tool_copy, from_bodies, new_bodies, false);
      if(new_bodies.size() > 1)
      {
        DLIList<BodySM*> final_bodies;
        unite(new_bodies, final_bodies, false); 
        newBody1 = final_bodies.get();
      }
      else
        newBody1 = new_bodies.get();
    }
*/
    stat = result_1_imprint(BodyPtr1, tool_copy, newBody1);
    if(stat == CUBIT_FAILURE)
      result = 2;
  }
  else if(result == 3)
  {
    TopoDS_Face* face = face_list.pop();
    int j = OCCQueryEngine::instance()->OCCMap->Find(*face);
    OCCSurface* cut_face = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(j))->second; 
    OCCBody* tool_body = cut_face->my_body();
    stat = result_3_imprint(BodyPtr1,tool_body, newBody1); 
    if(stat == CUBIT_FAILURE)
      result = 2;
  }

  if(result && keep_old)
  {
    delete shape1;
    shape1 = NULL;
    PRINT_INFO("There's no imprint on the first body.\n");
    newBody1 = BodyPtr1;
  }

  tbs.clean_out();
  result = imprint_toposhapes(shape2, shape1, face_list);
  if(result == 0)
  {
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape2);
    newBody2 = CAST_TO(tbs.get(),BodySM);     
  }
  
  else if(result == 1)
  {
    //for cylinder overlapping cases, doing boolean operation is necessary for 
    //now. 10/25/11
    BodySM* tool_copy = copy_body(BodyPtr1);
/*
    DLIList<BodySM*> tools;
    tools.append(tool_copy);
    DLIList<BodySM*> from_bodies ;
    BodySM* from_copy = copy_body(BodyPtr2);
    from_bodies.append(from_copy);
    DLIList<BodySM*> new_bodies;
    stat = subtract(tools, from_bodies, new_bodies, false, false);
    if (stat)
    {
      tool_copy = copy_body(BodyPtr1);
      from_bodies.clean_out();
      from_bodies.append(BodyPtr2);
      intersect(tool_copy, from_bodies, new_bodies, false);
      if(new_bodies.size() > 1)
      {
        DLIList<BodySM*> final_bodies;
        unite(new_bodies, final_bodies, false);       
        newBody2 = final_bodies.get();
      }
      else
        newBody2 = new_bodies.get();
    }
*/
    stat = result_1_imprint(BodyPtr2, tool_copy, newBody2);
    if(stat == CUBIT_FAILURE)
      result = 2; 
  }

  else if(result == 3)
  {
    TopoDS_Face* face = face_list.pop();
    int j = OCCQueryEngine::instance()->OCCMap->Find(*face);
    OCCSurface* cut_face = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(j))->second; 
    OCCBody* tool_body = cut_face->my_body();
    stat = result_3_imprint(BodyPtr2,tool_body, newBody2); 
    if(stat == CUBIT_FAILURE)
      result = 2;
  }

  if(result && keep_old)
  {
    delete shape2;
    shape2 = NULL;
    PRINT_INFO("There's no imprint on the second body.\n");
    newBody2 = BodyPtr2;
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_shape_list
// Member Type: PRIVATE
// Description: get the TopoDS_Shape list for imprinting use. 
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus OCCModifyEngine::get_shape_list(DLIList<BodySM*>& BodySM_list, 
                                         DLIList<TopoDS_Shape*>& shape_list,
                                         DLIList<CubitBoolean>& is_volume,
                                         bool  keep_old,
                                         DLIList<CubitBox*>* b_boxes) const
{
  OCCBody* occ_body = NULL;
  shape_list.clean_out();
  is_volume.clean_out();
  CubitStatus stat = CUBIT_SUCCESS;
  for(int i = 0; i <BodySM_list.size(); i++)
  {
    occ_body = CAST_TO(BodySM_list.get_and_step(), OCCBody);
    if (!occ_body)
      continue;

    if(b_boxes)
    {
      CubitBox *tool_box = new CubitBox(occ_body->get_bounding_box());
      b_boxes->append(tool_box);
    }

    TopoDS_Compound* shape = occ_body->get_TopoDS_Shape();
    if( shape && !shape->IsNull())
    {
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*shape);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*shape);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append_unique(shape);
      TopExp_Explorer Ex(*shape, TopAbs_SOLID);
      if(Ex.More())
        is_volume.append( CUBIT_TRUE);
      else
        is_volume.append( CUBIT_FALSE);
      continue;
    }

    DLIList<OCCSurface*> surfaces;
    DLIList<OCCShell*>   shells;
    DLIList<OCCCurve*>   curves;    

    surfaces = occ_body->my_sheet_surfaces();
    shells = occ_body->shells();
    if(surfaces.size() + shells.size() > 1)
    {
      PRINT_ERROR("Can't do boolean operation on multiple-volume body.\n");
      stat = CUBIT_FAILURE;
      break;
    }
    is_volume.append( CUBIT_TRUE);

    if(surfaces.size() == 1)
    {
      TopoDS_Face* topo_face = surfaces.get()->get_TopoDS_Face();
      if(!topo_face)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*topo_face);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_face);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(topo_face);
      is_volume.last();
      is_volume.change_to( CUBIT_FALSE);
    }
    else if(shells.size() == 1)
    {
      TopoDS_Shell* topo_shell = shells.get()->get_TopoDS_Shell();
      if(!topo_shell)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*topo_shell);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_shell);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(topo_shell);
      is_volume.change_to( CUBIT_FALSE);
    }

    else 
    {
      DLIList<Lump*> lumps = occ_body->lumps();
      if (lumps.size() > 1)
      {
        PRINT_ERROR("Can't do boolean operation on multi-volumes types. \n");
        stat = CUBIT_FAILURE;
        break;
      }

      TopoDS_Solid* solid = CAST_TO(lumps.get(), OCCLump)->get_TopoDS_Solid();
      if(!solid)
      {
        stat = CUBIT_FAILURE;
        break;
      }
      if(keep_old)
      {
        BRepBuilderAPI_Copy api_copy(*solid);
        TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
        TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
        shape_list.append(Shape1);
      }
      else
        shape_list.append(solid);
    }
  }
  if(!stat)
  {   
    for (int i = 0; keep_old && i < shape_list.size(); i++)
    {
          TopoDS_Shape* shape = shape_list.get_and_step();
          delete shape;
          shape = NULL;
    }
    shape_list.clean_out();
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : imprint multiple bodies at once
// Member Type: PUBLIC
// Description: imprint boolean operation on OCC-based bodies
// Author     : Jane HU
// Date       : 04/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint(DLIList<BodySM*> &from_body_list ,
                                     DLIList<BodySM*> &new_from_body_list,
                                     bool keep_old,
                                     DLIList<TopologyBridge*>* new_tbs,
                                     DLIList<TopologyBridge*>* att_tbs) const
{
  CubitStatus success = CUBIT_SUCCESS;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_vo;

  //keep record of all vertices and edges and faces in the from_body_list,
  //for comparison to generated new_tbs and att_tbs.
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCCurve*> curves;
  DLIList<OCCPoint*> points;

  CubitStatus stat = get_shape_list(from_body_list, shape_list, is_vo,keep_old);

  if(!stat)
    return stat;
 
  int size = shape_list.size();
  // total number of imprints to be done
  int total_imprints = size * (size -1);

  if( size > 2 )
  {
     char message[128];
     sprintf(message, "Imprinting %d OCC Bodies", from_body_list.size() ); 
     AppUtil::instance()->progress_tool()->start(0, total_imprints, message);
  }
  
  std::map<OCCSurface*, std::pair<CubitVector, int> > surf_property_map;
  std::map<OCCCurve*, std::pair<CubitVector, int> > curve_property_map;

  for(int i = 0; i < size; i++)
  {
    TopoDS_Shape* shape1 = shape_list[i];
    CubitBoolean modified = CUBIT_FALSE;

    DLIList<TopoDS_Face*> face_list;
    for(int j = i+1; j < size+i; j ++)
    {
       if (AppUtil::instance()->interrupt())
       {
          success = CUBIT_FAILURE;
          break;
       }

       TopoDS_Shape* shape2 = shape_list[j%size];
       DLIList<TopologyBridge*> tbs;
       int result = imprint_toposhapes(shape1, shape2, face_list);
       if(result == 0)
       {
          DLIList<TopologyBridge*> tbs;
          tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape1);
          from_body_list[i] = CAST_TO(tbs.get(),BodySM);
          modified = CUBIT_TRUE;
       }
       else if(result == 1)
       {
         //for cylinder overlapping cases, doing boolean operation is necessary for 
         //now. 10/25/11
         BodySM* tool_copy = copy_body(from_body_list[j%size]);
         BodySM* newBody = NULL;

         stat = result_1_imprint(from_body_list[i], tool_copy, newBody);
         if (stat && newBody)
         {
           modified = CUBIT_TRUE;
           DLIList<TopoDS_Shape*> shapes;
           DLIList<CubitBoolean> is_volume;
           DLIList<BodySM*> new_bodies;
           new_bodies.append(newBody);
           get_shape_list(new_bodies, shapes, is_volume,false);
           shape1 = shapes.get();
           //shape_list[i] = shape1;
           from_body_list[i] = newBody;
         }
       }
       else if(result == 3)
       {
         BodySM* newBody = NULL;
         BodySM* oldBody = from_body_list[i];
         TopExp_Explorer Ex;
         Ex.Init(*shape1, TopAbs_SOLID);
         int nSolid = 0;
         for(; Ex.More(); Ex.Next())
           nSolid++;
         do
         { 
           TopoDS_Face* face = face_list.pop();
           int k = OCCQueryEngine::instance()->OCCMap->Find(*face);
           OCCSurface* cut_face = (OCCSurface*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
           OCCBody* tool_body = cut_face->my_body();

           stat = result_3_imprint(oldBody, tool_body, newBody);
           if (newBody)
           {
             modified = CUBIT_TRUE;
             DLIList<TopoDS_Shape*> shapes;
             DLIList<CubitBoolean> is_volume;
             DLIList<BodySM*> new_bodies;
             new_bodies.append(newBody);
             get_shape_list(new_bodies, shapes, is_volume,false);
             shape1 = shapes.get();
             //shape_list[i] = shape1;
             from_body_list[i] = newBody;
           }
           if(nSolid == 1)
             break;
           else
             nSolid --;
           result = imprint_toposhapes(shape1, shape2, face_list); 
           oldBody = newBody;
         }while (result == 3);
       }
    }
    if(modified)
      new_from_body_list.append(from_body_list[i]);
    
    shape_list.reset();
    if( size > 2 )
      AppUtil::instance()->progress_tool()->step();
  }

  if( size > 2 )
    AppUtil::instance()->progress_tool()->end();

  if( AppUtil::instance()->interrupt() )
        PRINT_INFO("Imprint aborted.\n");

  return success;
}

CubitStatus OCCModifyEngine::result_3_imprint(BodySM* from_body, 
                                              BodySM* tool_body, 
                                              BodySM*& newBody)const
{
  DLIList<BodySM*> from_bodies;
  from_bodies.append(from_body); 
  DLIList<BodySM*> results_list;
  DLIList<BodySM*> neighbor_imprint_list;
  CubitStatus stat = webcut(from_bodies, tool_body, neighbor_imprint_list, 
                            results_list);
  if(results_list.size() > 1)
  {
    DLIList<BodySM*> bodies;
    unite(results_list, bodies, false); 
    newBody = bodies.get();
  }
  else if (results_list.size() == 1)
    newBody = results_list.get();
 
  return stat;

}

CubitStatus OCCModifyEngine::result_1_imprint(BodySM* from_body,  
                                              BodySM* tool_body,
                                              BodySM*& newBody)const
{
  //for compound of solids in from_body, try to determine if subtracting or
  //intersecting will keep the from_body un-deleted by checking subtracting
  //first, if it kept the from_body, don't take any further risk, do intersect
  //first, then subtract.

  DLIList<TopoDS_Shape*> from_bodies_shapes;
  DLIList<BodySM*> from_bodies ;
  DLIList<CubitBoolean> is_volume;
  from_bodies.append(from_body);
  //get the from_bodies underling shapes
  CubitStatus stat = get_shape_list(from_bodies, from_bodies_shapes, is_volume, CUBIT_FALSE);

  TopoDS_Shape* from_shape = from_bodies_shapes.get();
  int num_solids = 0;
  if(from_shape->ShapeType() == TopAbs_COMPOUND)
  {
    TopExp_Explorer Ex;
    for (Ex.Init(*from_shape, TopAbs_SOLID);Ex.More(); Ex.Next())
      num_solids++;
  }
  BodySM* tool_copy = copy_body(tool_body);
  DLIList<BodySM*> tools;
  tools.append(tool_copy);
  from_bodies.clean_out() ;
  BodySM* from_copy = copy_body(from_body);
  from_bodies.append(from_copy);
  DLIList<BodySM*> new_bodies;
  stat = subtract(tools, from_bodies, new_bodies, false, false);
  if (stat)
  {
    if (num_solids > 1 && from_copy == new_bodies.get())
    {
      //double check if subtract kept the original from_bodies pointer.
      //if so, do intersect first, then subtract.
      tool_copy = copy_body(tool_body);
      from_bodies.clean_out(); 
      from_copy = copy_body(from_body);
      from_bodies.append(from_copy);
      OCCQueryEngine::instance()->delete_solid_model_entities(new_bodies);
      new_bodies.clean_out();
      intersect(tool_copy, from_bodies, new_bodies, false);

      tool_copy = copy_body(tool_body);
      tools.clean_out();
      tools.append(tool_copy);
      from_bodies.clean_out();
      from_bodies.append(from_body);
      subtract(tools, from_bodies, new_bodies, false, false);
    }
    else
    { 
      tool_copy = copy_body(tool_body);
      from_bodies.clean_out();
      from_bodies.append(from_body);
      intersect(tool_copy, from_bodies, new_bodies, false);
    }
    
    //make sure the first body in new_bodies is the one we want to keep.
    new_bodies.reverse();
    if(new_bodies.size() > 1)
    {
       DLIList<BodySM*> final_bodies;
       unite(new_bodies, final_bodies, false);
       new_bodies.clean_out();
       new_bodies = final_bodies;
    }
    newBody = new_bodies.get();
  }
  return stat;
} 

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint curves onto body_list
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                           DLIList<Curve*> &ref_edge_list,
                                           DLIList<BodySM*>& new_body_list,
                                           DLIList<TopologyBridge*> &temp_bridges,
                                           bool keep_old,
                                           bool show_messages) const
{
  CubitStatus success = CUBIT_SUCCESS;
  DLIList<TopoDS_Shape*> shape_list, tool_shapes;
  DLIList<CubitBoolean> is_vo;
  CubitStatus stat = get_shape_list(body_list, shape_list, is_vo, keep_old);
  if (!stat)
    return stat;

  int size = ref_edge_list.size();
  // total number of imprints to be done

  if( size > 2 && show_messages)
  {
     char message[128];
     sprintf(message, "Imprinting %d OCC Bodies", body_list.size() );
     AppUtil::instance()->progress_tool()->start(0, size, message);
  }
  for (int i = 0; i < ref_edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(ref_edge_list.get_and_step(), OCCCurve) ;
    if (!curve)
      continue;

    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
    if (edge->IsNull())
      continue;
    
    if (AppUtil::instance()->interrupt())
    {
       success = CUBIT_FAILURE;
       break;
    }
    DLIList<TopoDS_Face*> face_list;
    for(int j = 0; j < shape_list.size(); j ++)
    {
      TopoDS_Shape* shape = shape_list.get();
        
      int result = imprint_toposhapes(shape, (TopoDS_Shape*)edge, face_list);
      if (result == 0)
        shape_list.change_to(shape);
      shape_list.step();
      body_list.step();
    }

    if( size > 2 )
      AppUtil::instance()->progress_tool()->step();
  }   

  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
    new_body_list.append(CAST_TO(tbs.get(),BodySM));
  }

  if( size > 2 )
    AppUtil::instance()->progress_tool()->end();

  if( AppUtil::instance()->interrupt() )
        PRINT_INFO("Imprint aborted.\n"); 
  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: to be consistante with Acis imprint.
//              The surfaces must be part of a body, but the curves 
//              just have to be valid OCC edge.
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint( DLIList<Surface*> &ref_face_list,
                                      DLIList<Curve*> &edge_list,
                                      DLIList<TopologyBridge*> &temp_bridges,
                                      DLIList<BodySM*>& new_body_list,
                                      bool keep_old) const
{
  DLIList<TopoDS_Face*> face_list;
  DLIList<TopoDS_Shape*> shape_list;
 
  face_edge_imprint(ref_face_list, edge_list, face_list, shape_list, keep_old);

  TopExp_Explorer Ex;
  int num_face = 0;
  TopoDS_Shape* shape = shape_list.get_and_step();
  for (Ex.Init(*shape, TopAbs_FACE); Ex.More(); Ex.Next())
    { 
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      num_face++;
    }

  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    if (!shape || shape->IsNull())
      continue;
    if (shape->ShapeType() == TopAbs_COMPOUND)
    {
      if(!OCCQueryEngine::instance()->OCCMap->IsBound(*shape)) 
      {
        TopExp_Explorer Ex;
        for (Ex.Init(*shape, TopAbs_SOLID);Ex.More(); Ex.Next())
        {
          TopoDS_Shape subshape = Ex.Current();
          tbs += OCCQueryEngine::instance()->populate_topology_bridge(subshape);
          new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
        }
      }
    }
    else
    {
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
      new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
    }
  }

/*
  if (keep_old)
  {
    for(int i = 0; i < face_list.size(); i++)
    {
      TopoDS_Face* face = face_list.get();
      face->Nullify();
      delete face;
      face =NULL;
    }
  }
*/
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : face_edge_imprint
// Member Type: PRIVATE
// Description: to be consistante with Acis imprint.
//              The surfaces must be part of a body, but the curves
//              just have to be valid OCC edge.
// Author     : Jane Hu
// Date       : 05/08
//===============================================================================
CubitStatus 
OCCModifyEngine::face_edge_imprint( DLIList<Surface*> &ref_face_list,
                                    DLIList<Curve*> &edge_list,
                                    DLIList<TopoDS_Face*>& face_list,
                                    DLIList<TopoDS_Shape*>& shape_list,
                                    bool keep_old ) const
{
  for(int i = 0; i <ref_face_list.size(); i++)
  {
    OCCSurface* surface = CAST_TO(ref_face_list.get_and_step(), OCCSurface);
    if(!surface)
      continue;

    TopoDS_Face* topo_face = surface->get_TopoDS_Face();
    face_list.append(topo_face);

    OCCBody* occ_body = NULL;
    OCCShell* shell = surface->my_shell();
    if(shell && shell->my_body())
      occ_body = shell->my_body();
    else
    {
      DLIList<OCCBody*> bodies;
      surface->get_bodies(bodies);
      if(bodies.size() != 1)
      {
        PRINT_ERROR("Can't find the corresponding manifold solid body.\n");
        return CUBIT_FAILURE;
      }
      occ_body = bodies.get();
    }
    TopoDS_Shape *shape ; 
    occ_body->get_TopoDS_Shape(shape);

    if(shape && !shape->IsNull())
      shape_list.append_unique(shape);
  }

  if(keep_old)
  {
    for(int i = 0; i < shape_list.size(); i++)
    {
      TopoDS_Shape* shape = shape_list.get();
      BRepBuilderAPI_Copy api_copy(*shape);
      TopoDS_Shape newShape = api_copy.ModifiedShape(*shape);
      TopoDS_Shape* Shape1 = new TopoDS_Shape(newShape);
      for(int j = 0; j < face_list.size(); j++)
      {
        TopoDS_Face* face = face_list.get();
        TopExp_Explorer Ex, Ex_new;
        for (Ex.Init(*shape, TopAbs_FACE); Ex.More(); Ex.Next())
        {
          if(face->IsSame(Ex.Current()))
          {
            face = new TopoDS_Face(TopoDS::Face(api_copy.ModifiedShape(*face)));
            face_list.change_to(face);
          }
        }
        face_list.step();
      }
      shape_list.change_to(Shape1);
      shape_list.step();
    }
  }

  for (int i = 0; i < edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(edge_list.get_and_step(), OCCCurve) ;
    if (!curve)
      continue;

    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
    if (edge->IsNull())
      continue;

    for(int j = 0; j < shape_list.size(); j ++)
    {
      TopoDS_Shape* shape = shape_list.get_and_step();
      DLIList<TopoDS_Face*> record_faces;
      for (int e = 0; e < face_list.size(); e++)
        record_faces.append(face_list.get_and_step()); 

      imprint_toposhapes(shape, (TopoDS_Shape*)edge, face_list);
      for (int e = 0; e < record_faces.size()&& keep_old; e++)  
      {
        TopoDS_Face* test_face = record_faces.get_and_step();
        if(!face_list.move_to(test_face))
        {
          test_face->Nullify();
          delete test_face;
        }
      }
    }
  }

  for(int j = 0; keep_old && j < face_list.size(); j++)
  {
    TopoDS_Face* face = face_list.get_and_step(); 
    face->Nullify();
    delete face;
  }

  TopExp_Explorer Ex;
  int num_face = 0;
  TopoDS_Shape* shape = shape_list.get_and_step();
  for (Ex.Init(*shape, TopAbs_FACE); Ex.More(); Ex.Next())
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      num_face++;
    }

  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: To be consistent with  AcisModifyEngine, althought it's hard 
//              to have a GUI interface for users to input. All surface must
//              on the same body. 
// Author     : Jane HU 
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint( DLIList<Surface*>& surface_list,
                                   DLIList<DLIList<Curve*>*>& curve_lists_list,
                                   BodySM*& new_body,
                                   bool keep_old,
                                   bool expand,
                                   DLIList<TopologyBridge*> *new_tbs,
                                   DLIList<TopologyBridge*> *att_tbs  ) const
{
  DLIList<TopoDS_Face*> face_list;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<TopoDS_Shape*> shape_list_all;
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCCurve*> curves;
  DLIList<OCCPoint*> points;
 
  assert (surface_list.size() == curve_lists_list.size());
  DLIList<OCCBody*> bodies;

  std::map<OCCSurface*, std::pair<CubitVector, double> > surf_property_map;
  std::map<OCCCurve*, std::pair<CubitVector, double> > curve_property_map;

  for(int j = 0; j < surface_list.size(); j++)
  {
    Surface* surface = surface_list.get_and_step();
    
    //keep record of old bodies, surfaces, curves and points
    OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
    occ_surface->get_bodies(bodies);
    if(j == 0)
    {
      bodies.get()->get_all_surfaces(surfaces);
      bodies.get()->get_all_curves(curves);
      bodies.get()->get_all_points(points);
 
      //save surface with its area and center info in the map
      for(int i = 0; i < surfaces.size(); i++)
      {
        OCCSurface* surf = CAST_TO(surfaces.get(), OCCSurface);
        double d = surf->measure();
        CubitVector center = surf->center_point();
        surf_property_map.insert(valType
             (surf, std::pair<CubitVector, double>(center,d)));
      }

      //save curve with its length and center info in the map
      for(int i = 0; i < curves.size(); i++)
      {
        OCCCurve* curve = CAST_TO(curves.get(), OCCCurve);
        double d = curve->measure();
        CubitVector center = curve->center_point();
        curve_property_map.insert(valType2
             (curve, std::pair<CubitVector, double>(center,d)));
      }
    }
    DLIList<Surface*> ref_face_list;
    ref_face_list.append(surface);
    DLIList<Curve*> *edge_list = curve_lists_list.get_and_step();
    face_edge_imprint(ref_face_list, *edge_list, face_list, shape_list, keep_old);

    for(int i = 0; i < shape_list.size(); i++)
    {
      TopoDS_Shape* shape = shape_list.get_and_step();
      shape_list_all.append_unique(shape);
    }
    shape_list.clean_out();

    if (keep_old)
    {
      for(int i = 0; i < face_list.size(); i++)
      {
        TopoDS_Face* face = face_list.get();
        face->Nullify();
        delete face;
        face = NULL;
      }
    }

    face_list.clean_out();
  }

  assert (bodies.size() == 1);

  DLIList<BodySM*> new_body_list;
  shape_to_bodySM(shape_list, new_body_list);
  
  if (new_body_list.size() == 1)
  {
    new_body = new_body_list.get();
    //find new_tbs and att_tbs;
    DLIList<OCCSurface*> new_surfs;
    DLIList<OCCCurve*> new_curves;
    DLIList<OCCPoint*> new_points;
    if(new_tbs || att_tbs) 
    {
      OCCBody* occ_body = CAST_TO(new_body, OCCBody);
      occ_body->get_all_surfaces(new_surfs);
      occ_body->get_all_curves(new_curves);
      occ_body->get_all_points(new_points);
    } 
    if(new_tbs)
      get_new_tbs(surf_property_map, curve_property_map, points, new_surfs, 
                  new_curves, new_points, new_tbs);
    if(att_tbs)
      get_att_tbs(new_surfs, new_curves, new_points, "COMPOSITE_GEOM",
                  att_tbs);

    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : shape_to_bodySM
// Member Type: PRIVATE
// Description: After imprint, update shape list to bodySM_list
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
void OCCModifyEngine::shape_to_bodySM( DLIList<TopoDS_Shape*> shape_list,
                                       DLIList<BodySM*>& new_body_list)const
{
  for(int j = 0; j < shape_list.size(); j ++)
  {
    DLIList<TopologyBridge*> tbs;
    TopoDS_Shape* shape = shape_list.get_and_step();
    if (shape->ShapeType() == TopAbs_COMPOUND)
    {
      if(!OCCQueryEngine::instance()->OCCMap->IsBound(*shape))
      {
        TopExp_Explorer Ex;
        for (Ex.Init(*shape, TopAbs_SOLID);Ex.More(); Ex.Next())
        {
          TopoDS_Shape subshape = Ex.Current();
          tbs += OCCQueryEngine::instance()->populate_topology_bridge(subshape);
          new_body_list.append_unique(CAST_TO(tbs.get(),BodySM));
        }
      }
    }
    else
    {
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape);
      BodySM* body = CAST_TO(tbs.get(),BodySM);
      if(body != NULL)
        new_body_list.append_unique(body);
    }
  }
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: Imprints locations to bodies (for splitting curves, there's
//              no known ways to put hard points on surfaces in OCC, so I just
//              add free_vertex on OCCSurface definition, mesh should look on
//              this structure).   
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                          DLIList<CubitVector> &vector_list,
                                          DLIList<BodySM*>& new_body_list,
                                          bool keep_old,
                                          DLIList<TopologyBridge*>* new_tbs,
                                          DLIList<TopologyBridge*>* att_tbs,
                                          double *tol_in ,
                                          bool clean_up_slivers ) const
{
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_vo;
  double tol = 0.1;
  if(tol_in)
    tol = *tol_in;
  CubitStatus stat = get_shape_list(body_list, shape_list, is_vo, keep_old);
  if(!stat)
    return stat;

  for (int i = 0; i < body_list.size(); i++)
  {
    OCCBody* body = CAST_TO(body_list.get_and_step(), OCCBody);
    TopoDS_Shape* from_shape = shape_list.get();
    if (!body)
      continue;
    DLIList<OCCSurface*> surfaces;
    DLIList<OCCCurve*> curves;

    body->get_all_surfaces(surfaces);
    body->get_all_curves(curves);
    
    CubitBoolean on_vertex = CUBIT_FALSE;
    CubitBoolean on_curve = CUBIT_FALSE;
    for (int j = 0; j < vector_list.size(); j ++)
    {
      CubitVector& v = vector_list[j];
      for (int k = 0;  k < curves.size(); k ++)
      {
         OCCCurve* curve = curves.get_and_step();
         CubitPointContainment pc = curve->point_containment(v);
         if(pc == CUBIT_PNT_BOUNDARY)
         {
           on_vertex = CUBIT_TRUE;
           break;
         }

         else if( pc == CUBIT_PNT_INSIDE)
         {
           on_curve = CUBIT_TRUE;
           //first make sure it won't generate a sliver curve
           //with respect to tol.
           if(clean_up_slivers)
           {
             double u = curve->u_from_position(v);
             double u_min, u_max;
             curve->get_param_range(u_min, u_max);
             double l1 = curve->length_from_u(u_min, u);
             double l2 = curve->measure() - l1;
             if(l1 <= tol || l2 <= tol)
               break;
           }

           const CubitVector location = v;
           DLIList<Curve*> created_curves;
           stat = split_shape_by_location(from_shape, (Curve*)curve, location, 
                                          created_curves);
           shape_list.change_to(from_shape);
           if(new_tbs)
             for(int ic=0; ic < created_curves.size(); ic++)
               new_tbs->append(created_curves.get_and_step());

           curves.remove(curve);
           for(int ic = 0; ic < created_curves.size(); ic++)
             curves.append(CAST_TO(created_curves.get_and_step(), OCCCurve));
           break;
         }  
       } 
       if(on_vertex || on_curve)
         continue;

       //check possible on surface
       for(int n = 0; n < surfaces.size(); n ++)
       {
          OCCSurface* surface = surfaces.get_and_step();
          if(!surface->is_position_on(v))
            continue;
           
          CubitPointContainment ps = surface->point_containment(v);
          if(ps == CUBIT_PNT_INSIDE)
          {
             TBPoint* p = make_Point(v);
             if(p)
             {
               surface->add_hardpoint(CAST_TO(p, OCCPoint));
               if(new_tbs)
                 new_tbs->append(p);
             }
             break;
          }
       }
    }
    shape_list.step();
  }       

  shape_to_bodySM(shape_list, new_body_list);
  
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCCurve*>   curves;
  DLIList<OCCPoint*>   points;
  for (int i = 0; i < new_body_list.size(); i++)
  {
    OCCBody * body = CAST_TO(new_body_list.get_and_step(), OCCBody);
    body->get_all_surfaces(surfaces);
    body->get_all_curves(curves);
    body->get_all_points(points);
  }

  get_att_tbs(surfaces, curves, points, "COMPOSITE_GEOM", att_tbs);

  return stat;
}


//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: Projects a list of Curves on to a list of Surfaces
//              and imprint the faces with the new Curves
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     
OCCModifyEngine::imprint_projected_edges( DLIList<Surface*> &ref_face_list,
                                          DLIList<Curve*> &ref_edge_list,
                                          DLIList<BodySM*>& new_body_list,
                                          DLIList<Curve*>& kept_free_edges,
                                          bool keep_old_body,
                                          bool keep_free_edges) const
{
  DLIList<Curve*> projected_curves;
  CubitStatus 
     stat = project_edges(ref_face_list, ref_edge_list, projected_curves);
  if(!stat)
    return stat;

  // imprint Surface with curves
  DLIList<TopologyBridge*> temp_bridges;
  stat = imprint(ref_face_list, projected_curves, temp_bridges, new_body_list, keep_old_body );

  if(keep_free_edges)
  {
     kept_free_edges += projected_curves;
     return  stat;
  }

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< projected_curves.size();i++)
  {
    // Now delete this Curve and its underlying solid model entities

    Curve* curve = projected_curves.get_and_step();
    stat = OCCQueryEngine::instance()->delete_solid_model_entities( curve );
    if (stat == CUBIT_FAILURE)
    {
       PRINT_ERROR("In OCCQueryEngine::delete_geometry\n"
                 "       Could not delete OCCCurve.\n"
                 "       The Model database is likely corrupted "
                 "due to\n       this unsuccessful deletion.\n" );
    }
  } 
  return stat;
}
//===============================================================================
// Function   : project_edges
// Member Type: PUBLIC
// Description: Projects a list of Curves on to a list of Surfaces
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::project_edges( DLIList<Surface*> &ref_face_list,
                                            DLIList<Curve*> &ref_edge_list,
                                            DLIList<Curve*> &projected_curves,
                                            bool print_error ) const

{
  CubitVector* v = NULL;
  Curve* projected_curve = NULL;
  DLIList<TBPoint*> points;
  //project curves onto surfaces.
  for(int i = 0; i < ref_edge_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(ref_edge_list.get_and_step(), OCCCurve);
    if(!curve)
       continue;

    for (int j = 0; j < ref_face_list.size(); j++)
    {
      OCCSurface* surface = CAST_TO(ref_face_list.get_and_step(), OCCSurface); 
      if(!surface)
        continue;
      if(surface->is_closed_in_U() || surface->is_closed_in_V())
      {
        if(print_error)
          PRINT_ERROR("This function can't project curves on closed surfaces.\n");
        return CUBIT_FAILURE;
      }
      
      projected_curve = NULL;
      projected_curve = curve->project_curve(surface, points, CUBIT_FALSE, v);
      if(projected_curve)
        projected_curves.append_unique(projected_curve);
    }
  }
  while(points.size() > 0)
    delete points.pop();
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: Projects a list of curves on to a list of surfaces
//              and imprint the bodies with the new curves
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus 
OCCModifyEngine::imprint_projected_edges(DLIList<Surface*> &ref_face_list,
                                         DLIList<BodySM*> &body_list,
                                         DLIList<Curve*> &ref_edge_list,
                                         DLIList<BodySM*>& new_body_list,
                                         bool keep_old,
                                         bool keep_free_edges) const
{
  DLIList<Curve*> projected_curves;
  CubitStatus
     stat = project_edges(ref_face_list, ref_edge_list, projected_curves);
  if(!stat)
    return stat; 
  return CUBIT_FAILURE;

  // imprint bodies with curves
  DLIList<TopologyBridge*> temp_bridges;
  stat = imprint(body_list,projected_curves, new_body_list, temp_bridges, keep_old);

  if (keep_free_edges)
        return  stat;

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< projected_curves.size();i++)
  {
    // Now delete this Curve 
    Curve* curve = projected_curves.get_and_step();
    stat = OCCQueryEngine::instance()->
          delete_solid_model_entities( curve );
    if (stat == CUBIT_FAILURE)
    {
       PRINT_ERROR("In OCCModifyEngine::delete_geometry\n"
                   "       Could not delete Curve.\n"
                   "       The Model database is likely corrupted "
                   "due to\n       this unsuccessful deletion.\n" );
    }
  }
  return stat; 
}

//===============================================================================
// Function   : intersect
// Member Type: PUBLIC
// Description: intersect boolean operation of body with list of bodies.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::intersect(BodySM*  tool_body_ptr,
                                       DLIList<BodySM*>  &from_bodies,
                                       DLIList<BodySM*>  &new_bodies,
                                       bool  keep_old,
                                       bool preview) const
{
  DLIList<BodySM*> tool_bodies;
  DLIList<TopoDS_Shape*> tool_shapes;
  DLIList<CubitBoolean> is_tool_volume, is_volume;
  
  tool_bodies.append(tool_body_ptr);
  //get tool_body's underlying shape, copy it, so boolean wouldn't touch it.
  CubitStatus stat = 
       get_shape_list(tool_bodies, tool_shapes, is_tool_volume, CUBIT_TRUE); 
  if(!stat)
    return stat;

  DLIList<TopoDS_Shape*> shape_list;
  stat =  get_shape_list(from_bodies, shape_list, is_volume, keep_old);
  if(!stat)
  {
    for (int i = 0; i < tool_shapes.size(); i++)
    {
       TopoDS_Shape* shape = tool_shapes.get_and_step();
       delete shape;
       shape = NULL; 
    }
    tool_shapes.clean_out();
    return CUBIT_FAILURE;
  }

  TopoDS_Shape* tool_shape = tool_shapes.get();
  CubitBoolean has_changed;
  DLIList<TopologyBridge*> tbs;
  DLIList<TopoDS_Shape*> preview_list;
  for (int i = 0; i < shape_list.size(); i++)
  { 
    TopoDS_Shape* from_shape = shape_list.get_and_step();
    BodySM* from_body = from_bodies.get_and_step();
    BRepAlgoAPI_Common intersector(*from_shape, *tool_shape);

    TopTools_ListOfShape shapes;
    shapes.Assign(intersector.Modified(*tool_shape));
    
    if ( shapes.Extent() > 1)
    {
      PRINT_ERROR("Tool has multiple intersection with the shape, make it simpler. \n");
      continue;
    }
    TopoDS_Shape common_shape;
    if (shapes.IsEmpty())
      common_shape = intersector.Shape();
    else 
      common_shape = shapes.First();

    if (!common_shape.IsNull())
    {
      TopExp_Explorer Ex(common_shape, TopAbs_SOLID);
      if(!Ex.More() && is_volume[i] == CUBIT_TRUE && 
         common_shape.ShapeType() > TopAbs_SOLID)
      {
        has_changed = CUBIT_TRUE;
        *from_shape = common_shape;
      }
      
      else if (is_volume[i] == CUBIT_FALSE && 
               common_shape.ShapeType() > TopAbs_FACE) 
      {
        Ex.Init(common_shape, TopAbs_SHELL, TopAbs_SOLID);
        if(!Ex.More())
          Ex.Init(common_shape, TopAbs_FACE);
        if(!Ex.More() && is_volume[i] == CUBIT_FALSE)
        {
          has_changed = CUBIT_TRUE;
          *from_shape = common_shape;
        }
      }
      else
      {
        double after_mass = 0.0;
        GProp_GProps myProps;
        if(is_volume[i])
          BRepGProp::VolumeProperties(common_shape, myProps);
        
        else
          BRepGProp::SurfaceProperties(common_shape, myProps);
        after_mass = myProps.Mass();
        if(after_mass > TOL)
          check_operation(common_shape, from_shape, is_volume[i], has_changed,
                    &intersector, keep_old);
        else
        {
          if(!keep_old)
            OCCQueryEngine::instance()->delete_solid_model_entities(from_body); 
          from_shape = NULL;
        }
      }
    }

    if(from_shape == NULL || from_shape->IsNull() )
    {
      PRINT_INFO("The %d body did not have common part with the tool_body.\n", i+1);
    }
    else
    {
      if(preview)
      {
        TopoDS_Shape* p_shape = new TopoDS_Shape(common_shape);
        preview_list.append(p_shape);
      }
      else
        tbs += OCCQueryEngine::instance()->populate_topology_bridge(*from_shape);
    }
  }

  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < preview_list.size(); i++)
    {
      TopoDS_Shape* new_shape = preview_list.get_and_step();
      // Draw this face
      OCCDrawTool::instance()->draw_TopoDS_Shape( new_shape, CUBIT_BLUE_INDEX,
                                                  CUBIT_FALSE, CUBIT_TRUE );
      delete new_shape;
    }
  }
  
  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }
  
  //if(tbs.size() == 0)
  //  stat = CUBIT_FAILURE;
    
  //ok, we're done with all cuts, delete unnecessaries.
  if(!keep_old)
    OCCQueryEngine::instance()->delete_solid_model_entities(tool_body_ptr);   

  for(int i = 0; i < tool_shapes.size(); i++)
  {
    TopoDS_Shape* shape = tool_shapes.get_and_step();
    shape->Nullify();
    delete shape;
    shape = NULL;
  }

  if(keep_old)
  {
    int size  = shape_list.size();
    for (int i = 0; i < size; i++)
    {
      TopoDS_Shape* shape = shape_list.pop();
      shape->Nullify();
      delete shape;
      shape = NULL;
    }
  }
  if(!stat)
    return stat;
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : check_operation
// Member Type: PRIVATE
// Description: check and update the from_shape according to type of the body.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
void OCCModifyEngine::check_operation(TopoDS_Shape& cut_shape,
                                      TopoDS_Shape*& from_shape, //output
                                      CubitBoolean  is_volume,
                                      CubitBoolean& has_changed, //output
                                      BRepAlgoAPI_BooleanOperation* op,
                                      CubitBoolean keep_old) const
{
   //compare to see if the from_shape has gotten cut.
   if(is_volume)
   {
     GProp_GProps myProps;
     BRepGProp::VolumeProperties(*from_shape, myProps);
     double orig_mass = myProps.Mass();
     TopTools_IndexedMapOfShape M;
     TopExp::MapShapes(cut_shape, TopAbs_SOLID, M);
     double after_mass = 0.0;
     CubitBoolean no_volume = CUBIT_FALSE;
     if(M.Extent() > 0)
     {
       BRepGProp::VolumeProperties(cut_shape, myProps);
       after_mass = myProps.Mass();
     }
     else
       no_volume = CUBIT_TRUE;
 
     if(fabs(-after_mass + orig_mass) <= TOL)
     {
        has_changed= CUBIT_FALSE; //common is itself
        return;
     }

     //got cut. Update the entities
     if(after_mass < TOL || no_volume) //no common section
       cut_shape.Nullify();
     has_changed = CUBIT_TRUE;
     TopExp_Explorer Ex;
     int num_solid = 0;
     Ex.Init(*from_shape, TopAbs_SOLID);
     TopoDS_Solid old_solid;
     for(; Ex.More(); Ex.Next())
     {
       num_solid ++;
       old_solid = TopoDS::Solid(Ex.Current());
     }
     if(num_solid == 1)
       OCCLump::update_OCC_entity(old_solid , cut_shape, op);
     
     else if(num_solid > 1)
       OCCBody::update_OCC_entity(*from_shape, cut_shape, op);
   }
   else
   {
     GProp_GProps myProps;
     BRepGProp::SurfaceProperties(*from_shape, myProps);
     double orig_mass = myProps.Mass();
     BRepGProp::SurfaceProperties(cut_shape, myProps);
     double after_mass = myProps.Mass();
     if(fabs(-after_mass + orig_mass) <= TOL)
     {
       has_changed= CUBIT_FALSE; //common is itself, or not cut
       return;
     }
     //got cut. Update the entities
     if(after_mass < TOL)//no common section
       cut_shape.Nullify();
     has_changed = CUBIT_TRUE;
     if(from_shape->ShapeType() == TopAbs_SHELL)
     {
       TopoDS_Shell old_shell = TopoDS::Shell(*from_shape);
       OCCShell::update_OCC_entity(old_shell,cut_shape, op);
     }
     else
     {
       TopoDS_Face old_face = TopoDS::Face(*from_shape);
       OCCSurface::update_OCC_entity(old_face,cut_shape, op);
     }
  }
//  if(keep_old) - Must not be deleted, causes random failure on OSX
//    delete from_shape;
  from_shape = new TopoDS_Shape(cut_shape);
}

//===============================================================================
// Function   : chop
// Member Type: PUBLIC
// Description: chop boolean operation between OCC-based bodies
//              bodies has a size() = 2, a blank body and a tool body.
//              chops the blank with the  tool, returing the body formed
//              by subtracting the tool from the blank, and the body formed
//              by intersecting the tool with the blank, simultaneously.
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus  OCCModifyEngine::chop(DLIList<BodySM*>& bodies, 
                                   DLIList<BodySM*> &intersectBodies, 
                                   DLIList<BodySM*> &outsideBodies,
                                   BodySM*& leftoversBody,
                                   bool keep_old ,
                                   bool nonreg) const
{
  //according to Acis chop function, leftoverBody = 0;
  leftoversBody = 0;

  //there's no effect of nonreg. keep_old mean if to keep the tool_body
  if(bodies.size() != 2)
  {
    PRINT_WARNING("Chop operation works only on two bodies. Nothing modified\n");  
    return CUBIT_FAILURE; 
  }
  
  //outsideBodies keeps the surface, curve ids if keep_old is false.
  BodySM* blank_body = bodies.get();
  
  //copy blank_body for intersect operation, because it will get changed.
  DLIList<BodySM*> tool_bodies, from_bodies;
  from_bodies.append(blank_body);
  BodySM* tool_body = bodies.step_and_get();
  tool_bodies.append(tool_body);
  
  CubitStatus stat = intersect(tool_body, from_bodies, 
                               intersectBodies, CUBIT_TRUE);

  if(!stat)
    return CUBIT_FAILURE;

  stat = subtract(tool_bodies, from_bodies, outsideBodies, 
                  CUBIT_FALSE, keep_old);
  
  return stat;
}

//===============================================================================
// Function   : unite
// Member Type: PUBLIC
// Description: unite boolean operation between OCC-based bodies
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus     OCCModifyEngine::unite(DLIList<BodySM*> &bodies, 
                                       DLIList<BodySM*> &newBodies,
                                       bool keep_old) const
{
  if(bodies.size() < 2)
  {
    newBodies = bodies;
    return CUBIT_SUCCESS;
  }

  //In order to distinguish bodies who are not intersecting each other to 
  //avoid doing the boolean, check the minimum distance of the bodies first.
  DLIList<BodySM*> revised_bodies;
  DLIList<BodySM*> overlaped_bodies;
  DLIList<TopoDS_Shape*> overlap_shapes;
  DLIList<TopoDS_Shape*> shape_list;
  DLIList<TopoDS_Shape*> revised_shapes;
  DLIList<CubitBoolean> is_volume;
  CubitStatus stat =
        get_shape_list(bodies, shape_list, is_volume, keep_old); 

  if( !stat )
    return CUBIT_FAILURE;

  while(bodies.size() > 0)
  {
    TopoDS_Shape *shape1 = shape_list.pop();
    BodySM* first_body = bodies.pop();
    CubitBoolean intersect = false; 
    shape_list.reset();
    int size = shape_list.size();

    overlaped_bodies.clean_out();
    overlap_shapes.clean_out();
    for (int k = 0 ; k < size; k++)
    {
      TopoDS_Shape *shape2 = shape_list.get();
      BodySM* sec_body = bodies.get();
      BRepExtrema_DistShapeShape dist(*shape1, *shape2);
      //dist.Perform();
      if(dist.IsDone() && dist.Value() < TOL)
      {
        overlaped_bodies.append(sec_body);
        overlap_shapes.append(shape2);
        bodies.change_to(NULL);
        shape_list.change_to(NULL);
        intersect = true;
      }
      bodies.step();
      shape_list.step();
    }
    bodies.remove_all_with_value(NULL);
    shape_list.remove_all_with_value(NULL);
    if (intersect != true)
    {
      revised_bodies.append(first_body);
      revised_shapes.append(shape1);
    }
    else
    {
      overlaped_bodies.append(first_body);
      overlap_shapes.append(shape1); 

      //find a non-sheet body to be the first shape
      TopoDS_Shape* first_shape;
      TopoDS_Shape* second_shape;
      CubitBoolean first_is_volume = false;
      for(int i = 0; i < overlap_shapes.size(); i++)
      {
        
        first_shape = overlap_shapes.get();
        if(first_shape->ShapeType() < TopAbs_SHELL)
        {
          first_is_volume = true;
          break;
        }
        overlap_shapes.step();
        overlaped_bodies.step();
      }
      overlap_shapes.remove(first_shape);
      BodySM* body1 = overlaped_bodies.remove();
      BodySM* body2 = NULL;
      overlap_shapes.reset();
      overlaped_bodies.reset();
      for(int i = 0; i < overlap_shapes.size(); i++)
      {
        second_shape = overlap_shapes.get_and_step();
        body2 = overlaped_bodies.get_and_step(); 
        BRepAlgoAPI_Fuse fuser(*first_shape, *second_shape);
        TopoDS_Shape new_shape = fuser.Shape();

        CubitBoolean has_changed;
        TopTools_IndexedMapOfShape M1, M2, M_new;
        TopExp::MapShapes(*first_shape, TopAbs_SOLID, M1); 
        TopExp::MapShapes(*second_shape, TopAbs_SOLID, M2);
        TopExp::MapShapes(new_shape, TopAbs_SOLID, M_new);
        if(M_new.Extent() == 1 && M1.Extent() > 1 && M2.Extent() == 1)
        {
          check_operation(new_shape,second_shape, is_volume[i], has_changed, &fuser, keep_old);
          check_operation(new_shape, first_shape, first_is_volume, has_changed, &fuser, keep_old);
        }
 
        else if(M_new.Extent() == 1 && M1.Extent() > 1 && M2.Extent() > 1)
        //two compound bodies unite into one solid lump body
        {
          OCCQueryEngine::instance()->copy_attributes(*first_shape, new_shape);
          OCCQueryEngine::instance()->copy_attributes(*second_shape, new_shape);
          first_shape = new TopoDS_Shape(new_shape);
          OCCQueryEngine::instance()->delete_solid_model_entities(body1);
          OCCQueryEngine::instance()->delete_solid_model_entities(body2);
        } 
        else
        {
          check_operation(new_shape, first_shape, first_is_volume, has_changed, &fuser, keep_old);

          check_operation(new_shape,second_shape, is_volume[i], has_changed, &fuser, keep_old);
        }
      }
      //ok, we're done with all unites, construct new Body'
      DLIList<TopologyBridge*> tbs;
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*first_shape);

      BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
      if (bodysm)
      {
        CAST_TO(bodysm, OCCBody)->get_TopoDS_Shape(first_shape);
        revised_bodies.append(bodysm);
        revised_shapes.append(first_shape);
      }
    } 
  }

  if (revised_bodies.size() > 1)
  {
    revised_shapes.clean_out();
    is_volume.clean_out();
    stat =
        get_shape_list(revised_bodies, revised_shapes, is_volume, true);

    //simply make all bodies into a compound
    TopoDS_Compound Co ;
    BRep_Builder B;
    B.MakeCompound(Co);
    for(int i = 0; i < revised_shapes.size(); i++)
      B.Add(Co, *revised_shapes.get_and_step());

    BodySM* body = OCCQueryEngine::instance()->populate_topology_bridge(Co);
    if(body)
      newBodies.append(body);

    //delete all other individial bodies
    if(keep_old == false)
      for(int i = 0; i < revised_bodies.size(); i++)
        OCCQueryEngine::instance()->
          delete_solid_model_entities(revised_bodies.get_and_step());
    revised_bodies.clean_out();
  }

  else if(revised_bodies.size() == 1)
    newBodies = revised_bodies;
  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::thicken( DLIList<BodySM*>& bodies,
                                      DLIList<BodySM*>& new_bodies,
                                      double depth,
                                      CubitBoolean both) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : hollow
// Member Type: PUBLIC
// Description: Hollow existing solid body by remove one or several surfaces 
//              Can only take one body at a time.
//              depth > 0, thick body going outside bodies
//              depth < 0, thick body going inside bodies
// Author     : Jane Hu 
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine::hollow( DLIList<BodySM*>& bodies, 
                                     DLIList<Surface*>& surfs_to_remove,
                                     DLIList<BodySM*>& new_bodies,
                                     double depth) const
{
  if(bodies.size() != 1 || surfs_to_remove.size() < 1)
  {
    PRINT_ERROR("Making thick solid in OCC will take one body and at least one surface at a time.\n"); 
    return CUBIT_FAILURE;
  }

  DLIList<TopoDS_Shape*> shape_list;
  DLIList<CubitBoolean> is_volume;
  CubitStatus stat = get_shape_list(bodies, shape_list, is_volume, CUBIT_FALSE);

  if(!stat)
    return stat;

  if(!is_volume.get())//sheet body
  {
    PRINT_ERROR("Making thick solid in OCC needs an initial solid body to hollow with.\n");
    return CUBIT_FAILURE;
  }

  //make sure the body to be hollowed has only one lump
  OCCBody* occ_body = CAST_TO(bodies.get(), OCCBody);
  DLIList<Lump*> lumps;
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCShell*> shells;
  surfaces = occ_body->my_sheet_surfaces();
  shells = occ_body->shells();
  lumps = occ_body->lumps();
  if(lumps.size()!=1 || surfaces.size() != 0 || shells.size() != 0)
  {
    PRINT_ERROR("bodies with more than one lump can't be hollowed to make a thick body.\n");
    return CUBIT_FAILURE;
  }

  //make sure surfs_to_remove are all in bodies
  TopTools_ListOfShape face_shapes;
  occ_body->get_all_surfaces(surfaces);
  for(int i = 0; i < surfs_to_remove.size(); i++)
  {
    OCCSurface* occ_surf = CAST_TO(surfs_to_remove.get(), OCCSurface);
    if(!occ_surf)
      continue;
    if(!surfaces.is_in_list(occ_surf))
      continue;
    TopoDS_Face * face = occ_surf->get_TopoDS_Face();
    face_shapes.Append(*face); 
  }

  if(face_shapes.IsEmpty())
  {
    PRINT_ERROR("The surfaces provided should be from the body to be hollowed.\n");
    return CUBIT_FAILURE;
  }
  
  double dTOL = 1.e-3; //hard coded for now, can be changed by application
  TopoDS_Shape* solid = shape_list.get();
  BRepOffsetAPI_MakeThickSolid hollower(*solid, face_shapes, depth, dTOL,
                                        BRepOffset_Skin, Standard_False,
                                        Standard_False, GeomAbs_Intersection);
  TopoDS_Shape new_shape = hollower.Shape();
  if(solid->ShapeType() == TopAbs_SOLID)
  {
    TopoDS_Solid old_solid = TopoDS::Solid(*solid);
    OCCLump::update_OCC_entity(old_solid , new_shape, &hollower); 
  }
  else if(solid->ShapeType() == TopAbs_COMPOUND)
    OCCBody::update_OCC_entity(*solid, new_shape, &hollower);
 
  //ok, we're done with all hollowing, construct new Body'
  DLIList<TopologyBridge*> tbs;
  tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }

  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : flip_normals
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 06/08
//===============================================================================
CubitStatus OCCModifyEngine :: flip_normals( DLIList<Surface*>& face_list ) const
{
  return CUBIT_FAILURE;
/*
  DLIList<Surface*> surface_list;
  while (face_list.size())
  {
    OCCSurface* occ_surface = CAST_TO(face_list.pop(), OCCSurface);
    OCCShell* occ_shell = occ_surface->my_shell();
    DLIList<OCCSurface*> surfaces;
    surfaces.append(occ_surface);
    if(occ_shell) //find all surfaces in face_list that belong to this shell
    {
      int size = face_list.size();
      for ( int i = 0; i < size; i++)
      {
        occ_surface = CAST_TO(face_list.get(), OCCSurface); 
        if(occ_shell == occ_surface->my_shell())
          surfaces.append(CAST_TO(face_list.remove(),OCCSurface));
        else
          face_list.step();
      } 
      
      if (!occ_shell->is_sheet())
      {
        DLIList<OCCSurface*> memberSurfs = occ_shell->getMemberSurfaces();
        for (int i = 0; i < memberSurfs.size(); ++i)
        {
          occ_surface = memberSurfs.get_and_step();
          if(surfaces.is_in_list(occ_surface))
          {
            // would need to implement flipping the normal here
            surface_list.append(occ_surface);
          }
        }
      }
    }        
    if(!occ_shell || occ_shell->is_sheet()) //sheet body 
    {
      // would need to implement flipping the normal here
      surface_list.append(occ_surface);
      PRINT_INFO( "Modified volume\n" );
    }
  }
  face_list = surface_list;
  return CUBIT_SUCCESS;
*/
}


//===============================================================================
// Function   : sweep_translational
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 09/08
//===============================================================================
CubitStatus OCCModifyEngine::sweep_translational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& sweep_vector,
  double draft_angle, //in Radius
  int draft_type, //RightCorner=1 or RoundCorner =2
  bool switchside,//not used, shell and surfaces are one sided, not like Acis
  bool rigid, //not used here, in Acis, it means whether the end surface is
              // parallel to the starting surface, or perpendicular to the path
  bool anchor_entity, //not used in OCC
  bool keep_old) const
{
  //in OCC, there's no sweep surface with draft option, this can be done by
  //creating draft shell then make solid to achieve.
  //while if draft_angle is 0, directly use sweep functions.

  gp_Dir adir(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());
  gp_Vec aVec(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());
 
  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, &sweep_vector, toposhape);   
      if(!stat)
        continue;
    }
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
    }
  
    DLIList<TopologyBridge*> tbs;
    //create the draft or the sweep
    BodySM* bodysm = NULL;
    if( draft_angle == 0.)
    {
      BRepSweep_Prism swept(toposhape, aVec);
      TopoDS_Shape new_shape = swept.Shape();

      tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
      assert(tbs.size() == 1);

      bodysm = CAST_TO(tbs.get(), BodySM); 
    }

    else
    {
      BRepOffsetAPI_MakeDraft draft(toposhape, adir, draft_angle);
      BRepBuilderAPI_TransitionMode Cornertype;
      if(draft_type == 1)
        Cornertype = BRepBuilderAPI_RightCorner;
      else 
        Cornertype = BRepBuilderAPI_RoundCorner;

      draft.SetOptions(Cornertype);
      draft.Perform(sweep_vector.length());
      TopoDS_Shape new_shape = draft.Shape();

      tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);

      assert(tbs.size() == 1);

      bodysm = CAST_TO(tbs.get(), BodySM);
      if(bodysm && surface != NULL) //only gets swept side and original surfaces
      {
         //get surfaces from the shell body and create a top surface to
         //make a swept solid.
         DLIList<OCCShell*> shells = CAST_TO(bodysm, OCCBody)->shells();
         if(shells.size() == 0)
         {
           PRINT_WARNING("Sweep surface failed inside OCC engine.\n");
           return CUBIT_FAILURE;
         }
         assert(shells.size() == 1);
         DLIList<OCCSurface*> memberSurfaces =
             shells.get()->getMemberSurfaces();
         DLIList<Surface*> surface_list;
         for (int i = 0; i < memberSurfaces.size(); ++i)
           surface_list.append(memberSurfaces.get_and_step());

         //create the top surface from edges.
         DLIList<OCCCoEdge*> coedges;
         for(int i = 0; i < surface_list.size(); i++)
           CAST_TO(surface_list.get_and_step(), OCCSurface)->get_coedges(coedges);
         for(int i = 0; i < coedges.size(); i++)
         {
           OCCCoEdge* coedge = coedges[i];
           if(coedge == NULL)
             continue;
           for(int j = i+1; j < coedges.size(); j++)
           {
              OCCCoEdge* temp_coedge = coedges[j];
              if(temp_coedge == NULL)
                continue; 
              if(coedge->curve() == temp_coedge->curve()) //Since the shell 
              // is open, the sense of curve can be either the same or opposite. 
              {
                coedges.move_to(coedge);
                coedges.change_to((OCCCoEdge*)NULL);
                coedges.move_to(temp_coedge);
                coedges.change_to((OCCCoEdge*)NULL);
              }
           }
         } 
         coedges.remove_all_with_value(NULL);
         assert(coedges.size() > 0);
         DLIList<Curve*> curves;
         for(int i = 0; i < coedges.size(); i++)
           curves.append(coedges.get_and_step()->curve());

         Surface* surf = make_Surface(PLANE_SURFACE_TYPE, curves);
         if(!surf)
           surf = make_Surface(BEST_FIT_SURFACE_TYPE, curves);
         if(!surf)
         {
           PRINT_ERROR("Can't calculate for the top surface.\n");
           continue;
         }
         surface_list.append(surf);
         DLIList<BodySM*> bodies;
         create_solid_bodies_from_surfs(surface_list, bodies);

         if(bodies.size() == 1)
           bodysm = bodies.get();
         else
         {
           PRINT_WARNING("Sweep surface failed in creating solid.\n");
           return CUBIT_FAILURE;
         }
       }
    }
    if(bodysm && !keep_old && surface != NULL)
    {
      //have to unite the new geometry with the old one.
      DLIList<BodySM*> bodies;
      DLIList<OCCBody*> occ_bodies;
      surface->get_bodies(occ_bodies);
      if(occ_bodies.size() == 1)
      {
        OCCBody* old_body = occ_bodies.get();
        //delete sheet body if surface is standalong.
        if (old_body->is_sheet_body())
        {
          OCCQueryEngine::instance()->delete_solid_model_entities(old_body);
          result_body_list.append(bodysm);
        }
        else
        {
          bodies.append(CAST_TO(old_body, BodySM));
          bodies.append(bodysm);
          DLIList<BodySM*> newBodies;
          bool keep_old = CUBIT_FALSE;
          CubitStatus stat = unite(bodies, newBodies, keep_old);
          if(stat)
            result_body_list.append(newBodies.get());
        }
      }
    }
    else if (bodysm)
      result_body_list.append(bodysm);
  }
  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::get_sweepable_toposhape(OCCCurve*& curve,
                                              TopoDS_Shape& toposhape)const
{
  TopoDS_Edge *edge = curve->get_TopoDS_Edge( );
  BRepBuilderAPI_Copy api_copy(*edge);
  toposhape = api_copy.ModifiedShape(*edge);
  TopoDS_Edge new_edge = TopoDS::Edge(toposhape);
  toposhape = BRepBuilderAPI_MakeWire(new_edge);

  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::get_sweepable_toposhape(OCCSurface*& surface,
                                                  const CubitVector* sweep_v_p,
                                                  TopoDS_Shape& toposhape)const
{
  GeometryEntity* ref_ent = surface;

  if(surface != NULL)
  {
    //Make copy of the surface's topo_shape.
    TopoDS_Shape* toposhape_prt =
          OCCQueryEngine::instance()->get_TopoDS_Shape_of_entity(ref_ent);

    if(!toposhape_prt)
    {
      PRINT_WARNING("GeometryEntity without TopoDS_Shape found.\n");
      return CUBIT_FAILURE;
    }

    BRepBuilderAPI_Copy api_copy(*toposhape_prt);
    toposhape = api_copy.ModifiedShape(*toposhape_prt);
    if(sweep_v_p)
    {
      CubitVector center = surface->center_point();
      CubitVector normal;
      surface->closest_point(center,NULL,&normal);
      CubitVector sweep_vector = *sweep_v_p;
      // TODO: Determine where there might be need to have the normal vector
      // in the opposite direction if normal % sweep_vector > 0
      if (normal % sweep_vector == 0)
      {
        PRINT_ERROR("Sweeping direction should not be on the surface.\n");
        return CUBIT_FAILURE;
      }
    }
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine::sweep_perpendicular(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  double distance,
  double draft_angle,
  int draft_type,
  bool switchside, //has no effect
  bool rigid, //has no effect
  bool anchor_entity, //not used in OCC
  bool keep_old) const
{
  //find the vector perpendicular to the ref_ent normal, and sweep_translate
  //the 'distance' along this vector
  DLIList<GeometryEntity*> edge_list;
  CubitVector vec;
  for(int i = 0; i < ref_ent_list.size(); i++)
  {
     GeometryEntity *ref_ent = ref_ent_list.get_and_step();
     Surface *face = CAST_TO(ref_ent, Surface);
     Curve* edge = CAST_TO(ref_ent, Curve);
     DLIList<GeometryEntity*> face_list;
     if(face != NULL)
     {
        OCCSurface* occ_face = CAST_TO(face, OCCSurface);
        CubitVector center = occ_face->center_point();
        CubitVector closest_p, unit_normal;
        CubitStatus stat = 
                    occ_face->closest_point(center, &closest_p, &unit_normal);
        if(stat)
        {
          vec = distance * unit_normal;
          face_list.append(ref_ent);
          stat = sweep_translational(face_list, result_body_list, vec, 
                                     draft_angle, draft_type, switchside,
                                     rigid, anchor_entity, keep_old);
       }
     }
     else if (edge != NULL)
     {
        edge_list.append(ref_ent);
     }
  }
  if(edge_list.size())
    PRINT_ERROR("Curves cannot be swept perpendicularly, please use the vector sweep.\n");

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : sweep_rotational
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu 
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_rotational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& point,
  const CubitVector& direction,
  double angle, //in radians
  int steps,  //not used
  double draft_angle, //not used
  int draft_type,  //not used
  bool switchside, //not used
  bool make_solid,
  bool rigid,  //not used
  bool anchor_entity, //not used
  bool keep_old ) const  
{
  gp_Dir adir(direction.x(), direction.y(), direction.z()); 
  gp_Pnt pt = gp_Pnt( point.x(), point.y(), point.z());
  gp_Ax1 axis = gp_Ax1(pt, adir);

  gp_Lin line = gp_Lin(axis);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(line);
  OCCCurve* acurve = CAST_TO(OCCQueryEngine::instance()->populate_topology_bridge(edge, CUBIT_TRUE), OCCCurve);
  assert(acurve);

  CubitVector start;
  CubitVector end;

  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface or curve for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, (CubitVector*)NULL, toposhape);
      if(!stat)
        continue;
      //only non-intersecting of surface and axis can be swept.
      DLIList<CubitVector> intersect_pts;
      OCCQueryEngine::instance()->get_intersections(acurve, surface,
                                    intersect_pts, CUBIT_TRUE);
      if(intersect_pts.size() > 0)
      { 
        PRINT_ERROR("Only surfaces with no intersection point with the axis can be revolve-swept.\n");
        continue;
      } 
    }
    else if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
      //closed curve can't intersect with the axis, while open curve can only
      //intersect the axis at the end points. 
      //only curve not intersecting with axis in curve's middle locations
      //can be revolved
      DLIList<CubitVector> intersect_pts;
      OCCQueryEngine::instance()->get_intersections(curve, acurve,
                                  intersect_pts, CUBIT_TRUE, CUBIT_TRUE);
      if(!toposhape.Closed())
      {
        //get start and end points
        DLIList<OCCPoint*> point_list;
        curve->get_points(point_list);
        assert(2 == point_list.size());
        GeometryType type = curve->geometry_type();
        start = point_list.get_and_step()->coordinates();
        end = point_list.get()->coordinates();
        CubitBoolean start_int = CUBIT_FALSE;
        CubitBoolean end_int = CUBIT_FALSE;
        if(intersect_pts.size() > 0)
        {
          CubitBoolean non_int = CUBIT_FALSE;
          for(int i = 0; i < intersect_pts.size(); i++)
          {
             CubitVector &prt = intersect_pts.get_and_step();
             if(prt.distance_between(start) > TOL &&
                prt.distance_between(end) > TOL)
             {
                non_int = CUBIT_TRUE;
                PRINT_ERROR("Only curves with no intersection point with the axis can be revolve-swept.\n");
                break;
             }
             else if(prt.distance_between(start) <= TOL)
                start_int = CUBIT_TRUE;
             else if(prt.distance_between(end) <= TOL)
                end_int = CUBIT_TRUE;
          }
          if(non_int)
            continue;
          if(start_int && end_int && type == STRAIGHT_CURVE_TYPE)
          {
            PRINT_ERROR("Sweep along curve itself is not allowed.\n");
            continue;
          } 
        }
      }
      else
      {
        if(intersect_pts.size() > 0)
        {
          PRINT_ERROR("Only curves with no intersection point with the axis can be revolve-swept.\n");
          continue;
        }  
      }
    } 
    else
    {
      PRINT_ERROR("Only surface or curve can be revolve-swept.\n");
      continue;
    }
    TopoDS_Shape new_shape;
    DLIList<TopologyBridge*> tbs;
    if(make_solid && curve != NULL )
    //giving an open wire and want a solid
    {
      if(!toposhape.Closed())
      {
        //project the start and end points onto the axis
        CubitBoolean start_closed = CUBIT_FALSE;
        CubitBoolean end_closed = CUBIT_FALSE;
        if(acurve->point_containment(start) != CUBIT_PNT_OFF)
          start_closed = CUBIT_TRUE;
        if(acurve->point_containment(end) != CUBIT_PNT_OFF)
          end_closed = CUBIT_TRUE; 
        CubitVector start_proj, end_proj;
        TopoDS_Edge edge1, edge2;
        BRepBuilderAPI_MakeWire m_wire;
        if(!start_closed)
        {
          acurve->closest_point(start, start_proj);
          gp_Pnt pt1 = gp_Pnt( start.x(), start.y(), start.z()); 
          gp_Pnt pt2 = gp_Pnt( start_proj.x(), start_proj.y(), start_proj.z());
          edge1 = BRepBuilderAPI_MakeEdge(pt1, pt2);
          m_wire.Add(edge1);
          m_wire.Add(TopoDS::Wire(toposhape));
        }
        else
        {
          m_wire.Add(TopoDS::Wire(toposhape));
          start_proj = start;
        }
 
        if(!end_closed)
        {
          acurve->closest_point(end,end_proj);
          gp_Pnt pt1 = gp_Pnt( end.x(), end.y(), end.z());
          gp_Pnt pt2 = gp_Pnt( end_proj.x(), end_proj.y(), end_proj.z());
          edge2 = BRepBuilderAPI_MakeEdge(pt1, pt2);
          m_wire.Add(edge2);
        }
      
        else
          end_proj = end;
        
        gp_Pnt pt1 = gp_Pnt( end_proj.x(), end_proj.y(), end_proj.z());
        gp_Pnt pt2 = gp_Pnt( start_proj.x(), start_proj.y(), start_proj.z());
        TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(pt2, pt1);
        m_wire.Add(edge3);
      
        TopoDS_Wire wire = m_wire.Wire();
        toposhape = BRepBuilderAPI_MakeFace(wire);
      }
      else //closed
      {
        TopoDS_Wire wire = TopoDS::Wire(toposhape);
        toposhape = BRepBuilderAPI_MakeFace(wire);
      }
    }
    BRepSweep_Revol revol(toposhape, axis, angle);
    new_shape = revol.Shape();

    tbs += OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
    assert(tbs.size() == 1);

    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);

    if(bodysm && !keep_old && surface != NULL)
    {
      //have to unite the new geometry with the old one.
      DLIList<BodySM*> bodies; 
      DLIList<OCCBody*> occ_bodies;
      surface->get_bodies(occ_bodies);
      if(occ_bodies.size() == 1)
      {
        OCCBody* old_body = occ_bodies.get();
        //delete sheet body if surface is standalong.
        if (old_body->is_sheet_body())
        {
          OCCQueryEngine::instance()->delete_solid_model_entities(old_body);
          result_body_list.append(bodysm);
        }
        else
        {
          bodies.append(CAST_TO(old_body, BodySM));
          bodies.append(bodysm);
          DLIList<BodySM*> newBodies;
          CubitStatus stat = unite(bodies, newBodies, CUBIT_FALSE);
          if(stat)
            result_body_list.append(newBodies.get());
        }
      }
    }
    else if (bodysm)
      result_body_list.append(bodysm);
    continue;
  }
  OCCQueryEngine::instance()->delete_solid_model_entities(acurve);
  if(result_body_list.size()>0)
    return CUBIT_SUCCESS;
  else 
    return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_along_curve
// Member Type: PUBLIC
// Description: The ref_edge_list must provide a list of curves which are
//              connected, and making G1 continuous wire.
// Author     : Jane Hu
// Date       : 10/08
//===============================================================================
CubitStatus OCCModifyEngine::sweep_along_curve(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  DLIList<Curve*>& ref_edge_list,
  double draft_angle, //only used for straight curve case
  int draft_type, //only used for straight curve case
  bool rigid, //not used
  bool anchor_entity, //not used
  bool keep_old) const 
{
  //make wire out of ref_edge_list
  BRepBuilderAPI_MakeWire awire;
  TopTools_ListOfShape L;
  OCCCurve* occ_curve = NULL;
  GeometryType type = UNDEFINED_CURVE_TYPE;
  int num_curve = 0;
  for(int i = 0; i < ref_edge_list.size(); i++)
  {
    Curve* curve = ref_edge_list.get_and_step();
    occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    TopoDS_Edge* topoedge = occ_curve->get_TopoDS_Edge( );
    BRepBuilderAPI_Copy api_copy(*topoedge);
    TopoDS_Shape newShape = api_copy.ModifiedShape(*topoedge);
    L.Append(newShape);
    type = occ_curve->geometry_type();
    num_curve++;
  }
  if(L.IsEmpty())
  {
    PRINT_ERROR("There's no valid sweeping path.\n");
    return CUBIT_FAILURE;
  }
  
  if(num_curve == 1 && type == STRAIGHT_CURVE_TYPE && draft_angle != 0.0)
  {
    DLIList<OCCPoint*> point_list;
    occ_curve->get_points(point_list);
    CubitVector v1 = point_list.get_and_step()->coordinates();
    CubitVector v2 = point_list.get()->coordinates();
    CubitVector sweep_vector = v2-v1;
    return sweep_translational(ref_ent_list,result_body_list,sweep_vector,
                               draft_angle, draft_type, CUBIT_FALSE, 
                               rigid, anchor_entity, keep_old); 
  }
  awire.Add(L);
  TopoDS_Wire wire;
  wire = awire.Wire();

  BRepTools_WireExplorer it(wire);
  int num_edges = 0;
  for(; it.More(); it.Next())
    num_edges++; 
  
  BRepLib_FuseEdges fuser(wire);
  fuser.SetConcatBSpl();
  fuser.Perform();
  TopoDS_Shape  spline = fuser.Shape();
  wire = TopoDS::Wire(spline);

  DLIList<TopologyBridge*> tbs;
  for (int i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    //Make copy of the surface or curve for later to build solid.
    OCCSurface* surface = CAST_TO(ref_ent, OCCSurface);
    OCCCurve* curve = CAST_TO(ref_ent, OCCCurve);
    TopoDS_Shape toposhape ;
    if(surface != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(surface, (CubitVector*)NULL, toposhape);
      if(!stat)
        continue;
    } 
    else if(curve != NULL)
    {
      CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
      if(!stat)
        continue;
    }

    //sweep along the wire
    BRepOffsetAPI_MakePipe maker(wire, toposhape);
    if(!maker.IsDone())
    {
      PRINT_ERROR("Can't sweep along the provided curve(s).\n");
      continue;
    }
    TopoDS_Shape newShape = maker.Shape();
    
    tbs += OCCQueryEngine::instance()->populate_topology_bridge(newShape);
    assert(tbs.size() == 1);

    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);

    if(bodysm && !keep_old && surface != NULL)
    {
      //have to unite the new geometry with the old one.
      DLIList<BodySM*> bodies;
      DLIList<OCCBody*> occ_bodies;
      surface->get_bodies(occ_bodies);
      if(occ_bodies.size() == 1)
      {
        OCCBody* old_body = occ_bodies.get();
        //delete sheet body if surface is standalong.
        if (old_body->is_sheet_body())
        {
          OCCQueryEngine::instance()->delete_solid_model_entities(old_body);
          result_body_list.append(bodysm);
        }
        else
        {
          bodies.append(CAST_TO(old_body, BodySM));
          bodies.append(bodysm);
          DLIList<BodySM*> newBodies;
          CubitStatus stat = unite(bodies, newBodies, CUBIT_FALSE);
          if(stat)
            result_body_list.append(newBodies.get());
        }
      }
    }
    else if (bodysm)
      result_body_list.append(bodysm);
  }
  return CUBIT_SUCCESS;
}

CubitStatus OCCModifyEngine::sweep_to_body(
                                   DLIList<Curve*> curve_list,
                                   BodySM *target_body,
                                   CubitVector sweep_vector,
                                   DLIList<BodySM*> &new_bodies,
                                   bool unite) const
{
  TopoDS_Shape *stop_shape = NULL;
  OCCBody* occ_body = CAST_TO(target_body, OCCBody);
  occ_body->get_TopoDS_Shape(stop_shape);

  gp_Dir adir(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());

  DLIList<BodySM*> swept_bodies;
  for(int i = 0; i < curve_list.size(); i++)
  {
    OCCCurve* curve = CAST_TO(curve_list.get_and_step(), OCCCurve);
    TopoDS_Shape toposhape ;
    CubitStatus stat = get_sweepable_toposhape(curve, toposhape);
    if(!stat)
      continue;
    BRepOffsetAPI_MakeDraft draft(toposhape, adir, 0.0);
    draft.Perform(*stop_shape);
    TopoDS_Shape new_shape = draft.Shape();
    DLIList<TopologyBridge*> tbs;
    tbs = OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
    assert(tbs.size() == 1); 
    BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
    if (bodysm)
      swept_bodies.append(bodysm);
  }
  if(unite)
    return this->unite(swept_bodies, new_bodies, CUBIT_FALSE);
  new_bodies += swept_bodies;
  return CUBIT_SUCCESS; 
}

CubitStatus OCCModifyEngine::sweep_to_body(
                                   Surface  *source_surface,
                                   BodySM *target_body,
                                   CubitVector sweep_vector,
                                   DLIList<BodySM*> &new_bodies ) const
{
  TopoDS_Shape *stop_shape = NULL;
  OCCBody* occ_body = CAST_TO(target_body, OCCBody);
  occ_body->get_TopoDS_Shape(stop_shape);

  gp_Dir adir(sweep_vector.x(), sweep_vector.y(), sweep_vector.z());

  DLIList<BodySM*> swept_bodies;
  OCCSurface* surf = CAST_TO(source_surface, OCCSurface);
  TopoDS_Shape toposhape ;
  CubitStatus stat = get_sweepable_toposhape(surf, &sweep_vector, toposhape); 
  if(!stat)
    return CUBIT_FAILURE;
  BRepOffsetAPI_MakeDraft draft(toposhape, adir, 0.0);
  draft.Perform(*stop_shape);
  TopoDS_Shape new_shape = draft.Shape();
  DLIList<TopologyBridge*> tbs;
  tbs = OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
  assert(tbs.size() == 1);
  BodySM* bodysm = CAST_TO(tbs.get(), BodySM);
  //get surfaces from the shell body and create a top surface to
  //make a swept solid.
  DLIList<OCCShell*> shells = CAST_TO(bodysm, OCCBody)->shells();
  if(shells.size() == 0)
  {
    PRINT_WARNING("Sweep surface failed inside OCC engine.\n");
    return CUBIT_FAILURE;
  }
  assert(shells.size() == 1);
  DLIList<OCCSurface*> memberSurfaces = shells.get()->getMemberSurfaces();
  DLIList<Surface*> surface_list;
  for (int i = 0; i < memberSurfaces.size(); ++i)
    surface_list.append(memberSurfaces.get_and_step());

  //create the top surface from edges.
  DLIList<OCCCoEdge*> coedges;
  for(int i = 0; i < surface_list.size(); i++)
    CAST_TO(surface_list.get_and_step(), OCCSurface)->get_coedges(coedges);
  for(int i = 0; i < coedges.size(); i++)
  {
    OCCCoEdge* coedge = coedges[i];
    if(coedge == NULL)
      continue;
    for(int j = i+1; j < coedges.size(); j++)
    {
      OCCCoEdge* temp_coedge = coedges[j];
      if(temp_coedge == NULL)
        continue;
      if(coedge->curve() == temp_coedge->curve()) //Since the shell
      // is open, the sense of curve can be either the same or opposite.
      {
         coedges.move_to(coedge);
         coedges.change_to((OCCCoEdge*)NULL);
         coedges.move_to(temp_coedge);
         coedges.change_to((OCCCoEdge*)NULL);
       }
    }
  }
  coedges.remove_all_with_value(NULL);
  assert(coedges.size() > 0);
  DLIList<Curve*> curves;
  for(int i = 0; i < coedges.size(); i++)
    curves.append(coedges.get_and_step()->curve());

  Surface* surface = make_Surface(PLANE_SURFACE_TYPE, curves);
  if(!surface)
    surface = make_Surface(BEST_FIT_SURFACE_TYPE, curves);
  if(!surface)
  {
    PRINT_ERROR("Can't calculate for the top surface.\n");
    return CUBIT_FAILURE;
  }
  surface_list.append(surface);
  DLIList<BodySM*> bodies;
    create_solid_bodies_from_surfs(surface_list, bodies);

  if(bodies.size() == 1)
    bodysm = bodies.get();
  else
  {
    PRINT_WARNING("Sweep surface failed in creating solid.\n");
    return CUBIT_FAILURE;
  }
  if (bodysm)
    new_bodies.append(bodysm);
  return CUBIT_SUCCESS;
}

//HEADER- Webcut-related functions

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                              const CubitVector &v1,
                              const CubitVector &v2,
                              const CubitVector &v3,
                              DLIList<BodySM*>& neighbor_imprint_list,
                              DLIList<BodySM*>& results_list,
                              ImprintType imprint_type,
                              bool preview ) const
{
  CubitStatus stat;
  DLIList<BodySM*> new_BodySMs;

  if(preview)
  {
    TopoDS_Face* p_face;
    stat = get_3_point_plane(v1, v2, v3, p_face);
    if(!stat)
      return stat;
    GfxPreview::clear();
    OCCDrawTool::instance()->draw_TopoDS_Shape(p_face, CUBIT_BLUE_INDEX, CUBIT_TRUE);
    delete p_face;
    p_face = NULL;
    return CUBIT_SUCCESS;
  }

  stat = OCCModifyEngine::instance()->section(webcut_body_list, v1, v2, v3, new_BodySMs, true, true,false);
  if(stat == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't webcut the bodies using a plane determined by 3 points.\n");
    return stat;
  }
  
  stat = OCCModifyEngine::instance()->section(webcut_body_list, v1, v2, v3, new_BodySMs, false, false, false);
  if(stat == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't webcut the bodies using a plane determined by 3 points.\n");
    return stat;
  }

  BodySM* new_body1, *new_body2;
  if(imprint_type > NO_IMPRINT)
  {
    for(int i = 0; i < new_BodySMs.size()-1; i ++)
    {
      BodySM* body1 = new_BodySMs[i];
      for(int j = i+1; j < new_BodySMs.size(); j++)
      {
        BodySM* body2 = new_BodySMs[j];
        stat =  this->imprint( body1, body2, new_body1, new_body2, false);
        if(new_body1 && body1 != new_body1)
          new_BodySMs[i] = new_body1;
        if(new_body2 && body2 != new_body2)
          new_BodySMs[j] = new_body2; 
      }
    }
  }

  results_list = new_BodySMs;

  //now imprint with the neighbors
  if( imprint_type == INCLUDE_NEIGHBORS  )
  {
    // Loop over all the neighboring Bodies
    DLIList<TopoDS_Shape*> shape_list1, shape_list2;
    DLIList<CubitBoolean> is_volume;
    stat = get_shape_list(neighbor_imprint_list,shape_list1,is_volume,false);
    if (!stat)
    {
      	PRINT_WARNING("Can't imprint using neighouring bodies.\n");
        return CUBIT_SUCCESS;
    }    
    is_volume.clean_out();

    stat = get_shape_list(new_BodySMs,shape_list2,is_volume,false);
    { 
        PRINT_WARNING("Can't imprint using neighouring bodies.\n");
        return CUBIT_SUCCESS;
    }
    // Loop over all the neighboring shapes
    DLIList<TopoDS_Face*> face_list;
    for (int i=shape_list1.size(); i--;)
    {
      TopoDS_Shape *neighbor_shape = shape_list1.get_and_step() ;
      for(int j = 0; j < shape_list2.size(); j ++)
      {
        TopoDS_Shape* shape2 = shape_list2[j];
        int result =  this->imprint_toposhapes( shape2, neighbor_shape,
                      face_list);
        if(result == 0)
          shape_list2[j] = shape2;   
      }
    }

    DLIList<TopologyBridge*> tbs;
    results_list.clean_out();
    for(int j = 0; j < shape_list2.size(); j ++)
    {
      TopoDS_Shape* shape2 = shape_list2[j];
        
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape2);
      BodySM* newBody1 = CAST_TO(tbs.get(),BodySM);
      results_list.append(newBody1);    
    }
  }

  return CUBIT_SUCCESS;  
}

//===============================================================================
// Function   : webcuts a list of bodies using another Body as the tool.
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 01/09
//===============================================================================
CubitStatus    OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                                BodySM const* tool_body,
                                DLIList<BodySM*>& neighbor_imprint_list,
                                DLIList<BodySM*>& results_list,
                                ImprintType imprint_type,
                                bool preview ) const
{
  // tool_body and webct_body_list will be kept and webcut result is in 
  //results_list.
  //tool_body is a const pointer points to varible BodySM object
  //here trying to create a non-const pointer points to the same BodySM object.

  BodySM *body = const_cast<BodySM*>(tool_body);
  CubitStatus stat;
  DLIList<BodySM*> tool_bodies;
  if(preview)
    tool_bodies.append(body);

  BodySM *copy_tool = copy_body(body);
 
  if(!preview)
    tool_bodies.append(copy_tool);

  // preview the tool
  if (preview)
  {
    DLIList<TopoDS_Shape*> shape_list;
    DLIList<CubitBoolean> is_volume;
    stat = get_shape_list(tool_bodies,shape_list,is_volume,false);
    if(!stat)
    {
        PRINT_WARNING("tool_body is not an OCC body.\n");
        return CUBIT_FAILURE;
    }
    GfxPreview::clear();
    OCCDrawTool::instance()->draw_TopoDS_Shape(shape_list.get(), 
                             CUBIT_BLUE_INDEX, CUBIT_TRUE);

    return CUBIT_SUCCESS;
  }

  //if the tool_body is a volume, use intersect & subtract,
  //if it's a shell or face body, use section then.
  OCCBody* occ_body = CAST_TO(body, OCCBody);
  DLIList<Lump*> lumps;
  lumps = occ_body->lumps();
  DLIList<OCCShell*> shells;
  shells = occ_body->shells();
  DLIList<OCCSurface*> surfaces;
  surfaces = occ_body->my_sheet_surfaces();

  if(lumps.size() == 0 && (shells.size()==1 || surfaces.size() == 1))
  {
    OCCSurface * surface = NULL;
    TopoDS_Shell* topo_shell = NULL;
    TopoDS_Face *topo_face = NULL;
    if(shells.size() == 1)
    {
      OCCShell* shell = shells.get();
      surface = shell->getMemberSurfaces().get();
      topo_shell = shell->get_TopoDS_Shell();
    }
    else
    {
      surface = surfaces.get();
      topo_face = surface->get_TopoDS_Face();
    }
    CubitVector point_1, normal ;
    CubitStatus rsl = surface->get_point_normal(point_1, normal);
    assert(rsl);

    //compare the tool body's location with the webcut_body_list, to see if
    //one side of the normal will totally cover one or more volumes
    //while the other side cut through one volume, with no whole volume
    //covered. Make sure the tool body's normal direct to single cut-
    //through volume to make sure the boolean operation won't destory
    //webcut_body_list and create a new body.
    CubitBoolean reversed = CUBIT_FALSE;
    if(webcut_body_list.size() == 1)
    {
      DLIList<OCCSurface*> cut_surfaces;
      OCCSurface* cut_face = NULL;
      if(surfaces.size() == 1)
      {
        cut_surfaces = occ_body->my_sheet_surfaces();
        cut_face = cut_surfaces.get();
      }
      else
      {
        DLIList<OCCShell*> shells;
        shells = occ_body->shells();
        cut_face = shells.get()->my_surface();
      }
      if(cut_face != NULL)
      {
        CubitVector center = cut_face->center_point();
        CubitVector normal;
        cut_face->get_point_normal(center, normal);
        DLIList<Lump*> lumps;
        BodySM *from_body = webcut_body_list.get();
        from_body->lumps(lumps);
        int right = 0, left = 0;
        for(int kk = 0; kk < lumps.size(); kk++)
        {
          CubitBox box = lumps.get_and_step()->bounding_box();
          if(box <= center) //center outside of box
          {
            CubitVector v = box.center() - center;
            if(v%normal > 0)
              left ++;
            else
              right++;
          }
        }
        if(right > 0 && left == 0)
          reversed = CUBIT_TRUE;
      }
    }
    gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
    gp_Dir normal_dir(normal.x(), normal.y(), normal.z());
    gp_Vec vec(normal_dir);
    if(reversed)
      vec = -vec;
    pt =  pt.Translated(vec);

    TopoDS_Solid solid;
    if(shells.size() == 1)
      solid = BRepPrimAPI_MakeHalfSpace(*topo_shell, pt);
    else
      solid = BRepPrimAPI_MakeHalfSpace(*topo_face,pt);

    DLIList<CubitBoolean> is_tool_volume;
    is_tool_volume.append(CUBIT_TRUE);
    DLIList<CubitBox*> tool_boxes ;
    Bnd_Box box;
    BRepBndLib::Add(solid, box);
    double min[3], max[3];
    box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);
    CubitBox* cBox = new CubitBox(min, max);

    tool_boxes.append(cBox);
    DLIList<TopoDS_Shape*> solids;
    solids.append(&solid);
    rsl = do_subtract(webcut_body_list, solids, is_tool_volume,
                     &tool_boxes, results_list, CUBIT_TRUE) ;
    delete cBox;
    if(!rsl)
    {
      PRINT_ERROR("Failed to webcut the bodies.\n");
      return CUBIT_FAILURE;
    }
    gp_Vec vec2 =  -2 * vec;
    pt = pt.Translated(vec2);
    if(shells.size() == 1)
      solid = BRepPrimAPI_MakeHalfSpace(*topo_shell, pt);
    else
      solid = BRepPrimAPI_MakeHalfSpace(*topo_face,pt);
    Bnd_Box box2;
    BRepBndLib::Add(solid, box2);
    box2.Get(min[0], min[1], min[2], max[0], max[1], max[2]);
    cBox = new CubitBox(min, max);
    tool_boxes.clean_out();
    tool_boxes.append(cBox);
    solids.clean_out();
    solids.append(&solid);
    rsl = do_subtract(webcut_body_list, solids, is_tool_volume,
                     &tool_boxes, results_list, CUBIT_FALSE) ;
    //make sure the original body's pointer is the first.
    results_list.reverse();
    delete cBox;
    return rsl;
  }

  stat = intersect(body, webcut_body_list, results_list,
                               CUBIT_TRUE);
 
  if(!stat)
  { 
    PRINT_ERROR("Failed to webcut the bodies.\n"); 
    return CUBIT_FAILURE;
  }

  bool imprint = CUBIT_TRUE;
  if(imprint_type == NO_IMPRINT)
    imprint = CUBIT_FALSE;

  stat = subtract(tool_bodies, webcut_body_list, results_list, imprint, 
                  CUBIT_FALSE);

  //intersect doesn't have to imprint option, so first do this imprint.
  BodySM* new_body1, *new_body2;
  if(imprint_type > NO_IMPRINT)
  {
    for(int i = 0; i < results_list.size()-1; i ++)
    {
      BodySM* body1 = results_list[i];
      for(int j = i+1; j < results_list.size(); j++)
      {
        BodySM* body2 = results_list[j];
        stat =  this->imprint( body1, body2, new_body1, new_body2, false);
        if(new_body1 && body1 != new_body1)
          results_list[i] = new_body1;
        if(new_body2 && body2 != new_body2)
          results_list[j] = new_body2;
      }
    }
  }

  //now imprint with the neighbors
  if( imprint_type == INCLUDE_NEIGHBORS  )
  {
    // Loop over all the neighboring Bodies
    DLIList<TopoDS_Shape*> shape_list1, shape_list2;
    DLIList<CubitBoolean> is_volume;
    stat = get_shape_list(neighbor_imprint_list,shape_list1,is_volume,false);
    if (!stat)
    {
        PRINT_WARNING("Can't imprint using neighouring bodies.\n");
        return CUBIT_SUCCESS;
    }
    is_volume.clean_out();

    stat = get_shape_list(results_list,shape_list2,is_volume,false);
    {
        PRINT_WARNING("Can't imprint using neighouring bodies.\n");
        return CUBIT_SUCCESS;
    }
    // Loop over all the neighboring shapes
    DLIList<TopoDS_Face*> face_list;
    for (int i=shape_list1.size(); i--;)
    {
      TopoDS_Shape *neighbor_shape = shape_list1.get_and_step() ;
      for(int j = 0; j < shape_list2.size(); j ++)
      {
        TopoDS_Shape* shape2 = shape_list2[j];
        int result =  this->imprint_toposhapes( shape2, neighbor_shape, face_list);
        if(result == 0)
          shape_list2[j] = shape2;
      }
    }

    DLIList<TopologyBridge*> tbs;
    results_list.clean_out();
    for(int j = 0; j < shape_list2.size(); j ++)
    {
      TopoDS_Shape* shape2 = shape_list2[j];

      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*shape2);
      BodySM* newBody1 = CAST_TO(tbs.get(),BodySM);
      results_list.append(newBody1);
    }
  } 
  return stat;
}

//Dummy functions to fulfill the pure virtural functions in parents, should
//not be called.       ---Jane Hu, 02/17/2012 
CubitStatus OCCModifyEngine::webcut_with_sheet(
                                    DLIList<BodySM*>& webcut_body_list,
                                    BodySM *sheet_body,
                                    DLIList<BodySM*>& neighbor_imprint_list,
                                    DLIList<BodySM*> &new_bodies,
                                    ImprintType imprint_type ,
                                    bool preview )
{
  return CUBIT_FAILURE;                     
}

CubitStatus OCCModifyEngine::webcut_with_extended_sheet(
                                    DLIList<BodySM*> &webcut_body_list,
                                    DLIList<Surface*> &surface_list,
                                    DLIList<BodySM*>& neighbor_imprint_list,
                                    DLIList<BodySM*> &new_bodies,
                                    int &num_cut,
                                    ImprintType imprint_type ,
                                    bool preview )
{
  PRINT_ERROR("This feature is not implemented.\n");
  return CUBIT_FAILURE;
} 

CubitStatus OCCModifyEngine::webcut_with_sweep_curves(
                            DLIList<BodySM*> &blank_bodies,
                            DLIList<Curve*> &curves,
                            const CubitVector& sweep_vector,
                            bool through_all,
                            Surface *stop_surf,
                            Curve *curve_to_sweep_along,
                            DLIList<BodySM*>& neighbor_imprint_list,
                            DLIList<BodySM*> &results_list,
                            ImprintType imprint_type, 
                            CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_curves_rotated(
                            DLIList<BodySM*> &blank_bodies,
                            DLIList<Curve*> &curves,
                            const CubitVector &point,
                            const CubitVector &sweep_axis,
                            double angle,
                            Surface *stop_surf,
                            DLIList<BodySM*>& neighbor_imprint_list,
                            DLIList<BodySM*> &results_list,
                            ImprintType imprint_type,
                            CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces_rotated(
                            DLIList<BodySM*> &blank_bodies,
                            DLIList<Surface*> &surfaces,
                            const CubitVector &point,
                            const CubitVector &sweep_axis,
                            double angle,
                            Surface *stop_surf,
                            bool up_to_next,
                            DLIList<BodySM*>& neighbor_imprint_list,
                            DLIList<BodySM*> &results_list,
                            ImprintType imprint_type ,
                            CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces(
                            DLIList<BodySM*> &blank_bodies,
                            DLIList<Surface*> &surfaces,
                            const CubitVector& sweep_vector,
                            bool sweep_perp,
                            bool through_all,
                            bool outward,
                            bool up_to_next,
                            Surface *stop_surf,
                            Curve *curve_to_sweep_along,
                            DLIList<BodySM*>& neighbor_imprint_list,
                            DLIList<BodySM*> &results_list,
                            ImprintType imprint_type ,
                            CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_cylinder(
                                        DLIList<BodySM*> &webcut_body_list,
                                        double radius,
                                        const CubitVector &axis,
                                        const CubitVector &center,
                                        DLIList<BodySM*>& neighbor_imprint_list,
                                        DLIList<BodySM*>& results_list,
                                        ImprintType imprint_type ,
                                        CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_brick(
                                     DLIList<BodySM*>& webcut_body_list,
                                     const CubitVector &center,
                                     const CubitVector axes[3],
                                     const CubitVector &extension,
                                     DLIList<BodySM*>& neighbor_imprint_list,
                                     DLIList<BodySM*> &results_list,
                                     ImprintType imprint_type ,
                                     CubitBoolean preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_planar_sheet(
                                            DLIList<BodySM*>& webcut_body_list,
                                            const CubitVector &center,
                                            const CubitVector axes[2],
                                            double width, double height,
                                            DLIList<BodySM*>& neighbor_imprint_list,
                                            DLIList<BodySM*> &results_list,
                                            ImprintType imprint_type ,
                                            bool preview )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_curve_loop(
                                         DLIList<BodySM*> &webcut_body_list,
                                         DLIList<Curve*> &ref_edge_list,
                                         DLIList<BodySM*>& neighbor_imprint_list,
                                         DLIList<BodySM*>& results_list,
                                         ImprintType imprint_type ,
                                         bool preview )
{
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : separate_surfaces
// Member Type: PUBLIC
// Description: separate surfaces in shells or compounds into independent ones,
//              still keep the original shell.
// Author     : Jane Hu
// Date       : 02/11
//===============================================================================
CubitStatus OCCModifyEngine::separate_surfaces( DLIList<Surface*> &surf_list,
                                         DLIList<BodySM*> &new_bodies )
{
  DLIList<OCCSurface*> surf_need_work;
  //find out the complexity that the surf_list involves.
  for (int i = 0; i < surf_list.size(); i++)
  {
    Surface* surf = surf_list.get_and_step();
    OCCSurface *surface = CAST_TO(surf, OCCSurface);
    OCCBody* body = surface->my_body();
    if(body != NULL) //either a single surface body or a compound body
    {
      DLIList<Lump*> lumps;
      body->lumps(lumps);
      if (lumps.size() > 0) //compound body
      {
        surf_need_work.append(surface);
        continue;
      }
      DLIList< OCCShell*> shells;
      body->shells(shells);
      if(shells.size() > 0) //compound body
        surf_need_work.append(surface);
      continue;
    }

    OCCShell* shell;
    shell = surface->my_shell();
    if(shell != NULL && shell->my_surface() == NULL) //sheet body
    {
      surf_need_work.append(surface);
      continue;
    }
  }
  
  //create sheet body for surf_need_work list.
  for(int i = 0; i < surf_need_work.size(); i++)
  {
    Surface* copy_surf = make_Surface(surf_need_work.get_and_step());
    if (copy_surf == NULL)
    {
       PRINT_ERROR("Cannot create an OCC sheet bodySM from the given bodySM.\n");
       return CUBIT_FAILURE;
    }

    OCCSurface* occ_surf = CAST_TO(copy_surf, OCCSurface);
    if(occ_surf != NULL)
     new_bodies.append( occ_surf->my_body() );
  }
  return CUBIT_SUCCESS;
}
//===============================================================================
// Function   : section
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 11/08
//===============================================================================
CubitStatus OCCModifyEngine::section( DLIList<BodySM*> &section_body_list,
                                      const CubitVector &point_1,
                                      const CubitVector &point_2,
                                      const CubitVector &point_3,
                                      DLIList<BodySM*>& new_body_list,
                                      bool keep_normal_side,
                                      bool keep_old,
                                      bool keep_both_sides)
{
  if (keep_both_sides == CUBIT_TRUE )
  {
     PRINT_ERROR("Please use webcut to for keep both sides option.\n");
     return CUBIT_FAILURE;
  }
 
  //Calculate normal of the section plan
  CubitVector v1, v2, normal;
  v1 = point_2 - point_1;
  v2 = point_3 - point_1; 
  normal = ~(v1 * v2); 
  if(fabs(normal.length() - 1) > TOL)
  {
     PRINT_ERROR("The three points are co-linear, and can't be used as a cutting plane.\n");
     return CUBIT_FAILURE;
  }
  
  if(keep_normal_side)
    normal *= -1;

  gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
  gp_Dir normal_dir(normal.x(), normal.y(), normal.z()); 
  gp_Pln plane(pt, normal_dir);
  gp_Vec vec(normal_dir);
  pt =  pt.Translated(vec);

  TopoDS_Face face = BRepBuilderAPI_MakeFace(plane);
  TopoDS_Solid solid = BRepPrimAPI_MakeHalfSpace(face, pt);
   
  DLIList<CubitBoolean> is_tool_volume;
  is_tool_volume.append(CUBIT_TRUE);
  DLIList<CubitBox*> tool_boxes ;
  Bnd_Box box;
  BRepBndLib::Add(solid, box);
  double min[3], max[3];
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);
  CubitBox* cBox = new CubitBox(min, max);
  
  tool_boxes.append(cBox);
  DLIList<TopoDS_Shape*> solids;
  solids.append(&solid);
  CubitStatus stat = do_subtract(section_body_list, solids, is_tool_volume,
                     &tool_boxes, new_body_list, keep_old) ;
  delete cBox;
  return stat;
}

//===============================================================================
// Function   : split_body
// Member Type: PUBLIC
// Description: Splits multiple lumps in one body into separate bodies
// Author     : Jane Hu 
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::split_body( BodySM *body_ptr,
                                         DLIList<BodySM*> &new_bodies )
{
  OCCBody* occ_body = CAST_TO(body_ptr, OCCBody);
  if(!occ_body)
  {
     PRINT_ERROR("This is not an OCC body to be split.\n");
     return CUBIT_FAILURE;
  }
  DLIList<Lump*> lumps = occ_body->lumps();
  if(lumps.size() == 1) 
  {
    new_bodies.append(body_ptr);
    return CUBIT_SUCCESS;
  }
  for(int i = 0; i < lumps.size(); i++)
  {
    Lump* lump = lumps.get_and_step();
    OCCLump* occ_lump = CAST_TO(lump, OCCLump);
    OCCSurface* occ_surface = occ_lump->my_sheet_surface();
    //first delete the body which bounds all the stuff.
    if (i == 0)
      OCCQueryEngine::instance()->unhook_BodySM_from_OCC(body_ptr, CUBIT_FALSE);
    if(occ_surface) 
    {
      TopoDS_Face* face = occ_surface->get_TopoDS_Face();
      Surface* surface = 
        OCCQueryEngine::instance()->populate_topology_bridge(*face,
                                                             CUBIT_TRUE);
      new_bodies.append(CAST_TO(surface, OCCSurface)->my_body()); 
      continue;
    }
    OCCShell* occ_shell = occ_lump->my_shell();
    if(occ_shell) 
    {
      TopoDS_Shell* shell = occ_shell->get_TopoDS_Shell();
      OCCShell* ashell = 
         OCCQueryEngine::instance()->populate_topology_bridge(*shell,
                                                              CUBIT_TRUE);
      new_bodies.append(ashell->my_body());
      continue;
    }
    else
    {
      TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();
      Lump* alump = 
        OCCQueryEngine::instance()->populate_topology_bridge(*solid,
                                                             CUBIT_TRUE);
      new_bodies.append(CAST_TO(alump,OCCLump)->get_body());
      continue;
    }
  }
  
  if(lumps.size() > 1)
    delete body_ptr;
  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : reverse_body
// Member Type: PUBLIC
// Description: Turn body inside-out
// Author     : Jane Hu
// Date       : 03/03/09
//===============================================================================
CubitStatus OCCModifyEngine::reverse_body( BodySM* body_ptr )
{
  OCCBody* occ_body = CAST_TO(body_ptr, OCCBody);
  if (!occ_body)
  {
     PRINT_ERROR("Cannot reverse a non-OCC bodySM .\n"
                 "Possible incompatible geometry engines.\n");
     return CUBIT_FAILURE;
  }

  TopoDS_Shape* orig_S;
  TopoDS_Shape S;
  BRep_Builder B;
  occ_body->get_TopoDS_Shape(orig_S);
  S = orig_S->EmptyCopied();
  TopoDS_Iterator it(*orig_S);
  while (it.More()) {
    B.Add(S,it.Value().Reversed());
    it.Next();
  } 
  occ_body->set_TopoDS_Shape(TopoDS::Compound(S));  
  //Bind the new shape and its underlining sub-shapes.
  TopExp_Explorer Ex_orig, Ex;
  int k = -1;
  Ex.Init(S, TopAbs_COMPOUND);
  Ex_orig.Init(*orig_S, TopAbs_COMPOUND);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    if(OCCQueryEngine::instance()->OCCMap->IsBound(Ex.Current()))
    {
      k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());   
      OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
      OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
      TopExp_Explorer Ex_old_solid, Ex_solid;
      Ex_old_solid.Init(*orig_S,TopAbs_SOLID);
      Ex_solid.Init(S, TopAbs_SOLID);
      DLIList<Lump*> lumps = occ_body->lumps();
      for (; Ex_old_solid.More(), Ex_solid.More(); Ex_old_solid.Next(), Ex_solid.Next())
      {
        k = OCCQueryEngine::instance()->OCCMap->Find(Ex_old_solid.Current());
        OCCQueryEngine::instance()->OCCMap->UnBind(Ex_old_solid.Current());
        OCCQueryEngine::instance()->OCCMap->Bind(Ex_solid.Current(), k);
        OCCLump* occ_lump = CAST_TO(lumps.get_and_step(), OCCLump);
        occ_lump->set_TopoDS_Solid(TopoDS::Solid(Ex_solid.Current()));
      }
    } 
    
    else
    {
      Lump *lump = occ_body->lumps().get();
      OCCLump* occ_lump = CAST_TO(lump, OCCLump);
      TopoDS_Solid solid = *(occ_lump->get_TopoDS_Solid());
      k = OCCQueryEngine::instance()->OCCMap->Find(solid);
      OCCQueryEngine::instance()->OCCMap->UnBind(solid);
      TopExp_Explorer Ex_local;
      Ex_local.Init(S, TopAbs_SOLID);
      OCCQueryEngine::instance()->OCCMap->Bind(Ex_local.Current(), k);
      occ_lump->set_TopoDS_Solid(TopoDS::Solid(Ex_local.Current())); 
    }
  }  
      
  Ex.Init(S, TopAbs_SHELL);
  Ex_orig.Init(*orig_S, TopAbs_SHELL);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCShell *shell = (OCCShell*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    shell->set_TopoDS_Shell(TopoDS::Shell(Ex.Current())); 
  }

  Ex.Init(S, TopAbs_FACE);
  Ex_orig.Init(*orig_S, TopAbs_FACE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCSurface *surface = (OCCSurface *)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    TopoDS_Face face = TopoDS::Face(Ex.Current());
    surface->set_TopoDS_Face(face);
  }

  Ex.Init(S, TopAbs_WIRE);
  Ex_orig.Init(*orig_S, TopAbs_WIRE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCLoop* wire = (OCCLoop*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    wire->set_TopoDS_Wire(TopoDS::Wire(Ex.Current()));
  }

  Ex.Init(S, TopAbs_EDGE);
  Ex_orig.Init(*orig_S, TopAbs_EDGE);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCCurve* edge = (OCCCurve*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    edge->set_TopoDS_Edge(TopoDS::Edge(Ex.Current()));
  }

  Ex.Init(S, TopAbs_VERTEX);
  Ex_orig.Init(*orig_S, TopAbs_VERTEX);
  for (; Ex_orig.More(), Ex.More(); Ex_orig.Next(), Ex.Next())
  {
    k = OCCQueryEngine::instance()->OCCMap->Find(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->UnBind(Ex_orig.Current());
    OCCQueryEngine::instance()->OCCMap->Bind(Ex.Current(), k);
    OCCPoint* point = (OCCPoint*)(OCCQueryEngine::instance()->OccToCGM->find(k))->second;
    point->set_TopoDS_Vertex(TopoDS::Vertex(Ex.Current()));
  }
  return CUBIT_SUCCESS;
}
    


//===============================================================================
// Function   : split_periodic
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::split_periodic( BodySM * /*body_ptr*/,
                                               BodySM *& /*new_body*/ )
{
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_body
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::regularize_body( BodySM * /*body_ptr*/,
                                                   BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_refentity
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus  OCCModifyEngine::regularize_entity( GeometryEntity * /*old_entity_ptr*/,  
                                                      BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : offset_curves
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
CubitStatus OCCModifyEngine::offset_curves( DLIList<Curve*>& curves, 
                                            DLIList<Curve*>& new_curves,
                                            double offset_distance,
                                            const CubitVector& offset_direction, 
                                            int gap_type )
{
  //gap_type has no effect here.
  if( curves.size() == 1 && (offset_direction.x() ||
        offset_direction.y() || offset_direction.z()) )
  {
    gp_Dir offset(offset_direction.x(), offset_direction.y(), offset_direction.z());
    //check if the curve is straight, error out for non-straight curve.
    Curve* curve = curves.get();
    OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
    if(occ_curve->geometry_type() == STRAIGHT_CURVE_TYPE)
    {
      //based on the definitions of the offset curve on OCC, calculate the
      //offset direction to make it consistant with ACIS
      double u1, u2;
      occ_curve->get_param_range(u1, u2);
      CubitVector p1, p2, tangent;
      occ_curve->position_from_u(u1, p1);
      occ_curve->position_from_u(u2, p2);
      tangent = p2 -p1;
      CubitVector cal_offset_dir =  offset_direction * tangent;
      gp_Dir offset(cal_offset_dir.x(), cal_offset_dir.y(), cal_offset_dir.z()); 
      TopoDS_Edge * edge = occ_curve->get_TopoDS_Edge();
      Standard_Real first;
      Standard_Real last;
      Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last);
      Handle(Geom_OffsetCurve) new_curve =
        new Geom_OffsetCurve(myCurve,offset_distance, offset);
      if(!new_curve)
      {
        TopologyEntity *entity = curve->topology_entity();
        BasicTopologyEntity *bte = CAST_TO(entity, BasicTopologyEntity);
        PRINT_ERROR("Can't create offset curve for curve %d.\n", bte->id());
        return CUBIT_FAILURE;
      }
      TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(new_curve, u1, u2);
      Curve* offset_curve = OCCQueryEngine::instance()->populate_topology_bridge(new_edge, CUBIT_TRUE);
      new_curves.append(offset_curve);
      return CUBIT_SUCCESS;
    }
    else
      PRINT_WARNING( "Direction qualifier ignored - only valid for one straight curve\n" );
  }

  else if( offset_direction.x() || offset_direction.y() || offset_direction.z() )
      PRINT_WARNING( "Direction qualifier ignored - only valid for one straight curve\n" );

  else if (offset_direction.x() == 0.0 && offset_direction.y() == 0.0 &&
           offset_direction.z() == 0.0)
  {
    for(int i = 0 ; i < curves.size(); i++)
    {
      Curve* curve = curves.get_and_step();
      OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
      if(occ_curve->geometry_type() == STRAIGHT_CURVE_TYPE)
      {
        // TODO: Remove this condition?
        // A planar nonlinear wire consisting of straight curves can be offset
        // but it isn't clear that the intent was to offset a wire instead
        // of each curve individually.
        PRINT_ERROR("Must have an offset direction for any straight curve.\n");
        return CUBIT_FAILURE;
      }
    }
  }
  //make wire out of ref_edge_list
  BRepBuilderAPI_MakeWire awire;
  TopTools_ListOfShape L;
  OCCCurve* occ_curve = NULL;
  BRepBuilderAPI_Copy shapeCopier;
  for(int i = 0 ; i < curves.size(); i++)
  {
    Curve* curve = curves.get_and_step();
    occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    shapeCopier.Perform(*(occ_curve->get_TopoDS_Edge()));
    L.Append(shapeCopier.Shape());
  }

  awire.Add(L);
  if(awire.Error())
  {
    PRINT_ERROR("Curves must form a planar chain in order for offset to work.\n");
    return CUBIT_FAILURE;
  }
  TopoDS_Wire wire;
  wire = awire.Wire();

#if OCC_VERSION_MAJOR < 6 || (OCC_VERSION_MAJOR == 6 && OCC_VERSION_MINOR < 8)
  CubitBoolean closed = wire.Closed();
  BRepOffsetAPI_MakeOffset offCurveBuilder(wire);
  offCurveBuilder.Perform(offset_distance);
  offCurveBuilder.Build();

  // Prior to OCC version 6.8 / OCE version 0.17
  // For an open wire:
  // It returns a closed wire, with two connecting curves and a positive offset
  // curve list and a negative offset curve list.
  // For a closed wire:
  // It returns a closed wire outside the original one if offset_distance > 0
  // It returns a closed wire inside the original one if offset_distance < 0

  TopoDS_Shape offWireShape = offCurveBuilder.Shape();
  wire = TopoDS::Wire(offWireShape);

  if (!closed)
  {
    int num_curves = countEdges(wire);
    if (num_curves != (2 + 2*curves.size()) )
    {
      PRINT_ERROR("Opencascade can't calculate one offset curve for each"
                  " input curve. It's a limitation, try to use another type"
                  " of curve to offset.\n");
      return CUBIT_FAILURE;
    }

    BRepTools_WireExplorer Ex(wire);
    Ex.Next(); // omit the initial cap/connecting curve
    for (int i = 0; i < curves.size(); ++i)
    {
      TopoDS_Edge new_edge = Ex.Current();
      Curve* curve = curves.get_and_step();
      Curve* offset_curve = OCCQueryEngine::instance()->
          populate_topology_bridge(new_edge, CUBIT_TRUE);
      // double check here to make sure we get the correct offset curve
      // which is if the offset_distance is positive, the new curve should
      // be longer than the original curve.
      if (i == 0)
      {
        double d_offset = offset_curve->measure();
        double d_orig = curve->measure();
        if ((offset_distance > 0 && d_offset < d_orig) ||
            (offset_distance < 0 && d_offset > d_orig))
        {
          for (int j = 0; j < curves.size(); ++j)
          {
            Ex.Next();
          }
          Ex.Next();
        }
        new_edge = Ex.Current();
        offset_curve = OCCQueryEngine::instance()->
            populate_topology_bridge(new_edge, CUBIT_TRUE);
      }
      new_curves.append(offset_curve);
      Ex.Next();
    }

    // successful offset of open wire
    return CUBIT_SUCCESS;
  }

  // if wire is closed
  for (BRepTools_WireExplorer Ex(wire); Ex.More(); Ex.Next())
  {
    TopoDS_Edge new_edge = Ex.Current();
    Curve* offset_curve = OCCQueryEngine::instance()->
        populate_topology_bridge(new_edge, CUBIT_TRUE);
    new_curves.append(offset_curve);
  }

  // successful offset of closed wire
  return CUBIT_SUCCESS;

#else

  // The behavior changed as of OCC 6.8 and OCE 0.17.
  // Not all cases have been tested.
  // Documentation in BRepFill_OffsetWire, which is used by
  // BRepOffsetAPI_MakeOffset, says that an offset will be "to the left"
  // of the spine wire, with each point at the specified offset from the
  // the closest point on the original wire.
  // Negative offsets are "to the right."
  // "Left" and "right" imply a normal direction on the face that contains
  // the wire.
  // The contract for this method in CGM is not documented well, but the
  // previous implementation suggests that a negative offset is intended
  // to return the shorter of the two open offset curves.  That may not be
  // well defined, and is hard to interpret when there are multiple
  // edges in the wire.
  //
  // When testing the offset functionality in a recent release version of
  // Cubit, it seemed that the direction of the offset depends on
  // how the wire turns.  If a closed wire has a positive turning number or
  // an open wire has what might be called a positive turning angle, then
  // the positive offset is to the right.  If the turning angle was
  // negative, then the positive offset was consistently to the left.
  // This behavior is independent of normal direction, since what appears
  // to be right and positive turning angle from one normal direction becomes
  // left and negative turning angle when the normal direction is flipped.

  // The current implementation here returns the result only if the number of
  // offset curves matches the number of input curves for both the
  // positive and negative offset.  In that situation, both the positive
  // and negative offsets are computed, the length of each wire is
  // computed, and the shorter result is returned if the offset is
  // negative.  The longer result is returned if the offset is positive.
  // If the lengths are equal then the curve(s) that OCC computed for the
  // offset with the same sign as the input offset are used.  This is
  // quite different from the behavior of Cubit, but passes the previously
  // existing offset_curves test and agrees with the local in-method
  // documentation that existed prior to OCC 6.8.

  // The WireExplorer is used in this method, which may miss some edges
  // if the result wire does not have all edges connected end to end.

  double positiveOffset = fabs(offset_distance);
  BRepOffsetAPI_MakeOffset posOffCurveBldr(wire);
  // Using an open offset, as in the following line, seems to produce results
  // closer to the original intent of this method, but the original behavior is
  // preserved better with the default Standard_False closed offset.
//  BRepOffsetAPI_MakeOffset posOffCurveBldr(wire, GeomAbs_Arc, Standard_True);
  posOffCurveBldr.Perform(positiveOffset);
  TopoDS_Shape posOffShape = posOffCurveBldr.Shape();
  TopoDS_Wire posOffWire = TopoDS::Wire(posOffShape);
  if (countEdges(posOffWire) != curves.size())
  {
      PRINT_ERROR("The number of positive offset curves computed by the OCE"
                  " geometry engine does not equal the number of curves in"
                  " the input.\n");
      // Remark: It appears to be possible to match curves in the output
      // to curves in the input by calling BRepOffsetAPI_MakeOffset::Generated.
      // One could also directly use BRepFill_OffsetWire.
      // Perhaps compounds curves could be constructed to make the
      // result curves match the input.
      return CUBIT_FAILURE;
  }

  double negativeOffset = -1.0*positiveOffset;
  BRepOffsetAPI_MakeOffset negOffCurveBldr(wire);
  // Using an open offset, as in the following line, seems to produce results
  // closer to the original intent of this method, but the original behavior is
  // preserved better with the default Standard_False closed offset.
//  BRepOffsetAPI_MakeOffset negOffCurveBldr(wire, GeomAbs_Arc, Standard_True);
  negOffCurveBldr.Perform(negativeOffset);
  TopoDS_Shape negOffShape = negOffCurveBldr.Shape();
  TopoDS_Wire negOffWire = TopoDS::Wire(negOffShape);
  if (countEdges(negOffWire) != curves.size())
  {
      PRINT_ERROR("The number of negative offset curves computed by the OCE"
                  " geometry engine does not equal the number of curves in"
                  " the input.\n");
      return CUBIT_FAILURE;
  }

  // compute wire lengths and choose which wire to return
  TopoDS_Wire* wirePtrs[2];
  double wireLengths[2];
  wirePtrs[0] = &posOffWire;
  wirePtrs[1] = &negOffWire;
  wireLengths[0] = wireLengths[1] = 0.0;

  for (int wireIndex = 0; wireIndex < 2; ++wireIndex)
  {
    TopoDS_Wire* wirePtr = wirePtrs[wireIndex];
    double wireLen = 0.0;
    for (BRepTools_WireExplorer Ex(*wirePtr); Ex.More(); Ex.Next())
    {
      GProp_GProps edgeProps;
      BRepGProp::LinearProperties(Ex.Current(), edgeProps);
      wireLen += edgeProps.Mass();
    }
    wireLengths[wireIndex] = wireLen;
  }

  TopoDS_Wire* chosenWirePtr = wirePtrs[0];
  if (offset_distance >= 0)
  {
    if (wireLengths[1] > wireLengths[0])
      chosenWirePtr = wirePtrs[1];
  }
  else if (wireLengths[1] <= wireLengths[0])
    chosenWirePtr = wirePtrs[1];

  for (BRepTools_WireExplorer Ex(*chosenWirePtr); Ex.More(); Ex.Next())
  {
    Curve* offset_curve = OCCQueryEngine::instance()->
        populate_topology_bridge(Ex.Current(), CUBIT_TRUE);
    new_curves.append(offset_curve);
  }
  return CUBIT_SUCCESS;
#endif
}

//=============================================================================
// Function   : countEdges
// Member Type: private
// Description: Count the number of edges encountered by a
//     BRepTools_WireExplorer exploration of the specified wire.
//=============================================================================
int OCCModifyEngine::countEdges(TopoDS_Wire &wire)
{
  int numEdges = 0;
  for (BRepTools_WireExplorer Ex(wire); Ex.More(); Ex.Next())
    ++numEdges;
  return numEdges;
}

//===============================================================================
// Function   : trim_curve
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::trim_curve( Curve* trim_curve, 
                                    const CubitVector& trim_vector,
                                    const CubitVector& keep_vector,
                                    bool keep )
{
  OCCCurve* occ_crv = CAST_TO(trim_curve, OCCCurve);
  if(!occ_crv) 
  {
    PRINT_ERROR("This is not a OCC curve to be trimmed.\n");
    return (Curve*)NULL;
  }
  
  //Determine the trimmed curve's parameter range.
  double u1, u2;
  occ_crv->get_param_range(u1, u2);
  double trim_u = occ_crv->u_from_position(trim_vector);
  double keep_u = occ_crv->u_from_position(keep_vector);
  if(trim_u > u2+TOL || trim_u < u1 - TOL)
  {
    PRINT_ERROR("The trim_vector is outside of the curve range.\n");
    return (Curve*)NULL;
  }
 
  if(keep_u > trim_u )
     u1 =  trim_u;
  else if(keep_u < trim_u)
     u2 = trim_u;

  else
  {
    PRINT_ERROR("Can't determine which part of the curve to be kept.\n");
    return (Curve*)NULL;
  }
  //get the Geom_Curve of the OCCCurve
  TopoDS_Edge * edge = occ_crv->get_TopoDS_Edge();
  Standard_Real first;
  Standard_Real last;
  Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last);
 
  //Trim the curve
  TopoDS_Edge t_edge = BRepBuilderAPI_MakeEdge(myCurve, u1, u2);
  Curve* t_curve = OCCQueryEngine::instance()->populate_topology_bridge(t_edge, CUBIT_TRUE);  
  if(!keep)
  {
    DLIList<OCCLoop*> loops = occ_crv->loops();
    if(loops.size() == 0)
      OCCQueryEngine::instance()->delete_solid_model_entities(trim_curve);
  }
  return t_curve;
}

//===============================================================================
// Function   : create_body_from_surfs
// Member Type: PUBLIC
// Description:
// Author     : Jane Hu
// Date       : 4/22/08
//===============================================================================
CubitStatus OCCModifyEngine::create_solid_bodies_from_surfs(DLIList<Surface*> & ref_face_list,
                                          DLIList<BodySM*>& new_bodies,
                                          bool keep_old,
                                          bool heal,
                                          bool sheet) const
{
  //keep_old and heal are ignored, always delete old.
  //all surfaces should be stand along surface bodies or shell bodies' surface
  CubitStatus stat = CUBIT_SUCCESS;
  if(!sheet)
  {
    Lump* lump = make_Lump(ref_face_list);
    if (!lump)
      return CUBIT_FAILURE;
  
    new_bodies.append(CAST_TO(lump, OCCLump)->get_body());
  }
 
  else
  {
    DLIList<BodySM*> bodies;
    for(int i = 0 ; i < ref_face_list.size(); i++)
    {
      OCCSurface* surf = CAST_TO(ref_face_list.get_and_step(), OCCSurface);
      BodySM* body = surf->my_body();
      if(body)
        bodies.append_unique(body);
    }
    stat = unite(bodies, new_bodies, keep_old); 
  }
  return stat;
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( TBPoint* pt1, 
                                          TBPoint* pt2,
                                          TBPoint* pt3, 
                                          bool full,
                                          bool preview )
{ 
  Curve* new_curve = NULL;
  if(!full)
  {
    CubitVector v2(pt2->coordinates());
    new_curve = const_cast<OCCModifyEngine*> (this)->
                  make_Curve(ARC_CURVE_TYPE,pt1,pt3, &v2);
  }
  else
  {
    CubitVector v1(pt1->coordinates());
    CubitVector v2(pt2->coordinates());
    CubitVector v3(pt3->coordinates());

    gp_Pnt gp_pt1(v1.x(),v1.y(), v1.z());
    gp_Pnt gp_pt2(v2.x(),v2.y(), v2.z());
    gp_Pnt gp_pt3(v3.x(),v3.y(), v3.z());

    Handle(Geom_Circle) curve_ptr;
    curve_ptr = GC_MakeCircle(gp_pt1,gp_pt2,gp_pt3); 

    OCCPoint* occ_pt1 = CAST_TO(const_cast<TBPoint*>(pt1),OCCPoint);
    TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt1);
    new_curve = OCCQueryEngine::instance()->populate_topology_bridge(new_edge,
                                                                   CUBIT_TRUE);
  }
  if(preview)
  {
    GfxPreview::clear();
    OCCCurve* occ_curve = CAST_TO(new_curve, OCCCurve);
    TopoDS_Edge* h_edge = occ_curve->get_TopoDS_Edge();
    // Draw this edge
    OCCDrawTool::instance()->draw_EDGE( h_edge, CUBIT_BLUE_INDEX, CUBIT_TRUE );
     
    OCCQueryEngine::instance()->delete_solid_model_entities(new_curve);
    return (Curve*) NULL;
  }
  else
    return new_curve;
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Curve* curve1, 
                                          Curve* curve2,
                                          Curve* curve3, 
                                          bool full, 
                                          bool preview  )
{ 
  OCCCurve* occ_crv1 = CAST_TO(curve1, OCCCurve);
  OCCCurve* occ_crv2 = CAST_TO(curve2, OCCCurve);
  OCCCurve* occ_crv3 = CAST_TO(curve3, OCCCurve);
  GeometryType type1 = occ_crv1->geometry_type();
  GeometryType type2 = occ_crv2->geometry_type();
  GeometryType type3 = occ_crv3->geometry_type();
  if(type1 != STRAIGHT_CURVE_TYPE || type2 != STRAIGHT_CURVE_TYPE ||
     type3 != STRAIGHT_CURVE_TYPE)
  {
    PRINT_WARNING("Need three straight curves to calculate incenter.\n");
    return (Curve*) NULL;
  } 
  
  //0.check that non of the curves are parallel of each other.
  DLIList<CubitVector> intscts;
  CubitVector vt1, vt2, vt3;
  CubitBoolean none = CUBIT_FALSE;
  OCCQueryEngine::instance()->get_intersections(curve1, curve2, intscts,none,none);
  vt1 = intscts.get();
  intscts.clean_out();
  OCCQueryEngine::instance()->get_intersections(curve2,curve3, intscts,none,none);
  vt2 = intscts.get();
  intscts.clean_out();
  OCCQueryEngine::instance()->get_intersections(curve3, curve1, intscts,none,none);
  vt3 = intscts.get();

  double u11, u12, u21, u22, u31, u32;
  occ_crv1->get_param_range(u11, u12);
  occ_crv2->get_param_range(u21, u22);
  occ_crv3->get_param_range(u31, u32);

  CubitVector tangent1, tangent2, tangent3;
  occ_crv1->get_tangent(vt1, tangent1);
  occ_crv2->get_tangent(vt2, tangent2);
  occ_crv3->get_tangent(vt3, tangent3); 

  CubitVector normal1 = tangent1 * tangent2;
  CubitVector normal2 = tangent2 * tangent3;
  CubitVector normal3 = tangent3 * tangent1;
  if( normal1.length()< TOL || normal2.length()< TOL ||
      normal3.length() < TOL )
  {
    PRINT_WARNING("Three curves must be able to form a triangle.\n");
    return (Curve*) NULL;
  }

  //normals must parallel to each other, meaning all curves must be on
  //the same plane.
  normal1.normalize();
  normal2.normalize();
  normal3.normalize();
  
  CubitVector parallel1 = normal1 * normal2;
  CubitVector parallel2 = normal2 * normal3;
  CubitVector parallel3 = normal3 * normal1;
  if(parallel1.length() > TOL || parallel2.length() > TOL ||
     parallel3.length() > TOL)
  {
    PRINT_WARNING("Three curves must be able to form a triangle.\n");
    return (Curve*) NULL;
  }
  //1.find the angle between each of the two curves
  double angle1, angle2, angle3;
  angle1 = tangent1.interior_angle(tangent2);
  angle2 = tangent2.interior_angle(tangent3);
  angle3 = tangent3.interior_angle(tangent1);

  //2.create curves to bisection each of the angle passing through the
  // vertices of the triangle
  CubitVector t_curve11 = 
         vectorRotate(angle1/2.0, normal1, tangent1);  
  t_curve11.normalize();
  CubitVector p11 = vt1+t_curve11;

  CubitVector t_curve12 = 
         vectorRotate(90.0 - angle1/2.0, -normal1, tangent1);
  t_curve12.normalize();
  CubitVector p12 = vt1 + t_curve12;

  CubitVector t_curve21 =
         vectorRotate(angle2/2.0, normal2, tangent2);
  t_curve21.normalize();
  CubitVector p21 = vt2 + t_curve21;

  CubitVector t_curve22 = 
         vectorRotate(90.0 - angle2/2.0, -normal2, tangent2);
  t_curve22.normalize();
  CubitVector p22 = vt2 + t_curve22;

  CubitVector t_curve31 = 
         vectorRotate(angle3/2.0, normal3, tangent3);
  t_curve31.normalize();

  CubitVector t_curve32 =
         vectorRotate(90.0 - angle3/2.0, -normal3, tangent3);
  t_curve32.normalize();

  //3. find the three intersection points which when connected with the vertices,
  //intersect at same point.
  CubitVector line_p[4], c_p[4], c_ptemp;
  double sc, tc;

  IntersectionTool int_tool;
  for(int i = 0; i < 4; i++)
  {
    if( i == 0)
    {
      line_p[0] = vt1;
      line_p[1] = p11;
      line_p[2] = vt2;
      line_p[3] = p21;
    } 
    else if(i == 1)
      line_p[3] = p22;
    else if(i == 2)
      line_p[1] = p12;
    else
      line_p[3] = p21;  
    int_tool.closest_points_on_segments(line_p[0], line_p[1], line_p[2],
                    line_p[3], c_ptemp, c_p[i], sc, tc);
        
    //check if the closeset point from c_p[i] to three curves are on three curves
    CubitVector closest1, closest2, closest3;
    occ_crv1->closest_point(c_p[i], closest1);
    double u = occ_crv1->u_from_position(closest1);
    if(u > u11-TOL && u < u12 + TOL)
    {
      occ_crv2->closest_point(c_p[i], closest2);
      u = occ_crv2->u_from_position(closest2);
      if(u > u21 - TOL && u < u22 + TOL)
      {
        occ_crv3->closest_point(c_p[i], closest3);
        u = occ_crv3->u_from_position(closest3);
        if(u > u31 - TOL && u < u32 + TOL)
        {
        //4. use the 3 intersection points to find the arc or circle.
          OCCPoint occ_p1(closest1);
          OCCPoint occ_p2(closest2);     
          OCCPoint occ_p3(closest3);
          return create_arc_three(&occ_p1, &occ_p2, &occ_p3, full, preview);
        } 
      }
    }
  }
  PRINT_ERROR("Can't find the tangent points to create circle.\n");
  return (Curve*) NULL;
} 

//===============================================================================
// Function   : create_arc_center_edge
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 12/08
//===============================================================================
Curve* OCCModifyEngine::create_arc_center_edge( TBPoint* pt1, 
                                                TBPoint* pt2,
                                                TBPoint* pt3,
                                                const CubitVector& normal, 
                                                double radius,
                                                bool full,
                                                bool preview ) 
{ 
  CubitVector vec1 = pt1->coordinates(); // Center of arc
  CubitVector vec2 = pt2->coordinates(); // Position on arc
  CubitVector vec3 = pt3->coordinates(); // Position on arc

  CubitVector dir1( vec1, vec2 );
  CubitVector dir2( vec1, vec3 );
  // Re-adjust vec2, vec3 if radius was given
  if( radius != CUBIT_DBL_MAX )
  {
     CubitVector vec;
     vec1.next_point( dir1, radius, vec );
     if(vec.distance_between(vec2) > TOL)
     {
       vec2 = vec;
       pt2 = new OCCPoint(vec);
     }
     vec1.next_point( dir2, radius, vec );
     if(vec.distance_between(vec3) > TOL)
     {
       vec3 = vec;
       pt3 = new OCCPoint(vec);
     }
  }
 
  else
  {
    radius = vec1.distance_between(vec2);
    CubitVector vec;
    vec1.next_point( dir2, radius, vec );
    if(vec.distance_between(vec3) > TOL)
    {
       vec3 = vec;
       pt3 = new OCCPoint(vec);
    }
  }
 
  CubitVector normal_dir = normal;
  if(normal_dir.length() > TOL)
  {
    normal_dir.normalize();
    //verify sense
    if((dir1 * dir2) % normal_dir < 0.0)
    {
      TBPoint* p = pt2;
      pt2 = pt3;
      pt3 = p;
    }
    else if((dir1 * dir2) % normal_dir == 0.0)
    {
      PRINT_ERROR("Normal can't be on the plan of the arc.\n");
      return (Curve*) NULL;
    }
  } 
  else 
  {
    normal_dir = dir1 * dir2;
    normal_dir.normalize();
  }
 
  Handle(Geom_Circle) curve_ptr;
  gp_Dir norm(normal_dir.x(), normal_dir.y(), normal_dir.z());
  gp_Pnt center = gp_Pnt( vec1.x(), vec1.y(), vec1.z());
  curve_ptr = GC_MakeCircle(center,norm,radius);

  OCCPoint* occ_pt1 = CAST_TO(const_cast<TBPoint*>(pt2),OCCPoint);
  TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
  TopoDS_Edge new_edge;
  if(full)
    new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt1);
  
  else
  {
    Handle(Geom_TrimmedCurve) arc;
    gp_Pnt on_arc1 = gp_Pnt( vec2.x(), vec2.y(), vec2.z());
    gp_Pnt on_arc2 = gp_Pnt( vec3.x(), vec3.y(), vec3.z());
    arc = GC_MakeArcOfCircle(curve_ptr->Circ(), on_arc1, on_arc2, Standard_True);
    OCCPoint* occ_pt2 = CAST_TO(const_cast<TBPoint*>(pt3),OCCPoint); 
    TopoDS_Vertex * vt2 = occ_pt2->get_TopoDS_Vertex();
    new_edge = BRepBuilderAPI_MakeEdge(arc, *vt1, *vt2);
  } 

  if(preview)
  {
    GfxPreview::clear();
    // Draw this edge
    OCCDrawTool::instance()->draw_EDGE( &new_edge, CUBIT_BLUE_INDEX, CUBIT_TRUE );
    return (Curve*) NULL;
  }
  return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
}

//===============================================================================
// Function   : create_curve_combine
// Member Type: PUBLIC
// Description: create a curve of combination of several curves.  
// Author     : Jane Hu 
// Date       : 03/09
//===============================================================================

CubitStatus 
OCCModifyEngine::create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr )
{
  int i;

  DLIList<OCCCurve*> occ_curves(curve_list.size());
  CAST_LIST( curve_list, occ_curves, OCCCurve );
  if (curve_list.size() != occ_curves.size())
  {
    PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       Not all input curves are OCC Curves.\n");
    return CUBIT_FAILURE;
  }
  
  BRepBuilderAPI_MakeWire aWire(*(occ_curves.get_and_step()->get_TopoDS_Edge()));
  for(i =1 ; i < curve_list.size(); i++)
  {
    OCCCurve* curve = occ_curves.get_and_step();
    TopoDS_Edge* edge = curve->get_TopoDS_Edge();
      aWire.Add(*edge);
    if(!aWire.IsDone())
    {
      PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The curves are not all connected.\n");
      return CUBIT_FAILURE;
    }
  }
  TopoDS_Wire wire = aWire.Wire(); 
  BRepAdaptor_CompCurve comp_curve(wire);
  GeomAbs_CurveType type = comp_curve.GetType();
  GeomAbs_Shape cont = comp_curve.Continuity();
  if(cont < GeomAbs_G1)
  {
    PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The combined curve is not G1 continued.\n");
    return CUBIT_FAILURE;
  }

  //find the start/end vertices for the combined curve.
  double first_u = comp_curve.FirstParameter();
  double last_u = comp_curve.LastParameter();
  gp_Pnt first = comp_curve.Value(first_u);
  gp_Pnt last = comp_curve.Value(last_u);
  CubitVector first_v(first.X(), first.Y(), first.Z());
  CubitVector last_v(last.X(), last.Y(), last.Z());
  if(first_v.about_equal(last_v))
    comp_curve.SetPeriodic(Standard_True);

  DLIList<CubitVector> v_list;
  v_list.append(first_v);
  v_list.append(last_v);
  v_list.reset();

  DLIList<gp_Pnt*> V_list;
  for(int j = 0; j < 2; j++)
  {
    DLIList<TopologyBridge*> children;
    if (j == 0)
      occ_curves.reset();
    else
      occ_curves.last();
    occ_curves.get()->get_children_virt(children);
    CubitVector v = v_list.get_and_step();
    for(i = 0 ; i < children.size(); i++)
    {
      OCCPoint* vertex = CAST_TO(children.get_and_step(), OCCPoint); 
      CubitVector xyz = vertex->coordinates();
      if(xyz.about_equal(v))
      {
        gp_Pnt p ( v.x(), v.y(), v.z());
        V_list.append(&p);
        break;
      }
    }
  }
   
  V_list.reset();
  TopoDS_Edge topo_edge; 
  gp_Lin line;
  gp_Circ circle;
  gp_Elips ellip;
  gp_Hypr hypr;
  gp_Parab parab;
  Handle_Geom_BezierCurve bezier;
  Handle_Geom_BSplineCurve spline;
  switch(type)
  {
    case GeomAbs_Line:
      line = comp_curve.Line();
      topo_edge = BRepBuilderAPI_MakeEdge(line,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Circle:
      circle = comp_curve.Circle();
      topo_edge = BRepBuilderAPI_MakeEdge(circle,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Ellipse:
      ellip = comp_curve.Ellipse();
      topo_edge = BRepBuilderAPI_MakeEdge(ellip,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Hyperbola:
      hypr = comp_curve.Hyperbola();
      topo_edge = BRepBuilderAPI_MakeEdge(hypr,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_Parabola:
      parab = comp_curve.Parabola();
      topo_edge = BRepBuilderAPI_MakeEdge(parab,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_BezierCurve:
      bezier = comp_curve.Bezier();
      topo_edge = BRepBuilderAPI_MakeEdge(bezier,*V_list.get_and_step(), *V_list.get() );
      break;
    case GeomAbs_BSplineCurve:
      spline = comp_curve.BSpline();
      topo_edge = BRepBuilderAPI_MakeEdge(spline,*V_list.get_and_step(), *V_list.get() );
      break;
    default:
      PRINT_ERROR("In OCCModifyEngine::create_curve_combine\n"
                "       The combined curve is not G1 continued.\n");
      return CUBIT_FAILURE;
  }
  TopoDS_Edge *topo_edge_ptr = new TopoDS_Edge(topo_edge);
  OCCCurve* occ_c = new OCCCurve(topo_edge_ptr);
  new_curve_ptr = occ_c;
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_gqe
// Member Type: PUBLIC
// Description: get the facet geometry query engince instance pointer
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
GeometryQueryEngine *OCCModifyEngine::get_gqe()
{
  return OCCQueryEngine::instance();
}

//===============================================================================
// Function   : is_modify_engine
// Member Type: PUBLIC
// Description: return CUBIT_TRUE if the tb_ptr belongs to this modify engine
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitBoolean OCCModifyEngine::is_modify_engine(const TopologyBridge *tb_ptr) const 
{
  return tb_ptr->get_geometry_query_engine() == OCCQueryEngine::instance();
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* curve1, 
                                              Curve* curve2,
                                              DLIList<CubitVector>& out_list,
                                              double offset,
                                              CubitBoolean ext_first ) 
{
  //offset the curve1 in both directions of normal direction of two curves at
  //center points.
    
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* curve1, 
                                           Surface* face_ptr,
                                           DLIList<CubitVector> & out_list,
                                           double offset,
                                           CubitBoolean ext_surf )
{
  Surface* new_surface = face_ptr;
  if(ext_surf)
    new_surface = make_Surface(face_ptr, CUBIT_TRUE);

  BodySM* bodysm = NULL;
  CubitStatus status = CUBIT_SUCCESS;
  status = create_offset_surface(new_surface, bodysm, offset); 
  if(status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't offset surface. \n");
    return status;
  }
  DLIList<OCCSurface*> surfaces = CAST_TO(bodysm, OCCBody)->my_sheet_surfaces();
  OCCSurface* surface = surfaces.get();
  status = OCCQueryEngine::instance()->get_intersections(curve1, surface, out_list);
  
  if(ext_surf || offset)
    OCCQueryEngine::instance()->delete_solid_model_entities(surface);
 
  //offset surface in opposite direction
  if(!offset)
    return status;

  status = create_offset_surface(new_surface, bodysm, -offset);
  if(status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't offset surface. \n");
    return status;
  }
  surfaces.clean_out();
  surfaces = CAST_TO(bodysm, OCCBody)->my_sheet_surfaces();
  surface = surfaces.get();
  status = OCCQueryEngine::instance()->get_intersections(curve1, surface, out_list);
  OCCQueryEngine::instance()->delete_solid_model_entities(surface);
  if(ext_surf)
    OCCQueryEngine::instance()->delete_solid_model_entities(new_surface);
  return status;
}

//===============================================================================
// Function   : surface_intersection
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu  
// Date       : 03/09
//===============================================================================
CubitStatus OCCModifyEngine::surface_intersection( Surface * surface1,
                                                   Surface * surface2,
                                                   DLIList<Curve*> &intscts,
                                                   const double tol) const
{
  OCCSurface *occ_surface1 =  CAST_TO(surface1, OCCSurface);
  if (occ_surface1 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCSurface *occ_surface2 =  CAST_TO(surface2, OCCSurface);
  if (occ_surface2 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }
  
  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(*(occ_surface1->get_TopoDS_Face()),
                                            *(occ_surface2->get_TopoDS_Face()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curve and surface.\n");
      return CUBIT_FAILURE;
    }

  if (distShapeShape.Value() < tol) //intersect
    {
      int numSol = distShapeShape.NbSolution();
      for (int i = 1; i <= numSol; i++)
        {
          TopoDS_Shape shape = distShapeShape.SupportOnShape1(i);
          if(shape.ShapeType() != TopAbs_EDGE)
            continue;

          TopoDS_Edge* edge = new TopoDS_Edge(TopoDS::Edge(shape));
          OCCCurve* cv = new OCCCurve(edge);
          intscts.append(cv);
        }
    }

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_3_point_plane
// Member Type: PUBLIC
// Description:
// Author     : Jane Hu
// Date       : 01/11
//===============================================================================
CubitStatus OCCModifyEngine::get_3_point_plane( const CubitVector & point_1,
                                            const CubitVector & point_2,
                                            const CubitVector & point_3,
                                            TopoDS_Face*& p_face)const
{
  //Calculate normal of the plane
  CubitVector v1, v2, normal;
  v1 = point_2 - point_1;
  v2 = point_3 - point_1;
  normal = ~(v1 * v2);
  if(fabs(normal.length() - 1) > TOL)
  {
     PRINT_ERROR("The three points are co-linear, and can't be used as a cutting plane.\n");
     return CUBIT_FAILURE;
  }

  gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
  gp_Dir normal_dir(normal.x(), normal.y(), normal.z());
  gp_Pln plane(pt, normal_dir);

  TopoDS_Face face = BRepBuilderAPI_MakeFace(plane);
  p_face = new TopoDS_Face(face);
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : get_mid_plane
// Member Type: PUBLIC
// Description: 
// Author     : Jane Hu
// Date       : 01/09
//===============================================================================
CubitStatus OCCModifyEngine::get_mid_plane( const CubitVector & point_1,
                                            const CubitVector & point_2,
                                            const CubitVector & point_3,
                                            BodySM * body_to_trim_to,
                                            BodySM *& result_body ) const
{
  //Calculate normal of the mid  plane
  CubitVector v1, v2, normal;
  v1 = point_2 - point_1;
  v2 = point_3 - point_1;
  normal = ~(v1 * v2);
  if(fabs(normal.length() - 1) > TOL)
  {
     PRINT_ERROR("The three points are co-linear, and can't be used as a cutting plane.\n");
     return CUBIT_FAILURE;
  }
 
  gp_Pnt pt = gp_Pnt( point_1.x(), point_1.y(), point_1.z());
  gp_Dir normal_dir(normal.x(), normal.y(), normal.z());
  gp_Pln plane(pt, normal_dir);

  TopoDS_Face face = BRepBuilderAPI_MakeFace(plane);
  Surface *surf = OCCQueryEngine::instance()->populate_topology_bridge(face,
                                               CUBIT_TRUE);
  if(!surf)
  {
    PRINT_ERROR("Can't create cutting plane.\n");
    return CUBIT_FAILURE;
  }
  
  BodySM* tool = CAST_TO(surf, OCCSurface)->my_body();
  DLIList<BodySM*> from_bodies;
  from_bodies.append(body_to_trim_to);
  DLIList<BodySM*> midplane_bodies;
  CubitStatus stat = intersect(tool, from_bodies, midplane_bodies, 
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);
  if (midplane_bodies.size() == 1)
    result_body = midplane_bodies.get();
  else {
    for (int i = 0; i < midplane_bodies.size(); ++i)
      OCCQueryEngine::instance()->delete_solid_model_entities(midplane_bodies.get_and_step());
    stat = CUBIT_FAILURE;
  }
  return stat;
}

//=============================================================================
// Function   : get_spheric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 spheric surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_spheric_mid_surface( Surface* surface_ptr1,
                                    Surface* surface_ptr2,
                                    BodySM* body_to_trim_to,
                                    BodySM*& result_body ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != SPHERE_SURFACE_TYPE ||
     occ_surf2->geometry_type() != SPHERE_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be sphere type.\n");
    return CUBIT_FAILURE;
  }
 
  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  gp_Sphere sphere1 = asurface1.Sphere();
  gp_Sphere sphere2 = asurface2.Sphere();

  gp_Pnt center1 = sphere1.Location();
  gp_Pnt center2 = sphere2.Location();

  if(!center1.IsEqual(center2, TOL))
  {
    PRINT_ERROR( "Spheres need to have the same center.\n");
    return CUBIT_FAILURE;
  }

  double radius = sphere1.Radius()/2.0 + sphere2.Radius()/2.0;
  BodySM* tool = sphere(radius);
  CubitVector center(center1.X(), center1.Y(), center1.Z());
  OCCQueryEngine::instance()->translate(tool, center);

  //get the tool surfaces as the tool
  DLIList<Surface*> surfaces;
  tool->surfaces(surfaces);
  assert (surfaces.size() == 1);
  tool = make_BodySM(surfaces.get());

  DLIList<BodySM*> from_bodies, midsurface_bodies;
  from_bodies.append(body_to_trim_to);

  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);
  if (midsurface_bodies.size() == 1)
    result_body = midsurface_bodies.get();
  else {
    for (int i = 0; i < midsurface_bodies.size(); ++i)
      OCCQueryEngine::instance()->delete_solid_model_entities(midsurface_bodies.get_and_step());
    stat = CUBIT_FAILURE;
  }

  return stat;  
}

//=============================================================================
// Function   : get_conic_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 conic surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_conic_mid_surface( Surface* surface_ptr1,
                                    Surface* surface_ptr2,
                                    BodySM* body_to_trim_to,
                                    BodySM*& result_body ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != CONE_SURFACE_TYPE ||
     occ_surf2->geometry_type() != CONE_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be conic type.\n");
    return CUBIT_FAILURE;
  }  

  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  GeomAbs_SurfaceType  type1 = asurface1.GetType();
  GeomAbs_SurfaceType  type2 = asurface2.GetType();
  if(type1 != type2)
  {
    PRINT_ERROR( "Both surfaces provided should be both cylinder or cone type.\n");
    return CUBIT_FAILURE;
  }

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  double height = (bounding_box.diagonal()).length();
  OCCBody* body = CAST_TO(body_to_trim_to, OCCBody);
  CubitVector centroid;
  double volume;
  body->mass_properties(centroid, volume);
  BodySM* tool;
  if(type1 == GeomAbs_Cylinder)
  {
    gp_Cylinder cyl1 = asurface1.Cylinder();
    gp_Cylinder cyl2 = asurface2.Cylinder(); 
    gp_Ax1  axis1 = cyl1.Axis();
    gp_Ax1  axis2 = cyl2.Axis();
    if(!axis1.IsCoaxial(axis2, 0.001, TOL))
    {
      PRINT_ERROR( "Cylinders need to have the same axis of symmetry.\n");
      return CUBIT_FAILURE;
    }
    double radius = cyl1.Radius()/2.0 + cyl2.Radius()/2.0; 
    gp_Ax2 axis;
    axis.SetAxis(axis1);
    TopoDS_Solid S = BRepPrimAPI_MakeCylinder(axis, radius, height);

    Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

    if (lump == NULL)
    {
      PRINT_ERROR("In OCCModifyEngine::get_conic_mid_surface\n"
                  "   Cannot create a cylinder for given radius.\n");
      return CUBIT_FAILURE;
    }

    tool = CAST_TO(lump, OCCLump)->get_body();
    double z = centroid.z();
    z -= height/2.0;
    centroid.z(z);
    OCCQueryEngine::instance()->translate(tool, centroid);  
  }

  else  //GeomAbs_Cone
  {
    gp_Cone cone1 = asurface1.Cone();
    gp_Cone cone2 = asurface2.Cone();
    double angle1 = cone1.SemiAngle();
    double angle2 = cone2.SemiAngle();
    if(fabs(angle1 - angle2) > 0.001)
    {
      PRINT_ERROR( "Cones do not have the same semi-angle.\n");
      return CUBIT_FAILURE;
    }
    gp_Ax1  axis1 = cone1.Axis();
    gp_Ax1  axis2 = cone2.Axis();
    if(!axis1.IsCoaxial(axis2, 0.001, TOL))
    {
      PRINT_ERROR( "Cones need to have the same axis of symmetry.\n");
      return CUBIT_FAILURE;
    } 
    if(axis1.IsOpposite(axis2, 0.001))
    {
      PRINT_ERROR( "Cones need to have the same orientation of axis.\n");
      return CUBIT_FAILURE;
    }
    double r1 = cone1.RefRadius()/2.0 + cone2.RefRadius()/2.0; 
    gp_Ax3 axis;
    axis.SetAxis(axis1);
    gp_Cone cone(axis, angle1, r1);
    TopoDS_Face face = BRepBuilderAPI_MakeFace(cone);
    Surface* surface = 
      OCCQueryEngine::instance()->populate_topology_bridge(face, CUBIT_TRUE);
    tool = CAST_TO(surface,OCCSurface)->my_body();  
  } 
  DLIList<BodySM*> from_bodies, midsurface_bodies;
  from_bodies.append(body_to_trim_to);

  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);
  if (midsurface_bodies.size() == 1)
    result_body = midsurface_bodies.get();
  else {
    for (int i = 0; i < midsurface_bodies.size(); ++i)
      OCCQueryEngine::instance()->delete_solid_model_entities(midsurface_bodies.get_and_step());
    stat = CUBIT_FAILURE;
  }

  return stat;
}

//=============================================================================
// Function   : get_toric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 toric surfaces.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::get_toric_mid_surface( Surface* surface_ptr1,
                                     Surface* surface_ptr2,
                                     BodySM* body_to_trim_to,
                                     BodySM*& result_body ) const
{
  OCCSurface* occ_surf1 = CAST_TO(surface_ptr1, OCCSurface);
  OCCSurface* occ_surf2 = CAST_TO(surface_ptr2, OCCSurface);
  if(occ_surf1->geometry_type() != TORUS_SURFACE_TYPE ||
     occ_surf2->geometry_type() != TORUS_SURFACE_TYPE)
  {
    PRINT_ERROR( "Both surfaces provided should be toric type.\n");
    return CUBIT_FAILURE;
  }

  BRepAdaptor_Surface asurface1(*occ_surf1->get_TopoDS_Face());
  BRepAdaptor_Surface asurface2(*occ_surf2->get_TopoDS_Face());

  gp_Torus torus1 = asurface1.Torus();
  gp_Torus torus2 = asurface2.Torus();

  gp_Pnt center1 = torus1.Location();
  gp_Pnt center2 = torus2.Location();

  if(!center1.IsEqual(center2, TOL))
  {
    PRINT_ERROR( "Torii need to have the same center.\n");
    return CUBIT_FAILURE;
  }

  double major_r1 = torus1.MajorRadius();
  double major_r2 = torus2.MajorRadius();
  if(fabs(major_r1 - major_r2) > TOL)
  {
    PRINT_ERROR( "Torii need to have the same major radius.\n");
    return CUBIT_FAILURE;
  }

  gp_Ax1 axis1 = torus1.Axis();
  gp_Ax1 axis2 = torus2.Axis();
  if(!axis1.IsCoaxial(axis2, 0.001, TOL))
  {
    PRINT_ERROR( "Torii need to have the same axis of symmetry.\n");
    return CUBIT_FAILURE;
  }

  double radius = torus1.MinorRadius()/2.0 + torus2.MinorRadius()/2.0;
  gp_Ax2 axis;
  axis.SetAxis(axis1);
  TopoDS_Solid S = BRepPrimAPI_MakeTorus(axis, major_r1, radius);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::get_toric_mid_surface\n"
                "   Cannot create a torus for given radii.\n");
    return CUBIT_FAILURE;
  }

  BodySM* tool = CAST_TO(lump, OCCLump)->get_body();

  DLIList<BodySM*> from_bodies, midsurface_bodies;
  from_bodies.append(body_to_trim_to);

  CubitStatus stat = intersect(tool, from_bodies, midsurface_bodies,
                               CUBIT_TRUE);
  OCCQueryEngine::instance()->delete_solid_model_entities(tool);
  if (midsurface_bodies.size() == 1)
    result_body = midsurface_bodies.get();
  else {
    for (int i = 0; i < midsurface_bodies.size(); ++i)
      OCCQueryEngine::instance()->delete_solid_model_entities(midsurface_bodies.get_and_step());
    stat = CUBIT_FAILURE;
  }

  return stat; 
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : Jane Hu
// Date       : 03/2009
//=============================================================================
CubitStatus OCCModifyEngine::tweak_chamfer( DLIList<Curve*> & curve_list, 
                                            double left_offset,
                                            DLIList<BodySM*> & new_bodysm_list,
                                            double right_offset,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
  CubitStatus stat;
  int count = 0;
  if(right_offset <= 0.0)
    right_offset = left_offset;

  for(int i = 0; i < curve_list.size(); i++)
  {
    BodySM * new_bodysm_ptr = NULL;
    stat = tweak_fillet(curve_list.get_and_step(), left_offset, right_offset,
                      new_bodysm_ptr , keep_old_body, CUBIT_FALSE, CUBIT_FALSE);
    if(stat && new_bodysm_ptr)
    {
      new_bodysm_list.append_unique(new_bodysm_ptr);
      count = new_bodysm_list.size();
    }
    else
      break;
  }

  if(count == 0)
    return CUBIT_FAILURE;

  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape; 
      CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape(modified_shape);
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
      }
    }
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body 
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If 
//              multiple vertices are supplied, only one offset value is 
//              allowed and the edges are not used.
// Author     : Jane Hu
// Date       : 03/09
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_chamfer( DLIList<TBPoint*> & point_list, 
                                double offset1,
                                DLIList<BodySM*> & new_bodysm_list,
                                Curve * edge1,
                                double offset2,
                                Curve * edge2,
                                double offset3,
                                Curve * edge3,
                                CubitBoolean keep_old_body,
                                CubitBoolean preview ) const
{
  // Sort out vertices between sheet and solid bodies
  DLIList<TBPoint*> solid_points, sheet_points;
  DLIList<OCCSurface*> s_list;
  DLIList<OCCBody*> bodies;
  if( sort_points_by_body_type( point_list, solid_points, sheet_points, 
                                s_list, bodies ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if(solid_points.size() > 0 && solid_points.size() != bodies.size())
  {
    PRINT_ERROR( "cannot find bodies corresponding to the points.\n" );
    return CUBIT_FAILURE;
  }

  if(sheet_points.size() > 0 && sheet_points.size() != s_list.size())
  {
    PRINT_ERROR( "cannot find surfaces corresponding to the points.\n" );
    return CUBIT_FAILURE;
  }

  // Do simple forms
  if( edge1 == NULL || offset2 <= 0.0 )
  {
    if( tweak_chamfer_solid( solid_points, bodies, offset1, new_bodysm_list,
      keep_old_body, preview )== CUBIT_FAILURE )
      return CUBIT_FAILURE;
    return tweak_fillet_chamfer_sheet( sheet_points, s_list, offset1, 
           CUBIT_FALSE, new_bodysm_list, keep_old_body, preview );
  }

  if( solid_points.size() > 1 || sheet_points.size() > 1 )
  {
    PRINT_ERROR( "cannot chamfer multiple vertices with a variable radius.\n" );
    return CUBIT_FAILURE;
  }

  if( solid_points.size() )
  {
    TBPoint *point_ptr = solid_points.get();

    if( tweak_chamfer_solid( point_ptr, bodies.get(), offset1, edge1, 
        offset2, edge2, offset3, edge3,
        new_bodysm_list, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }

  if( sheet_points.size() )
  {
    TBPoint *point_ptr = sheet_points.get();

    if( tweak_chamfer_sheet( point_ptr, s_list.get(), offset1, edge1, offset2, 
        edge2, new_bodysm_list, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }
  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::tweak_chamfer_solid( DLIList<TBPoint*> &point_list,
                                    DLIList<OCCBody*> &bodies,
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )const
{
  //if point_list.size() > 1, after the first chamfer operation, the rest of 
  //the points in point_list maybe deleted and recreated, so can't use the 
  //points' pointers any more.
  DLIList<CubitVector> p_locs;
  if(point_list.size() > 1)
  {
    for (int i = 1; i < point_list.size(); i++)
     p_locs.append(point_list[i]->coordinates());
  }
  for(int i = 0; i < point_list.size(); i++)
  {
    DLIList<TopologyBridge*> parents;
    TBPoint* point = point_list.get_and_step();
    OCCBody* body = bodies.get_and_step();
    DLIList<OCCPoint*> new_p_list;
    if( i > 0)
    {
      body->get_all_points(new_p_list);
      if(!new_p_list.move_to((OCCPoint*)point)) 
      {
        CubitVector p_loc = p_locs[i-1];
        for(int j = 0; j < new_p_list.size(); j++)
        {
          CubitVector test_loc = new_p_list.step_and_get()->coordinates();
          if(test_loc == p_loc)
          {
            point = new_p_list.get(); 
            break;
          }
        }
      }
    }
    OCCPoint* occ_point = CAST_TO(point, OCCPoint);
    if(occ_point != NULL)
      occ_point->get_parents_virt(parents); //OCCCurves
    assert(parents.size() == 3);
    DLIList<Curve*> curves;
    for(int j = 0; j < 3; j++)
    {
      OCCCurve* occ_curve = CAST_TO(parents.get_and_step(), OCCCurve); 
      curves.append(occ_curve);
    }
    CubitStatus stat;
    stat = tweak_chamfer_solid(point, body,radius, curves.pop(), radius, 
                        curves.pop(), radius, curves.pop(), new_bodysm_list,
                        keep_old_body, CUBIT_FALSE );
    if(!stat)
      return CUBIT_FAILURE;
  }
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape ;
      CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape(modified_shape);
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
      }
    }
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }

  return CUBIT_SUCCESS;

}

CubitStatus
OCCModifyEngine::tweak_chamfer_solid( TBPoint* point_ptr,
                                    OCCBody* body,
                                    double r1,
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    double r3,
                                    Curve *c3,
                                    DLIList<BodySM *> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )const
{
  if(r1 <= 0.0 || r2 <= 0.0 || r3 <= 0.0)
  {
    PRINT_ERROR( "Chamfer radii must be greater than zero.\n" );
    return CUBIT_FAILURE;
  }
  
  DLIList<Curve*> curves;
  curves.append(c1);
  curves.append(c2);
  curves.append(c3);
  
  DLIList<double> radii;
  radii.append(r1);
  radii.append(r2);
  radii.append(r3);

  //check point on curves
  OCCPoint* occ_point = CAST_TO(point_ptr, OCCPoint);
  CubitVector position = occ_point->coordinates();
  DLIList<CubitVector> locations;
  for(int i = 0; i < 3; i++)
  {
    OCCCurve *occ_curve = NULL;
    occ_curve = CAST_TO(curves.get_and_step(), OCCCurve);
    double length = occ_curve->measure();

    DLIList<OCCPoint*> point_list;
    occ_curve->get_points(point_list);
    CubitBoolean in = point_list.is_in_list(occ_point); 
    if(!in)
    {
      PRINT_ERROR( "Point is not on curve.\n" );
      return CUBIT_FAILURE;
    }
    //find cutting points on curves
    double u, u1, u2;
    occ_curve->get_param_range(u1,u2);
    u = occ_curve->u_from_position(position); 
    if(fabs(u-u1) < TOL)
      u = occ_curve->u_from_arc_length(u1, radii[i]);
    else
      u = occ_curve->u_from_arc_length(u1, length-radii[i]);
    CubitVector c_p;
    occ_curve->position_from_u(u, c_p);
    locations.append(c_p);
  }

  //decide normal
  CubitVector v1, v2, normal;
  CubitVector point_1 = locations.pop();
  CubitVector point_2 = locations.pop();
  CubitVector point_3 = locations.pop();
  v1 = point_2 - point_1;
  v2 = point_1 - point_3;
  normal = ~(v1 * v2); 
  CubitVector center;
  double volume;
  body->mass_properties(center, volume); 
  CubitVector dir = ~(center - position);
  if(normal % dir > 0.0)//1, 3, 2 order
  {
    CubitVector v = point_2;
    point_2 = point_3;
    point_3 = v;
  }

  DLIList<BodySM*> bodies;
  BodySM* new_body;
  if(keep_old_body || preview)
    new_body = copy_body(body); 
  else
    new_body = body;

  bodies.append(new_body);
  const CubitVector p1 = point_1;
  const CubitVector p2 = point_2;
  const CubitVector p3 = point_3;
  CubitStatus status = const_cast<OCCModifyEngine*> (this)->
                               section(bodies, p1, p2, p3, 
                               new_bodysm_list, true,false, false);    
  if(!status)
    return CUBIT_FAILURE;

  if(!preview)
    return CUBIT_SUCCESS;

  GfxPreview::clear();

  for(int i = 0; i < new_bodysm_list.size(); i++)
  {
     BodySM* new_bodysm = new_bodysm_list.get_and_step();
     DLIList<OCCSurface*> surfs = CAST_TO(new_bodysm, OCCBody)->my_sheet_surfaces();
     for(int j = 0 ; j < surfs.size(); j++)
     {
       TopoDS_Face* modified_shape = surfs.get_and_step()->get_TopoDS_Face();
       // Draw this face
       OCCDrawTool::instance()->draw_FACE( modified_shape, CUBIT_BLUE_INDEX, CUBIT_TRUE );
     }
  }
  OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
  new_bodysm_list.clean_out();

  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::sort_points_by_body_type( DLIList<TBPoint*> &point_list,
                                         DLIList<TBPoint*> &solid_points,
                                         DLIList<TBPoint*> &sheet_points,
                                         DLIList<OCCSurface*> &s_list,
                                         DLIList<OCCBody*> &bodies )const
{
  for (int i = 0; i < point_list.size(); i++)
  {
    DLIList<TopologyBridge*> parents;
    OCCPoint* point = CAST_TO(point_list.get_and_step(), OCCPoint);
    int curve_size = 0;
    if(point != NULL)
    {
      point->get_parents_virt(parents); //OCCCurves
      if(parents.size() < 2)
      {
        PRINT_ERROR( "Vertex found not attached to any surfaces.\n" );
        return CUBIT_FAILURE;
      } 
      else if(parents.size() > 3)
      {
        PRINT_ERROR( "Vertex found attached to multiple bodies.\n" );
        return CUBIT_FAILURE;
      }
      curve_size = parents.size();
    }
    
    OCCCurve* occ_curve = CAST_TO(parents.get(), OCCCurve);
    parents.clean_out();
    occ_curve->get_parents_virt(parents); //OCCCoEdges
    if(parents.size() == 0)
    {
      PRINT_ERROR( "Vertex found not attached to any surfaces.\n" );
      return CUBIT_FAILURE;
    }
    OCCCoEdge* coedge = CAST_TO(parents.get(), OCCCoEdge);
    parents.clean_out();
    coedge->get_parents_virt(parents);  //OCCLoops
    assert(parents.size() > 0);
    OCCLoop* loop = CAST_TO(parents.get(), OCCLoop);
    parents.clean_out();
    loop->get_parents_virt(parents); //OCCSurface
    assert(parents.size() > 0);
    OCCSurface* s = CAST_TO(parents.get(), OCCSurface); 
    if(s->my_body() != NULL && curve_size == 2) //sheet body
    {
      s_list.append(s);
      sheet_points.append(point);
    }
    else if(s->my_body() != NULL && curve_size == 3) //shell body
    {
      PRINT_ERROR( "Vertex found attached to multiple surfaces but not on bodies.\n" );
      return CUBIT_FAILURE;
    }
    else
    {
      solid_points.append(point);
      DLIList<OCCBody*> solid_bodies;
      s->get_bodies(solid_bodies);
      assert(solid_bodies.size() == 1);
      bodies += solid_bodies;
    }
  }
  return CUBIT_SUCCESS;
}
//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid 
//              bodies.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( DLIList<Curve*> & curve_list, 
                                           double radius,
                                           DLIList<BodySM*> & new_bodysm_list,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview )const 
{
  CubitStatus stat;
  int count = 0;
  for(int i = 0; i < curve_list.size(); i++) 
  {
    BodySM * new_bodysm_ptr = NULL;
    stat = tweak_fillet(curve_list.get_and_step(), radius, radius, 
                        new_bodysm_ptr , keep_old_body, CUBIT_FALSE);
    if(stat && new_bodysm_ptr)
    {
      new_bodysm_list.append_unique(new_bodysm_ptr);
      count = new_bodysm_list.size();
    }
    else
      break;
  }

  if(count == 0) 
    return CUBIT_FAILURE;
 
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape ; 
      CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape(modified_shape);
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
      }
    }
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);  
    new_bodysm_list.clean_out();
  }
  
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid
//              body.  The fillet/chamfer has a variable radius from the
//              start to the end of the curve.
// Author     : Jane Hu
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * curve_ptr,
                                           double start_radius,
                                           double end_radius,
                                           BodySM *& new_bodysm_ptr,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview )const
{
  return tweak_fillet(curve_ptr, start_radius, end_radius, new_bodysm_ptr,
                      keep_old_body, preview, CUBIT_TRUE);
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: private
// Description: Create a round fillet (or blend) at the given curves on a solid 
//              body.  The fillet/chamfer has a variable radius from the 
//              start to the end of the curve.
// Author     : Jane Hu 
// Date       : 01/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * curve_ptr, 
                                           double start_radius,
                                           double end_radius,
                                           BodySM *& new_bodysm_ptr,
                                           CubitBoolean keep_old_body,
                                           CubitBoolean preview,
                                           CubitBoolean if_fillet ) const
{
  //check if this id is valid 
  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  DLIList <OCCBody* > *bodies = oqe->BodyList;  
  DLIList<OCCCurve*> curves;
  for(int j = 0; j <  bodies->size(); j++)
  {
    OCCBody* body = bodies->get_and_step();
    body->get_all_curves(curves);
  }

  bool curve_alive = false;
  for(int j = 0; j <  curves.size(); j++)
  {
    if(curve_ptr == curves.get_and_step())
    {
      curve_alive = true; 
      break;
    }
  }

  if(!curve_alive)
  {
    PRINT_ERROR("This curve is not valid in the current model.\n");
    return CUBIT_FAILURE;
  }

  OCCCurve *occ_curve = CAST_TO(curve_ptr, OCCCurve);
  TopoDS_Edge* topo_edge = occ_curve->get_TopoDS_Edge();

  TopTools_IndexedDataMapOfShapeListOfShape M;
  DLIList<TopoDS_Shape*> shape_list;
  TopoDS_Face* s = NULL;
  for(int j = 0; j <  bodies->size(); j++)
  {
    OCCBody* body = bodies->get_and_step();
    TopExp_Explorer Ex;
    TopoDS_Shape* pshape ; 
    body->get_TopoDS_Shape(pshape);
    
    if (pshape && !pshape->IsNull())
    {
      M.Clear();
      TopExp::MapShapesAndAncestors(*pshape, TopAbs_EDGE, TopAbs_COMPOUND, M);
      if(!M.Contains(*topo_edge))
         continue;

      shape_list.append_unique(pshape);
    }   

    if(!if_fillet) //for chamfer, need to know the face for the curve.
    {
      DLIList<TopologyBridge*> parents;
      occ_curve->get_parents_virt(parents); //OCCCoEdges
      assert(parents.size() > 1);
      OCCCoEdge* coedge = CAST_TO(parents.get(), OCCCoEdge);
      parents.clean_out();
      coedge->get_parents_virt(parents);  //OCCLoops
      assert(parents.size() > 0);
      OCCLoop* loop = CAST_TO(parents.get(), OCCLoop);
      parents.clean_out(); 
      loop->get_parents_virt(parents); //OCCSurface
      assert(parents.size() > 0);
      s = CAST_TO(parents.get(), OCCSurface)->get_TopoDS_Face();
    }
  }
  if(shape_list.size() != 1)
  {
    PRINT_ERROR("Fillets must be created on solids.\n");
    return CUBIT_FAILURE;
  }

  TopoDS_Shape newShape;
  TopoDS_Shape* shape = shape_list.get();
  if(keep_old_body)
  {
    BRepBuilderAPI_Copy api_copy(*shape);
    newShape = api_copy.ModifiedShape(*shape);
  }
  else
    newShape = *shape;

  BRepBuilderAPI_MakeShape* fillet;
  if(if_fillet)
  {
    fillet = new BRepFilletAPI_MakeFillet(newShape);
    dynamic_cast<BRepFilletAPI_MakeFillet*>(fillet)->Add(start_radius, end_radius, *topo_edge);
  }
  else
  {
    fillet = new BRepFilletAPI_MakeChamfer(newShape);
    dynamic_cast<BRepFilletAPI_MakeChamfer*>(fillet)->Add(start_radius, end_radius, *topo_edge, *s);
  }
  fillet->Build();

  if(!fillet->IsDone())
  {
    PRINT_ERROR("Can't create fillet on given curve.\n");
    return CUBIT_FAILURE;
  } 
  TopoDS_Shape modified_shape = fillet->Shape();

  if( !preview )
  {
    TopExp_Explorer Ex;
    Ex.Init(newShape, TopAbs_SOLID);
    TopoDS_Solid old_solid = TopoDS::Solid(Ex.Current());
    OCCLump::update_OCC_entity(old_solid , modified_shape, fillet);     
    DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->populate_topology_bridge(modified_shape);
    new_bodysm_ptr = CAST_TO(tbs.get(), BodySM);  
  }
  else
  {
    GfxPreview::clear();

    TopExp_Explorer Ex;
    Ex.Init(modified_shape, TopAbs_FACE); 
    for( ; Ex.More(); Ex.Next() )
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      // Draw this face
      OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
    }
  }
  delete fillet;
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Jane Hu
// Date       : 01/09 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_fillet( DLIList<TBPoint*> & ref_vertex_list, 
                               double radius,
                               DLIList<BodySM*> & new_bodysm_list,
                               CubitBoolean keep_old_body,
                               CubitBoolean preview ) const
{
  DLIList<OCCSurface*> s_list;
  return tweak_fillet_chamfer_sheet(ref_vertex_list, s_list, radius, CUBIT_TRUE,
         new_bodysm_list, keep_old_body, preview);
}

CubitStatus
OCCModifyEngine::tweak_fillet_chamfer_sheet( DLIList<TBPoint*> & ref_vertex_list,
                               DLIList<OCCSurface*> faces,
                               double radius,
                               CubitBoolean is_fillet,
                               DLIList<BodySM*> & new_bodysm_list,
                               CubitBoolean keep_old_body,
                               CubitBoolean preview )const
{
  TopTools_IndexedDataMapOfShapeListOfShape M;

  for(int i = 0; i < ref_vertex_list.size(); i ++)
  {
    TBPoint* pnt = ref_vertex_list.get_and_step();
    OCCPoint* occ_pnt = CAST_TO(pnt, OCCPoint);
    TopoDS_Vertex* vertex = occ_pnt->get_TopoDS_Vertex();
    OCCSurface* face = NULL;

    if( faces.size() == 0)
    {
      OCCQueryEngine* oqe = OCCQueryEngine::instance();
      //make sure the vertex is on sheet body, not on a volume.
      DLIList <OCCBody*> *bodies = oqe->BodyList;
      for(int k =0 ; k < bodies->size(); k++)
      {
        OCCBody* occ_body = bodies->get_and_step();
        TopExp_Explorer Ex;
        TopoDS_Shape* pshape ; 
        occ_body->get_TopoDS_Shape(pshape);
        if(pshape && !pshape->IsNull())
        {
          M.Clear();
          TopExp::MapShapesAndAncestors(*pshape, TopAbs_VERTEX, TopAbs_COMPOUND, M);
          if(M.Contains(*vertex) && !occ_body->is_sheet_body())
          {
            PRINT_ERROR("Fillet on vertex can only be done on sheet body.\n");
            return CUBIT_FAILURE;
          }
          else if(M.Contains(*vertex) && occ_body->is_sheet_body())
          {
            DLIList<OCCSurface*> surfaces = occ_body->my_sheet_surfaces();
            face = surfaces.get();
          }
        }
      } 
      //find corresponding faces.
      if(face == NULL) //in case there's a bug in code.
      {
        DLIList <OCCSurface* > *surfaces = oqe->SurfaceList;
        for(int k =0 ; k < surfaces->size(); k++)
        {
          OCCSurface* occ_face = surfaces->get_and_step();
          TopoDS_Face* topo_face = occ_face->get_TopoDS_Face();
          TopExp_Explorer Ex;
          M.Clear();
          TopExp::MapShapesAndAncestors(*topo_face, TopAbs_VERTEX, TopAbs_FACE, M);
          if(!M.Contains(*vertex))
            continue;
          face = occ_face; 
          break;
        }
      }
    }
    else
      face = faces.get_and_step();

    if(face == NULL)
    {
      PRINT_ERROR("Can't find corresponding surface for the vertex.\n");
      return CUBIT_FAILURE;
    }

    if(!is_fillet)
    {
      //find the two edges sharing the vertex.
      DLIList<OCCCurve*> curves;
      face->get_curves(curves);
      int size = curves.size();
      for(int j = 0; j < size; j ++)
      {
        DLIList<OCCPoint*> point_list;   
        OCCCurve *curve = curves.get();
        curve->get_points(point_list);
        if(!point_list.is_in_list(occ_pnt))
          curves.remove();
        else
          curves.step();
      }
      assert (curves.size()==2);
      tweak_chamfer_sheet(pnt, face, radius, curves.pop(), radius, curves.pop(),
                        new_bodysm_list, keep_old_body, CUBIT_FALSE); 
    }

    else
    {
      TopoDS_Face *shape = face->get_TopoDS_Face();
      TopoDS_Face newShape;
      if(keep_old_body)
      {
        BRepBuilderAPI_Copy api_copy(*shape);
        newShape = TopoDS::Face(api_copy.ModifiedShape(*shape));
      }
      else
        newShape = *shape;

      BRepFilletAPI_MakeFillet2d fillet(newShape);
      TopoDS_Edge fillet_edge;
      fillet_edge = fillet.AddFillet(*vertex, radius);
      fillet.Build();
      if(fillet.Status() != ChFi2d_IsDone)
      {
        PRINT_ERROR("Can't create fillet on given curve.\n");
        return CUBIT_FAILURE;
      }
      TopoDS_Shape modified_shape = fillet.Shape();

      TopExp_Explorer Ex;
      OCCSurface::update_OCC_entity(newShape , modified_shape, &fillet, vertex);
      TopoDS_Face modified_face = TopoDS::Face(modified_shape);
      Surface* surf = OCCQueryEngine::instance()->populate_topology_bridge(modified_face, CUBIT_TRUE);
      BodySM* new_bodysm_ptr = CAST_TO(surf,OCCSurface)->my_body();
      new_bodysm_list.append_unique(new_bodysm_ptr);
    }
  }

  if(!preview )
    return CUBIT_SUCCESS;

  GfxPreview::clear();

  for(int i = 0; i < new_bodysm_list.size(); i++)
  {
     BodySM* new_bodysm = new_bodysm_list.get_and_step();
     DLIList<OCCSurface*> surfs = CAST_TO(new_bodysm, OCCBody)->my_sheet_surfaces();
     for(int j = 0; j < surfs.size(); j ++)
     {
       TopoDS_Face* modified_shape = surfs.get_and_step()->get_TopoDS_Face();
       // Draw this face
       OCCDrawTool::instance()->draw_FACE( modified_shape, CUBIT_BLUE_INDEX, CUBIT_TRUE );
     }
  }
  OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
  new_bodysm_list.clean_out();

  return CUBIT_SUCCESS;
}

CubitStatus
OCCModifyEngine::tweak_chamfer_sheet(TBPoint* pnt,
                                     OCCSurface* face,
                                     double d1,
                                     Curve* edge1,
                                     double d2,
                                     Curve* edge2,
                                     DLIList<BodySM*> & new_bodysm_list,
                                     CubitBoolean keep_old_body,
                                     CubitBoolean preview ) const
{
  TopoDS_Face *shape = face->get_TopoDS_Face();
  TopoDS_Face newShape;
  if(keep_old_body)
  {
    BRepBuilderAPI_Copy api_copy(*shape);
    newShape = TopoDS::Face(api_copy.ModifiedShape(*shape));
  }
  else
    newShape = *shape;

  BRepFilletAPI_MakeFillet2d fillet(newShape);

  TopoDS_Edge fillet_edge;
  if(edge1 == NULL || edge2 == NULL)
  {
    PRINT_ERROR("Cannot find the two edges for the vertex.\n");
    return CUBIT_FAILURE;
  }
  TopoDS_Edge* topo_e1 = CAST_TO(edge1, OCCCurve)->get_TopoDS_Edge();
  TopoDS_Edge* topo_e2 = CAST_TO(edge2, OCCCurve)->get_TopoDS_Edge();
  TopoDS_Vertex common_v;
  TopExp::CommonVertex(*topo_e1, *topo_e2, common_v);
  fillet_edge = fillet.AddChamfer( *topo_e1, *topo_e2, d1, d2);

  fillet.Build() ;
  if(fillet.Status() != ChFi2d_IsDone)
  {
    PRINT_ERROR("Can't create chamfer on given vertex.\n");
    return CUBIT_FAILURE;
  }
  TopoDS_Shape modified_shape = fillet.Shape();

  if( !preview )
  {
    TopExp_Explorer Ex;
    Ex.Init(newShape, TopAbs_FACE);
    TopoDS_Face old_face = TopoDS::Face(Ex.Current());
    OCCSurface::update_OCC_entity(old_face , modified_shape, &fillet, &common_v);
    DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->populate_topology_bridge(modified_shape);
    BodySM* new_bodysm_ptr = CAST_TO(tbs.get(), BodySM);
    new_bodysm_list.append_unique(new_bodysm_ptr);
  }
  else
  {
    GfxPreview::clear();

    TopExp_Explorer Ex;
    Ex.Init(modified_shape, TopAbs_FACE);
    for( ; Ex.More(); Ex.Next() )
    {
      TopoDS_Face face = TopoDS::Face(Ex.Current());
      // Draw this face
      OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
    }
  }
  return CUBIT_SUCCESS;
}
//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Jane Hu
// Date       : 04/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Surface*> & surface_list, 
                                         const CubitVector & delta,
                                         DLIList<BodySM*> & new_bodysm_list, 
                                         CubitBoolean keep_old_body ,
                                         CubitBoolean preview) const
{
  CubitStatus stat;
  for(int i = 0 ; i < surface_list.size(); i++)
  {
    Surface* surf = surface_list.get_and_step();
    OCCSurface* occ_surf = CAST_TO(surf, OCCSurface);
    if(!occ_surf)
      continue;

    BodySM* original_body = occ_surf->my_body();
    if(original_body == NULL)
    {
      DLIList<OCCBody*> original_bodies;
      occ_surf->get_bodies(original_bodies);
      if(original_bodies.size() > 1)
      {
        PRINT_ERROR( "Cannot tweak move the surface in non-mainfold solids. \n");
        return CUBIT_FAILURE;
      }
      else if(original_bodies.size() == 0)
      {
        PRINT_ERROR( "Interal error: Can't find associated solid. \n");
        return CUBIT_FAILURE;
      }
      original_body = original_bodies.get();
      assert(original_body != NULL);
    }

    //check to make sure that the surf is not on a sheet body. 
    OCCLump* lump = occ_surf->my_lump(); 
    if(lump != NULL && (lump->my_sheet_surface() || lump->my_shell()))
    {
      PRINT_ERROR( "Cannot tweak move surfaces that are not in a solid\n");
      return CUBIT_FAILURE;
    }
    DLIList<GeometryEntity*> ref_ent_list;
    ref_ent_list.append(occ_surf);
    DLIList<BodySM*> result_bodies;
    stat = sweep_translational(ref_ent_list, result_bodies, delta, 0.0, 1,
                               false, true, false, true); 
    if(stat == CUBIT_FAILURE)
    {
      PRINT_ERROR( "Cannot tweak move the surface. \n");
      return CUBIT_FAILURE;
    }
    assert(result_bodies.size() == 1);
    
    //determine if the delta is to trim the existing body or extend it.
    CubitVector center_point;
    center_point = occ_surf->center_point();
    center_point += ~delta;
     
    CubitBoolean trim = CUBIT_FALSE;
    if(original_body->point_containment(center_point) == CUBIT_PNT_INSIDE)
       trim = CUBIT_TRUE;

    //subtract or unite the two bodies
    if(!trim)
    {
      result_bodies.insert_first(original_body);
      stat = unite(result_bodies, new_bodysm_list, keep_old_body);
    }
    else
    {
      DLIList<BodySM*> from_bodies;
      from_bodies.append(original_body);
      stat = subtract(result_bodies, from_bodies, new_bodysm_list, CUBIT_FALSE,
                      keep_old_body);
    }
    if(stat == CUBIT_FAILURE)
      return CUBIT_FAILURE;
  } 
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape ;
      CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape(modified_shape);
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
      }
    }
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Jane Hu
// Date       : 04/09
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Curve*> & curves,
                                         const CubitVector & delta,
                                         DLIList<BodySM*> & new_bodysm_list, 
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview ) const
{
  gp_Dir offset_dir(delta.x(), delta.y(), delta.z());
  double length = delta.length();
  
  for(int i = 0 ; i < curves.size(); i++)
  {
    Curve* curve = curves.get_and_step();
    OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
    if(!occ_curve)
      continue;
    //check to make sure that the curve is on a sheet body.
    DLIList<OCCLoop*> loops;
    loops = occ_curve->loops();
    if(loops.size() == 0)
    {
      PRINT_ERROR( "Cannot tweak move curves that are free\n");
      return CUBIT_FAILURE;
    }
    else if(loops.size() != 1)
    {
      PRINT_ERROR( "Can only tweak move curves attached to one surface\n");
      return CUBIT_FAILURE;
    }
    //determine if the delta is to trim the existing surface or extend it.
    double u_low, u_upper;
    CubitVector a_point;
    occ_curve->get_param_range(u_low, u_upper);
    occ_curve->position_from_u((u_low + u_upper)/2, a_point);
    a_point += ~delta;
    DLIList<TopologyBridge*> parents;
    loops.get()->get_parents_virt(parents);
    OCCSurface* surface = CAST_TO(parents.get(), OCCSurface);
    BodySM* original_body = surface->my_body();
    if(!original_body)
      original_body = surface->my_shell()->my_body();
    if(!original_body)
    {
      PRINT_ERROR("Can't tweak move curves on volumes.\n");
      return CUBIT_FAILURE;
    }
    CubitBoolean trim = CUBIT_FALSE;
    if(surface->point_containment(a_point) == CUBIT_PNT_INSIDE)
       trim = CUBIT_TRUE;

    TopoDS_Edge * edge = occ_curve->get_TopoDS_Edge();
    Standard_Real first;
    Standard_Real last;
    Handle(Geom_Curve) myCurve = BRep_Tool::Curve(*edge, first, last);
    Handle(Geom_SurfaceOfLinearExtrusion) new_surface = 
               new Geom_SurfaceOfLinearExtrusion(myCurve, offset_dir); 
    Handle(Geom_RectangularTrimmedSurface) trimmed_surface = 
               new  Geom_RectangularTrimmedSurface(new_surface, first, last,
                                                   0, length);    
    if(trimmed_surface == NULL)
    { 
      PRINT_ERROR( "Can not tweak move the %dth curve\n", i);
      return CUBIT_FAILURE;
    }
#if OCC_VERSION_MINOR > 5
    TopoDS_Face FACE = BRepBuilderAPI_MakeFace(trimmed_surface, TOL);
#else
  #if OCC_VERSION_MAINTENANCE < 2
    TopoDS_Face FACE = BRepBuilderAPI_MakeFace(trimmed_surface);
  #else
    TopoDS_Face FACE = BRepBuilderAPI_MakeFace(trimmed_surface, TOL);
  #endif
#endif
    Surface*  extrude_surf= OCCQueryEngine::instance()->populate_topology_bridge(FACE, CUBIT_TRUE);
    BodySM* body = CAST_TO(extrude_surf, OCCSurface)->my_body();
    //subtract or unite the two surfaces
    DLIList<BodySM*> bodies;
    bodies.append(original_body);
    CubitStatus stat;
    if(!trim)
    {
      bodies.append(body);
      stat = unite(bodies, new_bodysm_list, keep_old_body);
    }
    else 
    {
      DLIList<BodySM*> tool_bodies;
      tool_bodies.append(body);
      stat = subtract(tool_bodies, bodies, new_bodysm_list, CUBIT_FALSE, 
                      keep_old_body);
    } 
    if(stat == CUBIT_FAILURE)
      return CUBIT_FAILURE;    
  }
  if(preview)
  {
    GfxPreview::clear();
    for(int i = 0; i < new_bodysm_list.size(); i++)
    {
      BodySM* new_bodysm = new_bodysm_list.get_and_step();
      TopoDS_Shape* modified_shape ;
      CAST_TO(new_bodysm, OCCBody)->get_TopoDS_Shape(modified_shape);
      TopExp_Explorer Ex;
      Ex.Init(*modified_shape, TopAbs_FACE);
      for( ; Ex.More(); Ex.Next() )
      {
        TopoDS_Face face = TopoDS::Face(Ex.Current());
        // Draw this face
        OCCDrawTool::instance()->draw_FACE( &face, CUBIT_BLUE_INDEX, CUBIT_TRUE );
      }
    }
    OCCQueryEngine::instance()->delete_solid_model_entities(new_bodysm_list);
    new_bodysm_list.clean_out();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Surface*> & /*surface_list*/, 
                                           double /*offset_distance*/,
                                           DLIList<Surface*> *add_surface_list_ptr,
                                           DLIList<double>*,
                                           DLIList<BodySM*> & /*new_bodysm_list*/,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Curve*> & /*curve_list*/,  
                                             double /*offset_distance*/,
                                             DLIList<Curve*>*,
                                             DLIList<double>*,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the 
//              remaining surfaces to fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Surface*> & /*surface_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*extend_adjoining*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove curves from a sheet body and then extend the 
//              remaining curves or fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Curve*> & /*curve_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a target 
//              surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Surface*> & /*surface_list*/,
                                           DLIList<Surface*> & ,
                                           DLIList<BodySM*> & /*new_bodysm_list*/,
                                           CubitBoolean extend_flg ,
                                           CubitPlane *limit_plane ,
                                           CubitBoolean /*reverse_flg*/,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Surface*> & /*target_surfs*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean extend_flg ,
                                           CubitPlane *limit_plane ,
                                           CubitBoolean ,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview,
                                           double max_area_increase ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a target
//              curve that is part of a sheet body.  The target is a surface 
//              created by thickening the owning surface of the target curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Curve*> & /*target_curve_ptr*/, 
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean extend_flg,
                                           CubitPlane *limit_plane ,
                                           CubitBoolean,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview,
                                           double max_area_increase ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::tweak_target( TBPoint *point_ptr,
                                    DLIList<Surface*> &modify_surface_list,
                                    CubitVector &target_loc,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body ,
                                    CubitBoolean preview  ) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: split curve at split_location, upper level geometry gets 
//              updated as well.
// Author     : Jane Hu
// Date       : 02/11
//================================================================================
CubitStatus  OCCModifyEngine::split_curve( Curve* curve_to_split,
                                    const CubitVector& split_location,
                                    DLIList<Curve*>& created_curves ) 
{
  //find if the curve is stand-along or in a body.
  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  DLIList <OCCBody* > *bodies = oqe->BodyList;
  DLIList<OCCLoop*> *loops = oqe->WireList;
  OCCCurve* occ_curve = CAST_TO(curve_to_split, OCCCurve);
  TopoDS_Edge* edge = occ_curve->get_TopoDS_Edge();
 
  LocOpe_SplitShape splitor;
  TopoDS_Shape from_shape;
  CubitBoolean found = false;
  //bodies consists compounds or solids.
  for(int i = 0; i < bodies->size(); i ++)
  {
    OCCBody* body = bodies->get_and_step();
    TopoDS_Shape* topo_shape;
    body->get_TopoDS_Shape(topo_shape);
    from_shape = *topo_shape;
    TopTools_IndexedDataMapOfShapeListOfShape M;

    TopExp::MapShapesAndAncestors(from_shape, TopAbs_EDGE, TopAbs_SHAPE, M);
    
    if(M.Contains(*edge) )   
    {
      found = true;
      break;
    }
  }
  
  for(int i = 0 ; found == false && i < loops->size(); i++)
  {
    OCCLoop* loop = loops->get_and_step();
    TopoDS_Wire* topo_loop = loop->get_TopoDS_Wire();
    from_shape = *topo_loop;
    TopTools_IndexedDataMapOfShapeListOfShape M;
    TopExp::MapShapesAndAncestors(from_shape, TopAbs_EDGE, TopAbs_WIRE, M);  
    if(M.Contains(*edge))
    {
      found = true;
      break;
    }
  }
  
  if(!found)
    from_shape = *edge;
  
  TopoDS_Shape* p_shape = &from_shape;
  CubitStatus status = split_shape_by_location(p_shape, curve_to_split, 
                                   split_location, created_curves);

  DLIList<TopoDS_Shape*> shape_list;
  DLIList<BodySM*> new_body_list;
  shape_list.append(p_shape);
  shape_to_bodySM(shape_list, new_body_list);
  return status;
}

//================================================================================
// Function   : split_shape_by_location
// Description: split curve at split_location, upper level geometry (from_shape)
//              gets updated as well.
// Author     : Jane Hu
// Date       : 02/11
//================================================================================
CubitStatus OCCModifyEngine::split_shape_by_location(TopoDS_Shape *&from_shape,
                                           Curve* curve_to_split,
                                           const CubitVector& split_location,
                                           DLIList<Curve*>& created_curves)const
{
  LocOpe_SplitShape splitor(*from_shape);
  CubitBoolean is_edge = (from_shape->ShapeType() == TopAbs_EDGE);
  TopoDS_Edge edge = *CAST_TO(curve_to_split, OCCCurve)->get_TopoDS_Edge();
  gp_Pnt pt = gp_Pnt(split_location.x(), split_location.y(), split_location.z());
  TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(pt);
  double param = curve_to_split->u_from_position(split_location);
  //double check that the split location is in curve
  double u_max, u_min;
  curve_to_split->get_param_range(u_min, u_max);
  if( param >= u_max || param <= u_min)
    return CUBIT_SUCCESS;

  splitor.Add(vertex, param, edge);
  
  //update the curve_list
  TopTools_ListOfShape edge_shapes;
  edge_shapes.Assign(splitor.DescendantShapes(edge));
  while(edge_shapes.Extent())
  {
     TopoDS_Shape edge_shape = edge_shapes.First();
     TopoDS_Edge occ_edge = TopoDS::Edge(edge_shape);
     OCCCurve* test_curve;
     if(!OCCQueryEngine::instance()->OCCMap->IsBound(occ_edge))
     {
       if(is_edge)
         test_curve = CAST_TO(OCCQueryEngine::instance()->populate_topology_bridge(occ_edge, CUBIT_TRUE), OCCCurve); 
       else
         test_curve = CAST_TO(OCCQueryEngine::instance()->populate_topology_bridge(occ_edge), OCCCurve);
       DLIList<OCCPoint*> points;
       test_curve->get_points(points);
       for(int i = 0 ; i <  points.size(); i++)
         points.get_and_step()->remove_curve(test_curve);
     }
     else
       test_curve = CAST_TO(OCCQueryEngine::instance()->populate_topology_bridge(occ_edge), OCCCurve);

     if(test_curve)
     {
       created_curves.append_unique(test_curve);
       //remove the points' curvelist of curve_to_split
       DLIList<OCCPoint*> points;
       test_curve->get_points(points);
       for(int i = 0; i < points.size(); i++)
       {
         OCCPoint* occ_p = points.get_and_step();
         occ_p->remove_curve(CAST_TO(curve_to_split, OCCCurve));
         if(from_shape->ShapeType() > TopAbs_EDGE)
           occ_p->remove_curve(test_curve);
       } 
     }
     edge_shapes.RemoveFirst();
  }

  TopTools_ListOfShape shapes;
  shapes.Assign(splitor.DescendantShapes(*from_shape));
  if(from_shape->ShapeType() ==TopAbs_COMPOUND)
    OCCBody::update_OCC_entity(*from_shape, shapes.First(),
                   (BRepBuilderAPI_MakeShape*) NULL, &splitor);

  else if(shapes.First().ShapeType() == TopAbs_SOLID)
    OCCLump::update_OCC_entity(TopoDS::Solid(*from_shape),
                   shapes.First(), (BRepBuilderAPI_MakeShape*) NULL, &splitor);

  else if(shapes.First().ShapeType() == TopAbs_SHELL)
    OCCShell::update_OCC_entity(TopoDS::Shell(*from_shape),
                   shapes.First(), (BRepBuilderAPI_MakeShape*) NULL, &splitor);

  else if(shapes.First().ShapeType() == TopAbs_FACE)
    OCCSurface::update_OCC_entity(TopoDS::Face(*from_shape),
              shapes.First(), (BRepBuilderAPI_MakeShape*) NULL, NULL, &splitor);

  else if(shapes.First().ShapeType() == TopAbs_WIRE)
    OCCLoop::update_OCC_entity(TopoDS::Wire(*from_shape), &splitor);
 
  from_shape = new TopoDS_Shape(shapes.First());
  return CUBIT_SUCCESS;
}

CubitStatus OCCModifyEngine::remove_curve_slivers( BodySM* body,
                                                   double lengthlimit ) const
{
  DLIList<CubitBoolean> is_volume;
  DLIList<BodySM*> bodies;
  DLIList<TopoDS_Shape*> shapes;
  bodies.append(body);
  CubitStatus status = get_shape_list(bodies, shapes, is_volume, CUBIT_FALSE);
  if(!status)
  {
    PRINT_ERROR("Can't find underlying TopoDS_Shape for this body.\n");
    return CUBIT_FAILURE;
  }
  Handle(ShapeBuild_ReShape) context;
  TopoDS_Shape new_shape = ShapeFix::RemoveSmallEdges(*shapes.get(), 
                           lengthlimit, context);
  new_shape = context->Apply(new_shape, TopAbs_COMPOUND);
  if(context->Status(ShapeExtend_OK))
  {
    PRINT_INFO("There's no small edges on this body.\n");
    return CUBIT_SUCCESS;
  }
     
  else if(context->Status(ShapeExtend_FAIL))
  {
    PRINT_ERROR("Small edges can't be removed from this body.\n");
    return CUBIT_FAILURE;
  }

  OCCQueryEngine::instance()->delete_solid_model_entities(body);
  DLIList<TopologyBridge*>tbs = OCCQueryEngine::instance()->populate_topology_bridge(new_shape);
  
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Surface*>& /*ref_face_list*/, 
                                                   BodySM *& /*new_body*/,
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_u*/, 
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_v*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Curve*>& /*u_curves*/, 
                                                   DLIList<Curve*>& /*v_curves*/,
                                                   BodySM *& /*new_body*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset surface.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_offset_surface( Surface* face_ptr, 
                                                    BodySM*& new_body, 
                                                    double offset ) const
{
  //create offset surface from its center along center normal of distance 
  //"offset"
  OCCSurface *occ_surface =  CAST_TO(face_ptr, OCCSurface);
  if (occ_surface == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    } 

  Surface* c_surface = NULL;
  c_surface = make_Surface(occ_surface);
  if (c_surface == NULL)
  {
    PRINT_ERROR("Cannot copy surface in sweep_translational.\n");
    return CUBIT_FAILURE;
  }
  occ_surface = CAST_TO(c_surface, OCCSurface);

  CubitVector center = occ_surface->center_point();
  CubitVector normal;
  occ_surface->closest_point(center,NULL,&normal); 
  CubitVector v = normal * offset;
  OCCQueryEngine::instance()->translate(occ_surface, v);
  new_body = occ_surface->my_body();
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Creates an offset body.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_offset_body( BodySM* body_ptr, 
                                                 BodySM*& new_bodysm, 
                                                 double offset_distance ) const
{
  PRINT_ERROR("Function not implemented because offset_distance \n"
               "doesn't show offset direction.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a skin surface.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
CubitStatus OCCModifyEngine::create_skin_surface( DLIList<Curve*>& curves, 
                                                  BodySM*& new_body,
                                                  DLIList<Curve*>&) const
{
   new_body = NULL;
   Surface* surf = make_Surface(BEST_FIT_SURFACE_TYPE, curves);
   if(surf)
   {
     new_body = CAST_TO(surf, OCCSurface)->my_body();
     return CUBIT_SUCCESS;
   } 
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a shell body from lofting surfaces.
// Author     : Jane Hu
// Date       : 01/09
//================================================================================
//CubitStatus OCCModifyEngine::loft_surfaces( Surface * face1, 
//                                            const double & /*takeoff1*/,
//                                            Surface * face2, 
//                                            const double & /*takeoff2*/,
//                                           BodySM*& new_body,
//                                            CubitBoolean /*arc_length_option*/, 
//                                            CubitBoolean /*twist_option*/,
//                                            CubitBoolean /*align_direction*/, 
//                                            CubitBoolean /*perpendicular*/,
//                                            CubitBoolean /*simplify_option*/ ) const
//{
//   BRepOffsetAPI_ThruSections loft(CUBIT_FALSE);
//   CubitStatus stat = do_loft(loft, face1, face2);
//   if(!stat)
//     return CUBIT_FAILURE;
//
//   TopoDS_Shape shape = loft.Shape();
//   TopoDS_Shell shell = TopoDS::Shell(shape);
//   OCCShell* occ_shell = OCCQueryEngine::instance()->populate_topology_bridge(shell, CUBIT_TRUE);
//   if (occ_shell == NULL)
//   {
//     PRINT_ERROR("In OCCModifyEngine::loft_surfaces\n"
//                 "   Cannot create a loft surface for given surfaces.\n");
//     return CUBIT_FAILURE;
//   }
//   new_body = occ_shell->my_body();
//   return CUBIT_SUCCESS;
//}

//================================================================================
// Description: Creates a solid body by lofting surfaces between surfaces
// Author     : Jane Hu
// Date       : 02/12
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces_to_body( DLIList<Surface*> &surfaces,
                             DLIList<double> &takeoff_factor_list, //not used
                             DLIList<Surface*> &takeoff_vector_surface_list, //not used
                             DLIList<CubitVector> &surface_takeoff_vector_list, //not used
                             DLIList<Curve*> &takeoff_vector_curve_list, //not used
                             DLIList<CubitVector> &curve_takeoff_vector_list, //not used
                             DLIList<Curve*> &guides, 
                             DLIList<TBPoint*> &match_vertices_list, //not used
                             BodySM*& new_body,
                             CubitBoolean global_guides,
                             CubitBoolean closed, //not used
                             CubitBoolean show_matching_curves,
                             CubitBoolean preview)const
{
   BRepOffsetAPI_ThruSections loft(CUBIT_TRUE);
   if(global_guides || guides.size() > 0)
     loft.Init(CUBIT_TRUE, CUBIT_TRUE);
   else
     loft.SetSmoothing(CUBIT_TRUE);
   CubitStatus stat = do_loft(loft, surfaces) ;
   if(!stat)
     return CUBIT_FAILURE;

   TopoDS_Shape shape = loft.Shape();
   if(preview && show_matching_curves)
   {
     PRINT_WARNING("Can't show matching curves in OCC.\n");
   }
   else if(preview)
   {
     GfxPreview::clear();
     // Draw this topoDS_shape
     OCCDrawTool::instance()->draw_TopoDS_Shape( &shape, CUBIT_BLUE_INDEX, CUBIT_FALSE, CUBIT_TRUE );

     return CUBIT_SUCCESS;
   } 
   TopoDS_Solid solid = TopoDS::Solid(shape);
   Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(solid, CUBIT_TRUE);
   if (lump == NULL)
   {
     PRINT_ERROR("In OCCModifyEngine::loft_surfaces_to_body\n"
                 "   Cannot create a loft body for given surfaces.\n");
     return CUBIT_FAILURE;
   }
   new_body = CAST_TO(lump, OCCLump)->get_body();
   return CUBIT_SUCCESS;   
}
 
CubitStatus OCCModifyEngine::do_loft(BRepOffsetAPI_ThruSections& loft,
                                  DLIList<Surface*> &surfaces) const
{
   for(int i = 0; i < surfaces.size(); i++)
   {
     OCCSurface* surf = CAST_TO(surfaces.get_and_step(), OCCSurface);
     if(!surf)
     {
       PRINT_ERROR("Surfaces are not OCC type.\n");
       return CUBIT_FAILURE;
     }
     TopoDS_Face* topo_face = surf->get_TopoDS_Face();
     TopExp_Explorer Ex;
     Ex.Init(*topo_face, TopAbs_WIRE);
     TopoDS_Wire wire = TopoDS::Wire(Ex.Current());
     if(Ex.More())
     {
       Ex.Next();
       if(Ex.Current().ShapeType() == TopAbs_WIRE)
       {
         PRINT_ERROR("Surface must have only one loop.\n");
         return CUBIT_FAILURE;
       }
     }
     loft.AddWire(wire);
   }
   loft.Build();
   if(!loft.IsDone())
   {
     PRINT_ERROR("Surfaces can't be loft into a body.\n");
     return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
}  

//================================================================================
// Description: Creates a surface using a list of vectors, project to surface if
//              given. Those points are connected into a closed wire and then
//              a surface is created on it.
// Author     : Jane Hu
// Date       : 02/11
//================================================================================
CubitStatus OCCModifyEngine::create_surface( DLIList<CubitVector*>& vec_list, 
                                             BodySM *& new_body, 
                                             Surface * surface_ptr,
                                             CubitBoolean project_points) const
{
  int i;
  CubitStatus stat;
  DLIList<TBPoint*> new_points; 
  DLIList<CubitVector*> new_vec_list;
  if (surface_ptr)
  {
    // Check the project_points option and do the necessary checks or projections.
    if (project_points)
    {
      // Create a new list of points that are projected to the surface
      vec_list.reset();
      CubitVector *vec_ptr, new_vec;
      for( i=0; i<vec_list.size(); i++ )
      {
        vec_ptr = vec_list.get_and_step();
        stat = surface_ptr->closest_point( *vec_ptr, &new_vec );
        if(stat)
          new_vec_list.append(new CubitVector(new_vec) );
        else
        {
          PRINT_ERROR("Can't project the %dth point onto surface \n", i+1); 
          for (int j = 0; j < new_vec_list.size(); j++)
            delete new_vec_list.get_and_step();
          return CUBIT_FAILURE;
        }         
      }
    } 

    else
    {
      // Make sure the points lie on the surface
      vec_list.reset();
      CubitVector loc_on_surf;
      CubitVector *vec_ptr;
      for( i=0; i<vec_list.size(); i++ )
      {
        vec_ptr = vec_list.get_and_step();
        surface_ptr->closest_point( *vec_ptr, &loc_on_surf );

        if (!vec_ptr->within_tolerance(loc_on_surf, GEOMETRY_RESABS))
        {
          PRINT_ERROR("all locations must lie on Surface\n" );
          return CUBIT_FAILURE;
        }
      }
      new_vec_list = vec_list;
    }
  }

  // Make the surface in the solid modeller
  else
    new_vec_list = vec_list;

  for( i=0; i<new_vec_list.size(); i++ )
  {
     CubitVector* vec_ptr = new_vec_list.get_and_step();
     OCCPoint* point = new OCCPoint(*vec_ptr);
     new_points.append( (TBPoint*)point );
  }

  stat =  create_surface( new_points, new_body, NULL);

  if(project_points && surface_ptr)
  {
    for( i=0; i<new_vec_list.size(); i++ )
    {
      CubitVector *vec_ptr = new_vec_list.get_and_step();
      delete vec_ptr;
    }
  }
  return stat;
}

//================================================================================
// Description: Creates a surface using a list of points.
//              Those points are connected into a closed wire and then
//              a surface is created on it.
// Author     : Jane Hu
// Date       : 02/11
//================================================================================
CubitStatus
OCCModifyEngine::create_surface( DLIList<TBPoint*>& points,
                               BodySM *&new_body,
                               Surface * /*on_surface*/ )const
{
  TBPoint *start_point = points.get_and_step();
  DLIList<Curve*> curve_list;
  for(int i = 0; i < points.size(); i++)
  {
    CubitVector coord1 = start_point->coordinates();
    TBPoint *end_point = points.get_and_step();
    CubitVector coord2 = end_point->coordinates();
    if(coord1.within_tolerance( coord2, GEOMETRY_RESABS ) )
    {
       PRINT_ERROR( "Attempt to create a line between coincident points at (%f, %f, %f)\n",
          coord1.x(), coord1.y(), coord1.z() );
       for(int  j=0; j<curve_list.size(); j++ )
       {
         Curve* curve_ptr = curve_list.get_and_step();
         OCCQueryEngine::instance()->delete_solid_model_entities( curve_ptr );
       }
       return CUBIT_FAILURE;
    }

    Curve* curve = make_Curve(start_point, end_point);
    curve_list.append(curve);
    start_point = end_point;
  }

  //make surface out of curves.
  Surface* surf = make_Surface(BEST_FIT_SURFACE_TYPE, curve_list); 
  OCCSurface* occ_surf;
  if(surf)
    occ_surf = CAST_TO(surf, OCCSurface);
  else
  {
    PRINT_ERROR("Failed to create a surface from given points. \n");
    return CUBIT_FAILURE;
  }
  new_body = occ_surf->my_body(); 
  return CUBIT_SUCCESS; 
}


//================================================================================
// Description: Creates a weld surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_weld_surface( CubitVector & /*root*/,
                                                    Surface * /*ref_face1*/, 
                                                    double /*leg1*/, 
                                                    Surface * /*ref_face2*/, 
                                                    double /*leg2*/,
                                                    BodySM *& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::scale( BodySM *&body, const CubitVector& factors )
{
  return OCCQueryEngine::instance()->scale( body, factors );
}

//================================================================================
// Description: According to what I read in AcisModifyEngine, 
//              bridges_in are already the invalid_tbs, and we want to 
//              get all vertices, edges and surfaces below bridges_in into
//              bridges_out.
// Author     : Jane Hu 
// Date       : 02/11
//================================================================================
void OCCModifyEngine::get_possible_invalid_tbs(
                             DLIList<TopologyBridge*> &bridges_in,
                             DLIList<TopologyBridge*> &bridges_out)
{
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCCurve*> curves;
  DLIList<OCCPoint*> points;

  for (int i = 0; i < bridges_in.size(); i++)
  {
     TopologyBridge* tb_in = bridges_in.get_and_step();
     if(! tb_in)
       continue;

     OCCBody *occ_body = NULL;
     if(OCCLump * lump_ptr = CAST_TO( tb_in,OCCLump))
     {
       BodySM* body = lump_ptr->get_body();
       assert (body);
       occ_body = CAST_TO(body, OCCBody);
     }
 
     else if( occ_body || (occ_body = CAST_TO( tb_in, OCCBody)))
     {
       occ_body->get_all_surfaces(surfaces);
       occ_body->get_all_curves(curves);
       occ_body->get_all_points(points);
     }

     else if(OCCSurface *surface_ptr = CAST_TO( tb_in, OCCSurface))
     {
       surface_ptr->get_curves(curves);
       surface_ptr->get_points(points);
       surfaces.append(surface_ptr);
     }
     else if( OCCCurve *curve_ptr = CAST_TO( tb_in, OCCCurve))
     {
       curves.append(curve_ptr);
       curve_ptr->get_points(points);
     }
     else if( OCCPoint *point_ptr = CAST_TO( tb_in, OCCPoint))
     {
       points.append(point_ptr);
     }
  }
  for (int i = 0; i < surfaces.size(); i++)
    bridges_out.append((TopologyBridge*)surfaces.get_and_step());
  for (int i = 0; i < curves.size(); i++)
    bridges_out.append((TopologyBridge*)curves.get_and_step());
  for (int i = 0; i < points.size(); i++)
    bridges_out.append((TopologyBridge*)points.get_and_step());
}
CubitStatus OCCModifyEngine::curve_surface_intersection( Surface *surface,
                                            Curve* curve,
                                            DLIList<Curve*> &new_curves ) const
{
  return CUBIT_FAILURE;
}

BodySM* OCCModifyEngine::make_extended_sheet( DLIList<Surface*> &surface_list,
                                              CubitBox *clip_box_ptr,
                                              bool preview ) const
{
  PRINT_ERROR("This feature is not implemented.\n");
  return (BodySM*) NULL;
}

CubitStatus OCCModifyEngine::remove_topology(DLIList<Curve*> &curves_to_remove,
                                       DLIList<Surface*> &surfs_to_remove,
                                       double backoff_distance,
                                       double small_edge_size,
                                       DLIList<BodySM*> &new_bodysm_list,
                                       CubitBoolean preview) const
{
  PRINT_ERROR("This feature is not implemented.\n");
  return CUBIT_FAILURE;;
}

CubitStatus OCCModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                               DLIList<BodySM*> &new_bodies,
                                               DLIList<TopologyBridge*>*,
                                               DLIList<TopologyBridge*>* )  const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::tolerant_imprint(DLIList<Surface*> &surfs_in,
                             DLIList<BodySM*> &new_bodysm_list) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::tolerant_imprint_surface_with_curves(
                                  Surface *surface_to_imprint,
                                  DLIList<Curve*> &curves,
                                  DLIList<TopologyBridge*> &temporary_bridges,
                                  BodySM *&new_body,
                                  DLIList<TopologyBridge*> *new_tbs,
                                  DLIList<TopologyBridge*> *att_tbs) const
{
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_new_tbs
// Member Type: PUBLIC
// Description: given old list of surfaces, curves and points, and new list
//              of surfaces, curves and points, using physical comparison to
//              find out what entities are newly generated by imprinting.
//              And add imprint feature on the new entities.
// Author     : Jane HU
// Date       : 02/11
//===============================================================================
void OCCModifyEngine::get_new_tbs( 
      std::map<OCCSurface*, std::pair<CubitVector, double> >& surf_property_map,
      std::map<OCCCurve*, std::pair<CubitVector, double> >& curve_property_map,
      DLIList<OCCPoint*> &points, 
      DLIList<OCCSurface*> &new_surfaces,
      DLIList<OCCCurve*> &new_curves, 
      DLIList<OCCPoint*> &new_points, 
      DLIList<TopologyBridge*> *new_tbs)const
{
/*
  CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );

  CubitBoolean found = false;
  for (int i = 0; i < new_surfaces.size() ; i++)
  {
    OCCSurface* new_surf = new_surfaces.get_and_step(); 
    double new_d = new_surf->measure();
    CubitVector new_center = new_surf->center_point();
    std::map<OCCSurface*, std::pair<CubitVector, int> >::iterator it = 
            surf_property_map.find( new_surf);
    if(it != surf_property_map.end())
    {
      std::pair<CubitVector, int> pair;
      pair = it->second;
      double d = pair.second;
      CubitVector center = pair.first;
      surf_property_map.erase(it);
      if(center.about_equal(new_center) && fabs(d- new_d) <= TOL)
        found = true;
    }
    if(!found)
    {
      new_tbs->append(new_surf); 
      //Add an imprint feature to the topo.
      TopoDS_Face* new_shape = new_surf->get_TopoDS_Face(); 
      OCCAttribSet::append_attribute(tmp_attrib, *new_shape); 
    }
    found = false;
  }

  for (int i = 0; i < new_curves.size() ; i++)
  {
    OCCCurve* new_curve = new_curves.get_and_step();
    double new_d = new_curve->measure();
    CubitVector new_center;
    new_curve->position_from_u(0.5, new_center);
    std::map<OCCCurve*, std::pair<CubitVector, int> >::iterator it = 
      curve_property_map.find(new_curve);
    if( it != curve_property_map.end())
    {
      std::pair<CubitVector, int> pair;
      pair = it->second;
      double d = pair.second;
      CubitVector center = pair.first;
      curve_property_map.erase(it);
      if(center.about_equal(new_center) && fabs(d- new_d) <= TOL)
        found = true;
    }
    if(!found)
    {
      new_tbs->append(new_curve);
      //Add an imprint feature to the topo.
      TopoDS_Edge* new_shape = new_curve->get_TopoDS_Edge();
      OCCAttribSet::append_attribute(tmp_attrib, *new_shape);
    }
    found = false;
  }

  int num_new_points = new_points.size() -points.size(); 
  if(num_new_points > 0)
  {
    int  new_vertices_count = 0;
    CubitBoolean found = false; 
    for( int i = 0; i < new_points.size() && new_vertices_count < num_new_points; i++)
    {
      OCCPoint* new_point = new_points.get_and_step();
      for(int j = 0; j < points.size(); j ++)
      {
        OCCPoint* point = points.get(); 
        if(point == new_point)
        {
          found = true;
          points.remove();
          break;
        }
        points.step();
      }
      if(!found)
      {
        new_tbs->append(new_point);
        new_vertices_count ++;
        //Add an imprint feature to the topo.
        TopoDS_Vertex* new_shape = new_point->get_TopoDS_Vertex();
        OCCAttribSet::append_attribute(tmp_attrib, *new_shape);
      }
      found = false;
    }
  }
  delete tmp_attrib;
*/
}

//===============================================================================
// Function   : get_att_tbs
// Member Type: PUBLIC
// Description: given new list of surfaces, curves and points for modified
//              bodies, find out what entities are "name" entities
// Author     : Jane HU
// Date       : 02/11
//===============================================================================
void OCCModifyEngine::get_att_tbs(DLIList<OCCSurface*> &new_surfaces,
                                  DLIList<OCCCurve*> &new_curves,
                                  DLIList<OCCPoint*> &new_points,
                                  const CubitString& name,
                                  DLIList<TopologyBridge*> *att_tbs)const
{
  DLIList<CubitSimpleAttrib> csa_list;
  for (int i = 0 ; i < new_surfaces.size(); i++)
  {
    OCCSurface* occ_surf = new_surfaces.get_and_step();
    csa_list.clean_out();
    occ_surf->get_simple_attribute(name, csa_list);
    if(csa_list.size() > 0)
      att_tbs->append(occ_surf);    
  }

  for (int i = 0 ; i < new_curves.size(); i++)
  {
    OCCCurve* occ_curve = new_curves.get_and_step();
    csa_list.clean_out();
    occ_curve->get_simple_attribute(name, csa_list);
    if(csa_list.size() > 0)
      att_tbs->append(occ_curve);
  }
 
  for (int i = 0 ; i < new_points.size(); i++)
  {
    OCCPoint* occ_point = new_points.get_and_step();
    csa_list.clean_out(); 
    occ_point->get_simple_attribute(name, csa_list);
    if(csa_list.size() > 0)
      att_tbs->append(occ_point);
  }
}

Curve* OCCModifyEngine::create_arc_radius(const CubitVector &center,
                                          TBPoint* ref_vertex_start,
                                          TBPoint* ref_vertex_end,
                                          const CubitVector &normal,
                                          double radius,
                                          bool full ,
                                          CubitBoolean preview  )
{
  CubitVector v1 = ref_vertex_start->coordinates();
  CubitVector v2 = ref_vertex_end->coordinates();  
  if ((2 * radius) < v1.distance_between(v2))
  {
    PRINT_ERROR("Unable to create acr from given radius and vertices. \n");
    return (Curve*) NULL;
  }
  OCCPoint* center_point = new OCCPoint(center);
  Curve* curve =  create_arc_center_edge(center_point, ref_vertex_start, 
                                ref_vertex_end,
                                normal, radius, full , preview  );
  delete center_point;
  
  return curve;
}

Curve* OCCModifyEngine::create_arc(const CubitVector& position,
                                   double radius,
                                   double start_angle,//in degree
                                   double end_angle,//in degree
                                   CubitVector plane,
                                   bool preview )
{
  gp_Pnt P(position.x(), position.y(), position.z());
  gp_Dir normal(plane.x(), plane.y(), plane.z());  
  gp_Ax1 axis(P, normal);
  Handle(Geom_Circle) curve_ptr;
  curve_ptr = GC_MakeCircle( axis, radius);
  gp_Circ circle = curve_ptr->Circ();
  start_angle *=  CUBIT_PI/180;
  end_angle *=  CUBIT_PI/180; 
  CubitBoolean sense = CUBIT_TRUE ;

  Handle_Geom_TrimmedCurve new_curve_ptr;
  new_curve_ptr =  GC_MakeArcOfCircle(circle, start_angle, end_angle, sense);
  TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(new_curve_ptr);
  Curve*  new_curve = OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
  if(preview)
  {
    GfxPreview::clear();
    OCCCurve* occ_curve = CAST_TO(new_curve, OCCCurve);
    TopoDS_Edge* h_edge = occ_curve->get_TopoDS_Edge();
    // Draw this edge
    OCCDrawTool::instance()->draw_EDGE( h_edge, CUBIT_BLUE_INDEX, CUBIT_TRUE );

    OCCQueryEngine::instance()->delete_solid_model_entities(new_curve);
    return (Curve*) NULL;
  }
  else
    return new_curve;
}

CubitStatus OCCModifyEngine::webcut_across_translate( 
                                     DLIList<BodySM*>& body_list,
                                     Surface* plane_surf1,
                                     Surface* plane_surf2,
                                     DLIList<BodySM*>& neighbor_imprint_list,
                                     DLIList<BodySM*>& results_list,
                                     ImprintType imprint_type ,
                                     bool preview ) const
{
  PRINT_WARNING("Not implemented yet. \n");
  return CUBIT_FAILURE;
}
//void OCCModifyEngine::start_tracking_history( DLIList<TopologyBridge*> &bridges,
//                                              OCCHistory &history_object,
//                                              bool ignore_parents )
//{
//}

//void OCCModifyEngine::stop_tracking_history( DLIList<BodySM*> &new_bodies,
//                                             OCCHistory &history_object )
//{
//  DLIList<TopologyBridge*> tbs;
//  CAST_LIST( new_bodies, tbs, TopologyBridge );

  //stop_tracking_history( tbs, history_object );
//}
 
//Following codes are OCC engine only, implemented by Boyd Tidwell
CubitStatus OCCModifyEngine::create_rectangle_surface( double width,
            double height,
            CubitVector plane,
            BodySM *&sheet_body) const
{
  //create points at the 4 corners
  CubitVector pos1( width*0.5, height*0.5, 0 );
  CubitVector pos2( -width*0.5, height*0.5, 0 );
  CubitVector pos3( -width*0.5, -height*0.5, 0 );
  CubitVector pos4( width*0.5, -height*0.5, 0 );

  DLIList<CubitVector*> vec_list(4);
  vec_list.append( &pos1 );
  vec_list.append( &pos2 );
  vec_list.append( &pos3 );
  vec_list.append( &pos4 );

  CubitStatus status = create_surface(vec_list, sheet_body, 0, CUBIT_FALSE );

  //rotate it so that it is aligned with the plane defined by vector 'plane'

  CubitVector current_plane(0,0,1);
  if (!plane.about_equal(current_plane))
  {
    CubitVector axis_of_rotation = current_plane * plane;
    double angle_of_rotation = axis_of_rotation.vector_angle(
         current_plane, plane );

    OCCQueryEngine::instance()->rotate( sheet_body, axis_of_rotation,
                                        RADIANS_TO_DEGREES(angle_of_rotation) );
  }

  return status;
}

CubitStatus OCCModifyEngine::create_ellipse_surface( TBPoint *pt1,
                                                     TBPoint *pt3,
                                                     CubitVector center_point,
                                                     BodySM *&sheet_body) const
{
  // make ellipse by creating two ellipse arcs in each sense direction
  Curve *ellipse1 = make_elliptical_Curve( pt1, pt3, center_point, 0, 360, CUBIT_FORWARD );
  if( NULL == ellipse1 )
    return CUBIT_FAILURE;
  DLIList<Curve*> curve_list(2);
  curve_list.append( ellipse1 );

  Curve *ellipse2 = make_elliptical_Curve( pt1, pt3, center_point, 0, 360, CUBIT_REVERSED );
  if( NULL == ellipse2 )
    return CUBIT_FAILURE;
  curve_list.append( ellipse2 );

  Surface *tmp_surface = make_Surface( PLANE_SURFACE_TYPE, curve_list, NULL, false );
  if( NULL == tmp_surface )
    return CUBIT_FAILURE;
  sheet_body = make_BodySM( tmp_surface );
  if( NULL == sheet_body)
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitStatus OCCModifyEngine::create_ellipse_surface( double major_radius,
                                                     double minor_radius,
                                                     CubitVector plane,
                                                     BodySM *&sheet_body) const
{
  //create the points
  CubitVector tmp_pt( major_radius, 0, 0 );
  TBPoint *pt1 = make_Point( tmp_pt );

  tmp_pt.set( 0, minor_radius, 0 );
  TBPoint *pt2 = make_Point( tmp_pt );

  CubitVector center_point(0,0,0);

  // make ellipse by creating two ellipse arcs in each sense direction
  Curve *ellipse1 = make_elliptical_Curve( pt1, pt2, center_point, 0, 360, CUBIT_FORWARD );
  if( NULL == ellipse1 )
    return CUBIT_FAILURE;
  DLIList<Curve*> curve_list(2);
  curve_list.append( ellipse1 );

  Curve *ellipse2 = make_elliptical_Curve( pt1, pt2, center_point, 0, 360, CUBIT_REVERSED );
  if( NULL == ellipse2 )
    return CUBIT_FAILURE;
  curve_list.append( ellipse2 );

  Surface *tmp_surface = make_Surface(PLANE_SURFACE_TYPE, curve_list, NULL, false );
  if( NULL == tmp_surface )
    return CUBIT_FAILURE;

  sheet_body = make_BodySM( tmp_surface );
  if( NULL == sheet_body)
    return CUBIT_FAILURE;

  //rotate it so that it is aligned with the plane defined by vector 'plane'
  CubitVector current_plane(0,0,1);
  if (!plane.about_equal(current_plane))
  {
    CubitVector axis_of_rotation = current_plane * plane;
    double angle_of_rotation = axis_of_rotation.vector_angle(current_plane, plane );

    OCCQueryEngine::instance()->rotate( sheet_body, axis_of_rotation,
                                        RADIANS_TO_DEGREES(angle_of_rotation) );
  }

  return CUBIT_SUCCESS;
}

Curve* OCCModifyEngine::make_elliptical_Curve( TBPoint const* point1,
                                               TBPoint const* point2,
                                               CubitVector &center_point,
                                               double start_angle,
                                               double end_angle,
                                               CubitSense sense) const
{
  GeometryType curve_type = ELLIPSE_CURVE_TYPE;
  Curve *curve = NULL;
  if(sense ==CUBIT_FORWARD)
    curve = make_Curve(curve_type, point1, point2, &center_point );
  else
    curve = make_Curve(curve_type, point2, point1, &center_point );
  if ( curve == NULL )
  {
    PRINT_ERROR("In OCCModifyEngine::make_elliptical_Curve\n"
                "       Cannot make Curve object.\n");
    return (Curve *)NULL;
  }
  return curve;
}

//circle surface passes three points.
CubitStatus OCCModifyEngine::create_circle_surface( TBPoint *pt1,
                                                    CubitVector center_point,
                                                    TBPoint *pt3,
                                                    BodySM *&sheet_body)const
{
  CubitVector * center = new CubitVector(center_point);
  Curve *circle = make_Curve( ARC_CURVE_TYPE, pt1, pt3, center);
  delete center;

  DLIList<Curve*> curve_list(1);
  curve_list.append( circle );
  Surface *tmp_surface = make_Surface(PLANE_SURFACE_TYPE, curve_list, NULL, false );

  if( NULL == tmp_surface )
    return CUBIT_FAILURE;

  sheet_body = make_BodySM( tmp_surface );
  if( NULL == sheet_body)
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//circle surface passes pt1, pt3, and centered at center_point.
CubitStatus OCCModifyEngine::create_circle_surface( TBPoint *pt1,
                                                    TBPoint *pt3,
                                                    CubitVector center_point,
                                                    BodySM *&sheet_body)const
{
  CubitVector v1 = pt1->coordinates();
  CubitVector v3 = pt3->coordinates();
  CubitVector normal(0,0,0);
  CubitVector dir1 = v1 - center_point;
  CubitVector dir2 = v3 - center_point;

  normal = dir1 * dir2;
  if(normal.length() < TOL)
  {
    PRINT_ERROR("The given points can't decide the surface plane. \n");
    return CUBIT_FAILURE;
  }
  Handle(Geom_Circle) curve_ptr;
  gp_Dir norm(normal.x(), normal.y(), normal.z());
  gp_Pnt center = gp_Pnt( center_point.x(), center_point.y(), center_point.z());

  double radius = v1.distance_between(center_point);
  curve_ptr = GC_MakeCircle(center,norm,radius);

  OCCPoint* occ_pt1 = CAST_TO(const_cast<TBPoint*>(pt1),OCCPoint);
  TopoDS_Vertex * vt1 = occ_pt1->get_TopoDS_Vertex();
  TopoDS_Edge new_edge;
  new_edge = BRepBuilderAPI_MakeEdge(curve_ptr, *vt1, *vt1);

  Curve* circle = OCCQueryEngine::instance()->
                   populate_topology_bridge(new_edge, CUBIT_TRUE);
  DLIList<Curve*> curve_list(1);
  curve_list.append( circle );
  Surface *tmp_surface = make_Surface(PLANE_SURFACE_TYPE, curve_list, NULL,
                                      false );
  if( NULL == tmp_surface )
    return CUBIT_FAILURE;

  sheet_body = OCCModifyEngine::instance()->make_BodySM( tmp_surface );
  if( NULL == sheet_body)
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitStatus OCCModifyEngine::create_circle_surface( double radius,
                                                    CubitVector plane,
                                                    BodySM *&sheet_body) const
{
  CubitVector center_pt;
  center_pt.set(0, 0, 0);
  CubitVector pt2, pt3;
  if (plane.x() > 0 )
  {
    pt2.set( 0, radius, 0 );
    pt3.set( 0, 0, radius);
  }
  else if (plane.y() > 0)
  {
    pt2.set(radius, 0, 0);
    pt3.set(0, 0, -radius);
  }
  else if (plane.z() > 0)
  {
    pt2.set(radius, 0, 0);
    pt3.set(0, radius, 0);
  }
  else
  {
    PRINT_ERROR("In OCCModifyEngine::create_circle_surface\n"
                "       Invalid plane specified.\n");
    return CUBIT_FAILURE;
  }

  TBPoint *tbpt2 = make_Point( pt2 );
  TBPoint *tbpt3 = make_Point( pt3 );

  CubitStatus stat = create_circle_surface(tbpt2, tbpt3, center_pt, sheet_body); 
  return stat;

}

CubitStatus OCCModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                        DLIList<BodySM*> &new_bodies,
                                        double overlap_tol,
                                        double imprint_tol,
										DLIList<TopologyBridge*> *new_tbs,
										DLIList<TopologyBridge*> *att_tbs) const
{
	return CUBIT_FAILURE;
}

#ifdef CGM_KCM  
CubitStatus OCCModifyEngine::mesh2brep(std::vector<double> &xvals,
                        std::vector<double> &yvals,
                        std::vector<double> &zvals,
                        std::vector<unsigned int> &tri_connectivity,
                        DLIList<BodySM*> &new_body_sms) const
{
	return CUBIT_FAILURE;
}
#endif

// EOF
