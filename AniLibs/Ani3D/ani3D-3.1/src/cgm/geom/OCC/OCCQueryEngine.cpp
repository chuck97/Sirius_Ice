//-------------------------------------------------------------------------
// Filename      : OCCQueryEngine.cpp
//
// Purpose       : Implementation of the OCCQueryEngine class.
//                 This class provides OCC-based implementations
//                 of various virtual functions in the GeometryQueryEngine
//                 hierarchy.
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 7/17/00
//
//-------------------------------------------------------------------------
#include <Standard_Version.hxx>
#include <Standard_Stream.hxx>
//#include <Standard_SStream.hxx>
//#include <Standard_String.hxx>
//#include <stringbuf>
#if OCC_VERSION_MINOR < 5
  #include "TDataStd_Shape.hxx"
  typedef TDataStd_Shape TDataXtd_Shape;
  typedef Handle_TDataStd_Shape Handle_TDataXtd_Shape;
#else
  #include "TDataXtd_Shape.hxx"
#endif

#include "BRep_Tool.hxx"
#include "gp_Pnt.hxx"
#include "gp_Ax1.hxx"
#include "gp_Ax2.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "Interface_Static.hxx"
#include "BRepBuilderAPI.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepBuilderAPI_GTransform.hxx"
#include "BRepBuilderAPI_ModifyShape.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "OCCShapeAttributeSet.hpp"
//#include "OCCBinToolsShapeSet.hpp"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "Poly_Triangle.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "Handle_Poly_Triangulation.hxx"
#include "GCPnts_TangentialDeflection.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopExp.hxx"
#ifdef HAVE_OCC_STEP
#  include "STEPControl_Reader.hxx"
#  include "STEPControl_Writer.hxx"
#  include "STEPControl_StepModelType.hxx"
#endif
#ifdef HAVE_OCC_IGES
#  include "IGESControl_Reader.hxx"
#  include "IGESControl_Writer.hxx"
#endif
#include "IFSelect_ReturnStatus.hxx"
#include "BndLib_Add3dCurve.hxx"
#include "Poly_Polygon3D.hxx"
#include "Handle_Poly_Polygon3D.hxx"
#include "BRepMesh_FastDiscret.hxx"
#include "OCCQueryEngine.hpp"
#include "OCCModifyEngine.hpp"
#include "Poly_Triangulation.hxx"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "GeometryEntity.hpp"
#include "DLIList.hpp"
#include "CubitBox.hpp"
#include "CubitString.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCAttribSet.hpp"
#include "GMem.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitObserver.hpp"
#include "GfxDebug.hpp"
#include <stdio.h>
#include <errno.h>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include "TopoDS_Compound.hxx"
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <Precision.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <TDF_ChildIterator.hxx>
#include <BinTools_ShapeSet.hxx>
#include "Standard_Boolean.hxx"

#include "TDF_Label.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "TDocStd_Document.hxx"
#include "TCollection_ExtendedString.hxx"
#include "gp_Lin.hxx"
using namespace NCubitFile;

#include "CGMEngineDynamicLoader.hpp"


CGM_ENGINE_EXPORT_CREATE_GQE(OpenCascade)
{
  return OCCQueryEngine::instance();
}

OCCQueryEngine* OCCQueryEngine::instance_ = NULL;

typedef std::map<int, TopologyBridge*>::value_type valType;
typedef std::map<int, TDF_Label>::value_type labType;
int OCCQueryEngine::iTotalTBCreated = 0;
int OCCQueryEngine::total_coedges = 0;
#define NUM_PTS_UV 30
//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
OCCQueryEngine* OCCQueryEngine::instance()
{
  if (instance_ == NULL ) {
    instance_ = new OCCQueryEngine;
  }
  return instance_;
}

//================================================================================
// Description:  default constructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::OCCQueryEngine()
{
  GeometryQueryTool::instance()->add_gqe( this );
  OCCMap = new TopTools_DataMapOfShapeInteger;
  OccToCGM = new std::map<int, TopologyBridge*>;
  Shape_Label_Map = new std::map<int, TDF_Label>;
  BodyList = new DLIList<OCCBody*>;
  WireList = new DLIList<OCCLoop*>;
  SurfaceList = new DLIList<OCCSurface*>;
  CurveList = new DLIList<OCCCurve*>;
  CubitString name("Doc");
  TCollection_ExtendedString xString((Standard_CString)name.c_str(), CUBIT_TRUE);
  MyDF = new TDocStd_Document(xString);
  mainLabel = MyDF->Main();
  EXPORT_ATTRIB = CUBIT_TRUE;
}

//================================================================================
// Description:  destructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::~OCCQueryEngine()
{
  instance_ = NULL;
  delete OCCMap;
  delete OccToCGM;
  delete Shape_Label_Map;
  delete BodyList;
  delete WireList;
  delete SurfaceList;
  delete CurveList;
}

int OCCQueryEngine::get_major_version()
{
  return OCC_VERSION_MAJOR;
}

int OCCQueryEngine::get_minor_version()
{
  return OCC_VERSION_MINOR;
}

int OCCQueryEngine::get_subminor_version()
{
  return OCC_VERSION_MAINTENANCE;
}

CubitString OCCQueryEngine::get_engine_version_string()
{
  return CubitString("OpenCascade ") + OCC_VERSION_STRING;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
const type_info& OCCQueryEngine::entity_type_info() const
{
  return typeid(OCCQueryEngine);
}

//================================================================================
// Description: reflect body about an axis
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
CubitStatus OCCQueryEngine::reflect( BodySM *bodysm,
                                     const CubitVector& axis)
{
  OCCBody *body = CAST_TO(bodysm, OCCBody);
  if (!body)
    {
      PRINT_ERROR("Attempt to reflect OCC-based geometry Body.  This body is not OCCBody.");
      return CUBIT_FAILURE;
    }

  body->reflect( axis.x(), axis.y(), axis.z() );

  return CUBIT_SUCCESS;
}


//================================================================================
// Description:  This function queries OCC for the necessary facets
//               information needed in facetting a RefFace.  This
//               information is stored and output in gMem.  The
//               number of triangles, points and facets are also
//               output.
//               normal_tolerance is in degree.
//               max_edge_length is not considered in getting graphics.
// Author     :  Jane Hu
// Date       :  10/25/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Surface* surface_ptr,
                                          GMem* g_mem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length) const
{
  // Because this may be unnecessarily called twice,
  // say there is one triangle.
  if (!g_mem)
      return CUBIT_SUCCESS;

  if(max_edge_length > get_sme_resabs_tolerance())
  {
    PRINT_WARNING("OCC surface's tessilation doesn't consider edge_length.\n");
    PRINT_WARNING("max_edge_length argument is ignored. \n");
  }

  OCCSurface *occ_surface_ptr = CAST_TO(surface_ptr, OCCSurface);
  TopoDS_Face * Topo_Face = occ_surface_ptr->get_TopoDS_Face();

  return get_graphics( Topo_Face, g_mem, normal_tolerance, distance_tolerance, max_edge_length );
}

CubitStatus OCCQueryEngine::get_graphics( TopoDS_Face* face_ptr,
                                          GMem* g_mem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                         double max_edge_length) const
{
  if (!face_ptr)
    return CUBIT_FAILURE;

  TopLoc_Location L;
  Handle_Poly_Triangulation facets = BRep_Tool::Triangulation(*face_ptr, L);

  gp_Trsf tf = L.Transformation();

  if(facets.IsNull() || facets->NbTriangles() == 0)
  {
    //do triangulation
    if(distance_tolerance <= 0.0)
      distance_tolerance = 0.01;
    double angle  = CUBIT_PI * normal_tolerance/180;
    BRepAdaptor_Surface asurface(*face_ptr);
    Bnd_Box aBox;
    BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
    BRepMesh_FastDiscret *myMesh =
#if OCC_VERSION_MAJOR <= 6 && OCC_VERSION_MINOR < 8
      new BRepMesh_FastDiscret(distance_tolerance, *face_ptr, aBox, angle, Standard_True, Standard_True);
#else
      new BRepMesh_FastDiscret(*face_ptr, distance_tolerance, angle, aBox, Standard_True, Standard_True);
#endif
    if (myMesh != NULL) delete myMesh;
    facets = BRep_Tool::Triangulation(*face_ptr, L);
    if(facets.IsNull() || facets->NbTriangles() == 0)
    {
      PRINT_ERROR("Can't get triangulation representation for this surface.\n");
      return CUBIT_FAILURE;
    }
  }
  //if necessary, the face tolerance can be returned. now, no use.
  //double tol = BRep_Tool::Tolerance(*Topo_Face);   

  int  number_points = facets->NbNodes();
  int  number_triangles = facets->NbTriangles();
  int  number_facets = 4 * number_triangles; 
  
  Poly_Array1OfTriangle triangles(0, number_triangles-1);
  triangles.Assign( facets->Triangles() );
  int *facetList =  new int[number_facets];
  //needs to test that N1, N2, N3 index are starting from 0 to number_points-1
  //otherwise needs to update either facetList or gPnts to make consistent.
  //It's possible also that N's starting from 1.
  int minN = 100;
  for (int i = 0; i < triangles.Length(); i++)
    {
      Poly_Triangle triangle = triangles.Value( i );
      int N1, N2, N3;
      triangle.Get(N1, N2, N3); 
      facetList[4 * i] = 3;
      facetList[4 * i + 1] = N1;
      minN = (minN < N1 ? minN : N1);
      facetList[4 * i + 2] = N2;
      minN = (minN < N2 ? minN : N2);
      facetList[4 * i + 3] = N3;
      minN = (minN < N3 ? minN : N3);
    } 
  if(minN != 0)
  {
    //subtract the minN from the facetList for all i+1, i+2, i+3 points
    for (int i = 0; i < triangles.Length(); i++)
    {
      facetList[4 * i + 1] -= minN;
      facetList[4 * i + 2] -= minN;
      facetList[4 * i + 3] -= minN;
    }
  }
  g_mem->replace_facet_list( facetList, number_facets, number_facets); 

  TColgp_Array1OfPnt points(0,  number_points-1);
  points.Assign(facets->Nodes());
  GPoint *gPnts= new GPoint[number_points];
  for (int i = 0; i < number_points ; i ++)
    {
      gp_Pnt gp_pnt = points.Value(i);
      if( !L.IsIdentity())
        gp_pnt.Transform(tf);

      GPoint gPnt;
      gPnt.x = gp_pnt.X();
      gPnt.y = gp_pnt.Y();
      gPnt.z = gp_pnt.Z();
      gPnts[i] = gPnt;
    }
  g_mem->replace_point_list( gPnts, number_points, number_points );

  return CUBIT_SUCCESS;
}

//================================================================================
// Description: This function queries OCC for the edge information
//              needed in facetting a RefEdge.  This information is
//              stored and output in g_mem.
// Author     : Jane Hu
// Date       : 10/26/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Curve* curve_ptr,
                                          GMem* gMem,
                                          double angle_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length ) const
{
  //  get the OCCCurve.
  OCCCurve *occ_curve_ptr = CAST_TO(curve_ptr,OCCCurve);
  assert (gMem);

  if(max_edge_length > get_sme_resabs_tolerance())
  {
    PRINT_WARNING("OCC surface's tessilation doesn't consider edge_length.\n");
    PRINT_WARNING("max_edge_length argument is ignored. \n");
  }
    
  TopoDS_Edge *edge_ptr = occ_curve_ptr->get_TopoDS_Edge();

  return get_graphics( edge_ptr, gMem, angle_tolerance, distance_tolerance, max_edge_length );  
}



CubitStatus OCCQueryEngine::get_graphics( TopoDS_Edge* edge_ptr,
                                          GMem* gMem,
                                          double angle_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length ) const
{
  if (!edge_ptr)
    return CUBIT_FAILURE;

  //do tessellation
  if(distance_tolerance <= 0.0)
    distance_tolerance = 0.2;
  if(angle_tolerance == 0.0)
    angle_tolerance = 5;
  double angle  = CUBIT_PI * angle_tolerance/180;
  BRepAdaptor_Curve acurve(*edge_ptr);
  GCPnts_TangentialDeflection *myMesh = 
        new GCPnts_TangentialDeflection(acurve, angle, 
                                        distance_tolerance);
  if (myMesh == NULL) 
  {
    PRINT_ERROR("Can't tessellate for this curve.\n");
    return CUBIT_FAILURE;
  }
  int num_points = myMesh->NbPoints();

  //! Note: If the polygon is closed, the point of closure is 
  //! repeated at the end of its table of nodes. Thus, on a closed 
  //! triangle the function NbNodes returns 4? 
  GPoint *gPnts= new GPoint[num_points];
  for (int i = 1; i <= num_points ; i ++)
    {
      gp_Pnt gp_pnt = myMesh->Value(i);
      GPoint gPnt;
      gPnt.x = gp_pnt.X();
      gPnt.y = gp_pnt.Y();
      gPnt.z = gp_pnt.Z();
      gPnts[i-1] = gPnt;
    }
  gMem->replace_point_list( gPnts, num_points, num_points );
 
  delete myMesh;
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Given surface and number of point on u and v parametric
//              direction, find the 3-d point locations
// Author     : Jane Hu
// Date       : 10/22/07
//================================================================================
CubitStatus OCCQueryEngine::get_isoparametric_points(Surface* surface,
                                                     int &nu, int &nv,
                                                     GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  if (Tops_face == NULL)
  {
    PRINT_ERROR("This surface is not OCCSurface.");
    return CUBIT_FAILURE;
  }
  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(*Tops_face);

  double u1, u2, v1, v2;
  BRepTools::UVBounds(*Tops_face, u1, u2, v1, v2);

  assert (nu > 1 && nv > 1);
  if(nu <= 1)
    nu = NUM_PTS_UV;
  if(nv <= 1)
    nv = NUM_PTS_UV;
  g_mem = new GMem[nu];
 

  // calculate the nu curves
  double interval_u = (u2-u1)/double(nu-1);
  double interval_v = (v2 - v1)/(double)(nv - 1);

  for (int i = 0; i < nu; i++)
  {
    Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u1 + i * interval_u); 
    g_mem[i].allocate_polylines(nv-1);
    for (int j = 0; j <  nv; j++)
	  {
	    gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval_v);
	    g_mem[i].point_list()[j].x = pnt.X();
	    g_mem[i].point_list()[j].y = pnt.Y();
	    g_mem[i].point_list()[j].z = pnt.Z();
	  }
    g_mem[i].pointListCount = nv;
  }

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_u_isoparametric_points(Surface* surface,
                                                       double v, int&n,
                                                       GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
    {
      PRINT_ERROR("This surface is not OCCSurface.");
      return CUBIT_FAILURE;
    }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  //n must be given to calculate the points.
  if (n <= 1)
    n = NUM_PTS_UV;
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (u2 - u1)/(n -1); 
  
  Handle_Geom_Curve HGeom_curve = HGeom_surface->VIso(v);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
    {
      gp_Pnt pnt = HGeom_curve->Value(u1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
    }
  g_mem->pointListCount = n;

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_v_isoparametric_points(Surface* surface,
                                                       double u, int&n,
                                                       GMem*&g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
    {
      PRINT_ERROR("This surface is not OCCSurface.");
      return CUBIT_FAILURE;
    }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  if (n <= 1)
    n = NUM_PTS_UV;
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (v2 - v1)/(n -1);

  Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
    {
      gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
    }
  g_mem->pointListCount = n;

  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus OCCQueryEngine::transform_vec_position( CubitVector const& ,
                                                    BodySM *,
                                                    CubitVector & ) const
{
  //Nobody is using this function in ACIS yet.
  PRINT_ERROR("Option not supported for OCC based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description:  Calculate for intersection points between a curve and 
//               a segment defined by two points or
//               between two curves or between a curve and a surface.
// Author     :  Jane Hu
// Date       :  10/15/07
//================================================================================
CubitStatus OCCQueryEngine::get_intersections( Curve* curve, 
                                               CubitVector& point1,
                                               CubitVector& point2,
                                               DLIList<CubitVector>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCPoint *pt1 = new OCCPoint(point1);
  OCCPoint *pt2 = new OCCPoint(point2);
  Curve *curve2 = 
    OCCModifyEngine::instance()->make_Curve(pt1, pt2);
  if (curve2 == NULL)
    {
      PRINT_ERROR( "Unable to create OCC EDGE from points\n" );
      return CUBIT_FAILURE;
    }
  
  OCCCurve *occ_curve2 = CAST_TO(curve2, OCCCurve);
  CubitStatus stat = get_intersections(occ_curve, occ_curve2, intscts, bounded, closest);
  delete_solid_model_entities(occ_curve2);
  return stat;
}

CubitStatus OCCQueryEngine::get_intersections( Curve* curve1, 
                                               Curve* curve2,
                                               DLIList<CubitVector>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  //If this function has shortcomes in using BRepExtrema_DistShapeShape,
  //look also at IntTools_EdgeEdge.
  OCCCurve *occ_curve1 =  CAST_TO(curve1, OCCCurve);
  if (occ_curve1 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCCurve *occ_curve2 =  CAST_TO(curve2, OCCCurve);
  if (occ_curve2 == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(
                                            *(occ_curve1->get_TopoDS_Edge()),
                                            *(occ_curve2->get_TopoDS_Edge()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curves.\n");
      return CUBIT_FAILURE;
    }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
    {
      int numPnt = distShapeShape.NbSolution();
      for (int i = 1; i <= numPnt; i++)
	{
	  gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
	  CubitVector cv(aPoint.X(), aPoint.Y(), aPoint.Z());
	  intscts.append(cv);
	}
    }
  return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::get_intersections( Curve* curve, Surface* surface,
                                   DLIList<CubitVector>& intscts,
                                   CubitBoolean bounded )
{
  // There's no effect of bounded =  false. 
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }

  OCCSurface *occ_surface =  CAST_TO(surface, OCCSurface);
  if (occ_surface == NULL)
    {
      PRINT_ERROR("Option not supported for non-occ based geometry.\n");
      return CUBIT_FAILURE;
    }
   
  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(*(occ_curve->get_TopoDS_Edge()),
                                            *(occ_surface->get_TopoDS_Face()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR("Cannot calculate the intersection points for the input curve and surface.\n");
      return CUBIT_FAILURE;
    }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
    {
      int numPnt = distShapeShape.NbSolution();
      for (int i = 1; i <= numPnt; i++)
	{
	  gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
	  CubitVector cv(aPoint.X(), aPoint.Y(), aPoint.Z());
	  intscts.append(cv);
	}
    }
 
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     : Jane Hu
// Date       : 10/30/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_extrema( DLIList<GeometryEntity*> &ref_entity_list,
				const CubitVector *dir1,
				const CubitVector *dir2,
				const CubitVector *dir3,
				CubitVector &extrema,
				GeometryEntity *&extrema_entity_ptr )
{
  //in Acis, the api_entity_extrema is used to calculate "possible 
  //self-intersecting sweeping and to align lofting sections"
  PRINT_ERROR("There's no such call in OCC correponding to Acis call."); 
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     : Jane Hu
// Date       : 10/19/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_entity_distance( GeometryEntity *entity1,
                                        GeometryEntity *entity2,
                                        CubitVector &pos1, CubitVector &pos2,
                                        double &distance )
{
  TopoDS_Shape * shape1;
  TopoDS_Shape * shape2;
  if ((shape1 = get_TopoDS_Shape_of_entity(entity1)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  if( (shape2 = get_TopoDS_Shape_of_entity( entity2 )) == NULL )
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n");
      return CUBIT_FAILURE;
    }

  BRepExtrema_DistShapeShape distShapeShape(*shape1, *shape2);
  //distShapeShape.Perform();
  
  if (!distShapeShape.IsDone())
    {
      PRINT_ERROR( "problem occured getting distance between OCC entities.\n"
		   "       Aborting.\n");
      return CUBIT_FAILURE;
    }

  distance = distShapeShape.Value();
  gp_Pnt pnt1 = distShapeShape.PointOnShape1(1);
  gp_Pnt pnt2 = distShapeShape.PointOnShape2(1);
  pos1 = CubitVector(pnt1.X(), pnt1.Y(), pnt1.Z());
  pos2 = CubitVector(pnt2.X(), pnt2.Y(), pnt2.Z());
  return CUBIT_SUCCESS;
}

TopoDS_Shape* OCCQueryEngine::get_TopoDS_Shape_of_entity(TopologyBridge * entity_ptr)
{
  if (OCCBody * occ_body = CAST_TO( entity_ptr, OCCBody))
    {
      TopoDS_Shape* theShape;
      occ_body->get_TopoDS_Shape(theShape);
      return theShape;
    }

  else if (OCCLump * lump_ptr = CAST_TO( entity_ptr,OCCLump))
    {
      TopoDS_Solid * theSolid = lump_ptr->get_TopoDS_Solid();
      if(theSolid)
	return (TopoDS_Shape*) theSolid; 
      else
	{
	  PRINT_ERROR("OCCLump without TopoDS_Solid at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}
    }

  else if( OCCSurface * surface_ptr = CAST_TO( entity_ptr, OCCSurface))
    {
      TopoDS_Face *theFace = surface_ptr->get_TopoDS_Face();
      if(!theFace)
	{
	  PRINT_ERROR("OCCSurface without TopoDS_Face at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) theFace;
    }

  else if( OCCCurve * curve_ptr = CAST_TO( entity_ptr, OCCCurve))
    {
      TopoDS_Edge *theEdge = curve_ptr->get_TopoDS_Edge();
      if (!theEdge)
	{
	  PRINT_ERROR("OCCCurve without TopoDS_Edge at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) theEdge;
    }

  else if( OCCPoint * point_ptr = CAST_TO( entity_ptr, OCCPoint))
    {
      TopoDS_Vertex *thePoint = point_ptr->get_TopoDS_Vertex(); 
      if (!thePoint)
	{
	  PRINT_ERROR("OCCPoint without TopoDS_Point at %s:%d.\n", __FILE__, __LINE__ );
	  return NULL;
	}

      return (TopoDS_Shape*) thePoint;
    }
  
  PRINT_ERROR("Non-OCC TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return NULL;

}
//===========================================================================
//Function Name: save_temp_geom_file
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary FACET file
//              If file is created, CubitString 'created_file' and
//              'create_file_type' are set
//Author:       Corey Ernst
//Date:         1/18/2003
//===========================================================================

CubitStatus OCCQueryEngine::save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
						 const char *file_name,
						 const CubitString &cubit_version,
						 CubitString &created_file,
						 CubitString &created_file_type)
{
  int size_before = ref_entity_list.size();
  CubitString temp_filename(file_name);
  temp_filename += ".occ";

  ModelExportOptions M_O;
  if( export_solid_model( ref_entity_list, temp_filename.c_str(), OCC_TYPE,
                          cubit_version, M_O ) == CUBIT_FAILURE )
    {
      PRINT_ERROR( "Error occured while trying to save temporary OCC_BASED_GEOMETRY file\n");
      return CUBIT_FAILURE;
    }

  int size_after = ref_entity_list.size();

  if( size_before > size_after )
    {
      created_file +=  temp_filename;
      created_file_type += "OCC";
    }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:export_solid_model
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
						const char* file_name,
						Model_File_Type file_type,
						const CubitString &,
						ModelExportOptions &)
{
  if(  file_type != OCC_TYPE  && 
       file_type != STEP_TYPE &&
       file_type != IGES_TYPE)
    {
      //PRINT_ERROR("The specified file type, %s, is not supported!\n", filetype );
      return CUBIT_FAILURE;
    }
 
/*
  char* name = "write.iges.unit";
  Standard_CString orig_unit;
  char* unit = "M";
  if(strcmp( file_type, "IGES") == 0 && unit != NULL)
  {
    orig_unit = Interface_Static::CVal(name);
    Interface_Static::SetCVal (name, unit); 
  }
*/
  DLIList<OCCBody*>    OCC_bodies;
  DLIList<OCCSurface*> OCC_surfaces;
  DLIList<OCCCurve*>   OCC_curves;
  DLIList<OCCPoint*>   OCC_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
    {
      TopologyBridge* ref_entity_ptr = ref_entity_list.get();
      CubitBoolean handled = CUBIT_TRUE;

      //if it is a Vertex
      if( OCCPoint* pt = CAST_TO( ref_entity_ptr, OCCPoint) )
	OCC_points.append( pt );

      //if it is a Curve
      else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
	OCC_curves.append( curve );
    
      //if it is a surface
      else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
	OCC_surfaces.append( surf );
   
      //if it is a Body
      else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
	OCC_bodies.append( body );

      else
	handled = CUBIT_FALSE;

      if( handled == CUBIT_TRUE )
	{
	  ref_entities_handled.append( ref_entity_ptr );
	  ref_entity_list.change_to(NULL);
	}

      ref_entity_list.step();
    }

  ref_entity_list.remove_all_with_value(NULL);

  int free_body_count = OCC_bodies.size();
  int free_curve_count = OCC_curves.size();
  int free_point_count = OCC_points.size();
  int free_surface_count = OCC_surfaces.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_surface_count == 0 && 
      free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //save the facets (geometry info )
  CubitStatus status;

  //write out topology and attributes
  status = write_topology( file_name, file_type,
                           OCC_bodies, OCC_surfaces,
                           OCC_curves, OCC_points );
/*
  //set the unit back.
  if(strcmp( file_type, "IGES") == 0 && unit != NULL) 
    Interface_Static::SetCVal (name,orig_unit);
*/
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_surface_count || 
      free_curve_count || free_point_count )
    PRINT_INFO( "\nExported:" );

  int flg = 0;

  if( free_body_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_body_count == 1 )
         PRINT_INFO( "%4d OCC Body to %s\n", free_body_count, file_name );
      else
         PRINT_INFO( "%4d OCC Bodies to %s\n", free_body_count, file_name );
    }

  if( free_surface_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_surface_count == 1 )
	PRINT_INFO( "%4d OCC Surface to %s\n", free_surface_count, file_name );
      else
	PRINT_INFO( "%4d OCC Surface to %s\n", free_surface_count, file_name );
    }

  if( free_curve_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
	PRINT_INFO( "%4d OCC Curve to %s\n", free_curve_count, file_name );
      else
	PRINT_INFO( "%4d OCC Curves to %s\n", free_curve_count, file_name );
    }

  if( free_point_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
	PRINT_INFO( "%4d OCC Point to %s\n", free_point_count, file_name );
      else
	PRINT_INFO( "%4d OCC Points to %s\n", free_point_count, file_name );
    }
  PRINT_INFO( "\n" );

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
						char*& p_buffer,
						int& n_buffer_size,
						bool b_export_buffer)
{
  DLIList<OCCBody*>    OCC_bodies;
  DLIList<OCCSurface*> OCC_surfaces;
  DLIList<OCCCurve*>   OCC_curves;
  DLIList<OCCPoint*>   OCC_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
    {
      TopologyBridge* ref_entity_ptr = ref_entity_list.get();
      CubitBoolean handled = CUBIT_TRUE;

      //if it is a Vertex
      if( OCCPoint* pt = CAST_TO( ref_entity_ptr, OCCPoint) )
	OCC_points.append( pt );

      //if it is a Curve
      else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
	OCC_curves.append( curve );
    
      //if it is a surface
      else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
	OCC_surfaces.append( surf );
   
      //if it is a Body
      else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
	OCC_bodies.append( body );

      else
	handled = CUBIT_FALSE;

      if( handled == CUBIT_TRUE )
	{
	  ref_entities_handled.append( ref_entity_ptr );
	  ref_entity_list.change_to(NULL);
	}

      ref_entity_list.step();
    }

  ref_entity_list.remove_all_with_value(NULL);

  int free_body_count = OCC_bodies.size();
  int free_curve_count = OCC_curves.size();
  int free_point_count = OCC_points.size();
  int free_surface_count = OCC_surfaces.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_surface_count == 0 && 
      free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //save the facets (geometry info )
  CubitStatus status;

  //write out topology and attributes
  status = write_topology( p_buffer, n_buffer_size,
			   b_export_buffer,
			   OCC_bodies, OCC_surfaces,
			   OCC_curves, OCC_points);
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_surface_count || 
      free_curve_count || free_point_count )
  {
    if (b_export_buffer) PRINT_INFO( "\nExported:" );
    else PRINT_INFO( "\nSize checked:" );
  }
  int flg = 0;

  if( free_body_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_body_count == 1 )
         PRINT_INFO( "%4d OCC Body to buffer\n", free_body_count );
      else
         PRINT_INFO( "%4d OCC Bodies to buffer\n", free_body_count );
    }

  if( free_surface_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_surface_count == 1 )
	PRINT_INFO( "%4d OCC Surface to buffer\n", free_surface_count );
      else
	PRINT_INFO( "%4d OCC Surface to buffer\n", free_surface_count );
    }

  if( free_curve_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
	PRINT_INFO( "%4d OCC Curve to buffer\n", free_curve_count );
      else
	PRINT_INFO( "%4d OCC Curves to buffer\n", free_curve_count );
    }

  if( free_point_count )
    {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
	PRINT_INFO( "%4d OCC Point to buffer\n", free_point_count );
      else
	PRINT_INFO( "%4d OCC Points to buffer\n", free_point_count );
    }
  PRINT_INFO( "\n" );

  return CUBIT_SUCCESS;
}

void OCCQueryEngine::body_attributes_for_writing(DLIList<OCCBody*> &OCC_bodies,
                                  BRep_Builder &B,
                                  TopoDS_Compound &Co,
                                  DLIList<OCCLump*> &single_lumps,
                                  DLIList< DLIList<CubitSimpleAttrib>*> &lists)
{
  //Add every shape to the compound
  OCCLump* lump = NULL;
  DLIList<CubitSimpleAttrib> body_csa_list;

  for (int i = 0; i < OCC_bodies.size(); i++)
  {
    OCCBody* body = OCC_bodies.get_and_step();
    TopoDS_Compound *shape = body->get_TopoDS_Shape();
    if (shape == NULL || shape->IsNull()) //single lump or sheet or shell body
    {
       DLIList<OCCSurface*> surfaces = body->my_sheet_surfaces();
       DLIList<OCCShell*> shells = body->shells();
       DLIList<Lump*> lumps = body->lumps();
       if(surfaces.size() == 1)
         B.Add(Co,*(surfaces.get()->get_TopoDS_Face()));
       else if (shells.size() == 1)
         B.Add(Co,*(shells.get()->get_TopoDS_Shell()));
       else
       {
         lump = CAST_TO(lumps.get(), OCCLump);
         B.Add(Co, *(lump->get_TopoDS_Solid()));
         //if body has attributes, add them to the solid.
         DLIList<CubitSimpleAttrib> csa_list;
         body->get_simple_attribute(csa_list);
         body_csa_list.clean_out();
         for(int i = 0; i < csa_list.size(); i++)
         {
           CubitSimpleAttrib body_csa = csa_list.get_and_step();
           CubitString pre_fix("#SINGLELUMP%");
           pre_fix += CubitString::number(i);

           body_csa.string_data_list().insert(body_csa.string_data_list().begin(), pre_fix);
           lump->append_simple_attribute_virt(body_csa);
           body_csa_list.append(body_csa);
         }
         if(csa_list.size() > 0)
         {
           single_lumps.append(lump);
           lists.append(new DLIList<CubitSimpleAttrib>(body_csa_list));
         }
       }
       continue;
    }
    B.Add(Co, *shape);
  }
} 
//===========================================================================
//Function Name:write_topology
//Member Type:  PUBLIC
//Description:  function called for write out temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus
OCCQueryEngine::write_topology( const char* file_name,
                                Model_File_Type file_type,
                                DLIList<OCCBody*> &OCC_bodies,
                                DLIList<OCCSurface*> &OCC_surfaces,
                                DLIList<OCCCurve*> &OCC_curves,
                                DLIList<OCCPoint*> &OCC_points )
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  DLIList<OCCLump*> single_lumps;
  DLIList< DLIList<CubitSimpleAttrib>*> lists;
  OCCLump* lump = NULL;

  if(OCC_bodies.size() > 0)
    body_attributes_for_writing(OCC_bodies, B, Co, single_lumps, lists);

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  if(file_type == OCC_TYPE)
  {
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    CubitBoolean result = Write(Co, const_cast<char*>(file_name),label);
    //remove the body attributes from lump
    for (int i = 0; i < single_lumps.size(); i++)
    {
      lump = single_lumps.get_and_step();
      DLIList<CubitSimpleAttrib>* p_csas = lists.get_and_step();
      for(int j = 0 ; j < p_csas->size(); j ++)
      {
        const CubitSimpleAttrib& csa = p_csas->get_and_step();
        lump->remove_simple_attribute_virt(csa);  
      }
      delete p_csas;
    }
    if(!result)
      return CUBIT_FAILURE;
  } 

#ifdef HAVE_OCC_STEP
  else if(file_type == STEP_TYPE)
  {
    STEPControl_Writer writer;
    writer.Model( Standard_True);
    writer.Transfer(Co, STEPControl_AsIs );
    IFSelect_ReturnStatus stat = writer.Write( (char*) file_name);
    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    }
  }
#endif
#ifdef HAVE_OCC_IGES
 else if (file_type == IGES_TYPE) // IGES file
  {
    IGESControl_Writer writer;
    writer.AddShape(Co);
    writer.ComputeModel();
    Standard_Boolean  stat = writer.Write( (char*) file_name);
    if (!stat )
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    }
  }
#endif
  else {
    PRINT_ERROR("File formats other than OCC, STEP and IGES are  not supported by OCC.\n");
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::write_topology( char*& p_buffer,
				int& n_buffer_size,
				bool b_export_buffer,
				DLIList<OCCBody*> &OCC_bodies,
				DLIList<OCCSurface*> &OCC_surfaces,
				DLIList<OCCCurve*> &OCC_curves,
				DLIList<OCCPoint*> &OCC_points)
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  DLIList<OCCLump*> single_lumps;
  DLIList< DLIList<CubitSimpleAttrib>*> lists;
  OCCLump* lump = NULL;

  if(OCC_bodies.size() > 0)
    body_attributes_for_writing(OCC_bodies, B, Co, single_lumps, lists);

  for (i = 0; i < OCC_surfaces.size(); i++)
    {
      TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
      B.Add(Co, *face);
    }

  //Add standalone wires to the export BRep file
  for (i = 0; i < WireList->size(); i++)
    {
      TopoDS_Wire *wire = WireList->get_and_step()->get_TopoDS_Wire();
      B.Add(Co, *wire);
    }

  for (i = 0; i < OCC_curves.size(); i++)
    {
      TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
      B.Add(Co, *edge);
    }

  for (i = 0; i < OCC_points.size(); i++)
    {
      TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
      B.Add(Co, *vertex);
    }
 
  //if(strcmp(file_type, "OCC") == 0)
  //{
    TDF_Label label;
    if(EXPORT_ATTRIB)
      label = mainLabel;

    if(!Write(Co, p_buffer, n_buffer_size, b_export_buffer, label))
      return CUBIT_FAILURE;

    //remove the body attributes from lump
    for (int i = 0; i < single_lumps.size(); i++)
    {
      lump = single_lumps.get_and_step();
      DLIList<CubitSimpleAttrib>* p_csas = lists.get_and_step();
      for(int j = 0 ; j < p_csas->size(); j ++)
      {
        const CubitSimpleAttrib& csa = p_csas->get_and_step();
        lump->remove_simple_attribute_virt(csa);
      }
      delete p_csas;
    }

  return CUBIT_SUCCESS;
}

CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
                                   const Standard_CString File,
                                   TDF_Label label) 
{
  ofstream os;
  os.open(File, ios::out);
  if (!os.rdbuf()->is_open()) return Standard_False;
  
  CubitBoolean isGood = (os.good() && !os.eof());
  if(!isGood)
    return isGood;

  OCCShapeAttributeSet SS;
  SS.Add(Sh);

  os << "DBRep_DrawableShape\n";  // for easy Draw read
  SS.Write(os);
  isGood = os.good();
  if(isGood )
    SS.Write(Sh,os,&label);
  os.flush();
  isGood = os.good();
  os.close();
  isGood = os.good() && isGood;

  return isGood;
}

CubitBoolean OCCQueryEngine::Write(const TopoDS_Shape& Sh,
				   char*& pBuffer,
				   int& n_buffer_size,
				   bool b_write_buffer,
                                   TDF_Label label)
{
  // make buffer as ouput stream
  std::stringbuf sb;
  std::iostream os(&sb);
  OCCShapeAttributeSet SS;
  
  // write to output stream
  SS.Add(Sh);
  os << "DBRep_DrawableShape\n";  // for easy Draw read
  SS.Write(os);
  CubitBoolean isGood = os.good();
  if (!isGood) return isGood;
  SS.Write(Sh,os,&label);
  isGood = os.good();
  if (!isGood) return isGood;
  
  n_buffer_size = os.rdbuf()->pubseekoff(0, std::ios_base::end, std::ios::out);

  // get real geometries from output stream to buffer
  if (b_write_buffer) os.read(pBuffer, n_buffer_size);
  
  return CUBIT_TRUE;
}
                                   
CubitBoolean OCCQueryEngine::Read(TopoDS_Shape& Sh,
                                  const Standard_CString File,
                                  TDF_Label label)
{
  ifstream in( File );
  if (in.fail()) {
    PRINT_INFO("%s: Cannot open file", File );
    return CUBIT_FAILURE;
  }

  OCCShapeAttributeSet SS;
  SS.Read(in, CUBIT_TRUE);
  int nbshapes = SS.NbShapes();
  if(!nbshapes) return CUBIT_FALSE;
  SS.Read(Sh,in,nbshapes, &label);
  return CUBIT_TRUE;
}
                                   
CubitBoolean OCCQueryEngine::Read(TopoDS_Shape& Sh,
				  const char* pBuffer,
				  const int n_buffer_size,
                                  TDF_Label label)
{
  // make buffer as input stream
  std::stringbuf sb;
  std::iostream is(&sb);
  is.write(pBuffer, n_buffer_size);
  
  // read from input stream
  OCCShapeAttributeSet SS;
  SS.Read(is, false);
  int nbshapes = SS.NbShapes();
  if (!nbshapes) return CUBIT_FALSE;
  SS.Read(Sh, is, nbshapes, &label);

  return CUBIT_TRUE;
}

CubitStatus
OCCQueryEngine::import_temp_geom_file(FILE* file_ptr,
                                      const char* file_name,
                                      Model_File_Type file_type,
                                      DLIList<TopologyBridge*> &bridge_list )
{
  ModelImportOptions M_O;
  return import_solid_model( file_name, file_type, bridge_list, M_O );
}

//===========================================================================
//Function Name:import_solid_model
//Member Type:  PUBLIC
//Description:  function called for read in temporary Brep file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus OCCQueryEngine::import_solid_model(
					       const char* file_name ,
					       Model_File_Type file_type,
					       DLIList<TopologyBridge*> &imported_entities,
                                               ModelImportOptions &import_options)
{
  TopoDS_Shape *aShape = new TopoDS_Shape;
  
  //BRep_Builder aBuilder;
  if(file_type == OCC_TYPE)
  {
    Standard_Boolean result = Read(*aShape, (char*) file_name, mainLabel);
    if (result==0) return CUBIT_FAILURE;
  }
#ifdef HAVE_OCC_STEP 
  else if (file_type == STEP_TYPE)
  {
    STEPControl_Reader reader;
/*
    char* name = "xstep.cascade.unit";
    char* unit = "M"; 
    Standard_CString orig_unit = Interface_Static::CVal(name);
    Interface_Static::SetCVal (name, unit); //this set is good for both step
                                            // and iges files.
*/
    IFSelect_ReturnStatus stat = reader.ReadFile( (char*) file_name);
/*
    //set the unit back.
    Interface_Static::SetCVal (name,orig_unit);
*/
    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    } 
    reader.TransferRoots();
    *aShape = reader.OneShape(); 
  }
#endif
#ifdef HAVE_OCC_IGES
  else if(file_type == IGES_TYPE)
  {
    IGESControl_Reader reader;
/*
    char* name = "xstep.cascade.unit";
    char* unit = "M";
    Interface_Static::SetCVal (name, unit); //this set is good for both step
                                            // and iges files.
*/
    const Standard_CString string1 = file_name;
    IFSelect_ReturnStatus stat = reader.ReadFile( string1);

//    Interface_Static::SetCVal (name, "MM");

    if (stat  != IFSelect_RetDone)
    {
       PRINT_INFO("%s: Cannot open file", file_name );
       return CUBIT_FAILURE;
    } 
    reader.TransferRoots(); 
    *aShape = reader.OneShape();
  } 
#endif
  else 
  {
    PRINT_ERROR("File formats other than OCC, STEP and IGES are not supported by OCC\n");
    return CUBIT_FAILURE;
  }
    
  //All read in shapes are wrapped inside a compound shape. Ignore this one
  TopoDS_Iterator it(*aShape);
  if(aShape->ShapeType() != TopAbs_COMPOUND)
  {
    imported_entities += populate_topology_bridge(*aShape);
    return CUBIT_SUCCESS;
  }

  for(;it.More();it.Next())
  {
    TopoDS_Shape shape = it.Value();
    imported_entities += populate_topology_bridge(shape);
  }

  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::import_solid_model(DLIList<TopologyBridge*> &imported_entities,
					       const char* pBuffer,
					       const int n_buffer_size)
{
  TopoDS_Shape *aShape = new TopoDS_Shape;
  Standard_Boolean result = Read(*aShape, pBuffer, n_buffer_size, mainLabel);
  if (result==0) return CUBIT_FAILURE;
  
  //All read in shapes are wrapped inside a compound shape. Ignore this one
  TopoDS_Iterator it(*aShape);
  if(aShape->ShapeType() != TopAbs_COMPOUND)
  {
    imported_entities += populate_topology_bridge(*aShape);
    return CUBIT_SUCCESS;
  }

  for(;it.More();it.Next())
  {
    TopoDS_Shape shape = it.Value();
    imported_entities += populate_topology_bridge(shape);
  }

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:populate_topology_bridge
//Member Type:  PUBLIC
//Description:  function called for populating topology bridge for OCC entity 
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

DLIList<TopologyBridge*> OCCQueryEngine::populate_topology_bridge(TopoDS_Shape& aShape)
{
  DLIList<TopologyBridge*> tblist;
  // suitable to populate for a  TopoDS_Compound shape.
  if ( aShape.ShapeType() == TopAbs_COMPOUND)
    tblist.append(populate_topology_bridge(TopoDS::Compound(aShape)));

  else if(aShape.ShapeType() == TopAbs_SOLID)
  {
    Lump *lump = 
    populate_topology_bridge(TopoDS::Solid(aShape), CUBIT_TRUE);
    tblist.append(CAST_TO(lump, OCCLump)->get_body());
  }

  else if(aShape.ShapeType() == TopAbs_SHELL)
  {
    OCCShell* shell =
       populate_topology_bridge(TopoDS::Shell(aShape), CUBIT_TRUE);
    tblist.append(shell->my_body());
  }

  else if(aShape.ShapeType() == TopAbs_FACE)
  {
    Surface* face =
      populate_topology_bridge(TopoDS::Face(aShape), CUBIT_TRUE);
    if(face)
      tblist.append(CAST_TO(face, OCCSurface)->my_body());
  }

  else if(aShape.ShapeType() == TopAbs_WIRE)
    populate_topology_bridge(TopoDS::Wire(aShape), CUBIT_TRUE);

  else if(aShape.ShapeType() == TopAbs_EDGE)
    tblist.append(populate_topology_bridge(TopoDS::Edge(aShape),CUBIT_TRUE));

  else if(aShape.ShapeType() == TopAbs_VERTEX)
    tblist.append(populate_topology_bridge(TopoDS::Vertex(aShape), CUBIT_TRUE));
  else
    PRINT_ERROR("Wrong topology type is given. \n");
  tblist.remove_all_with_value(NULL);
  return tblist;
}

BodySM* OCCQueryEngine::populate_topology_bridge(const TopoDS_Compound& aShape)
{
  if(aShape.IsNull())
    return (BodySM*)NULL;
  OCCBody *body = (OCCBody*)NULL;
  if (!OCCMap->IsBound(aShape) || 
       OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
    {
      //check to see if this compound has only one lump which is already in 
      //in another body. Unite operation will return a one lump compound.
      //Check also if this compound has only shells. which is or has faces that
      //are already in another body. Unite faces will return such shell.
      TopExp_Explorer Ex;
      int num_lumps = 0, num_shells = 0, num_faces = 0;
      TopoDS_Solid solid;
      for (Ex.Init(aShape, TopAbs_SOLID); Ex.More(); Ex.Next()) 
      {
        num_lumps++;
        solid = TopoDS::Solid(Ex.Current());
      }

      TopoDS_Shell shell;
      for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
      {
        num_shells++;
        shell = TopoDS::Shell(Ex.Current());
      }

      TopoDS_Face face;
      for (Ex.Init(aShape, TopAbs_FACE, TopAbs_SHELL);  Ex.More(); Ex.Next())
      {
        num_faces++;
        face = TopoDS::Face(Ex.Current());
      }

      if(num_faces + num_shells + num_lumps == 1)
      {
        if (num_faces  == 1 && (!OCCMap->IsBound(face) ||
                      OccToCGM->find(OCCMap->Find(face)) == OccToCGM->end()))
        {
          Surface* surface = populate_topology_bridge(face, CUBIT_TRUE);
          return CAST_TO(surface, OCCSurface)->my_body();
        }
        else if (num_shells == 1 && (!OCCMap->IsBound(shell)||
                  OccToCGM->find(OCCMap->Find(shell)) == OccToCGM->end()))
        {
          OCCShell* occ_shell = populate_topology_bridge(shell, CUBIT_TRUE);
          return occ_shell->my_body();
        }
        else if( num_lumps == 1 && (!OCCMap->IsBound(solid) ||
                  OccToCGM->find(OCCMap->Find(solid)) == OccToCGM->end()))
        {
          Lump* lump= populate_topology_bridge(solid, CUBIT_TRUE);
          return CAST_TO(lump, OCCLump)->get_body();
        }
        else //find existing body
        {
          int k;
          if(num_lumps == 1)
          {
            k = OCCMap->Find(solid);
            OCCLump* lump = (OCCLump*)(OccToCGM->find(k))->second;
            body = CAST_TO(lump->get_body(), OCCBody);
          }
          else if (num_shells == 1)
          {
            k = OCCMap->Find(shell);
            OCCShell* occ_shell = (OCCShell*)(OccToCGM->find(k))->second;
            body = occ_shell->my_body();
          }
          else
          {
            k = OCCMap->Find(face);
            OCCSurface* occ_surface = (OCCSurface*) (OccToCGM->find(k))->second;
            body = occ_surface->my_body();
          }
        }
      } 

      else
      {
        TopoDS_Compound *comsolid = new TopoDS_Compound;
        *comsolid = aShape;
        body = new OCCBody(comsolid);
        int current_id;
        if(!OCCMap->IsBound(aShape))
        {
          iTotalTBCreated++;
          current_id = iTotalTBCreated;
          OCCMap->Bind(aShape, current_id);
        }

        else
          current_id = OCCMap->Find(aShape);

        OccToCGM->insert(valType(current_id,
                             (TopologyBridge*)body));
        BodyList->append(body);
      }
    }
    else
    {
      int k = OCCMap->Find(aShape);
      body = (OCCBody*)(OccToCGM->find(k))->second;
      TopoDS_Compound compound = aShape;
      body->set_TopoDS_Shape(compound);
    }

  TopExp_Explorer Ex;
  DLIList<Lump*> lumps;
  for (Ex.Init(aShape, TopAbs_SOLID); Ex.More(); Ex.Next())
  {
     TopoDS_Shape sh = Ex.Current();
     Lump* lump = populate_topology_bridge(TopoDS::Solid(sh));
     lumps.append(lump);
     CAST_TO(lump, OCCLump)->add_body(body);
     //add the solid shape into map
     TopoDS_Shape parent = aShape;
     int current_id;
     add_shape_to_map(sh, parent, current_id);
     if(!OCCMap->IsBound(sh) ||
        OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
        OccToCGM->insert(valType(current_id, (TopologyBridge*)lump));
  }
  body->lumps(lumps);

  DLIList<OCCShell*> shells;
  for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID);Ex.More(); Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();
    OCCShell * shell = populate_topology_bridge(TopoDS::Shell(sh)); 
    OCCLump* lump = shell->my_lump();
    if(lump == (OCCLump*)NULL)
      lump = new OCCLump(NULL, NULL, shell);
    lumps.append(lump);
    lump->add_body(body);
    shell->set_body(body);
    shell->set_lump(lump);  
    shells.append(shell);
    TopoDS_Shape parent = aShape;
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)shell));
  }
  body->shells(shells);
  
  DLIList<OCCSurface*> surfaces;
  for (Ex.Init(aShape, TopAbs_FACE, TopAbs_SHELL);Ex.More(); Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();
    Surface* face = populate_topology_bridge(TopoDS::Face(sh));
    OCCSurface *surface = CAST_TO(face, OCCSurface);
    OCCShell* shell = surface->my_shell();
    if (shell == (OCCShell*) NULL)
      shell = new OCCShell(NULL, surface);
    OCCLump* lump = surface->my_lump();
    if(lump == (OCCLump*) NULL)
      lump = new OCCLump(NULL, surface);
    lumps.append(lump);
    lump->add_body(body);
    surface->set_body(body);
    surface->set_lump(lump);
    surface->set_shell(shell);
    shell->set_body(body);
    shell->set_lump(lump);
    TopoDS_Shape parent = aShape;
    surfaces.append(surface);
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)face));

  } 
  body->set_sheet_surfaces(surfaces);
  return body;
}

Lump* OCCQueryEngine::populate_topology_bridge(const TopoDS_Solid& aShape,
		 			       CubitBoolean build_body)
{
  if(aShape.IsNull())
    return (Lump*)NULL;

  OCCLump *lump = NULL;
  OCCBody *body = NULL;
  int current_lump_number = 0;
  if (!OCCMap->IsBound(aShape) || 
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
  {
    TopoDS_Solid *posolid =  new TopoDS_Solid;
    *posolid = aShape;
    lump = new OCCLump(posolid);
    if (build_body)
    {
      if(OCCMap->IsBound(aShape))
        current_lump_number = OCCMap->Find(aShape);
      else
      {
        iTotalTBCreated++;
        current_lump_number = iTotalTBCreated;
      }
      body = new OCCBody(NULL, NULL, NULL, lump);
      DLIList<CubitSimpleAttrib> csa_list;
      lump->get_simple_attribute(csa_list);
      //if there's body attribute, append it to body and delete it from lump.
      for(int i = 0; i < csa_list.size(); i++)
      {
        const CubitSimpleAttrib& csa = csa_list.get_and_step();
        CubitString type = csa.string_data_list()[0];
        CubitString subtype = type.substr(0,12);
        if(subtype == "#SINGLELUMP%")
        {  
          lump->remove_simple_attribute_virt(csa);
          CubitSimpleAttrib copy = csa;
          copy.string_data_list().erase(copy.string_data_list().begin());
          body->append_simple_attribute_virt(copy);
        }
      }
      BodyList->append(body);
      lump->add_body(body);
    }
  }
  else 
  {
    int k = OCCMap->Find(aShape);
    lump = (OCCLump*)(OccToCGM->find(k))->second;
    lump->set_TopoDS_Solid(aShape);
    body = static_cast<OCCBody*>(lump->get_body());
    TopoDS_Shape *b_shape = NULL;
    if(body)
      body->get_TopoDS_Shape(b_shape);
    if (!body || b_shape == NULL)
    { 
      if(body)
        BodyList->remove(body);
      body = new OCCBody(NULL, NULL, NULL, lump);
      BodyList->append(body);
      lump->add_body(body);
    }
  }

  TopoDS_Compound *shape = NULL;
  if(body)
    shape = body->get_TopoDS_Shape();

  if(build_body && OCCMap->IsBound(aShape) && shape && !shape->IsNull())
  {
    PRINT_ERROR("Single lump body shouldn't have Compound shape.\n");
    return (Lump*) NULL;
  }

  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_SHELL); Ex.More(); Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();
    OCCShell* shell = populate_topology_bridge(TopoDS::Shell(sh));
    shell->set_lump(lump);
    shell->set_body(body);
    TopoDS_Shape parent = aShape;
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)shell));
  } 

  if(build_body && (!OCCMap->IsBound(aShape)||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end()))
  {
    if (!OCCMap->IsBound(aShape))
      OCCMap->Bind(aShape, current_lump_number);
    OccToCGM->insert(valType(current_lump_number,
                       (TopologyBridge*)lump));
  }
  return lump;
}

void OCCQueryEngine::add_shape_to_map(TopoDS_Shape& sh,
                                      TopoDS_Shape& aShape, /*In, parent */
                                      int &current_id /*Out*/)
{
  //add the shape into map
  if(sh.IsNull())
    return;

  if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
  {
     DLIList<TopoDS_Shape*> list;
     //find the sh shape without aShape's location.
     TopLoc_Location L = aShape.Location();
     TopoDS_Iterator it(aShape, Standard_True, Standard_False);
     TopoDS_Shape bare_shape;
     for(; it.More(); it.Next())
     {
       TopoDS_Shape test_shape = it.Value();
       test_shape.Move(L);
       if(sh.IsEqual(test_shape))
       {
         bare_shape = it.Value();
         break;
       }
     }
     if(!OCCMap->IsBound(sh))
     {
       CubitBoolean shape_found = false;
       if(OCCMap->IsBound(bare_shape))
       {
         current_id = OCCMap->Find(bare_shape);
         //There are two possiblities when coming here:
         //1. After OCCAttribute binds the bare_shape but not binds the topo.
         //2. The bare_shape is bound because it binds to a different topo. 
         if( OccToCGM->find(current_id) == OccToCGM->end() )        
         {
           OCCMap->UnBind(bare_shape);
           std::map<int, TDF_Label>::iterator it = 
              Shape_Label_Map->find(current_id);
           if(it != Shape_Label_Map->end())
           {
             TDF_Label aLabel = (*it).second;
             Handle_TDataXtd_Shape attr_shape = TDataXtd_Shape::Set(aLabel, sh);
           }
           shape_found = true;
         }
       }
       if(!shape_found)
       {
         iTotalTBCreated++;
         current_id = iTotalTBCreated;
       }
       OCCMap->Bind(sh, current_id);
    }
    else
      current_id = OCCMap->Find(sh);
  }
  else
    current_id = OCCMap->Find(sh);
}

OCCShell* OCCQueryEngine::populate_topology_bridge(const TopoDS_Shell& aShape,
						   CubitBoolean standalone)
{
  if(aShape.IsNull())
    return (OCCShell*)NULL;
  OCCShell *shell ;
  if (!OCCMap->IsBound(aShape) ||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
  {
    if(standalone)
    {
      //check if just has one Face,if so, don't make new shell.
      TopExp_Explorer Ex;
      int num_faces = 0;
      TopoDS_Face topo_face;
      for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next())
      {
        topo_face = TopoDS::Face(Ex.Current());
        num_faces++;
      }
      if(num_faces == 1)
      {
        Surface* face = populate_topology_bridge(topo_face, standalone);
        return CAST_TO(face, OCCSurface)->my_shell();
      }
    }
    TopoDS_Shell *poshell = new TopoDS_Shell;
    *poshell = aShape;
    shell = new OCCShell(poshell);
    shell->set_body(NULL);
    shell->set_lump(NULL); 
    if(standalone)
    {
      OCCLump* lump = new OCCLump(NULL, NULL, shell);
      OCCBody* body = new OCCBody(NULL, NULL, shell);
      BodyList->append(body);
      shell->set_body(body);
      shell->set_lump(lump);
      int k;
      if(!OCCMap->IsBound(aShape))
      {
        iTotalTBCreated++;
        k = iTotalTBCreated;
        OCCMap->Bind(*poshell, k);
      }
      else
        k = OCCMap->Find(aShape); 
      OccToCGM->insert(valType(k, (TopologyBridge*)shell));
    }
  }
  else
  {
    int k = OCCMap->Find(aShape);
    shell = (OCCShell*)(OccToCGM->find(k))->second;
    shell->set_TopoDS_Shell(aShape);
  }

  DLIList<OCCSurface*> memberSurfaces;
  TopExp_Explorer Ex;
  for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();
    TopoDS_Face topo_face = TopoDS::Face(sh);
    Surface* face =
      populate_topology_bridge(topo_face, CUBIT_FALSE);
    
    TopoDS_Shape parent = aShape; 
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)face));

    if(!face)
      continue;
    OCCSurface *occ_surface = CAST_TO(face, OCCSurface);
    //check if surface was a sheet surface, delete it if so
    if(occ_surface->my_shell() != NULL && occ_surface->my_shell() != shell)
    {
       //if Sheet_body, delete this sheet body
       OCCBody* occ_body = occ_surface->my_body();
       if(occ_body != NULL)
       {
          delete_solid_model_entities(occ_body);
          face =
            populate_topology_bridge(topo_face, CUBIT_FALSE);
   
          if(!face)
            continue;
          occ_surface = CAST_TO(face, OCCSurface);
          shell->set_sheet_surface(occ_surface);
       }
    }

    memberSurfaces.append(occ_surface);
    occ_surface->set_shell(shell);
  }

  shell->setMemberSurfaces(memberSurfaces);
  return shell;
}

Surface* OCCQueryEngine::populate_topology_bridge(const TopoDS_Face& aShape,
                                                  CubitBoolean build_body)
{
  if(aShape.IsNull())
    return (Surface*)NULL;
  OCCSurface *surface = NULL;
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(aShape, myProps);
  double area = myProps.Mass();
  double tol = get_sme_resabs_tolerance();
  if(area < tol * tol && area > 0.)
    PRINT_WARNING("Generated a sliver surface. \n");
  
  else if (area < 0.0)
    PRINT_WARNING("Generated a negative area surface. \n");

  if (!OCCMap->IsBound(aShape) ||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
  {
    TopoDS_Face *poface = new TopoDS_Face;
    *poface = aShape;
    surface = new OCCSurface(poface);
    SurfaceList->append(surface);
    surface->set_body(NULL);
    surface->set_lump(NULL);
    surface->set_shell(NULL);
    if(build_body)
    {
      OCCShell* shell = new OCCShell(NULL, surface);
      OCCLump* lump = new OCCLump(NULL, surface);
      OCCBody* body = new OCCBody(NULL, surface);
      surface->set_body(body);
      surface->set_lump(lump);
      surface->set_shell(shell);
      shell->set_body(body);
      shell->set_lump(lump);
      BodyList->append(body);
      int k;
      if(!OCCMap->IsBound(aShape))
      {
        iTotalTBCreated++;
        k = iTotalTBCreated;
        OCCMap->Bind(*poface, iTotalTBCreated);
      }
      else
        k = OCCMap->Find(aShape);
      OccToCGM->insert(valType(k, (TopologyBridge*)surface));
    }
  } 

  else 
  {
    int k = OCCMap->Find(aShape);
    surface = (OCCSurface*)(OccToCGM->find(k))->second;
    TopoDS_Face aFace(aShape);
    surface->set_TopoDS_Face(aFace);
  }

  TopExp_Explorer Ex;
  for (Ex.Init(aShape.Oriented(TopAbs_FORWARD), TopAbs_WIRE); Ex.More();
      Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();

    OCCLoop* loop = populate_topology_bridge(TopoDS::Wire(sh));
/*
    if( aShape.Orientation() == TopAbs_REVERSED )
    {
      DLIList<OCCCoEdge*> coedges_new = loop->coedges();
      coedges_new.reverse();
      //Reverse all coedges' senses.
      for (int i = 0; i < coedges_new.size(); i ++)
      {
        OCCCoEdge* coedge = coedges_new.get_and_step();
        CubitSense sense = (coedge->sense() == CUBIT_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
        coedge->set_sense(sense);
      }
      loop->coedges(coedges_new);
    } 
*/
    TopoDS_Shape parent = aShape;
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)loop));
  } 

  return surface;
}

OCCLoop* OCCQueryEngine::populate_topology_bridge(const TopoDS_Wire& aShape,
						  CubitBoolean standalone)
{
  if(aShape.IsNull())
    return (OCCLoop*)NULL;

  BRepTools_WireExplorer Ex;

  OCCLoop *loop ;
  if (!OCCMap->IsBound(aShape) ||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
  {
    TopoDS_Wire *powire = new TopoDS_Wire;
    *powire = aShape;
    loop = new OCCLoop(powire);
    if(standalone)
    {
      int k;
      if(!OCCMap->IsBound(aShape))
      {
        iTotalTBCreated++;
        k = iTotalTBCreated;
        OCCMap->Bind(aShape, iTotalTBCreated);
      }
      else
        k = OCCMap->Find(aShape);
      OccToCGM->insert(valType(k, (TopologyBridge*)loop));
      WireList->append(loop);
    }
  }
  else
  {
    int k = OCCMap->Find(aShape);
    loop = (OCCLoop*)(OccToCGM->find(k))->second;
    loop->set_TopoDS_Wire(aShape);
  }

  CubitVector v;

  DLIList <OCCCoEdge*> coedges_old, coedges_new;
  coedges_old = loop->coedges();

  for (Ex.Init(aShape); Ex.More(); Ex.Next())
  {
    const TopoDS_Edge& anEdgeForPop = Ex.Current();
    TopoDS_Shape crv = anEdgeForPop;
    Curve* curve = populate_topology_bridge(anEdgeForPop);
    if(!curve)
      continue;
    TopoDS_Shape parent = aShape;
    int current_id;
    add_shape_to_map(crv, parent, current_id);
    if(!OCCMap->IsBound(crv) ||
     OccToCGM->find(OCCMap->Find(crv)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)curve));

    OCCCurve *occ_curve = CAST_TO(curve, OCCCurve);
    DLIList<OCCLoop*> loops = occ_curve->loops();
    CubitBoolean exist = CUBIT_FALSE;
    OCCCoEdge * coedge = NULL;
    int size = coedges_old.size();
    CubitSense sense = (crv.Orientation() == TopAbs_FORWARD ?
        CUBIT_FORWARD : CUBIT_REVERSED);

    //for the cylinder side face, there are 4 coedges, 2 of them are seam
    //edges and should have opposite sense.
    for(int i = 0; i < coedges_new.size(); i++)
    {
      coedge =  coedges_new.get_and_step();
      Curve* test_c = coedge->curve();
      if(test_c == curve)
      {       
        if(sense == coedge->sense())
          sense = (sense == CUBIT_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
        
        break;    
      }
    }

    for(int i = 0; i < size; i++)
    {
      coedge = coedges_old.get_and_step();
      if(coedge->curve() == curve && coedge->sense() ==  sense)
      {
        coedge->set_mark(1);
        exist = CUBIT_TRUE;
        coedge->set_sense(sense);
        coedges_new.append(coedge);
        break;
      }
    }   
    
    if(!exist )
    {
      //search through the curve loops
      for(int i = 0; i < loops.size() ; i++)
      {
        OCCLoop* occ_loop = loops.get_and_step();
        TopoDS_Wire* wire = occ_loop->get_TopoDS_Wire();
        if (size > 0 && (wire->IsNull() || !OCCMap->IsBound(*wire)))
        { 
          DLIList<OCCCoEdge*> coedge_list = occ_loop->coedges();
          for(int j = 0; j < coedge_list.size(); j++)
          {
            OCCCoEdge * test_coedge = coedge_list.get_and_step();
            occ_loop->remove_coedge(test_coedge);
            occ_curve->remove_loop(occ_loop);
            delete test_coedge;
          }
          if(!wire->IsNull())
            wire->Nullify();
          WireList->remove(occ_loop);
        }
      }

      if(occ_curve->loops().size() == 2)
      {
        // If there are some assumptions made about manifold geometry . . .
        //there must be a loop which doesn't have this curve anymore
        //this is been found in subtract cases while one solid becomes
        //2 solid, and one face becomes two faces. One face uses/updates
        //the old face while the other face generates face and wire from
        //new. However, in order to uses the curves, the old curve is kept
        //as possible, so curve's loops get kept, but since it's going to 
        //associate with new loop, the old loop should be removed from the 
        //loop list.
        DLIList<OCCLoop*> old_loops = occ_curve->loops();
        for (int i = 0; i < 2; i++)
        {
          int found = 0;
          OCCLoop* old_loop = old_loops.get_and_step();
          TopoDS_Wire* thewire = old_loop->get_TopoDS_Wire();
          if (thewire < (void*) 0x1000 || thewire->IsNull())
          {
            found = 1;
            break;
          }
          DLIList<OCCCoEdge*> test_coedges = old_loop->coedges();
          for(int j = 0; j < test_coedges.size() ; j++)
          {
            if(test_coedges.get()->curve() != curve)
              test_coedges.step();
            else
            {
              found = 1;
              break;
            }
          }
          if(!found)
            occ_curve->remove_loop(old_loop); 
        }  
      }
      //for unite case, it's possible that the a curve has 3 coedges because
      //opencascade do not perform unite on surfaces.
      coedge = new OCCCoEdge( curve, loop, sense);
      coedges_new.append(coedge);
      occ_curve->add_loop(loop);
    }
  }

  //Double check edge sense to make sure it consists with loop direction
  //coedges size = 2 case is checked in the face level so the face normal will
  //be considered.
  if(coedges_new.size() > 2)
  {
    OCCCoEdge* coedge = coedges_new.get_and_step();
    double    d1, d1_, d2, d2_;
    OCCCoEdge* next_coedge = coedges_new.get();
    //because of tolerance issue, now check current loop end vertices distance
    //compared with reversed loop's end vertices distance, which ever shorter
    //will be the direction.
    //current loop's parameter not using " _ ", reversed loop's parameter
    //uses " _ ".
    if (coedge->sense() == CUBIT_FORWARD)
    { 
      d1 = coedge->curve()->end_param();
      d1_ = coedge->curve()->start_param();
    }
    else
    {
      d1 = coedge->curve()->start_param();
      d1_ = coedge->curve()->end_param();
    }
    if( next_coedge->sense() == CUBIT_FORWARD)
    {
      d2_ = next_coedge->curve()->end_param();
      d2 = next_coedge->curve()->start_param();
    }
    else
    {
      d2 = next_coedge->curve()->end_param();
      d2_ = next_coedge->curve()->start_param();
    } 
    CubitVector v1, v1_, v2, v2_;
    coedge->curve()->position_from_u(d1, v1);
    coedge->curve()->position_from_u(d1_, v1_);
    next_coedge->curve()->position_from_u(d2, v2);
    next_coedge->curve()->position_from_u(d2_, v2_);
    if(v1.distance_between(v2) > v1_.distance_between(v2_))
    {
      //Reverse all coedges' senses.
      for (int i = 0; i < coedges_new.size(); i ++)
      {
        coedge = coedges_new.get_and_step();
        CubitSense sense = (coedge->sense() == CUBIT_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
        coedge->set_sense(sense);
      }
    }
  }

/*
  if(aShape.Orientation() == TopAbs_REVERSED)
    coedges_new.reverse();
  //Reverse all coedges' senses.
  for (int i = 0; i < coedges_new.size(); i ++)
  {
    OCCCoEdge* coedge = coedges_new.get_and_step();
    CubitSense sense = (coedge->sense() == CUBIT_FORWARD ? CUBIT_REVERSED : CUBIT_FORWARD);
    coedge->set_sense(sense);
  }
*/
  loop->coedges(coedges_new);

  //remove unused coedges
  for(int i = 0; i < coedges_old.size(); i++)
  {
    OCCCoEdge *coedge = coedges_old.get_and_step();
    if(coedge->get_mark() == 0)
      delete coedge;
    else
      coedge->set_mark(0);
  }
  return loop;
}

Curve* OCCQueryEngine::populate_topology_bridge(const TopoDS_Edge& aShape,
                                                CubitBoolean stand_along )
{
  if(aShape.IsNull())
    return (Curve*)NULL;
  Curve *curve;
  GProp_GProps myProps;
  BRepGProp::LinearProperties(aShape, myProps);
  double length =  myProps.Mass();
  TopExp_Explorer Ex;
  if(length < get_sme_resabs_tolerance())
  {
    //check if the two vertices are acctually the same point.
    CubitVector v[2];
    int i = 0;
    for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next())
    {
      TopoDS_Vertex vt = TopoDS::Vertex(Ex.Current());
      gp_Pnt pt = BRep_Tool::Pnt(vt);
      v[i] = CubitVector(pt.X(), pt.Y(), pt.Z());
      i++;
    }  
      
    if(v[0] == v[1])
      return (Curve*) NULL;
    else  
      PRINT_WARNING("Generated a sliver curve. \n");
  }

  if (!OCCMap->IsBound(aShape) ||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end())
  {
    TopoDS_Edge *poedge = new TopoDS_Edge;
    *poedge = aShape;
    curve = new OCCCurve(poedge);
    CurveList->append((OCCCurve*)curve);
    if(stand_along)
    {
      int k;
      if(!OCCMap->IsBound(aShape)) 
      {
        iTotalTBCreated++;
        k = iTotalTBCreated;
        OCCMap->Bind(*poedge, iTotalTBCreated);
      }
      else 
        k = OCCMap->Find(aShape);
      OccToCGM->insert(valType(k, (TopologyBridge*)curve));
    }
  }
  else 
  {
    int i = OCCMap->Find(aShape);
    curve = (Curve*)(OccToCGM->find(i))->second;
    CAST_TO(curve, OCCCurve)->set_TopoDS_Edge(aShape);
  }

  for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next())
  {
    TopoDS_Shape sh = Ex.Current();

   /* bool alreadyWrapped=false;
    if(OCCMap->IsBound(sh))
    {
      alreadyWrapped=true;
    }*/

    TBPoint* point = populate_topology_bridge(TopoDS::Vertex(Ex.Current()));
    CAST_TO(point, OCCPoint)->add_curve(CAST_TO(curve, OCCCurve));
    TopoDS_Shape parent = aShape;
    int current_id;
    add_shape_to_map(sh, parent, current_id);
    if(!OCCMap->IsBound(sh) ||
     OccToCGM->find(OCCMap->Find(sh)) == OccToCGM->end())
       OccToCGM->insert(valType(current_id, (TopologyBridge*)point));
     
  
   /* if(alreadyWrapped)
       PRINT_INFO("Vertex Already Wrapped:     "); 
    else
       PRINT_INFO("Vertex Wrapped:     ");  
    PRINT_INFO("     Shape ID = %d",  OCCMap->Find(sh) ); 
    PRINT_INFO("     TBPoint Address: %p", point );  
    PRINT_INFO("     Shape Orientation: %s\n", sh.Orientation()==TopAbs_FORWARD ? "Forward" : "Reversed" ); 
*/

  }

  return curve;
}

TBPoint* OCCQueryEngine::populate_topology_bridge(const TopoDS_Vertex& aShape,
                                                CubitBoolean stand_along)
{
  if(aShape.IsNull())
    return (TBPoint*)NULL;
  OCCPoint *point;
  if (iTotalTBCreated == 0 || !OCCMap->IsBound(aShape) ||
      OccToCGM->find(OCCMap->Find(aShape)) == OccToCGM->end()) 
  {
    TopoDS_Vertex *povertex = new TopoDS_Vertex;
    *povertex = aShape;
    point = new OCCPoint(povertex);
    if(stand_along)
    { 
      int k;
      if(!OCCMap->IsBound(aShape)) 
      {
        iTotalTBCreated++;
        k = iTotalTBCreated;
        OCCMap->Bind(*povertex, iTotalTBCreated);
      }
      else 
        k = OCCMap->Find(aShape);
      OccToCGM->insert(valType(k, (TopologyBridge*)point));
    }

  } 
  else 
  {
    int i = OCCMap->Find(aShape);
    point = (OCCPoint*)(OccToCGM->find(i))->second;
    point->set_TopoDS_Vertex(aShape);
  }
  return point;
}

TopologyBridge* OCCQueryEngine::occ_to_cgm(const TopoDS_Shape& shape)
{
  if(!OCCMap->IsBound(shape))
    return (TopologyBridge*) NULL;

  int k = OCCMap->Find(shape);
  return (OccToCGM->find(k))->second;
}	

//-----------------------------------------------------------------------
// Purpose       : Deletes all solid model entities associated with the
//                 Bodies in the input list.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 7/23/11
//-------------------------------------------------------------------------
void OCCQueryEngine::delete_solid_model_entities(DLIList<BodySM*>&bodyList)const
{
  delete_bodies(bodyList, true);
}
 
//-------------------------------------------------------------------------
// Purpose       : Deletes all solid model entities associated with the
//                 Bodies in the input list depending on remove_lower_entities
//                 flag.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 4/23/01
//-------------------------------------------------------------------------
void OCCQueryEngine::delete_bodies(DLIList<BodySM*>&bodyList,
                                   bool remove_lower_entities ) const
{
  BodySM* BodyPtr = NULL;
  for (int i = 0; i < bodyList.size(); i++ )
    {
      BodyPtr = bodyList.get_and_step();
      this->delete_body(BodyPtr, remove_lower_entities);
    }

  return;
}

CubitStatus OCCQueryEngine::delete_solid_model_entities(
    GeometryEntity* ref_entity_ptr,
    bool remove_lower_entities) const
{
     //Lump
   Lump* lump = CAST_TO(ref_entity_ptr, Lump);
   if(lump != NULL)
   {
     BodySM* body = CAST_TO(lump, OCCLump)->get_body();
     DLIList<Lump*> lumps = CAST_TO(body, OCCBody)->lumps();
 
     if (remove_lower_entities)
       return delete_solid_model_entities(body);

     DLIList<TopologyBridge*> children;
     for(int i = 0; i < lumps.size(); i++)
     {
       lump = lumps.get_and_step();
       CAST_TO(lump, OCCLump)->get_children_virt(children);
     }

     CubitStatus stat = this->unhook_BodySM_from_OCC(body); 
     if(stat)
     {
       while (children.size())
          delete children.pop();
       while(lumps.size())
          delete lumps.pop();
       delete body;
     }
     return stat;
   }

     // Surface
   Surface* ref_face_ptr = CAST_TO(ref_entity_ptr, Surface);
   if (ref_face_ptr != NULL)
   {
     if (remove_lower_entities)
       return ( this->delete_solid_model_entities(ref_face_ptr) );
     CubitStatus stat = this->unhook_Surface_from_OCC(ref_face_ptr);
     if(stat)
       delete ref_face_ptr;
     return stat;
   }

     // Curve
   Curve* ref_edge_ptr = CAST_TO(ref_entity_ptr, Curve);
   if (ref_edge_ptr != NULL)
   {
      if (remove_lower_entities)
        return ( this->delete_solid_model_entities(ref_edge_ptr));
      CubitStatus stat = this->unhook_Curve_from_OCC(ref_edge_ptr);
      if(stat)
        delete ref_edge_ptr;
      return stat;
   }

     // Point
   TBPoint* ref_vertex_ptr = CAST_TO(ref_entity_ptr, TBPoint);
   if (ref_vertex_ptr != NULL)
   {
      return ( this->delete_solid_model_entities(ref_vertex_ptr) );
   }

     // Oops!
   PRINT_ERROR("In AcisQueryEngine::delete_solid_model_entities\n"
               "       Can only delete solid model entities underlying \n"
               "RefFaces, RefEdges and RefVertices.\n");
   return CUBIT_FAILURE;

}
//-------------------------------------------------------------------------
// Purpose       : Delete a OCCBody and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( BodySM* bodysm) const
{
  return delete_body(bodysm, true);
}
//-------------------------------------------------------------------------
// Purpose       : Delete a OCCBody and child entities depending on
//                 remove_lower_entities flag.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_body( BodySM* bodysm,
                             bool remove_lower_entities) const
{
  OCCBody* occ_body = static_cast<OCCBody*>(bodysm);
  if (!occ_body)
    return CUBIT_FAILURE;

  DLIList<Lump*> lumps;
  DLIList<ShellSM*> shell_list;
  DLIList<OCCSurface*> surfaces = occ_body->my_sheet_surfaces();
  for(int i = 0;  i <surfaces.size(); i++)
  { 
    OCCSurface* occ_surface = surfaces.get_and_step();
    occ_surface->set_body((OCCBody*)NULL);
    if(remove_lower_entities)
    {
      delete occ_surface->my_lump();
      OCCShell* shell = occ_surface->my_shell();
      delete_solid_model_entities(occ_surface);
      delete shell;
    }
  }

  DLIList<OCCShell*> shells = occ_body->shells();
  for(int i = 0;  i <shells.size(); i++)
  {
    OCCShell* occ_shell = shells.get_and_step();
    occ_shell->set_body((OCCBody*)NULL);
    if(remove_lower_entities)
    {
      delete occ_shell->my_lump();
      DLIList<TopologyBridge*> tb_surfaces;
      occ_shell->get_children_virt(tb_surfaces);
      unhook_ShellSM_from_OCC(occ_shell);
      for(int k = 0; k < tb_surfaces.size(); k++)
        delete_solid_model_entities(CAST_TO(tb_surfaces.get_and_step(), Surface));
      delete occ_shell;
    }
  }

  DLIList<TopologyBridge*> children;
  lumps = occ_body->lumps();
  int size = lumps.size();

  for(int i =0; i < size; i++)
  {
    Lump* lump = lumps.get_and_step();
    OCCLump* occ_lump = CAST_TO(lump, OCCLump);
    if (!occ_lump)
       continue;
    occ_lump->remove_body();
    if(remove_lower_entities)
    {
      children.clean_out();
      occ_lump->get_children_virt(children);
      for(int j = 0; j < children.size(); j++)
      {
        ShellSM* shell = CAST_TO(children.get_and_step(), ShellSM);
      
        if (shell)
          shell_list.append(shell);
        DLIList<TopologyBridge*> tb_surfaces;
        shell->get_children_virt(tb_surfaces);
        for(int k = 0; k < tb_surfaces.size(); k++)
          delete_solid_model_entities(CAST_TO(tb_surfaces.get_and_step(), Surface));
      }
    }
  }
  CubitStatus stat = CUBIT_SUCCESS;
  stat = unhook_BodySM_from_OCC(bodysm, remove_lower_entities);

  if(remove_lower_entities)
  {
    for(int j = 0; j < shell_list.size(); j++)
       delete shell_list.get_and_step();

    for(int i =0; i < lumps.size(); i++)
       delete lumps.get_and_step(); 
  }

  BodyList->remove(occ_body);
  delete bodysm;
  return stat;
}

CubitStatus
OCCQueryEngine::unhook_BodySM_from_OCC( BodySM* bodysm ,
                                        bool remove_lower_entities)const
{
  OCCBody* occ_body = dynamic_cast<OCCBody*>(bodysm);
  if (!occ_body)
    return CUBIT_FAILURE;

  TopoDS_Shape* shape = occ_body->get_TopoDS_Shape();

  if (shape && !shape->IsNull())
  {
    //remove the entry from label tree
    OCCAttribSet::remove_attribute(*shape) ;

    //remove the entry from the map
    int k;
    if(shape && !shape->IsNull() && OCCMap->IsBound(*shape))
    {
        k = OCCMap->Find(*shape);
        OCCMap->UnBind(*shape);

        if(!OccToCGM->erase(k))
          PRINT_ERROR("The OccBody and iCreatedTotal %i pair is not in the map!", k);
    }
  }

  DLIList<Lump*> lumps = occ_body->lumps();
  for(int i =0; i < lumps.size()&& remove_lower_entities; i++)
  {
     Lump* lump = lumps.get_and_step();
     //OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     //if(occ_lump)
     //  occ_lump->remove_body();

     unhook_Lump_from_OCC(lump);
  }

  if (shape && !shape->IsNull())
    shape->Nullify();
  return CUBIT_SUCCESS;
} 
  
//-------------------------------------------------------------------------
// Purpose       : unhook a OCClumps and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 11/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Lump_from_OCC( Lump* lump ) const
{
  if (lump == NULL)
    return CUBIT_FAILURE;

  OCCLump* occ_lump = dynamic_cast<OCCLump*>(lump);
  if (!occ_lump)
    return CUBIT_FAILURE;

  TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();

  if(!solid)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*solid) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*solid))
    {
      k = OCCMap->Find(*solid);

      OCCMap->UnBind(*solid);

      if(!OccToCGM->erase(k))
	PRINT_ERROR("The OccLump and iCreatedTotal pair is not in the map!");
    }
  
  DLIList<TopologyBridge*> children;
  occ_lump->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     ShellSM* shell = CAST_TO(children.get_and_step(), ShellSM); 
     unhook_ShellSM_from_OCC(shell);
  }
  if (occ_lump->get_body() != NULL)
    BodyList->remove(CAST_TO(occ_lump->get_body(), OCCBody));

  if(!solid->IsNull())
     solid->Nullify();
  return CUBIT_SUCCESS;
} 

//-------------------------------------------------------------------------
// Purpose       : unhook a ShellSM from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_ShellSM_from_OCC( ShellSM* shell) const
{
  OCCShell* occ_shell = dynamic_cast<OCCShell*>(shell);
  if (!occ_shell)
    return CUBIT_FAILURE;

  TopoDS_Shell* Shell = occ_shell->get_TopoDS_Shell();

  if(!Shell)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*Shell))
    {
      k = OCCMap->Find(*Shell);

      OCCMap->UnBind(*Shell);

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccShell and iCreatedTotal pair is not in the map!");
    }

  if(!Shell->IsNull())
    Shell->Nullify();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCSurface and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Surface* surface)const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf)
    return CUBIT_FAILURE;

  OCCBody* sheet_body = fsurf->my_body();
  if(sheet_body != NULL)
    BodyList->remove(sheet_body);

  double d = fsurf->measure();
  if(d < 0.0)
    PRINT_WARNING("Negative area surface is being deleted. \n");

  DLIList<TopologyBridge*> children;
  fsurf->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     LoopSM* loop = CAST_TO(children.get_and_step(), LoopSM);
     delete_loop(loop);
  }
  CubitStatus stat = unhook_Surface_from_OCC(surface);
  if (stat)
    delete surface;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Surface from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Surface_from_OCC( Surface* surface) const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf)
    return CUBIT_FAILURE;

  TopoDS_Face *face = fsurf->get_TopoDS_Face();

  if(!face)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*face) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*face))
    {
      k = OCCMap->Find(*face);

      OCCMap->UnBind(*face);

      if(!OccToCGM->erase(k))
        PRINT_WARNING("The OccSurface and iCreatedTotal pair is not in the map!");
    }
  SurfaceList->remove(fsurf);
  if(!face->IsNull())
    face->Nullify();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCLoop and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_loop( LoopSM* loopsm)const
{
  OCCLoop* occ_loop = dynamic_cast<OCCLoop*>(loopsm);
  if (!occ_loop)
    return CUBIT_FAILURE;

  DLIList<OCCCoEdge*> children;
  children = occ_loop->coedges();
  int size = children.size();
  DLIList<Curve*> curves;
  while(size > 0)
  {
     OCCCoEdge* coedge = children.pop();
     Curve* curve = coedge->curve();
     curves.append_unique(curve);
     size = children.size();
  }
   
  CubitStatus status = unhook_LoopSM_from_OCC(loopsm);

  for(int i = 0; i < curves.size(); i ++)
  {
    OCCCurve* occ_curve = CAST_TO(curves.get_and_step(), OCCCurve);
    unhook_coedges_of_a_curve(occ_curve, occ_loop);
    occ_curve->remove_loop(occ_loop);
    delete_solid_model_entities(occ_curve);
  }

  if (status)
  {
    WireList->remove(occ_loop);
    delete loopsm;
  }
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a LoopSM from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_LoopSM_from_OCC( LoopSM* loopsm) const
{
  OCCLoop* occ_loop = dynamic_cast<OCCLoop*>(loopsm);
  if (!occ_loop)
    return CUBIT_FAILURE;

  TopoDS_Wire* wire = occ_loop->get_TopoDS_Wire();

  if(!wire)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*wire))
    {
      k = OCCMap->Find(*wire);

      OCCMap->UnBind(*wire);

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccLoop and iCreatedTotal pair is not in the map!");
    }

  if(!wire->IsNull())
    wire->Nullify(); 
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose      : unhook a list of OCCCoEdges from their underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_CoEdges_from_OCC( DLIList<OCCCoEdge*>& coedges) const
{
  int size = coedges.size();
  while(size > 0)
  {
     OCCCoEdge* coedge = coedges.pop();

     LoopSM* loopsm = coedge->loop();
     OCCLoop* loop = CAST_TO(loopsm, OCCLoop);
     assert(loop);
     loop->remove_coedge(coedge);
 
     delete coedge;

     size = coedges.size();
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCCurve and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Curve* curve)const
{
  OCCCurve* fcurve = CAST_TO(curve, OCCCurve);
  if (!fcurve )
    return CUBIT_FAILURE;

  if(fcurve->loops().size() > 0 && fcurve->loops().get() != NULL)
    return CUBIT_FAILURE;

  DLIList<TopologyBridge*> children;
  fcurve->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     TBPoint* point = CAST_TO(children.get_and_step(), TBPoint);
     CAST_TO(point, OCCPoint)->remove_curve(fcurve);
     delete_solid_model_entities(point);
  }
  
  CubitStatus stat = unhook_Curve_from_OCC(curve);
  if (!stat)
    return CUBIT_FAILURE;

  CurveList->remove(fcurve);
  delete fcurve;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Curve from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Curve_from_OCC( Curve* curve ) const
{
  OCCCurve* fcurve = dynamic_cast<OCCCurve*>(curve);
  if (!fcurve )
    return CUBIT_FAILURE;

  DLIList<TopologyBridge*> children;
  fcurve->get_children_virt(children);
  for(int i = 0; i < children.size(); i++)
  {
     TBPoint* point = CAST_TO(children.get_and_step(), TBPoint);
     CAST_TO(point, OCCPoint)->remove_curve(fcurve);
  }

  unhook_coedges_of_a_curve(fcurve, NULL);
    
  fcurve->clean_loops();
  TopoDS_Edge* edge = fcurve->get_TopoDS_Edge();
  if (!edge)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*edge) ;
  
  //remove the entry from the map
  int k;
  if(edge && !edge->IsNull() && OCCMap->IsBound(*edge))
    {
      k = OCCMap->Find(*edge);

      OCCMap->UnBind(*edge);

      if(!OccToCGM->erase(k))
        PRINT_WARNING("The OccCurve and iCreatedTotal pair is not in the map!");
    }
  CurveList->remove(fcurve); 
  if(!edge->IsNull())
    edge->Nullify();
  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------------
// Purpose       : Delete a OCCPoint and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( TBPoint* point) const
{
  OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);
  if (!fpoint)
    return CUBIT_FAILURE;

  if(fpoint->num_curves() > 0)
    return CUBIT_FAILURE;

  CubitStatus stat = unhook_Point_from_OCC(point);
  
  if(stat)
    delete point;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : unhook a Point from its underlining OCC entity.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::unhook_Point_from_OCC( TBPoint* point) const 
{
  OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);
  if (!fpoint)
    return CUBIT_FAILURE;

  TopoDS_Vertex* vertex = fpoint->get_TopoDS_Vertex();
  if (!vertex)
    return CUBIT_FAILURE;

  //remove the entry from label tree
  OCCAttribSet::remove_attribute(*vertex) ;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*vertex))
    {
      k = OCCMap->Find(*vertex);

      OCCMap->UnBind(*vertex);

      if(!OccToCGM->erase(k))
        PRINT_ERROR("The OccPoint and iCreatedTotal pair is not in the map!");
    }
  if(!vertex->IsNull())
    vertex->Nullify();
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : fire a ray at the specified body, returning the entities hit.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 12/12/07
//-------------------------------------------------------------------------
CubitStatus OCCQueryEngine::fire_ray(  CubitVector &origin,
                        CubitVector &direction,
                        DLIList<TopologyBridge*> &at_entity_list,
                        DLIList<double> &ray_params,
                        int max_hits ,
                        double ray_radius ,
                        DLIList<TopologyBridge*> *hit_entity_list_ptr )const
{
  CubitStatus status = CUBIT_SUCCESS;

  //- fire a ray at the specified body, returning the entities hit and
  //- the parameters along the ray; return CUBIT_FAILURE if error
  // - line body intersection. 
  gp_Pnt p(origin.x(), origin.y(), origin.z());
  gp_Dir dir(direction.x(), direction.y(), direction.z());
  gp_Lin L(p, dir);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(L); 
  
  for(int i = 0; i < at_entity_list.size(); i++)
  {
    TopologyBridge* tb = at_entity_list.get_and_step();
    TopoDS_Shape *shape;
    OCCBody *occBody = CAST_TO(tb, OCCBody);
    if (occBody )
    {
      occBody->get_TopoDS_Shape(shape);
    
      BRepExtrema_DistShapeShape distShapeShape(edge, *shape);
      //distShapeShape.Perform();
      if (!distShapeShape.IsDone())
      {
        PRINT_ERROR("Cannot calculate the intersection points for the input body.\n");
        return CUBIT_FAILURE;
      }

      if (distShapeShape.Value() < get_sme_resabs_tolerance())
      {
        int numPnt = distShapeShape.NbSolution();
        for (int i = 1; i <= numPnt; i++)
        {
	  double para;
	  distShapeShape.ParOnEdgeS1(i , para);
	  ray_params.append(para);
        }
        hit_entity_list_ptr->append(tb);
      } 
    }
  }
  return status;
}

double OCCQueryEngine::get_sme_resabs_tolerance() const
{
  return Precision::Confusion(); 
}
// Gets solid modeler's resolution absolute tolerance

double OCCQueryEngine::set_sme_resabs_tolerance( double p)
{
  double old_p = get_sme_resabs_tolerance();
  BRepBuilderAPI::Precision(p);
  return old_p;
}

CubitStatus OCCQueryEngine::set_int_option( const char* , int )
{
  PRINT_ERROR("OCCQueryEngine::set_int_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_dbl_option( const char* , double )
{
  PRINT_ERROR("OCCQueryEngine::set_dbl_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_str_option( const char* , const char* )
{
  PRINT_ERROR("OCCQueryEngine::set_str_option not yet implemented.\n");
  return CUBIT_FAILURE;
}
//- Set solid modeler options


//===========================================================================
//Function Name: ensure_is_ascii_stl_file
//Member Type:
//Description: returns CUBIT_TRUE in is_ascii if fp is an ascii stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================
CubitStatus OCCQueryEngine::ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii)
{

  char line[128]="";

  if (fgets(line, 128, fp)==NULL)
    {
      return CUBIT_FAILURE;
    }
  if (fgets(line, 128, fp)==NULL)
    {
      return CUBIT_FAILURE;
    }
  if (strlen(line)==127)
    {
      if (fgets(line, 128, fp)==NULL)
	{
	  return CUBIT_FAILURE;
	}
    }


  unsigned int dummy_int=0;

  while (isspace(line[dummy_int])&& dummy_int<strlen(line)) dummy_int++;

  if (strlen(line)-dummy_int>5)
    {
      if (tolower(line[dummy_int++])=='f' &&
	  tolower(line[dummy_int++])=='a' &&
	  tolower(line[dummy_int++])=='c' &&
	  tolower(line[dummy_int++])=='e' &&
	  tolower(line[dummy_int])=='t')
	{
	  if (fgets(line, 128, fp)==NULL)
	    {
	      return CUBIT_FAILURE;
	    }
	  dummy_int=0;
	  while (isspace(line[dummy_int])&& dummy_int<strlen(line))
	    {
	      dummy_int++;
	    }
	  if (strlen(line)-dummy_int>5)
	    {
	      if (tolower(line[dummy_int++])=='o' &&
		  tolower(line[dummy_int++])=='u' &&
		  tolower(line[dummy_int++])=='t' &&
		  tolower(line[dummy_int++])=='e' &&
		  tolower(line[dummy_int])=='r')
		{
		  if (fgets(line, 128, fp)==NULL)
		    {
		      return CUBIT_FAILURE;
		    }
		  dummy_int=0;
		  while (isspace(line[dummy_int])&& dummy_int<strlen(line)) {
		    dummy_int++;
		  }
		  if (strlen(line)-dummy_int>6)
		    {
		      if (tolower(line[dummy_int++])=='v' &&
			  tolower(line[dummy_int++])=='e' &&
			  tolower(line[dummy_int++])=='r' &&
			  tolower(line[dummy_int++])=='t' &&
			  tolower(line[dummy_int++])=='e'	&&
			  tolower(line[dummy_int])=='x')
			{
			  is_ascii=CUBIT_TRUE;
			}
		    }
		}
	    }
	}
    }
  return CUBIT_SUCCESS;
}


//=============================================================================
//Function:   create_super_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a list of BodySMs
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus OCCQueryEngine::create_super_bounding_box(
						      DLIList<BodySM*>& body_list,
						      CubitBox& super_box )
{
  BodySM *bodySM;
  int i;
  CubitStatus status = CUBIT_SUCCESS;

  body_list.reset();
  for ( i = 0; i < body_list.size(); i++ ) {
    bodySM = body_list.get_and_step();  
    OCCBody* occBody = CAST_TO(bodySM, OCCBody);
    super_box |= occBody->get_bounding_box();
  }

  return status;
}

CubitStatus OCCQueryEngine::restore_transform( BodySM* body )
{
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  OCCBody* theBody = dynamic_cast<OCCBody*>(body);
  return theBody ? theBody->move( d.x(), d.y(), d.z() ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::rotate( BodySM* body, const CubitVector& v, double a )
{
  // a is in degree.
  OCCBody* occ_bod = dynamic_cast<OCCBody*>(body);
  a *= CUBIT_PI/180;
  return occ_bod ? occ_bod->rotate( v.x(), v.y(), v.z(), a ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, double factor )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( factor ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, const CubitVector& f )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( f.x(), f.y(), f.z() ) : CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Solid, Surface, Curve, or Vertex
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/23/07
//-------------------------------------------------------------------------
CubitStatus OCCQueryEngine::translate( GeometryEntity* entity,
                                       const CubitVector& v )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Vec aVec(v.x(), v.y(),v.z());
  gp_Trsf aTrsf;  
  aTrsf.SetTranslation(aVec);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::update_entity_shape(GeometryEntity* entity_ptr,
					BRepBuilderAPI_ModifyShape* aBRepTrsf,
                                        BRepAlgoAPI_BooleanOperation *op)
{
  if (OCCBody *body_ptr = CAST_TO( entity_ptr, OCCBody))
    {
      body_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  else if( OCCSurface *surface_ptr = CAST_TO( entity_ptr, OCCSurface))
    {
      surface_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  else if( OCCCurve *curve_ptr = CAST_TO( entity_ptr, OCCCurve))
    {
       curve_ptr->update_OCC_entity(aBRepTrsf, op); 
       return CUBIT_SUCCESS;
    }

  else if( OCCPoint *point_ptr = CAST_TO( entity_ptr, OCCPoint))
    {
      point_ptr->update_OCC_entity(aBRepTrsf, op);
      return CUBIT_SUCCESS;
    }

  PRINT_ERROR("Non-OCC TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return CUBIT_FAILURE;
}

//a is angular value of rotation in radians
CubitStatus OCCQueryEngine::rotate( GeometryEntity* entity,
                                    const CubitVector& v, double a )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax1 anAxis(aOrigin, aDir);

  //a is angular value of rotation in radians 
  gp_Trsf aTrsf;
  aTrsf.SetRotation(anAxis, a);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* entity, double f )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_Trsf aTrsf;
  aTrsf.SetScaleFactor(f);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* entity, 
                                   const CubitVector& v )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
                   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

  gp_GTrsf gTrsf;
  gTrsf.SetValue(1,1, v.x());
  gTrsf.SetValue(2,2, v.y());
  gTrsf.SetValue(3,3, v.z());

  BRepBuilderAPI_GTransform gBRepTrsf(gTrsf);
  gBRepTrsf.Perform(*shape);
  update_entity_shape(entity, &gBRepTrsf);
  return CUBIT_SUCCESS;
}

// like ACIS, here v is the normal of symmetric plane.
CubitStatus OCCQueryEngine::reflect( GeometryEntity* entity, 
                                     const CubitVector&  v)
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
    {
      PRINT_ERROR( "problem occured getting OCC entity.\n"
		   "       Aborting.\n" );
      return CUBIT_FAILURE;
    }
 
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax2 anAx2(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAx2);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  update_entity_shape(entity, &aBRepTrsf);
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : bodies_overlap
// Member Type: PUBLIC
// Description: determine if OCC-based bodies overlap
// Author     : Jane Hu 
// Date       : 10/07
//===============================================================================
CubitBoolean OCCQueryEngine::bodies_overlap (BodySM * body_ptr_1,
                                             BodySM * body_ptr_2 ) const
{
  OCCBody *occ_body1 = CAST_TO(body_ptr_1, OCCBody);
  if (!occ_body1)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
      return CUBIT_FALSE;
    }
 
  OCCBody *occ_body2 = CAST_TO(body_ptr_2, OCCBody);
  if (!occ_body2)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
      return CUBIT_FALSE;
    }

  CubitBox box_1 = occ_body1->get_bounding_box();
  CubitBox box_2 = occ_body2->get_bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1;
  occ_body1->get_TopoDS_Shape(shape1);
  TopoDS_Shape *shape2;
  occ_body2->get_TopoDS_Shape(shape2); 
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_SOLID); Ex1.More(); Ex1.Next())
    {
      TopoDS_Solid *posolid1 =  new TopoDS_Solid;
      *posolid1 = TopoDS::Solid(Ex1.Current());
      OCCLump * lump1 = new OCCLump(posolid1); 
      for (Ex2.Init(*shape2, TopAbs_SOLID); Ex2.More(); Ex2.Next())
	{
	  TopoDS_Solid *posolid2 =  new TopoDS_Solid;
	  *posolid2 = TopoDS::Solid(Ex2.Current());
	  OCCLump * lump2 = new OCCLump(posolid2);
	  CubitBoolean is_overlap = volumes_overlap(lump1, lump2);
	  if(is_overlap)
	    return CUBIT_TRUE;
	}
    }
  return CUBIT_FALSE;
}

CubitBoolean OCCQueryEngine::volumes_overlap (Lump *lump1, Lump *lump2 ) const
{
  OCCLump *occ_lump1 = CAST_TO(lump1, OCCLump);
  if (!occ_lump1)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
      return CUBIT_FALSE;
    }

  OCCLump *occ_lump2 = CAST_TO(lump2, OCCLump);
  if (!occ_lump2)
    {
      PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
      return CUBIT_FALSE;
    }

  CubitBox box_1 = occ_lump1->bounding_box();
  CubitBox box_2 = occ_lump2->bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1 = (TopoDS_Shape*)(occ_lump1->get_TopoDS_Solid());
  TopoDS_Shape *shape2 = (TopoDS_Shape*)(occ_lump2->get_TopoDS_Solid());
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_FACE); Ex1.More(); Ex1.Next())  
    {
      for (Ex2.Init(*shape2, TopAbs_FACE); Ex2.More(); Ex2.Next()) 
	{
	  BRepAlgoAPI_Section section(Ex1.Current(), Ex2.Current());
	  if (section.HasGenerated())
	    return CUBIT_TRUE;
	}
    }
  return CUBIT_FALSE;
}

void OCCQueryEngine::copy_attributes(TopoDS_Shape& old_shape,
                                     TopoDS_Shape& new_shape)
{
  if(new_shape.IsNull())
    return;

  //update the attribute label tree
  DLIList<CubitSimpleAttrib> list;
  OCCAttribSet::get_attributes(old_shape, list);

  for(int i = 0; i < list.size(); i ++)
  {
    const CubitSimpleAttrib& s_attr = list.get_and_step();
    TopAbs_ShapeEnum type = old_shape.ShapeType();
    if(new_shape.ShapeType() < type)
    {
      TopTools_IndexedMapOfShape M;
      TopExp::MapShapes(new_shape, type, M); 
      for(int j = 1; j <= M.Extent() ; j++ )
      {
        TopoDS_Shape sub_shape = M(j);
        OCCAttribSet::append_attribute(s_attr, sub_shape); 
      }
    }
    else    
      OCCAttribSet::append_attribute(s_attr, new_shape);
  }
}

int OCCQueryEngine::update_OCC_map(TopoDS_Shape& old_shape, 
                                   TopoDS_Shape& new_shape)
{
  if (old_shape.IsNull() || !OCCMap->IsBound(old_shape) || 
      old_shape.IsEqual(new_shape))
    return -1;

  //update the attribute label tree
  int current_id = OCCMap->Find(old_shape);
  std::map<int, TDF_Label>::iterator it_lab =
          Shape_Label_Map->find(current_id);
  CubitBoolean newlyBound = CUBIT_FALSE;
  TopTools_IndexedMapOfShape M;
  TopoDS_Shape new_subshape;
  new_subshape.Nullify();

  if(old_shape.ShapeType() == TopAbs_SOLID)
    TopExp::MapShapes(new_shape, TopAbs_SOLID, M);
  else if (old_shape.ShapeType() == TopAbs_SHELL)
    TopExp::MapShapes(new_shape, TopAbs_SHELL, M);
  else if (old_shape.ShapeType() == TopAbs_FACE)
    TopExp::MapShapes(new_shape, TopAbs_FACE, M);

  if(it_lab != Shape_Label_Map->end())
  {
     CubitBoolean isNewShapeBound = CUBIT_FALSE;
     if(old_shape.ShapeType() > TopAbs_COMPOUND && !new_shape.IsNull() &&
        new_shape.ShapeType() == TopAbs_COMPOUND && M.Extent() == 1)
     {
       new_subshape = M(1); 
       isNewShapeBound = OCCMap->IsBound(new_subshape);
       copy_attributes(old_shape, new_subshape);
       Shape_Label_Map->erase(current_id);
       if(isNewShapeBound != OCCMap->IsBound(new_subshape)) 
         newlyBound = CUBIT_TRUE;
     }
     else
     {
       isNewShapeBound = OCCMap->IsBound(new_shape);
       copy_attributes(old_shape, new_shape);
       Shape_Label_Map->erase(current_id);
       if(isNewShapeBound != OCCMap->IsBound(new_shape))
         newlyBound = CUBIT_TRUE; 
     }
  }

  //update CGM-OCC map
  int k = current_id;
  assert (k > 0 && k <= iTotalTBCreated);

  std::map<int, TopologyBridge*>::iterator it = OccToCGM->find(k);
  TopologyBridge* tb = NULL;
  if (it != OccToCGM->end()) 
     tb = (*it).second;

  //unless just changing location, if the TShape is going to change, remove
  //old curve_list .
  CubitBoolean curve_removed = CUBIT_FALSE;
  if(old_shape.ShapeType() == TopAbs_VERTEX && !new_shape.IsPartner(old_shape)) 
  {
    //remove the curve list associated with the vertex too.
    GeometryEntity* ge =  CAST_TO(tb, GeometryEntity);
    if (ge)
    {
      OCCPoint* test_p = CAST_TO(ge, OCCPoint);
      if(test_p)
      {
        test_p->clear_curves();
        curve_removed = CUBIT_TRUE;
      }
    }
  }

  OCCMap->UnBind(old_shape);
  if(new_subshape.IsNull())
    new_subshape = new_shape;

  if (tb && TopAbs_SOLID == old_shape.ShapeType() && !new_shape.IsNull() && 
       TopAbs_COMPOUND == new_shape.ShapeType() && M.Extent() > 1)
  {
    OccToCGM->erase(k);
    GeometryEntity* ge =  CAST_TO(tb, GeometryEntity);
    if(ge)
      delete_solid_model_entities( ge, CUBIT_TRUE);
    return k;
  }

  else if(tb && TopAbs_FACE == old_shape.ShapeType() && !new_shape.IsNull() &&
       TopAbs_COMPOUND == new_shape.ShapeType() && M.Extent() > 1)
  {
    GeometryEntity* ge =  CAST_TO(tb, GeometryEntity);
    if(ge)
    {
      OCCSurface *face = CAST_TO(ge, OCCSurface);
      OCCShell* shell = face->my_shell();
      if (shell)
      {
        TopoDS_Shell* shape = shell->get_TopoDS_Shell();
        assert(!shape );
        if (shape) {
          PRINT_ERROR("Unexpected non-NULL TopoDS_Shell pointer.\n");
          return -1;
        }
      }
      OCCLump* lump = face->my_lump();
      if(lump)
      {
        TopoDS_Solid* shape = lump->get_TopoDS_Solid();
        assert(!shape );
        if (shape) {
          PRINT_ERROR("Unexpected non-NULL TopoDS_Solid pointer.\n");
          return -1;
        }
        delete lump;
      }
      OCCBody* body = face->my_body();
      if (body)
        delete body;
      
      delete_solid_model_entities(ge, CUBIT_TRUE); 
      if(shell)
        delete shell;
    }
    return k;
  }

  else if (tb && ((!new_subshape.IsNull() && !old_shape.IsSame(new_subshape)&&
        OCCMap->IsBound(new_subshape) && 
        OccToCGM->find(OCCMap->Find(new_subshape))!= OccToCGM->end()) || 
        new_subshape.IsNull()))
  //already has a TB built on new_shape
  {
    //delete the second TB corresponding to old_shape
    OccToCGM->erase(k);
    GeometryEntity* ge =  CAST_TO(tb, GeometryEntity);
    if(ge)
    {
      //PRINT_INFO("TB: %p\n",ge);
      Lump* lump = CAST_TO(ge, Lump);
      if(lump)
      {
        BodySM* body = CAST_TO(lump,OCCLump)->get_body();
        if(body)
        {
          //OCCBody* occ_body = CAST_TO(body, OCCBody);
          //TopoDS_Compound* pshape = occ_body->get_TopoDS_Shape(); 
          //if(!pshape || pshape->IsNull())
          delete_solid_model_entities(body);
        }
      }

      else
      {
        OCCPoint* test_p = CAST_TO(ge, OCCPoint);
        if(test_p)
          test_p->clear_curves();
        delete_solid_model_entities( ge, CUBIT_FALSE );
      }
    }

    else
    {
      ShellSM * shell = CAST_TO(tb, ShellSM);
      if(shell)
      {
        DLIList<OCCSurface*> memberSurfaces;
        OCCShell* occShell = CAST_TO(shell, OCCShell);
        memberSurfaces = occShell->getMemberSurfaces();
        while (memberSurfaces.size())
        {
          OCCSurface* memberSurf = memberSurfaces.pop();
          if (SurfaceList->is_in_list(memberSurf))
            // if the surface has not been unhooked and deleted
            memberSurf->set_shell((OCCShell*)NULL);
        }

        OCCLump* lump = CAST_TO(shell, OCCShell)->my_lump();
        if(lump && (lump->get_TopoDS_Solid() == NULL || 
           !OCCMap->IsBound(*(lump->get_TopoDS_Solid()))))
        {
          delete CAST_TO(shell, OCCShell)->my_body();
          delete lump;
        }
        unhook_ShellSM_from_OCC(shell);
        delete shell;
        return k;
      }
      LoopSM* loop = CAST_TO(tb, LoopSM);
      if(loop)
      {
         DLIList<OCCCoEdge*> children;
         children = CAST_TO(loop, OCCLoop)->coedges();
         while(children.size())
         {
           OCCCoEdge* coedge = children.pop();
           CAST_TO(coedge->curve(), OCCCurve)->remove_loop(CAST_TO(loop, OCCLoop));
           delete (OCCCoEdge*)coedge;
         }
         unhook_LoopSM_from_OCC(loop);
         delete loop;
         return k;
      }
    }
  }

  else
  {   
    //if the new_shape is bounded in copy_attribute, unbind it and rebind to 
    //the old k
    if(newlyBound)
    {
      int new_k = OCCMap->Find(new_subshape);
      TDF_Label aLabel;
      CubitBoolean found = CUBIT_FALSE;
      OCCAttribSet::FindShape(new_subshape, aLabel, found);
      assert(found);
      Shape_Label_Map->erase(new_k);
      Shape_Label_Map->insert(labType(k, aLabel)); 
      OCCMap->UnBind(new_subshape);
    }      
    if(!OCCMap->IsBound(new_subshape))
      OCCMap->Bind(new_subshape, k);
    if(tb && !curve_removed)
      set_TopoDS_Shape(tb, new_subshape);
  }
  return k;
}

void OCCQueryEngine::unhook_coedges_of_a_curve(OCCCurve* curve,
                                               OCCLoop* loop)const 
{
  DLIList<OCCCoEdge*> coedges;
  DLIList<OCCCoEdge *> children ;
  if (loop != NULL)
    children = loop->coedges();
  else
  {
    DLIList<OCCLoop*> loops;
    loops = curve->loops();
    for (int i = 0; i < loops.size(); i ++)
      children += loops.get_and_step()->coedges();
  } 
  for(int j = 0; j < children.size(); j++)
  {
     OCCCoEdge* coedge = children.get_and_step();
     if (coedge->curve() == curve)
       coedges.append(coedge);
  }
  unhook_CoEdges_from_OCC(coedges);
}

void OCCQueryEngine::set_TopoDS_Shape(TopologyBridge* tb,
                                      TopoDS_Shape& shape)
{
  BodySM* body = CAST_TO(tb, BodySM);
  if(body)
    return CAST_TO(body, OCCBody)->set_TopoDS_Shape(TopoDS::Compound(shape));

  Lump* lump = CAST_TO(tb, Lump);
  if(lump)
    return CAST_TO(lump, OCCLump)->set_TopoDS_Solid(TopoDS::Solid(shape));

  ShellSM* shell = CAST_TO(tb, ShellSM);
  if(shell)
    return CAST_TO(shell, OCCShell)->set_TopoDS_Shell(TopoDS::Shell(shape));

  Surface* surface = CAST_TO(tb, Surface);
  if (surface)
    return CAST_TO(surface, OCCSurface)->set_TopoDS_Face(TopoDS::Face(shape));

  LoopSM* loop =  CAST_TO(tb, LoopSM);
  if(loop)
    return CAST_TO(loop, OCCLoop)->set_TopoDS_Wire(TopoDS::Wire(shape));

  Curve* curve = CAST_TO(tb, Curve);
  if (curve)
    return CAST_TO(curve, OCCCurve)->set_TopoDS_Edge(TopoDS::Edge(shape));

  TBPoint* point = CAST_TO(tb, TBPoint);
  if(point)
    return CAST_TO(point, OCCPoint)->set_TopoDS_Vertex(TopoDS::Vertex(shape));
  
}

//Added by Jane Hu on 06/17/11
void OCCQueryEngine::bound_TopoDS_Shape(const TopoDS_Shape & aShape)
{
  (iTotalTBCreated)++;
  switch (aShape.ShapeType())
  {
    case TopAbs_COMPOUND:
      OCCMap->Bind(TopoDS::Compound(aShape), iTotalTBCreated); 
      break;
    case TopAbs_SOLID:
      OCCMap->Bind(TopoDS::Solid(aShape), iTotalTBCreated);
      break;
    case TopAbs_SHELL:
      OCCMap->Bind(TopoDS::Shell(aShape), iTotalTBCreated);
      break;
    case TopAbs_FACE: 
      OCCMap->Bind(TopoDS::Face(aShape), iTotalTBCreated);
      break;
    case TopAbs_WIRE:
      OCCMap->Bind(TopoDS::Wire(aShape), iTotalTBCreated);
      break;
    case TopAbs_EDGE:
      OCCMap->Bind(TopoDS::Edge(aShape), iTotalTBCreated);
      break;
    case TopAbs_VERTEX:
      OCCMap->Bind(TopoDS::Vertex(aShape), iTotalTBCreated);
      break;

    default:
      (iTotalTBCreated)--;
  }
}
//EOF
