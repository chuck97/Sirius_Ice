/**
 * \file makept.cpp
 *
 * \brief makept, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * and performs varies checks for bodies, surfaces, curves and vertices.
 */

#undef NDEBUG
#include <cassert>
#include <Standard_Version.hxx>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CastTo.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "OCCSurface.hpp"
#include "OCCCurve.hpp"
#include "OCCShell.hpp"
#include "TopoDS_Shape.hxx"
#include "InitCGMA.cpp"
#include "CubitCompat.hpp"

#ifndef SRCDIR
# define SRCDIR .
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#define SRCPATH STRINGIFY(SRCDIR) "/"

// forward declare some functions used and defined later
CubitStatus read_geometry(int, const char **, bool local = false);
CubitStatus make_Point();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("=======================================\n");


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus status = InitCGMA::initialize_cgma("OCC");
  if (CUBIT_SUCCESS != status) return 1;

  //Do make point.
  status = make_Point();
  if (status == CUBIT_FAILURE) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val != 7 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  else
    ret_val = 0;

  return ret_val;
  
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
CubitStatus read_geometry(int num_files, const char **argv, bool local)
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
  PRINT_SEPARATOR;

  for (i = 0; i < num_files; i++) {
    std::string filename( local ? "./" : SRCPATH );
    filename += argv[i];
    status = CubitCompat_import_solid_model(filename.c_str(), "OCC");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", filename.c_str());
      abort();
    }
  }
  PRINT_SEPARATOR;

  return CUBIT_SUCCESS;
}

CubitStatus make_Point()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  OCCModifyEngine* ome = OCCModifyEngine::instance();

  //test for tweak fillet and chamfer
  Body* body = gmti->brick(10, 10, 10);
  //Added for code development testing, can be removed later.
  DLIList<RefFace*> test_ref_faces;
  DLIList<RefFace*> test_ref_faces2;
  body->ref_faces(test_ref_faces);
  test_ref_faces2.append(test_ref_faces.pop());
  CubitVector sweep_direction(0,0,6);
  DLIList<Body*> all_bodies;
  CubitStatus sta = gmti->tweak_move(test_ref_faces2, sweep_direction, all_bodies);
  if (CUBIT_SUCCESS != sta) {
    PRINT_ERROR("GeometryModifyTool::tweak_move failed.\n");
    return sta;
  }
  assert(all_bodies.size() == 1);
  double volume = all_bodies.get()->measure();
  assert(fabs(volume-1600) <= 0.00001);
  int n = all_bodies.get()->num_ref_faces();
  assert(n == 10);

  test_ref_faces.clean_out();
  test_ref_faces2.clean_out();
  all_bodies.get()->ref_faces(test_ref_faces);
  for(int i = 0 ; i < test_ref_faces.size(); i++)
  {
    double max = test_ref_faces.get()->bounding_box().max_z();
    double min = test_ref_faces.get()->bounding_box().min_z();
    if(max-11. < 0.000001 && max-11. > -0.000001 &&
       min-11. < 0.000001 && min-11. > -0.000001)
      test_ref_faces2.append(test_ref_faces.get());
    test_ref_faces.step();
  }
  sweep_direction.set(0,0,-6);
  all_bodies.clean_out();
  sta = gmti->tweak_move(test_ref_faces2, sweep_direction, all_bodies);
  assert(all_bodies.size() == 1);
  volume = all_bodies.get()->measure();
  assert(fabs(volume-1000) <= 0.00001);
  n = all_bodies.get()->num_ref_faces();
  assert(n == 6);
  body = all_bodies.get();
  
  Body* body2 = gmti->brick(10, 10, 10);
  CubitVector scale(1, 2, 5);
  CubitVector p1(0, 0, 0);
  gmti->scale(body2, p1, scale);
  assert(fabs(body2->measure()-10000)< 0.0001);
  CubitBox box1;
  box1 = body2->bounding_box();
  assert(fabs(box1.max_x()-5)<0.001 && fabs(box1.min_x() + 5) < 0.001);
  assert(fabs(box1.max_y()-10)<0.001 && fabs(box1.min_y() + 10) < 0.001);
  assert(fabs(box1.max_z()-25)<0.001 && fabs(box1.min_z() + 25) < 0.001);
  scale.z(0.2);
  scale.y(0.5);
  gmti->scale(body2, p1, scale);
  assert(fabs(body2->measure()-1000)< 0.0001);

  CubitVector away(11, 0 , 0);
  DLIList<Body*> bods;
  bods.append(body2);
  gti->translate(bods,away);
  DLIList<RefEdge*> ref_edges;
  body->ref_edges(ref_edges);
  DLIList<Body*> new_bodies;
  int isize=ref_edges.size();
  //isize = 12
  for(int i = 2; i < isize; i++)
    ref_edges.pop();

  DLIList<RefVertex*> ref_vertices, ref_vertices2;
  body->ref_vertices(ref_vertices);
  RefVertex* vtx1 ;
  ref_vertices2.append(ref_vertices.pop());
  ref_vertices2.append(ref_vertices.pop());

  CubitStatus status = gmti->tweak_fillet(ref_vertices2, 1, new_bodies);
  assert(status == CUBIT_FAILURE);

  //make a sheet surfaces to do fillet on vertices
  DLIList<RefFace*> ref_faces;
  body->ref_faces(ref_faces);
  RefFace* face = gmti->make_RefFace(ref_faces.get());
  RefFace* face2 = gmti->make_RefFace(ref_faces.get());
  RefFace* tweak_face = gmti->make_RefFace(ref_faces.get());
  DLIList<RefEdge*> tweak_edges;
  tweak_face->ref_edges(tweak_edges);
  tweak_edges.pop();
  tweak_edges.pop();
  tweak_edges.pop();
  Body* tweak_body = tweak_face->body();
  away *= 3;
  DLIList<Body*> tweak_bods;
  bods.append(tweak_body);
  gti->translate(tweak_bods, away);
  away /= 3;
  CubitVector delta(-1,0, 0);
  DLIList<Body*> tweak_bodies;
  gmti->tweak_move( tweak_edges, delta, tweak_bodies, CUBIT_FALSE,
                    CUBIT_FALSE);
  assert(tweak_bodies.size() == 1);
  CubitBoolean is_sheet = tweak_bodies.get()->is_sheet_body();
  assert(is_sheet == CUBIT_FALSE);
  DLIList<OCCSurface*> _surfaces;
  OCCBody *occ_body;
  occ_body = CAST_TO(tweak_bodies.get()->get_body_sm_ptr(), OCCBody);
  occ_body->get_all_surfaces(_surfaces);
  double area = _surfaces.get_and_step()->measure() + _surfaces.get()->measure();
  assert(fabs(area - 110) <= 0.00001);

  Body* sheet_body = face->body();
  Body* sheet_body2 = face2->body();
  away *= 5;
  DLIList<Body*> sheet_bods;
  sheet_bods.append(sheet_body);
  gti->translate(sheet_bods, away);
  away /= 2;
  DLIList<Body*> sheet_bods2;
  sheet_bods2.append(sheet_body2);
  gti->translate(sheet_bods2, away);

  ref_vertices2.clean_out();
  ref_vertices.clean_out();
  face->ref_vertices(ref_vertices2);
  isize=ref_vertices2.size();
  for(int i = 2; i < isize; i++)
    ref_vertices.append(ref_vertices2.pop());
  
  status = gmti->tweak_fillet(ref_vertices2, 1, new_bodies);
  DLIList<RefEdge*> ref_edges_check;
  new_bodies.get()->ref_edges(ref_edges_check);
  isize=ref_edges_check.size();
  assert(isize == 6);

  status = gmti->tweak_chamfer(ref_vertices, 1, new_bodies);
  ref_edges_check.clean_out();
  new_bodies.get()->ref_edges(ref_edges_check);
  isize=ref_edges_check.size();
  assert(isize == 8);

  ref_vertices.clean_out();
  face2->ref_vertices(ref_vertices);
  ref_vertices.pop();
  ref_vertices.pop();
  ref_vertices.pop();
  ref_edges_check.clean_out();
  ref_vertices.get()->ref_edges(ref_edges_check);
  new_bodies.clean_out();
  status = gmti->tweak_chamfer(ref_vertices, 1, new_bodies,ref_edges_check.pop(), 0.5, ref_edges_check.pop());
 
  ref_edges_check.clean_out();
  new_bodies.get()->ref_edges(ref_edges_check);
  isize=ref_edges_check.size();
  assert(isize == 5);

  new_bodies.clean_out();
  gmti->tweak_fillet(ref_edges, 1, new_bodies, CUBIT_FALSE, CUBIT_FALSE);

  ref_edges.clean_out();
  new_bodies.get()->ref_edges(ref_edges);
  isize=ref_edges.size();
  assert(isize == 21);

  ref_edges.clean_out();
  new_bodies.clean_out(); 

  body2->ref_edges(ref_edges);
  isize=ref_edges.size();
  assert(isize == 12);
  for(int i = 2; i < isize; i++)
    ref_edges.pop(); 
  
  gmti->tweak_chamfer(ref_edges, 1, new_bodies, 1);

  ref_edges.clean_out();
  new_bodies.get()->ref_edges(ref_edges);
  isize=ref_edges.size();
  assert(isize == 17);

  ref_vertices.clean_out();
  ref_vertices2.clean_out();
  new_bodies.get()->ref_vertices(ref_vertices);
  
  for(int i = 0; i < 3; i++)
  {
    vtx1 = ref_vertices.pop();
    ref_vertices2.append(vtx1);
  }

  new_bodies.clean_out();
  ref_edges.clean_out();
  ref_vertices2.get()->ref_edges(ref_edges);
  status = gmti->tweak_chamfer(ref_vertices2, 1, new_bodies, ref_edges.pop(), 
                         0.7, ref_edges.pop(), 0.5, ref_edges.pop());
  assert(status == CUBIT_FAILURE);

  ref_vertices.clean_out();
  ref_vertices2.pop();
  ref_vertices2.pop();
  ref_vertices2.get()->ref_edges(ref_edges);
  gmti->tweak_chamfer(ref_vertices2, 2, new_bodies, ref_edges.pop(), 1,
                      ref_edges.pop(), 0.5, ref_edges.pop());

  new_bodies.get()->ref_vertices(ref_vertices);
  new_bodies.clean_out();
  int v_size = ref_vertices.size();
  for (int i=0; i < v_size; i++)
  {
    RefVertex* vertex = ref_vertices.get_and_step();
    if(vertex->coordinates().x() != 16)
    {
      ref_vertices.remove(vertex); 
      continue;
    }
    if((vertex->coordinates().y()== 5 && vertex->coordinates().z() == 5)
   || (vertex->coordinates().y()== -5 && vertex->coordinates().z() == -5))
    continue;
    else
      ref_vertices.remove(vertex);
  }
   
  assert (ref_vertices.size() ==2);
  gmti->tweak_chamfer(ref_vertices, 1, new_bodies);

  ref_edges.clean_out();
  new_bodies.get()->ref_edges(ref_edges);
  isize=ref_edges.size();
  assert(isize == 26);
 
  CubitStatus rsl = CUBIT_SUCCESS;
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "fillet.occ";
  const char * filetype = "OCC";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 5);
  remove(filename);

  DLIList<Body*> bodies;
  DLIList<RefEntity*>  free_entities;
  gti->bodies(bodies);
  gti->get_free_ref_entities(free_entities);

  //delete all entities
  gti->delete_Body(bodies);

  //bug fixing for creating arc curve.
   CubitVector vti(6.0, 0.0, 0.0);
   RefVertex* vert1 = GeometryModifyTool::instance()->make_RefVertex(vti);
   CubitVector vtii(6.0, 8.0, 0.0);
   RefVertex* vert2 = GeometryModifyTool::instance()->make_RefVertex(vtii);
   CubitVector v3(0, 8.0, 0.0);
   RefEdge* edge6 = gmti->make_RefEdge(ARC_CURVE_TYPE, vert1, vert2, &v3 );
   double d = edge6->measure();
   assert(d - 22.14297 < 0.00001);

   DLIList<RefEntity*> listVertexEdge;
   edge6->get_child_ref_entities(listVertexEdge);
   int numVertexEdge = listVertexEdge.size();
   RefVertex* tmpVertex = NULL;
   CubitVector coord;
   for (int i = 0; i < numVertexEdge; i++)
   {
       tmpVertex = CAST_TO(listVertexEdge.get_and_step(), RefVertex);
       coord = tmpVertex->coordinates();
   }
   int numTotalVertex  = GeometryQueryTool::instance()->num_ref_vertices();
   DLIList<RefVertex*> ref_v;
   GeometryQueryTool::instance()->ref_vertices(ref_v);
   for (int i = 0; i < numTotalVertex; i++)
   {
       tmpVertex = ref_v.get_and_step();
       coord = tmpVertex->coordinates();
   }

  bodies.clean_out();
  free_entities.clean_out();
  gti->bodies(bodies);
  gti->get_free_ref_entities(free_entities);

  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  //other tests
  body = gmti->brick(10, 10, 10);
  CubitVector v(15,0,0);
  BodySM* bodysm = body->get_body_sm_ptr();
  DLIList<OCCSurface*> occ_surfaces;
  CAST_TO(bodysm, OCCBody)->get_all_surfaces(occ_surfaces);  
  ref_faces.clean_out();
  ref_edges.clean_out();
  body->ref_faces(ref_faces);

  DLIList<RefFace*> faces_to_stitch;
  RefFace* test_face = NULL;
  for(int i = 0 ; i < ref_faces.size(); i++)
  {
    RefFace* refface = gmti->make_RefFace(ref_faces.get_and_step());
    if(refface)
	faces_to_stitch.append(refface);
  }

  DLIList<Surface*> surface_list;
  test_face = gmti->make_RefFace(ref_faces.get());
  surface_list.append(test_face->get_surface_ptr());
  CubitVector v_test(-5,0,0);
  CubitVector normal, test_normal(-1,0,0);
  normal = test_face->normal_at(v_test); //(-1,0,0)
  assert(normal == test_normal);
  ome->flip_normals(surface_list);
  normal = test_face->normal_at(v_test); //(1,0, 0)
  assert(normal == -test_normal);
  bodies.clean_out();
  test_face->bodies(bodies);
  gti->delete_Body(bodies);
  surface_list.clean_out();
  bodies.clean_out(); 
  gti->delete_Body(body);

  DLIList<BodySM*> bodysm_list;
  DLIList<RefFace*> face_list;
  DLIList<RefVertex*> vertices;
  for(int i = 0; i < faces_to_stitch.size(); i++)
  {
    //move each refface by (15,0,0)
    RefFace* refface = faces_to_stitch.get_and_step();
    refface->ref_vertices(vertices);
    Body* body = refface->ref_volume()->get_body_ptr();
    bodysm = body->get_body_sm_ptr();
    Surface* surface = refface->get_surface_ptr();
    surface_list.append(surface);
    bodysm_list.append(bodysm);
    bods.clean_out();
    bods.append(body);
    gti->translate(bods, v);
    body->ref_faces(face_list);
    occ_surfaces.clean_out();
    CAST_TO(bodysm, OCCBody)->get_all_surfaces(occ_surfaces);
  }

  //create solid from surfaces 
  //GeometryModifyEngine *gme = gmti->get_engine(occ_surfaces.get());
  new_bodies.clean_out();
  gmti->create_solid_bodies_from_surfs(face_list, new_bodies);

  ref_entity_list.clean_out();
  num_ents_exported=0;
  filename = "stitch.occ";

  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 1);

  gti->bodies(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  // Read in the geometry from files specified on the command line
  const char *argv = "stitch.occ";
  status = read_geometry(1, &argv, true);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 1 volume.
  remove(filename);

  bodies.clean_out();
  free_entities.clean_out();
  gti->bodies(bodies); 
  gti->get_free_ref_entities(free_entities);

  //delete all entities
  gti->delete_Body(bodies);

  for (int j = free_entities.size(); j--;)
    {
      gti->delete_RefEntity( free_entities.get_and_step());
    }

  //test for cylinder making
  Body* from_body = gmti->cylinder(10, 4, 3, 2);
  Body* from_body2 = gmti->cylinder(8, 4, 2, 0);
  Body* tool_body = gmti->cylinder(10, 1, 1, 1);
  d = from_body->measure(); //d = 219.91
  assert(d - 219.91 < 0.01 && d > 219.91);
  bodysm = from_body->get_body_sm_ptr();
  CubitVector center1;
  CAST_TO(bodysm, OCCBody)->mass_properties(center1, d);
  d = from_body2->measure();//d = 67.02
  assert(d - 67.02 < 0.01&& d > 67.02);
  d = tool_body->measure(); //d = 31.41
  assert(d - 31.41 < 0.01 && d > 31.41);
  
  // test for sphere making
  Body* test_body = gmti->sphere(5);
  d = test_body->measure(); //d = 523.598 
  assert(d - 523.598 < 0.01 && d >  523.598);
  //test for prism making
  test_body = gmti->prism(10, 4, 4,2);
  d =  test_body->measure(); //d = 320
  assert(-d + 320 < 0.1 && d < 320);
  //test for pyramid making
  test_body = gmti->pyramid(10, 4, 5, 2, 3);
  d =  test_body->measure(); //d = 130.6666
  assert(d - 130.666 < 0.001&& d > 130.666);

  //test for torus making
  test_body =  gmti->torus(10,5);
  d =  test_body->measure(); //d = 4934.8
  assert(d - 4934.8 < 0.1&& d > 4934.8);

  CubitVector factors(1,1,0.5);
  gmti->scale(test_body, p1, factors);  
  d = test_body->measure();
  assert( d > 2475 && d < 2475.5);

  //test for planar_sheet making
  CubitVector p2(0, 0, 1);
  CubitVector p3(0, 1, 1);
  CubitVector p4(0, 1, 0);
  test_body =  gmti->planar_sheet(p1, p2, p3, p4);
  d =  test_body->measure(); //d = 1
  assert(1-d < 0.001 && d < 1);
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 

  //test for subtract
  from_body = gmti->brick(10, 10, 10);
  int width = 10; //we can also test for width < 10
  from_body2 = gmti->brick(4, width, 4);
  tool_body = gmti->brick(10, 10, 10);  
  CubitVector v_move(3,width/2.0,1);
  CubitVector v_movei(5,5,6);
  CubitVector v_moveii(5,5,5);
  DLIList<Body*> from_bods2;
  from_bods2.append(from_body2);
  gti->translate(from_bods2,v_move);
  DLIList<Body*> tool_bods;
  tool_bods.append(tool_body);
  gti->translate(tool_bods, v_movei);
  DLIList<Body*> from_bods;
  from_bods.append(from_body);
  gti->translate(from_bods, v_moveii);
  DLIList<Body*> from_bodies;
  from_bodies.append(from_body);
  new_bodies.clean_out();
  rsl = gmti->subtract(from_body2, from_bodies, new_bodies, 
                       CUBIT_TRUE, CUBIT_FALSE);
  //new bodies has one body, new body has 10 ref-faces, 5 of them are remaining
  //with old id, 5 of them are new faces.
  assert(new_bodies.size() == 1);
  assert(new_bodies.get()->num_ref_faces() == 10);

  from_bodies=new_bodies;
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body,from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);
  d = new_bodies.step_and_get()->measure(); //d = 50
  v = new_bodies.get()->center_point(); //v = (5, 5, .05)
  n = new_bodies.get()->num_ref_faces();
  assert(new_bodies.size() == 1);
  assert( n == 12);
  assert (60 -d < 0.0001 && d <= 60);
  CubitVector test_v(5, 5, .5);
  assert (v == test_v);
  //new bodies has 1 body with 2 volumes, a volume = 10 and the other  
  //volume = 50; each of them has 6 ref_faces, of which 3 are new and 3 are
  //remaining (unchanged or modified).

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies); 
  
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //there shouldn't be any free_entites.

  from_body = gmti->brick(10, 10, 10);
  tool_body = gmti->brick(11, 1, 1);
  CubitVector v_move4(5,6,4);
  CubitVector v_move4i(5.5,0.5,0.5);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods,v_move4);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move4i);
  Body* cp_from_body = gmti->copy_body(from_body);
  Body* cp_from_body2 = gmti->copy_body(cp_from_body);

  //test edge imprint on body
  ref_edges.clean_out();
  tool_body->ref_edges(ref_edges);
  from_bodies.clean_out();
  from_bodies.append(cp_from_body);
  new_bodies.clean_out();
  CubitStatus stat = gmti->imprint(from_bodies, ref_edges, new_bodies, CUBIT_FALSE, CUBIT_TRUE );
  assert(new_bodies.get()->num_ref_faces() == 8);

  //test edge imprint on surface
  CubitVector vv(5,1,4.0);
  face_list.clean_out();
  cp_from_body2->ref_faces(face_list);
  int size = face_list.size();
  DLIList<RefFace*> unimprint_faces;
  for(int i = 0; i < size; i++)
  {
    CubitVector v = face_list.get()->center_point();
    if(!v.about_equal(vv))
      unimprint_faces.append(face_list.remove());
    else
      face_list.step();
  }
  assert(face_list.size() == 1);
  ref_edges.step_and_get()->ref_vertices(vertices);
  new_bodies.clean_out();
  stat = gmti->imprint(unimprint_faces, ref_edges, new_bodies, CUBIT_FALSE);
  assert(new_bodies.get()->num_ref_faces() == 6);

  new_bodies.clean_out();
  face_list.clean_out();

  //imprint a point on an edge, split it
  CubitVector pp1(10,1,8);
  CubitVector pp2(10,5,9);
  CubitVector pp3(10,1,6);
  DLIList<CubitVector> vectors;
  vectors.append(pp1);
  vectors.append(pp2);
  vectors.append(pp3);
  from_bodies.clean_out();
  from_bodies.append(cp_from_body2);
  new_bodies.clean_out(); 
  stat = gmti->imprint(from_bodies, vectors, new_bodies, CUBIT_FALSE);

  n = new_bodies.get()->num_ref_edges();//n = 17
  assert( n == 17);

  new_bodies.get()->ref_faces(face_list);
  for(int i = 0; i < size; i++)
  {
    CubitVector v = face_list.get()->center_point();
    if(!v.about_equal(vv))
      face_list.remove();
    else
      face_list.step();
  }
  new_bodies.clean_out();
  stat = gmti->imprint(face_list, ref_edges, new_bodies, CUBIT_FALSE);

  n = new_bodies.get()->num_ref_edges();//n = 21
  assert( n == 21);

  //test for multi-cut imprint for subtract.
  from_bodies.clean_out();
  from_bodies.append(from_body);
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_FALSE);

  n = new_bodies.get()->num_ref_faces();
  assert(n == 8);
  n = new_bodies.get()->num_ref_edges();
  assert(n == 18);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //there shouldn't be any free_entites.

  //test for shell body subtract.
  tool_body = gmti->brick(1, 1, 1);
  CubitVector new_move(0.5,0.5,0.5);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,new_move);
  //just need two surfaces.
  DLIList<RefFace*> reffaces;
  tool_body->ref_faces(reffaces);
  CubitVector v1(1, 0.5, 0.5);
  CubitVector v2(0.5, 0.5, 1);
  DLIList<Surface*> surfaces;
  for (int i = 0; i < tool_body->num_ref_faces(); i++)
  {
    RefFace* face = reffaces.get_and_step();
    CubitVector v = face->center_point();
    if (v.about_equal(v1) || v.about_equal(v2))
      surfaces.append(face->get_surface_ptr());
  } 
  assert(surfaces.size() == 2);
  DLIList<BodySM*> body_list;
  for (int i = 0; i < surfaces.size(); i++)
  {
    Surface* surface = surfaces.get_and_step();
    surface = ome->make_Surface(surface);
    body_list.append(CAST_TO(surface,OCCSurface)->my_body());
  } 
  bodysm = NULL;

  //test stitch surfaces operation
  DLIList<BodySM*> new_surf_bodies;
  ome->stitch(body_list, new_surf_bodies, CUBIT_FALSE, 1E-6);
  assert(new_surf_bodies.size() == 1);
 
  bodies.clean_out();
  bodies.append(tool_body);
  gti->delete_Body(bodies);

  //test flip_normal for a shell body.
  bodysm = new_surf_bodies.get();
  from_body2 = gti->make_Body(bodysm); 
  test_body = gmti->copy_body(from_body2);
  is_sheet = test_body->is_sheet_body();
  assert(is_sheet == CUBIT_FALSE); //2 faces shell

  RefEntity* ref_body =  CAST_TO(test_body, RefEntity);
  DLIList<RefEntity*> ref_entities;
  ref_faces.clean_out();
  ref_body->get_all_child_ref_entities(ref_entities);
  for(int i = 0; i < ref_entities.size(); i++) 
  {
    RefFace *face = CAST_TO(ref_entities.get_and_step(), RefFace);
    if (face != NULL)
      ref_faces.append(face);
  }
  CubitVector test_normal0(1,0,0);
  CubitVector test_normal2(0,0,1);
  normal = ref_faces.get()->normal_at(v1); //(1,0,0)
  if (normal != test_normal0)
  {
    assert( normal ==test_normal2);
    ref_faces.step();
  }
  normal = ref_faces.step_and_get()->normal_at(v2); //(0,0,1)
  assert(normal == test_normal2);
  surfaces.clean_out();
  surfaces.append(ref_faces.step_and_get()->get_surface_ptr());
  ome->flip_normals(surfaces);
  normal = ref_faces.get()->normal_at(v1); //(-1,0,0)
  assert(normal == -test_normal0);
  normal = ref_faces.step_and_get()->normal_at(v2); //(0,0,1)
  assert(normal == test_normal2);
  gti->delete_Body(test_body);
  surfaces.clean_out();
   

  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move3(2,3,2);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move3);
  Body* copy_tool_body = gmti->copy_body(tool_body);
  Body* copy_tool_body2 = gmti->copy_body(tool_body);

  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(copy_tool_body);
  from_bodies.append(from_body2);
  stat = gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE); 
  assert(new_bodies.size() == 1 && from_bodies.size() == 2);
  assert(new_bodies.get()->num_ref_faces() == 7);

  //test body imprinted by curves.
  ref_edges.clean_out();
  from_body2->ref_edges(ref_edges);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(copy_tool_body2);
  stat = gmti->imprint(from_bodies, ref_edges, new_bodies, CUBIT_FALSE, CUBIT_TRUE );
  assert(new_bodies.get()->num_ref_faces() == 6);

  //test body cutting a shell, one surface got cut as the result. 
  CubitVector v_move6(1,-1,0);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move6);
  from_bodies.clean_out();
  from_bodies.append(from_body2);
  new_bodies.clean_out();
  rsl = gmti->subtract(tool_body, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);

  d = new_bodies.step_and_get()->measure();
  v = new_bodies.get()->center_point();
  n = new_bodies.get()->num_ref_faces();
  assert( n == 1);
  assert(1-d<0.00001 && d < 1);
  CubitVector test_v1(0.5,0.5,1);
  assert(v == test_v1);
 
  from_bodies.clean_out();
  from_bodies.append(tool_body);
  
  //test a shell cutting a body, failed operation with an Error message.
  rsl = gmti->subtract(from_body2, from_bodies, new_bodies,
                       CUBIT_TRUE, CUBIT_TRUE);
  //WARNING: Surfaces or Shells can't be used to cut a body.
  //ERROR: Subtract FAILED
  assert (rsl == CUBIT_FAILURE);

  //test solid solid imprint
  tool_body  = gmti->brick(4, 4, 4);
  CubitVector v_move5(2,2.5,2);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move5);
  from_body = gmti->brick(1,1,1);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  stat =  gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE);
  //one body gets 4 cuts, the other 3.
  //one body gets 10 faces, the other 9.
  assert(new_bodies.size() == 2);
  assert(new_bodies.get()->num_ref_faces() == 10);
  assert(new_bodies.step_and_get()->num_ref_faces() == 9);
 
  //test imprint projected edges
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);  
  assert(free_entities.size() == 0);
 
  from_body  = gmti->brick(4, 4, 4);
  CubitVector v_move9(2,5,2);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods,v_move9);
  ref_faces.clean_out();
  from_body->ref_faces(ref_faces);
  tool_body = gmti->brick(11, 1, 1);  
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods, v_move4i);
  ref_edges.clean_out();
  tool_body->ref_edges(ref_edges);
  new_bodies.clean_out();
  gmti->imprint_projected_edges(ref_faces,ref_edges, new_bodies, CUBIT_TRUE,
       CUBIT_FALSE); 
  if(new_bodies.size())
  {
    n = new_bodies.get()->num_ref_faces();//n = 8, new_bodies.size() == 1
    assert(n == 8 && new_bodies.size() == 1);
  }
  //delete all bodies.
  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //test body-body intersect.
  //1. from body is the commom body, no update
  tool_body  = gmti->brick(4, 4, 4);
  from_body = gmti->brick(1,1,1);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move);
  CubitVector new_move2(2, 2, 2);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods, new_move2);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 1
  assert(1-d < 0.00001 && d <= 1); 

  //2. common body is part of from body, update the correponding face
  tool_body  = gmti->brick(4, 4, 4);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move5);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  new_bodies.clean_out();
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 0.5
  assert(0.5-d < 0.0001 && d < 0.5);

  //3. there's no common body, from body is deleted or kept depending on
  //keep-old flag. 
  tool_body  = gmti->brick(4, 4, 4);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move5);
  CubitVector v_m(0, 0.5, 0);
  gti->translate(tool_bods,v_m);
  from_bodies.clean_out();
  new_bodies.clean_out();
  from_bodies.append(from_body);
  //nothing changed 
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_TRUE);
  //"The 1 body did not have common part with the tool_body."
  assert(0 == new_bodies.size());

  //from_body get deleted
  stat =  gmti->intersect(tool_body, from_bodies, new_bodies, CUBIT_FALSE);
  //"The 1 body did not have common part with the tool_body."
  assert(0 == new_bodies.size());
  //Destroyed volume(s): 44, 46

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);

  //test chop operation
  tool_body  = gmti->brick(4, 4, 4);
  from_body = gmti->brick(1,1,1);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move5);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();

  //new_bodies = intersect bodies;bodies = outside bodies; body is dummy 
  stat =  gmti->chop(from_bodies, new_bodies, bodies, body, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 0.5
  assert(0.5-d<0.00001 && d < 0.5);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);
  //test chop 2
  from_body = gmti->brick(4, 4, 4);
  tool_body = gmti->brick(1,1,1);
  CubitVector v_move7(0.5, 1, 0.5);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move7);

  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move2);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  bodies.clean_out();
  stat =  gmti->chop(from_bodies, new_bodies, bodies, body, CUBIT_FALSE);
  d = bodies.get()->measure();//d = 63
  assert(d-63<0.0001 && d>63); 

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  assert(free_entities.size() == 0);

  //test unite 1
  tool_body = gmti->brick(1,1,1);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move7);
  gti->translate(tool_bods,v_m);
  from_body = gmti->brick(1,1,1);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move);
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 2
  n = new_bodies.get()->num_ref_faces(); //n = 10  
  assert(fabs(d-2)<0.0001 && n == 10); 

  filename = "glue.occ";
  filetype = "OCC";
  ref_entity_list.clean_out();
 
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                num_ents_exported, cubit_version);
  remove(filename);

  bodies.clean_out();
  gti->bodies(bodies); //bodies.size() = 1
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //delete all entities
  gti->delete_Body(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //test unite 2
  tool_body = gmti->brick(4, 4, 4);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move5);
  from_body = gmti->brick(1,1,1);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, new_move);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 64.5
  n = new_bodies.get()->num_ref_faces(); 
  assert(fabs(d-64.5)<0.00001 && n == 13);
  //n = 13. because fusion keeps the imprintings.

  bodies.clean_out();
  gti->bodies(bodies); //bodies.size() = 1
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //delete all entities
  gti->delete_Body(bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);

  //test making thick body.
  from_body = gmti->cylinder(10, 4, 4, 4); 
  tool_body = gmti->cylinder(5,1,1,1);
  CubitVector v_move8(0,0,12.5);
  CubitVector v_move8i(0 , 0 , 5);
  tool_bods.clean_out();
  tool_bods.append(tool_body);
  gti->translate(tool_bods,v_move8);
  from_bods.clean_out();
  from_bods.append(from_body);
  gti->translate(from_bods, v_move8i);
  from_bodies.clean_out();
  from_bodies.append(from_body);
  from_bodies.append(tool_body);
  new_bodies.clean_out();
  stat = gmti->unite(from_bodies, new_bodies, CUBIT_FALSE);
  d = new_bodies.get()->measure(); //d = 518.3627
  n = new_bodies.get()->num_ref_faces(); //n = 5
  assert(d-518.3627<0.0001 && d > 518.3627 && n==5);

  //find the top most surface as the opening of the thick body.
  ref_faces.clean_out();
  new_bodies.get()->ref_faces(ref_faces);
  CubitVector center(0,0,15);
  for(int i = 0; i < n; i++)
  {
    if(ref_faces.step_and_get()->is_planar() && 
       ref_faces.get()->center_point() == center )
      break;
  }
  RefFace* sweep_face = gmti->make_RefFace(ref_faces.get());
  DLIList<RefFace*> faces_to_remove;
  faces_to_remove.append(ref_faces.get());
  from_bodies = new_bodies;
  new_bodies.clean_out();
  stat = gmti->hollow(from_bodies, faces_to_remove, new_bodies, -.2);
  n = new_bodies.get()->num_ref_faces(); //n = 9
  d = new_bodies.get()->measure(); //d = 72.4074
  assert(n == 9 && d-72.4074<0.0001 && d>72.4074);

  bodies.clean_out();
  gti->bodies(bodies); //bodies.size() = 2
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);
  //delete hollowed body and keep sweep_face. 
  gti->delete_Body(new_bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities); //free_entities.size() = 0
  assert(free_entities.size() == 0);

  DLIList<RefEntity*> refentities;
  refentities.append(sweep_face);
  RefFace* draft_face = gmti->make_RefFace(sweep_face);
  RefFace* draft_face2 = gmti->make_RefFace(sweep_face);
  RefFace* perp_face = gmti->make_RefFace(sweep_face);
  RefFace* rotate_face = gmti->make_RefFace(sweep_face);
 
  new_bodies.clean_out();
  gmti->sweep_translational(refentities, v_move8, 0.0, 1, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, new_bodies);  
  body = new_bodies.get();
  d = body->measure();
  assert(d - 39.2699 < 0.0001 && d >39.2699);

  refentities.clean_out();
  refentities.append(draft_face);
  CubitVector v_move8ii(0,0,10);
  new_bodies.clean_out();
  gmti->sweep_translational(refentities, v_move8ii, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, new_bodies); 
  body = new_bodies.get();
  d = body->measure();
  //d = 66.3676  theoretical calculation is 66.7833, error 0.62%
  assert(d - 66.3676 < 0.0001 && d > 66.3676);

  v_move8ii.z(-10);
  refentities.clean_out();
  refentities.append(draft_face2);
  new_bodies.clean_out();
  gmti->sweep_translational(refentities, v_move8ii, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, new_bodies);
  body = new_bodies.get();
  d = body->measure();
  //d = 66.3676  theoretical calculation is 66.7833, error 0.62%
  assert(d - 66.3676 < 0.0001 && d > 66.3676);

  v_move8ii.z(10);
  DLIList<RefEdge*> edges;
  body->ref_edges(edges);
  refentities.clean_out();
  for (int i = 0; i < edges.size(); i++)
  {
    RefEdge* edge = edges.get();
    double length = edge->measure();
    if (length > 6.28  && length < 6.3)
    {
      refentities.append(edge);
      break;
    }
    edges.step();
  }
  new_bodies.clean_out();
  gmti->sweep_translational(refentities, v_move8ii, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, new_bodies);
  body = new_bodies.get();
  d = body->measure();
  //d = area = 90.1292 theoretica calculation is 90.5754, error 0.49%
  assert(d - 90.1292 < 0.0001 && d > 90.1292);

  refentities.clean_out();
  refentities.append(perp_face);
  new_bodies.clean_out();
  gmti->sweep_perpendicular(refentities, 10.0, 0.087, 1, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, new_bodies);
  body = new_bodies.get();
  d = body->measure();
  //d = 66.3676  theoretical calculation is 66.7833, error 0.62%
  assert(d - 66.3676 < 0.0001 && d > 66.3676);
  
  //Testing for sweep_rotational function
  //1. Negative testing: surface rotates along an intersection axis, fails.
  refentities.clean_out();
  refentities.append(rotate_face); 
  DLIList<RefEdge*> rotate_edges;
  rotate_face->ref_edges(rotate_edges);
  CubitVector sweep_axis(1,0,0);   
  CubitVector point(0, 2, 15);
  double angle = 1.57;
  new_bodies.clean_out();
  stat =gmti->sweep_rotational(refentities, center, sweep_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0.0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  assert(stat == CUBIT_FAILURE);

  //2. surface rotates along a non-intersection axis
  refentities.append(rotate_face);
  new_bodies.clean_out();
  gmti->sweep_rotational(refentities, point, sweep_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(d - 9.8646 < 0.0001 && d > 9.8646);

  //3. closed curve rotates along an intersecting axis, fails
  refentities.clean_out();
  refentities.append(rotate_edges.get());
  new_bodies.clean_out();
  stat =gmti->sweep_rotational(refentities, center, sweep_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE); 
  assert(stat == CUBIT_FAILURE);

  //4. closed curve rotates along a non-intersecting axis with make_solid
  //   option set to be true.
  refentities.append(rotate_edges.get());
  gmti->sweep_rotational(refentities, point, sweep_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(d - 9.8646 < 0.0001 && d > 9.8646);

  //5. closed curve rotates along a non-intersecting axis with make_solid
  //   option set to be false.
  refentities.clean_out();
  refentities.append(rotate_edges.get());
  new_bodies.clean_out();
  gmti->sweep_rotational(refentities, point, sweep_axis, -angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(d - 19.7292 < 0.0001 && d> 19.7292);

  //6. open curve rotates along an intersecting axis with make_solid
  //   option set to be true. Failed
  CubitVector pt1(0,1,15);
  CubitVector pt2(0,-1,15);
  CubitVector pt3(-1,0,15);
  RefVertex* vt1 = gmti->make_RefVertex(pt1);
  RefVertex* vt2 = gmti->make_RefVertex(pt2);
  RefEdge* edge1 = gmti->make_RefEdge(ARC_CURVE_TYPE, vt1, vt2, &pt3); 
  CubitVector apoint(0,-0.5,15);
  refentities.clean_out();
  refentities.append(edge1);
  new_bodies.clean_out();
  stat=gmti->sweep_rotational(refentities, apoint, sweep_axis, angle, 
                              new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                              CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  assert(stat == CUBIT_FAILURE);

  //7. open curve rotates along an intersecting axis at end points with 
  //   make_solid option set to be true.
  refentities.append(edge1);
  CubitVector rotate_axis(0, 1, 0);
  gmti->sweep_rotational(refentities, center, rotate_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(1.046667-d < 0.0001 && d < 1.046667);

  //8. open curve rotates along an intersecting axis at end points with
  //   make_solid option set to be false.
  refentities.clean_out();
  refentities.append(edge1);
  new_bodies.clean_out();
  gmti->sweep_rotational(refentities, center, rotate_axis, -angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert( d - 3.14 > -0.0001 && d < 3.14);

  //9. staight line rotate about itself fails.
  RefEdge* edge2 = gmti->make_RefEdge(vt1, vt2, rotate_face);
  refentities.clean_out();
  refentities.append(edge2);
  new_bodies.clean_out();
  stat=gmti->sweep_rotational(refentities, center, rotate_axis, -angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE); 
  assert(stat == CUBIT_FAILURE);

  //10. open curve rotates along a non-intersect axis with make_solid option
  //    set to be true.
  refentities.clean_out();
  refentities.append(edge1);
  CubitVector off_center(1,0,15);
  gmti->sweep_rotational(refentities, off_center, rotate_axis, angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(d - 5.0828 < 0.0001 && d > 5.0828);

  //11. open curve rotates along a non-intersect axis with make_solid option
  //    set to be false.
  refentities.clean_out();
  refentities.append(edge1);
  new_bodies.clean_out();
  gmti->sweep_rotational(refentities, off_center, rotate_axis, -angle, 
                         new_bodies, CUBIT_FALSE, CUBIT_TRUE, 0, 0, 1,
                         CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE);
  body = new_bodies.get();
  d = body->measure();
  assert(d - 8.07227 < 0.001 && d > 8.07227);

  //test sweep_along_curve function
  //sweep along a straight curve with draft
  refentities.clean_out();
  refentities.append(rotate_face);
  CubitVector pt4(0,1,0);
  RefVertex* vt4 = gmti->make_RefVertex(pt4);
  RefEdge* edge3 = gmti->make_RefEdge(STRAIGHT_CURVE_TYPE, vt1, vt4);
  edges.clean_out();
  edges.append(edge3);
  new_bodies.clean_out();
  gmti->sweep_along_curve(refentities, edges, new_bodies, CUBIT_FALSE, 
                          CUBIT_TRUE,0.087,1);
  body = new_bodies.get();
  d = body->measure();
  //d = 134.7152, theoretical calculation is 135.66363, error 0.67%
  assert(d -134.7152 < 0.0001 && d > 134.7152);

  //sweep along two curves which make a G1 continous wire, no draft is performed
  refentities.clean_out();
  refentities.append(rotate_face);
  CubitVector pt5(0,2,7.5);
  CubitVector pt6(0,0,-7.5);
  CubitVector pt7(0,1,-15);
  RefVertex* vt5 = gmti->make_RefVertex(pt7);
  DLIList<CubitVector*> vector_list;
  CubitVector coord1 = vt1->coordinates();
  CubitVector coord2 = vt4->coordinates();
  vector_list.append(&coord1);
  vector_list.append(&pt5);
  vector_list.append(&coord2);
  
  RefEdge* edge4 = gmti->make_RefEdge( SPLINE_CURVE_TYPE, vt1, vt4, vector_list);
  vector_list.clean_out();
  vector_list.append(&coord2);
  vector_list.append(&pt6);
  CubitVector coord3 = vt5->coordinates();
  vector_list.append(&coord3);
  RefEdge* edge5 = gmti->make_RefEdge( SPLINE_CURVE_TYPE, vt4, vt5, vector_list);  
  edges.clean_out();
  edges.append(edge4);
  edges.append(edge5);
  new_bodies.clean_out();
  gmti->sweep_along_curve(refentities, edges, new_bodies, CUBIT_FALSE, 
                          CUBIT_TRUE, 0.087,1);
  body = new_bodies.get();
  d = body->measure();
  //d = 93.697, no effect of draft angle.
#if OCC_VERSION_MINOR > 5
  assert(d - 92.1335 < 0.001 && d > 92.1335);
#else
  OCCQueryEngine* oqe = OCCQueryEngine::instance();
  if(oqe->get_subminor_version() >= 2)
    assert(d - 92.1335 < 0.001 && d > 92.1335);
  else
    assert(d - 93.697 < 0.001 && d > 93.697);
#endif

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  for(int i = 0; i < free_entities.size(); i++)
    gti->delete_RefEntity(free_entities.get_and_step());

  //So far created volume 80.

  // Read in the geometry from files specified on the command line
  argv = "Cylinder_1.brep";
  status = read_geometry(1, &argv);
  argv = "Cylinder_2.brep";
  status = read_geometry(1, &argv);
  if (status == CUBIT_FAILURE) exit(1);
  //Read in 2 volumes.

  from_bodies.clean_out();
  new_bodies.clean_out();
  gti->bodies(from_bodies);
  stat = gmti->imprint(from_bodies, new_bodies, CUBIT_FALSE);

  num_ents_exported=0;
  filename = "imprint.brep";
  ref_entity_list.clean_out();
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 2);
  assert (CAST_TO(ref_entity_list.get(), Body)->num_ref_faces() == 7);
  assert (CAST_TO(ref_entity_list.get_and_step(), Body)->num_ref_edges() == 11);
  assert (CAST_TO(ref_entity_list.get(), Body)->num_ref_volumes() == 1);
  remove(filename);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  for(int i = 0; i < free_entities.size(); i++)
    gti->delete_RefEntity(free_entities.get_and_step());

  // Read in the geometry from files specified on the command line
  argv = "webcut.brep";
  read_geometry(1, &argv, false);

  DLIList<Body*> old_bodies, junk;
  new_bodies.clean_out();
  gti->bodies(old_bodies);
  CubitVector center_l(0, 0, 0);
  CubitVector axis(0.,0.,1.);
  double radius = 0.5;
  rsl= gmti->webcut_with_cylinder(old_bodies,radius,axis, center_l,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  num_ents_exported=0;
  filename = "after_webcut.brep";
  ref_entity_list.clean_out();
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 2);
  remove(filename);

  bodies.clean_out();
  gti->bodies(bodies);
  //delete all entities
  gti->delete_Body(bodies);

  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  for(int i = 0; i < free_entities.size(); i++)
    gti->delete_RefEntity(free_entities.get_and_step());

  //test for webcut with curve sweep option
  argv = "webcut.brep";
  read_geometry(1, &argv, false);
  argv = "Line_1.brep";
  read_geometry(1, &argv, false);
  argv = "Line_2.brep";
  read_geometry(1, &argv, false);
  
  old_bodies.clean_out();
  new_bodies.clean_out();
  gti->bodies(old_bodies);
  free_entities.clean_out();
  gti->get_free_ref_entities(free_entities);
  DLIList<RefEdge*> curves;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* free_edge = CAST_TO(free_entities.get_and_step(), RefEdge);
    curves.append(free_edge);
  }

  axis.set(0.,1.,0.);

  rsl= gmti->webcut_with_sweep_curves(old_bodies,curves,axis,true,NULL,NULL,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;
  
  num_ents_exported=0;
  filename = "webcut_sweep.brep";
  ref_entity_list.clean_out();
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);
  assert(num_ents_exported == 4);
  remove(filename);
  return CUBIT_SUCCESS;
}
