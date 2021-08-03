//-------------------------------------------------------------------------
// Filename      : OCCQueryEngine.hpp
//
// Purpose       : OCC geometry engine.
//
//                 This class is implemented as a Singleton pattern. Only
//                 one instance is created and it is accessed through the 
//                 {instance()} static member function.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/07
//
//-------------------------------------------------------------------------

#ifndef OCC_GEOMETRY_ENGINE_HPP
#define OCC_GEOMETRY_ENGINE_HPP

// ********** BEGIN STANDARD INCLUDES         **********

#include <typeinfo>
#if !defined(NT) && !defined(CANT_USE_STD)
using std::type_info;
#endif

// ********** END STANDARD INCLUDES           **********
#include "TDF_Label.hxx"
// ********** BEGIN CUBIT INCLUDES            **********
#include "CubitFileIOWrapper.hpp"
#include "GeometryQueryEngine.hpp"
#include "Handle_TDocStd_Document.hxx"
#include <map>

// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS
class BRep_Builder;
class TopologyEntity;
class TopologyBridge;
class RefEntity;
class Body;
class Shell;
class ShellSM;
class Loop;
class Chain;
class LoopSM;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class TBPoint;
class Curve;
class Surface;
class Lump;
class BodySM;

class GeometryEntity;
class BodySM;
class ShellSM;
class Surface;
class Curve;
class CoEdgeSM;
class TopologyEntity;
class CubitBox;
class CubitString;
class CubitSimpleAttrib;

class OtherSolidModelingEntity;
class OCCLump;
class OCCShell;
class OCCLoop;
class OCCSurface;
class OCCBody;
class OCCCoEdge;
class OCCCurve;
class OCCPoint;
 
class TopTools_DataMapOfShapeInteger;
class BRepAlgoAPI_BooleanOperation;
class TopTools_IndexedMapOfShape;
class BRepBuilderAPI_ModifyShape;
class TopoDS_Vertex;
class TopoDS_Edge;
class TopoDS_Shape;
class TopoDS_Wire;
class TopoDS_Face;
class TopoDS_Solid;
class TopoDS_Shell;
class TopoDS_Compound;
// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********

// ********** END ENUM DEFINITIONS            **********
class OCCQueryEngine : public GeometryQueryEngine
{
public:
// ********** BEGIN FRIEND DECLARATIONS        **********
  friend class OCCSurface;
  
// ********** END FRIEND DECLARATIONS        **********

//HEADER- Constructor and Destructor
  
  static OCCQueryEngine* instance();
    //- Singleton pattern
    //- Controlled access and creation of the sole instance of this class.

  CubitBoolean EXPORT_ATTRIB;

  void copy_attributes(TopoDS_Shape& old_shape,
                       TopoDS_Shape& new_shape);

  int update_OCC_map(TopoDS_Shape& old_shape, TopoDS_Shape& new_shape);

  virtual ~OCCQueryEngine();
  
  const char* modeler_type()
     { return "OCC"; }

  TopoDS_Shape *get_TopoDS_Shape_of_entity( TopologyBridge * entity);
  
  void body_attributes_for_writing(DLIList<OCCBody*> &OCC_bodies, //input
                                 BRep_Builder &B, //input
                                 TopoDS_Compound &Co, //input and output,
                                 DLIList<OCCLump*> &single_lumps, //output
                                 DLIList< DLIList<CubitSimpleAttrib>*> &lists);

  int get_major_version();

  int get_minor_version();

  int get_subminor_version();

  CubitString get_engine_version_string();

//HEADER- RTTI and safe casting functions.
  
  virtual const type_info& entity_type_info() const ;
    //R- The geometric modeler type
    //- This function returns the type of the geometric modeler.
  
  virtual CubitBoolean is_solid_modeler_type() const 
    {return CUBIT_FALSE;}
    //R CubitBoolean
    //R- This  is not a solid modeling engine.
//HEADER- Functions for importing and exporting solid models.

  virtual CubitStatus restore_transform( BodySM* body );

  using GeometryQueryEngine::get_graphics;

  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                    GMem* gMem,
                                    unsigned short normal_tolerance = 15,
                                    double distance_tolerance = 0,
                                    double longest_edge = 0) const;
    CubitStatus get_graphics( TopoDS_Face* face_ptr,
                            GMem* g_mem,                                          
                            unsigned short normal_tolerance = 15,
                            double distance_tolerance = 0,
                            double max_edge_length = 0) const;


  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    GMem* gMem = NULL,
                                    double angle_tolerance=0,
                                    double distance_tolerance=0,
                                    double max_edge_length = 0.0 ) const;



  CubitStatus get_graphics( TopoDS_Edge* edge_ptr,
    GMem* gMem = NULL,
    double angle_tolerance = 0,
    double distance_tolerance = 0,
    double max_edge_length = 0 ) const;


    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_edge_ptr
    //I- The RefEdge for which hoops facetting information will be
    //I- gathered.
    //O numSteps
    //O- The number of points in resulting polyline.
    //O gMem
    //O- The storage place for edges, involved in facetting.
    //I tolerance
    //I- The tolerance deviation used when facetting the curve (optional
    //I- and currently IGNORED by this engine).
    //- This function gathers and outputs ACIS edge information for
    //- hoops involved in facetting a RefEdge.  If all goes well,
    //- CUBIT_SUCCESS is returned.  Otherwise, CUBIT_FAILURE is
    //- returned.
 
  virtual CubitStatus fire_ray( CubitVector &origin,
                                CubitVector &direction,
                                DLIList<TopologyBridge*> &at_entity_list,
                                DLIList<double> &ray_params,
                                int max_hits = 0,
                                double ray_radius = 0.0,
                                DLIList<TopologyBridge*> *hit_entity_list=0 ) const ;
    //- Fire a ray at specified entities, returning the parameters (distances)
    //- along the ray and optionally the entities hit.  Returned lists are
    //- appended to.  Input entities can be any of bodies, volumes, faces,
    //- edges or vertices.  Optionally you can specify the maximum number of
    //- hits to return (default = 0 = unlimited), and the ray radius to use for
    //- intersecting the entities (default = 0.0 = use modeller default).

  virtual CubitStatus get_isoparametric_points(Surface* ,
                                               int&, int&,
                                               GMem*&) const;
  
  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const;
  
  virtual CubitStatus transform_vec_position( 
    CubitVector const& ,
    BodySM *,
    CubitVector & ) const;
  
  virtual CubitStatus get_intersections( Curve*, CubitVector& point1,
                                         CubitVector&,
                                         DLIList<CubitVector>& ,
                                         CubitBoolean,
                                         CubitBoolean );

  virtual CubitStatus get_intersections( Curve* , Curve* ,
                                         DLIList<CubitVector>& ,
                                         CubitBoolean,
                                         CubitBoolean );
  virtual CubitStatus get_intersections( Curve* ref_edge, Surface* ref_face,
                                        DLIList<CubitVector>& intersection_list,
                                        CubitBoolean bounded = CUBIT_FALSE );

  virtual CubitStatus entity_extrema( DLIList<GeometryEntity*> &ref_entity_list, 
                                      const CubitVector *dir1, 
                                      const CubitVector *dir2,
                                      const CubitVector *dir3, 
                                      CubitVector &extrema,
                                      GeometryEntity *&extrema_entity_ptr );
  //- Gets the extrema position along the first given direction. If there 
  //- is more than one extrema position, the other directions will be used 
  //- to determine a unique position.  Directions 2 and 3 can be NULL.
  //- Entities supported include bodies, volumes, surfaces, curves and
  //- vertices.  The entity the extrema is found on is also returned.

  virtual CubitStatus entity_entity_distance( GeometryEntity *ref_entity_ptr1,
                                              GeometryEntity *ref_entity_ptr2,
                                              CubitVector &pos1, CubitVector &pos2,
                                              double &distance );
  //- Gets the minimum distance between two entities and the closest positions 
  //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                          const char* file_name,
                                          Model_File_Type  file_type,
                                          const CubitString &cubit_version,
                                          ModelExportOptions &export_options );


  // write shapes to buffer as binary format
  virtual CubitStatus export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
					  char*& p_buffer,
					  int& n_buffer_size,
					  bool b_export_buffer);

  virtual CubitStatus save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
                                          const char *file_name,
                                          const CubitString &cubit_version,
                                          CubitString &created_file,
                                          CubitString &created_file_type ); 

 virtual CubitStatus import_temp_geom_file(FILE* file_ptr,
                                      const char* file_name,
                                      Model_File_Type file_type,
                                      DLIList<TopologyBridge*> &bridge_list );

 virtual CubitStatus import_solid_model(const char* file_name,
                                        Model_File_Type file_type,
                                        DLIList<TopologyBridge*>& imported_entities,
                                        ModelImportOptions& options);
  virtual CubitStatus import_solid_model(DLIList<TopologyBridge*> &imported_entities,
					 const char* pBuffer,
					 const int n_buffer_size);
    
  CubitStatus unhook_BodySM_from_OCC( BodySM* bodysm,
                                    bool remove_lower_entities=CUBIT_TRUE)const;
  CubitStatus unhook_Surface_from_OCC( Surface* surface) const;
  CubitStatus unhook_Curve_from_OCC( Curve* curve) const;
  CubitStatus unhook_Point_from_OCC( TBPoint* point) const;
  void bound_TopoDS_Shape(const TopoDS_Shape & aShape);

private:
  CubitStatus import_solid_model(FILE *file_ptr,
                                 const char* /*file_type*/,
                                 DLIList<TopologyBridge*> &imported_entities,
                                 CubitBoolean print_results = CUBIT_TRUE,
                                 const char* logfile_name = NULL,
                                 CubitBoolean heal_step = CUBIT_TRUE,
                                 CubitBoolean import_bodies = CUBIT_TRUE,
                                 CubitBoolean import_surfaces = CUBIT_TRUE,
                                 CubitBoolean import_curves = CUBIT_TRUE,
                                 CubitBoolean import_vertices = CUBIT_TRUE,
                                 CubitBoolean free_surfaces = CUBIT_TRUE);

  CubitStatus unhook_Lump_from_OCC( Lump* lump ) const;
  CubitStatus unhook_ShellSM_from_OCC( ShellSM* shell ) const;
  CubitStatus unhook_CoEdges_from_OCC( DLIList<OCCCoEdge*>& coedges) const;
  CubitStatus unhook_LoopSM_from_OCC( LoopSM* loopsm) const;
  CubitStatus delete_loop( LoopSM* loopsm)const;
  void unhook_coedges_of_a_curve(OCCCurve* curve,
                                 OCCLoop*  loop)const;

  void add_shape_to_map(TopoDS_Shape& sh,
                        TopoDS_Shape& aShape, /*In, parent shape*/
                        int &current_id /*Out*/);
public:
  virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list)const;
    //- Deletes all solid model entities associated with the Bodies in 
    //- the input list. 
  void delete_bodies(DLIList<BodySM*>& body_list,
                     bool remove_lower_entities =CUBIT_TRUE) const;
  using GeometryQueryEngine::delete_solid_model_entities;
  virtual CubitStatus delete_solid_model_entities(
          GeometryEntity* ref_entity_ptr,
          bool remove_lower_entities) const;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr)const; 
  CubitStatus delete_body(BodySM* body_ptr ,
                          bool remove_lower_entities =CUBIT_TRUE) const;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr)const;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr)const; 
  virtual CubitStatus delete_solid_model_entities( TBPoint* point_ptr)const;

  virtual double get_sme_resabs_tolerance() const; // Gets solid modeler's resolution absolute tolerance
  virtual double set_sme_resabs_tolerance( double new_resabs );

  virtual CubitStatus set_int_option( const char* opt_name, int val );
  virtual CubitStatus set_dbl_option( const char* opt_name, double val );
  virtual CubitStatus set_str_option( const char* opt_name, const char* val );
    //- Set solid modeler options

  CubitStatus ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii);
  //- returns true in is_ascii if fp points to an ascii stl file

  CubitStatus create_super_bounding_box(
                                DLIList<BodySM*>& body_list,
                                CubitBox& super_box );

  //ModifyShape refers to only Transform and GTransform for now (1/10/11)
  CubitStatus update_entity_shape(GeometryEntity* entity_ptr,
                                  BRepBuilderAPI_ModifyShape* aTranf,
                                  BRepAlgoAPI_BooleanOperation *op = NULL);

  void set_TopoDS_Shape(TopologyBridge* tb, TopoDS_Shape& new_shape);
  CubitStatus translate( BodySM* body, const CubitVector& offset );
  CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle );
  CubitStatus scale    ( BodySM* body, double factor );
  CubitStatus scale    ( BodySM* body, const CubitVector& factors );
  CubitStatus reflect  ( BodySM* body, const CubitVector& axis );

  CubitStatus translate( GeometryEntity* ent, const CubitVector& offset );
  CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees );
  CubitStatus scale    ( GeometryEntity* ent, double factor );
  CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors );
  CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis );

  CubitStatus get_connected_patch( DLIList<OCCSurface*>& remaining_surfs,
                                   DLIList<OCCSurface*>& output_patch );
  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const;
  //R CubitBoolean
  //R- CUBIT_TRUE if the two bodies overlap, CUBIT_FALSE if they don't
  //R- overlap.  If the bodies are touching the function
  //R- should return CUBIT_FALSE.
  //I body_ptr_1, body_ptr_2
  //I- The two body pointers that are being tested for overlap.
  //-  The function uses the intersect call to test if the bodies
  //-  are overlaping.  The full intersect Boolean is needed to see if
  //-  the bodies actually overlap and don't just touch.

  TopologyBridge* occ_to_cgm(const TopoDS_Shape& shape);
  
  virtual CubitBoolean volumes_overlap (Lump *lump1, Lump *lump2 ) const ;

  DLIList<TopologyBridge*> populate_topology_bridge(TopoDS_Shape& aShape);
  BodySM* populate_topology_bridge(const TopoDS_Compound& aShape);
  Lump* populate_topology_bridge(const TopoDS_Solid& aShape, 
				   CubitBoolean build_body = CUBIT_FALSE);
  Surface* populate_topology_bridge(const TopoDS_Face& aShape,
                                    CubitBoolean build_body = CUBIT_FALSE);
  Curve* populate_topology_bridge(const TopoDS_Edge& aShape,
                                  CubitBoolean stand_along = CUBIT_FALSE );
  TBPoint* populate_topology_bridge(const TopoDS_Vertex& aShape,
                                  CubitBoolean stand_along = CUBIT_FALSE);

  OCCShell* populate_topology_bridge(const TopoDS_Shell& aShape,
                                     CubitBoolean standalone = CUBIT_FALSE );
  DLIList<OCCBody*> *BodyList ;
  DLIList<OCCSurface*> *SurfaceList ;
  DLIList<OCCLoop*> *WireList; //standalone wire list
  DLIList<OCCCurve*> *CurveList ;
  Handle(TDocStd_Document) MyDF;
  TDF_Label mainLabel;
  TopTools_DataMapOfShapeInteger* OCCMap;
  std::map<int, TopologyBridge*>* OccToCGM;
  std::map<int, TDF_Label>* Shape_Label_Map;
  static int iTotalTBCreated ;
  static int total_coedges;
protected:
  
  OCCQueryEngine();
  
private:

  OCCLoop* populate_topology_bridge(const TopoDS_Wire& aShape,
				    CubitBoolean standalone = CUBIT_FALSE);  

  CubitStatus write_topology( const char* file_name, 
                              Model_File_Type file_type,
                              DLIList<OCCBody*> &facet_bodies,
                              DLIList<OCCSurface*> &facet_surfaces,
                              DLIList<OCCCurve*> &facet_curves,
                              DLIList<OCCPoint*> &facet_points );

  CubitStatus write_topology( char*& p_buffer,
			      int& n_buffer_size,
			      bool b_export_buffer,
			      DLIList<OCCBody*> &OCC_bodies,
			      DLIList<OCCSurface*> &OCC_surfaces,
			      DLIList<OCCCurve*> &OCC_curves,
			      DLIList<OCCPoint*> &OCC_points);
  
  CubitBoolean Write(const TopoDS_Shape& Sh,
                     const Standard_CString File,
                     TDF_Label label);

  CubitBoolean Write(const TopoDS_Shape& Sh,
		     char*& p_buffer,
		     int& n_buffer_size,
		     bool b_export_buffer,
                     TDF_Label label);
  
  CubitBoolean Read(TopoDS_Shape& Shapes,
                    const Standard_CString File,
                    TDF_Label label);

  CubitBoolean Read(TopoDS_Shape& Sh,
		    const char* pBuffer,
		    const int n_buffer_size,
                    TDF_Label label);

  static OCCQueryEngine* instance_;
    //- static pointer to unique instance of this class
};

// ********** BEGIN INLINE FUNCTIONS          **********
// ********** END INLINE FUNCTIONS            **********

// ********** BEGIN FRIEND FUNCTIONS          **********
// ********** END FRIEND FUNCTIONS            **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END HELPER CLASS DECLARATIONS   **********

#endif
