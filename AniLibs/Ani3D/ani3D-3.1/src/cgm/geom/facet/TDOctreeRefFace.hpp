//- Class:       TDOctreeRefFace
//- Description: Tool data for storing additional information related to surface faceting
//- Owner:       W. R. Quadors
//- Checked by:
//- Version:

#ifndef TD_SKELETON_ref_FACE

#define TD_SKELETON_ref_FACE

#include "CubitDefines.h"
#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "CastTo.hpp"

class RefFace;
class CubitPoint;
class CubitFacet;
class CubitFacetEdge;
class CubitTransformMatrix;

  
class TDOctreeRefFace : public ToolData
{
private:

    //static MemoryManager memoryManager;
    //- memory management object

  DLIList<CubitPoint *> *ptrCubitPointList;

  DLIList<CubitFacetEdge *> *ptrCubitFacetEdgeList;

  DLIList<CubitFacet *> *ptrCubitFacetList;

  RefFace *refFace;

  DLIList<int> loopIndex;

  int lastCurveID;
	
  CubitBoolean visit;
    // This stores the starting curve id of every loop in a ref_face
    // Using this adjacency relation between the curves can be found
    // loopIndexLastCurveID stores the last curve id of all loops

    // If bad surface facetting is detected, don't create 2dmat
  CubitBoolean create_2dmat;
  
public:

  TDOctreeRefFace();
    //- constructor

  ~TDOctreeRefFace();
    //- destructor
 
  inline CubitBoolean get_visit( void )const{ return visit; }
  inline void set_visit( CubitBoolean type ){ visit = type; }

  void set_create_2dmat(const CubitBoolean new_val) {create_2dmat = new_val;}
  CubitBoolean get_create_2dmat() {return create_2dmat;}

  void set_ref_face( RefFace *ptr_ref_face ){ refFace = ptr_ref_face; } 
 
  void append_list_item( int loop_id ){ loopIndex.append(loop_id); }

  static int is_td_octree_ref_face(const ToolData* td)
      {return (CAST_TO(const_cast<ToolData*>(td), TDOctreeRefFace) != NULL);}
		
  void set_last_curve_id( int id ){ lastCurveID = id; }
    //- this is used along with loopIndex to determine non-adj curves

  DLIList<CubitPoint *> *get_ptr_cubit_point_list(){ return ptrCubitPointList; }
  void set_ptr_cubit_point_list( DLIList<CubitPoint *> *ptr_cubit_point_list ){ ptrCubitPointList = ptr_cubit_point_list; } 
 
  DLIList<CubitFacetEdge *> *get_ptr_cubit_facet_edge_list(){ return ptrCubitFacetEdgeList; }
  void set_ptr_cubit_facet_edge_list( DLIList<CubitFacetEdge *> *ptr_cubit_facet_edge_list ){ ptrCubitFacetEdgeList = ptr_cubit_facet_edge_list; } 
 
  DLIList<CubitFacet *> *get_ptr_cubit_facet_list(){ return ptrCubitFacetList; }
  void set_ptr_cubit_facet_list( DLIList<CubitFacet *> *ptr_cubit_facet_list ){ ptrCubitFacetList = ptr_cubit_facet_list; } 
 
    //- get and set pointers to DLIList of CubitPoint, CubitFacetEdge, CubitFacet 
		
  
  CubitBoolean is_adj_curves( int id1, int id2 );
    //- returns true if the two curves are adjacent

    //SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
    
    // static void set_memory_allocation_increment(int increment = 0){memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment
  
    // static void destroy_memory(){memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

  CubitFacet* split_facet_into_two_facets( CubitFacet *target_facet, CubitPoint* edge1_pt,
                                           CubitPoint* edge2_pt, 
                                           CubitFacetEdge *boundary_edge,
                                           CubitPoint* new_pt,
                                           DLIList<CubitFacet*> &facet_list,
                                           DLIList<CubitFacetEdge*> &facet_edge_list,
                                           DLIList<CubitPoint*> &point_list);

  CubitBoolean split_facet_locally_along_edge( CubitFacet *target_facet, CubitPoint* edge1_pt,
                                               CubitPoint* edge2_pt, 
                                               CubitFacetEdge *boundary_edge,
                                               int num_of_segments,
                                               DLIList<CubitFacet*> &facet_list,
                                               DLIList<CubitFacetEdge*> &facet_edge_list,
                                               DLIList<CubitPoint*> &point_list); 
  
  CubitBoolean split_facet_type_00( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_01( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_02( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_03( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_12( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_13( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_23( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
  CubitBoolean split_facet_type_33( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list );
	  

  static CubitStatus add_td( RefFace *ref_face );
  static TDOctreeRefFace* get_td( RefFace *ref_face );

  CubitBoolean check_valid_facets(CubitBoolean disable_if_bad);

};
    

#endif // TD_SKELETON_ref_FACE


//EOF

