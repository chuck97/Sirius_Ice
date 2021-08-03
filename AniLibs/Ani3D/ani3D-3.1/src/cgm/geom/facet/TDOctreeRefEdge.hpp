//- Class:       TDOctreeRefEdge
//- Description: Tool data for storing additional information related to faceting edges
//- Owner:       W. R. Quadors
//- Checked by:
//- Version:

#ifndef TD_SKELETON_MREF_EDGE

#define TD_SKELETON_MREF_EDGE

#include "CubitDefines.h"
#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "CastTo.hpp"
#include "GMem.hpp"

struct GPoint;
class RefEdge;
class CubitVector;
  
class TDOctreeRefEdge : public ToolData
{
private:
    //static MemoryManager memoryManager;
    //- memory management object

  GMem *resultPointData;
  RefEdge *refEdge;
  CubitBoolean visit;

public:

  TDOctreeRefEdge();
    //- constructor

  ~TDOctreeRefEdge();
    //- destructor
 
  inline CubitBoolean get_visit( void ){ return visit; }
  inline void set_visit( CubitBoolean type ){ visit = type; }

    // deletes all Gpoints in the list
  void clean_gpoint_list( void );

  void set_mref_edge( RefEdge *ptr_mref_edge ){ refEdge = ptr_mref_edge; } 
 
  static int is_td_skl_mref_edge(const ToolData* td){return (CAST_TO(const_cast<ToolData*>(td), TDOctreeRefEdge) != NULL);}

  GPoint *get_head_gpoint(){ return resultPointData->point_list(); }
	
  int get_gpoint_list_length( void ){ return resultPointData->pointListCount; }

    //SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
    
    // static void set_memory_allocation_increment(int increment = 0){memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment
  
    // static void destroy_memory(){memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

  static CubitStatus add_td( RefEdge *mref_edge );
  static TDOctreeRefEdge* get_td( RefEdge *mref_edge );

  static void decimate_curve_points_for_source_entity(GMem &point_data, double angle_tolerance, GMem& result_point_data);
  static CubitBoolean find_curve_curvature_using_three_points(CubitVector point_a, CubitVector point_b, CubitVector point_c, CubitVector &curvature );

  CubitBoolean generate_gpoint_list( double decimation_ang );
};
    

#endif // TD_SKELETON_MREF_EDGE


//EOF

