//-------------------------------------------------------------------------
// Filename      : CubitOctree.hpp
//
// Purpose       : Defines the grid node of the CubitOctree model	
//
// Creator       : William Roshan Quadros
//
// Creation Date :  01/01/2003
//
// Owner         : 
//-------------------------------------------------------------------------
#ifndef CUBITOCTREENODE_H
#define CUBITOCTREENODE_H
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "PriorityQueue.hpp"
#include <queue>
#include "CubitDefines.h"
//#include "MemoryManager.hpp"
#include "OctreeIntersectionData.hpp"
#include "CubitOctreeConstants.hpp"

class RefFace;
class CubitOctree;
class CubitOctreeCell;
class CubitFacet;
class CubitOctreeNode;
class CubitOctreeGenerator;


class CubitOctreeNode{
  
public:
  
  CubitOctreeNode( const CubitVector &cen, CubitOctreeCell *parent_cell, const int i, const int j, const int k );
  CubitOctreeNode( const double &x, const double &y, const double &z );
  void initialize_constructor( const double &x, const double &y, const double &z );
    //- Constructor
  
  ~CubitOctreeNode();
    //- Distructor
  
  void display( OctreeNodeConstant type = NODE_SIZE, float draw_size = -1 ); 
    //- Displays the grid nodes
  
  inline int get_color(){ return color; }
    //- Returns color

  inline void set_color (const int new_color) {color = new_color;}
  
  inline int get_num(){ return num; }
    //- Returns the serial number of the grid node

  inline CubitVector& get_coord(){ return coord; }
    //- Returns the cartesian coordinates
  
  double get_size( OctreeSourceEntityType type ) const;
    //-	Returns size of the node

  int id(){ return num; }

  void set_size( double s, int type );
    //-	Sets size at the node
  
  void set_adj_node( enum OctreePosition type, CubitOctreeNode *ptr_grid_node );
    //- sets the adjacent nodes
  
  void set_adj_node_distance( enum OctreePosition type, int dist );
    //- sets the adjacent nodes distance
  
  double manhattan_distance_adj_node( int index );
  inline double manhattan_distance_adj_node( CubitOctreeNode *ptr_adj_node );  

  inline double get_distance(){ return distance; }
  
  inline double x(){ return coord.x(); }
  inline double y(){ return coord.y(); }
  inline double z(){ return coord.z(); }
    //- Returns cartesian co-ordinates 
  
  inline void x( double value ){ coord.x( value ); }
  inline void y( double value ){ coord.y( value ); }
  inline void z( double value ){ coord.z( value ); }
    //- Sets cartician co-ordinates of node
  
  CubitOctreeNode *get_adj_node( int num );
    //- Returns one of the six adjacent nodes
  
  int get_adj_node_distance( enum OctreePosition type );
    //- returns the distance with the adj node in the direction of type
  
  void append_list_item( CubitVector facet_normal );
  
   void calc_facet_patch_distance_normal(DLIList<OctreeIntersectionData*> &idatas, int num_use_idatas, double &patch_distance, CubitVector &patch_normal, CubitBoolean sort, CubitBoolean set_Refface);
  void SAT_find_face_distance_average_normal ();

  inline CubitOctreeCell* get_min_depth_cell(){ return minDepthCell; }
    //- returns minimum depth cell

  inline int get_cell_depth_difference(){ return cellDepthDifference; }
  

  int get_counter();
    //- uses static counter to number the cell

  static inline void reset_counter() {mCounter = -1;}

  inline void update_adj_cell( CubitOctreeCell *ptr_cell, int i, int j, int k ) 
      {
        adjCell[i][j][k] = ptr_cell;
      }
    //- updates the cell at i, j, k location
  
  int find_min_depth_cell_and_depth_difference( void  );
    //- finds the maximum depth adjacent cell and difference in adj cells' depth
  
  void calculate_size_based_on_cell_dimension( double bbox_dimension );
    //- calculates size using 6 adj nodes distance
  
  inline CubitOctreeCell *get_adj_cell( const int i, const int j, const int k ) { return adjCell[i][j][k]; }
  
    // MARK is used in establishing connection between skl point and grid nodes contained in it using DFS
  inline void set_mark( CubitBoolean type ){ mark = type; }
  inline CubitBoolean get_mark(){ return mark; }
  
  inline void set_visit( CubitBoolean type ){ visit = type;	}
    //- sets the visit with type ( CUBIT_FALSE or CUBIT_TRUE )
    
  inline CubitBoolean get_visit() const {return visit;}
    //- Returns visit
    
  inline CubitVector get_normal() const {return mNormal;}


  
  int find_half_space( CubitFacet *ptr_facet );

    //- returns half space direction
  
  inline CubitBoolean get_halfspace_direction(){ return halfspaceDirection; }
  inline void set_halfspace_direction( const CubitBoolean type ){ halfspaceDirection = type ; }

  void find_intersection_with_facet( CubitOctreeType type, RefFace *ptr_face, CubitFacet *ptr_facet,  DLIList<CubitOctreeNode*> &boundary_white_node_list);
    //- checks intesection between the lines joining grid node and adjacent nodes with the facet.
    //- queue_mat_generation stores the BLACK grinodes that contains distanace and normal for MAT generation.
  
  CubitBoolean find_intersection_point( int axis, CubitVector grid_node0, CubitVector grid_node1, CubitVector &facet_normal, CubitVector facet_vert0, CubitVector facet_vert1, CubitVector facet_vert2, CubitVector &int_point, double &para );
    //- returns true if the intersection between the line segment and facet takes place
    //- para stores the parameter of intersection point int_point

  static CubitBoolean is_same_side(const CubitVector &p1, const CubitVector &p2, const CubitVector &a, const CubitVector &b);
  static CubitBoolean is_intersection_point_contained_inside_facet( const CubitVector &int_point, const CubitVector &facet_vert0, const CubitVector &facet_vert1, const CubitVector &facet_vert2 );
    //- returns true if intersection point is contained inside the facet
    //- interior angle at the intersecton point should add up to 360 deg.
  
     
    //- Appends the source point

  inline void append_list_item( OctreeIntersectionData *ptr_int_data ) 
      {
        octreeIntersectionDataList.push( ptr_int_data );
      }
  
    //- Appends intesection data like, face, int_point, facet etc.

  DLIList<OctreeIntersectionData*>* get_idata_list() {return &octreeIntersectionDataList;}
  
  CubitBoolean bfs_to_find_contained_nodes(int max_onode_per_src);
  CubitBoolean dfs_to_find_contained_nodes(int max_onode_per_src);
    // Breadth first search is used to find nodes contained inside the box of skeleton source point  

    
    //- sizing function due to each type of source point is evaluated separately.
    //- then these function can be blended to find final sizing function.
    //- curently the minimum of all function is selected.

  static CubitBoolean compare_function( CubitOctreeNode *&a, CubitOctreeNode *&b );

  void find_distance_at_adj_node(PriorityQueue<CubitOctreeNode *> *heap );
  
  CubitBoolean find_size_using_adj_node();
  
  
    //static void set_scope_box_limits( double max_x, double min_x, double max_y, double min_y, double max_z, double min_z);
    //static void reset_scope_box_limits();
    //- To avoid the overflow of stack during linking source point and grid nodes
    //static double scpBoxMaxX, scpBoxMinX, scpBoxMaxY, scpBoxMinY, scpBoxMaxZ, scpBoxMinZ;
  //SetDynamicMemoryAllocation(memoryManager)
  //- class specific new and delete operators
  
private:
  

  
  
private:
  int num;
    //- Serial number of node

  static int mCounter;
  
  CubitBoolean visit;
  CubitBoolean mark;
  
  short color; 
    // in CubitOctree based approach
    //- GREY is default ( all nodes outside WHITE or Solid )
    //- BLACK (node inside solid), WHITE (node at the boundary during facet CubitOctree intersection)

  CubitBoolean halfspaceDirection;
  
  double size;

    //- Size at the grid node is initially set to 0.0

  
   // OPT: keep pointer to list and clear list after intersection
  DLIList<OctreeIntersectionData *> octreeIntersectionDataList;
  
  // OPT: use cubit hash to temp hold this data
  int cellDepthDifference;
    //- used in splitting CubitOctree cells during smooth transition
  
  
  CubitOctreeCell *minDepthCell;
    //- keeps the pointer to minimum depth adjacent cell
    
  double distance;
  
  CubitVector coord; 
  // Co-ordinates of the grid node
  
  CubitOctreeNode *adjGridNode[6]; 
    //0 to 5 LEFT, RIGHT, BOTTOM, TOP, BACK, FRONT
  
  int adjNodeDistance[6];  ///???
    //- stores the lengths between adjacent nodes
  
  CubitOctreeCell *adjCell[2][2][2];
  
 
  RefFace *refFace;
    //- Face associated with the grid node used in wave propagation
  
  CubitVector mNormal;
    //- Normal at the grid_node is used in MAT generation.

  
  //static MemoryManager memoryManager;
    //- memory management object

};

#endif

//EOF

