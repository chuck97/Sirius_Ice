//-------------------------------------------------------------------------
// Filename      : CubitOctreeCell.hpp
//
// Purpose       : A cell of the CubitOctree model is defined 	
//
// Creator       : William Roshan Quadros
//
// Creation Date :  01/01/2003
//
// Owner         : 
//-------------------------------------------------------------------------
#ifndef CUBITOCTREECELL_H
#define CUBITOCTREECELL_H 

#include "DLIList.hpp"
#include "CubitVector.hpp"

class CubitOctreeNode;
class CubitOctree;
class RefFace;
class CubitBox;
class CubitFacet;
class CubitFacetEdge; 
class CubitPoint;
class OctreeIntersectionData;
class OctreeFacetPointData;


class CubitOctreeCell
{
  
public:
  CubitOctreeCell(){}
  CubitOctreeCell( CubitVector center, double dimension, int level, CubitOctreeCell *parent_cell );
    //- Constructor
  
  virtual ~CubitOctreeCell();
    //- Distructor
  
  void display( CubitOctree *ptr_octree, int opt);
    //- Display of cells
    
  void display_color_wireframe( CubitOctree *ptr_octree);
    //- diaplay cell in wireframe by color coding size at node
    
  void display_octreefacetpointdata( void );
    //-  Displays the facet point data contained in the list

  void release_octreefacetpointdata( void );
    //- Releases the facet point data stored at cells 

  void set_oct_grid_node(const int i, const int j, const int k, CubitOctreeNode *ptr_node );
  
  inline int get_depth( void ){ return depth; }
    //- returns depth
  
  inline double get_dimension( void ){ return dimension; }
    //- returns the dimension
  
  void set_child( int i, int j, int k, CubitOctreeCell *ptr_child_cell );
    //- updates the child of the cell
  
  inline CubitVector get_center( void ){ return mCenter; }
    //- returns the center of the cell
  
  inline CubitOctreeNode *get_octree_grid_node( const int i, const int j, const int k ){ return cubitOctreeNode[i][j][k]; }
    //- returns the node of the cell
  
  inline int id( void ){ return num; }
    //- returns the numbers
  
  CubitOctreeCell * get_child( const int i, const int j, const int k );
    //- returns the child at i,j,k
  
  inline void set_leaf( bool type ){ leaf = type; }
  
  CubitOctreeCell * find_leaf_octree_cell(const CubitVector &point );
    //- Find the leaf cell 
  
  void subdivide_octree_based_on_point( CubitPoint *ptr_point );
    //- subdivides the CubitOctree based on point ( either original point or centroid )
  
  bool is_leaf( void ){ return leaf; }
    //- returns the status of leaf 
  
  void distribute_facet_points_among_children( void );
    // distributes the facet points and associated face among the children
  
  int num_of_facet_point_data( void ){ return octreeFacetPointDataList.size(); }
    //- returns number of facet points and centroids
  
  CubitBoolean append_list_item( OctreeFacetPointData *ptr_facet_point_data );
    //- appends the arguments into the respective list

  inline DLIList<OctreeFacetPointData*>* get_facet_point_data_list() 
      {
        return &octreeFacetPointDataList;
      }
  
  
  inline void set_mark( CubitBoolean type ){ mark = type; }
  inline CubitBoolean get_mark( void ){ return mark; }

  inline void set_visit( CubitBoolean type ){ visit = type; }
  inline CubitBoolean get_visit( void ){ return visit; }
  
  CubitBoolean add_adjacent_unmarked_cells( DLIList< CubitOctreeCell *> &queue );
  
  CubitBoolean is_intersects_box(const CubitBox &box );
  
  CubitBoolean does_facet_intersect_octreecell( CubitFacet *ptr_facet );


  void set_color_and_intersection_datas(CubitFacet *ptr_facet, RefFace *ptr_Ref_face,
#ifdef USE_octree_BGMESH
                                        DLIList<CubitOctreeCell*> *greyCellList,
#endif
                                        CubitSense surf_sense=CUBIT_FORWARD);
    //- if cell intersects supplied facet call this function to color cell grey and attached OctreeIntersectionDatas to its nodes
    //- function assumes that grid nodes of cell still have halfspace settings corresponding to the supplied facet
  
    // The function assumes that nodes have half space information
  CubitBoolean does_contain_positive_and_negative_nodes();

  CubitBoolean is_facet_point_data_present( const CubitVector &coord );
  CubitBoolean is_facet_point_data_present( OctreeFacetPointData *new_facet_point_data );
    //- Currently CubitVector is passed and distance between existing points and coord is compared
    //- with SKL_EPSILON in determining if the coord already exists.  
    //- To avoid SKL_EPSILON CubitPoint srl number can be used.  This function is currently not used.

  CubitBoolean interpolate_grey_octreecell_node( void );
    //- Finds size at WHITE nodes of GREY cells
  
  void coloring( /*DLIList<CubitOctreeCell *> &grey_cell_list,*/ DLIList<CubitOctreeCell*> &black_cell_list);
    //- Colors the cells with WHITE (cell outside solid), BLACK (cell inside solid), and GREY (cell at the boundary)

    // Trilienar interpolation	
    // if number of nodes with zero size exceeds  MAX_NUM_ZERO_SIZE_NODES
    // then average size is used instead of zero
    // Even GREY grid nodes are consided 
    // If a node has zero size then that node is not considered
    // Therefore color is not important but the size
  double trilinear_interpolation(const CubitVector &point );
  
  double inverse_distance_interpolation(const CubitVector &point );
  double min_size_interpolation(const CubitVector &point );
  double min_distance_interpolation(const CubitVector &point );
  
  CubitStatus find_indices_in_parent( int *index );

#ifndef NDEBUG
  void write_octreecell_sizing_info_file( FILE *pof, DLIList<CubitOctreeCell*> &stack );
#endif

    /*
      void set_node( CubitOctreeNode *r_t_n_, CubitOctreeNode *r_t_f_, CubitOctreeNode *r_b_n_, CubitOctreeNode *r_b_f_, CubitOctreeNode *l_t_n_, CubitOctreeNode *l_t_f_, CubitOctreeNode *l_b_n_, CubitOctreeNode *l_b_f_ );
        //- Sets the eight nodes of the cell
    
        void set_octreecell( CubitOctreeCell *left_cell, CubitOctreeCell *right_cell, CubitOctreeCell *bottom_cell, CubitOctreeCell *top_cell, CubitOctreeCell *far_cell, CubitOctreeCell *near_cell );
          //- Sets the six adjacent cells
          */    
  inline void mark_color( int type ){ color = type; }
    //- Cell contained completely inside the solid is call BLACK
    //- Cell contained completely outside the solid is calle WHITE
    //- Cell at the boundary with some nodes inside and others outside
    //	are called GREY

  inline int get_color() {return color;}
    //- Returns the color of cell

private:
  
  CubitVector mCenter;
    //- Center of cell
  
  double dimension;
    //- size of the CubitOctree cell
  
  int depth;
    //- depth of the node
  
  CubitBoolean leaf;
    //- to check leaf 
  
  int num;
    //- Serial number of cell
  
  int color;	
    //- WHITE (outside solid), BLACK (inside solid), and GREY (on boundary)
  
  CubitBoolean mark;
    //- Mark used for visiting cell;

  CubitBoolean visit;
  
    //CubitOctreeNode* (*CubitOctreeNode)[2][2];
  CubitOctreeNode* cubitOctreeNode[2][2][2];
  
    //- Eight corner nodes of cell
  
    //CubitOctree *CubitOctree;
    //- CubitOctree is currently not used
  
  CubitOctreeCell *parent;
  CubitOctreeCell *children[2][2][2];
    //- eight children
  
  DLIList< OctreeFacetPointData *> octreeFacetPointDataList;
    //- Stores the facet points 
  
  DLIList<OctreeIntersectionData*> *myFacetList;
  
};

#endif

//EOF

