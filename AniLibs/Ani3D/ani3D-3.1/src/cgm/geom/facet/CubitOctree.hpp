//-------------------------------------------------------------------------
// Filename      : CubitOctree.hpp
//CubitOctree
// Purpose       : CubitOctree model of the input solids/surfaces is generated. 
//                 
//
// Creator       : William Roshan Quadros
//
// Creation Date :  09/01/2013
//
// Owner         : 
//-------------------------------------------------------------------------
#ifndef CUBITOCTREE_H
#define CUBITOCTREE_H

#include "CubitVector.hpp"
#include "CubitOctreeConstants.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"


class RefVolume;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;


class CubitOctreeNode;
class CubitOctreeCell;
class FaceAdjacency;
class VertexFaceIncident;
class CubitFacet;
class CubitFacetEdge;
class FacetEntity;
class CubitPoint;
class CubitBox;
class TDSkeletonCubitPoint;
class CubitOctreeGenerator;
class OctreeFacetPointData;


// Generation of CubitOctree model of the CAD model:
// Graphics Facets are used to detect curvature changes in surface and curves.
// Vertices and Centroid of the Facets are used to generate the CubitOctree.
// Bounding Box of the facets, Edge Length etc can be used.
// Each Cell is made to contain only one point.

class CubitOctree
{
  
public:
  	 
  CubitOctree( CubitOctreeGenerator *lattice_gen );
    //- Constructor
  
  virtual ~CubitOctree();
    //- Distructor

  void release_memory( void );
    //- clears the memory for assembly

  inline int get_max_depth( void){ return maxDepth;}
  inline int get_min_depth( void){ return minDepth;}
  inline void set_max_depth( int max_depth ) { maxDepth = max_depth; }
  inline void set_min_depth( int min_depth ) { minDepth = min_depth; }
  
  inline DLIList<CubitOctreeNode*> *get_whitenodelist( void ){ return &boundaryWhiteGridNodeList; }

  inline DLIList<CubitOctreeCell*> *get_greycelllist( void ){ return &greyCellList; }
  
  inline DLIList<CubitOctreeCell*> *get_blackcelllist() {return &blackCellList;}

  inline DLIList<CubitOctreeNode*> *get_gridnodevector( void ){ return &gridNodeVector; }

  inline void set_status_octree_coloring( CubitBoolean type ){ statusCubitOctreeColoring = type; }
  
  // Methods related to CubitOctree Generation 
  inline CubitOctreeCell *get_root( void ){ return root; }

  void color_octreecell( void );
  //- NOTE: during CubitOctree generation only the points at the facets are marked white and black
  //- during mat generation whole interior points are colored BLACK and the later the cells' are colored
  //- Identifies boundary (GREY), interior (BLACK) and exterior cells (WHITE)
  
  
  CubitOctreeCell *find_octreecell( const CubitVector &pnt );
  //- Returns a cell containing point pnt
  
  CubitBoolean initialize_octree_generation( void );
  //- Initializes the root cell which is equal to bbox
  
  CubitBoolean build_octree_till_min_depth( CubitOctreeCell *ptr_cell );
  //- Builds the CubitOctree till the minimum depth
  
  CubitBoolean subdivide_octree_cell( CubitOctreeCell *ptr_cell, DLIList<CubitOctreeNode*> *ptr_queue, int max_depth);
  //- subdivides the cell
  
  CubitBoolean subdivide_octree_based_on_facet_point( OctreeFacetPointData *ptr_facet_point_data, int max_depth );
  //- use the facet center and vertices to subdivide the CubitOctree
  
  double size_at_a_point( const CubitVector &point, int size_type );

  CubitPointContainment point_containment( CubitVector tmp_vec, double tolerance );
  
  
  inline double get_max_size_grid_node( void ){ return maxSizeGridNode; }
  //- Returns maximum size
  
  inline double get_min_size_grid_node( void ) const{ return minSizeGridNode; }
  //- Returns minimum size
  
  void display( const int mode, int opt = 0 );
  //- Display of CubitOctree
  
  
  // For better 3D medial resolution near adjacent surfaces with a large dihedral angle (outward normals) between them,
  // CubitOctree cells intersecting with the shared curve can be refined to the maximum allowable depth.
  // First call find_cells... to build the list of cells to refine, then call refine_cells_to_target_depth on them
  CubitBoolean find_cells_based_on_surface_angle(DLIList<CubitOctreeCell*> &refine_cell_list, DLIList<CubitFacetEdge*> &crease_edges,
                                                 DLIList<CubitFacet*> &crease_facets, RefFace *one, RefFace *two,
                                                 const CubitBoolean draw_facets=CUBIT_FALSE, double dihedral_angle_thresh=CUBIT_PI/2.0);
  
  
  CubitBoolean find_octree_cells_contained_inside_bbox( CubitFacet *ptr_facet, DLIList<CubitOctreeCell*> &CubitOctree_cell_list );
  
  //- output: CubitOctree cells which intersect with the bbox are stored in CubitOctree_cell_list
  
  CubitBoolean mark_positive_and_negative_octree_grid_nodes( CubitFacet *ptr_facet, DLIList<CubitOctreeCell*> &CubitOctree_cell_list, DLIList< CubitOctreeNode *> &CubitOctree_grid_node_list );
  
  //- output: the nodes of cells that intersect with the bbox are stored in octee_cell_list
  
  CubitBoolean find_intersection_between_grid_edges_and_facet( CubitOctreeType type, RefFace *ptr_face, CubitFacet *ptr_facet, DLIList<CubitOctreeNode *> &CubitOctree_grid_node_list  );
  //- finds intersection between line segment (edge ) incident at every grid_node contained inside bbox and triangle ( facet)
  
  
  void find_cell_list_intersecting_line_segment(const CubitVector &p0, const CubitVector &p1, DLIList<CubitOctreeCell*> &cell_list);
  //- returns the list of CubitOctree cells intersecting with line segment defined by end points
  
  double get_scaled_from_wrld_size_grid_node( double wrld_size );

  
  static double capsule_distance_to_facet(const CubitVector &point, CubitFacet *lp_facet, CubitVector &int_point, CubitBoolean use_projection_only = CUBIT_FALSE);
  
  CubitBoolean establish_smooth_transition_of_cells( int max_depth ); 

private:
  
  inline CubitBoolean apd_vector_item( CubitOctreeNode *ptr_grid_node );
    //- appends the ptr_grid_node in the CubitOctreeNodeList
  
  double calculate_depth_based_on_size( double size );

  inline double get_bbox_dimension( void ){ return bboxDimension; }

  
  
  void refine_cells_to_target_depth(DLIList<CubitOctreeCell*> &refine_cell_list, const int target_depth);

    // To optimize intersection between facet and CubitOctree cells, all the cells not intersecting with the plance containing facet are removed
    // The nodes are marked with +ve or -ve half space info before calling this function
  CubitBoolean mark_cells_that_intersect_facet_plane( /*CubitFacet *ptr_facet,*/ DLIList<CubitOctreeCell*> &CubitOctree_cell_list );
  
    // Stores the CubitOctree cells intersecting with the facet
  CubitBoolean unmark_octree_cells_not_intersecting_with_facet( CubitFacet *ptr_facet, DLIList<CubitOctreeCell *> &CubitOctree_cell_list );

    //- CubitOctree is subdivided till skeleton max depth if size or scope of source entity is less than the dimention of
    //- the cell containing it.

  CubitBoolean subdivide_octree_from_leaf_cell( CubitOctreeCell *ptr_cell, OctreeFacetPointData *ptr_facet_point_data, int max_depth );
 
    //- subdivides the CubitOctree starting from ptr_cell such that the ptr_new_point is contained in a uniue cell
    //- if the maxDepth is not reached.

  void gather_initial_grid_nodes_for_smooth_transition( DLIList<CubitOctreeNode *> &queue );
    //- queue hold potential nodes which has more than or equal to 2 cell depth difference
 
  
    //- establishes smooth transition by keeping adjacent cells differing by 1 level
  
    //CubitBoolean find_intersection_between_octree_and_facets( DLIList<CubitOctreeNode *> &queue_for_mat_generation );
    //- Finds the intersection between the CubitOctree line segments and the facets

  void reset_internal_node_size( float value, int type );
    //- Resets the node size with value  
	
  void find_radius_at_boundary_node( void );
	
 

  void find_max_min_size_grid_node_for_scaling( void );
    //- Finds maximum and minimum size used for scaling.


  
   static CubitStatus circumcenter(CubitVector &a, CubitVector &b, CubitVector &c, CubitVector &center, double &radius);
    /*  
        void build_boundary_white_node_list( void );
          //- Uses grey CubitOctrees to find white boundary grid nodes;
          */


  void prepare_for_background_mesh();
    // performs necessary actions on CubitOctree to make way for background mesh

#ifndef NDEBUG 

  void write_sizing_info_file( const char *file_name );
    //- Writes sizing function at grid nodes to a file

  void write_octree_sizing_info_file( const char *file_name );
    //- Writes sizing function stored in the CubitOctree on to a file.

  void write_matlab_sizing_info_file( const char *file_name );
    //- Writes sizing function stored in the CubitOctree for visualization in matlab

#endif
  
private:
  	
  CubitBoolean statusCubitOctreeColoring;
  

    //- Maximum and minimum radius of source points and grid node
    //- this is used in scaling the source points between 1 and 10 for display purpose
		
  double bboxMinY; // minY
  double bboxMaxY; // maxY
  double bboxMaxZ; // maxZ
  double bboxMinZ; // minZ 
  double bboxMaxX; // maxX 
  double bboxMinX; // minX
    //- Limits of the bounding box of the volume
 
  CubitOctreeCell *root;
    //- Root of the CubitOctree
  
  CubitVector bboxCenter;
    //- Center of the bounding box
  
  double bboxDimension;
    //- Dimension of the bounding box
  
  double epsilonBBox;
    //- Actual bounding box is enlarged by epsilon_bbox in all directions
   
  DLIList<CubitOctreeCell *> greyCellList;
    //- List of boundary CubitOctrees

  DLIList<CubitOctreeCell*> blackCellList;
    //- List of internal CubitOctree cells
  
  DLIList<CubitOctreeNode *> boundaryWhiteGridNodeList;
    //- List of white nodes of grey cells
  
  DLIList <CubitOctreeNode *> gridNodeVector;
    //- Grid Nodes of CubitOctree

  DLIList <CubitOctreeNode*> greyGridNodeVector;
  
    
  CubitOctreeGenerator *cubitOctreeGenerator;


  double maxSizeGridNode;
  double minSizeGridNode;
  
  int minDepth;
  int maxDepth;

};

#endif
    
//EOF


