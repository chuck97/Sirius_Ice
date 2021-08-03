#include "CubitOctree.hpp"
#include "RefFace.hpp" 
#include "CubitOctreeConstants.hpp"

#include "CubitOctreeCell.hpp"
#include "CubitOctreeNode.hpp"  
#include "CubitFacet.hpp" 
#include "CubitFacetEdge.hpp" 
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "TDOctreeRefFace.hpp"
#include "FacetDataUtil.hpp"
#include "CubitOctreeGenerator.hpp"
#include "OctreeFacetPointData.hpp"



/* -------------------- Methods of CubitOctree ---------------- */
//      Takes the following parameters
//		min_depth : Minimum depth of the CubitOctree (avoids large cubes).
//      max_depht : Maximum depth of the CubitOctree (restricted by memory). 
	
CubitOctree::CubitOctree( CubitOctreeGenerator *lattice_gen ){

  CubitOctreeNode::reset_counter();
  
  cubitOctreeGenerator = lattice_gen;


  maxSizeGridNode = -CUBIT_DBL_MAX; 
  minSizeGridNode = CUBIT_DBL_MAX;

  CubitVector min_corner; 
  CubitVector max_corner;  
  
  lattice_gen->get_bounding_box( min_corner, max_corner );

    // note (vvyas 7/2005): perturbing bbox unsymmetrically to help avoid degenerate cases
    // in SAT intersection code
  CubitVector min_perturbation(-.0001,-.0002,-.0003);
  CubitVector max_perturbation(.0005,.00012,.00028);

  min_corner += min_perturbation;
  max_corner += max_perturbation;
  
  CubitBox bbox(min_corner, max_corner );
  
  double x_range = bbox.x_range();
  double y_range = bbox.y_range();
  double z_range = bbox.z_range();
  
  double max_range;
  max_range = ( x_range > y_range ) ? x_range : y_range;
  max_range = ( z_range > max_range ) ? z_range : max_range;
  
  bboxCenter.x(( min_corner.x() + max_corner.x() ) / 2.0);
  bboxCenter.y(( min_corner.y() + max_corner.y() ) / 2.0);
  bboxCenter.z(( min_corner.z() + max_corner.z() ) / 2.0);
  
  epsilonBBox = max_range / pow(2.0, OCTREE_MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL] ) * 0.15;
  bboxDimension = max_range + 2 * epsilonBBox;
  
  
  bboxMaxX = bboxCenter.x() + bboxDimension / 2.0;
  bboxMinX = bboxCenter.x() - bboxDimension / 2.0;
  bboxMaxY = bboxCenter.y() + bboxDimension / 2.0;
  bboxMinY = bboxCenter.y() - bboxDimension / 2.0;
  bboxMaxZ = bboxCenter.z() + bboxDimension / 2.0;
  bboxMinZ = bboxCenter.z() - bboxDimension / 2.0;
  
    // Data structure flags are used to know if a perticualr data structure is already established or not. 

  statusCubitOctreeColoring = CUBIT_FALSE;
  
  minDepth = OCTREE_MIN_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL]; 
  maxDepth = OCTREE_MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL]; 
  
  root = NULL;
}

CubitOctree::~CubitOctree(){
  
  greyCellList.clean_out();
  
  boundaryWhiteGridNodeList.clean_out();

  int i;
  int num_elm;
  
    //PRINT_INFO("BEFORE: Deleted grid nodes \n");
    //system("pause");
  num_elm = gridNodeVector.size();
  CubitOctreeNode *ognode;
  for( i = 0; i < num_elm; i++ ){
    ognode = gridNodeVector.pop();
    if (ognode) {delete ognode;}
  }
  for (i=0; i < greyGridNodeVector.size(); ++i)
  {
    delete (greyGridNodeVector.get_and_step());
  }
  
    //PRINT_INFO("AFTER: Deleted grid nodes \n");
    //system("pause");
  delete root;

    //PRINT_INFO("AFTER: Deleted CubitOctree cells \n");
    // system("pause");
}


double CubitOctree::calculate_depth_based_on_size( double size ){

  return ceil(  log((double)get_bbox_dimension() / (double)size ) / log(2.0)  );
}



CubitBoolean CubitOctree::initialize_octree_generation( void ){
  
  
    // Normalize the center and size of the CubitOctree cell
    // bbox dimension is 1.0 with center at 1/2, 1/2, 1/2
  CubitVector norm_center;
  int level = 0;
  int i,j, k, l, m, n;
  double norm_dim;
  
    //norm_center.x(100.0 / 2);
    //norm_center.y(100.0 / 2);
    //norm_center.z(100.0 / 2);
    //norm_dim = 100.0;
  
  norm_center = bboxCenter;
  norm_dim = bboxDimension;
  
    //VERBOSE PRINT_INFO(" BBOX Center = [ %f %f %f ] Dim = %f \n",norm_center.x(), norm_center.y(), norm_center.z(), norm_dim );
  
  root = new CubitOctreeCell( norm_center, norm_dim, level, NULL );
  CubitVector shift;
  
  
  
  CubitOctreeNode *ptr_oct_node[2][2][2];
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if( i == 0 )
          shift.x( -norm_dim / 2.0 );
        else
          shift.x( norm_dim / 2.0 );
	      
        if( j == 0 )
          shift.y( -norm_dim / 2.0 );
        else
          shift.y( norm_dim / 2.0 );		
	      
        if( k == 0 )
          shift.z( -norm_dim / 2.0 );
        else
          shift.z( norm_dim / 2.0 );
	      
        l = ( i == 1 )?0:1;  // fliped 0 -> 1;  1 -> 0
        m = ( j == 1 )?0:1;
        n = ( k == 1 )?0:1;
	      
        ptr_oct_node[i][j][k] = new CubitOctreeNode( norm_center + shift, root, l, m, n ); 
        apd_vector_item( ptr_oct_node[i][j][k] );
      }
    }
  }
  
    // Connecting CubitOctree Grid Nodes
    //-	adj_node[0] = O_FRONT Node
    //	adj_node[1] = O_BACK Node
    //	adj_node[2] = O_RIGHT Node
    //	adj_node[3]	= O_LEFT Node
    //	adj_node[4] = O_TOP Node
    //	adj_node[5] = O_BOTTOM Node
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
	
        if( i == 0 ){
          ptr_oct_node[i][j][k]->set_adj_node( O_FRONT, ptr_oct_node[1][j][k] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_FRONT, 0 );
        }
        else{
          ptr_oct_node[i][j][k]->set_adj_node( O_BACK, ptr_oct_node[0][j][k] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_BACK, 0 );
        }
	      
        if( j == 0 ){
          ptr_oct_node[i][j][k]->set_adj_node( O_RIGHT, ptr_oct_node[i][1][k] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_RIGHT, 0 );
        }
        else{
          ptr_oct_node[i][j][k]->set_adj_node( O_LEFT, ptr_oct_node[i][0][k] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_LEFT, 0 );
        }
	      
        if( k == 0 ){
          ptr_oct_node[i][j][k]->set_adj_node( O_TOP, ptr_oct_node[i][j][1] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_TOP, 0 );
        }
        else{
          ptr_oct_node[i][j][k]->set_adj_node( O_BOTTOM, ptr_oct_node[i][j][0] );
          ptr_oct_node[i][j][k]->set_adj_node_distance( O_BOTTOM, 0 );
        }
	
      }
    }
  }
  
  
    // Update the nodes of the root
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        root->set_oct_grid_node( i, j, k, ptr_oct_node[i][j][k] );
	
      }
    }
  }
  
    //root->display( this );
  
  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::build_octree_till_min_depth( CubitOctreeCell *ptr_cell ){
  
  int i, j, k;
  
  if( ptr_cell->get_depth() < minDepth ){
    
    subdivide_octree_cell( ptr_cell, NULL, minDepth );
    
      // PRINT_INFO(" establish connection between new nodes and pointed old nodes \n");
    
      // Recursively call all the childen		
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
            // PRINT_INFO("Child : %d %d %d \n",i,j,k );
          build_octree_till_min_depth( ptr_cell->get_child(i, j, k) );
        }
      }
    }

      // PRINT_INFO(" finished calling all recursively all children \n");
  
  } 
  return CUBIT_TRUE;
  
}



CubitBoolean CubitOctree::subdivide_octree_based_on_facet_point( OctreeFacetPointData *ptr_facet_point_data, int max_depth  ){
  
  CubitOctreeCell *ptr_leaf_cell;

  ptr_leaf_cell = root->find_leaf_octree_cell( ptr_facet_point_data->coordinates() );
  
  CubitBoolean result = subdivide_octree_from_leaf_cell( ptr_leaf_cell, ptr_facet_point_data, max_depth );
  if( result == CUBIT_FALSE ){
      // Simillar point exist; delete this facet_point_data
    delete ptr_facet_point_data;
    return CUBIT_FALSE;
  }

  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::find_cells_based_on_surface_angle(DLIList<CubitOctreeCell*> &refine_cell_list, DLIList<CubitFacetEdge*> &crease_edges, DLIList<CubitFacet*> &crease_facets,
                                                       RefFace *one, RefFace *two, const CubitBoolean draw_facets, double dihedral_angle_thresh)
{
  int i, j, k;
  
  CubitFacetEdge *edge;
  CubitBoolean water_tight = CUBIT_FALSE;
  
  DLIList<DLIList<CubitFacet*>*> shell_list;
  DLIList<CubitQuadFacet*> dummy;
  DLIList<CubitFacet*> mark1, mark2, facets, facet1, facet2, draw_facet_list;
  DLIList<CubitFacetEdge*> new_edge1, new_edge2;
  DLIList<CubitPoint*> new_points;

  CubitFacet *temp_facet = NULL;
  
  if (one == NULL || two == NULL)
  {
    return CUBIT_FALSE;
  }

  dihedral_angle_thresh = cos(dihedral_angle_thresh);

  TDOctreeRefFace *td_skl_one = TDOctreeRefFace::get_td(one),
      *td_skl_two = TDOctreeRefFace::get_td(two);

  if (td_skl_one == NULL || td_skl_two == NULL ||
      td_skl_one->get_ptr_cubit_facet_list() == NULL || td_skl_two->get_ptr_cubit_facet_list() == NULL)
  {
    PRINT_DEBUG_157("CubitOctree::find_cells_based_on_surface_angle: null TDs or null facet list ptr\n");
    PRINT_DEBUG_157("RefFace one->id()=%d, td_skl=%p\ntwo->id()=%d, td_skl=%p\n", one->id(), (void*)td_skl_one,
                    two->id(), (void*)td_skl_two);
    PRINT_DEBUG_157("This may happen when one volume sizing function \'steals\'\na surface from an existing sizing function.");
    return CUBIT_FALSE;
  }

  if (!td_skl_one->get_create_2dmat() || !td_skl_two->get_create_2dmat()) {return CUBIT_FALSE;}
  

  FacetDataUtil::copy_facets(*TDOctreeRefFace::get_td(one)->get_ptr_cubit_facet_list(), facet1, new_points, new_edge1);
  FacetDataUtil::copy_facets(*TDOctreeRefFace::get_td(two)->get_ptr_cubit_facet_list(), facet2, new_points, new_edge2);

  for (i=0; i < new_edge1.size(); ++i)
  {
    edge = new_edge1.get_and_step();
    if (edge->number_tris() == 1)
    {
      temp_facet = edge->adj_facet(0);
      if (!temp_facet->marked())
      {
        mark1.append(edge->adj_facet(0));
        temp_facet->marked(CUBIT_TRUE);
      }
        //facets.append(edge->adj_facet(0));
    }
  }

  for (i=0; i < new_edge2.size(); ++i)
  {
    edge = new_edge2.get_and_step();
    if (edge->number_tris() == 1)
    {
      temp_facet = edge->adj_facet(0);
      if (!temp_facet->marked())
      {
        mark2.append(edge->adj_facet(0));
        temp_facet->marked(CUBIT_TRUE);
      }
        //facets.append(edge->adj_facet(0));
    }
  }

  for (i=0; i < mark1.size(); ++i)
  {
    mark1.get_and_step()->marked(CUBIT_FALSE);
  }
  
  for (i=0; i < mark2.size(); ++i)
  {
    mark2.get_and_step()->marked(CUBIT_FALSE);
  }
  
  facets += facet1;
  facets += facet2;
  
  facet1.clean_out(); //delete facet1;
  facet2.clean_out(); //delete facet2;
  new_edge1.clean_out(); //delete new_edge1;
  new_edge2.clean_out(); //delete new_edge2;
  new_points.clean_out(); //delete new_points;
  

  CubitStatus rv = FacetDataUtil::split_into_shells(/*facet1*/facets, dummy, shell_list, water_tight);
  if (!rv || shell_list.size() == 0)
  {
    PRINT_DEBUG_157("CubitOctree::find_cells_based_on_surface_angle: could not form shells from facet list!\n");
    PRINT_DEBUG_157("return value was %d, number of shells is %d\n", rv, shell_list.size());
      //SVDrawTool::clear_non_retained();
      //SVDrawTool::draw_facets(facets, CUBIT_RED_INDEX);
      //SVDrawTool::mouse_xforms();
    return CUBIT_FALSE;
  }

      
    /*
  int num_null_shells = 0;
  for (i=0; i < shell_list.size(); ++i)
  {
    if (shell_list.get_and_step() == NULL)
    {
      ++num_null_shells;
    }
  }
  if (num_null_shells)
  {
    PRINT_INFO("CubitOctree::find_cells_based_on_surface_angle: %d of %d shells formed from facets are NULL!\n", num_null_shells, shell_list.size());
  }
    */
      
////////////////  
  rv = FacetDataUtil::stitch_facets(shell_list, 1e-4, water_tight, CUBIT_FALSE);
  if (!rv)
  {
    PRINT_DEBUG_157("CubitOctree::find_cells_based_on_surface_angle: could not stitch facets!\n");
      //SVDrawTool::clear_non_retained();
      //SVDrawTool::draw_facets(facets, CUBIT_RED_INDEX);
      //SVDrawTool::mouse_xforms();
    return CUBIT_FALSE;
  }
////////////////////
  
    //TDSkeletonCubitFacet::stitch_facets(facet1,facet2);
      
  for (i=0; i < mark1.size(); ++i)
  {
    CubitFacet *facet = mark1.get_and_step();
    
    if (facet->marked()) {continue;}
    facet->marked(CUBIT_TRUE);
    
    for (j=0; j < mark2.size(); ++j)
    {
      for (k=0; k < 3; ++k)
      {
        CubitFacet *other = facet->shared_facet(facet->point(k), facet->point((k+1)%3));
        if (other != NULL && other == mark2[j] && !mark2[j]->marked() && (facet->normal() % other->normal() < dihedral_angle_thresh))
        {
            //crease_edges.append(facet->point(k)->shared_edge(facet->point((k+1)%3)));
          find_cell_list_intersecting_line_segment(facet->point(k)->coordinates(), facet->point((k+1)%3)->coordinates(), refine_cell_list);
          mark2[j]->marked(CUBIT_TRUE);
          
          if (draw_facets)
          {
            draw_facet_list.append(facet);
            draw_facet_list.append(other);
          }
        }
      }  
    }
  }
  
  if (draw_facets) {
    //SVDrawTool::draw_facets(draw_facet_list, CUBIT_MAGENTA_INDEX);
  }

  if (refine_cell_list.size() > 0)
  {
      //crease_facets += facets;
    FacetDataUtil::delete_facets(facets);
    return CUBIT_TRUE;
  }

  FacetDataUtil::delete_facets(facets);
  return CUBIT_FALSE;
}

void CubitOctree::refine_cells_to_target_depth(DLIList<CubitOctreeCell*> &refine_cell_list, const int target_depth)
{
  int i, j, k;
  CubitOctreeCell *current_cell = NULL;
  
  while (refine_cell_list.size() > 0)
  {
    current_cell = refine_cell_list.get();
    if (subdivide_octree_cell(current_cell, NULL, target_depth))
    {
      for (i=0; i < 2; ++i)
      {
        for (j=0; j < 2; ++j)
        {
          for (k=0; k < 2; ++k)
          {
            if (current_cell->get_child(i,j,k) != NULL)
            {
              refine_cell_list.append(current_cell->get_child(i,j,k));
            }
          }
        }
      }
    }
    refine_cell_list.extract();
  }
}

// Checks for only marked cells and unmarks it if it doesn't interesect
CubitBoolean CubitOctree::unmark_octree_cells_not_intersecting_with_facet( CubitFacet *ptr_facet, DLIList<CubitOctreeCell *> &CubitOctree_cell_list){

  int i;
  CubitOctreeCell *ptr_cell;

  for( i = 0; i < CubitOctree_cell_list.size(); i++ ){
    ptr_cell = CubitOctree_cell_list.get_and_step();
    if( ptr_cell->get_mark() ==  CUBIT_TRUE ){
      if( ptr_cell->does_facet_intersect_octreecell( ptr_facet ) == CUBIT_FALSE ){
        ptr_cell->set_mark( CUBIT_FALSE );
      }
    }
  }
  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::mark_cells_that_intersect_facet_plane( /*CubitFacet *ptr_facet,*/ DLIList<CubitOctreeCell*> &CubitOctree_cell_list )
{
  int i;
  CubitOctreeCell *ptr_cell;
  
  for( i = 0; i < CubitOctree_cell_list.size(); i++ ){
    ptr_cell = CubitOctree_cell_list.get_and_step();
    if( ptr_cell->does_contain_positive_and_negative_nodes() == CUBIT_TRUE ){
      ptr_cell->set_mark( CUBIT_TRUE );
    }
  }
   
  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::subdivide_octree_from_leaf_cell( CubitOctreeCell *ptr_cell, OctreeFacetPointData *facet_point_data, int max_depth ){
  
  CubitVector cell_center = ptr_cell->get_center();
  int i, j, k;
  
  CubitVector coord;

    // Check if any other coord already present in the cell is close to new facet_point_data
    // At the common curves the facet points of adjacent surfaces can come very close to one another
    // Therefore id will not work

  if( ptr_cell->num_of_facet_point_data() != 0 ){
    if( ptr_cell->is_facet_point_data_present( facet_point_data ) ){
      return CUBIT_FALSE;
    }
  }
  

    // If CubitOctree cell contains on fact points just insert the point don't subdivide (Ref. PR-CubitOctree property)
    // 
  if( ptr_cell->num_of_facet_point_data() == 0 || ptr_cell->get_depth() >= max_depth ){
      
    ptr_cell->append_list_item( facet_point_data );
    return CUBIT_TRUE;
    
  }

  coord = facet_point_data->coordinates();  
  bool outside_tolerance = CUBIT_TRUE;
  DLIList<OctreeFacetPointData*> *fpdata = ptr_cell->get_facet_point_data_list();
  for (i=0; i < fpdata->size(); ++i)
  {
    if ((fpdata->get_and_step()->coordinates()-coord).length_squared() <= OCTREE_EPSILON * OCTREE_EPSILON)
    {
      outside_tolerance = CUBIT_FALSE;
    }
  }
  
  if (!outside_tolerance)
  {
    ptr_cell->append_list_item(facet_point_data);
    return CUBIT_TRUE;
    
  }
  else{
  
    subdivide_octree_cell( ptr_cell, NULL, max_depth );
    
    
      // find the leaf cell containing the point
      // call the subdivision recursively from leaf cell containing the point.
    i = j = k = 1;
    
    if( coord.x() < cell_center.x() )
      i = 0; 
    
    if( coord.y() < cell_center.y() )
      j = 0;
    
    if( coord.z() < cell_center.z() )
      k = 0; 
    
    subdivide_octree_from_leaf_cell( ptr_cell->get_child(i, j, k), facet_point_data, max_depth );
    
  }
  return CUBIT_TRUE;
  
}

CubitBoolean CubitOctree::find_octree_cells_contained_inside_bbox( CubitFacet *ptr_facet, DLIList<CubitOctreeCell*> &CubitOctree_cell_list ){
  
  CubitPoint *p0, *p1, *p2;
  CubitBox facet_bbox;
  DLIList<CubitOctreeCell*> queue;
  DLIList<CubitOctreeCell*> marked_list;
  CubitOctreeCell *CubitOctree_cell1, *CubitOctree_cell2, *CubitOctree_cell3;
  CubitOctreeCell *ptr_cell;
  CubitBoolean result;
  int i;
  
  queue.clean_out();
  
  facet_bbox = ptr_facet->bounding_box();
    //PRINT_DEBUG_157(" BBox min : %f %f %f max : %f %f %f ", facet_bbox.minimum().x(), facet_bbox.minimum().y(), facet_bbox.minimum().z(), facet_bbox.maximum().x(),  facet_bbox.maximum().y(),  facet_bbox.maximum().z() );
  
  ptr_facet->points( p0, p1, p2 );
  
  
  CubitOctree_cell1 = find_octreecell( p0->coordinates() );
    //CubitOctree_cell_list.push( CubitOctree_cell1 );
  queue.push( CubitOctree_cell1 );
  CubitOctree_cell1->set_mark( CUBIT_TRUE );
  
  CubitOctree_cell2 = find_octreecell( p1->coordinates() );
    //CubitOctree_cell_list.push( CubitOctree_cell2 );
  if( CubitOctree_cell2->get_mark() == CUBIT_FALSE ){
    queue.push( CubitOctree_cell2 );
    CubitOctree_cell2->set_mark( CUBIT_TRUE );
  }
  
  CubitOctree_cell3 = find_octreecell( p2->coordinates() );
    //CubitOctree_cell_list.push( CubitOctree_cell3 );
  if( CubitOctree_cell3->get_mark() == CUBIT_FALSE ){
    queue.push( CubitOctree_cell3 );
    CubitOctree_cell3->set_mark( CUBIT_TRUE );
  }
  
  while( queue.size() > 0 ){
    
    ptr_cell = queue.pop();
    marked_list.push( ptr_cell );
    
    result = ptr_cell->is_intersects_box( facet_bbox ); 	
    
    if( result == CUBIT_TRUE ){
      CubitOctree_cell_list.push( ptr_cell );
        // add adjacent unmarked CubitOctree cells to queue
      ptr_cell->add_adjacent_unmarked_cells( queue );
    }
    
  }
  
  for( i = 0; i < marked_list.size(); i++ ){		
    ptr_cell = marked_list.get_and_step();
    ptr_cell->set_mark( CUBIT_FALSE );
      //ptr_cell->display(this);
  }
  
  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::mark_positive_and_negative_octree_grid_nodes( CubitFacet *ptr_facet, DLIList<CubitOctreeCell*> &CubitOctree_cell_list, DLIList<CubitOctreeNode *> &CubitOctree_grid_node_list ){
  
  CubitOctreeCell *ptr_cell;
  int i, l, m, n;
  CubitOctreeNode *ptr_grid_node;
  
  for( i = 0; i < CubitOctree_cell_list.size(); i++ ){
    ptr_cell = CubitOctree_cell_list.get_and_step();
    
    for( l = 0; l < 2; l++ ){
      for( m = 0; m < 2; m++ ){
        for( n = 0; n < 2; n++ ){
          ptr_grid_node = ptr_cell->get_octree_grid_node( l, m, n );
	  
          if( ptr_grid_node->get_mark() == CUBIT_FALSE ){
            ptr_grid_node->find_half_space( ptr_facet );
            ptr_grid_node->set_mark( CUBIT_TRUE );
            CubitOctree_grid_node_list.push( ptr_grid_node );
          }
        }
      }
    }
  }
  
  return CUBIT_TRUE;
}


CubitBoolean CubitOctree::find_intersection_between_grid_edges_and_facet(  CubitOctreeType type, RefFace *ptr_face, CubitFacet *ptr_facet, DLIList<CubitOctreeNode *> &CubitOctree_grid_node_list  ){
  
  int i;
  CubitOctreeNode *ptr_grid_node;
  
  for( i = 0; i < CubitOctree_grid_node_list.size(); i++ ){
    ptr_grid_node = CubitOctree_grid_node_list.get_and_step();
    if( ptr_grid_node->get_halfspace_direction() == OCTREE_POSITIVE ){
      ptr_grid_node->find_intersection_with_facet( type,  ptr_face, ptr_facet, boundaryWhiteGridNodeList );
    }
  }
  
    // reset all +ve grid  nodes to -ve
  for( i = 0; i < CubitOctree_grid_node_list.size(); i++ ){
    ptr_grid_node = CubitOctree_grid_node_list.get_and_step();
    ptr_grid_node->set_mark( CUBIT_FALSE );	
    if( ptr_grid_node->get_halfspace_direction() == OCTREE_POSITIVE ){
      ptr_grid_node->set_halfspace_direction( OCTREE_NEGATIVE );
    }
  }
  
  return CUBIT_TRUE;
}


void CubitOctree::gather_initial_grid_nodes_for_smooth_transition( DLIList<CubitOctreeNode*> &queue ){
  int i; 
  CubitOctreeNode *ptr_grid_node;
  int depth_difference;

  for( i = 0; i < gridNodeVector.size(); i++ ){
    ptr_grid_node = gridNodeVector.get_and_step();
    
    depth_difference = ptr_grid_node->find_min_depth_cell_and_depth_difference();
      // By default mark is maintained CUBIT_FALSE
    if( depth_difference >= 2 ){
      queue.push( ptr_grid_node );
      ptr_grid_node->set_mark( CUBIT_TRUE );  
    }
  }

}


// Visit every black node and check the depth of the adjacent cells.
// If the cells differ more than 1 then split the maximum cell.
CubitBoolean CubitOctree::establish_smooth_transition_of_cells( int max_depth ){
  
  CubitOctreeNode *ptr_grid_node;
  int depth_difference;
  DLIList<CubitOctreeNode *> queue;

  gather_initial_grid_nodes_for_smooth_transition( queue );

  while( queue.size() > 0 ){
    ptr_grid_node = queue.remove();
    ptr_grid_node->set_mark( CUBIT_FALSE );
      // WARNING: Why Cubit Black only
      // if( ptr_grid_node->get_color() == CUBIT_BLACK_INDEX ){
    depth_difference = ptr_grid_node->get_cell_depth_difference();
    if( depth_difference >= 2 ){
      subdivide_octree_cell( ptr_grid_node->get_min_depth_cell(), &queue, max_depth );
    }
      //}
  }

  return CUBIT_TRUE;
}

CubitBoolean CubitOctree::subdivide_octree_cell( CubitOctreeCell *ptr_cell, DLIList<CubitOctreeNode*> *ptr_queue, int max_depth ){

  CubitVector cell_center, shift;
  int new_depth, i, j, k, l, m, n, p, q, r;
  double norm_dim, total_dist, new_half_point;
  CubitOctreeCell *ptr_child_cell;
  CubitOctreeNode *local_node_matrix[3][3][3], *ptr_node, *ptr_new_oct_grid_node;
  CubitBoolean new_node_matrix[3][3][3];
  CubitOctreeCell *ptr_new_cell;

  if( ptr_cell->get_depth() < max_depth && ptr_cell->is_leaf() )
  {
    ptr_cell->set_leaf( CUBIT_FALSE );
    norm_dim = ptr_cell->get_dimension();
    new_depth = ptr_cell->get_depth() + 1;
    cell_center = ptr_cell->get_center();
      // PRINT_INFO("Num = %d Depth = %d  Coord = %2.12e %2.12e %2.12e\n", ptr_cell->get_num(), depth, cell_center.x(), cell_center.y(), cell_center.z() );
        
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
	        
          if( i == 0 )
            shift.x( -norm_dim / 4.0 );
          else
            shift.x( norm_dim / 4.0 );
	        
          if( j == 0 )
            shift.y( -norm_dim / 4.0 );
          else
            shift.y( norm_dim / 4.0 );		
	        
          if( k == 0 )
            shift.z( -norm_dim / 4.0 );
          else
            shift.z( norm_dim / 4.0 );
	        
	        
          ptr_child_cell = new CubitOctreeCell( cell_center  + shift, norm_dim / 2.0, new_depth, ptr_cell );
          ptr_cell->set_child( i, j, k, ptr_child_cell );
	        
        }
      }
    }
    
      // PRINT_INFO(" 8 children created \n");
      // Update the wire frame model 
    
    for( i = 0; i < 3; i++ ){
      for( j = 0; j < 3; j++ ){
        for( k = 0; k < 3; k++ ){
          new_node_matrix[i][j][k] = CUBIT_FALSE;
        }
      }
    }
    
      // find the eight corner nodes of the local matrix
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          if( i == 1 )
            l = 2;
          else
            l = 0;
	        
          if( j == 1 )
            m = 2;
          else
            m = 0;
	        
          if( k == 1 )
            n = 2;
          else
            n = 0;
	        
          local_node_matrix[l][m][n] = ptr_cell->get_octree_grid_node( i, j, k );
        }
      }
    }
    
      // PRINT_INFO(" find 8 corner nodes of local 3 x 3 matrix \n");

    new_half_point = 1.0 / pow(2.0, new_depth);
    
      // To find 8
    ptr_node = local_node_matrix[2][0][2];
    total_dist = 0.0;
    while(  total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BOTTOM ) );
      ptr_node = ptr_node->get_adj_node( O_BOTTOM );			
    }
    
    if( total_dist == new_half_point ){
      local_node_matrix[2][0][1] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][2])->x(), (local_node_matrix[2][0][2])->y(), (local_node_matrix[2][0][2])->z() - norm_dim /2 );
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[2][0][1] = ptr_new_oct_grid_node;
      new_node_matrix[2][0][1] = CUBIT_TRUE;
    }
    
    
      // To find 9
    ptr_node = local_node_matrix[2][2][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BOTTOM ) );
      ptr_node = ptr_node->get_adj_node( O_BOTTOM );
    }
    
    if( total_dist == new_half_point ){
      local_node_matrix[2][2][1] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][2][2])->x(), (local_node_matrix[2][2][2])->y(), (local_node_matrix[2][2][2])->z() - norm_dim /2 );
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[2][2][1] = ptr_new_oct_grid_node;
      new_node_matrix[2][2][1] = CUBIT_TRUE;
    }
    
      // To find 10
    ptr_node = local_node_matrix[0][2][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BOTTOM ) );
      ptr_node = ptr_node->get_adj_node( O_BOTTOM );
    }
    
    if( total_dist == new_half_point ){
      local_node_matrix[0][2][1] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[0][2][2])->x(), (local_node_matrix[0][2][2])->y(), (local_node_matrix[0][2][2])->z() - norm_dim /2 );
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[0][2][1] = ptr_new_oct_grid_node;
      new_node_matrix[0][2][1] = CUBIT_TRUE;
    }
    
      // To find 11
    ptr_node = local_node_matrix[0][0][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BOTTOM ) );
      ptr_node = ptr_node->get_adj_node( O_BOTTOM );
    }		
    if( total_dist == new_half_point ){
      local_node_matrix[0][0][1] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[0][0][2])->x(), (local_node_matrix[0][0][2])->y(), (local_node_matrix[0][0][2])->z() - norm_dim /2 );
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[0][0][1] = ptr_new_oct_grid_node;
      new_node_matrix[0][0][1] = CUBIT_TRUE;
    }
    
      // To find 12
    ptr_node = local_node_matrix[2][0][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
      ptr_node = ptr_node->get_adj_node( O_RIGHT );
    }		
    if( total_dist == new_half_point ){
      local_node_matrix[2][1][2] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][2])->x(), (local_node_matrix[2][0][2])->y() + norm_dim / 2.0, (local_node_matrix[2][0][2])->z() ); 
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[2][1][2] = ptr_new_oct_grid_node;
      new_node_matrix[2][1][2] = CUBIT_TRUE;
    }
    
    
      // To find 13
    ptr_node = local_node_matrix[2][0][0];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
      ptr_node = ptr_node->get_adj_node( O_RIGHT );
    }		
    if( total_dist == new_half_point ){
      local_node_matrix[2][1][0] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][0])->x(), (local_node_matrix[2][0][0])->y() + norm_dim / 2.0, (local_node_matrix[2][0][0])->z() ); 
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[2][1][0] = ptr_new_oct_grid_node;
      new_node_matrix[2][1][0] = CUBIT_TRUE;
    }
    
      // To find 14
    ptr_node = local_node_matrix[0][0][0];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
      ptr_node = ptr_node->get_adj_node( O_RIGHT );
    }		
    if( total_dist == new_half_point ){
      local_node_matrix[0][1][0] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[0][0][0])->x(), (local_node_matrix[0][0][0])->y() + norm_dim / 2.0, (local_node_matrix[0][0][0])->z() ); 
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[0][1][0] = ptr_new_oct_grid_node;
      new_node_matrix[0][1][0] = CUBIT_TRUE;
    }
    
    
      // To find 15
    ptr_node = local_node_matrix[0][0][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
      ptr_node = ptr_node->get_adj_node( O_RIGHT );
    }		
    if( total_dist == new_half_point ){
      local_node_matrix[0][1][2] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[0][0][2])->x(), (local_node_matrix[0][0][2])->y() + norm_dim / 2.0, (local_node_matrix[0][0][2])->z() ); 
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[0][1][2] = ptr_new_oct_grid_node;
      new_node_matrix[0][1][2] = CUBIT_TRUE;
    }
    
    
      // To find 16
    ptr_node = local_node_matrix[2][0][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
      ptr_node = ptr_node->get_adj_node( O_BACK );
    }
    
    if( total_dist == new_half_point ){
      local_node_matrix[1][0][2] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][2])->x() - norm_dim / 2.0, (local_node_matrix[2][0][2])->y(), (local_node_matrix[2][0][2])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][0][2] = ptr_new_oct_grid_node;
      new_node_matrix[1][0][2] = CUBIT_TRUE;
    }
    
    
      // To find 17
    ptr_node = local_node_matrix[2][0][0];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
      ptr_node = ptr_node->get_adj_node( O_BACK );
    }
    if( total_dist == new_half_point ){
      local_node_matrix[1][0][0] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][0])->x() - norm_dim / 2.0, (local_node_matrix[2][0][0])->y(), (local_node_matrix[2][0][0])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][0][0] = ptr_new_oct_grid_node;
      new_node_matrix[1][0][0] = CUBIT_TRUE;
    }
    
    
      // To find 18
    ptr_node = local_node_matrix[2][2][0];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
      ptr_node = ptr_node->get_adj_node( O_BACK );
    }
    if( total_dist == new_half_point ){
      local_node_matrix[1][2][0] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][2][0])->x() - norm_dim / 2.0, (local_node_matrix[2][2][0])->y(), (local_node_matrix[2][2][0])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][2][0] = ptr_new_oct_grid_node;
      new_node_matrix[1][2][0] = CUBIT_TRUE;
    }
    
    
      // To find 19
    ptr_node = local_node_matrix[2][2][2];
    total_dist = 0.0;
    while( total_dist < new_half_point ){
      total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
      ptr_node = ptr_node->get_adj_node( O_BACK );
    }
    if( total_dist == new_half_point ){
      local_node_matrix[1][2][2] = ptr_node;
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][2][2])->x() - norm_dim / 2.0, (local_node_matrix[2][2][2])->y(), (local_node_matrix[2][2][2])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][2][2] = ptr_new_oct_grid_node;
      new_node_matrix[1][2][2] = CUBIT_TRUE;
    }
    
      // PRINT_INFO(" found nodes on edges of local 3 x 3 matrix \n");
      // To find 20
    ptr_node = local_node_matrix[2][0][1];
    if( new_node_matrix[2][0][1] == CUBIT_FALSE && ptr_node->get_adj_node( O_RIGHT ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
        ptr_node = ptr_node->get_adj_node( O_RIGHT );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[2][1][1] = ptr_node;
          // PRINT_INFO(" old node 20 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
      
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][1])->x(), (local_node_matrix[2][0][1])->y() + norm_dim / 2.0, (local_node_matrix[2][0][1])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[2][1][1] = ptr_new_oct_grid_node;
      new_node_matrix[2][1][1] = CUBIT_TRUE;
        // PRINT_INFO(" new node 20 created \n");
    }
    
    
    
      // To find 21
    ptr_node = local_node_matrix[0][0][1];
    if( new_node_matrix[0][0][1] == CUBIT_FALSE && ptr_node->get_adj_node( O_RIGHT ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
        ptr_node = ptr_node->get_adj_node( O_RIGHT );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[0][1][1] = ptr_node;
          // PRINT_INFO(" old node 21 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[0][0][1])->x(), (local_node_matrix[0][0][1])->y() + norm_dim / 2.0, (local_node_matrix[0][0][1])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[0][1][1] = ptr_new_oct_grid_node;
      new_node_matrix[0][1][1] = CUBIT_TRUE;
        // PRINT_INFO(" new node 21 created \n");
    }
    
      // To find 22	
    ptr_node = local_node_matrix[2][0][1];
    if( new_node_matrix[2][0][1] == CUBIT_FALSE && ptr_node->get_adj_node( O_BACK ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
        ptr_node = ptr_node->get_adj_node( O_BACK );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[1][0][1] = ptr_node;
          // PRINT_INFO(" old node 22 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][0][1])->x() - norm_dim / 2.0, (local_node_matrix[2][0][1])->y(), (local_node_matrix[2][0][1])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][0][1] = ptr_new_oct_grid_node;
      new_node_matrix[1][0][1] = CUBIT_TRUE;
        // PRINT_INFO(" new node 22 created \n");
    }
    
      // To find 23	
    ptr_node = local_node_matrix[2][2][1];
    if( new_node_matrix[2][2][1] == CUBIT_FALSE && ptr_node->get_adj_node( O_BACK ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_BACK ) );
        ptr_node = ptr_node->get_adj_node( O_BACK );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[1][2][1] = ptr_node;
          // PRINT_INFO(" old node 23 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[2][2][1])->x() - norm_dim / 2.0, (local_node_matrix[2][2][1])->y(), (local_node_matrix[2][2][1])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][2][1] = ptr_new_oct_grid_node;
      new_node_matrix[1][2][1] = CUBIT_TRUE;
        // PRINT_INFO(" new node 23 created \n");
    }
    
      // To find 24
    ptr_node = local_node_matrix[1][0][0];
    if( new_node_matrix[1][0][0] == CUBIT_FALSE && ptr_node->get_adj_node( O_RIGHT ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
        ptr_node = ptr_node->get_adj_node( O_RIGHT );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[1][1][0] = ptr_node;
          // PRINT_INFO(" old node 24 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[1][0][0])->x(), (local_node_matrix[1][0][0])->y() + norm_dim / 2.0, (local_node_matrix[1][0][0])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][1][0] = ptr_new_oct_grid_node;
      new_node_matrix[1][1][0] = CUBIT_TRUE;
        // PRINT_INFO(" new node 24 created \n");
    }
    
    
    
      // To find 25
    ptr_node = local_node_matrix[1][0][2];
    if( new_node_matrix[1][0][2] == CUBIT_FALSE && ptr_node->get_adj_node( O_RIGHT ) != NULL ){
      
      total_dist = 0.0;
      while( total_dist < new_half_point ){
        total_dist += 1.0 / pow( 2.0, ptr_node->get_adj_node_distance( O_RIGHT ) );
        ptr_node = ptr_node->get_adj_node( O_RIGHT );
      }
      if( total_dist == new_half_point ){
        local_node_matrix[1][1][2] = ptr_node;
          // PRINT_INFO(" old node 25 found \n");
      }
      else{
        PRINT_INFO("WARNING: Node not found \n");
      }
    }
    else{
      ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[1][0][2])->x(), (local_node_matrix[1][0][2])->y() + norm_dim / 2.0, (local_node_matrix[1][0][2])->z());
      apd_vector_item( ptr_new_oct_grid_node );
      local_node_matrix[1][1][2] = ptr_new_oct_grid_node;
      new_node_matrix[1][1][2] = CUBIT_TRUE;
        // PRINT_INFO(" new node 25 created \n");
    }
    
      // PRINT_INFO(" found nodes on faces of local 3 x 3 matrix \n");
    
      // To find 26
    ptr_new_oct_grid_node = new CubitOctreeNode( (local_node_matrix[1][1][2])->x(), (local_node_matrix[1][1][2])->y(), (local_node_matrix[1][1][2])->z() - norm_dim / 2.0);
    apd_vector_item( ptr_new_oct_grid_node );
    local_node_matrix[1][1][1] = ptr_new_oct_grid_node;
    new_node_matrix[1][1][1] = CUBIT_TRUE;
    
      // PRINT_INFO(" found nodes at center of  local 3 x 3 matrix \n");
    
      // update the eight nodes of new eight cells
    
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
	  
          ptr_new_cell = ptr_cell->get_child( i, j, k );
	        
          for( l = 0; l < 2; l++ ){
            for( m = 0; m < 2; m++ ){
              for( n = 0; n < 2; n++ ){
                ptr_new_cell->set_oct_grid_node( l, m, n, local_node_matrix[i+l][j+m][k+n] );
  	
                p = ( l == 1 )?0:1;  // fliped 0 -> 1;  1 -> 0
                q = ( m == 1 )?0:1;
                r = ( n == 1 )?0:1;
                  // update the eight cells of all 27 nodes
                local_node_matrix[i+l][j+m][k+n]->update_adj_cell( ptr_new_cell, p, q, r );
              }
            }
          }
        }
      }
    }

      // PRINT_INFO(" updates nodes of the 8 children and visa versa \n");
    
      // connnect the new nodes with the old nodes to update the wire frame
    
    for( i = 0; i < 3; i++ ){
      for( j = 0; j < 3; j++ ){
        for( k = 0; k < 3; k++ ){
          if( new_node_matrix[i][j][k] == CUBIT_TRUE ){
            if( i - 1 >= 0 ){
              local_node_matrix[i][j][k]->set_adj_node( O_BACK, local_node_matrix[i-1][j][k] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_BACK, new_depth );
              local_node_matrix[i-1][j][k]->set_adj_node( O_FRONT, local_node_matrix[i][j][k] );
              local_node_matrix[i-1][j][k]->set_adj_node_distance( O_FRONT, new_depth );
            }
            if( i + 1 <= 2 ){
              local_node_matrix[i][j][k]->set_adj_node( O_FRONT, local_node_matrix[i+1][j][k] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_FRONT, new_depth );
              local_node_matrix[i+1][j][k]->set_adj_node( O_BACK, local_node_matrix[i][j][k] );				
              local_node_matrix[i+1][j][k]->set_adj_node_distance( O_BACK, new_depth );
            }
	          
            if( j - 1 >= 0 ){
              local_node_matrix[i][j][k]->set_adj_node( O_LEFT, local_node_matrix[i][j-1][k] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_LEFT, new_depth );
              local_node_matrix[i][j-1][k]->set_adj_node( O_RIGHT, local_node_matrix[i][j][k] );
              local_node_matrix[i][j-1][k]->set_adj_node_distance( O_RIGHT, new_depth );
            }
            if( j + 1 <= 2 ){
              local_node_matrix[i][j][k]->set_adj_node( O_RIGHT, local_node_matrix[i][j+1][k] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_RIGHT, new_depth );
              local_node_matrix[i][j+1][k]->set_adj_node( O_LEFT, local_node_matrix[i][j][k] );
              local_node_matrix[i][j+1][k]->set_adj_node_distance( O_LEFT, new_depth );
            }
	          
            if( k - 1 >= 0 ){
              local_node_matrix[i][j][k]->set_adj_node( O_BOTTOM, local_node_matrix[i][j][k-1] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_BOTTOM,  new_depth );
              local_node_matrix[i][j][k-1]->set_adj_node( O_TOP, local_node_matrix[i][j][k] );
              local_node_matrix[i][j][k-1]->set_adj_node_distance( O_TOP, new_depth );
            }
            if( k + 1 <= 2 ){
              local_node_matrix[i][j][k]->set_adj_node( O_TOP, local_node_matrix[i][j][k+1] );
              local_node_matrix[i][j][k]->set_adj_node_distance( O_TOP, new_depth );
              local_node_matrix[i][j][k+1]->set_adj_node( O_BOTTOM, local_node_matrix[i][j][k] );
              local_node_matrix[i][j][k+1]->set_adj_node_distance( O_BOTTOM, new_depth );
            }
	          
          }
        }
      }
    }
    
      // distribute the facet_points and Ref_face among the children
    ptr_cell->distribute_facet_points_among_children();
    
    
    int depth_difference;
      // update the depthMaxDepthCell, depthMinDepthCell, cellDepthDifference, and minDepthCell 
    if( ptr_queue != NULL ){
      for( i = 0; i < 3; i++ ){
        for( j = 0; j < 3; j++ ){
          for( k = 0; k < 3; k++ ){
            depth_difference = local_node_matrix[i][j][k]->find_min_depth_cell_and_depth_difference( );        
              /*  RISKY TO OPTIMIZE ANY THING HERE 
                  NEEDS CAREFULL STUDY
                  if( new_node_matrix[i][j][k] == CUBIT_TRUE ){
                    // initialize the variables depthMaxDepthCell, depthMinDepthCell, cellDepthDifference, and minDepthCell  
                    local_node_matrix[i][j][k]->find_min_depth_cell_and_depth_difference(  );        
                    }
                    else{
                      // update the variables depthMaxDepthCell, depthMinDepthCell, cellDepthDifference, and minDepthCell 
                      max_depth = local_node_matrix[i][j][k]->get_depth_max_depth_cell();
                      min_depth = local_node_matrix[i][j][k]->get_depth_min_depth_cell();
                      if( new_depth < min_depth || new_depth > max_depth )

                      }
              */
            if( local_node_matrix[i][j][k]->get_mark() == CUBIT_FALSE  && depth_difference >= 2 ){
              ptr_queue->push( local_node_matrix[i][j][k] );
              local_node_matrix[i][j][k]->set_mark( CUBIT_TRUE );
            }
          }
        }
      }
    }
//    ptr_cell->delete_local_node_matrix();
    return CUBIT_TRUE;
  }
  else{
    return CUBIT_FALSE;
  }
}


void CubitOctree::color_octreecell( void ){
  
  root->coloring(/*greyCellList,*/ blackCellList);
  
    /* TESTING: Duplicates in grey cell list 
       CubitOctreeCell *ptr_cell;
       int i;
       for( i = 0; i < greyCellList.size(); i++ ){
       ptr_cell = greyCellList.get_and_step();
       PRINT_INFO("Grey cell ID: %d \n", ptr_cell->id() );
       }
    */
}

// p0 and p1 should be within the root's bbox
// note: this is a little sloppy, i'll clean it up when I fix the general neighbor traversal
void CubitOctree::find_cell_list_intersecting_line_segment (const CubitVector &p0, const CubitVector &p1, DLIList<CubitOctreeCell*> &cell_list){

//  plane order is like in enum OctreePosition
  
  CubitOctreeCell *start = find_octreecell(p0), *end = find_octreecell(p1);
  CubitOctreeCell *current = start, *previous = NULL;
  CubitVector curr_pt = p0;
  CubitVector dir = p1-p0, center;
    //CubitVector ii(1,0,0);
    //CubitVector jj(0,1,0);
    //CubitVector kk(0,0,1);

  int i;
  dir.normalize();
  double tx=0, ty=0, tz=0, h=0, a = dir.x(), b = dir.y(), c = dir.z();

  static unsigned int face_corners[6][4] =
      {
          {0,5,6,2},
          {4,3,1,7},
          {5,3,7,2},
          {0,4,1,6},
          {1,7,2,6},
          {4,3,5,0}
      };

  static unsigned int inwward_cells[6][4][3] =
      {
          {{0,0,0},{0,0,1},{1,0,0},{1,0,1}},
          {{0,1,0},{0,1,1},{1,1,0},{1,1,1}},
          {{0,1,0},{0,0,0},{1,0,0},{1,1,0}},
          {{0,1,1},{0,0,1},{1,0,1},{1,1,1}},
          {{0,0,0},{0,0,1},{0,1,1},{0,1,0}},
          {{1,0,0},{1,0,1},{1,1,1},{1,1,0}}
      };

  while (current != end)
  {
    if (current != previous)
    {
      cell_list.append(current);
      previous = current;
    }
      // do tn, pick face, march stuff here
    h = current->get_dimension();
    center = current->get_center();
    if (fabs(a) > CUBIT_DBL_MIN) {tx = h/(2*fabs(a)) - (curr_pt.x() - center.x())/a;}
    else {tx = CUBIT_DBL_MAX;}
    if (fabs(b) > CUBIT_DBL_MIN) {ty = h/(2*fabs(b)) - (curr_pt.y() - center.y())/b;}
    else {ty = CUBIT_DBL_MAX;}
    if (fabs(c) > CUBIT_DBL_MIN) {tz = h/(2*fabs(c)) - (curr_pt.z() - center.z())/c;}
    else {tz = CUBIT_DBL_MAX;}

    tx = fabs(tx);
    ty = fabs(ty);
    tz = fabs(tz);
    
    double t = 0.0;

    OctreePosition int_plane = O_UNSET;
    
    if (tx <= ty && tx <= tz)
    {
      if (a > 0) {t = ((center.x() + h/2) - curr_pt.x())/a; int_plane = O_FRONT;}
      else {t = ((center.x() - h/2) - curr_pt.x())/a; int_plane = O_BACK;}
    }
    
    else if (ty <= tx && ty <= tz)
    {
      if (b > 0) {t = ((center.y() + h/2) - curr_pt.y())/b; int_plane = O_RIGHT;}
      else {t = ((center.y() - h/2) - curr_pt.y())/b; int_plane = O_LEFT;}
    }
    
    else if (tz <= ty && tz <= tx)
    {
      if (c > 0) {t = ((center.z() + h/2) - curr_pt.z())/c; int_plane = O_TOP;}
      else {t = ((center.z() - h/2) - curr_pt.z())/c; int_plane = O_BOTTOM;}
    }
    
      // now we have determined the point of intersection on a face of the cell
      // based on the face, we should now check the grid nodes' adjacent cells and march if necessary to find the next cell to march in

    CubitOctreeNode* corners[] =
        {
            current->get_octree_grid_node(1,0,1),
            current->get_octree_grid_node(0,1,1),
            current->get_octree_grid_node(0,0,0),
            current->get_octree_grid_node(1,1,0),
          
            current->get_octree_grid_node(1,1,1),
            current->get_octree_grid_node(1,0,0),
            current->get_octree_grid_node(0,0,1),
            current->get_octree_grid_node(0,1,0),
        };
      /*
        CubitOctreeNode *cADF = current->get_octree_grid_node(1,0,1); //0
        CubitOctreeNode *cBCF = current->get_octree_grid_node(0,1,1); //1
        CubitOctreeNode *cCDE = current->get_octree_grid_node(0,0,0); //2
        CubitOctreeNode *cABE = current->get_octree_grid_node(1,1,0); //3

        CubitOctreeNode *cABF = current->get_octree_grid_node(1,1,1); //4
        CubitOctreeNode *cADE = current->get_octree_grid_node(1,0,0); //5
        CubitOctreeNode *cCDF = current->get_octree_grid_node(0,0,1); //6
        CubitOctreeNode *cBCE = current->get_octree_grid_node(0,1,0); //7
      */
      
  

      // each digit is an i,j,k index for the inward cell from a grid node

    CubitOctreeNode* int_face[4] = {corners[face_corners[int_plane][0]], corners[face_corners[int_plane][1]], corners[face_corners[int_plane][2]], corners[face_corners[int_plane][3]]};
    unsigned int* cells[4] = {inwward_cells[int_plane][0], inwward_cells[int_plane][1], inwward_cells[int_plane][2], inwward_cells[int_plane][3]};

    CubitOctreeCell* four_cells[4];
      //CubitOctreeCell *curr_cell = int_face[0]->get_adj_cell(cells[0][0],cells[0][1],cells[0][2]);
    bool one_face_neighbor = true;
    //bool same_depth = true;
    CubitOctreeCell *compare_to = NULL;
    
    for (i=0; i < 4; ++i)
    {
      four_cells[i] = int_face[i]->get_adj_cell(cells[i][0],cells[i][1],cells[i][2]);
      if (compare_to == NULL && four_cells[i] != NULL) {compare_to = four_cells[i];}
    }
    
    for (i=0; i < 4; ++i)
    {
      if (four_cells[i] == NULL)
      {
        one_face_neighbor = false;
        //same_depth = false;
        break;
      }
      if (four_cells[i] != NULL)
      {
        if (four_cells[i] != compare_to) {one_face_neighbor = false;}
        //if (four_cells[i]->get_depth() != compare_to->get_depth()) {same_depth = false;}
      }
    }

    curr_pt += t*dir;
    if (one_face_neighbor == true)
    {
      current = compare_to;//four_cells[0];//curr_cell;.
//      PRINT_INFO("Found one face neighbor for CubitOctree ray traversal\n");
    }

    else
    {
        //curr_pt += (t+.0001)*dir;
      int axes[6] = {1,1,2,2,0,0};
      int dirs[6] = {-1,1,-1,1,-1,1};
      double xyz[3];
      curr_pt.get_xyz(xyz);

      xyz[axes[int_plane]] += dirs[int_plane] * (get_bbox_dimension()/(2*pow((double)2.0,maxDepth)));

      CubitVector temp_pt;
      temp_pt.set(xyz);
        //    SVDrawTool::draw_point(temp_pt, CUBIT_WHITE_INDEX);
      
      current = find_octreecell(temp_pt);
        /* if (current == previous) {
           PRINT_INFO("Samet method failed!\n");
           }*/
      
        //PRINT_INFO("Did not find one face neighbor for CubitOctree ray traversal, following Samet method\n");
    }
      //SVDrawTool::mouse_xforms();
    if (current == end) {break;}
    
  }
  cell_list.append(end);
  
}



CubitOctreeCell *CubitOctree::find_octreecell( const CubitVector &pnt ){
  
  return( root->find_leaf_octree_cell( pnt ) );
  
}

CubitBoolean CubitOctree::apd_vector_item( CubitOctreeNode *ptr_node ){
  
  gridNodeVector.push( ptr_node);
  return CUBIT_TRUE;
} 


void CubitOctree::find_max_min_size_grid_node_for_scaling( void ){
  CubitOctreeNode *ptr_grid_node;
  int i;
  double size;

  for( i = 0; i < gridNodeVector.size(); i++ ){

    ptr_grid_node = gridNodeVector.get_and_step();
    if (ptr_grid_node)
    {
              
      size = ptr_grid_node->get_size(  OCTREE_SIZE_DEFAULT );
            
      if( size > 0.0 ){
          // TESTING
          //PRINT_DEBUG_157("size = %f\n", size);
              
        if( size < minSizeGridNode )
          minSizeGridNode = size;
              
        if( size > maxSizeGridNode )
          maxSizeGridNode = size;
      }
    }
            
  }

  PRINT_DEBUG_167("Max Size Grid Node =  %f\n", maxSizeGridNode );
  PRINT_DEBUG_167("Min Size Grid Node =  %f\n", minSizeGridNode );
  double average_size = (minSizeGridNode + maxSizeGridNode) / 2.0;
  minSizeGridNode -= .05*average_size;
  maxSizeGridNode += .05*average_size;
  PRINT_DEBUG_167("For drawing:\n");
  PRINT_DEBUG_167("Adjusted Max Size Grid Node =  %f\n", maxSizeGridNode );
  PRINT_DEBUG_167("Adjusted Min Size Grid Node =  %f\n", minSizeGridNode );

    //system("pause");
    
}



double CubitOctree::size_at_a_point( const CubitVector &point, int size_type ){
//PRINT_DEBUG_116("Entered the SizeAtPoint of SkeletonSizingFunction");
  static int mode = 1;
  double size = 0.0; 

//  SVDrawTool::draw_point(point, CUBIT_RED_INDEX);
  
  
  CubitOctreeCell *ptr_cell = find_octreecell( point );
  
  if( ptr_cell == NULL )
  {
    
    PRINT_INFO("WARNING: Point is outside the solid\n");
      //exit(1);
    if( mode == 0 )
      return size;

      /*
        mode = 0;
			
          // Search the adjacent cells and find the nearest sizing data
          CubitVector xyz[6];
          CubitVector mod_point;
          int i;
      
          for( i = 0; i < 5; i ++ ){
          xyz[i].x(0.0);
          xyz[i].y(0.0);
          xyz[i].z(0.0);
          }
      
          xyz[0].x( grid_size );
          xyz[1].x( -1 * grid_size );
          xyz[2].y( grid_size );
          xyz[3].y( -1 * grid_size );
          xyz[4].z( grid_size );
          xyz[5].z( -1 * grid_size );
      
          double total_size = 0;
          int count = 0;
          double size1;
      
          for( i = 0; i < 5; i ++ ){
          mod_point = point + xyz[i];
          size1 = size_at_point( mod_point, 0.0 );
      
          if( size1 != 0 ){
          total_size += size1;
          count++;
          }
          }
      
          size = total_size / count;
          if( size == 0.0 ){
          PRINT_INFO("WARNING: Point is outside the solid can't handle ( size set to average size )\n");
          size = 0.5;
          size = 2 * ( min_rad + ( size - 1 ) / (10 - 1) * ( max_rad - min_rad ) ) * grid_size * over_all_reduction_scale;
          }
      
          mode = 1;
      */
  }
  else
  {
    
    switch( DEFAULT_INTERPOLATION_INSIDE_CELL ) 
    {
    case TRILINEAR_INTERPOLATION_INSIDE_CELL:
      size = ptr_cell->trilinear_interpolation( point );
      break;
    
    case INVERSE_DISTANCE_INTERPOLATION_INSIDE_CELL:
      size = ptr_cell->inverse_distance_interpolation( point );
      break;

    case MIN_DISTANCE_INTERPOLATION_INSIDE_CELL:
      size = ptr_cell->min_distance_interpolation( point );
      break;

    case MIN_SIZE_INTERPOLATION_INSIDE_CELL:
      size = ptr_cell->min_size_interpolation( point );
      break;

    default:
#ifndef NDEBUG
        //PRINT_INFO(" WARNING: Interpolation method inside a cell is not set in SkeletonConstants.hpp\n");
#endif
      size = ptr_cell->trilinear_interpolation( point );
      break;
    }		
   
    if( size_type == SCALED_SIZE )
    {
      size = get_scaled_from_wrld_size_grid_node( size );
			
    }

		
      //	if( size_type == MESH_SIZE ){
      // Do nothing send size for meshing
      //}

  }
    // PRINT_INFO("Size at point <%f,%f,%f> = %f\n", point.x(), point.y(), point.z(), size);
  return size;
}


CubitPointContainment CubitOctree::point_containment( CubitVector point, double tolerance )
{
  
  CubitOctreeCell *ptr_cell = find_octreecell( point );
  
  if( ptr_cell == NULL )
  {
    // not known and needs full acis/mbg point containment
    return CUBIT_PNT_BOUNDARY;
  }
  else
  {
    switch( ptr_cell->get_color() )
    {
      case CUBIT_BLACK_INDEX:
        return CUBIT_PNT_INSIDE;
    
      case CUBIT_WHITE_INDEX:
        return CUBIT_PNT_OUTSIDE;
        
      case CUBIT_GREY_INDEX:
        return CUBIT_PNT_BOUNDARY;
    
      default:
        return CUBIT_PNT_BOUNDARY;
    }
  }
  return CUBIT_PNT_BOUNDARY;
}



double CubitOctree::get_scaled_from_wrld_size_grid_node( double wrld_size ){

  double scaled_size;
	
  scaled_size = 1.0 + ( ( wrld_size - minSizeGridNode ) / ( maxSizeGridNode - minSizeGridNode ) ) * 9.0;

    //PRINT_DEBUG_157("Scaled Size = %2.10f\n",scaled_size);
    //PRINT_INFO(" world = %f, scaled = %f\n", wrld_size, scaled_size );
  if( scaled_size > 10 ){
    PRINT_INFO("Testing ");
  }

  return scaled_size;
}





CubitStatus CubitOctree::circumcenter(CubitVector &a, CubitVector &b, CubitVector &c, CubitVector &center, double &radius)
{
  
    // Use coordinates relative to point `a' of the triangle. 
  CubitVector vec_ba(a, b);
  CubitVector vec_ca(a, c);

    // Squares of lengths of the edges incident to `a'. 
  double ba_length = vec_ba.length_squared();
  double ca_length = vec_ca.length_squared();

    // Cross product of these edges. 
    // (Take your chances with floating-point roundoff.)
  CubitVector cross_bc = vec_ba * vec_ca;

    // Calculate the denominator of the formulae. 
  double denominator = (cross_bc % cross_bc);
  if (denominator == 0.0) {return CUBIT_FAILURE;}
  assert(denominator != 0.0); // don't think we need this..

    // Calculate offset (from `a') of circumcenter. 
  center  = (ba_length * vec_ca - ca_length * vec_ba) * cross_bc * 0.5;
  center /= denominator;

    // radius is length from point `a' to center
  radius = center.length();

    // Add point `a' to get global coordinate of center
  center += a;

  return CUBIT_SUCCESS;
}

double CubitOctree::capsule_distance_to_facet(const CubitVector &point, CubitFacet *lp_facet, CubitVector &int_point, CubitBoolean use_projection_only)
{
    
  CubitVector c = lp_facet->center();
  CubitVector n = lp_facet->normal();
  double d = (point - c)%n;
  CubitVector proj = point - d*n;
  CubitVector points[3] = {lp_facet->point(0)->coordinates(), lp_facet->point(1)->coordinates(), lp_facet->point(2)->coordinates()};
  
  if (use_projection_only || CubitOctreeNode::is_intersection_point_contained_inside_facet(proj, points[0], points[1], points[2]))
  {
    int_point = proj;
    return fabs(d);
  }

  double dn[3] = {0,0,0};
  CubitVector pn[3];
  double t = 0;
  int i;

  for (i=0; i < 3; ++i)
  {
    CubitVector dir = points[(i+1)%3] - points[i];
    double length = dir.length();
    dir.normalize();
    t = ((point - points[i])%dir)/length;
    if (t >= 0 && t <= 1)
    {
      pn[i] = (length*t*dir + points[i]);
      dn[i] = (point - (t*length*dir + points[i])).length();
    }
    else
    {
      if (t < 0)
      {
        pn[i] = points[i];
        dn[i] = (point - points[i]).length();
      }
      else if (t > 1)
      {
        pn[i] = points[(i+1)%3];
        dn[i] = (point - points[(i+1)%3]).length();
      }
        //else {PRINT_INFO("Hmmm.... t is not within bounds: %f!!!!!!!!!!!!!!!!11\n", t);}
    }
  }
  
  if (dn[0] <= dn[1] && dn[0] <= dn[2]) {int_point = pn[0]; return dn[0];}
  if (dn[1] <= dn[0] && dn[1] <= dn[2]) {int_point = pn[1]; return dn[1];}
  if (dn[2] <= dn[1] && dn[2] <= dn[0]) {int_point = pn[2]; return dn[2];}

  else
  {
      //PRINT_INFO("capsule dist not found!!!!!!!!!!!!!1\n");
    return -1;
  }
}


#ifndef NDEBUG
void CubitOctree::write_sizing_info_file( const char *file_name ){
     
  int i,j,k;
  double size;
  
  FILE *pof = fopen( file_name, "w" );
  if( pof == NULL ){
    PRINT_INFO("ERROR: Opening Sizing Output File ");
    exit(0);
  }
     
  fprintf( pof, "TFD1\n" );
  int num_of_cells = (int)pow(2.0, maxDepth-1 );
  fprintf( pof, "%d %d %d\n", num_of_cells, num_of_cells, num_of_cells ); 
    //fprintf( pof, "%f %f %f %f %f %f\n", bboxMaxY, bboxMinZ, bboxMinX, bboxMinY, bboxMaxZ, bboxMaxX );
  fprintf( pof, "%f %f %f %f %f %f\n", bboxMinX, bboxMinY, bboxMinZ, bboxMaxX, bboxMaxY, bboxMaxZ );

  int l = num_of_cells + 1;
  int m = num_of_cells + 1;
  int n = num_of_cells + 1;
  CubitVector point;

  double grid_size = fabs( bboxMinY - bboxMaxY ) / num_of_cells; 
  
  for( k = 0; k < n; k++ )
  {
    for( j = 0; j < m; j++ )
    {
      for ( i = 0; i < l; i++ )
      {
        point.x( bboxMinX + grid_size * i );
        point.y( bboxMinY + grid_size * j );
        point.z( bboxMinZ + grid_size * k );
        
        size = size_at_a_point( point, OCTREE_SIZE_DEFAULT);	
        fprintf( pof, "%f %f %f %f %f %f %f %f %f %f %f %f\n",
                 1.0, 0.0, 0.0, size,
                 0.0, 1.0, 0.0, size,
                 0.0, 0.0, 1.0, size );
      }
    }
  }
  
  fclose( pof );
}

void CubitOctree::write_matlab_sizing_info_file( const char *file_name ){
     
  int i,j,k;
  double size;
  
  CubitOctreeNode *grid_node;

  for (i=0; i < gridNodeVector.size(); ++i)
  {
    grid_node = gridNodeVector.get_and_step();
    if (grid_node->get_color() == CUBIT_GREY_INDEX)
    {
      grid_node->set_size(0.0, OCTREE_SIZE_DEFAULT);
    }
  }
  
  FILE *pof = fopen( file_name, "w" );
  if( pof == NULL ){
    PRINT_INFO("ERROR: Opening Sizing Output File ");
    exit(0);
  }
     
  
  int num_of_cells = 50;//(int)pow(2.0, MATLAB_OUTPUT_MAX_DEPTH );
    //fprintf( pof, "%f %f %f %f %f %f\n", bboxMinX, bboxMinY, bboxMinZ, bboxMaxX, bboxMaxY, bboxMaxZ );

  int l = num_of_cells + 1;
  int m = num_of_cells + 1;
  int n = num_of_cells + 1;
  CubitVector point;

  double grid_size = fabs( bboxMinY - bboxMaxY ) / num_of_cells; 
  
  for( k = 0; k < n; k++ ){
    for( j = 0; j < m; j++ ){
      for ( i = 0; i < l; i++ ){			
        point.x( bboxMinX + grid_size * i );
        point.y( bboxMinY + grid_size * j );
        point.z( bboxMinZ + grid_size * k );
        
        size = size_at_a_point( point, OCTREE_SIZE_DEFAULT);	

        fprintf( pof, "%f %f %f %f\n",
                 point.x(), point.y(), point.z(), size );

      }
    }
  }
  
  fclose( pof );
} 

void CubitOctree::write_octree_sizing_info_file( const char *file_name ){
     
  int i;
  double size;
  CubitOctreeNode *ptr_onode;
  double FACTOR_FOR_BUBBLEMESH = 1 / 1.8;  // Based on 20% overlap of bubbles
  
  FILE *pof = fopen( file_name, "w" );
  if( pof == NULL ){
    PRINT_INFO("ERROR: Opening Sizing Output File ");
    exit(0);
  }
  
  DLIList<CubitOctreeCell*> stack;
  CubitOctreeCell *ptr_cell;

  fprintf( pof, "OCT 1\n" );
//   int num_of_cells = (int)pow(2.0, skeletonProxy->get_max_depth()-1 );
  fprintf( pof, "B %f %f %f %f %f %f\n", 
           bboxMinX, bboxMinY, bboxMinZ, bboxMaxX, bboxMaxY, bboxMaxZ );

  gridNodeVector.reset();
  for( i = 0; i < gridNodeVector.size(); i++ ){
    ptr_onode = gridNodeVector.get_and_step();
    size = ptr_onode->get_size(OCTREE_SIZE_DEFAULT);
		
    fprintf( pof, "V %d %f %f %f %f \n", 
             ptr_onode->get_num(), ptr_onode->x(), ptr_onode->y(), 
             ptr_onode->z(), size *  FACTOR_FOR_BUBBLEMESH  );
  }

  stack.append( root );

  while( stack.size() > 0 ){
    ptr_cell = stack.remove();
    ptr_cell->write_octreecell_sizing_info_file( pof, stack );
  }
	
  fprintf( pof, "X\n");
  fclose( pof );
}

#endif



//EOF
