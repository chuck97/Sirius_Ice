#include <map>
#include "CubitOctreeNode.hpp"
#include "CubitOctreeCell.hpp"
#include "CubitOctree.hpp"
#include "RefFace.hpp"
#include "CubitOctreeConstants.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "GfxDebug.hpp"
#include "PriorityQueue.hpp" 
#include "CubitOctreeGenerator.hpp"
#include "GMem.hpp"
#include "RTree.hpp"
#include "RefEdge.hpp"
#include "GeometryQueryTool.hpp"
#include "DLIList.hpp"
#include "RefFace.hpp"
#include "GeomMeasureTool.hpp"
#include "OctreeIntersectionData.hpp"


//#include "SVDrawTool.hpp"

int CubitOctreeNode::mCounter = -1;

int sort_by_size( CubitOctreeNode * &node1, CubitOctreeNode *&node2 )
{
  if (node1->get_size(OCTREE_SIZE_DEFAULT) <= node2->get_size(OCTREE_SIZE_DEFAULT)) {return -1;}
  else {return 1;}
}

 
/* --------------- Methods of CubitOctreeNode ------------- */
CubitOctreeNode::CubitOctreeNode( const CubitVector &cen, CubitOctreeCell *parent_cell, const int x_index, const int y_index, const int z_index )
{  
  initialize_constructor( cen.x(), cen.y(), cen.z() );
  adjCell[x_index][y_index][z_index] = parent_cell;
}

CubitOctreeNode::CubitOctreeNode( const double &x, const double &y, const double &z )
{
  initialize_constructor( x, y, z );
}

CubitOctreeNode::~CubitOctreeNode(void)
{
  int i;
  OctreeIntersectionData *idata;
  
  for (i=0; i < octreeIntersectionDataList.size(); ++i)
  {
    idata = octreeIntersectionDataList.get_and_step();
    if (idata != NULL)
    {
      delete idata;
    }
  }
  
  octreeIntersectionDataList.clean_out();
}

/// in CubitOctree outside nodes are GREY, inside nodes are BLACK.  And boundary nodes are WHITE.
void CubitOctreeNode::initialize_constructor( const double &x, const double &y, const double &z ){
  
  
  int i, j, k;
  color = CUBIT_GREY_INDEX;  // boundary nodes are colored white and inside nodes black
                       // during facet intersection the white nodes dominates black
  visit = CUBIT_FALSE; // mat generation
  mark = CUBIT_FALSE;
  
  halfspaceDirection = OCTREE_NEGATIVE;
  num = get_counter();
  
  size = 0.0;

   distance = CUBIT_DBL_MAX;
  
  cellDepthDifference = 0;
  minDepthCell = NULL;

  mNormal.x(0.0);
  mNormal.y(0.0);
  mNormal.z(0.0);

  coord.x(x);
  coord.y(y);
  coord.z(z);
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        adjCell[i][j][k] = NULL;
      }
    }
  }
  
    // LRBTFN
  for( i = 0; i < 6; i++ ){
    adjGridNode[i] = NULL;
    adjNodeDistance[i] = -1; // ?
  }
  
 refFace = NULL;
}

int CubitOctreeNode::get_counter()
{
    //static int counter = -1;
  mCounter++;
  return(mCounter);
}

void CubitOctreeNode::display( OctreeNodeConstant type, float draw_size )
{
  switch( type ){
    case NODE_SIZE:
        if (draw_size < 0)
        {
          draw_size = size;
        }
        //SkeletonDebug::draw_point( coord, draw_size );
        break;
    
    case NODE_DISTANCE:
        if( color == CUBIT_BLACK_INDEX ){
          //SkeletonDebug::draw_point( coord,  (float)draw_size ) ;
                   }
        break;
    
    case NODE_FACE_NUM:
        if( refFace != NULL ){
          //SkeletonDebug::draw_point( coord, (float)(( refFace->id() % 18) * 0.5)   );
        }
        else{
          PRINT_INFO("ERROR: Face Ptr Doesn't Exist in White Grid Node \n");
        }
        break;
    
    case NODE_NORMAL:
      if( color == CUBIT_BLACK_INDEX ){
        //SkeletonDebug::draw_line(coord, draw_size, coord+mNormal,draw_size);

      }
      break;
      

    
    default:
        //SkeletonDebug::draw_point( coord, 1 );
        break;
  }
}



double CubitOctreeNode::get_size( OctreeSourceEntityType type ) const
{
    return size;

}

void CubitOctreeNode::set_size( double s, int type ){
  
        size = s;
}


CubitOctreeNode* CubitOctreeNode::get_adj_node(int select)
{
  switch (select)
  {
    case O_LEFT:
        return adjGridNode[O_LEFT];
    
    case O_RIGHT:
        return adjGridNode[O_RIGHT];
    
    case O_BOTTOM:
        return adjGridNode[O_BOTTOM];
    
    case O_TOP:
        return adjGridNode[O_TOP];
    
    case O_BACK:
        return adjGridNode[O_BACK];
    
    case O_FRONT:
        return adjGridNode[O_FRONT];
    
    default:
        PRINT_INFO("ERROR: AdjNode index exceeded\n");
        return NULL;
  }
}

void CubitOctreeNode::set_adj_node(enum OctreePosition select, CubitOctreeNode *ptr_node)
{
  switch( select )
  {
    case O_LEFT:
        adjGridNode[O_LEFT] = ptr_node;
        break;
    
    case O_RIGHT:
        adjGridNode[O_RIGHT] = ptr_node;
        break;
    
    case O_BOTTOM:
        adjGridNode[O_BOTTOM] = ptr_node;
        break;
    
    case O_TOP:
        adjGridNode[O_TOP] = ptr_node; 
        break;
    
    case O_BACK:
        adjGridNode[O_BACK] = ptr_node;
        break;
    
    case O_FRONT:
        adjGridNode[O_FRONT] = ptr_node;
        break;
    
    default:
        PRINT_INFO("ERROR: AdjNode index exceeded");
    
  }
  
}


void CubitOctreeNode::set_adj_node_distance( enum OctreePosition select, int dist ){
  
  switch( select ){
    
    case O_LEFT:
        adjNodeDistance[O_LEFT] = dist;
        break;
    
    case O_RIGHT:
        adjNodeDistance[O_RIGHT] = dist;
        break;
    
    case O_BOTTOM:
        adjNodeDistance[O_BOTTOM] = dist;
        break;
    
    case O_TOP:
        adjNodeDistance[O_TOP] = dist; 
        break;
    
    case O_BACK:
        adjNodeDistance[O_BACK] = dist;
        break;
    
    case O_FRONT:
        adjNodeDistance[O_FRONT] = dist;
        break;
    
    default:
        PRINT_INFO(" ERROR: AdjNode index exceeded");
    
  }
  
}


double CubitOctreeNode::manhattan_distance_adj_node( int index ){
  
  switch( index ){

    case O_FRONT:
        return fabs( coord.x() - adjGridNode[index]->coord.x() );
    
    case O_BACK:
        return fabs( coord.x() - adjGridNode[index]->coord.x() );
    
    case O_RIGHT:
        return fabs( coord.y() - adjGridNode[index]->coord.y() );
    
    case O_LEFT:
        return fabs( coord.y() - adjGridNode[index]->coord.y() );
    
    case O_TOP:
        return fabs( coord.z() - adjGridNode[index]->coord.z() );
    
    case O_BOTTOM:
        return fabs( coord.z() - adjGridNode[index]->coord.z() );

    default:
        PRINT_INFO("WARNING: Adjacent grid node index doesn't exist");
        return 1.0;
  }

}

double CubitOctreeNode::manhattan_distance_adj_node( CubitOctreeNode *ptr_adj_node ){  
  return fabs( coord.x() - ptr_adj_node->coord.x() ) + fabs(coord.y() - ptr_adj_node->coord.y()) + fabs(coord.z() - ptr_adj_node->coord.z());
}

int CubitOctreeNode::get_adj_node_distance( enum OctreePosition select){
  
  switch( select ){
    
    case O_LEFT:
        return adjNodeDistance[O_LEFT];
    
    case O_RIGHT:
        return adjNodeDistance[O_RIGHT];
    
    case O_BOTTOM:
        return adjNodeDistance[O_BOTTOM];
    
    case O_TOP:
        return adjNodeDistance[O_TOP];     
    
    case O_BACK:
        return adjNodeDistance[O_BACK];
    
    case O_FRONT:
        return adjNodeDistance[O_FRONT];
    
    default:
        PRINT_INFO("ERROR: AdjNode index exceeded");
        return CUBIT_TRUE;
    
  }
  
}

int CubitOctreeNode::find_min_depth_cell_and_depth_difference( void ){
  
  int i, j, k;
  int depth_min_depth_cell = CUBIT_INT_MAX;
  int depth_max_detph_cell = -CUBIT_INT_MAX;
  int depth;

  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if( adjCell[i][j][k] != NULL ){
          depth = adjCell[i][j][k]->get_depth();
          if( depth < depth_min_depth_cell ){
            depth_min_depth_cell = depth;
            minDepthCell = adjCell[i][j][k];
          }
          if( depth > depth_max_detph_cell ){
            depth_max_detph_cell = depth;
          }
        }
      }
    }
  }
  
  if( depth_min_depth_cell != CUBIT_INT_MAX && depth_max_detph_cell != -CUBIT_INT_MAX ){
    cellDepthDifference = depth_max_detph_cell - depth_min_depth_cell;
  }
  else{
    cellDepthDifference = 0;
  }

  return cellDepthDifference;
}

void CubitOctreeNode::calculate_size_based_on_cell_dimension( double bbox_dimension )
{
  int i;
  int counter = 0;
  double sum = 0.0;
  
  for( i = 0; i < 6; i++ ){
    
    if( adjNodeDistance[i] >= 0  ){
      sum += 1.0 / pow( 2.0, adjNodeDistance[i] );
      counter++;
    }
  }
  
  size = sum / counter * bbox_dimension;
}

int CubitOctreeNode::find_half_space( CubitFacet *ptr_facet )
{
  if( (coord - ptr_facet->point(0)->coordinates()) % (ptr_facet->normal()) >= 0 )
  {
    halfspaceDirection = OCTREE_POSITIVE;
  }
  
  return halfspaceDirection;
}



void CubitOctreeNode::calc_facet_patch_distance_normal(DLIList<OctreeIntersectionData*> &idatas, int num_use_idatas, double &patch_distance, CubitVector &patch_normal, CubitBoolean sort, CubitBoolean set_Refface)
{
  int i, j;
  CubitBoolean duplicate = CUBIT_FALSE;
  for (i=0; i < idatas.size(); ++i)
  {
    idatas.reset();
    idatas.step(i);
    OctreeIntersectionData *first = idatas.get();

    
    for (j=i+1; j < idatas.size(); ++j)
    {
      idatas.reset();
      idatas.step(j);
      if (idatas.get() == first) {duplicate = CUBIT_TRUE;}
      break;
    }
    if (duplicate == CUBIT_TRUE) {break;}
    
  }
  if (duplicate) {PRINT_INFO("Duplicate idatas found for node %d\n", this->get_num());}
  
  if (sort) {idatas.sort(OctreeIntersectionData::compare_function);}

  idatas.reset();
  if (sort)
  {
    patch_distance = /*(coord - idatas.get()->get_int_point())%idatas.get()->get_facet_normal();*/ idatas.get()->get_length();
//    CubitFacet *ptr_facet = idatas.get()->get_facet_ptr();
//    CubitVector closest_point_on_facet = idatas.get()->get_int_point(), facet_normal = idatas.get()->get_facet_normal();
    
/*    if (fabs((coord - closest_point_on_facet)%facet_normal) != (coord-closest_point_on_facet).length())
      {
      SVDrawTool::clear_non_retained();
      SVDrawTool::draw_vector(coord, closest_point_on_facet, CUBIT_YELLOW_INDEX);
      PRINT_INFO("closest length is %lf\n", idatas.get()->get_length());
      DLIList<CubitFacet*> temp_facets;
      temp_facets.append(ptr_facet);
      
      SVDrawTool::draw_facets(temp_facets, CUBIT_YELLOW_INDEX);
      temp_facets.clean_out();
      for (i=0; i < idatas.size(); ++i)
      {
      CubitFacet *facet = idatas.get()->get_facet_ptr();
      if (facet != ptr_facet)
      {
      SVDrawTool::draw_vector(coord, idatas.get()->get_int_point(), CUBIT_RED_INDEX);
      PRINT_INFO("Other length is %lf\n", idatas.get()->get_length());
      }
        
      idatas.step();
      temp_facets.append(facet);
      }
      temp_facets.remove(ptr_facet);
      SVDrawTool::draw_facets(temp_facets, CUBIT_RED_INDEX);
      int j,k;

      for (i=0; i < 2; ++i)
      {
      for (j=0; j < 2; ++j)
      {

      for (k=0; k < 2; ++k)
      {
            
      CubitOctreeCell *CubitOctree_cell = adjCell[i][j][k];
      if (CubitOctree_cell == NULL) {continue;}
      double corners[3];
      double half_edge_length = CubitOctree_cell->get_dimension()/2.0;
      CubitVector center = CubitOctree_cell->get_center();
      center.get_xyz(corners);
      float box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
      corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
      SVDrawTool::draw_cube(box, CUBIT_GREEN_INDEX, SVDrawTool::WIRE);
      }
      }
      }
      
      SVDrawTool::mouse_xforms();
      
      }*/
  }
  else {patch_distance = -1;}

  if (set_Refface) {refFace = idatas.get()->get_face();}
  
  if (num_use_idatas == -1) {num_use_idatas = idatas.size();}
  
  patch_normal = CubitVector(0,0,0);
    
    /*for (i=0; i < num_use_idatas; ++i)
      {
      double len = idatas.get()->get_length();
      if (len < OCTREE_EPSILON) {len = OCTREE_EPSILON;}
      normal += idatas.get()->get_facet_normal() * 1/(len);
      idatas.step();
      }
      normal.normalize(); */
  
  OctreeIntersectionData *ptr_data;
    //double denominator = 0.0;
  double length;
        
    //distance = 0.0;
  for( i=0; i < num_use_idatas; i++ )
  {
    ptr_data = idatas.get_and_step();  // VED: important

    length = ptr_data->get_length();
    if( length < OCTREE_EPSILON ){
      patch_normal += ptr_data->get_normal() * (1/(OCTREE_EPSILON*OCTREE_EPSILON));        
        //distance += ( ptr_data->get_normal()%( coord - ptr_data->get_int_point() ) ) * (1/OCTREE_EPSILON)/denominator;
    }
    else{
      patch_normal += ptr_data->get_normal() * (1/(length*length));        
        //distance += ( ptr_data->get_normal()%( coord - ptr_data->get_int_point() ) ) * (1/length)/denominator;
    }
  }
  patch_normal.normalize();
}


void CubitOctreeNode::SAT_find_face_distance_average_normal ()
{
    // three cases:
    // 1) One facet (one RefFace) => use facet's normal and distance to facet, RefFace is set to facet's owning face
    // 2) Multiple facets, one refFace => use distance to closest facet and its owning RefFace, IDW normal from N closest facets
    // 3) Multiple facets, multiple RefFaces => use distance to closest facet and its owning RefFace, IDW normal from N closest facets

    //int i, j, k;

    // case: no intersection datas, this should not happen. I should prolly put an cassertere.
  if (octreeIntersectionDataList.size() == 0)
  {
    PRINT_ERROR("No OctreeIntersectionDatas attached to black grid node in queue for MAT generation!\n");
    return;
  }

    // case: one facet, just get distance, normal, and face from facet
  if (octreeIntersectionDataList.size() == 1)
  {
    distance = octreeIntersectionDataList.get()->get_length();
    mNormal = octreeIntersectionDataList.get()->get_facet_normal();
    refFace = octreeIntersectionDataList.get()->get_face();
    
   
      // remember to delete the idata - check if this is ok
    delete octreeIntersectionDataList.get();
    octreeIntersectionDataList.clean_out();
    return;
  }

  

    // case: multiple facets, one face
    // just initialize normal, distance, and Refface

    // case: two faces
    // determine if disconnected and if normals diverge => EWC case 1
    // else check angle 
  
  

    // now find N closest facets (keep list of idatas though)
    // choose closest one, use distance to it and use its RefFace
    // then use all and IDW to get normal

  int num_facets_to_use;
  if (N_CLOSEST_FACETS_FACTOR_FOR_FRONT_NORMALS == 0.00) {num_facets_to_use = 1;}
  else if (N_CLOSEST_FACETS_FACTOR_FOR_FRONT_NORMALS == 1.00) {num_facets_to_use = octreeIntersectionDataList.size();}
  else {num_facets_to_use = (int)(N_CLOSEST_FACETS_FACTOR_FOR_FRONT_NORMALS * octreeIntersectionDataList.size());}
  if (num_facets_to_use == 0) {num_facets_to_use = 1;}
  else if (num_facets_to_use > octreeIntersectionDataList.size()) {num_facets_to_use = octreeIntersectionDataList.size();}

    /*DLIList<OctreeIntersectionData*> n_closest_idatas;

    RTree<OctreeIntersectionData*> *rtree = new RTree<OctreeIntersectionData*>;
    double closest = CUBIT_DBL_MAX;
    for (i=0; i < octreeIntersectionDataList.size(); ++i)
    {
    rtree->add(octreeIntersectionDataList.get_and_step());
    }
    rtree->k_nearest_neighbor(coord, num_facets_to_use, closest, n_closest_idatas, OctreeIntersectionData::dist_sqr_to_vec);
    delete rtree;

    normal = CubitVector(0,0,0);
//  calc_facet_patch_distance_normal(octreeIntersectionDataList, num_facets_to_use, distance, normal, CUBIT_TRUE, CUBIT_TRUE);
  
  
refFace = n_closest_idatas.get()->get_face();


distance = n_closest_idatas.get()->get_length();

for (i=0; i < num_facets_to_use; ++i)
{
double len = n_closest_idatas.get()->get_length();
normal += n_closest_idatas.get_and_step()->get_facet_normal() * 1/(len);
}

normal.normalize();*/


    mNormal = CubitVector(0,0,0);
    calc_facet_patch_distance_normal(octreeIntersectionDataList, num_facets_to_use, distance, mNormal, CUBIT_TRUE, CUBIT_TRUE);
  
      /*
        octreeIntersectionDataList.sort(OctreeIntersectionData::compare_function);
        octreeIntersectionDataList.reset();
  

  
  
        refFace = octreeIntersectionDataList.get()->get_face();


        distance = octreeIntersectionDataList.get()->get_length();

        for (i=0; i < num_facets_to_use; ++i)
        {
        double len = octreeIntersectionDataList.get()->get_length();
        normal += octreeIntersectionDataList.get_and_step()->get_facet_normal() * 1/(len);
        }

        normal.normalize();
      */
      // now normal, distance to boundary, and Refface have been set
  
    }


// checks intesection between the lines joining grid node and adjacent nodes with the facet.

void CubitOctreeNode::find_intersection_with_facet( CubitOctreeType type, RefFace *ptr_face, CubitFacet *ptr_facet, DLIList<CubitOctreeNode*> &boundary_white_node_list )
{
  
  int i;
  CubitBoolean result;
  CubitVector int_point;
  CubitVector facet_normal;
  double para;
  
  facet_normal = ptr_facet->normal();
  
  for( i = 0; i < 6; i++ ){
    
    result = CUBIT_FALSE;
    
    if( adjGridNode[i] != NULL ){
      if( adjGridNode[i]->halfspaceDirection == OCTREE_NEGATIVE && adjGridNode[i]->mark == CUBIT_TRUE ){
          //-  adj_node[0] = O_FRONT Node
          //-  adj_node[1] = O_BACK Node
          //-  adj_node[2] = O_RIGHT Node
          //-  adj_node[3]  = O_LEFT Node
          //-  adj_node[4] = O_TOP Node
          //-  adj_node[5] = O_BOTTOM Node

        switch( i ){
          
          case O_FRONT:
              result = find_intersection_point( OCTREE_X, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          case O_BACK:
              result = find_intersection_point( OCTREE_X, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          
          case O_RIGHT:
              result = find_intersection_point( OCTREE_Y, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          
          case O_LEFT:
              result = find_intersection_point( OCTREE_Y, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          
          case O_TOP:
              result = find_intersection_point( OCTREE_Z, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          
          case O_BOTTOM:
              result = find_intersection_point( OCTREE_Z, coord, adjGridNode[i]->coord, facet_normal, ptr_facet->point(0)->coordinates(), ptr_facet->point(1)->coordinates(), ptr_facet->point(2)->coordinates(), int_point, para );
              break;
          
          default:
              break;
          
        }
      }
    }
   
    
    if( result == CUBIT_TRUE ){     
      OctreeIntersectionData * ptr_data;
            
      switch( type ){
        case CUBIT_OCTREE_VOLUME:                
              // update color of boundary node
            if( color != CUBIT_WHITE_INDEX ){
              if( color ==  CUBIT_BLACK_INDEX ){
                  // Mark the common cells incident on both end points of the edge as GREY
                  // To resolve the intersection problems I added the SAT intersection code -ved
              }
          
              boundary_white_node_list.push( this );
              refFace = ptr_face;
                //PRINT_DEBUG_157(" Testing:  Face Num = %d\n", RefFace->id() );
              color = CUBIT_WHITE_INDEX;
              distance = -1;
            }
            // WARNING: devided by zero.
             ptr_data = new OctreeIntersectionData( this, facet_normal * -1, int_point, (adjGridNode[i]->coord - int_point).length(), ptr_face );
            adjGridNode[i]->append_list_item( ptr_data );
            break;

        case CUBIT_OCTREE_FACE:
            refFace = ptr_face;
            color = CUBIT_BLACK_INDEX;
            distance = -1;
            adjGridNode[i]->color = CUBIT_BLACK_INDEX;
            break;

        default:
            PRINT_INFO("This case not yet implemented \n");
            break;
      }      
    }
  }
}





// returns true if the intersection between the line segment and facet takes place
// para stores the parameter of intersection point int_point
CubitBoolean CubitOctreeNode::find_intersection_point( int axis, CubitVector grid_node0, CubitVector grid_node1, CubitVector &facet_normal, CubitVector facet_vert0, CubitVector facet_vert1, CubitVector facet_vert2, CubitVector &int_point, double &para ){
  
  double A, B, C, D;
  double nominator, denominator;
  
  D = - ( facet_normal%(facet_vert1 - grid_node0) );
  A = facet_normal.x();
  B = facet_normal.y();
  C = facet_normal.z();
  
  switch( axis ){
    
    case OCTREE_X:
        if( fabs(A) < OCTREE_EPSILON ){
            //line parallel to facet 
            // both end points are intersection points
      
          return CUBIT_FALSE;
        }
        else{
          nominator = D;
          denominator = A * ( grid_node0.x() - grid_node1.x() );
          para = nominator / denominator;
        }
        break;
    
    
    case OCTREE_Y:
        if( fabs(B) < OCTREE_EPSILON ){
            //line parallel to facet 
            // both end points are intersection points
      
          return CUBIT_FALSE;
        }
        else{
          nominator = D;
          denominator = B * ( grid_node0.y() - grid_node1.y() );
          para = nominator / denominator;
        }
        break;
    
    
    case OCTREE_Z:
        if( fabs(C) < OCTREE_EPSILON ){
            //line parallel to facet 
            // both end points are intersection points
      
          return CUBIT_FALSE;
        }
        else{
          nominator = D;
          denominator = C * ( grid_node0.z() - grid_node1.z() );
          para = nominator / denominator;
        }
        break;
    
    default:
        break;
  }
  
  int_point = grid_node0 + para * ( grid_node1 - grid_node0 );
  return( is_intersection_point_contained_inside_facet( int_point, facet_vert0, facet_vert1, facet_vert2 ) );
  
  
}


CubitBoolean CubitOctreeNode::is_same_side(const CubitVector &p1, const CubitVector &p2, const CubitVector &a, const CubitVector &b)
{
  static CubitVector edge;
  edge = b-a;
    //cp1 = edge*(p1-a);
    //cp2 = edge*(p2-a);
  if ( (edge*(p1-a)) % (edge*(p2-a))  >= 0) {return CUBIT_TRUE;}
  return CUBIT_FALSE;
}


// returns true if intersection point is contained inside the facet
// interior angle at the intersecton point should add up to 360 deg.
CubitBoolean CubitOctreeNode::is_intersection_point_contained_inside_facet( const CubitVector &int_point, const CubitVector &facet_vert0, const CubitVector &facet_vert1, const CubitVector &facet_vert2 ){
  if ( fabs((int_point-facet_vert0)%( (facet_vert1-facet_vert0) * (facet_vert2-facet_vert0) )) > OCTREE_EPSILON)
  {
    return CUBIT_FALSE;
  }
  
  if (is_same_side(int_point,facet_vert0, facet_vert1,facet_vert2) && is_same_side(int_point,facet_vert1, facet_vert0,facet_vert2) && is_same_side(int_point,facet_vert2, facet_vert0,facet_vert1)) {return CUBIT_TRUE;}
  return CUBIT_FALSE;
  
    /*CubitVector line0, line1, line2;
      double ang0, ang1, ang2;
  
      line0 = facet_vert0 - int_point ;
      line0 = line0 / line0.length();
  
  
      line1 = facet_vert1 - int_point ;
      line1 = line1 / line1.length();
  
      line2 = facet_vert2 - int_point ;
      line2 = line2 / line2.length();
        //Added (M. Brewer) to check for an invalid number passed to acos
          //compute ang0
          double temp_double = line0 % line1;
          if(temp_double>1.0)
          temp_double=1.0;
          else if (temp_double<-1.0)
          temp_double=-1.0;
          ang0 = acos( temp_double  );
            //compute ang1
            temp_double = line1 % line2;
            if(temp_double>1.0)
            temp_double=1.0;
            else if (temp_double<-1.0)
            temp_double=-1.0;
            ang1 = acos( temp_double  );
              //compute ang2  
              temp_double = line2 % line0;
              if(temp_double>1.0)
              temp_double=1.0;
              else if (temp_double<-1.0)
              temp_double=-1.0;  
              ang2 = acos( temp_double  ); 
  
              if( fabs ( ( ang0 + ang1 + ang2 ) - CUBIT_PI * 2 ) < OCTREE_EPSILON ){
              return CUBIT_TRUE;
              }
              else{
              return CUBIT_FALSE;
              }*/
}





CubitBoolean CubitOctreeNode::find_size_using_adj_node()
{
  int i;
  double sum = 0.0;
  int count = 0;
  for( i = 0; i < 6; i ++ ){
    if( adjGridNode[i] != NULL ){
        // WARNING: Why we need to check for -1
        //if( adjGridNode[i]->size != -1 && adjGridNode[i]->size != CUBIT_DBL_MAX && adjGridNode[i]->size != 0){
        //if( adjGridNode[i]->size != CUBIT_DBL_MAX && adjGridNode[i]->size != 0.0 ){
      if( adjGridNode[i]->size != 0.0 ){
        sum += adjGridNode[i]->size;
        count++;
      }
    }
  }

  if( sum != 0.0 ){
    size = sum / count;
    return CUBIT_TRUE;
  }
  else
    return CUBIT_FALSE;
}



CubitBoolean CubitOctreeNode::compare_function( CubitOctreeNode *&a, CubitOctreeNode *&b ){

  if( a->get_distance() < b->get_distance() )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;

}

void CubitOctreeNode::find_distance_at_adj_node(PriorityQueue<CubitOctreeNode *> *heap )
{
  int i;
  CubitVector mat_pnt_center;
  double new_node_distance;
  
  // PRINT_INFO("Testing: Face id = %d\n", mrefFace->id());
  visit = CUBIT_TRUE;
  
  for( i = 0; i < 6; i++ ){
    
    // at the boundary of solid
    if( adjGridNode[i] == NULL )
      continue;
    
    
    //  *************** WARNING *********
    
    
    // if the adjGridNode is on boundary or outside don't do anything
    // Just testing white node is enough because we assume the facet and the grid
    // intersection done during octree generation is correct
    // At this point only at the boundary there is white and black nodes because of facet/octree intersection
    // and both at interior and outside we grey node
    //if( adjGridNode[i]->get_color() == CUBIT_WHITE_INDEX || adjGridNode[i]->get_color() == CUBIT_GREY_INDEX )
    if( adjGridNode[i]->get_color() == CUBIT_WHITE_INDEX || adjGridNode[i]->get_color() == CUBIT_YELLOW_INDEX)
      continue;
    
    // meeting of the front should be the major
    // criteria not the size value. check for null of
    // grid_node->face. then check for adjacency, black_white, and size...
    new_node_distance = distance + ( adjGridNode[i]->coord - coord )%mNormal;
    
    // WARNING: The below condition is very important.
    // Therefore I think the normal calculation is wrong
    // Front should always go forward.
    
    if(new_node_distance <= distance)
    {
      
      // OPEN IT
      // this happens at the concave region containing many surfaces.
      //PRINT_DEBUG_157(" Distance Calculation is wrong (discresing)\n");
      
      
      if (adjGridNode[i]->refFace != NULL)
      {
        //if (adjGridNode[i]->mrefFace != mrefFace)
        // {
        continue;
        //}
      }
      
      
        //if (adjGridNode[i]->mrefFace != mrefFace)
        //{
        adjGridNode[i]->color = CUBIT_BLACK_INDEX;
        //new_node_distance = distance;
        continue;
        //}
      
      
      //SVDrawTool::draw_vector(adjGridNode[i]->coord, closest, CUBIT_RED_INDEX);
      
      
    }
    //if( adjGridNode[i]->mrefFace != mrefFace )
    //{
    adjGridNode[i]->color = CUBIT_BLACK_INDEX;
    //}
    
    if( adjGridNode[i]->refFace == NULL ){
      // internal black nodes
      // not near skeleton
      adjGridNode[i]->distance = new_node_distance;
      adjGridNode[i]->refFace = refFace;
      adjGridNode[i]->mNormal = mNormal;
      heap->push( adjGridNode[i] );
    }
    else{
      
      
    }
  }
  
}





//EOF
