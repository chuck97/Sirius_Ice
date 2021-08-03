#include "CubitOctreeCell.hpp"

#include "CubitColor.hpp"
#include "CubitOctree.hpp"
#include "CubitOctreeConstants.hpp"
#include "CubitOctreeNode.hpp"
#include "CubitPoint.hpp"
#include "CubitBox.hpp" 
#include "GfxDebug.hpp" 
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "OctreeFacetPointData.hpp"
#include "OctreeIntersectionData.hpp"

// temporary includes
#include "CastTo.hpp"
#include "RefFace.hpp"
#include "RefFace.hpp"



/* -------------------- Methods of CubitOctreeCell ----------------- */

CubitOctreeCell::CubitOctreeCell( CubitVector cen, double dim, int level, CubitOctreeCell *parent_cell ){
  static int counter = 0; 
  int i, j, k; 
  
  leaf = CUBIT_TRUE;
  dimension = dim;
  mCenter = cen;
  depth = level;
  num = counter;
  counter ++;
  parent = parent_cell; 
  mark = CUBIT_FALSE;
  visit = CUBIT_FALSE;
  
  color = CUBIT_WHITE_INDEX;

//  CubitOctreeNode = new CubitOctreeNode* [2][2][2];

  myFacetList = NULL;
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        cubitOctreeNode[i][j][k] = NULL;
      }
    }
  }
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        children[i][j][k] = NULL;
      }
    }
  }
  
    //CubitOctree = NULL;
}

CubitOctreeCell::~CubitOctreeCell(){
  
  int i, j, k;
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if( children[i][j][k] != NULL ){
          delete children[i][j][k];
        }
      }
    }
  }

  int num_elm = octreeFacetPointDataList.size();
  for( i = 0; i < num_elm; i++ ){
    delete octreeFacetPointDataList.pop(); 
  }

    // delete_local_node_matrix();
  
    /*if (myFacetList)
      {
      for (i=0; i < myFacetList->size(); ++i)
      {
      delete myFacetList->get_and_step();
      }
      delete myFacetList;
      }*/ //uncomment later
  
}

void CubitOctreeCell::display_color_wireframe( CubitOctree *ptr_octree)
{
    
    //glColor3f( (cubitOctreeNode[0][0][0])->GetsizeColor(X), (cubitOctreeNode[0][0][0])->GetsizeColor(Y), (cubitOctreeNode[0][0][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][0][0])->GetsizeColor(X), (cubitOctreeNode[1][0][0])->GetsizeColor(Y), (cubitOctreeNode[1][0][0])->GetsizeColor(Z) );
    // GfxDebug::draw_line( cubitOctreeNode[0][0][0]->x(), cubitOctreeNode[0][0][0]->y(), cubitOctreeNode[0][0][0]->z(), cubitOctreeNode[1][0][0]->x(), cubitOctreeNode[1][0][0]->y(), cubitOctreeNode[1][0][0]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][0][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][0][0]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node( cubitOctreeNode[0][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node( cubitOctreeNode[1][0][0]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    
    //glColor3f( (cubitOctreeNode[0][1][0])->GetsizeColor(X), (cubitOctreeNode[0][1][0])->GetsizeColor(Y), (cubitOctreeNode[0][1][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][1][0])->GetsizeColor(X), (cubitOctreeNode[1][1][0])->GetsizeColor(Y), (cubitOctreeNode[1][1][0])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][1][0]->x(), cubitOctreeNode[0][1][0]->y(), cubitOctreeNode[0][1][0]->z(),cubitOctreeNode[1][1][0]->x(), cubitOctreeNode[1][1][0]->y(), cubitOctreeNode[1][1][0]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][1][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][1][0]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][1][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][1][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][1][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][0]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][0][0])->GetsizeColor(X), (cubitOctreeNode[0][0][0])->GetsizeColor(Y), (cubitOctreeNode[0][0][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[0][1][0])->GetsizeColor(X), (cubitOctreeNode[0][1][0])->GetsizeColor(Y), (cubitOctreeNode[0][1][0])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][0][0]->x(), cubitOctreeNode[0][0][0]->y(), cubitOctreeNode[0][0][0]->z(), cubitOctreeNode[0][1][0]->x(), cubitOctreeNode[0][1][0]->y(), cubitOctreeNode[0][1][0]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][0][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[0][1][0]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[0][1][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][1][0]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[1][0][0])->GetsizeColor(X), (cubitOctreeNode[1][0][0])->GetsizeColor(Y), (cubitOctreeNode[1][0][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][1][0])->GetsizeColor(X), (cubitOctreeNode[1][1][0])->GetsizeColor(Y), (cubitOctreeNode[1][1][0])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[1][0][0]->x(), cubitOctreeNode[1][0][0]->y(), cubitOctreeNode[1][0][0]->z(), cubitOctreeNode[1][1][0]->x(), cubitOctreeNode[1][1][0]->y(), cubitOctreeNode[1][1][0]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[1][0][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][1][0]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[1][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][1][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][0]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][0][1])->GetsizeColor(X), (cubitOctreeNode[0][0][1])->GetsizeColor(Y), (cubitOctreeNode[0][0][1])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][0][1])->GetsizeColor(X), (cubitOctreeNode[1][0][1])->GetsizeColor(Y), (cubitOctreeNode[1][0][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][0][1]->x(), cubitOctreeNode[0][0][1]->y(), cubitOctreeNode[0][0][1]->z(), cubitOctreeNode[1][0][1]->x(), cubitOctreeNode[1][0][1]->y(), cubitOctreeNode[1][0][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][0][1]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][0][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][1]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][0][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][1][1])->GetsizeColor(X), (cubitOctreeNode[0][1][1])->GetsizeColor(Y), (cubitOctreeNode[0][1][1])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][1][1])->GetsizeColor(X), (cubitOctreeNode[1][1][1])->GetsizeColor(Y), (cubitOctreeNode[1][1][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][1][1]->x(), cubitOctreeNode[0][1][1]->y(), cubitOctreeNode[0][1][1]->z(), cubitOctreeNode[1][1][1]->x(), cubitOctreeNode[1][1][1]->y(), cubitOctreeNode[1][1][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][1][1]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][1][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][1][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][1][1]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][1][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][0][1])->GetsizeColor(X), (cubitOctreeNode[0][0][1])->GetsizeColor(Y), (cubitOctreeNode[0][0][1])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[0][1][1])->GetsizeColor(X), (cubitOctreeNode[0][1][1])->GetsizeColor(Y), (cubitOctreeNode[0][1][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][0][1]->x(), cubitOctreeNode[0][0][1]->y(), cubitOctreeNode[0][0][1]->z(), cubitOctreeNode[0][1][1]->x(), cubitOctreeNode[0][1][1]->y(), cubitOctreeNode[0][1][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][0][1]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[0][1][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][1]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[0][1][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][1][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[1][0][1])->GetsizeColor(X), (cubitOctreeNode[1][0][1])->GetsizeColor(Y), (cubitOctreeNode[1][0][1])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][1][1])->GetsizeColor(X), (cubitOctreeNode[1][1][1])->GetsizeColor(Y), (cubitOctreeNode[1][1][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[1][0][1]->x(), cubitOctreeNode[1][0][1]->y(), cubitOctreeNode[1][0][1]->z(), cubitOctreeNode[1][1][1]->x(), cubitOctreeNode[1][1][1]->y(), cubitOctreeNode[1][1][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[1][0][1]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][1][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[1][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][0][1]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][1][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][0][0])->GetsizeColor(X), (cubitOctreeNode[0][0][0])->GetsizeColor(Y), (cubitOctreeNode[0][0][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[0][0][1])->GetsizeColor(X), (cubitOctreeNode[0][0][1])->GetsizeColor(Y), (cubitOctreeNode[0][0][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][0][0]->x(), cubitOctreeNode[0][0][0]->y(), cubitOctreeNode[0][0][0]->z(), cubitOctreeNode[0][0][1]->x(), cubitOctreeNode[0][0][1]->y(), cubitOctreeNode[0][0][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][0][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[0][0][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[0][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[0][1][0])->GetsizeColor(X), (cubitOctreeNode[0][1][0])->GetsizeColor(Y), (cubitOctreeNode[0][1][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[0][1][1])->GetsizeColor(X), (cubitOctreeNode[0][1][1])->GetsizeColor(Y), (cubitOctreeNode[0][1][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[0][1][0]->x(), cubitOctreeNode[0][1][0]->y(), cubitOctreeNode[0][1][0]->z(), cubitOctreeNode[0][1][1]->x(), cubitOctreeNode[0][1][1]->y(), cubitOctreeNode[0][1][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[0][1][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[0][1][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[0][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[1][0][0])->GetsizeColor(X), (cubitOctreeNode[1][0][0])->GetsizeColor(Y), (cubitOctreeNode[1][0][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][0][1])->GetsizeColor(X), (cubitOctreeNode[1][0][1])->GetsizeColor(Y), (cubitOctreeNode[1][0][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[1][0][0]->x(), cubitOctreeNode[1][0][0]->y(), cubitOctreeNode[1][0][0]->z(), cubitOctreeNode[1][0][1]->x(), cubitOctreeNode[1][0][1]->y(), cubitOctreeNode[1][0][1]->z(), CUBIT_BLUE_INDEX );
    //SkeletonDebug::draw_line( cubitOctreeNode[1][0][0]->get_coord(), cubitOctreeNode[1][0][0]->get_size( OCTREE_SIZE_DEFAULT), cubitOctreeNode[1][0][1]->get_coord(), cubitOctreeNode[1][0][1]->get_size( OCTREE_SIZE_DEFAULT) );
    if( cubitOctreeNode[1][0][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][0][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[0][0][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[0][0][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[0][0][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
    //glColor3f( (cubitOctreeNode[1][1][0])->GetsizeColor(X), (cubitOctreeNode[1][1][0])->GetsizeColor(Y), (cubitOctreeNode[1][1][0])->GetsizeColor(Z) );
    //glColor3f( (cubitOctreeNode[1][1][1])->GetsizeColor(X), (cubitOctreeNode[1][1][1])->GetsizeColor(Y), (cubitOctreeNode[1][1][1])->GetsizeColor(Z) );
    //GfxDebug::draw_line( cubitOctreeNode[1][1][0]->x(), cubitOctreeNode[1][1][0]->y(), cubitOctreeNode[1][1][0]->z(), cubitOctreeNode[1][1][1]->x(), cubitOctreeNode[1][1][1]->y(), cubitOctreeNode[1][1][1]->z(), CUBIT_BLUE_INDEX );
    if( cubitOctreeNode[1][1][0]->get_color() != CUBIT_GREY_INDEX && cubitOctreeNode[1][1][1]->get_color() != CUBIT_GREY_INDEX ){
        //SkeletonDebug::draw_line( cubitOctreeNode[1][1][0]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][0]->get_size( OCTREE_SIZE_DEFAULT)), cubitOctreeNode[1][1][1]->get_coord(), ptr_octree->get_scaled_from_wrld_size_grid_node(cubitOctreeNode[1][1][1]->get_size( OCTREE_SIZE_DEFAULT)) );
    }
    
}

void CubitOctreeCell::display( CubitOctree *ptr_octree, int opt ){
  
  int i, j, k;
  
  if( is_leaf() == CUBIT_FALSE ){
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
            //if( children[i][j][k] != NULL ){
          children[i][j][k]->display(ptr_octree, opt);
            //}
        }
      }
    }
  }
  else if (opt == 1)
  {
      // draw all cells in black wireframe mode
    //double corners[3];
    //double half_edge_length = get_dimension()/2.0;
    //CubitVector center = get_center();
    //center.get_xyz(corners);
    //double box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
                    //corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
    //SVDrawTool::draw_cube(box, CUBIT_BLACK_INDEX, SVDrawTool::WIRE);
    return;
  }
  else if (opt == 2)
  {
      // draw grey and black cells in black wireframe mode
      //(color == CUBIT_GREY_INDEX || color == CUBIT_BLACK_INDEX)
    if (color == CUBIT_BLACK_INDEX)
    {
      //double corners[3];
      //double half_edge_length = get_dimension()/2.0;
      //CubitVector center = get_center();
      //center.get_xyz(corners);
      //double box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
                      //corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
      //SVDrawTool::draw_cube(box, CUBIT_BLACK_INDEX, SVDrawTool::WIRE);
    }
    return;
  }
  else if (opt == 3)
  {
      // draw grey cells in white smoothshade mode
    if (color == CUBIT_GREY_INDEX)
    {
      //double corners[3];
      //double half_edge_length = get_dimension()/2.0;
      //CubitVector center = get_center();
      //center.get_xyz(corners);
      //double box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
                      //corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
      //SVDrawTool::draw_cube(box, CUBIT_WHITE_INDEX, SVDrawTool::SHADE);
      //SVDrawTool::draw_cube(box, CUBIT_BLACK_INDEX, SVDrawTool::WIRE);
    }
    return;
  }
  else if (opt == 4)
  {
    if (color != CUBIT_WHITE_INDEX)
    {
      //double corners[3];
      //double half_edge_length = get_dimension()/2.0;
      //CubitVector center = get_center();
      //center.get_xyz(corners);
      //double box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
                      //corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
      //SVDrawTool::draw_cube(box, CUBIT_WHITE_INDEX, SVDrawTool::SHADE);
      //SVDrawTool::draw_cube(box, CUBIT_BLACK_INDEX, SVDrawTool::WIRE);
    }
    return;
  }
  else if (opt == 5)
  {
      // draw grey and black cells in color wireframe mode
      if( color == CUBIT_GREY_INDEX || color == CUBIT_BLACK_INDEX )
          display_color_wireframe(ptr_octree);
      return;
  }
  else
    if( /*color == CUBIT_GREY_INDEX ||*/ color == CUBIT_BLACK_INDEX )
        display_color_wireframe(ptr_octree);
    return;
}

void CubitOctreeCell::release_octreefacetpointdata( void ){

  int i, j, k;
  
  if( is_leaf() == CUBIT_FALSE ){
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          children[i][j][k]->release_octreefacetpointdata( );
        }
      }
    }
  }
  else{
    int num_elm = octreeFacetPointDataList.size();
    for( i = 0; i < num_elm; i++ ){
      delete octreeFacetPointDataList.pop(); 
    }
  }
}

 

void CubitOctreeCell::display_octreefacetpointdata( void ){

  int i, j, k;
  
  if( is_leaf() == CUBIT_FALSE ){
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          children[i][j][k]->display_octreefacetpointdata( );
        }
      }
    }
  }
  else{
    for( i = 0; i < octreeFacetPointDataList.size(); i++ ){
      octreeFacetPointDataList.get_and_step()->display(); 
    }
  }
}

void CubitOctreeCell::set_oct_grid_node( const int i, const int j, const int k, CubitOctreeNode *ptr_node ){
  
  cubitOctreeNode[i][j][k] = ptr_node;
}

void CubitOctreeCell::set_child( int i, int j, int k, CubitOctreeCell *ptr_child_cell ){
  
  children[i][j][k] = ptr_child_cell;
  
}

CubitOctreeCell * CubitOctreeCell::get_child( const int i, const int j,const  int k ){
  
  return ( children[i][j][k] );
  
}

CubitOctreeCell * CubitOctreeCell::find_leaf_octree_cell(const CubitVector &point ){
  
  int i, j, k;
  
  if( is_leaf() == CUBIT_TRUE ){
    return this; 
  }
  else{
    
    i = j = k = 1;
    
    if( point.x() < mCenter.x() )
      i = 0; 
    
    if( point.y() < mCenter.y() )
      j = 0;
    
    if( point.z() < mCenter.z() )
      k = 0; 

    if (children[i][j][k] == NULL) {return NULL;}
    
    return children[i][j][k]->find_leaf_octree_cell( point );
    
  }
}

void CubitOctreeCell::distribute_facet_points_among_children( void ){
  int i;
  int l, m, n;
  
  OctreeFacetPointData *ptr_facet_point_data;
  
  octreeFacetPointDataList.reset();
  
  for( i = 0; i < octreeFacetPointDataList.size(); i++ ){
    
    ptr_facet_point_data = octreeFacetPointDataList.get_and_step();
 
    l = m = n = 1;
    
    if( ptr_facet_point_data->x() < mCenter.x() )
      l = 0; 
    
    if( ptr_facet_point_data->y() < mCenter.y() )
      m = 0;
    
    if( ptr_facet_point_data->z() < mCenter.z() )
      n = 0; 
    
    children[l][m][n]->append_list_item( ptr_facet_point_data );
    
  }
    // Clean the points and face list in the parent cell
  octreeFacetPointDataList.clean_out();
  
}


/*
  bool CubitOctreeCell::append_list_item( RefFace *ptr_Ref_face ){
  facetPointFaceList.push( ptr_Ref_face);
  return CUBIT_TRUE;
  }
*/
bool CubitOctreeCell::append_list_item( OctreeFacetPointData *ptr_facet_point_data ){
  octreeFacetPointDataList.push( ptr_facet_point_data );
  return CUBIT_TRUE;
}


bool CubitOctreeCell::add_adjacent_unmarked_cells( DLIList< CubitOctreeCell *> &queue ){
//optimized  
  static int i, j, k, l, m, n;
  static CubitOctreeNode *ptr_grid_node;
  static CubitOctreeCell *ptr_cell;
  static DLIList<CubitOctreeCell*> temp_list;
  
  for( i = 1; i >= 0; i-- ){
    for( j = 1; j >= 0; j-- ){
      for( k = 1; k >= 0; k-- ){
	
        ptr_grid_node = cubitOctreeNode[i][j][k];
	      
        for( l = 0; l < 2; l++ ){
          for( m = 0; m < 2; m++ ){
            for( n = 0; n < 2; n++ ){

              ptr_cell = ptr_grid_node->get_adj_cell( l, m, n );

              if( ptr_cell != NULL )
              {
                if( ptr_cell->mark == CUBIT_FALSE )
                {
                    //queue.append( ptr_cell );
                  temp_list.append(ptr_cell);
                    // ptr_cell->set_mark( CUBIT_TRUE );
                  ptr_cell->mark = CUBIT_TRUE;
                                  
                }
              }
              
            }
          }
        }
      }
    }
  }
  queue += temp_list;
  temp_list.clean_out();
  
  return CUBIT_TRUE;
}

bool CubitOctreeCell::is_intersects_box(const CubitBox &box ){
  
  double x_min, x_max, y_min, y_max, z_min, z_max;
  CubitVector box_min, box_max;
  int I, II, III, IV;
  
    //GfxDebug::draw_point( center.x(), center.y(), center.z(), CUBIT_RED_INDEX, itoa( num, buffer, 10 ) );
  
    // DEBUG:
    // changed from dimension -> dimenstion_2 
  double dimension_2 = dimension / 2.0;
  x_min = mCenter.x() - dimension_2;
  x_max = mCenter.x() + dimension_2;
  y_min = mCenter.y() - dimension_2;
  y_max = mCenter.y() + dimension_2;
  z_min = mCenter.z() - dimension_2;
  z_max = mCenter.z() + dimension_2;
  
  box_min = box.minimum();
  box_max = box.maximum();
  
    //PRINT_DEBUG_157("Testing: Bbox of facet MIN = [ %2.2lf %2.2lf %2.2lf ] MAX = [ %2.2lf %2.2lf %2.2lf ]\n",box_min.x(), box_min.y(), box_min.z(), box_max.x(), box_max.y(), box_max.z() );
    //PRINT_DEBUG_157("Testing:         Cell  MIN = [ %2.2lf %2.2lf %2.2lf ] MAX = [ %2.2lf %2.2lf %2.2lf ]\n",x_min, y_min, z_min, x_max, y_max, z_max );
  
    // A represent min of cell and B represent max of cell
    // C represent min of bbox and D represent max of bbox
    // A < C => I
    // B < C => II
    // A < D => III
    // B < D => IV
  
  I = II = III = IV = -1;
  if( x_max < box_min.x() ){
    II = 1;	
    I = 1;
  }
  else{
    if( x_min < box_min.x() ){
      I = 1;
    }
  }
  
  if( x_max <= box_max.x() ){
    III = 1;	
    IV = 1;
  }
  else{
    if( x_min <= box_max.x() ){
      III = 1;
    }
  }
  
  if( I * II == 1 && III * IV == 1 && I * III == 1 ){
      //PRINT_DEBUG_157(" X No Intersection \n");
    return CUBIT_FALSE;
  }
  
  
  
  I = II = III = IV = -1;
  if( y_max < box_min.y() ){
    II = 1;	
    I = 1;
  }
  else{
    if( y_min < box_min.y() ){
      I = 1;
    }
  }
  
  if( y_max <= box_max.y() ){
    III = 1;	
    IV = 1;
  }
  else{
    if( y_min <= box_max.y() ){
      III = 1;
    }
  }
  
  if( I * II == 1 && III * IV == 1 && I * III == 1 ){
      //PRINT_DEBUG_157(" Y No Intersection \n");
    return CUBIT_FALSE;
  }
  
  
  I = II = III = IV = -1;
  if( z_max < box_min.z() ){
    II = 1;	
    I = 1;
  }
  else{
    if( z_min < box_min.z() ){
      I = 1;
    }
  }
  
  if( z_max <= box_max.z() ){
    III = 1;	
    IV = 1;
  }
  else{
    if( z_min <= box_max.z() ){
      III = 1;
    }
  }
  
  if( I * II == 1 && III * IV == 1 && I * III == 1 ){
      //PRINT_DEBUG_157(" Z No Intersection \n");
    return CUBIT_FALSE;
  }
  
    //PRINT_DEBUG_157(" *********************** Intersection ************************* \n");
  return CUBIT_TRUE;
  
}


CubitBoolean CubitOctreeCell::is_facet_point_data_present( const CubitVector &coord ){
  
  OctreeFacetPointData *facet_point_data;
  int i;
  
  octreeFacetPointDataList.reset();
  for( i = 0; i < octreeFacetPointDataList.size(); i++ ){
    facet_point_data = octreeFacetPointDataList.get_and_step();
    if( (facet_point_data->coordinates() - coord ).length() < OCTREE_EPSILON ){
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
  
}

CubitBoolean CubitOctreeCell::is_facet_point_data_present( OctreeFacetPointData *new_facet_point_data ){
  
  OctreeFacetPointData *facet_point_data;
  int i;
  
  octreeFacetPointDataList.reset();
  for( i = 0; i < octreeFacetPointDataList.size(); i++ ){
    facet_point_data = octreeFacetPointDataList.get_and_step();
    if( facet_point_data->id() ==
        new_facet_point_data->id() ){
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
  
}

void CubitOctreeCell::coloring( /*DLIList<CubitOctreeCell *> &grey_cell_list ,*/ DLIList<CubitOctreeCell*> &black_cell_list )
{
  int i,j,k;
  int grey, black;

  grey = FALSE;
  black = TRUE;

  if( leaf == CUBIT_TRUE ){

    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          if( cubitOctreeNode[i][j][k]->get_color() != CUBIT_BLACK_INDEX )
          {
            black = FALSE;
          }
          else
          {
            grey = TRUE;
          }
        }
      }
    }

    if( black == TRUE )
    {
      color = CUBIT_BLACK_INDEX;
      black_cell_list.append(this);
    }
    else
      if( grey == TRUE ){
      color = CUBIT_GREY_INDEX;
      //grey_cell_list.push( this );
      }
      else
      {
        color = CUBIT_WHITE_INDEX; //constructor
      //return CUBIT_TRUE;
      }
    
  }
  else{
    for( i = 0; i < 2; i++ )
    {
      for( j = 0; j < 2; j++ )
      {
        for( k = 0; k < 2; k++ )
        {
          children[i][j][k]->coloring( /*grey_cell_list,*/ black_cell_list);
        }
      }
    }
  }

  //return CUBIT_FALSE;
  
}


CubitBoolean CubitOctreeCell::interpolate_grey_octreecell_node( void ){

  float avg_size;
  int counter;
  int i,j,k;

  avg_size = 0;
  counter = 0;
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
          //if( cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT ) != 0 && cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT ) != CUBIT_DBL_MAX  ){
          //if( cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT ) != 0.0 ){
        if (cubitOctreeNode[i][j][k]->get_color() == CUBIT_BLACK_INDEX)
        {
          avg_size += cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT );
          counter++;
        }
          //}
      }
    }
  }

  if( counter == 0 ){
    return CUBIT_FALSE;
  } 
  
  avg_size /= counter;

  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
          //if( cubitOctreeNode[i][j][k]->get_color() == CUBIT_WHITE_INDEX ){
        if( cubitOctreeNode[i][j][k]->get_color() != CUBIT_BLACK_INDEX ){
            //if( cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT ) == 0.0 ){
          cubitOctreeNode[i][j][k]->set_size( avg_size, OCTREE_SIZE_DEFAULT );
            //}
            //else{
            // cubitOctreeNode[i][j][k]->set_size( ( avg_size + cubitOctreeNode[i][j][k]->get_size( OCTREE_SIZE_DEFAULT ) ) / 2.0, OCTREE_SIZE_DEFAULT );
            //}					
        }
          //}
      }
    }
  }

  return CUBIT_TRUE;
}

// Trilienar interpolation	
// if number of nodes with zero size exceeds  MAX_NUM_ZERO_SIZE_NODES
// then average size is used instead of zero
// Even GREY grid nodes are consided 
// If a node has zero size then that node is not considered
// Therefore color is not important but the size
double CubitOctreeCell::trilinear_interpolation(const CubitVector &point ){
  double size;
  int i, j, k;
  CubitVector center, min_corner, max_corner;
  double dim, dim_2;
  double Nx[2], Ny[2], Nz[2];
  double size_at_node;
  CubitOctreeNode *ptr_node;

    // Find Shape functions at eight grid nodes of a cell
  center = get_center();
  dim = get_dimension();
  dim_2  =  dimension / 2.0;
		
  min_corner.x( center.x() - dim_2 );
  min_corner.y( center.y() - dim_2 );
  min_corner.z( center.z() - dim_2 );

  max_corner.x( center.x() + dim_2 );
  max_corner.y( center.y() + dim_2 );
  max_corner.z( center.z() + dim_2 );	
			
  Nx[0] = ( max_corner.x() - point.x() ) / dim;
  Ny[0] = ( max_corner.y() - point.y() ) / dim;
  Nz[0] = ( max_corner.z() - point.z() ) / dim;

  Nx[1] = ( point.x() - min_corner.x() ) / dim;
  Ny[1] = ( point.y() - min_corner.y() ) / dim;
  Nz[1] = ( point.z() - min_corner.z() ) / dim;

    // Trilinear interpolation	
    // if number of nodes with zero size exceeds  MAX_NUM_ZERO_SIZE_NODES
    // then average size is used instead of zero
	
		
    // Interpolate using non-zero nodes
  size = 0.0;

  int zero_node_counter = 0;

  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        ptr_node = get_octree_grid_node(i, j, k);
        size_at_node = ptr_node->get_size(OCTREE_SIZE_DEFAULT);
        if( size_at_node != 0.0 ){
          size +=  size_at_node * Nx[i] * Ny[j] * Nz[k];       	          
        }
        else{
          if( ptr_node->find_size_using_adj_node() == CUBIT_TRUE ){
            size += ptr_node->get_size( OCTREE_SIZE_DEFAULT ) * Nx[i] * Ny[j] * Nz[k];
          }
          else{
            zero_node_counter ++;
          }
        }
      }
    }
  }
  
    // For interior of volume avg_size will not be calculated
    // as no zero size node exist
  if( zero_node_counter > MAX_NUM_ZERO_SIZE_NODES ){

      // Find Average size
    double avg_size = 0.0;
    int counter = 0;

    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          size_at_node = get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT);
          if( size_at_node != 0.0 ){
            avg_size +=  size_at_node;
            counter ++;
          }					
        }
      }
    }

    if( counter > 0 ){
      avg_size /= counter;
			
        // use avg_size instead of zero during tri-linear interpolation
        // only at nodes with zero size
      for( i = 0; i < 2; i++ ){
        for( j = 0; j < 2; j++ ){
          for( k = 0; k < 2; k++ ){
            size_at_node = get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT);
            if( size_at_node == 0.0 ){
              size +=  avg_size * Nx[i] * Ny[j] * Nz[k];       	          
            }
          }
        }
      }

    }
    else{
      PRINT_DEBUG_157(" WARNING: Size of all nodes is zero \n" );
        //In size at a point default auto size will be used
        // autosize could be based on either average medial radius or
        // steve's autosize function
    }
  }
	
  return size;
}

  // Using Separating Axis Theorem (Eberly)
  // Mostly used Moller's article to implement this code

// returns true if intersection is present
CubitBoolean CubitOctreeCell::does_facet_intersect_octreecell(CubitFacet *ptr_facet)
{
  int i;
  
  CubitVector v0 = ptr_facet->point(0)->coordinates() - get_center(), v1 = ptr_facet->point(1)->coordinates() - get_center(),
      v2 = ptr_facet->point(2)->coordinates() - get_center();

  static CubitVector sep_axes[13];
  sep_axes[0] = ptr_facet->normal();
  sep_axes[1] = CubitVector(1,0,0);// note to self: optimize following trivial cross products
  sep_axes[2] = CubitVector(0,1,0);
  sep_axes[3] = CubitVector(0,0,1);
  sep_axes[4] = sep_axes[1]*(v1-v0);
  sep_axes[5] = sep_axes[1]*(v2-v1);
  sep_axes[6] = sep_axes[1]*(v0-v2);
  sep_axes[7] = sep_axes[2]*(v1-v0);
  sep_axes[8] = sep_axes[2]*(v2-v1);
  sep_axes[9] = sep_axes[2]*(v0-v2);
  sep_axes[10] = sep_axes[3]*(v1-v0);
  sep_axes[11] = sep_axes[3]*(v2-v1);
  sep_axes[12] = sep_axes[3]*(v0-v2);

  sep_axes[4].normalize();
  sep_axes[5].normalize();
  sep_axes[6].normalize();
  sep_axes[7].normalize();
  sep_axes[8].normalize();
  sep_axes[9].normalize();
  sep_axes[10].normalize();
  sep_axes[11].normalize();
  sep_axes[12].normalize();

  double rad_b, dot1, dot2, dot3, min_t, max_t;
  for (i=0; i < 13; ++i)
  {

      // rad_b is the "radius" of the projection of the cell onto the separating axis, centered about the projection of the cell's center
      // min_t and max_t are the extremes of the projection of the triangle onto the separating axis
    rad_b = (get_dimension() / 2.0) * (fabs(sep_axes[i].x()) + fabs(sep_axes[i].y()) + fabs(sep_axes[i].z()));
    dot1 = v0%sep_axes[i];
    dot2 = v1%sep_axes[i];
    dot3 = v2%sep_axes[i];
  
    min_t = CUBIT_MIN(CUBIT_MIN(dot1,dot2),dot3);
    max_t = CUBIT_MAX(CUBIT_MAX(dot1,dot2),dot3);

      // if the triangle's projection is on either side of the CubitOctree cell's projection, i.e. the projections don't overlap, then the two are disjoint
      //if ( min_t > rad_b || max_t < -rad_b)
    if ( (min_t - rad_b) > 1e-6 || (-rad_b -max_t) > 1e-6)
    {
      // found a separating axis so facet and cell do not intersect
      return CUBIT_FALSE;
    }
  }

    // OctreeIntersectionData code was here

  return CUBIT_TRUE;
}

void CubitOctreeCell::set_color_and_intersection_datas(CubitFacet *ptr_facet, RefFace *ptr_Ref_face,
#ifdef USE_octree_BGMESH
                                                  DLIList<CubitOctreeCell*> *greyCellList,
#endif
                                                  CubitSense surf_sense)
{
  int i, j, k;
  CubitVector facet_normal = ptr_facet->normal();
    //SVDrawTool::draw_vector(ptr_facet->center(),ptr_facet->center()+ptr_facet->normal(), CUBIT_MAGENTA_INDEX);
/*  if( ptr_facet->is_backwards() ) // Ved: note to self - check this
    {
    facet_normal = facet_normal * -1;
    }*/

    // if (!myFacetList) {myFacetList = new DLIList</*CubitFacet*/OctreeIntersectionData*>;} uncomment later

//    OctreeIntersectionData(CubitVector normal, int half_space, double len, RefFace *ptr_face, CubitVector closest_point_on_facet, CubitFacet *lp_facet);

/*  CubitVector dummy;
    OctreeIntersectionData *new_idata = new OctreeIntersectionData(dummy, 1, 0.0, ptr_Ref_face, dummy, ptr_facet);
    myFacetList->append(new_idata);*/ // uncomment later

    /* if (!myFaceList) {myFaceList = new DLIList<RefFace*>;}
       if (!myFaceList->where_is_item(ptr_Ref_face))
       {
       myFaceList->append(ptr_Ref_face);
       }*/

    // bool merged = (CAST_TO(ptr_Ref_face, RefFace)->num_parent_ref_entities() == 2);

  for (i=0; i < 2; ++i)
  {
    for (j=0; j < 2; ++j)
    {
      for (k=0; k < 2; ++k)
      {
        CubitVector closest_point_on_facet;
        double dist_to_facet;
        if (defaultDistanceMetric == CAPSULE_DIST) {dist_to_facet = CubitOctree::capsule_distance_to_facet(cubitOctreeNode[i][j][k]->get_coord(), ptr_facet, closest_point_on_facet);}
        else if (defaultDistanceMetric == PROJECTED_DIST) {dist_to_facet = CubitOctree::capsule_distance_to_facet(cubitOctreeNode[i][j][k]->get_coord(), ptr_facet, closest_point_on_facet, CUBIT_TRUE);}
       
          /*if (dist_to_facet < 0)
            {
              //PRINT_INFO("Invalid distance: %lf !!!!\n", dist_to_facet);
                //SVDrawTool::draw_point(cubitOctreeNode[i][j][k]->get_coord(), CUBIT_MAGENTA_INDEX);
          
                  //double test_dist = CubitOctree::capsule_distance_to_facet(cubitOctreeNode[i][j][k]->get_coord(), ptr_facet, closest_point_on_facet);
                  }*/

        OctreeIntersectionData *idata;
        
        if (/*merged == CUBIT_TRUE*/ surf_sense == CUBIT_UNKNOWN)
        {
          CubitVector temp_facet_normal = facet_normal;
          
          if (cubitOctreeNode[i][j][k]->get_halfspace_direction() == OCTREE_NEGATIVE)
          {
            temp_facet_normal *= -1.0;
          }

          
//          SVDrawTool::draw_point(cubitOctreeNode[i][j][k]->get_coord(), CUBIT_MAGENTA_INDEX);
          
          idata = new OctreeIntersectionData(temp_facet_normal, OCTREE_NEGATIVE, dist_to_facet, ptr_Ref_face, closest_point_on_facet, ptr_facet);
          idata->set_merged(1);
        }
        else
        {
          CubitVector temp_facet_normal = facet_normal;
          CubitBoolean temp_half_space = cubitOctreeNode[i][j][k]->get_halfspace_direction();
          if (surf_sense == CUBIT_FORWARD) {temp_facet_normal *= -1.0;}
          else if (surf_sense == CUBIT_REVERSED) {temp_half_space = !temp_half_space;}
          idata = new OctreeIntersectionData(temp_facet_normal, temp_half_space, dist_to_facet, ptr_Ref_face, closest_point_on_facet, ptr_facet);
        }
        
//        idata->set_facet_ptr(ptr_facet);
//        idata->set_int_point(closest_point_on_facet);
        cubitOctreeNode[i][j][k]->append_list_item(idata);
      }
    }
  }

  color = CUBIT_GREY_INDEX;
  
#ifdef USE_octree_BGMESH
  if (!visit)
  {
    greyCellList->append(this); // for background mesh generation
    visit = CUBIT_TRUE;
  }
#endif
  
}

CubitBoolean CubitOctreeCell::does_contain_positive_and_negative_nodes( void ){
  int i, j, k;
  CubitBoolean positive_space = CUBIT_FALSE;
  CubitBoolean negative_space = CUBIT_FALSE;
 
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if( cubitOctreeNode[i][j][k]->get_halfspace_direction() == OCTREE_POSITIVE )
          positive_space = CUBIT_TRUE;
        else
          negative_space = CUBIT_TRUE;
      }
    }
  }

  if( positive_space == CUBIT_TRUE && negative_space == CUBIT_TRUE ){
    return CUBIT_TRUE;
  }
  else{
    return CUBIT_FALSE;
  }
}

double CubitOctreeCell::inverse_distance_interpolation(const CubitVector &point ){
  double size;
  long double denominator=0;
  long double nominator=0;
  double weight;
  long double dist_func;
  CubitVector node_coord;
  double dist;
  int i, j, k;
  int type = DEFAULT_INVERSE_DISTANCE_SCHEME_INSIDE_CELL; 
  
  
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if(  get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT) != 0.0 ){
          node_coord = get_octree_grid_node(i, j, k)->get_coord();
	        
          dist = ( node_coord - point  ).length();
	        
          if( dist < OCTREE_EPSILON )
            dist = OCTREE_EPSILON;
	        
          switch( type ){
	          
            case INVERSE_LINEAR:
                dist_func = dist;
                break;
	          
            case INVERSE_QUADRATIC:
                dist_func = dist * dist;
                break;
	          
            case INVERSE_CUBIC:
                dist_func = dist * dist * dist;
                break;
	          
            default:
                dist_func = dist * dist;
          }
	        
          denominator += 1.0 / ( dist_func );
        }
      }
    }
  }
  
  size = 0.0;
  
  for( i = 0; i < 2; i ++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if(  get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT) != 0.0 ){
          node_coord = get_octree_grid_node(i,j,k)->get_coord();
	        
          dist = ( node_coord - point  ).length();
	        
          if( dist < OCTREE_EPSILON )
            dist = OCTREE_EPSILON;
	        
          switch( type ){
	          
            case INVERSE_LINEAR:
                dist_func = dist;
                break;
	          
            case INVERSE_QUADRATIC:
                dist_func = dist * dist;
                break;
	          
            case INVERSE_CUBIC:
                dist_func = dist * dist * dist;
                break;
            default:
                dist_func = dist * dist;
          }
	        
          nominator = 1.0 / ( dist_func );
	        
          weight = nominator / denominator ;
        
          size +=  get_octree_grid_node(i,j,k)->get_size( OCTREE_SIZE_DEFAULT ) * weight; 
            //PRINT_DEBUG_157(" size at grid node = %f\n",get_octree_grid_node(i,j,k)->get_size( OCTREE_SIZE_DEFAULT ));
        }	          
      }
    }
  } 

  return size;
}


double CubitOctreeCell::min_distance_interpolation(const CubitVector &point )
{
  double size;
  CubitVector node_coord;
  double dist;
  int i, j, k;
  
  double min_dist = CUBIT_DBL_MAX;
 
  size = 0.0;

  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if(  get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT) != 0.0 )
        {
          node_coord = get_octree_grid_node(i, j, k)->get_coord();
	        
          dist = ( node_coord - point  ).length();
	        
          if( dist < OCTREE_EPSILON )
            dist = OCTREE_EPSILON;
	        
          if( dist < min_dist )
          {
            min_dist = dist;
            size = get_octree_grid_node(i,j,k)->get_size( OCTREE_SIZE_DEFAULT );
          }          
        }
      }
    }
  }
  
  return size;
}



double CubitOctreeCell::min_size_interpolation(const CubitVector &point ){

  CubitVector node_coord;
  int i, j, k;
  
  double min_size = CUBIT_DBL_MAX;
  double size;
 
 

  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if(  get_octree_grid_node(i, j, k)->get_size(OCTREE_SIZE_DEFAULT) != 0.0 )
        {
          size = get_octree_grid_node(i,j,k)->get_size( OCTREE_SIZE_DEFAULT );
	        
          if( size < OCTREE_EPSILON )
            size = OCTREE_EPSILON;
	        
          if( size < min_size )
          {
            min_size = size;            
          }          
        }
      }
    }
  }
  
  return min_size;
}

CubitStatus CubitOctreeCell::find_indices_in_parent( int *index ){

  int i, j, k;
	
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        if( parent->get_child( i, j, k )->id() == num ){
          index[0] = i;
          index[1] = j;
          index[2] = k;
          return CUBIT_SUCCESS;
        }
      }
    }
  }

  return CUBIT_FAILURE;
}


#ifndef NDEBUG
void CubitOctreeCell::write_octreecell_sizing_info_file( FILE *pof, DLIList<CubitOctreeCell*> &stack ){
	
  int index[3];
  int i, j, k;

  if( parent == NULL ){
    fprintf( pof, "C %d -1 0 0 0 %d\n", num, !leaf );
  }
  else{
    find_indices_in_parent( index );
    fprintf( pof, "C %d %d %d %d %d %d \n", num, parent->id(), index[0], index[1], index[2], !leaf );
  }
	

  fprintf( pof, "V");
  for( i = 0; i < 2; i++ ){
    for( j = 0; j < 2; j++ ){
      for( k = 0; k < 2; k++ ){
        fprintf( pof," %d ", cubitOctreeNode[i][j][k]->get_num() );
      }
    }
  }	
  fprintf( pof, "\n");

  if( leaf == 0 ){  
    fprintf( pof, "H");
    for( i = 0; i < 2; i++ ){
      for( j = 0; j < 2; j++ ){
        for( k = 0; k < 2; k++ ){
          fprintf( pof," %d ", children[i][j][k]->id() );
          stack.append( children[i][j][k] );
        }
      }
    }	
    fprintf( pof, "\n");
  }
	
  fprintf( pof, "E\n");
}
#endif


//EOF

