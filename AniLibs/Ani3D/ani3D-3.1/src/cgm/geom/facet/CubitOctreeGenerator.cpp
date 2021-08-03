#include "CubitOctreeGenerator.hpp"
#include "CubitOctree.hpp"
#include "CubitOctreeConstants.hpp"
  
CubitOctreeGenerator::CubitOctreeGenerator( ){
  }

CubitOctreeGenerator::~CubitOctreeGenerator(){
    delete cubitOctree;
}
 
void CubitOctreeGenerator::color_lattice_cell( void ){
  cubitOctree->color_octreecell();
  cubitOctree->set_status_octree_coloring(CUBIT_TRUE);
}

double CubitOctreeGenerator::size_at_a_point( const CubitVector &point )
{
  double size = 0;
  
  size = size_at_point_in_octree( point, MESH_SIZE );
  
  return size;
}

double CubitOctreeGenerator::size_at_point_in_octree( const CubitVector &point, int type  ){
  return cubitOctree->size_at_a_point( point, type );
}


CubitPointContainment CubitOctreeGenerator::point_containment( CubitVector tmp_vec, double tolerance )
{
  return cubitOctree->point_containment( tmp_vec, tolerance );
}


// EOF

