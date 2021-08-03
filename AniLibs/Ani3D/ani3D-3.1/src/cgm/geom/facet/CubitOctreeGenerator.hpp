//-------------------------------------------------------------------------
// Filename      :  CubitOctreeGenerator.hpp
//
// Purpose       :  Base class for octree generation of surfaces, single volume or assembly
//
// Creator       :  William Roshan Quadros
//
// Creation Date :  08/28/2013
// Owner         : 
//-------------------------------------------------------------------------
#ifndef CUBIT_OCTREE_GENERATOR_H
#define CUBIT_OCTREE_GENERATOR_H

#include "CubitVector.hpp"

template <class X> class DLIList;

class CubitOctree;
class CubitOctreeNode;


// Base class for source point
class CubitOctreeGenerator{
  
public:
   
  
  CubitOctreeGenerator( );
    //- Constructor for octree
  
  virtual ~CubitOctreeGenerator();
    //- Distructor

  CubitOctree* get_octree_lattice( void ){ return cubitOctree;}
 
  double size_at_point_in_octree( const CubitVector &point, int type );
  double size_at_a_point( const CubitVector &point );

  void color_lattice_cell( void );
  
  virtual void get_bounding_box( CubitVector &min, CubitVector &max ) = 0;

 
  virtual CubitBoolean generate_lattice( void ) = 0;
  
  CubitPointContainment point_containment( CubitVector tmp_vec, double tolerance );

 
protected:


  CubitOctree *cubitOctree;
    //- Octree model of the solid;


private:

};



#endif

//EOF
