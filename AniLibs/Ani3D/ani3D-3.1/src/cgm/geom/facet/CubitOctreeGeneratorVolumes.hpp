//-------------------------------------------------------------------------
// Filename      :  CubitOctreeGeneratorVolumes.hpp
//
// Purpose       :  Upper level code for octree generation of  single volume or assembly
//
// Creator       :  William Roshan Quadros
//
// Creation Date :  08/28/2013
// Owner         : 
//-------------------------------------------------------------------------
#ifndef CUBIT_OCTREE_GENERATOR_VOLUMES_H
#define CUBIT_OCTREE_GENERATOR_VOLUMES_H

#include "CubitOctreeGenerator.hpp"
#include "CubitOctreeConstants.hpp"


template <class X> class DLIList;

class CubitOctree;
class CubitOctreeNode;
class RefVolume; 
class RefEntity;



class CubitOctreeGeneratorVolumes : public CubitOctreeGenerator
{
  
public:
  
    // Constructor for octree
  CubitOctreeGeneratorVolumes( DLIList<RefEntity*> &entity_list );
    
    //- Constructor for octree
  ~CubitOctreeGeneratorVolumes(){}
    //- Distructor
  
  CubitBoolean generate_lattice(void );  
  CubitBoolean find_intersection_between_octree_and_facets( DLIList<CubitOctreeNode *> &queue_for_mat_generation );
  CubitBoolean build_octree_till_skl_max_depth_based_on_facets(void);

  void get_bounding_box( CubitVector &min, CubitVector &max );


  // finds the optimal min and max depths of octree 
  static void find_optimal_min_and_max_octree_depths( DLIList< RefEntity *> &entity_list, int &min_depth, int &max_depth );
  static void find_optimal_min_and_max_depth( RefEntity *volume, int &local_min_depth, int &local_max_depth );

  
  void reset_td_octree_ref_face_visit( const CubitBoolean type );
  void reset_td_octree_ref_edge_visit( const CubitBoolean type );
  
  CubitBoolean build_td_mref_faces( void );
  CubitBoolean build_td_mref_edges( void );

  void color_octreenode_via_grassfire( DLIList<CubitOctreeNode *> &queue_for_mat_generation );
  CubitStatus generate_full_octree(void);
  
  
  void testing();
  
private:
  
  DLIList<RefEntity*> &entityList;
};

#endif

//EOF
