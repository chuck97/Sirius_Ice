//-------------------------------------------------------------------------
// Filename      : OctreeIntersectionData.hpp
//
// Purpose       : Data related to intersection between octree and facets
//
// Creator       : William Roshan Quadros
//
// Creation Date :  08/26/2013
//
// Owner         : 
//-------------------------------------------------------------------------
#ifndef OCTREEINTERSECTIONDATA_H
#define OCTREEINTERSECTIONDATA_H

#include "CubitVector.hpp"
#include "CubitDefines.h"


class RefFace;

class CubitFacet;
class CubitOctreeNode;

class OctreeIntersectionData{

public:
  OctreeIntersectionData(const CubitVector &normal, const CubitBoolean half_space, const double &len, RefFace *ptr_face, const CubitVector &closest_point_on_facet, CubitFacet *lp_facet); // SAT constructor
  OctreeIntersectionData( CubitOctreeNode *ptr_white_node, CubitVector normal, CubitVector int_point, double len, RefFace *ptr_face);
  ~OctreeIntersectionData(){}
  inline CubitVector get_normal(void){ return facetNormal; }
  inline CubitVector get_int_point(void){return intPoint; }
  inline double get_length(void){ return length; }
  inline RefFace * get_face(void){ return refFace; }
  inline CubitOctreeNode *get_white_node(void){ return whiteNode; }
  inline CubitFacet* get_facet_ptr() {return ptrFacet;} // for SAT code
  inline CubitBoolean get_halfspace() {return halfSpace;} // for SAT code
  inline CubitVector get_facet_normal() {return facetNormal;}

  inline static int compare_function(OctreeIntersectionData *&a, OctreeIntersectionData *&b)
  {
    if (a->get_length() < b->get_length()) {return -1;}
    else if (a->get_length() > b->get_length()) {return 1;}
    else {return 0;}
  }

  inline CubitBoolean is_merged() {return merged;}
  inline void set_merged(const CubitBoolean val) {merged = val;}

 private:
  CubitOctreeNode *whiteNode;
  CubitVector facetNormal;
  CubitVector intPoint;
  double length;
  RefFace *refFace;
  CubitFacet *ptrFacet;
  CubitBoolean halfSpace;
  CubitBoolean merged;

   
};
#endif

//EOF

