#include "OctreeIntersectionData.hpp"
#include "RefFace.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"



OctreeIntersectionData::OctreeIntersectionData(const CubitVector &normal, const CubitBoolean half_space, const double &len, RefFace *ptr_face, const CubitVector &closest_point_on_facet, CubitFacet *lp_facet)
{
  facetNormal = normal;
  halfSpace = half_space;
  length = len;
  refFace = ptr_face;
  intPoint = closest_point_on_facet;
  ptrFacet = lp_facet;
  merged = 0;
}

OctreeIntersectionData::OctreeIntersectionData( CubitOctreeNode *ptr_white_node, CubitVector normal, CubitVector int_point, double len, RefFace *ptr_face)
{
  facetNormal = normal;
  intPoint = int_point;
  refFace = ptr_face;
  length = len; 
  whiteNode = ptr_white_node;
}


//EOF
