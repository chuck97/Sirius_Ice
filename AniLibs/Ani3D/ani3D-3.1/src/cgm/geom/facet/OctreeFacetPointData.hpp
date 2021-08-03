//-------------------------------------------------------------------------
// Filename      : OctreeFacetPointData.hpp
//
// Purpose       : Holds Intersection Data	
//
// Creator       : William Roshan Quadros
//
// Creation Date :  01/01/2003
//
// Owner         : 
//-------------------------------------------------------------------------
#ifndef OCTREEFACETPOINTDATA_H
#define OCTREEFACETPOINTDATA_H 

#include "CubitVector.hpp"
#include "DLIList.hpp"


class CubitFacet;
class CubitFacetEdge; 
class CubitPoint;


enum OctreeFacetPointDataType { CUBIT_FACET_POINT_DATA_FACET, CUBIT_FACET_POINT_DATA_EDGE, CUBIT_FACET_POINT_DATA_POINT };

class OctreeFacetPointData
{

public:
  OctreeFacetPointData( CubitVector coord, CubitFacet *ptr_facet );
  OctreeFacetPointData( CubitVector coord, CubitFacetEdge *ptr_facet_edge );
  OctreeFacetPointData( CubitVector coord, CubitPoint *ptr_point );
  ~OctreeFacetPointData(){}

  CubitVector coordinates( void ){ return xyz; }
  double x( void ){ return xyz.x(); }
  double y( void ){ return xyz.y(); }
  double z( void ){ return xyz.z(); }

  int calculate_id( void );
  int id( void ){ return num; }

  void display( void );

  static CubitBoolean generate_facet_point_data_at_slender_facet( CubitFacet *ptr_facet, DLIList<OctreeFacetPointData *> &facet_point_data_list ); 
  static CubitBoolean generate_facet_point_data_based_on_curvature( CubitFacetEdge *ptr_facet_edge, /*double angle,*/ DLIList<OctreeFacetPointData *> &facet_point_data_list );


private:
  
  union
  {
    CubitFacet *facetPtr;
    CubitFacetEdge *facetEdgePtr;
    CubitPoint *facetPointPtr;
  };

  int num;
  CubitVector xyz;
  OctreeFacetPointDataType type;
  
};

#endif

//EOF

