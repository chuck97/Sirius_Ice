#include "OctreeFacetPointData.hpp"

#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
//#include "SVDrawTool.hpp"
#include "CubitColor.hpp"

/* -------------------- Methods of OctreeFacetPointData ----------------- */
int OctreeFacetPointData::calculate_id( void ){
  static int count = -1;
  count++;
  return count;
}

OctreeFacetPointData::OctreeFacetPointData( CubitVector coord, CubitFacet *ptr_facet ){
  xyz = coord;
  type = CUBIT_FACET_POINT_DATA_FACET;
  facetPtr = ptr_facet;
  num = calculate_id();
}

OctreeFacetPointData::OctreeFacetPointData( CubitVector coord, CubitFacetEdge *ptr_facet_edge ){
  xyz = coord;
  type = CUBIT_FACET_POINT_DATA_EDGE;
  facetEdgePtr = ptr_facet_edge;
  num = calculate_id();
}

OctreeFacetPointData::OctreeFacetPointData( CubitVector coord, CubitPoint *ptr_point ){
  xyz = coord;
  type = CUBIT_FACET_POINT_DATA_POINT;
  facetPointPtr = ptr_point;
  num = calculate_id();
}

CubitBoolean OctreeFacetPointData::generate_facet_point_data_based_on_curvature( CubitFacetEdge *ptr_facet_edge, /*double angle,*/ DLIList< OctreeFacetPointData *> &facet_point_data_list ){
  
  CubitFacet *ptr_adj_facet0 = ptr_facet_edge->adj_facet(0);
  CubitFacet *ptr_adj_facet1 = ptr_facet_edge->adj_facet(1);

    // Find minimum edge length
  int i;
  double edge_len, min_edge_len = CUBIT_DBL_MAX;
  for( i = 0; i < 3; i++ ){
    edge_len = ptr_adj_facet0->edge(i)->length();
    if( edge_len < min_edge_len ){
      min_edge_len = edge_len;
    }
  }
  for( i = 0; i < 3; i++ ){
    edge_len = ptr_adj_facet1->edge(i)->length();
    if( edge_len < min_edge_len ){
      min_edge_len = edge_len;
    }
  }
  
  int num_seg = (int)(ptr_facet_edge->length() / min_edge_len) + 1;

  // For triangles with high aspect ratio the number of segments can be enormous.  Limit
  // the number of segments.
  if(num_seg > 100)
    num_seg = 100;

  OctreeFacetPointData *facet_point_data;
  for( i = 1; i < num_seg; i++ ){
    facet_point_data = new OctreeFacetPointData( (ptr_facet_edge->point(0)->coordinates() * (num_seg - i) + ptr_facet_edge->point(1)->coordinates()  * (i)) / num_seg , ptr_facet_edge );
    facet_point_data_list.push( facet_point_data );
  }
  
  return CUBIT_TRUE;
}

void OctreeFacetPointData::display( void ){
  
  switch( type  ){
    case CUBIT_FACET_POINT_DATA_POINT:
        //SVDrawTool::draw_point( xyz, CUBIT_BLUE_INDEX );
        break;

    case CUBIT_FACET_POINT_DATA_EDGE:
        //SVDrawTool::draw_point( xyz, CUBIT_YELLOW_INDEX );
        break;
  
    case CUBIT_FACET_POINT_DATA_FACET:
        //SVDrawTool::draw_point( xyz, CUBIT_GREEN_INDEX );
        break;
    
    default:
        break;
  }

}



//EOF

