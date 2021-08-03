//- Class:       TDOctreeRefFace
//- Description: Tool data for storing additional information needed for MAT generation.
//- Owner:       W. R. Quadors
//- Checked by:
//- Version:
 
#include "TDOctreeRefFace.hpp"
#include "CubitDefines.h"
#include "ToolData.hpp"
#include "MemoryManager.hpp" 
#include "DLIList.hpp"
#include "CubitVector.hpp" 
#include "CastTo.hpp"

#include "CubitTransformMatrix.hpp"
#include "RefFace.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPointData.hpp"
#include "CubitOctreeConstants.hpp"
#include "FacetDataUtil.hpp"


// Constructor  
TDOctreeRefFace::TDOctreeRefFace(){

  ptrCubitPointList = NULL;
  ptrCubitFacetEdgeList = NULL;
  ptrCubitFacetList = NULL;
  visit = CUBIT_FALSE;
  create_2dmat = CUBIT_TRUE;
}

// Distructor
TDOctreeRefFace::~TDOctreeRefFace()
{
  if (ptrCubitFacetList != NULL)
  {
    FacetDataUtil::delete_facets(*ptrCubitFacetList);
    if( ptrCubitFacetList )
      delete ptrCubitFacetList;
    ptrCubitFacetList = NULL;

    if( ptrCubitFacetEdgeList )
      delete ptrCubitFacetEdgeList;
    ptrCubitFacetEdgeList = NULL;

    if( ptrCubitPointList )
      delete ptrCubitPointList;
    ptrCubitPointList = NULL;
  }
  
//PRINT_INFO("Inside ~TDOctreeRefFace\n");
}

//-------------------------------------------------------------------------
// Purpose       : To create TDOctreeRefFace to RefFace
//
// Special Notes :
//
// Creator       : W R Quadros
//
// Creation Date : 07/03
//------------------------------------------------------------------------- 
CubitStatus TDOctreeRefFace::add_td(  RefFace *ref_face )
{
  ToolData *td;
  td = ref_face->get_TD(&TDOctreeRefFace::is_td_octree_ref_face);
  if ( td == NULL )
  {
    TDOctreeRefFace *td_gm = new TDOctreeRefFace;
    ref_face->add_TD( td_gm);
    td_gm->set_ref_face( ref_face );
  }
  else
  {
    TDOctreeRefFace *td_gm = CAST_TO(td, TDOctreeRefFace);
    td_gm->set_ref_face( ref_face );

  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get the TDOctreeRefFace for a RefFace
//
// Special Notes :
//
// Creator       : W R Quadros
//
// Creation Date : 07/03
//------------------------------------------------------------------------- 
TDOctreeRefFace* TDOctreeRefFace::get_td( RefFace *ref_face )
{
  ToolData *td;
  td = ref_face->get_TD(&TDOctreeRefFace::is_td_octree_ref_face);
  if ( td != NULL )
  {
    TDOctreeRefFace *td_gm = CAST_TO(td, TDOctreeRefFace);
    return td_gm;
  }
  return (TDOctreeRefFace*) NULL;
} 
CubitBoolean TDOctreeRefFace::is_adj_curves( int id1, int id2 ){

  int loop_index1, loop_index2;
  int i;
  int first_curve_id, last_curve_id, start_curve_id;
		
  loopIndex.reset();
  loop_index1 = CUBIT_INT_MAX;
    // i indicates loop index
  for( i = 0; i < loopIndex.size(); i++ ){
    start_curve_id = loopIndex.get_and_step();
    if( id1 < start_curve_id ){
      loop_index1 = i - 1;
      break;
    }
  }
  if( loop_index1 == CUBIT_INT_MAX ){
    loop_index1 = i - 1;
  }
	
  loopIndex.reset();
  loop_index2 = CUBIT_INT_MAX;
    // i indicates loop index
  for( i = 0; i < loopIndex.size(); i++ ){
    start_curve_id = loopIndex.get_and_step();
    if( id2 < start_curve_id ){
      loop_index2 = i - 1;
      break;
    }
  }

  if( loop_index2 == CUBIT_INT_MAX ){
    loop_index2 = i - 1;
  }
	
  if( loop_index1 != loop_index2 ){
      // curves are not in same loop
    return CUBIT_FALSE;
  }
  else{
    if( abs(id1 - id2) == 1 ){
        // adj curves of same loop
      return CUBIT_TRUE;
    }
    loopIndex.reset();
    loopIndex.step( loop_index1 );
    first_curve_id = loopIndex.get_and_step();
    last_curve_id = loopIndex.get() - 1;

    if( last_curve_id < first_curve_id ){
//      PRINT_DEBUG_157(" id1 = %d and id2 = %d", id1, id2 );
//      PRINT_DEBUG_157(" Note: Last curve id ( %d ) is less than first curve id ( %d ) \n", last_curve_id, first_curve_id );
      last_curve_id = lastCurveID;  
//      PRINT_DEBUG_157(" last_curve_id is reset to %d \n", lastCurveID );      
    }

    if( abs(id1 - id2) == (last_curve_id - first_curve_id ) ){
        // firs and last curve of same loop
      return CUBIT_TRUE;
    }

      // part of same curve
    if( abs(id1 - id2) == 0 ){
      return CUBIT_TRUE;
    }
		
      // same loop but non-adj curves
    return CUBIT_FALSE;
  }
}

CubitBoolean TDOctreeRefFace::split_facet_type_03( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list ){

    // Do Nothing 
    // During CAT centroid will give the approximate branch point
  return CUBIT_TRUE;
}

CubitBoolean TDOctreeRefFace::split_facet_type_13( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list ){
  int i;
  CubitFacetEdge *boundary_edge = NULL, *ptr_edge;
  CubitPoint *boundary_edge_pnt0, *boundary_edge_pnt1, *other_boundary_point;
  double internal_angle;
  int num_of_segments;

  for( i = 0; i < 3; i++ ){
    ptr_edge = target_facet->edge(i);
    if( abs( ptr_edge->marked() ) == 1 ){
      boundary_edge = ptr_edge;
      break;
    }
  }

  boundary_edge_pnt0 = boundary_edge->point(0);
  boundary_edge_pnt1 = boundary_edge->point(1);

  for( i = 0; i < 3; i++ ){
    other_boundary_point = target_facet->point(i);
    if( other_boundary_point->id() != boundary_edge_pnt0->id() &&  other_boundary_point->id() != boundary_edge_pnt1->id() ){
      break;
    }
  }

  internal_angle = target_facet->angle( other_boundary_point );

  num_of_segments = (int)(floor( (internal_angle / (FACET_SPLITTING_INTERNAL_ANGLE )) + 0.5));

  if( num_of_segments > 1 )
    split_facet_locally_along_edge( target_facet, boundary_edge_pnt0, boundary_edge_pnt1, boundary_edge, num_of_segments, facet_list, facet_edge_list, point_list );

  return CUBIT_TRUE;
}

CubitBoolean TDOctreeRefFace::split_facet_type_23( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list ){

  int i;
  CubitFacetEdge *ptr_edge, *boundary_edge1, *boundary_edge2 = NULL,
      *boundary_split_edge;
  CubitPoint *other_boundary_point;
  double internal_angle;
  int num_of_segments;
  
  boundary_edge1 = NULL;
  for( i = 0; i < 3; i++ ){
    ptr_edge = target_facet->edge(i);
    if( abs( ptr_edge->marked() ) == 1 ){
      if( boundary_edge1 == NULL ){
        boundary_edge1 = ptr_edge;
      }
      else{
        boundary_edge2 = ptr_edge;
        break;
      }
    }
  }


  CubitFacetEdgeData *boundary_edge1_data = dynamic_cast<CubitFacetEdgeData*>( boundary_edge1 );
  CubitFacetEdgeData *boundary_edge2_data = dynamic_cast<CubitFacetEdgeData*>( boundary_edge2 );

  if( boundary_edge1_data->length() > boundary_edge2_data->length() ){
    boundary_split_edge = boundary_edge1;
  }
  else{
    boundary_split_edge = boundary_edge2;
  }

  for( i = 0; i < 3; i++ ){
    other_boundary_point = target_facet->point(i);
    if( other_boundary_point->id() != boundary_split_edge->point(0)->id() &&  other_boundary_point->id() != boundary_split_edge->point(1)->id() ){
      break;
    }
  }

  internal_angle = target_facet->angle( other_boundary_point );

  num_of_segments = (int)(floor( (internal_angle / (FACET_SPLITTING_INTERNAL_ANGLE )) + 0.5));

  if( num_of_segments > 1 )
    split_facet_locally_along_edge( target_facet, boundary_split_edge->point(0), boundary_split_edge->point(1), boundary_split_edge, num_of_segments, facet_list, facet_edge_list, point_list );

  return CUBIT_TRUE;
}

CubitBoolean TDOctreeRefFace::split_facet_type_33( CubitFacet *target_facet, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list ){

    // Do Nothing
    // 
    // WARNING: Degenerate Case
    // A triangualar surface patch is a facet.  
    // Skeleton should be calculated by taking circum/in center and radius.
  return CUBIT_TRUE;
}

CubitBoolean TDOctreeRefFace::split_facet_locally_along_edge( CubitFacet *target_facet, CubitPoint *edge1_pt, CubitPoint *edge2_pt, CubitFacetEdge *ptr_edge, int num_of_segments, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list){ 
  int i;
  CubitFacet *new_facet;
  CubitVector position, delta_length;
	
  delta_length = (edge2_pt->coordinates() - edge1_pt->coordinates()) / num_of_segments; 

  for( i = num_of_segments - 1; i >= 1; i-- ){
    position = edge1_pt->coordinates() + ( delta_length * i );

    CubitPointData* new_pt_data = new CubitPointData( position );
		
      // update TD and append the list
    CubitPoint *new_pt = CAST_TO(new_pt_data, CubitPoint);
    new_facet = split_facet_into_two_facets( target_facet, edge1_pt, edge2_pt, ptr_edge, new_pt, facet_list, facet_edge_list, point_list );
		
      // mark the newly created facets so that they are not checked for 
      // further splitting.  "this" facet is marked in Octree.cpp.
    new_facet->marked( -1 );
    edge2_pt = new_pt;
  }
  return CUBIT_TRUE; 
} 

CubitFacet *TDOctreeRefFace::split_facet_into_two_facets( CubitFacet *target_facet, CubitPoint* edge1_pt, CubitPoint* edge2_pt, CubitFacetEdge *edge, CubitPoint* new_pt, DLIList<CubitFacet*> &facet_list, DLIList<CubitFacetEdge*> &facet_edge_list, DLIList<CubitPoint*> &point_list){ 

    //add new_pt to point_list
  point_list.append( new_pt );
  new_pt->marked( abs( edge->marked() ));

    // split triangle
  CubitFacetData* facet_d = dynamic_cast<CubitFacetData*>(target_facet);
  assert(!!facet_d);

    // fix up existing facet
  int pt2_index = target_facet->point_index( edge2_pt );
  bool edge_reversed = ( edge1_pt == target_facet->point( (pt2_index+1) % 3 ) );
  int edge_index = (pt2_index + 1 + edge_reversed) % 3;
	
  edge2_pt->remove_facet( target_facet );
  facet_d->set_point( new_pt, pt2_index );
  new_pt->add_facet( target_facet );
  target_facet->update_plane();

    // make new facet and update facet list  
  CubitPoint* other_pt = target_facet->point( edge_index );
  CubitFacetData* new_facet;
  if( edge_reversed )
    new_facet = new CubitFacetData( other_pt, edge2_pt, new_pt );
  else
    new_facet = new CubitFacetData( other_pt, new_pt, edge2_pt );

    // facet constructor takes care of this
    // update the facets in points
    // this takes care of adj facets of new edges and 
    // conversely edges of two facets
    //edge2_pt->add_facet( new_facet ); 
    //new_pt->add_facet( new_facet );
    //other_pt->add_facet( new_facet );


    // Facets list in the points should be updated before creating new edge
    // the edge constructor generates adjacet facets.
    // split edge, if there is one
    //CubitFacetEdge* edge = edge1_pt->shared_edge( edge2_pt );
  CubitFacetEdgeData* new_edge = 0;
  if( edge ){
    CubitFacetEdgeData* edge_d = dynamic_cast<CubitFacetEdgeData*>(edge);
    assert(!!edge_d);
		
		
      // make sure new edge has same orientation as old edge
    new_edge = dynamic_cast<CubitFacetEdgeData*>(new_pt->shared_edge(edge2_pt));
    if( edge->point(0) == edge1_pt ){
      edge_d->set_point(new_pt, 1);
      if ( !new_edge )
        new_edge = new CubitFacetEdgeData( new_pt, edge2_pt );
      else if( new_edge->point(0) != new_pt )
        new_edge->flip();
    } 
    else {
      edge_d->set_point(new_pt, 0);
      if ( !new_edge )
        new_edge = new CubitFacetEdgeData( edge2_pt, new_pt );
      else if( new_edge->point(1) != new_pt )
        new_edge->flip();
    }
		
    new_edge->marked( abs( edge->marked() ) );
    facet_edge_list.append( new_edge );
  }
  else{
    PRINT_INFO("ERROR:Edge is doesn't exist. Splitting is not possible");
  }

  int sense;

    /*
      // constructor takes care of this
      if ( new_edge ) {
      assert(!new_facet->edge(0));
      new_facet->edge( new_edge, 0 );
      new_edge->add_facet( new_facet ); 
      sense = new_facet->point( 1 ) == new_edge->point(0) ? 1 : -1;
      new_facet->edge_use( sense, 0 );
      }
    */
    // facet_list appended
  facet_list.append( new_facet );

    // move other edge, if there is one
  int pt1_index = ( pt2_index + 2 - edge_reversed ) % 3;
  CubitFacetEdge* other_edge = target_facet->edge( pt1_index );
  int e_index;
  if( other_edge ){
    other_edge->remove_facet(target_facet);
      //target_facet->edge( 0, pt1_index );
    e_index = 1 + edge_reversed;
    assert(!new_facet->edge(e_index));
    new_facet->edge( other_edge, e_index );
    other_edge->add_facet( new_facet ); 
    sense = new_facet->point( ( e_index + 1 ) % 3 ) == other_edge->point(0) ? 1 : -1;
    new_facet->edge_use( sense, e_index );
  }

    // Add new mid edge in two facets

  CubitFacetEdgeData *new_mid_edge_data;
  CubitFacetEdge *new_mid_edge;
  if( edge_reversed )
    new_mid_edge_data = new CubitFacetEdgeData( other_pt, new_pt );
  else
    new_mid_edge_data = new CubitFacetEdgeData( new_pt, other_pt );

  new_mid_edge = CAST_TO( new_mid_edge_data, CubitFacetEdge );

  new_mid_edge->marked(2);

    // new_mid_edge is appended
  facet_edge_list.append( new_mid_edge );

    /* ////edge constructor takes care of this
       target_facet->edge( new_mid_edge, pt1_index ); 
	
       if( !edge_reversed ){  // new_mid_edge index 2
       new_facet->edge( new_mid_edge, ( e_index + 1 ) % 3 );
       sense = new_facet->point( ( e_index + 2 ) % 3 ) == new_mid_edge->point(0) ? 1 : -1;
       new_facet->edge_use( sense, (e_index + 1 ) % 3);
       }
       else{ // edge index 1
       new_facet->edge( new_mid_edge, ( e_index - 1 ) % 3 );
       sense = new_facet->point( ( e_index  ) % 3 ) == new_mid_edge->point(0) ? 1 : -1;
       new_facet->edge_use( sense, (e_index - 1 ) % 3);
       }
    */


#ifndef NDEBUG
    /* ---------------------- TESTING ----------------------- */
    // Testing the oerientation of faces w.r.t. the common edge
  
  int index0, index1;

  if( new_mid_edge->num_adj_facets() != 2 ){
    PRINT_INFO("ERROR:number of adjacent faces of new mid edge is not equal to 2 \n");
    assert(0);
  }

  CubitFacet *ptr_facet0 = new_mid_edge->adj_facet(0);
  CubitFacet *ptr_facet1 = new_mid_edge->adj_facet(1);

  for( index0 = 0; index0 < 3; index0++ ){
    if( ptr_facet0->edge(index0) == new_mid_edge )
      break;
  }

  for( index1 = 0; index1 < 3; index1++ ){
    if( ptr_facet1->edge(index1) == new_mid_edge )
      break;
  }

  if( index0 == 3 || index1 == 3 ){
    PRINT_INFO("ERROR: new edge doesn't exist at the adjacent facets \n");
    assert(0);
  }

  if( ptr_facet0->edge_use(index0) == ptr_facet1->edge_use(index1) ){
    PRINT_INFO("ERROR: The orientation of edges is adjacent facet is not proper \n");
    assert(0);
  }
  
#endif

    //PRINT_DEBUG_157("Passes\n");

  return CAST_TO( new_facet, CubitFacet );
}

// This only checks gfx facet edge valence right now (should be == 2)
CubitBoolean TDOctreeRefFace::check_valid_facets(CubitBoolean disable_if_bad)
{
  int i;
  CubitBoolean good_facets = CUBIT_TRUE;
  CubitFacetEdge *facet_edge;

  for (i=0; i < ptrCubitFacetEdgeList->size(); ++i)
  {
    facet_edge = ptrCubitFacetEdgeList->get_and_step();
    if (facet_edge->num_adj_facets() > 2)
    {
      good_facets = CUBIT_FALSE;
      break;
    }
  }

  if (disable_if_bad && !good_facets) {set_create_2dmat(CUBIT_FALSE);}
  return good_facets;
}

//EOF

