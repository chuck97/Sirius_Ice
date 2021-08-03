//-------------------------------------------------------------------------
// Filename      : TopologyEntity.cpp
//
// Purpose       : Implementation of the TopologyEntity class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/31/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitMessage.hpp"

#include "CastTo.hpp"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"

#include "ModelQueryEngine.hpp"
#include "GeometryQueryTool.hpp"

#include "Body.hpp"
#include "DLIList.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoVolume.hpp"
#include "CoFace.hpp"
#include "CoEdge.hpp"
#include "CoVertex.hpp"
#include "AppUtil.hpp"
#include "CubitEvent.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Casting query functions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
#define DECLARE_TOPO_ENT_QUERY_FUNC( TYPE, NAME, QUERYTYPE ) \
 CubitStatus TopologyEntity::NAME(DLIList<TYPE*>& list) \
 { \
   DLIList<TopologyEntity*> temp_list ; \
   CubitStatus result; \
   ModelQueryEngine *const mqe = ModelQueryEngine::instance(); \
   \
   result = mqe->query_model( *this, DagType::QUERYTYPE(), temp_list ); \
   if (result == CUBIT_FAILURE) \
   { \
      PRINT_ERROR("In TopologyEntity::" #NAME "\n"); \
      PRINT_ERROR("       Query failed for unknown reason.\n"); \
      return CUBIT_FAILURE; \
   } \
   \
   temp_list.reset(); \
   for (int i = temp_list.size(); i--; ) \
   { \
     TYPE* ptr = static_cast<TYPE*>(temp_list.get_and_step()); \
     assert(!!ptr); \
     list.append(ptr); \
   } \
   \
   return CUBIT_SUCCESS; \
}



//-------------------------------------------------------------------------
// Purpose       : Default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------

TopologyEntity::TopologyEntity() 
    : deactivatedStatus_(CUBIT_FALSE),
      encountered_(CUBIT_FALSE),
      bridgeMan()
{
  bridgeMan.set_entity(this);
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/10/96
//-------------------------------------------------------------------------

TopologyEntity::~TopologyEntity() 
{
  // Make sure that the TopologyEntity class does not have a pointer to this
  // TopologyEntity in its list.  If it does, then remove the instances
  // of the pointer
  DAG::instance()->remove(this);
}


GeometryQueryEngine*
TopologyEntity::get_geometry_query_engine() const
{
  TopologyBridge* bridge = 
    const_cast<TopologyEntity*>(this)->bridgeMan.topology_bridge();
  if (bridge)
    return bridge->get_geometry_query_engine();
  else
    return NULL;
}


DECLARE_TOPO_ENT_QUERY_FUNC(      Body,       bodies,       body_type )
DECLARE_TOPO_ENT_QUERY_FUNC(  CoVolume,   co_volumes,  co_volume_type )
DECLARE_TOPO_ENT_QUERY_FUNC( RefVolume,  ref_volumes, ref_volume_type )
DECLARE_TOPO_ENT_QUERY_FUNC(     Shell,       shells,      shell_type )
DECLARE_TOPO_ENT_QUERY_FUNC(    CoFace,     co_faces,    co_face_type )
DECLARE_TOPO_ENT_QUERY_FUNC(   RefFace,    ref_faces,   ref_face_type )
DECLARE_TOPO_ENT_QUERY_FUNC(      Loop,        loops,       loop_type )
DECLARE_TOPO_ENT_QUERY_FUNC(    CoEdge,     co_edges,    co_edge_type )
DECLARE_TOPO_ENT_QUERY_FUNC(   RefEdge,    ref_edges,   ref_edge_type )
DECLARE_TOPO_ENT_QUERY_FUNC(     Chain,       chains,      chain_type )
DECLARE_TOPO_ENT_QUERY_FUNC(  CoVertex,  co_vertices,  co_vertex_type )
DECLARE_TOPO_ENT_QUERY_FUNC( RefVertex, ref_vertices, ref_vertex_type )


RefVertex *TopologyEntity::ref_vertex()
{
  DLIList<RefVertex*> verts;
  ref_vertices(verts);
  if (verts.size() > 0) return verts.get();
  else return NULL;
}

RefEdge *TopologyEntity::ref_edge()
{
  DLIList<RefEdge*> edges;
  ref_edges(edges);
  if (edges.size() > 0) return edges.get();
  else return NULL;
}

RefFace *TopologyEntity::ref_face()
{
  DLIList<RefFace*> faces;
  ref_faces(faces);
  if (faces.size() > 0) return faces.get();
  else return NULL;
}

RefVolume *TopologyEntity::ref_volume()
{
  DLIList<RefVolume*> volumes;
  ref_volumes(volumes);
  if (volumes.size() > 0) return volumes.get();
  else return NULL;
}

Body *TopologyEntity::body()
{
  DLIList<Body*> these_bodies;
  this->bodies(these_bodies);
  if (these_bodies.size() > 0) return these_bodies.get();
  else return NULL;
}


CoEdge *TopologyEntity::co_edge() 
{
  DLIList<CoEdge*> these_co_edges;
  co_edges(these_co_edges);
  if (these_co_edges.size() > 0) return these_co_edges.get();
  else return NULL;
}

int TopologyEntity::num_loops() 
{
  DLIList<Loop*> these_loops;
  loops(these_loops);
  return these_loops.size();
}


int TopologyEntity::num_ref_volumes() 
{
  DLIList<RefVolume*> these_ref_volumes;
  ref_volumes(these_ref_volumes);
  return these_ref_volumes.size();
}

int TopologyEntity::num_ref_faces() 
{
  DLIList<RefFace*> these_ref_faces;
  ref_faces(these_ref_faces);
  return these_ref_faces.size();
}

int TopologyEntity::num_ref_edges() 
{
  DLIList<RefEdge*> these_ref_edges;
  ref_edges(these_ref_edges);
  return these_ref_edges.size();
}

int TopologyEntity::num_ref_vertices() 
{
  DLIList<RefVertex*> these_ref_vertices;
  ref_vertices(these_ref_vertices);
  return these_ref_vertices.size();
}


  //- return the number of connected entities of the given type

CubitBoolean TopologyEntity::is_directly_related( TopologyEntity* entity )
{
    // NULL?
  if ( !entity )
    return CUBIT_FALSE;
  
    // self?
  if (entity == this)
    return CUBIT_TRUE;  

    // same type but not self?
  if ( dag_type() == entity->dag_type() )  
    return CUBIT_FALSE;
  
    // Get the entities of the right type.
  DLIList<TopologyEntity*> kin_folk;
  CubitStatus result = ModelQueryEngine::instance()->
    query_model( *this, entity->dag_type(), kin_folk );
    //only fails if types are the same, caught above. 
  assert(result != CUBIT_FAILURE); 
  if (CUBIT_FAILURE == result) {
    PRINT_ERROR("ModelQueryEngine::query_model failed.\n");
    return CUBIT_FALSE;
  }

    // Search entities for passed in entity.
  TopologyEntity* model_entity = CAST_TO(entity, TopologyEntity) ;
  return kin_folk.is_in_list( model_entity );
}

CubitStatus TopologyEntity::set_topology_bridge(TopologyBridge* TB_ptr)
{
  assert (TB_ptr != NULL);
  CubitStatus rv = CUBIT_SUCCESS;
  
  BridgeManager* bridges = bridge_manager();
  bridges->remove_all_bridges();
  bridges->add_bridge(TB_ptr);

  return rv;
}

CubitStatus TopologyEntity::remove_from_DAG(CubitBoolean recurse_flag)
{
     // This counter will be used in the recursion to test whether
     // the call is from outside or from the function itself. When
     // the call comes from outside, the counter should always be
     // zero.

  CubitBoolean this_recurse = recurse_flag;
  if (recurse_flag == CUBIT_FALSE) recurse_flag = CUBIT_TRUE;

  DLIList<TopologyEntity*> childModEntList;

     // Check to see if there are no parents of this object.
   if ( get_parents() == 0 )
   {
      if (this_recurse == CUBIT_FALSE)
      {
         // Since we are not recursing, this is a top-level entity.
         // Notify the static observers that a top-level entity is being
         // destructed.  This must be done before children are disconnected.
       CubitObservable *top_level = CAST_TO(this, CubitObservable);
       AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOP_LEVEL_ENTITY_DESTRUCTED, dynamic_cast<RefEntity*>(top_level)));
      }

        // Go through all the children and remove their link to
        // the current object.

      TopologyEntity* tempModEntPtr = NULL ;
      TopologyEntity* childModEntPtr = NULL ;

      childModEntList.clean_out();
      disconnect_all_children(&childModEntList);

        // The following while conditional may not work...it depends on
        // what is_at_end does when you step over the end...CHECK THIS
      int i;
      for( i = 0 ; i < childModEntList.size() ; i++ )
      {
           // Get the next ModelEnti in the child list and make sure its
           // pointer to its parent is removed.
         tempModEntPtr = childModEntList.get_and_step();

           // Try remove() on the child ModEnt. If it comes back with
           // a success, then delete it.
         childModEntPtr = tempModEntPtr;

         if ( childModEntPtr->remove_from_DAG(recurse_flag) == CUBIT_SUCCESS )
         {
              // Now deactivate the child ModEnt
            childModEntPtr->deactivated(CUBIT_TRUE) ;

              // remove it from observables, just before we go to delete it
            CubitObservable *observable = CAST_TO(childModEntPtr, CubitObservable);
            if (observable)
            {
              AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_ENTITY_DESTRUCTED, dynamic_cast<RefEntity*>(observable)));
            }
         }
      }


        // If this is the top of the recursion, then clean out all the deactivated
        // entities.
      if (this_recurse == CUBIT_FALSE)
      {
         this->deactivated(CUBIT_TRUE) ;

         // remove it from observables, just before we go to delete it
         CubitObservable *observable = CAST_TO(childModEntPtr, CubitObservable);
         if (observable)
         {
           AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_ENTITY_DESTRUCTED, dynamic_cast<RefEntity*>(observable)));
         }
         GeometryQueryTool::instance()->cleanout_deactivated_geometry() ;
      }

      return CUBIT_SUCCESS ;
   }
   else
   {
      return CUBIT_FAILURE ;
   }
}

void TopologyEntity::disconnect_from_DAG()
{
    // disconnects this entity from any others
    // to which it is connected in the DAG; does not delete the DAGNode

  disconnect_all_children();
  disconnect_all_parents();
}

//-------------------------------------------------------------------------
// Purpose       : Get and set functions for the deactivated flag.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/02/96
//-------------------------------------------------------------------------

void TopologyEntity::deactivated(CubitBoolean flag)
{
   if (deactivatedStatus_ != flag)
   {
     deactivatedStatus_ = flag ;
     if (flag == CUBIT_TRUE)
     {
        DAG::instance()->add_deactivated_DAG_node(this) ;
     }
     else
     {
        DAG::instance()->remove_deactivated_DAG_node(this) ;
     }
   }
}

CubitBoolean TopologyEntity::deactivated() const
{
   return (CubitBoolean)deactivatedStatus_ ;
}




// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

