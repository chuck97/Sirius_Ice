
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "CubitDefines.h"
#include "AppUtil.hpp"
#include "ProgressTool.hpp"
#include "GeometryDefines.h"
#include "GeometryEntity.hpp"
#include "GeomMeasureTool.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "AnalyticGeometryTool.hpp"
#include "MergeTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeometryModifyEngine.hpp" // for creation of temporary construction entities
#include "DAG.hpp"
#include "TBOwnerSet.hpp"

#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "RefGroup.hpp"

#include "CoVertex.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "CoVolume.hpp"

#include "Chain.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Body.hpp"

#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "BodySM.hpp"
#include "ShellSM.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"

#include "CubitUtil.hpp"
#include "CubitAttrib.hpp"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"

#include "DLIList.hpp"
#include "GSaveOpen.hpp"

#include "CubitMessage.hpp"

#include "CastTo.hpp"
#include "CpuTimer.hpp"

#include "BridgeManager.hpp"
#include "TDUniqueId.hpp"
#include "CAMergePartner.hpp"
#include "CAActuateSet.hpp"
#include "CADeferredAttrib.hpp"
#include "CAUniqueId.hpp"

#include "SettingHandler.hpp"
#include "ModelQueryEngine.hpp"

#include "CubitTransformMatrix.hpp"
#include "CubitUndo.hpp"
#include "GMem.hpp"
#include "IntersectionTool.hpp"

#include "GfxPreview.hpp" //DJQ
#include "GfxDebug.hpp" //DJQ
#include "RefEntityName.hpp"

double GeometryQueryTool::geometryToleranceFactor = DEFAULT_GEOM_FACTOR;
GeometryQueryTool* GeometryQueryTool::instance_ = 0;
CubitBoolean GeometryQueryTool::useFacetBBox = CUBIT_FALSE;
CubitBoolean GeometryQueryTool::trackMergedAwayEnts = CUBIT_FALSE;
CubitBoolean GeometryQueryTool::importingSolidModel = CUBIT_FALSE;
CubitBoolean GeometryQueryTool::mergeGloballyOnImport = CUBIT_TRUE;
CubitBoolean GeometryQueryTool::clearUidMapBeforeImport = CUBIT_TRUE; 
DLIList<int> GeometryQueryTool::uidsOfImportingEnts;
int GeometryQueryTool::entitiesMergedAway = 0;

#ifndef CAT
  // Keep checking bounding box and use internal surf check for Sandia
  CubitBoolean GeometryQueryTool::bboxMergeTest = CUBIT_TRUE;
  int GeometryQueryTool::internalSurfaceMergeTest = 2; // 0=off, 1=all, 2=splines only
#else
//test
  // Cat prefers to avoid internal checks altogether as they cause problems with splines
  CubitBoolean GeometryQueryTool::bboxMergeTest = CUBIT_FALSE;
  int GeometryQueryTool::internalSurfaceMergeTest = 0; // 0=off, 1=all, 2=splines only
#endif // ndef CAT

double GeometryQueryTool::curveSliverCleanUpTolerance = geometryToleranceFactor*GEOMETRY_RESABS ;
double GeometryQueryTool::surfaceSliverCleanUpTolerance = -1.0;

// static CubitStatus import_actuate(DLIList<RefEntity*> &entity_list);

// ********** END STATIC DECLARATIONS      **********

#define PRINT(var) cout << #var << " = " << var << endl;

#define CUBIT_VERY_SMALL_NUMBER (CUBIT_RESABS * 1.0E-15)

// ********** BEGIN PUBLIC FUNCTIONS       **********
//-------------------------------------------------------------------------
// Purpose       : Controls access and creation of the sole instance of this
//                 class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------

GeometryQueryTool* GeometryQueryTool::instance(GeometryQueryEngine *GQEPtr)
{
     // Check to see if we have created an instance of the class
     // If proceed to create one.

   if (instance_ == 0)
   {
        // When creating the instance, we should always have a valid
        // GQEPtr. If not, complain.

     instance_ = new GeometryQueryTool (GQEPtr) ;

        // check to make sure there's a ref entity factory extant
     //RefEntityFactory *factory =
     RefEntityFactory::instance();
   }

     // If there is an existing instance of the class, check if there
     // was a request to set default solid modeling engine. If so, be nice
     // to the calling routine :) :) and kindly set the default solid
     // modeling engine.

   else if ( GQEPtr != NULL && !instance_->gqeList.move_to(GQEPtr)) {
       delete instance_->gqeList.remove();
       instance_->gqeList.insert(GQEPtr);
   }

     // Return the a pointer to the instance of the class.

   return instance_ ;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryQueryTool::~GeometryQueryTool ()
{

  //Kill the geometry query engine(s).
   int i;
   for (i = gqeList.size(); i > 0; i--)
   {
     delete gqeList.get_and_step();
   }
   gqeList.clean_out();

   instance_ = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Function to delete instance variable
//
// Special Notes :
//
// Creator       : Corey Ernst
//
// Creation Date : 12/31/07
//-------------------------------------------------------------------------
void GeometryQueryTool::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Creates temporary geometry files for save/restore
//
// Special Notes :
//
// Creator       : Corey Ernst
//
// Creation Date : 01/17/03
//-------------------------------------------------------------------------

CubitStatus GeometryQueryTool::save_temp_geom_files(DLIList<RefEntity*> &ref_entity_list,
                                             const char *base_filename,
                                             const CubitString &cubit_version,
                                             std::list<CubitString> &files_written,
                                             std::list<CubitString> &types_written)
{
  int i;
// clear all attributes
  if (ref_entity_list.size() == 0) {

      // All bodies are to be exported
    RefEntityFactory::instance()->ref_entity_list("Body",
                                                  ref_entity_list);

      // add free ref entities
    get_free_ref_entities(ref_entity_list);

  }

    // get all child entities
  DLIList<RefEntity*> child_list;
  RefEntity::get_all_child_ref_entities( ref_entity_list, child_list );

    // merge lists
  for(i = ref_entity_list.size(); i--; )
    ref_entity_list.get_and_step()->marked(1);
  for(i = child_list.size(); i--; )
    child_list.get_and_step()->marked(0);
  for(i = ref_entity_list.size(); i--; )
  {
    RefEntity* ent = ref_entity_list.get_and_step();
    if( ent->marked() )
    {
      ent->marked(0);
      child_list.append(ent);
    }
  }

    // now call auto update on this combined list; this will update both visible
    // and hidden entities; the combined list should be used here, but only the
    // export list should be exported (some of the hidden entities might be directly
    // related to other entities on the export list, and we want to avoid exporting
    // those entities twice)
  CubitAttribUser::auto_update_cubit_attrib(child_list);

    // Get list of TopologyBridges to save
  DLIList<TopologyBridge*> geometry_list(ref_entity_list.size()), ref_ent_bridges;
  ref_entity_list.reset();
  for ( i = ref_entity_list.size(); i--; )
  {
    RefEntity* ref_ent = ref_entity_list.get_and_step();
    TopologyEntity* topo_ent = dynamic_cast<TopologyEntity*>(ref_ent);
    if ( !topo_ent )
    {
      PRINT_ERROR("Attempt to save %s (%s %d) as geometry.\n",
        ref_ent->entity_name().c_str(), ref_ent->class_name(), ref_ent->id());
      continue;
    }

    ref_ent_bridges.clean_out();
    topo_ent->bridge_manager()->get_bridge_list(ref_ent_bridges);
    geometry_list += ref_ent_bridges;
  }

    // Save virtual.
  for (IGESet::reverse_iterator itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->export_geometry(geometry_list);

  //we need this list around so that we can remove the CSAs off
  // the entities under the virtual entities
  DLIList<TopologyBridge*> geometry_list2 = geometry_list;

  CubitStatus result = CUBIT_SUCCESS, temp_result;

  GeometryQueryEngine* gqe = NULL;
  //we're just going to mess with MBG stuff right now
  gqeList.reset();
  for( int k=0; k<gqeList.size(); k++)
  {
    gqe = gqeList.get_and_step();

    CubitString file_written;
    CubitString type_written;

    temp_result = gqe->save_temp_geom_file(geometry_list, base_filename,
                                           cubit_version, file_written,
                                           type_written );
    if( file_written.length() )
    {
      files_written.push_back( file_written );
      types_written.push_back( type_written );
    }

    if (CUBIT_SUCCESS == temp_result) result = temp_result;

  }


  //if there is still geometry left over....could not be handled
  if( geometry_list.size() )
    PRINT_ERROR("Not all geometry could be handled for save.\n");
  CubitAttribUser::clear_all_simple_attrib(child_list);

  //remove attributes off underlying entities of virtual geometry
  if( geometry_list2.size() )
    GeometryQueryTool::instance()->ige_remove_attributes( geometry_list2 );

  gqeList.reset();
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Export a list of solid model entities to a file.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 03/17/99
//-------------------------------------------------------------------------

CubitStatus GeometryQueryTool::export_solid_model(DLIList<RefEntity*>& ref_entity_list,
  char const* filename,
  Model_File_Type filetype,
  int &num_ents_exported,
  const CubitString &cubit_version,
  ModelExportOptions &export_options )
{
  if (0 == gqeList.size())
  {
    PRINT_WARNING("No active geometry engine.\n");
    return CUBIT_FAILURE;
  }

  int i;

  if (ref_entity_list.size() == 0) {

      // All bodies are to be exported
    RefEntityFactory::instance()->ref_entity_list("Body",
                                                  ref_entity_list);

      // add free ref entities
    get_free_ref_entities(ref_entity_list);

  }

    // Get TopologyBridges from RefEntities.
  DLIList<TopologyBridge*> bridge_list(ref_entity_list.size()),
                           parent_bridges, ref_ent_bridges;
  ref_entity_list.reset();
  for( i = ref_entity_list.size(); i--; )
  {
    ref_ent_bridges.clean_out();
    TopologyEntity* topo_ptr = dynamic_cast<TopologyEntity*>(ref_entity_list.get_and_step());
    if( topo_ptr )
      topo_ptr->bridge_manager()->get_bridge_list( ref_ent_bridges );
    bridge_list += ref_ent_bridges;
  }

    // Get all child RefEntities
  DLIList<RefEntity*> child_list;
  RefEntity::get_all_child_ref_entities( ref_entity_list, child_list );

    // Scan for free-but-merged entities in child list.
  child_list.reset();
  for (i = child_list.size(); i--; )
  {
    ref_ent_bridges.clean_out();
    TopologyEntity* topo_ptr = dynamic_cast<TopologyEntity*>(child_list.get_and_step());
    assert(!!topo_ptr);
    topo_ptr->bridge_manager()->get_bridge_list( ref_ent_bridges );
    ref_ent_bridges.reset();
    for (int j = ref_ent_bridges.size(); j--; )
    {
      TopologyBridge* bridge = ref_ent_bridges.get_and_step();
      parent_bridges.clean_out();
      bridge->get_parents(parent_bridges);
      if (parent_bridges.size() == 0)
        bridge_list.append_unique(bridge);
    }
  }

    // Merge lists so we have one big list of every
    // RefEntity to be saved.
  for(i = ref_entity_list.size(); i--; )
    ref_entity_list.get_and_step()->marked(1);

  for(i = child_list.size(); i--; )
    child_list.get_and_step()->marked(0);
  for(i = ref_entity_list.size(); i--; )
  {
    RefEntity* ent = ref_entity_list.get_and_step();
    if( ent->marked() )
    {
      ent->marked(0);
      child_list.append(ent);
    }
  }

  // Make a copy of the bridge list to be used below in
  // removing the virtual geometry attributes.  We can't
  // use the original bridge list because it gets emptied.
  DLIList<TopologyBridge*> copy_of_bridge_list = bridge_list;

  int num_input_entities = bridge_list.size();

    // now call auto update on this combined list; this will update both visible
    // and hidden entities; the combined list should be used here, but only the
    // export list should be exported (some of the hidden entities might be directly
    // related to other entities on the export list, and we want to avoid exporting
    // those entities twice)
  CubitAttribUser::auto_update_cubit_attrib(child_list);

    // Save virtual.
  IGESet::reverse_iterator itor;
  for (itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->export_geometry(bridge_list);

  CubitStatus result = CUBIT_SUCCESS, temp_result;

  if( !bridge_list.size() )
  {
    return CUBIT_SUCCESS;
  }
//   GeometryQueryEngine* gqe = NULL;
  //we're just going to mess with MBG stuff right now
  gqeList.reset();

  int num_ents_before = bridge_list.size();
  temp_result = gqeList.get()->export_solid_model(bridge_list, filename, filetype,
                                        cubit_version, export_options );
  if (temp_result == CUBIT_SUCCESS )
    result = temp_result;

  //if all geometry wasn't exported, warn user and print out
  //what wasn't
  if( bridge_list.size() != 0 )
  {
    if( bridge_list.size() == num_ents_before )
      PRINT_ERROR("No geometry exported.  Must set geometry engine to another type.\n");
    else
      PRINT_WARNING("Not all geometry could be handled for save.\n");

    PRINT_INFO("Set the geometry engine appropriately to export the following geometry:\n");
    int k;
    for(k=bridge_list.size(); k--; )
    {
      TopologyEntity *te = bridge_list.get()->topology_entity();
      if( te )
      {
        RefEntity *ref_ent = CAST_TO( te, RefEntity );
        if( ref_ent )
        {
          GeometryQueryEngine *gqe = bridge_list.get()->get_geometry_query_engine();
          PRINT_INFO("%s is of Geometry Engine type %s\n",
          ref_ent->entity_name().c_str(), gqe->modeler_type() );
        }
      }
      bridge_list.step();
    }
  }

    // clean off attributes off of underyling virtual geometry
  for (itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->remove_attributes(copy_of_bridge_list);

  CubitAttribUser::clear_all_simple_attrib(child_list);

  num_ents_exported = num_input_entities - bridge_list.size();
  gqeList.reset();
  return result;
}

// export entities to buffer
CubitStatus GeometryQueryTool::export_solid_model(DLIList<RefEntity*>& ref_entity_list,
                                                  char*& p_buffer,
                                                  int& n_buffer_size,
                                                  bool b_export_buffer)
{
  if (0 == gqeList.size())
  {
    PRINT_WARNING("No active geometry engine.\n");
    return CUBIT_FAILURE;
  }

  int i;

  if (ref_entity_list.size() == 0) {
    // All bodies are to be exported
    RefEntityFactory::instance()->ref_entity_list("Body",
                                                  ref_entity_list);

    // add free ref entities
    get_free_ref_entities(ref_entity_list);
  }

  // Get TopologyBridges from RefEntities.
  DLIList<TopologyBridge*> bridge_list(ref_entity_list.size()),
    parent_bridges, ref_ent_bridges;
  ref_entity_list.reset();
  for( i = ref_entity_list.size(); i--; ) {
    ref_ent_bridges.clean_out();
    TopologyEntity* topo_ptr = dynamic_cast<TopologyEntity*>(ref_entity_list.get_and_step());
    if( topo_ptr )
      topo_ptr->bridge_manager()->get_bridge_list( ref_ent_bridges );
    bridge_list += ref_ent_bridges;
  }

  // Get all child RefEntities
  DLIList<RefEntity*> child_list;
  RefEntity::get_all_child_ref_entities( ref_entity_list, child_list );

  // Scan for free-but-merged entities in child list.
  child_list.reset();
  for (i = child_list.size(); i--; )
  {
    ref_ent_bridges.clean_out();
    TopologyEntity* topo_ptr = dynamic_cast<TopologyEntity*>(child_list.get_and_step());
    assert(!!topo_ptr);
    topo_ptr->bridge_manager()->get_bridge_list( ref_ent_bridges );
    ref_ent_bridges.reset();
    for (int j = ref_ent_bridges.size(); j--; )
    {
      TopologyBridge* bridge = ref_ent_bridges.get_and_step();
      parent_bridges.clean_out();
      bridge->get_parents(parent_bridges);
      if (parent_bridges.size() == 0)
        bridge_list.append_unique(bridge);
    }
  }

  // Merge lists so we have one big list of every
  // RefEntity to be saved.
  for(i = ref_entity_list.size(); i--; )
    ref_entity_list.get_and_step()->marked(1);

  for(i = child_list.size(); i--; )
    child_list.get_and_step()->marked(0);
  for(i = ref_entity_list.size(); i--; )
  {
    RefEntity* ent = ref_entity_list.get_and_step();
    if( ent->marked() )
    {
      ent->marked(0);
      child_list.append(ent);
    }
  }

  // Make a copy of the bridge list to be used below in
  // removing the virtual geometry attributes.  We can't
  // use the original bridge list because it gets emptied.
  DLIList<TopologyBridge*> copy_of_bridge_list = bridge_list;

  /*int num_input_entities =*/ bridge_list.size();

    // now call auto update on this combined list; this will update both visible
    // and hidden entities; the combined list should be used here, but only the
    // export list should be exported (some of the hidden entities might be directly
    // related to other entities on the export list, and we want to avoid exporting
    // those entities twice)
  CubitAttribUser::auto_update_cubit_attrib(child_list);

    // Save virtual.
  IGESet::reverse_iterator itor;
  for (itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->export_geometry(bridge_list);

  CubitStatus result = CUBIT_SUCCESS, temp_result;

  if( !bridge_list.size() )
  {
    return CUBIT_SUCCESS;
  }
  //   GeometryQueryEngine* gqe = NULL;
  //we're just going to mess with MBG stuff right now
  gqeList.reset();

  int num_ents_before = bridge_list.size();
  temp_result = gqeList.get()->export_solid_model(bridge_list, p_buffer,
                                                  n_buffer_size, b_export_buffer);
  if (temp_result == CUBIT_SUCCESS )
    result = temp_result;

  //if all geometry wasn't exported, warn user and print out
  //what wasn't
  if( bridge_list.size() != 0 )
  {
    if( bridge_list.size() == num_ents_before )
      PRINT_ERROR("No geometry exported.  Must set geometry engine to another type.\n");
    else
      PRINT_WARNING("Not all geometry could be handled for save.\n");

    PRINT_INFO("Set the geometry engine appropriately to export the following geometry:\n");
    int k;
    for(k=bridge_list.size(); k--; )
    {
      TopologyEntity *te = bridge_list.get()->topology_entity();
      if( te )
      {
        RefEntity *ref_ent = CAST_TO( te, RefEntity );
        if( ref_ent )
        {
          GeometryQueryEngine *gqe = bridge_list.get()->get_geometry_query_engine();
          PRINT_INFO("%s is of Geometry Engine type %s\n",
          ref_ent->entity_name().c_str(), gqe->modeler_type() );
        }
      }
      bridge_list.step();
    }
  }

    // clean off attributes off of underyling virtual geometry
  for (itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->remove_attributes(copy_of_bridge_list);

  CubitAttribUser::clear_all_simple_attrib(child_list);

  //num_ents_exported = num_input_entities - bridge_list.size();
  gqeList.reset();
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Fire a ray at a list of entities and return the
//                 parameters along the ray (distance from origin of ray)
//                 where it hits the entities; optionally find the
//                 corresponding entities it hit.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 3/29/2007
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::fire_ray( CubitVector &origin,
                                         CubitVector &direction,
                                         DLIList<RefEntity*> &at_entity_list,
                                         DLIList<double> &ray_params,
                                         int max_hits,
                                         double ray_radius,
                                         DLIList<RefEntity*> *hit_entity_list_ptr )
{
  int i;

  // Do this in a way to account for the case if they happen to be from
  // different geometry engines (could easily occur with virtual geometry)

  DLIList<TopologyEntity*> te_fire_at_list;
  DLIList<TopologyEntity*> te_hit_entity_list;
  DLIList<TopologyEntity*> *te_hit_entity_list_ptr = 0;
  if( hit_entity_list_ptr )
    te_hit_entity_list_ptr = &te_hit_entity_list;

  GeometryQueryEngine *gqe = 0;
  GeometryQueryEngine *gqe_last = 0;
  RefEntity *ref_entity_ptr;
  TopologyEntity *topo_ptr;

  int loc_max_hits = max_hits;
  CubitBoolean hits_limited = CUBIT_FALSE;
  if( max_hits > 0 )
    hits_limited = CUBIT_TRUE;

  // Note we care about order
  at_entity_list.reset();
  for( i=at_entity_list.size(); i--; )
  {
    ref_entity_ptr = at_entity_list.get_and_step();

    topo_ptr = CAST_TO( ref_entity_ptr, TopologyEntity );
    if( !topo_ptr )
    {
      PRINT_ERROR( "Couldnt get topo_ptr\n" );
      continue;
    }

    gqe = topo_ptr->get_geometry_query_engine();
    if( !gqe )
    {
      PRINT_ERROR( "Unable to find geometry engine associated with an entity!\n" );
      return CUBIT_FAILURE;
    }

    if( !gqe_last )
    {
      gqe_last = gqe;
    }
    if( gqe != gqe_last )
    {
      if( hits_limited == CUBIT_TRUE )
      {
        loc_max_hits = max_hits - ray_params.size();
        if( loc_max_hits <= 0 )
          break;
      }

      if( fire_ray( origin, direction, te_fire_at_list, ray_params,
        loc_max_hits, ray_radius, te_hit_entity_list_ptr ) == CUBIT_FAILURE )
        return CUBIT_FAILURE;

      // Reset
      gqe_last = gqe;
      te_fire_at_list.clean_out();
    }

    te_fire_at_list.append( topo_ptr );
  }

  // Do the last ray fire, if necessary
  if( hits_limited == CUBIT_TRUE )
    loc_max_hits = max_hits - ray_params.size();
  if( hits_limited==CUBIT_FALSE || loc_max_hits>0 )
    if( fire_ray( origin, direction, te_fire_at_list, ray_params,
      loc_max_hits, ray_radius, te_hit_entity_list_ptr ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

  if( hit_entity_list_ptr )
  {
    // Find RefEntities from TopologyEntities
    te_hit_entity_list.reset();
    for( i=te_hit_entity_list.size(); i--; )
    {
      topo_ptr = te_hit_entity_list.get_and_step();
      if( !topo_ptr )
      {
        hit_entity_list_ptr->append( 0 );
        continue;
      }

      ref_entity_ptr = CAST_TO( topo_ptr, RefEntity );
      hit_entity_list_ptr->append( ref_entity_ptr );
    }
  }

  // Now, make sure we don't have more hits than asked for
  if( hits_limited == CUBIT_TRUE )
  {
    if( ray_params.size() <= max_hits )
      return CUBIT_SUCCESS;

    for( i=ray_params.size()-max_hits; i--; )
    {
      ray_params.last();
      ray_params.remove();
      if( hit_entity_list_ptr )
      {
        hit_entity_list_ptr->last();
        hit_entity_list_ptr->remove();
      }
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Fire a ray at a list of entities and return the
//                 parameters along the ray (distance from origin of ray)
//                 where it hits the entities; optionally find the
//                 corresponding entities it hit.
//
// Special Notes : All entities must be from the same geometry engine.
//
// Creator       : Steve Storm
//
// Creation Date : 5/19/2007
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::fire_ray( CubitVector &origin,
                                         CubitVector &direction,
                                         DLIList<TopologyEntity*> &at_entity_list,
                                         DLIList<double> &ray_params,
                                         int max_hits,
                                         double ray_radius,
                                         DLIList<TopologyEntity*> *hit_entity_list_ptr )
{
  if( !at_entity_list.size() )
    return CUBIT_SUCCESS;


  TopologyEntity *topo_ptr = at_entity_list.get();

  GeometryQueryEngine *gqe = topo_ptr->get_geometry_query_engine();
  if( !gqe )
  {
    PRINT_ERROR( "Unable to find geometry engine associated with an entity!\n" );
    return CUBIT_FAILURE;
  }


  DLIList<TopologyBridge*> tb_list;
  int i;
  at_entity_list.reset();
  for( i=at_entity_list.size(); i--; )
  {
    topo_ptr = at_entity_list.get_and_step();
    DLIList<TopologyBridge*> bridge_list;
    topo_ptr->bridge_manager()->get_bridge_list( bridge_list );
    bridge_list.reset();
    tb_list.append( bridge_list.get() );
  }

  // Setup temporary variables to pass
  DLIList<double> tmp_ray_params;

  DLIList<TopologyBridge*> tb_hit_list;
  DLIList<TopologyBridge*> *tb_hit_list_ptr = NULL;
  if( hit_entity_list_ptr )
    tb_hit_list_ptr = &tb_hit_list;

  // Do the ray fire.  Note we will sort the hits by distance and append to the output lists.
  if( gqe->fire_ray( origin, direction, tb_list, tmp_ray_params,
    max_hits, ray_radius, tb_hit_list_ptr ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  tmp_ray_params.reset();
  if( tb_hit_list_ptr) tb_hit_list_ptr->reset();

  // Append to output lists
  ray_params += tmp_ray_params;

  // Get the TE's from the TopologyBridges
  if( hit_entity_list_ptr )
  {
	  // First need to make sure that the entities hit are visible...entities hit
	  // may be hidden under virtual entities....we need to get the visible
	  // entities.
	  DLIList<TopologyBridge*> vis_tb_hit_list;
	  DLIList<CubitVector*> cv_list;

	  CubitVector *loc_ptr;
	  double param;
	  tmp_ray_params.reset();
	  for( i=tmp_ray_params.size(); i--; )
	  {
		  param = tmp_ray_params.get_and_step();

		  loc_ptr = new CubitVector;
		  origin.next_point( direction, param, *loc_ptr );
		  cv_list.append( loc_ptr );
	  }


	  TopologyBridge *bridge_ptr;
	  for( i=0; i<tb_hit_list_ptr->size(); i++ )
	  {
		  bridge_ptr = tb_hit_list_ptr->get_and_step();

		  TopologyBridge *visible_tb = NULL;
		  TBOwner* o2 = bridge_ptr->owner();

		  bool broke_early = false;
		  BridgeManager* bridge_manager2;
		  while (!(bridge_manager2 = dynamic_cast<BridgeManager*>(o2)))
		  {
			  if (TopologyBridge* bridge2 = dynamic_cast<TopologyBridge*>(o2))
			  {
				  GeometryQueryEngine* gqe2 = bridge2->get_geometry_query_engine();

				  //Let the VQE handle the work
				  visible_tb = gqe2->get_visible_entity_at_point(bridge_ptr, cv_list[i]);
				  if (visible_tb)
					  o2 = visible_tb->owner();
				  else
					  o2 = bridge2->owner();
			  }

			  else if(TBOwnerSet* set = dynamic_cast<TBOwnerSet*>(o2))
			  {
				  DLIList<TopologyBridge*> list2;
				  set->get_owners(list2);
				  list2.reset();

				  // This had better be the Virtual QE.
				  GeometryQueryEngine* gqe2 = list2.get()->get_geometry_query_engine();

				  //Let the VQE handle the work
				  visible_tb = gqe2->get_visible_entity_at_point(bridge_ptr, cv_list[i]);
				  if (visible_tb)
					  o2 = visible_tb->owner();
				  else
				  {
					  broke_early = true;
					  break;
				  }
			  }
			  else
			  {
				  broke_early = true;
				  break;
			  }
		  }

		  if (!broke_early)
			  visible_tb = bridge_manager2->topology_bridge();


		  if( visible_tb )
		  {
			  topo_ptr = visible_tb->topology_entity();
			  hit_entity_list_ptr->append( topo_ptr );
		  }
		  else
			  hit_entity_list_ptr->append( 0 );
	  }


	  // Free memory
	  while( cv_list.size() ) delete cv_list.pop();
  }

  // Sort ray_params (low to high) and sort hit_entity_list_ptr to match
  //  This will ensure entities are in order of who got hit first
  //  Do the sort by adding to a map (should auto-sort)
  std::map<double, TopologyEntity*> temp_map;
  for (i=0; i<ray_params.size(); i++)
	  temp_map.insert(std::map<double, 
    TopologyEntity*>::value_type( ray_params.get_and_step(), 
    hit_entity_list_ptr ? hit_entity_list_ptr->get_and_step() : 0 ) );

  // The map should be sorted, so iterate through it and add to the official lists
  ray_params.clean_out();
  if( hit_entity_list_ptr) hit_entity_list_ptr->clean_out();

  std::map<double, TopologyEntity*>::iterator iter;
  for (iter=temp_map.begin(); iter != temp_map.end(); iter++)
  {
	  ray_params.append(iter->first);
	  if( hit_entity_list_ptr) hit_entity_list_ptr->append(iter->second);
  }

	return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Reads in geometry and creates the necessary Reference
//                 entities associated with the input geometry. Has ability
//                 to read selectively read in either bodies, free surfaces,
//                 free curves or free vertices from the file.
//
//                 Valid file types are:
//                       "ACIS_SAT"    --  ACIS ASCII (SAT) file format
//                       "ACIS_SAB"    --  ACIS BINARY (SAB) file format
//                       "IGES"        --  IGES file
//                       "STEP"        --  STEP file
//                       "KCM"         --  KCM file
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 03/17/99
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::import_solid_model( const char* file_name,
  Model_File_Type file_type,
  ModelImportOptions &import_options,
  DLIList<RefEntity*> *imported_entities )
{
  if (0 == gqeList.size())
  {
    PRINT_WARNING("No active geometry engine.\n");
    return CUBIT_FAILURE;
  }

  importingSolidModel = CUBIT_TRUE;
  uidsOfImportingEnts.clean_out();
  mergeGloballyOnImport = import_options.merge_globally;
 
 if( clearUidMapBeforeImport )
    CAUniqueId::clear_out_old_to_new_map();

  if( GSaveOpen::performingUndo )
    mergeGloballyOnImport = CUBIT_TRUE;

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

    // reset the attribImporteds flags to facilitate attribute reporting
//  CubitAttrib::clear_attrib_importeds();

    // Use the default MQE to import a list of ToplogyBridges from the file.
  gqeList.reset();
  DLIList<TopologyBridge*> bridge_list;

  import_options.print_results = CUBIT_TRUE;

  // We will turn off the printing of name change warnings during import and just
  // keep track of how many happened. 
  bool prev_print_warning_flag = RefEntityName::instance()->print_name_change_warnings();
  RefEntityName::instance()->print_name_change_warnings(false);
  RefEntityName::instance()->num_name_change_warnings(0);
  
  CubitStatus status = gqeList.get()->import_solid_model( file_name,
      file_type, bridge_list, import_options );
      
  if( bridge_list.size() == 0 )
  {
    importingSolidModel = CUBIT_FALSE;
    mergeGloballyOnImport = CUBIT_FALSE; 
    return status;
  }

  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->import_geometry(bridge_list);

  bridge_list.reset();
  DLIList<RefEntity*> tmp_ent_list;
  status = construct_refentities(bridge_list, &tmp_ent_list);

  // Print out a summary warning message if names were changed.
  RefEntityName::instance()->print_name_change_warnings(prev_print_warning_flag);
  if(RefEntityName::instance()->num_name_change_warnings() > 0)
  {
    PRINT_WARNING("%d invalid names were found during import.\n\n",
            RefEntityName::instance()->num_name_change_warnings());
  }
    
  if( imported_entities )
    (*imported_entities) += tmp_ent_list;

  if( CubitUndo::get_undo_enabled() )
  {
    if( tmp_ent_list.size() && status == CUBIT_SUCCESS )
      CubitUndo::note_result_entities( tmp_ent_list );
    else
      CubitUndo::remove_last_undo();
  }

  //clear out all attributes that were set to actuate
  DLIList<RefEntity*> all_ents;
  RefEntity::get_all_child_ref_entities( tmp_ent_list, all_ents );
  all_ents += tmp_ent_list;
  CubitAttribUser::clear_all_simple_attrib_set_to_actuate( all_ents );

    // report attribs imported
//  CubitAttrib::report_attrib_importeds();

  importingSolidModel = CUBIT_FALSE;
  mergeGloballyOnImport = CUBIT_FALSE; 
    // now return
  return status;
}

// import entities to buffer
CubitStatus GeometryQueryTool::import_solid_model(DLIList<RefEntity*> *imported_entities,
                                                  const char* pBuffer,
                                                  const int n_buffer_size)
{
  if (0 == gqeList.size())
  {
    PRINT_WARNING("No active geometry engine.\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

    // reset the attribImporteds flags to facilitate attribute reporting
//  CubitAttrib::clear_attrib_importeds();

    // Use the default MQE to import a list of ToplogyBridges from the file.
  gqeList.reset();
  DLIList<TopologyBridge*> bridge_list;

  CubitStatus status = CUBIT_SUCCESS;
  for(int i = 0; i < gqeList.size(); i++)
  {
    status = gqeList.get_and_step()->import_solid_model( bridge_list, pBuffer, n_buffer_size );

    if( bridge_list.size() > 0 )
      break;
  }
  if(bridge_list.size() == 0)
    return status;

  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->import_geometry(bridge_list);

  bridge_list.reset();
  DLIList<RefEntity*> tmp_ent_list;
  status = construct_refentities(bridge_list, &tmp_ent_list);

  if( imported_entities )
    (*imported_entities) += tmp_ent_list;

  if( CubitUndo::get_undo_enabled() )
  {
    if( tmp_ent_list.size() && status == CUBIT_SUCCESS )
      CubitUndo::note_result_entities( tmp_ent_list );
    else
      CubitUndo::remove_last_undo();
  }

  //clear out all attributes that were set to actuate
//  DLIList<RefEntity*> all_ents;
//  RefEntity::get_all_child_ref_entities( tmp_ent_list, all_ents );
//  all_ents += tmp_ent_list;
//  CubitAttribUser::clear_all_simple_attrib_set_to_actuate( all_ents );

    // report attribs imported
//  CubitAttrib::report_attrib_importeds();


    // now return
  return status;
}

CubitStatus GeometryQueryTool::construct_refentities(DLIList<TopologyBridge*> &bridge_list,
                                                     DLIList<RefEntity*> *imported_entities)
{

    // Construct VGI Topology from the TopologyBridges.
  DLIList<Body*> body_list;
  DLIList<RefFace*> face_list;
  DLIList<RefEdge*> edge_list;
  DLIList<RefVertex*> vtx_list;
    // We are going to pop entities from bridge_list,
    // so reverse the list so that the entities get
    // created in the old order and the IDs remain
    // the same.
  bridge_list.reverse();
  bridge_list.reset();
  CubitStatus status = CUBIT_SUCCESS;

  bool do_progress = bridge_list.size() > 2 ? true : false;

  if(do_progress)
  {
    char message[128];
    sprintf(message, "Building %d CUBIT Entities", bridge_list.size() );
    AppUtil::instance()->progress_tool()->start(0, bridge_list.size(), "Progress",
        message, TRUE, TRUE);
  }


  while( bridge_list.size() && !AppUtil::instance()->interrupt() )
  {
    BodySM* bodysm_ptr;
    Surface* surface_ptr;
    Curve* curve_ptr;
    TBPoint* point_ptr;
    TopologyBridge* bridge_ptr = bridge_list.pop();
    if( (bodysm_ptr = dynamic_cast<BodySM*>(bridge_ptr) ) != NULL )
    {
      Body* body_ptr = make_Body( bodysm_ptr );
      if( body_ptr )
      {
        body_list.append( body_ptr );
      }
      else
      {
        PRINT_ERROR("\nError constructing VGI Body.\n");
        status = CUBIT_FAILURE;
      }
    }
    else if( (surface_ptr = dynamic_cast<Surface*>(bridge_ptr) ) != NULL )
    {
      bool is_free_face = true;
      RefFace* face_ptr = make_free_RefFace( surface_ptr, is_free_face );
      if( face_ptr )
      {
        face_list.append( face_ptr );
      }
      else
      {
        PRINT_ERROR("\nError constructing free VGI Surface.\n");
        status = CUBIT_FAILURE;
      }
    }
    else if( (curve_ptr = dynamic_cast<Curve*>(bridge_ptr) ) != NULL )
    {
      RefEdge* edge_ptr = make_free_RefEdge( curve_ptr );
      if( edge_ptr )
      {
        edge_list.append( edge_ptr );
      }
      else
      {
        PRINT_ERROR("\nError constructing free VGI Curve.\n");
        status = CUBIT_FAILURE;
      }
    }
    else if( (point_ptr = dynamic_cast<TBPoint*>(bridge_ptr) ) != NULL )
    {
      RefVertex* vtx_ptr = make_free_RefVertex( point_ptr );
      if( vtx_ptr )
      {
        vtx_list.append( vtx_ptr );
      }
      else
      {
        PRINT_ERROR("\nError constructing free VGI Vertex.\n");
        status = CUBIT_FAILURE;
      }
    }
    else
    {
      PRINT_ERROR("Unknown entity type imported by solid modeler.  Ignored.\n");
      status = CUBIT_FAILURE;
    }
    if(do_progress)
        AppUtil::instance()->progress_tool()->step();
  }

  if(do_progress)
    AppUtil::instance()->progress_tool()->end();

  PRINT_INFO("\n");

    // If the above loop was interrupted, bridge_list will not be
    // empty.  Free any bridges we did not make topology for.
  if( AppUtil::instance()->interrupt() )
  {
    PRINT_WARNING("Aborted construction of %d entities.\n", bridge_list.size());
  }
  while( bridge_list.size() )
  {
    TopologyBridge* bridge_ptr = bridge_list.pop();
    GeometryQueryEngine* gqe = bridge_ptr->get_geometry_query_engine();
    BodySM* body_ptr;
    Surface* surface_ptr;
    Curve* curve_ptr;
    TBPoint* point_ptr;
    if( (body_ptr = dynamic_cast<BodySM*>(bridge_ptr) ) != NULL )
      gqe->delete_solid_model_entities( body_ptr );
    else if( (surface_ptr = dynamic_cast<Surface*>(bridge_ptr) ) != NULL )
      gqe->delete_solid_model_entities( surface_ptr );
    else if( (curve_ptr = dynamic_cast<Curve*>(bridge_ptr) ) != NULL )
      gqe->delete_solid_model_entities( curve_ptr );
    else if( (point_ptr = dynamic_cast<TBPoint*>(bridge_ptr) ) != NULL )
      gqe->delete_solid_model_entities( point_ptr );
  }

/*
    // If any geometry was pre-merged during the
    // RefEntity construction, adjust the sense of
    // the entities to match what it was originally.
  DLIList<BasicTopologyEntity*> bte_list;
  DLIList<RefFace*> child_face_list;
  DLIList<RefEdge*> child_edge_list;
  for( int b = body_list.size(); b--; )
  {
    body_list.get()->ref_faces( child_face_list );
    CAST_LIST_TO_PARENT( child_face_list, bte_list );
    adjust_merge_sense( bte_list );
    body_list.get()->ref_edges( child_edge_list );
    CAST_LIST_TO_PARENT( child_edge_list, bte_list );
    adjust_merge_sense( bte_list );
    body_list.step();
  }
  CAST_LIST_TO_PARENT( face_list, bte_list );
  adjust_merge_sense( bte_list );
  for( int f = face_list.size(); f--; )
  {
    face_list.get_and_step()->ref_edges( child_edge_list );
    CAST_LIST_TO_PARENT( child_edge_list, bte_list );
    adjust_merge_sense( bte_list );
  }
  CAST_LIST_TO_PARENT( edge_list, bte_list );
  adjust_merge_sense( bte_list );
*/

    // Actuate Attributes
  DLIList<RefEntity*> ref_entity_list, temp_ref_list;

  CAST_LIST_TO_PARENT( body_list, temp_ref_list );
  ref_entity_list += temp_ref_list;
  CAST_LIST_TO_PARENT( face_list, temp_ref_list );
  ref_entity_list += temp_ref_list;
  CAST_LIST_TO_PARENT( edge_list, temp_ref_list );
  ref_entity_list += temp_ref_list;
  CAST_LIST_TO_PARENT(  vtx_list, temp_ref_list );
  ref_entity_list += temp_ref_list;

  import_actuate( ref_entity_list );

    // Pass back imported entities
  if( imported_entities )
  {
    *imported_entities = ref_entity_list;
  }

  char pre[100];
  if( body_list.size() )
  {
    DLIList<RefVolume*> temp_vols;
    for( int nv = body_list.size(); nv > 0; nv-- )
    {
      DLIList<RefVolume*> t2;
      body_list.get_and_step()->ref_volumes( t2 );
      temp_vols += t2;
    }

    DLIList<CubitEntity*> temp_list;
    if( DEBUG_FLAG( 153 ) )
    {
      CAST_LIST_TO_PARENT(body_list, temp_list);
      if( temp_list.size() == 1 )
         CubitUtil::list_entity_ids( "Constructed 1 Body: ", temp_list);
      else
      {
        sprintf( pre, "Constructed %d Bodies: ", temp_list.size() );
        CubitUtil::list_entity_ids( pre, temp_list );
      }
    }

    temp_list.clean_out();
    CAST_LIST_TO_PARENT( temp_vols, temp_list );
    if( temp_list.size() == 1 )
       CubitUtil::list_entity_ids( "Constructed 1 Volume: ", temp_list);
    else
    {
      sprintf( pre, "Constructed %d Volumes: ", temp_list.size() );
      CubitUtil::list_entity_ids( pre, temp_list );
    }
  }
  if( face_list.size() )
  {
    DLIList<CubitEntity*> temp_list;
    CAST_LIST_TO_PARENT(face_list, temp_list);
    if( temp_list.size() == 1 )
      CubitUtil::list_entity_ids( "Constructed 1 Free Surface: ", temp_list );
    else
    {
      sprintf( pre, "Constructed %d Free Surfaces: ", temp_list.size() );
      CubitUtil::list_entity_ids( pre, temp_list );
    }
  }
  if( edge_list.size() )
  {
    DLIList<CubitEntity*> temp_list;
    CAST_LIST_TO_PARENT(edge_list, temp_list);
    if( temp_list.size() == 1 )
      CubitUtil::list_entity_ids( "Constructed 1 Free Curve: ", temp_list );
    else
    {
      sprintf( pre, "Constructed %d Free Curves: ", temp_list.size() );
      CubitUtil::list_entity_ids( pre, temp_list );
    }
  }
  if( vtx_list.size() )
  {
    DLIList<CubitEntity*> temp_list;
    CAST_LIST_TO_PARENT(vtx_list, temp_list);
    if( temp_list.size() == 1 )
      CubitUtil::list_entity_ids( "Constructed 1 Free Vertex: ", temp_list );
    else
    {
      sprintf( pre, "Constructed %d Free Vertices: ", temp_list.size() );
      CubitUtil::list_entity_ids( pre, temp_list );
    }
  }

  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Build/update Body topology
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/21/02
//-------------------------------------------------------------------------
Body* GeometryQueryTool::make_Body(BodySM *bodysm_ptr) const
{
  int i;
  CpuTimer timer;
  bool body_created = false;
  bool body_modified = false;
  bool vol_modified = false, tmp_vol_modified;

  Body* body = CAST_TO( bodysm_ptr->topology_entity(), Body );
  if( ! body )
  {
    body = RefEntityFactory::instance()->construct_Body(bodysm_ptr);
    body_created = true;
  }
  if( ! body )
  {
    PRINT_ERROR("Body creation failed in GeometryQueryTool::make_Body\n");
    assert( !!body );
    return 0;
  }

  DLIList<Lump*> lumps;
  bodysm_ptr->lumps( lumps );

  DLIList<SenseEntity*> dead_covols, new_covols(lumps.size());
  lumps.reset();
  for( i = lumps.size(); i--; )
  {
    RefVolume* vol = make_RefVolume( lumps.get_and_step(), tmp_vol_modified );
    if( !vol )
    {
      PRINT_ERROR("RefVolume creation failed in GeometryQueryTool::make_Body\n");
      PRINT_ERROR("Aborting construction of Body %d\n", body->id() );
      assert( !!vol );
    }

    if (tmp_vol_modified)
      vol_modified = true;

    SenseEntity* covol = vol->find_sense_entity( body );
    if (!covol)
    {
      covol = new CoVolume( vol );
      body_modified = true;
    }

    new_covols.append(covol);
  }

  body->set_sense_entity_list( new_covols, dead_covols );

  while (dead_covols.size())
  {
    destroy_dead_entity(dead_covols.pop());
    body_modified = true;
  }

  PRINT_DEBUG_3("Constructed Body %d from BodySM in %f seconds.\n",
    body->id(), timer.cpu_secs() );

  if( body_created )
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::FREE_REF_ENTITY_GENERATED, body));
    CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_CREATED, body);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }
  else if( body_modified )
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_MODIFIED, body));
    CGMHistory::Event evt(CGMHistory::TOPOLOGY_CHANGED, body);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }
  else if( vol_modified )
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::GEOMETRY_MODIFIED, body));
    CGMHistory::Event evt(CGMHistory::GEOMETRY_CHANGED, body);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }

  PRINT_DEBUG_3("Updated graphics for Body %d in %f seconds.\n",
    body->id(), timer.cpu_secs() );

  return body;
}

//-------------------------------------------------------------------------
// Purpose       : Build/update RefVolume topology
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/21/02
//-------------------------------------------------------------------------
RefVolume* GeometryQueryTool::make_RefVolume( Lump* lump_ptr,
                                              bool& volume_modified ) const
{
  int i;

  DLIList<ShellSM*> shellsms;
  DLIList<GroupingEntity*> old_shells, new_shells;

  bool volume_created = false;
  volume_modified = false;
  RefVolume* volume = CAST_TO( lump_ptr->topology_entity(), RefVolume );
  if( ! volume )
  {
    volume = RefEntityFactory::instance()->construct_RefVolume(lump_ptr);
    volume_created = true;
  }
  if( ! volume )
  {
    PRINT_ERROR("Volume creation failed in GeometryQueryTool::make_RefVolume\n");
    assert( !!volume );
    return 0;
  }

  lump_ptr->shellsms( shellsms );
  shellsms.reset();
  for( i = shellsms.size(); i--; )
  {
    bool shell_modified;

      // Don't use a shell if it is in a different RefVolume.
      // Can't move shells or other RefVolume may appear  unmodified
      // in later call to make_RefVolume.
// This breaks things -- live without it for now
//    ShellSM* shellsm = shellsms.get_and_step();
//    Shell* shell = dynamic_cast<Shell*>(shellsm->topology_entity());
//    if (shell && shell->get_basic_topology_entity_ptr() &&
//                 shell->get_basic_topology_entity_ptr() != volume)
//      shell->bridge_manager()->remove_bridge(shellsm);

      // Build/update Shell
    Shell* shell = make_Shell( shellsms.get_and_step(), shell_modified );

    if( !shell )
    {
      PRINT_ERROR("Shell creation failed in GeometryQueryTool::make_RefVolume\n");
      PRINT_ERROR("Aborting construction of Volume %d\n", volume->id() );
      assert( !!shell );
      return 0;
    }
    if( shell_modified )
      volume_modified = true;

      // If the shell is new (not already part of RefVolume) then
      // we will be modifying the RefVolume when we add it later.
    BasicTopologyEntity* parent = shell->get_basic_topology_entity_ptr();
    if (parent != volume)
    {
      if (parent)
        parent->remove_grouping_entity(shell);
      volume_modified = true;
    }

    new_shells.append( shell );
  }

  volume->set_grouping_entity_list( new_shells, old_shells );
  while( old_shells.size() )
  {
    GroupingEntity* shell = old_shells.pop();
    destroy_dead_entity( shell );
    volume_modified = true;
  }

  if (volume->deactivated())
    volume_modified = false;
  else if(!volume_created && volume_modified )
  {

#ifndef ALPHA_LAYTRACKS3D
	AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_MODIFIED, volume));
	CGMHistory::Event evt(CGMHistory::TOPOLOGY_CHANGED, volume);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
#endif
  }

  return volume;
}

//-------------------------------------------------------------------------
// Purpose       : Build/update Shell topology
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/21/02
//-------------------------------------------------------------------------
Shell* GeometryQueryTool::make_Shell(ShellSM *shellsm_ptr, bool& shell_modified ) const
{
  Shell* shell = CAST_TO( shellsm_ptr->topology_entity(), Shell );
  shell_modified = false;
  if( !shell )
    shell = new Shell(shellsm_ptr);
  assert( shell != NULL );

  DLIList<CoFace*> face_cofaces(2);
  //shell->co_faces( old_cofaces );
  DLIList<SenseEntity*> new_cofaces, old_cofaces;

  DLIList<Surface*> surfaces;
  shellsm_ptr->surfaces( surfaces );

  int i;
  surfaces.reset();
  for( i = surfaces.size(); i--; )
  {
    Surface* surface = surfaces.get_and_step();
    RefFace* face = make_RefFace(surface);
    if( !face )
    {
      PRINT_ERROR("Surface creation failed in GeometryQueryTool::make_Shell\n");
      assert( !!face );
      return 0;
    }

    CubitSense sense = surface->get_shell_sense( shellsm_ptr );
    if( surface->bridge_sense() == CUBIT_REVERSED )
      sense = CubitUtil::opposite_sense( sense );

      // If sense is CUBIT_UNKNOWN, we have a two-sided
      // face and need to make two CoFaces.
    int num_cofaces = 1;
    if( sense == CUBIT_UNKNOWN )
    {
      DLIList<ShellSM*> surf_shells(2);
      surface->shellsms( surf_shells );
      assert( surf_shells.is_in_list(shellsm_ptr) );
      num_cofaces = 2;
      sense = CUBIT_FORWARD;
    }

      // If two CoFaces, loop once for each
    face_cofaces.clean_out();
    face->co_faces(face_cofaces);
    for( int k = 0; k < num_cofaces; k++ )
    {

      CoFace* coface = 0;
      for( int j = face_cofaces.size(); j--; )
      {
        CoFace* tmp_coface = face_cofaces.get();
        if( tmp_coface->get_sense() == sense &&
            tmp_coface->get_shell_ptr() == shell )
        {
          coface = tmp_coface;
          break;
        }
      }

      if (!coface)
      {
        coface = new CoFace();
        coface->attach_basic_topology_entity( face );
        shell_modified = true;
        if( coface->get_sense() != sense )
          coface->set_sense( sense );
      }

      new_cofaces.append(coface);

        // If we loop to create a second CoFace between
        // the same Surface-Shell pair, create it with
        // the opposite sense.
      sense = CubitUtil::opposite_sense(sense);
    }
  }

  shell->set_sense_entity_list( new_cofaces, old_cofaces );

    // Destroy any unwanted CoFaces (those remaining
    // in old_cofaces)
  while( old_cofaces.size() )
  {
    CoFace* dead_coface = dynamic_cast<CoFace*>(old_cofaces.pop());
    assert(!dead_coface->get_shell_ptr());
    destroy_dead_entity(dead_coface);
    shell_modified = true;
  }

  return shell;
}


//-------------------------------------------------------------------------
// Purpose       : Find the CoEdgeSM in the passed Surface that should
//                 be merged with the passed CoEdgeSM.
//
// Special Notes : Should be O(1) in the average case.  See comments below.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/21/02
//-------------------------------------------------------------------------
CoEdgeSM* GeometryQueryTool::find_merged_coedgesm( Surface* on_this_surf,
                                      CoEdgeSM* merged_with_this ) const
{
    // Get the relative sense of the merged surfaces
  DLIList<Surface*> coedge_surfaces(1);
  merged_with_this->surfaces( coedge_surfaces );
  if( coedge_surfaces.size() != 1 )
    return 0;
  bool surfaces_reversed =
    (coedge_surfaces.get()->bridge_sense() != on_this_surf->bridge_sense());

    // Get the Curve merge partner from the CoEdgeSM merge partner
  DLIList<Curve*> coedge_curves(1);
  merged_with_this->curves( coedge_curves );
  if( coedge_curves.size() != 1 )
  {
    assert(coedge_curves.size() == 1);
    return 0;
  }
  Curve* merged_curve = coedge_curves.get();

    // If the curve doesn't have a bridge manager, it is
    // not merged.
  if (!merged_curve->bridge_manager())
    return 0;

    // Get all curves merged with the merged_curve
  int num_bridges = merged_curve->bridge_manager()->number_of_bridges();
  DLIList<TopologyBridge*> curve_bridges( num_bridges );
  merged_curve->bridge_manager()->get_bridge_list( curve_bridges );

    // Loop through Curves, looking for the Curve with the
    // parent CoEdgeSM we want.  This loop should be
    // roughly O(1) for manifold geometry (and most non-
    // manifold geometry), where n is the number of curves
    // merged together in the RefEdge.  For manifold
    // geometry, there are never more than 2 CoEdgeSMs per
    // Curve.  For most models, n is <= 4, so this entire
    // function can be considered O(1) for most models.
  DLIList<CoEdgeSM*> curve_coedges(2);
  curve_bridges.reset();
  for( int i = curve_bridges.size(); i--; )
  {
    TopologyBridge* curve = curve_bridges.get_and_step();

    curve_coedges.clean_out();
    curve->coedgesms( curve_coedges );
    curve_coedges.reset();
    for( int j = curve_coedges.size(); j--; )
    {
      CoEdgeSM* coedgesm = curve_coedges.get_and_step();

        // Get parent Surface of CoEdgeSM
      coedge_surfaces.clean_out();
      coedgesm->surfaces( coedge_surfaces );
      assert( coedge_surfaces.size() == 1 );

        // Is the parent Surface the one we are looking for?
      if( coedge_surfaces.get() == on_this_surf )
      {
          // Curve may occur in surface twice for non-manifold
          // topology.  We need to make sure we have the CoEdgeSM
          // with the correct sense.

          // Are the curves reversed wrt to each other?
        bool curves_reversed =
          (merged_curve->bridge_sense() != curve->bridge_sense());
          // Should this coedgesm be reversed wrt to the merge partner?
        bool should_be_reversed = (curves_reversed != surfaces_reversed);
          // Are the coedgesm actually reversed wrt to each other?
        bool are_reversed = (merged_with_this->sense() != coedgesm->sense());

          // Do the coedges have the appropriate senses
        if( should_be_reversed == are_reversed )
        {
          return  coedgesm;
        }
      }
    } // for(j = curve_coedgesms)
  } // for(i = curve_coedges)

  return 0;
}

//-------------------------------------------------------------------------
// Purpose       : Build/update RefFace topology
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/21/02
//-------------------------------------------------------------------------
RefFace* GeometryQueryTool::make_RefFace( Surface* surface_ptr ) const
{
  int i;
  DLIList<LoopSM*> loopsms;
  DLIList<CoEdgeSM*> coedgesms;
  DLIList<Curve*> curves(1);

    // Get/construct RefFace for surface.
  bool face_created = false;
  bool face_modified = false;
  bool face_modified2 = false;
  RefFace* face = CAST_TO( surface_ptr->topology_entity(), RefFace );
  if( !face )
  {
    face = dynamic_cast<RefFace*>(check_mergeable_refentity(surface_ptr));

    //Including this call in here in case the surface a composite and
    //hidden curves need to be removed in the graphics.
    if( face )
    {
      face_modified2 = true;
    }
  }
  if( !face )
  {
    face = RefEntityFactory::instance()->construct_RefFace(surface_ptr);
    face_created = true;
  }
  assert(!!face);
  bool merged = face->bridge_manager()->number_of_bridges() > 1;

    // First construct all the RefEdges.  If the Curves
    // were merged, the RefEdges will be constructed merged.
    // We need the merged RefEdges to determined which
    // CoEdges and Loops need to be merged.  We could do this
    // with one loop by getting the curves directly from the
    // surface.  However, that might result in different
    // RefEdge IDs, so do it with nested loops like the old code.
  loopsms.clean_out();
  surface_ptr->loopsms( loopsms );
  loopsms.reset();
  for( i = loopsms.size(); i--; )
  {
    LoopSM* loopsm = loopsms.get_and_step();
    coedgesms.clean_out();
    loopsm->coedgesms(coedgesms);
    coedgesms.reset();
    for( int j = coedgesms.size(); j--; )
    {
      CoEdgeSM* coedgesm = coedgesms.get_and_step();
      curves.clean_out();
      coedgesm->curves(curves);
      if( curves.size() != 1 )
      {
        PRINT_ERROR("Invalid SolidModel topology encountered in "
                    "GeometryQueryTool::make_RefFace.  CoEdgeSM in "
                    "Surface %d has %d curves.\n", face->id(), curves.size());
        if( curves.size() == 0 )
          continue;
      }

      RefEdge* edge = make_RefEdge( curves.get() );
      if( !edge )
      {
        PRINT_ERROR("Failed to construct RefEdge in Surface %d\n",face->id());
      }
    }
  }

    // If the Surface was merged, make sure it can remain merged
  if (merged && !make_merged_RefFace(surface_ptr))
    merged = false;

    // If !merged, then if necessary unmerge the RefFace
  if( !merged && face->bridge_manager()->number_of_bridges() > 1 )
  {
    PRINT_WARNING("Unmerging Surface %d.\n", face->id() );

      // unmerge
    face->bridge_manager()->remove_bridge(surface_ptr);

      // reset bridge sense on unmerged surface
    if (surface_ptr->bridge_sense() == CUBIT_REVERSED)
      surface_ptr->reverse_bridge_sense();

      // if the first bridge in the old refface is now reversed,
      // flip the refface.
    if (face->get_surface_ptr()->bridge_sense() == CUBIT_REVERSED)
    {
      face->bridge_manager()->reverse_bridge_senses();
      face->reverse_topology();
    }

      // construct new refface for unmerged surface.
    face = RefEntityFactory::instance()->construct_RefFace(surface_ptr);
    face_created = true;
  }

    // If the sense of the Surface with respect to the RefFace is
    // reversed, need to construct Loops and CoEdges reversed.
  bool need_to_reverse = (surface_ptr->bridge_sense() == CUBIT_REVERSED);

    // Construct LoopSMs and CoEdgeSMs, and attach
  loopsms.clean_out();
  surface_ptr->loopsms( loopsms );
  loopsms.reset();
  DLIList<SenseEntity*> coedges, old_coedges;
  DLIList<GroupingEntity*> loops(loopsms.size()), dead_loops;
  for( i = loopsms.size(); i--; )
  {
      // Get/construct Loop
    LoopSM* loopsm = loopsms.get_and_step();
    Loop* loop = dynamic_cast<Loop*>(loopsm->topology_entity());
    if( loop && !merged && loop->bridge_manager()->number_of_bridges() > 1 )
    {
      loop->bridge_manager()->remove_bridge( loopsm );
      loop = new Loop( loopsm );
      face_modified = true;
    }
    else if( !loop )
    {
      loop = new Loop( loopsm );
      face_modified = true;
    }

      // If loop is already attached to some other RefFace
      // we cannot use it.  NOTE:  Removing the Loop from
      // the other RefFace is NOT an option, even if the
      // loop has only our loopsm.  If we were to remove it
      // we could potentially miss the fact that the other
      // RefFace was modified when doing make_RefFace for
      // RefFace.
    RefFace* loop_face = loop->get_ref_face_ptr();
    if( loop_face && loop_face != face )
    {
      loop->bridge_manager()->remove_bridge( loopsm );
      loop = new Loop( loopsm );
      face_modified = true;
    }
    loops.append(loop);

      // Get CoEdges in the appropriate order
    coedgesms.clean_out();
    loopsm->coedgesms( coedgesms );
    if( coedgesms.size() == 0 )
    {
      PRINT_ERROR("Encountered Loop with no CoEdges in Surface %d.\n",
        face->id() );
      continue;
    }
    if( need_to_reverse )
      coedgesms.reverse();
    coedgesms.reset();

      // Consturct/update each CoEdge
    int j;
    coedges.clean_out();
    for( j = coedgesms.size(); j--; )
    {
      CoEdgeSM* coedgesm = coedgesms.get_and_step();
      CoEdge* coedge = dynamic_cast<CoEdge*>(coedgesm->topology_entity());

        // If the CoEdge is merged but the surface isn't, unmerge
        // the CoEdge
      if( coedge && !merged && coedge->bridge_manager()->number_of_bridges() > 1 )
      {
        coedge->bridge_manager()->remove_bridge( coedgesm );
        coedge = 0;
      }

        // If CoEdge belongs to a loop in some other face, don't
        // use it.
      if( coedge && coedge->get_loop_ptr() && coedge->get_loop_ptr() != loop )
      {
        coedge->bridge_manager()->remove_bridge( coedgesm );
        coedge = 0;
      }

        // If the CoEdgeSM doesn't already have an owning CoEdge
        // (or we unmerged) create a new CoEdge
      if( !coedge )
      {
        coedge = new CoEdge();
        coedge->bridge_manager()->add_bridge( coedgesm );
        face_modified = true;
      }
      coedges.append( coedge );

        // Find the RefEdge to attach to
      curves.clean_out();
      coedgesm->curves( curves );
      if( curves.size() == 0 )
        continue; // error message should have been printed already

      Curve* curve = curves.get();
      RefEdge* edge = dynamic_cast<RefEdge*>(curve->topology_entity());
      assert( edge != NULL );

        // Figure out what the CoEdge sense should be.
        // The CoEdge sense should be the product ("product"
        // as in multiply) of three
        // senses: 1) the sense of the surface wrt its RefFace,
        // 2) the sense of the curve wrt its RefEdge, and 3)
        // the sense of the CoEdgeSM.
      bool curve_reversed = (CUBIT_REVERSED == curve->bridge_sense());
      bool reverse = (curve_reversed != need_to_reverse);
      CubitSense sense = coedgesm->sense();
      if( reverse )
        sense = CubitUtil::opposite_sense(sense);
      if( coedge->get_sense() != sense )
        coedge->set_sense( sense );

        // Attach RefEdge to CoEdge (if it isn't already)
      RefEdge* coedge_edge = coedge->get_ref_edge_ptr();
      if (coedge_edge != edge)
      {
        if (coedge_edge)
        {
          coedge_edge->remove_sense_entity(coedge);
          destroy_dead_entity(coedge_edge);
          face_modified = true;
        }
        coedge->attach_basic_topology_entity(edge);
      }
    }

      // Attach CoEdges to Loop
    old_coedges.clean_out();
    loop->set_sense_entity_list(coedges, old_coedges);

      // Clean up dead coedges
    while (old_coedges.size())
    {
      CoEdge* dead_coedge = dynamic_cast<CoEdge*>(old_coedges.pop());
      assert(!dead_coedge->get_loop_ptr());
      destroy_dead_entity(dead_coedge);
      face_modified = true;
    }
  } // end for( i == loopsms )

    // Attach loops
  face->set_grouping_entity_list( loops, dead_loops );

    // Remove any dead loops
  while (dead_loops.size())
  {
    GroupingEntity* loop = dead_loops.pop();
    destroy_dead_entity(loop);
    face_modified = true;
  }

  if( !face_created && face_modified && !face->deactivated() )
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::GEOMETRY_TOPOLOGY_MODIFIED, face));
    CGMHistory::Event evt2(CGMHistory::TOPOLOGY_CHANGED, face);
    const_cast<CGMHistory&>(mHistory).add_event(evt2);
    CGMHistory::Event evt(CGMHistory::GEOMETRY_CHANGED, face);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }
  else if(face_modified2)
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_MODIFIED, face));
    CGMHistory::Event evt(CGMHistory::TOPOLOGY_CHANGED, face);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }

  return face;
}


//-------------------------------------------------------------------------
// Purpose       : Check if Surface can remain merged.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/12/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::make_merged_RefFace(Surface* surface_ptr) const
{
    // Get the RefFace
  RefFace* face = dynamic_cast<RefFace*>(surface_ptr->topology_entity());
  if (!face)  // isn't merged
    return CUBIT_FAILURE;

    // Get a Surface to compare to.
  int num_bridges = face->bridge_manager()->number_of_bridges();
  if (num_bridges < 2) // isn't merged
    return CUBIT_FAILURE;

  // Get some other Surface in merged RefFace to compare
  // topology with.
  DLIList<TopologyBridge*> bridge_list(num_bridges);
  face->bridge_manager()->get_bridge_list(bridge_list);
  bridge_list.reset();
  if (bridge_list.get() == surface_ptr)
    bridge_list.step();
  TopologyBridge* merged_bridge = bridge_list.get();
  Surface* merged_surface = dynamic_cast<Surface*>(merged_bridge);
  assert(merged_surface && merged_surface != surface_ptr);

    // Get the relative sense of the Surfaces
  bool surfaces_reversed =
    (surface_ptr->bridge_sense() != merged_surface->bridge_sense());

    // Get LoopSMs for both surfaces
  DLIList<LoopSM*> merged_loops, loopsms;
  surface_ptr->loopsms(loopsms);
  merged_surface->loopsms(merged_loops);
  if( merged_loops.size() != loopsms.size() )
    return CUBIT_FAILURE;

    // For each LoopSM on the merged Surface, find the
    // corresponding LoopSM on this Surface and merge.
  DLIList<CoEdgeSM*> merged_coedges, coedgesms;
  DLIList<LoopSM*> coedge_loops(1);
  DLIList<Curve*> curves(1);
  merged_loops.reset();
  for( int i = merged_loops.size(); i--; )
  {
    LoopSM* merged_loop = merged_loops.get_and_step();

      // Find the LoopSM by choosing a CoEdgeSM from the
      // merged LoopSM, and finding corresponding CoEdgeSM
      // on the Surface we are building.

    merged_coedges.clean_out();
    merged_loop->coedgesms( merged_coedges );
    merged_coedges.reset();
    CoEdgeSM* merged_coedge = merged_coedges.get();

    CoEdgeSM* coedgesm = find_merged_coedgesm( surface_ptr, merged_coedge );
    if( !coedgesm )
      return CUBIT_FAILURE;

     // Get the LoopSM in surface_ptr from the CoEdgeSM
    coedge_loops.clean_out();
    coedgesm->loopsms( coedge_loops );
    assert( coedge_loops.size() == 1 );
    LoopSM* loopsm = coedge_loops.get();

      // Get all the CoEdgeSMs on the LoopSM we are going to merge.
    coedgesms.clean_out();
    loopsm->coedgesms( coedgesms );
    if( coedgesms.size() != merged_coedges.size() )
      return CUBIT_FAILURE;

      // Get the CoEdge list in the appropriate order and position
    if( surfaces_reversed )
      coedgesms.reverse();
    coedgesms.move_to(coedgesm);
    assert( coedgesms.get() == coedgesm );

      // Now check and merge CoEdgeSM pairs
    for( int j = coedgesms.size(); j--; )
    {
      coedgesm = coedgesms.get_and_step();
      merged_coedge = merged_coedges.get_and_step();

        // Get Curve in surface_ptr
      curves.clean_out();
      coedgesm->curves(curves);
      if( curves.size() != 1 )
      {
        PRINT_ERROR("Invalid SolidModel topology encountered in "
                    "GeometryQueryTool::make_RefFace.  CoEdgeSM in "
                    "has %d curves.\n", curves.size());
        return CUBIT_FAILURE;
      }
      Curve* curve = curves.get();

        // Get Curve in merged_surface
      curves.clean_out();
      merged_coedge->curves(curves);
      if( curves.size() != 1 )
      {
        PRINT_ERROR("Invalid SolidModel topology encountered in "
                    "GeometryQueryTool::make_RefFace.  CoEdgeSM in "
                    "has %d curves.\n", curves.size());
        return CUBIT_FAILURE;
      }
      Curve* merged_curve = curves.get_and_step();

        // Check if curves are merged
      if( merged_curve->bridge_manager() != curve->bridge_manager() )
        return CUBIT_FAILURE;

        // Check that relative sense of CoEdgeSMs is correct
      bool curves_reversed =
        (merged_curve->bridge_sense() != curve->bridge_sense());
        // Should this coedgesm be reversed wrt to the merge partner?
      bool should_be_reversed = (curves_reversed != surfaces_reversed);
        // Are the coedgesm actually reversed wrt to each other?
      bool are_reversed = (merged_coedge->sense() != coedgesm->sense());
      if( should_be_reversed != are_reversed )
        return CUBIT_FAILURE;

        // Merge the CoEdgeSMs
      CoEdge* coedge = dynamic_cast<CoEdge*>(coedgesm->topology_entity());
      CoEdge* merge  = dynamic_cast<CoEdge*>(merged_coedge->topology_entity());

        // If this coedgesm is not the primary coedgesm in the
        // CoEdge, check to make sure the primary CoEdgeSM's curve
        // is merged with our curve.
      if (coedge && coedge->bridge_manager()->number_of_bridges() > 1
           && coedge->bridge_manager()->topology_bridge() != coedgesm
           && coedge->bridge_manager()->topology_bridge() != merged_coedge )
      {
        curves.clean_out();
        coedge->bridge_manager()->topology_bridge()->curves(curves);
        assert(curves.size() == 1);
        if (curves.get()->owner() != curve->owner())
        {
          coedge->bridge_manager()->remove_bridge(coedgesm);
          assert(coedge->get_parents());
          coedge = 0;
        }
      }
      if (merge && merge->bridge_manager()->number_of_bridges() > 1
           && merge->bridge_manager()->topology_bridge() != coedgesm
           && merge->bridge_manager()->topology_bridge() != merged_coedge )
      {
        curves.clean_out();
        merge->bridge_manager()->topology_bridge()->curves(curves);
        assert(curves.size() == 1);
        if (curves.get()->owner() != curve->owner())
        {
          merge->bridge_manager()->remove_bridge(merged_coedge);
          assert(merge->get_parents());
          merge = 0;
        }
      }

        // If the two CoEdgeSM's aren't in the same CoEdge,
        // merge them.
      if (merge && merge != coedge)
      {
        if (coedge)
        {
          coedge->bridge_manager()->remove_bridge(coedgesm);
          destroy_dead_entity(coedge);
        }
        merge->bridge_manager()->add_bridge(coedgesm);
      }
        // If for some reason (Composite curve creatation does
        // this) neither CoEdgeSM has a CoEdge, create one and
        // re-merge them.
      else if (!coedge && !merge)
      {
        coedge = new CoEdge(coedgesm);
        coedge->bridge_manager()->add_bridge(merged_coedge);
      }
      else if(!coedge)
      {
        merge->bridge_manager()->add_bridge(coedgesm);
      }
      else if(!merge)
      {
        coedge->bridge_manager()->add_bridge(merged_coedge);
      }
    } // for( j = coedgesms )


    Loop* loop = dynamic_cast<Loop*>(loopsm->topology_entity());
    Loop* merged = dynamic_cast<Loop*>(merged_loop->topology_entity());
      // If the two LoopSMs aren't in the same Loop, merge them
    if( merged != 0 && merged != loop )
    {
      if( loop )
        loop->bridge_manager()->remove_bridge(loopsm);
      merged->bridge_manager()->add_bridge(loopsm);
    }
    else if( !loop && !merged )
    {
      loop = new Loop;
      loop->bridge_manager()->add_bridge(merged_loop);
    }
    else if( !loop )
    {
      merged->bridge_manager()->add_bridge(loopsm);
    }
    else if( !merged )
    {
      loop->bridge_manager()->add_bridge(merged_loop);
    }
  } // for( i = merged_loops )

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get/construct/update a RefEdge from a Curve
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/22/02
//-------------------------------------------------------------------------
RefEdge* GeometryQueryTool::make_RefEdge( Curve* curve_ptr ) const
{
  int i;

    // Get/construct RefEdge
  bool edge_created = false;
  bool edge_modified = false;
  RefEdge* edge = dynamic_cast<RefEdge*>(curve_ptr->topology_entity());
  if( !edge )
    edge = dynamic_cast<RefEdge*>(check_mergeable_refentity(curve_ptr));
  if( !edge )
  {
    edge = RefEntityFactory::instance()->construct_RefEdge(curve_ptr);
    edge_created = true;
  }
  assert(edge != NULL);
  bool reversed = (CUBIT_REVERSED == curve_ptr->bridge_sense());
  bool merged = edge->bridge_manager()->number_of_bridges() > 1;

    // Construct RefVertices
  DLIList<TBPoint*> points(2);
  curve_ptr->points( points );
  if( reversed )
    points.reverse();

  if( points.size() > 2 || points.size() < 1 )
  {
    PRINT_ERROR("Invalid SolidModel topology encountered in "
                "GeometryQueryTool::make_RefEdge.  Curve %d "
                "has %d vertices.\n", edge->id(), points.size() );
    assert( points.size() == 1 || points.size() == 2 );
    if( points.size() == 0 )
      return edge;
  }

  points.reset();
  RefVertex* start_vertex = make_RefVertex( points.get() );
  points.last();
  RefVertex* end_vertex = make_RefVertex( points.get() );

    // If curves are merged, make sure vertices are merged.
  if( merged )
  {
    int num_bridges = edge->bridge_manager()->number_of_bridges();
    DLIList<TopologyBridge*> bridge_list( num_bridges );
    edge->bridge_manager()->get_bridge_list( bridge_list );

    bridge_list.reset();
    for( i = bridge_list.size(); i--; )
      if( bridge_list.get() != curve_ptr )
        break;
      else
        bridge_list.step();

      // Get points on some other Curve that this curve is merged with.
    Curve* other_curve = dynamic_cast<Curve*>(bridge_list.get());
    assert( other_curve != curve_ptr) ;
    points.clean_out();
    other_curve->points(points);
    if( other_curve->bridge_sense() == CUBIT_REVERSED )
      points.reverse();

      // Check that points merged.
    points.reset();
    if( points.get()->topology_entity() != start_vertex )
      merged = false;
    points.last();
    if( points.get()->topology_entity() != end_vertex )
      merged = false;

    //perhaps the bridge sense just needs to be swapped 
    if( merged == false )
    {
      points.reverse();
      points.reset();
      TBPoint *point1 = points.get_and_step();
      TBPoint *point2 = points.get_and_step();
      if( point1->topology_entity() == start_vertex &&
          point2->topology_entity() == end_vertex )
      {
        merged = true;
        curve_ptr->reverse_bridge_sense();
      }
    }
  }
    // Unmerge the curve, if necessary.
  if( !merged && edge->bridge_manager()->number_of_bridges() > 1 )
  {
    PRINT_WARNING("Unmerging Curve %d\n", edge->id() );

      // unmerge
    edge->bridge_manager()->remove_bridge( curve_ptr );

      // reset bridge sense on un-merged curve
    if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    {
      curve_ptr->reverse_bridge_sense();
    }

      // if first curve in old edge now has reversed sense,
      // reverse the old edge also.
    //if (edge->get_curve_ptr()->bridge_sense() == CUBIT_REVERSED)
    //{
    //  edge->bridge_manager()->reverse_bridge_senses();
    //  edge->reverse_topology();
    //}

      // Create new RefEdge for un-merged curve
    edge = RefEntityFactory::instance()->construct_RefEdge( curve_ptr );
    PRINT_WARNING("Creating edge %d that was supposed to be merged\n", edge->id() );
    edge_created = true;

      // If we had swapped the vertex order because the curve
      // was reversed, switch them back because we unmerged so
      // we're going to create the new curve with a forward sense.
    if( reversed )
      std::swap(start_vertex,end_vertex);
  }

    // Get/construct Chain
  Chain* chain = edge->get_chain_ptr();
  if( !chain )
  {
    chain = new Chain();
    edge->add_grouping_entity( chain );
  }

    // Get any existing CoVertices
  DLIList<CoVertex*> co_vertices(2);
  chain->co_vertices( co_vertices );
  if( co_vertices.size() > 2 )
  {
    PRINT_ERROR("In GeometryQueryTool::make_RefEdge, encountered "
                "a Chain with %d CoVertices on Curve %d\n",
                co_vertices.size(), edge->id());
    assert( co_vertices.size() <= 2 );
  }

    // Now construct CoVertices if necessary, and attach
    // RefVertices.
  DLIList<RefVertex*> vertices(1);
  CoVertex* prev = 0;
    // First iteration is for start RefVertex
  RefVertex* vertex = start_vertex;
  co_vertices.reset();
    // Loop twice, once for each CoVertex.
  for( i = 0; i < 2; i++ )
  {
      // Construct CoVertex if one does not already exist
    CoVertex* cvtx = 0;
    if( co_vertices.size() > i )
    {
      cvtx = co_vertices.get_and_step();
    }
    else
    {
      cvtx = new CoVertex();
      chain->add_sense_entity( cvtx , prev );
    }
    prev = cvtx;

      // Get existing RefVertex from CoVertex, if it has one.
    vertices.clean_out();
    cvtx->ref_vertices( vertices );
    if( vertices.size() > 1 )
    {
      PRINT_ERROR("In GeometryQueryTool::make_RefEdge, encountered "
                  "a CoVertex with %d Vertices on Curve %d\n",
                   vertices.size(), edge->id());
      assert( vertices.size() <= 1 );
      vertices.reset();
    }
    RefVertex* cvtx_vertex = vertices.size() ? vertices.get() : 0;

      // Attach RefVertex to CoVertex if necessary
    if( cvtx_vertex != vertex )
    {
        // Wrong RefVertex.  Unhook from that one.
      if( cvtx_vertex )
      {
        cvtx_vertex->remove_sense_entity(cvtx);
        if (NULL == cvtx_vertex->get_point_ptr())
          destroy_dead_entity( cvtx_vertex );
      }

      cvtx->attach_basic_topology_entity( vertex );
      edge_modified = true;
    }

      // Second iteration is for end RefVertex
    vertex = end_vertex;
  }

  if( !edge_created && edge_modified && !edge->deactivated() )
  {
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::GEOMETRY_TOPOLOGY_MODIFIED, edge));
    CGMHistory::Event evt2(CGMHistory::TOPOLOGY_CHANGED, edge);
    const_cast<CGMHistory&>(mHistory).add_event(evt2);
    CGMHistory::Event evt(CGMHistory::GEOMETRY_CHANGED, edge);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }

  return edge;
}


RefEntity *GeometryQueryTool::check_mergeable_refentity(GeometryEntity *bridge) const
{
    // Disable for surfaces until we can correctly handle
    // surfaces that have virtual curves.
    // J.Kraftcheck  Apr.29,2002
//  if( dynamic_cast<Surface*>(bridge) )
//    return 0;

  if (CUBIT_FALSE == CGMApp::instance()->attrib_manager()->auto_actuate_flag(CA_MERGE_PARTNER))
    return NULL;

    // checks the topology bridge for a merge attribute, which might lead to another
    // refentity
  DLIList<CubitSimpleAttrib> csa_list;
  const char* name = CGMApp::instance()->attrib_manager()->att_internal_name(CA_MERGE_PARTNER);
  bridge->get_simple_attribute(name, csa_list);

  RefEntity *re_ptr = NULL;

  if(csa_list.size() > 0)
  {
    const CubitSimpleAttrib &csa_ptr = csa_list.size() ? csa_list.get() : CubitSimpleAttrib();

      // from the csa, call CAMP to get the partner, if one exists
      // (and return it)
    int merge_id = csa_ptr.int_data_list()[0];

    if( importingSolidModel && CUBIT_FALSE == GSaveOpen::performingUndo && 
        !mergeGloballyOnImport )
    {
      //if the merge id isn't already on some other entity already imported,
      //this entity doesn't have anybody to merge with...return NULL
      if( !uidsOfImportingEnts.is_in_list( merge_id) )
      {
        GeometryQueryTool::uidsOfImportingEnts.append( merge_id );
        return NULL;
      }
    }

    bool unique_append = false;
    if( GeometryQueryTool::trackMergedAwayEnts )
      unique_append = GeometryQueryTool::uidsOfImportingEnts.append_unique( merge_id );

    ToolDataUser *tdu = TDUniqueId::find_td_unique_id(merge_id);
    TopologyEntity *te = dynamic_cast<TopologyEntity*>(tdu);

    CubitSense sense = CAMergePartner::get_bridge_sense( csa_ptr );

    //We found the merge partner.....
    if (te != NULL)
    {
      //it's the first time you found this uid and you found an entity
      //
      if(GeometryQueryTool::trackMergedAwayEnts && unique_append )
        GeometryQueryTool::entitiesMergedAway++;

        // assume this merge attrib will be actuated, so remove the csa
      bridge->remove_simple_attribute_virt(csa_ptr);

        // now do the actual merge
      re_ptr = dynamic_cast<RefEntity*> (te);

      //compare relative sense.....reverse bridges as necessary
      Surface *surface = CAST_TO( bridge, Surface );
      Curve *curve= CAST_TO( bridge, Curve );
      CubitSense rel_sense = CUBIT_FORWARD;

      CubitBoolean is_survivor = CAMergePartner::is_survivor( csa_ptr );
      if( surface )
      {
        TopologyBridge *other_bridge = te->bridge_manager()->topology_bridge();
        Surface *other_surface = CAST_TO( other_bridge, Surface );
        rel_sense = relative_sense( surface, other_surface );
      }
      else if( curve )
      {
        TopologyBridge *other_bridge = te->bridge_manager()->topology_bridge();
        Curve *other_curve = CAST_TO( other_bridge, Curve );
        rel_sense = curve->relative_sense( other_curve );
      }

      //non-surviving merged entity is here....reverse
      //its bridge sense so that it will merge
      if( !is_survivor && sense != CUBIT_UNKNOWN && rel_sense == CUBIT_REVERSED )
        bridge->reverse_bridge_sense();
      //surviving merged entity is now being restored....
      //reverse the ref entity and make this entity the primary one
      if( is_survivor )
      {
        //If the passed in new bridge is the real survivor of the merge,
        //it should be the first bridge in the bridge manager.  The first
        //bridge in the bridge manager is always FORWARD wrt the TE.
        // If this new bridge is reversed wrt the current first bridge,
        // reverse all the other bridges so that they will all merge
        //successfully with this new first bridge.  Also add this new bridge
        //as the first bridge, since it was the original survivor.
        if( rel_sense == CUBIT_REVERSED )
        {
          te->reverse_topology();

          //reverse the sense of the bridges of this refentity
          DLIList<TopologyBridge*> bridge_list;
          te->bridge_manager()->get_bridge_list( bridge_list );

          int kk;
          for( kk=bridge_list.size(); kk--; )
          {
            TopologyBridge *tmp_bridge = bridge_list.get_and_step();
            tmp_bridge->reverse_bridge_sense();
          }
        } 

        //add this bridge as the primary bridge
        te->bridge_manager()->add_bridge_as_representation(bridge);
      }
      else
        te->bridge_manager()->add_bridge(bridge);
    }

      // Set the merge sense of the bridge, if it is saved in the
      // attribute.
    if( re_ptr && sense == CUBIT_UNKNOWN )
    {
      TopologyBridge* first = te->bridge_manager()->topology_bridge();
      Curve* curve;
      Surface* surface;
      if( (curve = dynamic_cast<Curve*>(bridge) ) != NULL )
      {
        Curve* first_curve = dynamic_cast<Curve*>(first);
        assert(first_curve != NULL);
        sense = first_curve->bridge_sense();
        if( first_curve->relative_sense(curve) == CUBIT_REVERSED )
          sense = CubitUtil::opposite_sense(sense);
      }
      else if( (surface = dynamic_cast<Surface*>(bridge) ) != NULL )
      {
        Surface* first_surf = dynamic_cast<Surface*>(first);
        assert(first_surf != NULL);
        sense = first_surf->bridge_sense();
        if( relative_sense( first_surf, surface ) == CUBIT_REVERSED )
          sense = CubitUtil::opposite_sense(sense);
      }
    }

    int id = CAMergePartner::get_saved_id( csa_ptr );
    if (id)
      bridge->set_saved_id( id );
  }

  return re_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Find relative senes of two Surfaces
//
// Special Notes : Used when CAMergePartner attrib is too old to
//                 contain sense information.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/22/02
//-------------------------------------------------------------------------
CubitSense GeometryQueryTool::relative_sense( Surface* surf1, Surface* surf2 )
{
  CubitVector point1 = surf1->bounding_box().center();
  CubitVector point2 = surf2->bounding_box().center();
  CubitVector point = 0.5 * (point1 + point2);

  CubitVector closest1, closest2, normal1, normal2;
  surf1->closest_point_trimmed( point, closest1 );
  surf2->closest_point( closest1, &closest2, &normal2 );
  surf1->closest_point( closest2, &closest1, &normal1 );

  return normal1 % normal2 < 0 ? CUBIT_REVERSED : CUBIT_FORWARD;
}

RefVertex* GeometryQueryTool::make_RefVertex(TBPoint* point) const
{
     // first, handle making the new RefVertex's; check for
     // existing chain afterwards
   RefVertex* vertex = 0;

     // First check to make sure we haven't already created a RefVertex
     // from this point
   TopologyEntity* topo = point->topology_entity();
   vertex = CAST_TO(topo, RefVertex);
  
   if ( vertex == NULL)
   {
     assert( topo == NULL );
     RefEntity *re_ptr = check_mergeable_refentity(point);
     if (re_ptr != NULL)
     {
       vertex = dynamic_cast<RefVertex*>(re_ptr);
     }
   }
 
   if (vertex == NULL)
   {
       // Create a RefVertex from this TBPoint.
     vertex = RefEntityFactory::instance()->construct_RefVertex(point);
   }
   assert(vertex != NULL);

   return vertex;
}



//-------------------------------------------------------------------------
// Purpose       : Build a free RefFace given a SurfaceSM.
//
// Special Notes : Purpose is to build a free surface.  Entities are
//                 added to the graphics.
//
// Creator       : Steve Storm
//
// Creation Date : 3/27/99
//-------------------------------------------------------------------------
RefFace* GeometryQueryTool::make_free_RefFace(Surface *surface_ptr,
                                              bool is_free_face ) const
{
   RefFace* ref_face_ptr = make_RefFace( surface_ptr );
   if( !ref_face_ptr )
   {
     PRINT_ERROR("Failed to construct free RefFace.\n");
     assert(!!ref_face_ptr);
     return 0;
   }

   // Add the new ref_face to the graphics
   if( is_free_face )
   {
     AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::FREE_REF_ENTITY_GENERATED, ref_face_ptr));
     
     CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_CREATED, ref_face_ptr);
     const_cast<CGMHistory&>(mHistory).add_event(evt);
   }

   return ref_face_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Build a free RefEdge given a CurveSM.
//
// Special Notes : Purpose is to build a free curve.  Entities are
//                 added to the graphics.
//
// Creator       : Steve Storm
//
// Creation Date : 3/27/99
//-------------------------------------------------------------------------
RefEdge* GeometryQueryTool::make_free_RefEdge(Curve *curve_ptr ) const
{
   RefEdge* ref_edge_ptr = make_RefEdge(curve_ptr);
   if( !ref_edge_ptr )
   {
     PRINT_ERROR("Failed to construct free RefEdge.\n");
     assert(!!ref_edge_ptr);
     return 0;
   }

   AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::FREE_REF_ENTITY_GENERATED, ref_edge_ptr));
   
   CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_CREATED, ref_edge_ptr);
   const_cast<CGMHistory&>(mHistory).add_event(evt);

   return ref_edge_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Build a free RefVertex given a PointSM.
//
// Special Notes : Purpose is to build a free vertex.  Vertex
//                 is added to the graphics.
//
// Creator       : Steve Storm
//
// Creation Date : 3/27/99
//-------------------------------------------------------------------------
RefVertex* GeometryQueryTool::make_free_RefVertex(TBPoint *point_ptr ) const
{
   RefVertex* ref_vertex_ptr = make_RefVertex( point_ptr );
   if( !ref_vertex_ptr )
   {
     PRINT_ERROR("Failed to construct free RefVertex.\n");
     assert(!!ref_vertex_ptr);
     return 0;
   }

   // Add the new ref_vertex to the graphics
   AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::FREE_REF_ENTITY_GENERATED, ref_vertex_ptr));
   
   CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_CREATED, ref_vertex_ptr);
   const_cast<CGMHistory&>(mHistory).add_event(evt);

   return ref_vertex_ptr;
}

void GeometryQueryTool::delete_Body(DLIList<Body*>& body_list)
{
  body_list.reset();
  for (int i = body_list.size(); i--; )
    delete_Body(body_list.get_and_step());
}

//-------------------------------------------------------------------------
// Purpose       : Delete a Body
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_Body( Body* body_ptr )
{
  BodySM* bodysm = body_ptr->get_body_sm_ptr();
  if (!bodysm)
  {
    PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n",body_ptr->id());
  }
  else
  {
      // Ask owning model engine to delete the TopologyBridges
    GeometryQueryEngine* gqe = bodysm->get_geometry_query_engine();
    gqe->delete_solid_model_entities(bodysm);

      // Clean up any virtual geometry that was on the body
    for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
      (*itor)->clean_out_deactivated_geometry();
  }

     // Check if body actually got deleted.
  return destroy_dead_entity( body_ptr );
}

//-------------------------------------------------------------------------
// Purpose       : Delete a Body
//
// Special Notes : Checks if child entities are merged and regenerates
//                 graphics for them.  This is necessary when 2 merged
//                 entities have been "force-merged", usually meaning they
//                 are not spatially equal.
//
// Creator       : Corey Ernst
//
// Creation Date : 08/25/04
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_single_Body( Body* body_ptr )
{
  BodySM* bodysm = body_ptr->get_body_sm_ptr();
  DLIList<RefEntity*> merged_children;
  DLIList<RefFace*> faces_to_reverse;
  DLIList<RefEdge*> edges_to_reverse;
  if (!bodysm)
  {
    PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n",body_ptr->id());
  }
  else
  {
    //get all the child entities that have been merged
    DLIList<RefEntity*> tmp_merged_children;
    MergeTool::instance()->contains_merged_children( body_ptr, tmp_merged_children );

    //get the owning bodies for each entity
    int i;
    for(i=tmp_merged_children.size(); i--;)
    {
      RefEntity *ref_ent = tmp_merged_children.get_and_step();
      TopologyEntity *tmp_entity = CAST_TO( ref_ent, TopologyEntity);
      DLIList<Body*> body_list;
      tmp_entity->bodies( body_list );
      //if 2 bodies own it, get body that is not "body_ptr"
      //for later graphics regeneration
      if( body_list.size() > 1 )
      {
        if( body_list.get() != body_ptr )
          merged_children.append( ref_ent );
        else if( body_list.step_and_get() != body_ptr )
          merged_children.append( ref_ent );
      }
    }

//    if( tmp_merged_children.size() )
//      MergeTool::instance()->unmerge( body_ptr );

    //fix up merged children -- some might need to be reversed
    for(i=merged_children.size(); i--; )
    {
      RefEntity *merged_child = merged_children.get_and_step();
      BasicTopologyEntity *bte = static_cast<BasicTopologyEntity*>(merged_child);

      //get the first bridge of the entity
      DLIList<TopologyBridge*> child_bridge_list;
      bte->bridge_manager()->get_bridge_list( child_bridge_list );
      child_bridge_list.reset();
      TopologyBridge *first_bridge = child_bridge_list.get_and_step();

      //if it is not the body we're deleting just continue
      BodySM *owning_body = first_bridge->bodysm();
      if( owning_body != bodysm )
        continue;

      RefFace *ref_face = CAST_TO( merged_child, RefFace );
      if( ref_face )
      {
        TopologyBridge *second_bridge = child_bridge_list.get_and_step();

        if( first_bridge->bridge_sense() != second_bridge->bridge_sense() )
          faces_to_reverse.append( ref_face );
        continue;
      }
      RefEdge *ref_edge = CAST_TO( merged_child, RefEdge );
      if( ref_edge )
      {
        //get merged_child's first topology bridge
        TopologyBridge *second_bridge = child_bridge_list.get_and_step();

        Curve *first_curve = CAST_TO( first_bridge, Curve );
        Curve *second_curve = CAST_TO( second_bridge, Curve );

        CubitSense relative_sense = first_curve->relative_sense( second_curve );

        if( relative_sense == CUBIT_REVERSED )
          edges_to_reverse.append( ref_edge );
        continue;
      }
    }

      // Ask owning model engine to delete the TopologyBridges
    GeometryQueryEngine* gqe = bodysm->get_geometry_query_engine();
    gqe->delete_solid_model_entities(bodysm);


      // Clean up any virtual geometry that was on the body
    for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
      (*itor)->clean_out_deactivated_geometry();
  }

  int i;
  for( i=faces_to_reverse.size(); i--; )
    faces_to_reverse.get_and_step()->reverse_normal();
  for( i=edges_to_reverse.size(); i--; )
    edges_to_reverse.get_and_step()->reverse_tangent();

    //regenerate graphics of merged entities
  for( i=merged_children.size(); i--; )
  {
    RefEntity* child = merged_children.get_and_step();
    AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::GEOMETRY_TOPOLOGY_MODIFIED, child));
    CGMHistory::Event evt2(CGMHistory::TOPOLOGY_CHANGED, child);
    const_cast<CGMHistory&>(mHistory).add_event(evt2);
    CGMHistory::Event evt(CGMHistory::GEOMETRY_CHANGED, child);
    const_cast<CGMHistory&>(mHistory).add_event(evt);
  }

  // Check if body actually got deleted.
  CubitStatus ret = destroy_dead_entity( body_ptr );

  // Send an unmerge event out for surfaces that were a part
  // of the body that was deleted.  Do this after the above
  // call to destroy_dead_entity() so that the Body
  // pointer associated with the deleted volume will no
  // longer be referenced by the merged surface.  This code
  // was added so that sidesets which are watching surface
  // modifications will know to update their data when 
  // the body is deleted. 
  for( i=merged_children.size(); i--; )
  {
    RefFace* child = CAST_TO(merged_children.get_and_step(), RefFace);
    if(child)
    {
      UnMergeEvent unmerge_event( child, child );
      AppUtil::instance()->send_event(unmerge_event );
    }
  }

  return ret;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a RefEntity
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_RefEntity( RefEntity* ref_entity_ptr )
{
  if (Body* body = dynamic_cast<Body*>(ref_entity_ptr))
    return delete_Body(body);

  if(RefFace* face = dynamic_cast<RefFace*>(ref_entity_ptr))
    return delete_RefFace(face);

  if(RefEdge* edge = dynamic_cast<RefEdge*>(ref_entity_ptr))
    return delete_RefEdge(edge);

  if(RefVertex* vtx = dynamic_cast<RefVertex*>(ref_entity_ptr))
    return delete_RefVertex(vtx);
  PRINT_ERROR("Cannot delete entity of type '%s'\n", ref_entity_ptr->class_name());
  return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a RefFace
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_RefFace( RefFace* ref_face )
{
    // Get the list of Surfaces owned by the RefFace
  BridgeManager* manager = ref_face->bridge_manager();
  DLIList<TopologyBridge*> bridge_list(manager->number_of_bridges());
  manager->get_bridge_list(bridge_list);

    // Ask each Surface's owning engine to destroy the surface.
  while( bridge_list.size() )
  {
    TopologyBridge* bridge = bridge_list.pop();
    Surface* surface = dynamic_cast<Surface*>(bridge);
    if (!surface)
    {
      PRINT_ERROR("RefFace %d is invalid -- attached TopologyBridge "
                  "is not a Surface.\n", ref_face->id());
      continue;
    }

    surface->get_geometry_query_engine()->
      delete_solid_model_entities( surface );
  }

    // Clean up any virtual geometry that was on the surface
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->clean_out_deactivated_geometry();

    // Check if all the Surfaces got deleted.
  return destroy_dead_entity( ref_face );
}

//-------------------------------------------------------------------------
// Purpose       : Delete a RefEdge
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_RefEdge( RefEdge* ref_edge )
{
    // Get the list of Curves owned by the RefEdge
  BridgeManager* manager = ref_edge->bridge_manager();
  DLIList<TopologyBridge*> bridge_list(manager->number_of_bridges());
  manager->get_bridge_list(bridge_list);

    // Ask each Curve's owning engine to destroy the surface.
  while( bridge_list.size() )
  {
    TopologyBridge* bridge = bridge_list.pop();
    Curve* curve = dynamic_cast<Curve*>(bridge);
    if (!curve)
    {
      PRINT_ERROR("RefEdge %d is invalid -- attached TopologyBridge "
                  "is not a Curve.\n", ref_edge->id());
      continue;
    }

    curve->get_geometry_query_engine()->
      delete_solid_model_entities( curve );
  }

    // Clean up any virtual geometry that was on the curve
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->clean_out_deactivated_geometry();

    // Check if all the curves got deleted.
  return destroy_dead_entity( ref_edge );
}

//-------------------------------------------------------------------------
// Purpose       : Delete a RefVertex
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::delete_RefVertex( RefVertex* ref_vertex )
{
    // Get the list of Points owned by the RefVertex
  BridgeManager* manager = ref_vertex->bridge_manager();
  DLIList<TopologyBridge*> bridge_list(manager->number_of_bridges());
  manager->get_bridge_list(bridge_list);

    // Ask each Curve's owning engine to destroy the surface.
  while( bridge_list.size() )
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBPoint* point = dynamic_cast<TBPoint*>(bridge);
    if (!point)
    {
      PRINT_ERROR("RefVertex %d is invalid -- attached TopologyBridge "
                  "is not a Point.\n", ref_vertex->id());
      continue;
    }

    point->get_geometry_query_engine()->
      delete_solid_model_entities( point );
  }

    // Clean up any virtual geometry that was on the Point
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->clean_out_deactivated_geometry();

    // Check if all the Points got deleted.
  return destroy_dead_entity( ref_vertex );
}

//Initialize all settings in this class
void GeometryQueryTool::initialize_settings()
{

  SettingHandler::instance()->add_setting("Merge Tolerance",
                                          GeometryQueryTool::set_geometry_factor,
					  GeometryQueryTool::get_geometry_factor);

  SettingHandler::instance()->add_setting("Merge Test BBox",
                                          GeometryQueryTool::set_merge_test_bbox,
					  GeometryQueryTool::get_merge_test_bbox);

  SettingHandler::instance()->add_setting("Merge Test InternalSurf",
					 GeometryQueryTool::set_merge_test_internal,
					 GeometryQueryTool::get_merge_test_internal);

  SettingHandler::instance()->add_setting("Facet BBox",
					 GeometryQueryTool::set_facet_bbox,
					 GeometryQueryTool::get_facet_bbox);
}

//-------------------------------------------------------------------------
// Purpose       : Constructor of the GeometryQueryTool class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryQueryTool::GeometryQueryTool(GeometryQueryEngine*gqe_ptr)
{
  if (gqe_ptr != NULL) add_gqe(gqe_ptr);
}

// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********


void GeometryQueryTool::geom_debug( DLIList<TopologyEntity*> topo_list )
{
     //This function was created March 1998 to assist in debugging
     //the unmerge command.  When sufficient time has passed and
     //it is no longer needed, this function can be deleted.
     //SR Jankovich
   for( int i = topo_list.size(); i > 0; i-- )
   {
      TopologyEntity *topo_entity = topo_list.get_and_step();
      if( !topo_entity )
          return;
      RefEntity *ref_entity = CAST_TO( topo_entity, RefEntity );
      if( !ref_entity )
          return;

      PRINT_INFO( "%s\n", ref_entity->entity_name().c_str() );
      DLIList<TopologyEntity*> next_topo_list;
      next_topo_list.clean_out();

      if( CAST_TO( ref_entity, Body ) )
      {

        DLIList<CoVolume*> co_vol_list;
        if( !topo_entity->co_volumes( co_vol_list ) )
           return;
        for( int j = co_vol_list.size(); j > 0; j-- )
        {
          CoVolume *co_vol = co_vol_list.get_and_step();
          PRINT_INFO( "   CoVolume %d (not id)\n",
                      co_vol_list.size() - j );
          TopologyEntity *temp_topo = CAST_TO( co_vol, TopologyEntity );
          if( !temp_topo )
             return;
          DLIList<RefVolume*> vol_list;
          if( !temp_topo->ref_volumes( vol_list ) )
             return;
          for( int k = vol_list.size(); k > 0; k-- )
          {
            RefVolume *vol = vol_list.get_and_step();
            PRINT_INFO( "      %s\n", vol->entity_name().c_str() );
            TopologyEntity *next_topo = CAST_TO( vol, TopologyEntity );
            next_topo_list.append( next_topo );
          }
        }
        GeometryQueryTool::geom_debug( next_topo_list );
      }
      else if( CAST_TO( ref_entity, RefVolume ) )
      {
        DLIList<Shell*> shell_list;
        if( !topo_entity->shells( shell_list ) )
           return;
        for( int m = shell_list.size(); m > 0; m-- )
        {
          Shell *shell = shell_list.get_and_step();
          PRINT_INFO( "   Shell %d (not id)\n", shell_list.size() - m );
          TopologyEntity *group_topo = CAST_TO( shell, TopologyEntity );
          if( !group_topo )
             return;
          DLIList<CoFace*> co_face_list;
          if( !group_topo->co_faces( co_face_list ) )
             return;
          for( int j = co_face_list.size(); j > 0; j-- )
          {
            CoFace *co_face = co_face_list.get_and_step();
            PRINT_INFO( "      CoFace %d (not id)\n",
                        co_face_list.size() - j );
            TopologyEntity *temp_topo = CAST_TO( co_face, TopologyEntity );
            if( !temp_topo )
               return;
            DLIList<RefFace*> face_list;
            if( !temp_topo->ref_faces( face_list ) )
               return;
            for( int k = face_list.size(); k > 0; k-- )
            {
              RefFace *face = face_list.get_and_step();
              PRINT_INFO( "         %s\n", face->entity_name().c_str() );
              TopologyEntity *next_topo = CAST_TO( face, TopologyEntity );
              next_topo_list.append( next_topo );
            }
          }
        }
        GeometryQueryTool::geom_debug( next_topo_list );
      }
      else if( CAST_TO( ref_entity, RefFace ) )
      {
        DLIList<Loop*> loop_list;
        if( !topo_entity->loops( loop_list ) )
           return;
        for( int m = loop_list.size(); m > 0; m-- )
        {
          Loop *loop = loop_list.get_and_step();
          PRINT_INFO( "   Loop %d (not id)\n", loop_list.size() - m );
          TopologyEntity *group_topo = CAST_TO( loop, TopologyEntity );
          if( !group_topo )
             return;
          DLIList<CoEdge*> co_edge_list;
          if( !group_topo->co_edges( co_edge_list ) )
             return;
          for( int j = co_edge_list.size(); j > 0; j-- )
          {
            CoEdge *co_edge = co_edge_list.get_and_step();
            PRINT_INFO( "      CoEdge %d (not id)\n",
                        co_edge_list.size() - j );
            TopologyEntity *temp_topo = CAST_TO( co_edge, TopologyEntity );
            if( !temp_topo )
               return;
            DLIList<RefEdge*> edge_list;
            if( !temp_topo->ref_edges( edge_list ) )
               return;
            for( int k = edge_list.size(); k > 0; k-- )
            {
              RefEdge *edge = edge_list.get_and_step();
              PRINT_INFO( "         %s\n", edge->entity_name().c_str() );
              TopologyEntity *next_topo = CAST_TO( edge, TopologyEntity );
              next_topo_list.append( next_topo );
            }
          }
        }
        GeometryQueryTool::geom_debug( next_topo_list );
      }
      else if( CAST_TO( ref_entity, RefEdge ) )
      {
        DLIList<Chain*> chain_list;
        if( !topo_entity->chains( chain_list ) )
           return;
        for( int m = chain_list.size(); m > 0; m-- )
        {
          Chain *chain = chain_list.get_and_step();
          PRINT_INFO( "   Chain %d (not id)\n", chain_list.size() - m );
          TopologyEntity *group_topo = CAST_TO( chain, TopologyEntity );
          if( !group_topo )
             return;
          DLIList<CoVertex*> co_vertex_list;
          if( !group_topo->co_vertices( co_vertex_list ) )
             return;
          for( int j = co_vertex_list.size(); j > 0; j-- )
          {
            CoVertex *co_vertex = co_vertex_list.get_and_step();
            PRINT_INFO( "      CoVertex %d (not id)\n",
                        co_vertex_list.size() - j );
            TopologyEntity *temp_topo = CAST_TO( co_vertex, TopologyEntity );
            if( !temp_topo )
               return;
            DLIList<RefVertex*> vertex_list;
            if( !temp_topo->ref_vertices( vertex_list ) )
               return;
            for( int k = vertex_list.size(); k > 0; k-- )
            {
              RefVertex *vertex = vertex_list.get_and_step();
              PRINT_INFO( "         %s\n",
                          vertex->entity_name().c_str() );
              TopologyEntity *next_topo = CAST_TO( vertex, TopologyEntity );
              next_topo_list.append( next_topo );
            }
          }
        }
        GeometryQueryTool::geom_debug( next_topo_list );
      }
      else if( CAST_TO( ref_entity, RefVertex ) )
         ;//Do nothing
      else
         PRINT_INFO( "UNKNOWN ENTITY TYPE!!!\n" );
   }
}

CubitBoolean GeometryQueryTool::about_spatially_equal (RefVertex* refVertex1,
                                                  RefVertex* refVertex2,
                                                  double tolerance_factor)
{
   return ( refVertex1->about_spatially_equal( refVertex2, tolerance_factor));
}

CubitBoolean GeometryQueryTool::about_spatially_equal (const CubitVector &vector1,
                                                  const CubitVector &vector2,
                                                  double tolerance_factor)
{
  double tol = GEOMETRY_RESABS * tolerance_factor;
  tol *= tol;
  double dist_sq = vector1.distance_between_squared(vector2);
  return !(dist_sq > tol);
  // within_tolerance() just checks coordinate aligned distances, not the actual distance.
//   return (vector1.within_tolerance( vector2, GEOMETRY_RESABS * tolerance_factor ));
}

double GeometryQueryTool::geometric_angle(RefEdge* ref_edge_1,
                                     RefEdge* ref_edge_2,
                                     RefFace* ref_face )
{
     // calculates internal surface angles given 2 refedges on the surface
   CoEdge *co_edge_1, *co_edge_2;
     //First get the two coedges that corrispond to the edges sent into
     //this function.  The two coedges are found from the ref_face's loop
     //where co_edge_1 is followed by co_edge_2 in the loop and co_edge_1 is
     //associated with ref_edge_1 while co_edge_2 is associated with
     //ref_edge_2.
   ref_edge_1->get_two_co_edges( ref_edge_2, ref_face, co_edge_1,
                                 co_edge_2 );

   return geometric_angle( co_edge_1, co_edge_2 );
}

double GeometryQueryTool::geometric_angle(CoEdge* co_edge_1,
                                     CoEdge* co_edge_2 )
{

  RefEdge *ref_edge_1 = co_edge_1->get_ref_edge_ptr();
  RefEdge *ref_edge_2 = co_edge_2->get_ref_edge_ptr();

    // return 2 pi for the tip of a hard line.
  if ( co_edge_1 != co_edge_2 && ref_edge_1 == ref_edge_2 )
    return 2.0 * CUBIT_PI;

  RefVertex *ref_vertex =
    co_edge_1->get_sense() == CUBIT_FORWARD ?
    ref_edge_1->end_vertex() :
    ref_edge_1->start_vertex();

  //RefVertex *junk_vertex1 = ref_edge_2->start_vertex();
  //RefVertex *junk_vertex2 = ref_edge_2->end_vertex();
  //CubitSense junk_sense = co_edge_2->get_sense();

  assert( ref_vertex ==
          ( co_edge_2->get_sense() == CUBIT_FORWARD ?
            ref_edge_2->start_vertex() :
            ref_edge_2->end_vertex() ) );

     // coordinates of common point
  CubitVector vertex_point = ref_vertex->coordinates();

    // Find normal to the the face at the common vertex of
    // the refedges. Use loop sense to artificially switch inner loops.
    // Use intrinsic normal.
  RefFace *ref_face = co_edge_1->get_ref_face();
  CubitVector normal = ref_face->normal_at(vertex_point, NULL);

    // Find directed tangents to determine interior angle
    // Use sense of edge with respect to this face's loop.
  CubitVector tangent_1, tangent_2;
  ref_edge_1->tangent( vertex_point, tangent_1 );
  ref_edge_2->tangent( vertex_point, tangent_2 );

  if ( co_edge_1->get_sense() == CUBIT_REVERSED )
    tangent_1 = -tangent_1;
  if ( co_edge_2->get_sense() == CUBIT_REVERSED )
    tangent_2 = -tangent_2;

    //  At this point we have the tangents going in the correct loop
    //  sense.
    //  Now get tangent pointing away from the center for the correct
    //  angle
  tangent_1 = -tangent_1;
    // Return angle from given tangents and normal to face
  double angle = normal.vector_angle( tangent_2, tangent_1 );
  if ( angle*180.0/CUBIT_PI >  360.0 - GEOMETRY_RESABS )
  {
      //try other points to make sure this is going the right way.
    CubitVector new_loc_1, new_loc_2;
    if ( ref_edge_1->start_vertex() == ref_vertex )
      ref_edge_1->position_from_fraction(0.01,
                                    new_loc_1);
    else
      ref_edge_1->position_from_fraction(0.99,
                                    new_loc_1);
    if ( ref_edge_2->start_vertex() == ref_vertex )
      ref_edge_2->position_from_fraction(0.01,
                                    new_loc_2);
    else
      ref_edge_2->position_from_fraction(0.99,
                                    new_loc_2);
      //Now just do the exact same thing as above but
      //use these new points...
    ref_edge_1->tangent( new_loc_1, tangent_1 );
    ref_edge_2->tangent( new_loc_2, tangent_2 );

    if ( co_edge_1->get_sense() == CUBIT_REVERSED )
      tangent_1 = -tangent_1;
    if ( co_edge_2->get_sense() == CUBIT_REVERSED )
      tangent_2 = -tangent_2;
    tangent_1 = -tangent_1;
      // Return angle from given tangents and normal to face
    angle = normal.vector_angle( tangent_2, tangent_1 );
    if ( angle < CUBIT_PI )
      angle = 0.0;
    else
      angle = 2.0*CUBIT_PI;
  }
  return angle;
}

CubitString GeometryQueryTool::get_engine_version_string()
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->get_engine_version_string();
  }
  else
  {
    return CubitString("No Active GeometryEngine");
  }
}

CubitBoolean
GeometryQueryTool::does_geom_contain_query_engine(DLIList<TopologyEntity*> &topo_list,
                                                  GeometryQueryEngine *engine) const
{
  GeometryQueryEngine *ge_ptr;
  for( int i=topo_list.size(); i--; )
  {
    ge_ptr = topo_list.get_and_step()->get_geometry_query_engine();
    if( ge_ptr == engine )
    {
      return CUBIT_TRUE;
    }
  }

  return CUBIT_FALSE;
}

CubitBoolean
GeometryQueryTool::does_geom_contain_query_engine(DLIList<RefEntity*> &ref_entity_list,
                                                  GeometryQueryEngine *engine,
                                                  CubitBoolean children_too) const
{
  DLIList<RefEntity*> complete_entity_list;

  // Check the check_children option and check all the children if necessary
  if (children_too)
  {
    //Make a complete list of all the RefEntitys and their children
    DLIList<RefEntity*> temp = ref_entity_list;
    RefEntity* ref_entity_ptr;
    int i;
    for( i=ref_entity_list.size(); i--; )
    {
      ref_entity_ptr = ref_entity_list.get_and_step();
      complete_entity_list.clean_out();
      ref_entity_ptr->get_all_child_ref_entities(complete_entity_list);
      temp += complete_entity_list;
    }
    complete_entity_list.clean_out();
    complete_entity_list.merge_unique(temp);
  }

  // Now check the RefEntities for the given geometry engine
  DLIList<TopologyEntity*> te_list;
  CAST_LIST(complete_entity_list, te_list, TopologyEntity);
  return does_geom_contain_query_engine(te_list, engine);
}

TopologyEntity* GeometryQueryTool::entity_from_bridge( TopologyBridge* bridge_ptr ) const
{
  if( !bridge_ptr )
    return NULL;
  TBOwner* owner = bridge_ptr->owner();
  BridgeManager* bridge_manager;
  while (!(bridge_manager = dynamic_cast<BridgeManager*>(owner)))
  {
    if (TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(owner))
      owner = bridge->owner();
    else if(TBOwnerSet* set = dynamic_cast<TBOwnerSet*>(owner))
    {
      DLIList<TopologyBridge*> list;
      set->get_owners(list);
      list.reset();
      owner = list.get()->owner();
    }
    else
      break;
  }

  return bridge_manager ? bridge_manager->topology_entity() : 0;
}


CubitStatus GeometryQueryTool::set_default_gqe(GeometryQueryEngine* gqe)
{
  int i;
  for (i = 0; i < gqeList.size(); i++)
  {
    if(gqe == gqeList.get())
      break;

    gqeList.step();
  }

  if(i == gqeList.size()) return CUBIT_FAILURE;

  GeometryQueryEngine* temp_ptr = gqeList.get();
  gqeList.remove();
  gqeList.insert_first(temp_ptr);
  PRINT_INFO("Geometry engine set to: %s\n", gqe->get_engine_version_string().c_str() );
  return CUBIT_SUCCESS;
}

CubitStatus GeometryQueryTool::set_export_allint_version(int version)
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->set_export_allint_version(version);
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return CUBIT_FAILURE;
  }
}

int GeometryQueryTool::get_allint_version()
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->get_allint_version();
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return 0;
  }
}

CubitStatus GeometryQueryTool::list_engine_versions(CubitString &versions)
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->list_engine_versions(versions);
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return CUBIT_FAILURE;
  }
}

double GeometryQueryTool::get_sme_resabs_tolerance()
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->get_sme_resabs_tolerance();
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return 0.0;
  }
}

double GeometryQueryTool::set_sme_resabs_tolerance( double new_resabs )
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->set_sme_resabs_tolerance( new_resabs );
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return 0.0;
  }
}

CubitStatus GeometryQueryTool::set_sme_int_option( const char* opt_name, int val )
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->set_int_option( opt_name, val );
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return CUBIT_FAILURE;
  }
}

CubitStatus GeometryQueryTool::set_sme_dbl_option( const char* opt_name, double val )
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->set_dbl_option( opt_name, val );
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return CUBIT_FAILURE;
  }
}

CubitStatus GeometryQueryTool::set_sme_str_option( const char* opt_name, const char* val )
{
  if (gqeList.size())
  {
    gqeList.reset();
    return gqeList.get()->set_str_option( opt_name, val );
  }
  else
  {
    PRINT_WARNING("No active geometry engine.");
    return CUBIT_FAILURE;
  }
}

CubitStatus GeometryQueryTool::get_intersections( RefEdge* ref_edge1,
                                                  CubitVector& point1,
                                                  CubitVector& point2,
                                                  DLIList<CubitVector>& intersection_list,
                                                  CubitBoolean bounded,
                                                  CubitBoolean closest)
{
  Curve* curve_ptr1 = ref_edge1->get_curve_ptr();
  if( curve_ptr1 == NULL )
  {
    if( curve_ptr1 == NULL )
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Curve %d\n"      "       This is a bug - please report it\n", ref_edge1->id() );
    return CUBIT_FAILURE;
  }

  if ( curve_ptr1->geometry_type() == STRAIGHT_CURVE_TYPE )
  {
    CubitVector dir = point2 - point1;
    return straightline_intersections(ref_edge1, point1, dir,
                             intersection_list, bounded, closest);
  }

  GeometryQueryEngine* GQE_ptr =
      curve_ptr1->get_geometry_query_engine();
  return GQE_ptr->get_intersections( curve_ptr1, point1, point2,
                                     intersection_list, bounded, closest);
}

CubitStatus GeometryQueryTool::straightline_intersections(RefEdge* ref_edge1,
                                 CubitVector & origin2,
                                 CubitVector & dir,
                                 DLIList<CubitVector> &intersection_list,
                                 CubitBoolean bounded ,
                                 CubitBoolean closest)
{
  Curve* curve_ptr1 = ref_edge1->get_curve_ptr();

  if( curve_ptr1 == NULL )
  {
    if( curve_ptr1 == NULL )
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Curve %d\n"      "       This is a bug - please report it\n", ref_edge1->id() );
    return CUBIT_FAILURE;
  }

  CubitVector dir2 = dir;
  dir2.normalize();
  assert( curve_ptr1->geometry_type() == STRAIGHT_CURVE_TYPE );
  // Proceed with the intersection calculation
  CubitVector origin1, dir1;
  double origin_pnt1[3], dir_vec1[3], origin_pnt2[3], dir_vec2[3];

  if( ref_edge1->get_point_direction( origin1, dir1 ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to get straight line information for Curve %d; aborting\n",
      ref_edge1->id() );
    return CUBIT_FAILURE;
  }

  origin1.get_xyz( origin_pnt1 ); origin2.get_xyz( origin_pnt2 );
  dir1.get_xyz( dir_vec1 ); dir2.get_xyz( dir_vec2 );

  AnalyticGeometryTool* agt = AnalyticGeometryTool::instance();

  int num_int;
  double int_pnt1[3], int_pnt2[3];
  num_int = agt->int_ln_ln( origin_pnt1, dir_vec1, origin_pnt2, dir_vec2,
    int_pnt1, int_pnt2 );
  if( num_int == 0 || (closest == CUBIT_FALSE && num_int == 2) )
  {
    if( num_int == 0 )
      PRINT_ERROR( "Curves %d and the straight line defined by position %g, %g, %g, and direction %g, %g, %g are parallel - no intersection exists\n",
 ref_edge1->id(),  origin2.x(), origin2.y(), origin2.z(), dir2.x(), dir2.y(), dir2.z());

    else
      PRINT_ERROR( "Curves %d and the straight line defined by position %g, %g, %g, and direction %g, %g, %g do not intersect\n", ref_edge1->id(),
      origin2.x(), origin2.y(), origin2.z(), dir2.x(), dir2.y(), dir2.z());
    return CUBIT_FAILURE;
  }

  if( bounded == CUBIT_TRUE )
  {
    CubitVector start1 = ref_edge1->start_vertex()->coordinates();
    CubitVector end1 = ref_edge1->end_vertex()->coordinates();
    CubitVector start2 = origin2;
    CubitVector end2 = origin2 + dir;

    double start_pnt1[3], end_pnt1[3], start_pnt2[3], end_pnt2[3];
    start1.get_xyz( start_pnt1 ); start2.get_xyz( start_pnt2 );
    end1.get_xyz( end_pnt1 ); end2.get_xyz( end_pnt2 );

    if( num_int == 1 )
    {
      // Vertex must be on both curves
      if( agt->is_pnt_on_ln_seg( int_pnt1, start_pnt1, end_pnt1 ) &&
        agt->is_pnt_on_ln_seg( int_pnt2, start_pnt2, end_pnt2 ) )
      {
        intersection_list.append( CubitVector(int_pnt1) );
      }
      else
      {
        PRINT_WARNING( "intersection point of Curves was not within bounds of both curves\n");
      }
    }
    else
    {
      // Only keep the vertices that are on the curve bounds
      if( agt->is_pnt_on_ln_seg( int_pnt1, start_pnt1, end_pnt1 ) )
      {
        intersection_list.append( CubitVector(int_pnt1) );
      }
      else
      {
        PRINT_WARNING( "intersection point on Curve %d was not within it's bounds\n",
          ref_edge1->id() );
      }
      if( agt->is_pnt_on_ln_seg( int_pnt2, start_pnt2, end_pnt2 ) )
      {
        intersection_list.append( CubitVector(int_pnt2) );
      }
      else
      {
        PRINT_WARNING( "intersection point on the straight line defined by position %g, %g, %g, and direction %g, %g, %g  was not within it's bounds\n",
        origin2.x(), origin2.y(), origin2.z(), dir2.x(), dir2.y(), dir2.z() );
      }
    }
    if( intersection_list.size() == 0 )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }
  else // Not bounded
  {
    intersection_list.append( CubitVector(int_pnt1) );
    if( num_int == 2 )
    {
      intersection_list.append( CubitVector(int_pnt2) );
    }
    return CUBIT_SUCCESS;
  }
}

CubitStatus GeometryQueryTool::get_intersections( RefEdge* ref_edge1,
                                                  RefEdge* ref_edge2,
                                                  DLIList<CubitVector>& intersection_list,
                                                  CubitBoolean bounded,
                                                  CubitBoolean closest)
{
  // If both curves are straight, compute their intersection; otherwise
  // use the geometry engine to do it.
  CubitStatus status = CUBIT_FAILURE;
  Curve* curve_ptr1 = ref_edge1->get_curve_ptr();
  Curve* curve_ptr2 = ref_edge2->get_curve_ptr();

  if( curve_ptr1 == NULL || curve_ptr2 == NULL )
  {
    if( curve_ptr1 == NULL )
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Curve %d\n"
      "       This is a bug - please report it\n", ref_edge1->id() );
    if( curve_ptr2 == NULL )
      PRINT_ERROR("Unable to retrieve underlying geometry entity of Curve %d\n"
      "       This is a bug - please report it\n", ref_edge2->id() );
    return CUBIT_FAILURE;
  }

  if( curve_ptr1->geometry_type() == STRAIGHT_CURVE_TYPE &&
    curve_ptr2->geometry_type() == STRAIGHT_CURVE_TYPE )
  {
    // Proceed with the intersection calculation
    CubitVector origin2, dir2;

    if( ref_edge2->get_point_direction( origin2, dir2 ) == CUBIT_FAILURE )
    {
      PRINT_ERROR( "Unable to get straight line information for Curve %d; aborting\n",
        ref_edge2->id() );
      return CUBIT_FAILURE;
    }

    if (bounded)
    {
       origin2 = ref_edge2->start_vertex()->coordinates();
       dir2 = ref_edge2->end_vertex()->coordinates() - origin2;
    }
    return straightline_intersections(ref_edge1, origin2, dir2,
                             intersection_list, bounded, closest);
  }

  if( closest == CUBIT_TRUE )
  {
    PRINT_ERROR( "'Near' option only works for straight lines\n" );
    return CUBIT_FAILURE;
  }

  // Use geometry engine to find intersections
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  DLIList<TopologyBridge*> curve_list1, curve_list2;
  entity_list.append(ref_edge1);
  entity_list.append(ref_edge2);

  //if there's virtual edge, find the underlying real curve.
  if (ref_edge1->get_geometry_query_engine())
    ref_edge1->get_geometry_query_engine()->get_underlying_curves(curve_ptr1, curve_list1);
  if (ref_edge2->get_geometry_query_engine())
    ref_edge2->get_geometry_query_engine()->get_underlying_curves(curve_ptr2, curve_list2);

  GeometryQueryEngine* gqe = common_query_engine( entity_list, bridge_list );

  //if they are both virtual...set gqe to NULL
  if( is_intermediate_geometry(curve_ptr1) && is_intermediate_geometry( curve_ptr2 ) )
    gqe = NULL;

  if( gqe == NULL && (curve_list1.size() > 0 || curve_list2.size() > 0) )
  {
    if (curve_list1.size() == 0)
      curve_list1.append (CAST_TO(curve_ptr1, TopologyBridge));
    if (curve_list2.size() == 0)
      curve_list2.append (CAST_TO(curve_ptr2, TopologyBridge));
    int i, j;
    for (i = 0; i <curve_list1.size(); i++)
    {
      TopologyBridge * tb1 = curve_list1.get_and_step();
      GeometryQueryEngine *gqe_ptr1 = tb1->get_geometry_query_engine();
      curve_ptr1 = CAST_TO(tb1, Curve);

      for(j = 0; j< curve_list2.size(); j++)
      {      
        TopologyBridge * tb2 = curve_list2.get_and_step();      
        GeometryQueryEngine *gqe_ptr2 = tb2->get_geometry_query_engine();
        if (gqe_ptr1 && gqe_ptr1 == gqe_ptr2 )
        {      
          curve_ptr2 = CAST_TO(tb2, Curve);
          status = gqe_ptr1->get_intersections(curve_ptr1, curve_ptr2,
            intersection_list, bounded, closest );
        }
      }
    }

    //remove duplicate intersections
    for( int k=0; k<intersection_list.size(); k++ )
    {
      intersection_list.reset();
      intersection_list.step(k+1);
      for( int s=k+1; s<intersection_list.size();)
      {
        if( intersection_list[k].distance_between( intersection_list[s] ) < GEOMETRY_RESABS )
        {
          intersection_list.remove();
        }
        else
        {
          s++;
          intersection_list.step();
        }
      }
    }

    return status;
  }

  else if (gqe == NULL)
  {
    PRINT_ERROR( "Curves %d and %d do not have the same underlying geometry modeling engine\n"
      "       For intersection calculations, they must be the same\n",
      ref_edge1->id(), ref_edge2->id() );
    return CUBIT_FAILURE;
  }

  bridge_list.reset();
  curve_ptr1 = dynamic_cast<Curve*>(bridge_list.next(0));
  curve_ptr2 = dynamic_cast<Curve*>(bridge_list.next(1));
  return gqe->get_intersections( curve_ptr1, curve_ptr2,
                                 intersection_list, bounded, closest );
}

CubitStatus
GeometryQueryTool::get_intersections( RefEdge* ref_edge, RefFace* ref_face,
                                      DLIList<CubitVector>& intersection_list,
                                      CubitBoolean bounded )
{
  // Use geometry engine to find intersections
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  entity_list.append(ref_edge);
  entity_list.append(ref_face);
  GeometryQueryEngine* gqe = common_query_engine( entity_list, bridge_list );

  if( gqe == NULL )
  {
    PRINT_ERROR( "Curve %d and Surface %d do not have the same underlying geometry query engine\n"
      "       For intersection calculations, they must be the same\n",
      ref_edge->id(), ref_face->id() );
    return CUBIT_FAILURE;
  }

  bridge_list.reset();
  Curve* curve = dynamic_cast<Curve*>(bridge_list.next(0));
  Surface* surf = dynamic_cast<Surface*>(bridge_list.next(1));
  return gqe->get_intersections( curve, surf, intersection_list, bounded );
}

CubitStatus
GeometryQueryTool::get_intersections(RefEdge* ref_edge, CubitPlane plane,
                                      DLIList<CubitVector> &intersection_list,
                                      CubitBoolean bounded, double extended_percent )
{
  CubitBox box = ref_edge->bounding_box();

  // create a Surface from the plane
  CubitVector p1, p2, p3, p4;
  AnalyticGeometryTool::instance()->min_pln_box_int_corners(plane, box, 1, extended_percent, p1, p2, p3, p4, true );

  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme = GeometryModifyTool::instance()->get_engine( ref_edge, &bridge );
  BodySM* sheet = gme->planar_sheet(p1, p2, p3, p4);
  if (!sheet)
  {
    PRINT_INFO("%s", "Unable to test for planar intersections\n");
    return CUBIT_FAILURE;
  }

  DLIList<Surface*> surfaces;
  sheet->surfaces( surfaces );

  Curve* curve = ref_edge->get_curve_ptr();
  GeometryQueryEngine* gqe = ref_edge->get_geometry_query_engine();
  CubitStatus status = gqe->get_intersections( curve, surfaces[0], intersection_list, bounded );

  gqe->delete_solid_model_entities(sheet);

  return status;
}

//===============================================================================
// Function   : entity_extrema
// Member Type: PUBLIC
// Description: Find extrema location on entity
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus
GeometryQueryTool::entity_extrema( RefEntity *ref_entity_ptr,
                                   const CubitVector *dir1,
                                   const CubitVector *dir2,
                                   const CubitVector *dir3,
                                   CubitVector &extrema,
                                   RefEntity *&extrema_entity_ptr )
{
  DLIList<RefEntity*> ref_entity_list;
  ref_entity_list.append( ref_entity_ptr );

  return entity_extrema( ref_entity_list, dir1, dir2, dir3, extrema,
                         extrema_entity_ptr );
}

//===============================================================================
// Function   : entity_extrema
// Member Type: PUBLIC
// Description: Find extrema location on a list of entities
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus
GeometryQueryTool::entity_extrema( DLIList<RefEntity*> &ref_entity_list,
                                   const CubitVector *dir1,
                                   const CubitVector *dir2,
                                   const CubitVector *dir3,
                                   CubitVector &extrema,
                                   RefEntity *&extrema_entity_ptr )
{
  if( ref_entity_list.size() == 0 )
  {
    PRINT_ERROR( "No entities found for extrema calculation.\n" );
    return CUBIT_FAILURE;
  }

  DLIList<TopologyEntity*> entity_list(ref_entity_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_entity_list.size());
  DLIList<RefVolume*> ref_vols;
    // Can only do BasicTopologyEntitys.  Relace Bodys with RefVolumes.
  ref_entity_list.reset();
  for (int i = ref_entity_list.size(); i--; )
  {
    RefEntity* entity = ref_entity_list.get_and_step();
    if (BasicTopologyEntity* bte = dynamic_cast<BasicTopologyEntity*>(entity))
    {
      entity_list.append(bte);
      continue;
    }

    if (Body* body = dynamic_cast<Body*>(entity))
    {
      ref_vols.clean_out();
      body->ref_volumes(ref_vols);
      ref_vols.reset();
      for (int j = ref_vols.size(); j--;)
        entity_list.append(ref_vols.get_and_step());
      continue;
    }

    PRINT_ERROR("Don't know how to handle entity of type '%s' in "
                "GQT::entity_extrema.\n", entity->class_name());
    return CUBIT_FAILURE;
  }


  GeometryQueryEngine* gqe = common_query_engine(entity_list, bridge_list);
  if( gqe == NULL )
  {
    PRINT_ERROR( "Entities specified for extrema do not have the same\n"
      "       underlying geometry query engine.\n"
      "       For extrema calculations, they must be the same.\n" );
    return CUBIT_FAILURE;
  }

    // Need GeometryEntities.
  DLIList<GeometryEntity*> geom_list(bridge_list.size());
  CAST_LIST(bridge_list, geom_list, GeometryEntity);

  // Check direction inputs
  if( dir1 == NULL )
  {
    PRINT_ERROR( "Direction not found for extrema calculation - it is required.\n" );
    return CUBIT_FAILURE;
  }
  if( dir2 == NULL && dir3 )
  {
    PRINT_ERROR( "Second direction not specified but last direction\n"
    "       was - this is not allowed for extrema calculation.\n" );
    return CUBIT_FAILURE;
  }

  GeometryEntity* extrema_geom = 0;
  CubitStatus status = gqe->entity_extrema( geom_list, dir1, dir2, dir3,
                                           extrema, extrema_geom );
  extrema_entity_ptr = dynamic_cast<RefEntity*>(entity_from_bridge(extrema_geom));
  return status;
}

CubitStatus
GeometryQueryTool::entity_entity_distance( RefEntity *ref_entity_ptr1,
                                           RefEntity *ref_entity_ptr2,
                                           CubitVector &pos1, CubitVector &pos2,
                                           double &distance )
{
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);

  BasicTopologyEntity* bte1 = dynamic_cast<BasicTopologyEntity*>(ref_entity_ptr1);
  BasicTopologyEntity* bte2 = dynamic_cast<BasicTopologyEntity*>(ref_entity_ptr2);
  if (!bte1 || !bte2)
  {
    const char* name = bte2 ? ref_entity_ptr1->class_name() : ref_entity_ptr2->class_name();
    PRINT_ERROR("Cannot calculate entity distance for entity of type '%s'\n", name);
    return CUBIT_FAILURE;
  }

  entity_list.append(bte1);
  entity_list.append(bte2);
  GeometryQueryEngine* gqe = common_query_engine(entity_list, bridge_list);
  if( gqe == NULL )
  {
    PRINT_ERROR( "%s and %s do not have the same underlying geometry query engine.\n"
      "       For distance calculations, they must be the same.\n",
      ref_entity_ptr1->entity_name().c_str(), ref_entity_ptr2->entity_name().c_str() );
    return CUBIT_FAILURE;
  }

  bridge_list.reset();
  GeometryEntity* geom1 = dynamic_cast<GeometryEntity*>(bridge_list.next(0));
  GeometryEntity* geom2 = dynamic_cast<GeometryEntity*>(bridge_list.next(1));
  return gqe->entity_entity_distance( geom1, geom2, pos1, pos2, distance );
}

CubitStatus
GeometryQueryTool::entity_entity_distance( GeometryEntity *ge1,
                                           GeometryEntity *ge2,
                                           CubitVector &pos1, CubitVector &pos2,
                                           double &distance )
{
  GeometryQueryEngine *gqe1 = ge1->get_geometry_query_engine();
  GeometryQueryEngine *gqe2 = ge2->get_geometry_query_engine();

  if(gqe1 != gqe2)
  {
    if(gqe1->is_intermediate_engine())
      return gqe1->entity_entity_distance(ge1, ge2, pos1, pos2, distance);
    else if(gqe2->is_intermediate_engine())
      return gqe2->entity_entity_distance(ge1, ge2, pos1, pos2, distance);
    else
    {
      PRINT_ERROR( "Entities do not have the same underlying geometry query engine.\n"
        "       For distance calculations, they must be the same.\n" );
      return CUBIT_FAILURE;
    }
  }
  else
  {
    return gqe1->entity_entity_distance(ge1, ge2, pos1, pos2, distance);
  }
}

CubitBox GeometryQueryTool::bounding_box_of_bodies()
{
   CubitBox total_bound;
   int i;

   if(GeometryQueryTool::instance()->num_bodies() > 0)
   {
     GeometryQueryTool::instance()->get_last_body();
        // Assign the first body's bounding box to total_bound
     total_bound = GeometryQueryTool::instance()->get_last_body()->bounding_box();

        // Create union of all remaining bounding boxes
      for (i = GeometryQueryTool::instance()->num_bodies(); --i; )
          total_bound |= GeometryQueryTool::instance()->get_next_body()->bounding_box();
   }
   else
   {
        // Set the box to [(0,0,0)(0,0,0)]
      CubitVector temp(0,0,0);
      total_bound.reset(temp);
   }
     // Return the box
   return total_bound;
}

void GeometryQueryTool::cleanout_deactivated_geometry ()
{


     // This routine removes from the model all RefEntities that have been
     // deactivated
   PRINT_DEBUG_17("\n\n...Cleaning out deactivated RefEntities "
               "and asociated Mesh Entities from the model\n");

     // First delete the meshes associated with the deactivated RefEntities
     // delete_meshes_of_deactivated_refEntities();

     // Now delete the deactivated RefEntities and remove all traces
     // of them from the DAG
   DAG::instance()->cleanout_deactivated_DAG_nodes();

     // Leave gracefully :-)
   PRINT_DEBUG_17(
               "Successfully cleaned out all deactivated RefEntities "
               "and associated Mesh Entities from the model.\n");

   return;
}

void GeometryQueryTool::delete_geometry ()
{

   CubitStatus status = CUBIT_FAILURE;

     // Delete all RefGroups
   RefGroup::delete_all_groups();
   Body* bodyPtr = NULL;
   //DLIList<Body*> bodyList;
   //GeometryQueryTool::instance()->bodies(bodyList);
   //bodyList.reset();
   //for ( i = bodyList.size(); i > 0; i--)
   //{
   //   bodyPtr = bodyList.get_and_step();

   while( (bodyPtr = GeometryQueryTool::instance()->get_last_body() ) != NULL )
	 {
        // Now delete this Body and its underlying solid model entities
      status = GeometryQueryTool::instance()->delete_Body( bodyPtr );
      if (status == CUBIT_FAILURE)
      {
         PRINT_ERROR("In GeometryQueryTool::delete_geometry\n"
                     "       Could not delete Body %d.\n"
                     "       The Model database is likely corrupted due"
                     " to\n       this unsuccessful deletion.\n",
                     bodyPtr->id() );
         assert( status != CUBIT_FAILURE) ;
         break;
      }
   }

     // Remove free-floating RefFaces
   RefFace* refFacePtr = NULL;
   //DLIList<RefFace*> refFaceList;
   //GeometryQueryTool::instance()->ref_faces(refFaceList);
   //refFaceList.reset();
   //for ( i = refFaceList.size(); i > 0; i--)
   //{
   //   refFacePtr = refFaceList.get_and_step();

   while( (refFacePtr = GeometryQueryTool::instance()->get_last_ref_face() ) != NULL )
	 {
        // NOTE-
        //   The following GeometryQueryTool call results in a call to
        //   the RefFace destructor which notifies the Model of
        //   the destruction event. Model, in turn, removes the RefFace
        //   pointer from refFaceList. Hence, the size of this list
        //   changes inside this "while" loop.
      status = GeometryQueryTool::instance()->delete_RefFace( refFacePtr );
      if (status == CUBIT_FAILURE)
      {
         PRINT_ERROR("In GeometryQueryTool::delete_geometry\n"
                     "       Could not delete RefFace %d.\n"
                     "       The Model database is likely corrupted "
                     "due to\n       this unsuccessful deletion.\n",
                     refFacePtr->id());
         assert( status != CUBIT_FAILURE) ;
         break;
      }
   }

     // Remove free-floating RefEdges
   RefEdge* refEdgePtr = NULL;
   //DLIList<RefEdge*> refEdgeList;
   //GeometryQueryTool::instance()->ref_edges(refEdgeList);
   //refEdgeList.reset();
   //for ( i = refEdgeList.size(); i > 0; i--)
   //{
	 //   refEdgePtr = refEdgeList.get_and_step();

   while( (refEdgePtr = GeometryQueryTool::instance()->get_last_ref_edge() ) != NULL )
	 {
        // NOTE-
        //   The following GeometryQueryTool call results in a call to
        //   the RefEdge destructor which notifies the Model of
        //   the destruction event. Model, in turn, removes the RefEdge
        //   pointer from refEdgeList. Hence, the size of this list
        //   changes inside this "while" loop.
      status = GeometryQueryTool::instance()->delete_RefEdge( refEdgePtr );
      if (status == CUBIT_FAILURE)
      {
         PRINT_ERROR("In GeometryQueryTool::delete_geometry\n"
                     "       Could not delete RefEdge %d.\n"
                     "       The Model database is likely corrupted "
                     "due to\n       this unsuccessful deletion.\n",
                     refEdgePtr->id());
         assert( status != CUBIT_FAILURE) ;
         break;
      }
   }

     // Remove free-floating RefVertex'es
   RefVertex* refVertexPtr = NULL;
   //DLIList<RefVertex*> refVertexList;
   //GeometryQueryTool::instance()->ref_vertices(refVertexList);
   //refVertexList.reset();
   //for ( i = refVertexList.size(); i > 0; i--)
   //{
   //   refVertexPtr = refVertexList.get_and_step();

   while( (refVertexPtr = GeometryQueryTool::instance()->get_last_ref_vertex() ) != NULL )
	 {
        // NOTE-
        //   The following GeometryQueryTool call results in a call to
        //   the RefVertex destructor which notifies the Model of
        //   the destruction event. Model, in turn, removes the RefVertex
        //   pointer from refVertexList. Hence, the size of this list
        //   changes inside this "while" loop.
      status = GeometryQueryTool::instance()->delete_RefVertex( refVertexPtr );
      if (status == CUBIT_FAILURE)
      {
         PRINT_ERROR("In GeometryQueryTool::delete_geometry\n"
                     "       Could not delete RefVertex %d.\n"
                     "       The Model database is likely corrupted "
                     "due to\n       this unsuccessful deletion.\n",
                     refVertexPtr->id());
         assert( status != CUBIT_FAILURE) ;
         break;
      }
   }


     // reset counters
   RefEntityFactory::instance()->reset_ids();

}

CubitStatus GeometryQueryTool::ref_entity_list(char const* keyword,
                                          DLIList<RefEntity*> &entity_list,
                                          const CubitBoolean print_errors)
{
  return RefEntityFactory::instance()->ref_entity_list(keyword, entity_list,
                                                       print_errors);
}

CubitBox GeometryQueryTool::model_bounding_box()
{
  int i;
  CubitBox result;
  DLIList<Body*> bodies;
  this->bodies(bodies);

  DLIList<RefEntity*> free_entity_list;
  this->get_free_ref_entities(free_entity_list);
  
  if (!(bodies.size() + free_entity_list.size()))
    return result;

  bodies.reset();
  free_entity_list.reset();
  
  if( bodies.size() )
  {
    result = bodies.get_and_step()->bounding_box();
    for ( i = bodies.size()-1; i--; )
      result |= bodies.get_and_step()->bounding_box();
    i = free_entity_list.size();
  }
  else
  {
    result = free_entity_list.get_and_step()->bounding_box();
    i = free_entity_list.size()-1;
  }

  for ( ; i>0 ; i-- )
    result |= free_entity_list.get_and_step()->bounding_box();

  return result;
}

void GeometryQueryTool::bodies(DLIList<Body*> &bodies)
{
  RefEntityFactory::instance()->bodies(bodies);
}

void GeometryQueryTool::ref_volumes(DLIList<RefVolume*> &ref_volumes)
{
  RefEntityFactory::instance()->ref_volumes(ref_volumes);
}

void GeometryQueryTool::ref_groups(DLIList<RefGroup*> &ref_groups)
{
  RefEntityFactory::instance()->ref_groups(ref_groups);
}

void GeometryQueryTool::ref_faces(DLIList<RefFace*> &ref_faces)
{
  RefEntityFactory::instance()->ref_faces(ref_faces);
}

void GeometryQueryTool::ref_edges(DLIList<RefEdge*> &ref_edges)
{
  RefEntityFactory::instance()->ref_edges(ref_edges);
}

void GeometryQueryTool::ref_vertices(DLIList<RefVertex*> &ref_vertices)
{
  RefEntityFactory::instance()->ref_vertices(ref_vertices);
}

void GeometryQueryTool::get_ordered_loops(RefFace* face, DLIList<Loop*> &loop_list)
{
  GeometryQueryEngine *gqe = face->get_geometry_query_engine();
  gqe->get_ordered_loops(face, loop_list);
}

#ifdef PROE
void GeometryQueryTool::ref_parts (DLIList<RefPart*> &ref_parts)
{
	RefEntityFactory::instance()->ref_parts(ref_parts);
}

void GeometryQueryTool::ref_assemblies (DLIList<RefAssembly*> &ref_assemblies)
{
	RefEntityFactory::instance()->ref_assemblies(ref_assemblies);
}
#endif

int GeometryQueryTool::num_bodies() const
{
  return RefEntityFactory::instance()->num_bodies();
}

int GeometryQueryTool::num_ref_volumes() const
{
  return RefEntityFactory::instance()->num_ref_volumes();
}

int GeometryQueryTool::num_ref_groups() const
{
  return RefEntityFactory::instance()->num_ref_groups();
}

int GeometryQueryTool::num_ref_faces() const
{
  return RefEntityFactory::instance()->num_ref_faces();
}

int GeometryQueryTool::num_ref_edges() const
{
  return RefEntityFactory::instance()->num_ref_edges();
}

int GeometryQueryTool::num_ref_vertices() const
{
  return RefEntityFactory::instance()->num_ref_vertices();
}

RefEntity *GeometryQueryTool::get_ref_entity(const char *type, int id)
{
  return RefEntityFactory::instance()->get_ref_entity(type, id);
}

RefEntity *GeometryQueryTool::get_ref_entity(const std::type_info& type, int id)
{
  return RefEntityFactory::instance()->get_ref_entity(type, id);
}

Body *GeometryQueryTool::get_body( int id)
{
  return RefEntityFactory::instance()->get_body(id);
}

RefVolume *GeometryQueryTool::get_ref_volume( int id)
{
  return RefEntityFactory::instance()->get_ref_volume(id);
}

RefGroup *GeometryQueryTool::get_ref_group( int id)
{
  return RefEntityFactory::instance()->get_ref_group(id);
}

RefFace *GeometryQueryTool::get_ref_face( int id)
{
  return RefEntityFactory::instance()->get_ref_face(id);
}

RefEdge *GeometryQueryTool::get_ref_edge( int id)
{
  return RefEntityFactory::instance()->get_ref_edge(id);
}

RefVertex *GeometryQueryTool::get_ref_vertex( int id)
{
  return RefEntityFactory::instance()->get_ref_vertex(id);
}

Body *GeometryQueryTool::get_first_body()
{
  return RefEntityFactory::instance()->get_first_body();
}

RefVolume *GeometryQueryTool::get_first_ref_volume()
{
  return RefEntityFactory::instance()->get_first_ref_volume();
}

RefGroup *GeometryQueryTool::get_first_ref_group()
{
  return RefEntityFactory::instance()->get_first_ref_group();
}

RefFace *GeometryQueryTool::get_first_ref_face()
{
  return RefEntityFactory::instance()->get_first_ref_face();
}

RefEdge *GeometryQueryTool::get_first_ref_edge()
{
  return RefEntityFactory::instance()->get_first_ref_edge();
}

RefVertex *GeometryQueryTool::get_first_ref_vertex()
{
  return RefEntityFactory::instance()->get_first_ref_vertex();
}

Body *GeometryQueryTool::get_next_body()
{
  return RefEntityFactory::instance()->get_next_body();
}

RefVolume *GeometryQueryTool::get_next_ref_volume()
{
  return RefEntityFactory::instance()->get_next_ref_volume();
}

RefGroup *GeometryQueryTool::get_next_ref_group()
{
  return RefEntityFactory::instance()->get_next_ref_group();
}

RefFace *GeometryQueryTool::get_next_ref_face()
{
  return RefEntityFactory::instance()->get_next_ref_face();
}

RefEdge *GeometryQueryTool::get_next_ref_edge()
{
  return RefEntityFactory::instance()->get_next_ref_edge();
}

RefVertex *GeometryQueryTool::get_next_ref_vertex()
{
  return RefEntityFactory::instance()->get_next_ref_vertex();
}


Body *GeometryQueryTool::get_last_body()
{
  return RefEntityFactory::instance()->get_last_body();
}

RefVolume *GeometryQueryTool::get_last_ref_volume()
{
  return RefEntityFactory::instance()->get_last_ref_volume();
}

RefGroup *GeometryQueryTool::get_last_ref_group()
{
  return RefEntityFactory::instance()->get_last_ref_group();
}

RefFace *GeometryQueryTool::get_last_ref_face()
{
  return RefEntityFactory::instance()->get_last_ref_face();
}

RefEdge *GeometryQueryTool::get_last_ref_edge()
{
  return RefEntityFactory::instance()->get_last_ref_edge();
}

RefVertex *GeometryQueryTool::get_last_ref_vertex()
{
  return RefEntityFactory::instance()->get_last_ref_vertex();
}

CubitStatus GeometryQueryTool::get_free_ref_entities(DLIList<RefEntity*> &free_entities)
{
    // go through the global entity lists looking for free entities

    // this algorithm works as follows: start by marking all entities in the model;
    // - bodies that don't have virtual children are NOT free, along with their children;
    //   unmark them all
    // - for all remaining entities:
    // . entities that are unmarked are not free since they were unmarked during the
    //   body check
    // . entities that are virtual are not free, but their underlying entities may be;
    //     check this, and unmark the underlying non-free entities, but not their children
    // . entities that don't have a TopologyBridge connected to a bodysm are free
    //
    // - Removed virtual stuff - no longer required - J.Kraftcheck 10-8-2003

  int i, j;
  RefEntityFactory *REF = RefEntityFactory::instance();

    // define a macro to mark all the children
#define MARK_CHILDREN(entity, index) \
   {ref_list.clean_out(); entity->get_all_child_ref_entities(ref_list); \
   for (index=ref_list.size(); index>0; index--) ref_list.get_and_step()->marked(0);}

    // mark all entities first
  for (i = REF->num_ref_volumes(); i > 0; i--)
    REF->get_next_ref_volume()->marked(1);
  for (i = REF->num_ref_faces(); i > 0; i--)
    REF->get_next_ref_face()->marked(1);
  for (i = REF->num_ref_edges(); i > 0; i--)
    REF->get_next_ref_edge()->marked(1);
  for (i = REF->num_ref_vertices(); i > 0; i--)
    REF->get_next_ref_vertex()->marked(1);

  DLIList<RefEntity*> ref_list;

    // first, mark all the children of bodies
  Body *body;
  for (i = REF->num_bodies(); i > 0; i--) {
    body = REF->get_next_body();
    MARK_CHILDREN(body, j);
  }

    // now go through them, checking for VG and free entities
  RefVolume *volume;
  for (i = REF->num_ref_volumes(); i > 0; i--) {
    volume = REF->get_next_ref_volume();
    if (volume->marked()) {
        // go through volume's children & unmark
      MARK_CHILDREN(volume, j);
      free_entities.append(volume);
      volume->marked(0);
    }
  }

  RefFace *face;
  for (i = REF->num_ref_faces(); i > 0; i--) {
    face = REF->get_next_ref_face();
    if (face->marked()) {
      MARK_CHILDREN(face, j);
      free_entities.append(face);
      face->marked(0);
    }
  }

  RefEdge *edge;
  for (i = REF->num_ref_edges(); i > 0; i--) {
    edge = REF->get_next_ref_edge();
    if (edge->marked()) {
      MARK_CHILDREN(edge, j);
      free_entities.append(edge);
      edge->marked(0);
    }
  }

  RefVertex *vertex;
  for (i = REF->num_ref_vertices(); i > 0; i--) {
    vertex = REF->get_next_ref_vertex();
    if (vertex->marked()) {
      free_entities.append(vertex);
      vertex->marked(0);
    }
  }

  return CUBIT_SUCCESS;
}


double GeometryQueryTool::surface_angle(RefFace *ref_face_1, RefFace *ref_face_2,
                                        RefEdge *ref_edge,
                                        RefVolume *ref_volume,
                                        double frac)
{

    // check for and supply missing ref edge if necessary
  if (ref_edge == NULL) {
    ref_edge = ref_face_1->common_ref_edge(ref_face_2);
    assert (ref_edge != NULL);
  }

  if (ref_volume == NULL) {
    ref_volume = ref_face_1->common_ref_volume(ref_face_2);
    assert(ref_volume != NULL);
  }

     //Find the dihedral angle for the ref_edge
   CubitVector mid_point;
   if ( frac != .5 )
     ref_edge->position_from_fraction(frac, mid_point);
   else
     mid_point = ref_edge->center_point();


     //Now find the normals of the surfaces at this point, wrt volume
   CubitVector surf_1_norm = ref_face_1->normal_at( mid_point, ref_volume );
   CubitVector surf_2_norm = ref_face_2->normal_at( mid_point, ref_volume );

     //Now we need to get the correct normal
   CubitVector tangent_vector;

     //This gets the correct tangent with respect for to the
     //ref_face_ptr.  Do the following check for non-manifold volumes.
     //This function will assert if the ref_edge has more than
     //one co_edge for the two_faces.  But this should be okay
     //since up above this function, that should be traped for...

    //Weed-out case where one edge is shared between more
    //than 2 surfaces of the same volume
    DLIList<RefFace*> tmp_faces;
    ref_edge->ref_faces( tmp_faces );

    if( tmp_faces.size() > 2 )
    {
      int kk;
      for(kk=tmp_faces.size(); kk--;)
      {
        if( !tmp_faces.get()->is_child( ref_volume ) )
          tmp_faces.change_to(NULL);
        tmp_faces.step();
      }
      tmp_faces.remove_all_with_value( NULL );
      if( tmp_faces.size() > 2 )
        //this isn't the type of surface we are looking for...
        return 0.0;
    }


   ref_edge->tangent( mid_point, tangent_vector, ref_face_1 );
   CubitSense face_1_sense = ref_face_1->sense( ref_volume );

   if( CUBIT_REVERSED != face_1_sense && CUBIT_FORWARD != face_1_sense )
      if(!(ref_volume->is_sheet()))
         return 0.0;

   if ( face_1_sense == CUBIT_REVERSED )
     tangent_vector = -tangent_vector;

   double angle = CUBIT_PI +
     tangent_vector.vector_angle(surf_2_norm, surf_1_norm );
     // note 1 and 2 switched.
     // tangent is cw, opposite of rhr

     // Range of above is pi/2 to 5pi/2, shift to 0 to 2pi.
   if ( angle >= 2.0* CUBIT_PI )
     angle -= 2.0 * CUBIT_PI;

     // if the angle is close to 0 or 2pi, use a center point check
     // to figure out if the interior angle is zero or 2pi
   const double min_tol = 0.1; // CUBIT_RESABS*100. too small
     // We could have a real overlap in some models, so keep
     // min_tol fairly big.
   const double max_tol = 2.0*CUBIT_PI - min_tol;

     // angle near zero - 2pi, make some checks at nearby locations
   if (angle <= min_tol || angle >= max_tol )
   {
       // try to get the inside-outness where the surfaces have a gap
       // between them

       // get a non-shared edge with a shared vertex between the surfaces

     int i;

       // get the next edge in the loop
     DLIList<DLIList<RefEdge*> > loops_1;
     DLIList<RefEdge*> *loop_1 = NULL;
     ref_face_1->ref_edge_loops( loops_1 );
     for ( i = loops_1.size(); i--; )
     {
       loop_1 = &loops_1.get_and_step();
       if ( loop_1->move_to( ref_edge ) )
         break;
     }
     assert( loop_1->get() == ref_edge );

     DLIList<DLIList<RefEdge*> > loops_2;
     DLIList<RefEdge*> *loop_2 = NULL;
     ref_face_2->ref_edge_loops( loops_2 );
     for ( i = loops_2.size(); i--; )
     {
       loop_2 = &loops_2.get_and_step();
       if ( loop_2->move_to( ref_edge ) )
         break;
     }
     assert( loop_2->get() == ref_edge );

     RefEdge *common_edge = ref_edge, *next_edge,
       *uncommon_edge_1 = NULL, *uncommon_edge_2 = NULL;
     RefVertex *common_vertex = NULL;

     for ( i = loop_1->size(); !common_vertex && i--; )
     {
       next_edge = loop_1->step_and_get();
       if ( loop_2->prev() == next_edge )
       {
         loop_2->back();
         common_edge = next_edge;
       }
       else if ( loop_2->next() == next_edge )
       {
         loop_2->step();
         common_edge = next_edge;
       }
       else
       {
         uncommon_edge_1 = next_edge;

           // if both vertices are shared, it doesn't matter which we chose.
         common_vertex = common_edge->common_ref_vertex( uncommon_edge_1 );

         if ( loop_2->next()->is_parent( common_vertex ) )
           uncommon_edge_2 = loop_2->next();
         else
         {
           assert( loop_2->prev()->is_parent( common_vertex ) );
           uncommon_edge_2 = loop_2->prev();
         }
       }
     }

     int too_far;

     CubitVector center_1, center_2, center_norm_1, center_norm_2;

       // we've found non-common edges with a common vertex
     int bad_angle = CUBIT_TRUE;
     if ( common_vertex )
     {

         // These two curves are only good geometrically if they
         // have a small angle between them: If the angle is too big,
         // then the closest pt will be the vertex, and we'll have to
         // start over with the surfaces or something.

         // reset midpoint and normals to common vertex...
       CubitVector vertex_coord = common_vertex->coordinates();

         // get a pair of close points on each surface
       center_1 = uncommon_edge_1->center_point();
       center_2 = uncommon_edge_2->center_point();

       double d_1 = (center_1 - vertex_coord).length_squared();
       double d_2 = (center_2 - vertex_coord).length_squared();
       int give_up = 0;
       i = 0;
       if ( d_1 <= d_2 )
       {
         do
         {
           center_2 = center_1;
           uncommon_edge_2->move_to_curve( center_2 );
           d_2 = (center_2 - vertex_coord).length_squared();
           bad_angle =  d_2 < d_1 * 0.2;
           give_up = ( d_1 < GEOMETRY_RESABS*10 || i++ > 10 );
           if ( give_up )
             break;
           if ( bad_angle )
           {
             center_1 += vertex_coord;
             center_1 /= 2.0;
             uncommon_edge_1->move_to_curve( center_1 );
             d_1 = (center_1 - vertex_coord).length_squared();
           }
         } while ( bad_angle );
       }
       else
       {
         do
         {
           center_1 = center_2;
           uncommon_edge_1->move_to_curve( center_1 );
           d_1 = (center_1 - vertex_coord).length_squared();
           bad_angle =  d_1 < d_2 * 0.2;
           give_up = ( d_2 < GEOMETRY_RESABS*10 || i++ > 10 );
           if ( give_up )
             break;
           if ( bad_angle )
           {
             center_2 += vertex_coord;
             center_2 /= 2.0;
             uncommon_edge_2->move_to_curve( center_2 );
             d_2 = (center_2 - vertex_coord).length_squared();
           }
         } while ( bad_angle );
       }
       if ( !bad_angle )
       {
         mid_point = vertex_coord;
         surf_1_norm = ref_face_1->normal_at( mid_point, ref_volume );
         surf_2_norm = ref_face_2->normal_at( mid_point, ref_volume );


         double best_d_1 = CUBIT_DBL_MAX;

         CubitVector test_center_1 = center_1;
         CubitVector test_center_norm_1;

         // may be too far away - make sure normal is roughly the same
         //  as at the midpoint

         too_far = CUBIT_TRUE;
         for ( i = 12; i-- && too_far; )
         {
           test_center_norm_1 =
             ref_face_1->normal_at( test_center_1, ref_volume );
           d_1 = test_center_norm_1 % surf_1_norm;
           if ( d_1 < best_d_1 )
           {
             center_1 = test_center_1;
             center_norm_1 = test_center_norm_1;
           }
           too_far =  d_1 < 0.2;
           if ( too_far && i )  // skip last time
           {
             test_center_1 += mid_point;
             test_center_1 /= 2.0;
             uncommon_edge_1->move_to_curve( test_center_1 );
           }
         }

         double best_d_2 = CUBIT_DBL_MAX;

         CubitVector test_center_2 = center_2;
         CubitVector test_center_norm_2;

           // may be too far away - make sure normal is roughly the same
           //  as at the surface midpoint
         too_far = CUBIT_TRUE;
         for ( i = 12; i-- && too_far; )
         {
           test_center_norm_2 =
             ref_face_2->normal_at( test_center_2, ref_volume );
           d_2 = test_center_norm_2 % surf_2_norm;
           if ( d_2 < best_d_2 )
           {
             center_2 = test_center_2;
             center_norm_2 = test_center_norm_2;
           }
           too_far =  d_2 < 0.2;
           if ( too_far && i )  // skip last time
           {
             test_center_2 += mid_point;
             test_center_2 /= 2.0;
             uncommon_edge_2->move_to_curve( test_center_2 );
           }
         }
       }
     }
       // surfaces share all edges! try the face center point
     if ( !common_vertex || bad_angle )
     {

         // get a pair of close points on each surface
       center_1 = ref_face_1->center_point();
       center_2 = ref_face_2->center_point();
       double d_1 = (center_1 - mid_point).length_squared();
       double d_2 = (center_2 - mid_point).length_squared();
       if ( d_1 <= d_2 )
       {
         center_2 = center_1;
         ref_face_2->move_to_surface( center_2 );
       }
       else
       {
         center_1 = center_2;
         ref_face_1->move_to_surface( center_1 );
       }

       double best_d_1 = CUBIT_DBL_MAX;

       CubitVector test_center_1 = center_1;
       CubitVector test_center_norm_1;

         // may be too far away - make sure normal is roughly the same
         //  as at the curve midpoint
       too_far = CUBIT_TRUE;
       for ( i = 12; i-- && too_far; )
       {
         test_center_norm_1 =
           ref_face_1->normal_at( test_center_1, ref_volume );
         d_1 = test_center_norm_1 % surf_1_norm;
         if ( d_1 < best_d_1 )
         {
           center_1 = test_center_1;
           center_norm_1 = test_center_norm_1;
         }
         too_far =  d_1 < 0.2;
         if ( too_far && i )  // skip last time
         {
           test_center_1 += mid_point;
           test_center_1 /= 2.0;
           ref_face_1->move_to_surface( test_center_1 );
         }
       }

         // surfaces share all edge! try the center point
       double best_d_2 = CUBIT_DBL_MAX;

       CubitVector test_center_2 = center_2;
       CubitVector test_center_norm_2;

         // may be too far away - make sure normal is roughly the same
         //  as at the curve midpoint
       too_far = CUBIT_TRUE;
       for ( i = 12; i-- && too_far; )
       {
         test_center_norm_2 =
           ref_face_2->normal_at( test_center_2, ref_volume );
         d_2 = test_center_norm_2 % surf_2_norm;
         if ( d_2 < best_d_2 )
         {
           center_2 = test_center_2;
           center_norm_2 = test_center_norm_2;
         }
         too_far =  d_2 < 0.2;
         if ( too_far && i )  // skip last time
         {
           test_center_2 += mid_point;
           test_center_2 /= 2.0;
           ref_face_2->move_to_surface( test_center_2 );
         }
       }
     }

       // gap vector from center_1 to center_2
     CubitVector gap_vector = center_2;
     gap_vector -= center_1;

     double gap_d = gap_vector % center_norm_1;
     gap_d += -gap_vector % center_norm_2;

     if ( gap_d >= 0 )
     {
       if ( angle < max_tol )
         angle = 2*CUBIT_PI;
         // else leave at its current nearly-2pi value
     }
     else
       if ( angle > min_tol )
         angle = 0.;
       // else leave at its current nearly-zero value

   } // if angle near zero - 2pi

   return angle;
}

//-------------------------------------------------------------------------
// Purpose       : Find a common geometry query engine
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/00
//-------------------------------------------------------------------------
GeometryQueryEngine* GeometryQueryTool::common_query_engine(
  DLIList<TopologyEntity*>& topology_list,
  DLIList<TopologyBridge*>& engine_bridges,
  CubitBoolean allow_default_engine ) const
{
  topology_list.reset();

  TopologyEntity* topo_ptr = topology_list.get_and_step();
  DLIList<TopologyBridge*> first_bridge_list;
  topo_ptr->bridge_manager()->get_bridge_list( first_bridge_list );

  first_bridge_list.reset();
  GeometryQueryEngine* gqe_ptr = 0;
  for( int i = first_bridge_list.size(); i > 0; i-- )
  {
    TopologyBridge* bridge_ptr = first_bridge_list.get_and_step();
    engine_bridges.clean_out();
    engine_bridges.append( bridge_ptr );
    gqe_ptr = bridge_ptr->get_geometry_query_engine();

    topology_list.reset();
    topology_list.step(); //skip first entry
    for( int j = topology_list.size(); j > 1; j-- )
    {
      topo_ptr = topology_list.get_and_step();
      bridge_ptr = topo_ptr->bridge_manager()->topology_bridge(gqe_ptr);
      if( bridge_ptr ) engine_bridges.append( bridge_ptr );
      else break;
    }

    if( engine_bridges.size() == topology_list.size() )
      break;

    gqe_ptr = 0;
  }

  if( !gqe_ptr )
  {
    engine_bridges.clean_out();

    if( allow_default_engine )
    {
      PRINT_WARNING("Entities do not belong to the same geometry "
        "engine.  Using the default  geometry engine.\n");
      gqe_ptr = default_gqe;
      topology_list.reset();
      for( int j = topology_list.size(); j > 0; j-- )
      {
        topo_ptr = topology_list.get_and_step();
        TopologyBridge* bridge_ptr =
          topo_ptr->bridge_manager()->topology_bridge( gqe_ptr );
        if( ! bridge_ptr )
          bridge_ptr = topo_ptr->bridge_manager()->topology_bridge();
        engine_bridges.append( bridge_ptr );
      }
    }
  }

  return gqe_ptr;
}


void GeometryQueryTool::get_connected_free_ref_entities(
  RefEntity *entity,
  const int merge_option,
  DLIList<Body*> &body_list,
  DLIList<RefFace*> &ref_face_list,
  DLIList<RefEdge*> &ref_edge_list,
  DLIList<RefVertex*> &ref_vertex_list )
{

  body_list.clean_out();
  ref_face_list.clean_out();
  ref_edge_list.clean_out();
  ref_vertex_list.clean_out();

  if ( !entity )
    return;

   // cases are different enough, just do something totally different
  if ( !merge_option  )
  {
      // get top level entities
    TopologyEntity *source_entity = CAST_TO(entity, TopologyEntity);

    Body *body = CAST_TO( source_entity, Body );
    if( body )
      body_list.append( body );
    else
        //  get attached bodies, if any
      source_entity->bodies( body_list );

    if(body_list.size() == 0)
    {
      RefEdge* ref_edge = CAST_TO(source_entity, RefEdge);
      if ( ref_edge )
        ref_edge_list.append( ref_edge );
      else
      {
        source_entity->ref_edges( ref_edge_list );
      }
      if ( ref_edge_list.size() == 0 )
      {
        RefVertex *vert = CAST_TO( source_entity, RefVertex );
        if ( vert )
          ref_vertex_list.append( vert );
      }
    }
    // this is the easy case, we're all done
    return;
  }


    // merge_option == 1
  DLIList <TopologyEntity*> source_list;
  DLIList <RefEntity*> source_entity_list;

  int i;

  DLIList <RefEntity*> free_list;

  DLIList<TopologyBridge*> bridge_list;

  TopologyEntity *topo_entity = CAST_TO( entity, TopologyEntity );
  source_list.append( topo_entity );
  DLIList<RefEntity*> entity_list;
  do
  {
      // get vertices (bottom level ref entities)
    source_entity_list.clean_out();    

    for (i = source_list.size(); i--; )
    { 
      TopologyEntity *source_entity = source_list.get_and_step();
      DLIList<RefVertex*> local_vert_list;
      DLIList<RefEntity*> tmp_source_entity_list;
      source_entity->ref_vertices( local_vert_list );
      if(local_vert_list.size() == 0)
      {
        DLIList<RefEdge*> local_edge_list;
        source_entity->ref_edges(local_edge_list);
        if(local_edge_list.size() == 0)
        {
          DLIList<RefFace*> local_face_list;
          source_entity->ref_faces(local_face_list);
          if(local_face_list.size() > 0)
          {
            CAST_LIST( local_face_list, tmp_source_entity_list, RefEntity);
            source_entity_list += tmp_source_entity_list;
          }
        }
        else
        {
          CAST_LIST( local_edge_list, tmp_source_entity_list, RefEntity);
          source_entity_list += tmp_source_entity_list;
        }
      }
      else
      {
        CAST_LIST( local_vert_list, tmp_source_entity_list, RefEntity);
        source_entity_list += tmp_source_entity_list;
      }
    }
    source_list.clean_out();

      // get top level entities
    for ( i = source_entity_list.size(); i--; )
    {
      RefEntity *source_entity = source_entity_list.get_and_step();
        // get all upwards related ref entities
      entity_list.clean_out();
        // get the bodies, too!
      source_entity->get_all_parent_ref_entities( entity_list, CUBIT_TRUE );
      entity_list.append( source_entity );

        // check each one to see if it is a new top level in the solid
        // modeller
      int j;
      for ( j = entity_list.size(); j--; )
      {
        // entity is top level if it has a topology bridge that has a bodysm
        RefEntity *ref_entity = entity_list.get_and_step(); 

        if ( !ref_entity->marked() )
        {
          bridge_list.clean_out();
          topo_entity = CAST_TO( ref_entity, TopologyEntity );
          topo_entity->bridge_manager()->get_bridge_list(bridge_list);
          int k;
          bool no_parents = true;
          for ( k = bridge_list.size(); k--; )
          {
            TopologyBridge* bridge = bridge_list.get_and_step();
            DLIList<TopologyBridge*> parents;
            bridge->get_parents( parents );            
            if( parents.size() )
            {
              no_parents = false;
              break;
            }
          }
          if( no_parents )
          {             
            ref_entity->marked( 1 );
            topo_entity = CAST_TO( ref_entity, TopologyEntity );
            source_list.append( topo_entity );
            free_list.append( ref_entity );
            continue;
          }
        }
      }
    }

  } while( source_list.size() );

    // xfer data, clean up marks
  for ( i = free_list.size(); i--; )
  {
    RefEntity *ref_entity = free_list.get_and_step();
    ref_entity->marked(0);
  }
  CAST_LIST( free_list, body_list, Body );
  CAST_LIST( free_list, ref_face_list, RefFace);
  CAST_LIST( free_list, ref_edge_list, RefEdge);
  CAST_LIST( free_list, ref_vertex_list, RefVertex);
}


CubitStatus GeometryQueryTool::register_intermediate_engine(
  IntermediateGeomEngine* engine_ptr )
{
  return (CubitStatus)igeSet.insert(engine_ptr).second;
}

void GeometryQueryTool::unregister_intermediate_engine(
  IntermediateGeomEngine* engine_ptr )
{
  igeSet.erase(engine_ptr);
}

void GeometryQueryTool::ige_export_geom( DLIList<TopologyBridge*> &geom_list )
{
  for (IGESet::reverse_iterator itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->export_geometry(geom_list);
}

void GeometryQueryTool::ige_push_imprint_attributes_before_modify
                                ( DLIList<BodySM*> &geom_list )
{
  for (IGESet::reverse_iterator itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->push_imprint_attributes_before_modify(geom_list);
}

void GeometryQueryTool::ige_push_named_attributes_to_curves_and_points
                                ( DLIList<TopologyBridge*> &tb_list, const char *name_in )
{
  for (IGESet::reverse_iterator itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->push_named_attributes_to_curves_and_points(tb_list, name_in);
}

void GeometryQueryTool::ige_remove_imprint_attributes_after_modify
                                ( DLIList<BodySM*> &old_sms,
                                  DLIList<BodySM*> &new_sms )
{
  for (IGESet::reverse_iterator itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
    (*itor)->remove_imprint_attributes_after_modify(old_sms, new_sms);
}

bool GeometryQueryTool::ige_is_composite(TBOwner *bridge_owner)
{
  bool ret = false;
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end() && ret != true; ++itor)
  {
    if((*itor)->is_composite(bridge_owner))
    {
      ret = true;
    }
  }
  return ret;
}

bool GeometryQueryTool::ige_is_composite(TopologyBridge *bridge)
{
  bool ret = false;
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end() && ret != true; ++itor)
  {
    if((*itor)->is_composite(bridge))
    {
      ret = true;
    }
  }
  return ret;
}


bool GeometryQueryTool::ige_is_partition(TBOwner *bridge_owner)
{
  bool ret = false;
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end() && ret != true; ++itor)
  {
    if((*itor)->is_partition(bridge_owner))
    {
      ret = true;
    }
  }
  return ret;
}

void GeometryQueryTool::ige_import_geom( DLIList<TopologyBridge*> &geom_list )
{
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->import_geometry(geom_list);
}

void GeometryQueryTool::get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                                            DLIList<TopologyBridge*> &tbs )
{
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->get_tbs_with_bridge_manager_as_owner( source_bridge, tbs );
}

void GeometryQueryTool::ige_attribute_after_imprinting(DLIList<TopologyBridge*> &tb_list,
                                                       DLIList<Body*> &old_bodies)
{
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->attribute_after_imprinting(tb_list, old_bodies);
}

void GeometryQueryTool::ige_remove_attributes( DLIList<TopologyBridge*> &geom_list )
{
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->remove_attributes(geom_list);
}

void GeometryQueryTool::ige_remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges)
{
  for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->remove_attributes_from_unmodifed_virtual(bridges);
}

bool GeometryQueryTool::contains_intermediate_geometry(DLIList<RefEntity*>& ref_entity_list) const
{
  // Unfortunately, the current implementation of partition creates virtual
  // bodies as well as virtual subentities.  Thus, we have to go up the tree
  // as well as down it.
  // TODO: Partitioning HAS got to change.  KGM 2/9/06
  //get the owning bodies
  DLIList<Body*> body_list;
  int j;
  for(j=ref_entity_list.size(); j--;)
  {
    TopologyEntity* te = dynamic_cast<TopologyEntity*>( ref_entity_list.get_and_step());
    te->bodies( body_list );
  }
  int i;
  for ( i = 0; i < body_list.size(); i++)
    if (GeometryQueryTool::instance()->
        is_intermediate_geometry(body_list.next(i)))
      return true;

  for ( i = 0; i < ref_entity_list.size(); i++)
    if (GeometryQueryTool::instance()->
        contains_intermediate_geometry(ref_entity_list.next(i)))
      return true;

  return false;
}

bool GeometryQueryTool::contains_intermediate_geometry(RefEntity* entity_ptr) const
{
  if (igeSet.empty())
    return false;

  DLIList<RefEntity*> children;
  entity_ptr->get_all_child_ref_entities(children);
  children.append(entity_ptr);

  while(children.size())
    if (is_intermediate_geometry(children.pop()))
      return true;

  return false;
}

bool GeometryQueryTool::is_intermediate_geometry(RefEntity* entity_ptr) const
{
  TopologyEntity* topo_ptr = CAST_TO(entity_ptr, TopologyEntity);
  if (!topo_ptr)
    return false;

  DLIList<TopologyBridge*> bridge_list;
  topo_ptr->bridge_manager()->get_bridge_list(bridge_list);
  while(bridge_list.size())
    if (is_intermediate_geometry(bridge_list.pop()))
      return true;

  return false;
}

bool GeometryQueryTool::is_intermediate_geometry( TopologyBridge* bridge ) const
{
  if (bridge->get_geometry_query_engine() == NULL )
    return true;
  else
    return bridge->get_geometry_query_engine()->is_intermediate_engine();
}


//-------------------------------------------------------------------------
// Purpose       : Destroy dead entity and dead children
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/08/03
//-------------------------------------------------------------------------
CubitStatus GeometryQueryTool::destroy_dead_entity(
                                            TopologyEntity* topo_ent,
                                            bool top ) const
{
  if (topo_ent->get_parents() || topo_ent->bridge_manager()->topology_bridge())
    return CUBIT_FAILURE;

  topo_ent->deactivated(CUBIT_TRUE);

  CubitObservable* ob = dynamic_cast<CubitObservable*>(topo_ent);
  if (ob)
  {
    if (dynamic_cast<RefEntity*>(topo_ent))
    {
        // "top" indicates if this call is the topmost call to this
        // recursive function.  It should be the case that if it is
        // the topmost call and the passed entity is a RefEntity, then
        // the entity was top-level (had no parent entities in the
        // topology graph.)  For cases where dead topology is cleaned
        // out and a dead RefEntity is not-top-most, it will have some
        // parent sense entity which this function will be called on and
        // that call will be the top-most one.
      if (top)
      {
        AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOP_LEVEL_ENTITY_DESTRUCTED, static_cast<RefEntity*>(ob)));
        CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_DELETED, static_cast<RefEntity*>(ob));
        const_cast<CGMHistory&>(mHistory).add_event(evt);
      }
      AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::TOPOLOGY_ENTITY_DESTRUCTED, static_cast<RefEntity*>(ob)));
    }
    else
      AppUtil::instance()->send_event(TopologyEvent(TopologyEvent::TOPOLOGY_ENTITY_DESTRUCTED, topo_ent));
  }

  DLIList<TopologyEntity*> child_list;
  topo_ent->disconnect_all_children(&child_list);
  child_list.reset();
  for (int i = child_list.size(); i--; )
  {
    TopologyEntity* child = child_list.get_and_step();
    TopologyEntity* child_topo_ent = dynamic_cast<TopologyEntity*>(child);
    if (child_topo_ent)
    {
      destroy_dead_entity(child_topo_ent, false);
      if (!child_topo_ent->deactivated() &&
           child_topo_ent->bridge_manager()->number_of_bridges() > 0 &&
           child_topo_ent->get_parents() == 0 &&
           NULL != (ob = dynamic_cast<CubitObservable*>(child_topo_ent)))
      {
        DLIList<TopologyBridge*> list1, list2;
        bool has_parents = false;
        child_topo_ent->bridge_manager()->get_bridge_list(list1);
        while (list1.size())
        {
          list2.clean_out();
          list1.pop()->get_parents(list2);
          if (list2.size())
            has_parents = true;
        }
        if (!has_parents)
        {
          AppUtil::instance()->send_event(GeometryEvent(GeometryEvent::FREE_REF_ENTITY_GENERATED, static_cast<RefEntity*>(ob)));
          CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_CREATED, static_cast<RefEntity*>(ob));
          const_cast<CGMHistory&>(mHistory).add_event(evt);
        }
      }
    }
  }

  if (top)
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();

  return CUBIT_SUCCESS;
}

// TODO - find out where this functionality belongs
CubitStatus GeometryQueryTool::import_actuate(DLIList<RefEntity*> &entity_list)
{
  int i;

    // given a Body list, actuates first the merge partner attribute
    // on all entities in the list, then actuates all other types of
    // attributes

  DLIList<TopologyEntity*> temp_list;
  DLIList<RefEntity*> refent_list;

  CubitBoolean auto_actuate_merge;
  auto_actuate_merge = CGMApp::instance()->attrib_manager()->auto_actuate_flag(CA_MERGE_PARTNER);


  if (auto_actuate_merge) {

      // Prevent MergeTool from destroying dead entities
      // after merging.  Otherwise we get stale pointers
      // in our lists.
      // Don't forget to turn this back on later!
    MergeTool::destroy_dead_geometry( false );

    DLIList<TopologyEntity*> me_list;
    CAST_LIST(entity_list, me_list, TopologyEntity);

    ModelQueryEngine *const mqe = ModelQueryEngine::instance();

    // vertices
    mqe->query_model(me_list, DagType::ref_vertex_type(), temp_list);
    CAST_LIST(temp_list, refent_list, RefEntity);
      // actuate merge attribute for vertices
    if (refent_list.size() > 0)
      refent_list.get()->actuate_cubit_attrib(refent_list, CA_MERGE_PARTNER);

    // edges
    temp_list.clean_out();
    refent_list.clean_out();
    mqe->query_model(me_list, DagType::ref_edge_type(), temp_list);
    CAST_LIST(temp_list, refent_list, RefEntity);
      // actuate merge attribute for edges
    if (refent_list.size() > 0)
      refent_list.get()->actuate_cubit_attrib(refent_list, CA_MERGE_PARTNER);

    // faces
    temp_list.clean_out();
    refent_list.clean_out();
    mqe->query_model(me_list, DagType::ref_face_type(), temp_list);
    CAST_LIST(temp_list, refent_list, RefEntity);
      // actuate merge attribute for faces
    if (refent_list.size() > 0)
      refent_list.get()->actuate_cubit_attrib(refent_list, CA_MERGE_PARTNER);

      // clean out entities destroyed during merge
    entity_list.reset();
    for( i = entity_list.size(); i--; )
    {
      TopologyEntity* me_ptr = dynamic_cast<TopologyEntity*>(entity_list.get());
      if( me_ptr && me_ptr->deactivated() )
        entity_list.extract();
      else
        entity_list.step();
    }

      // Restore merge tool setting, and clean up dead geometry
    MergeTool::destroy_dead_geometry( true );
    GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  }

    // now actuate other attributes
  RefEntity* entity_ptr;
  CAActuateSet actuate_set( entity_list );

    // actuate in increasing dimension
    // (vtx = 0, ..., volume = 3, body = 4)
  for( i = 0; i < 5; i++ )
  {
    actuate_set.set_current_dimension( i );
    while( (entity_ptr = actuate_set.remove_next() ) != NULL )
      entity_ptr->auto_actuate_cubit_attrib(CUBIT_FALSE);
  }

    // actuate deferred attribs
  CADeferredAttrib::cleanup_cadas(CUBIT_TRUE, CUBIT_TRUE);

    // finally, actuate the attributes that go after all other geometry changes
  for( i = 0; i < 4; i++ )
  {
    actuate_set.set_current_dimension( i );
    while( (entity_ptr = actuate_set.remove_next() ) != NULL )
      entity_ptr->auto_actuate_cubit_attrib(CUBIT_FALSE, CUBIT_TRUE);
  }

  return CUBIT_SUCCESS;
}

CubitBoolean GeometryQueryTool::okay_to_transform( Body* body ) const
{
  MergeTool* mt = MergeTool::instance();
  int i;

    // Check for merged vertices
  DLIList<RefVertex*> vertices;
  body->ref_vertices( vertices );
  for (i = vertices.size(); i--; )
    if ( mt->entity_merged( vertices.get_and_step() ) )
    {
      PRINT_ERROR("Cannot transform %s (Body %d) with merged geomtery.\n",
          body->entity_name().c_str(), body->id() );
      return CUBIT_FALSE;
    }

    // Need to check for merged surfaces because can have surfaces
    // w/out vertices.
  DLIList<RefFace*> surfaces;
  body->ref_faces( surfaces );
  for (i = surfaces.size(); i--; )
    if ( mt->entity_merged( surfaces.get_and_step() ) )
    {
      PRINT_ERROR("Cannot transform %s (Body %d) with merged geomtery.\n",
          body->entity_name().c_str(), body->id() );
      return CUBIT_FALSE;
    }

  return CUBIT_TRUE;
}

void GeometryQueryTool::translate( DLIList<RefEntity*> &entities_to_transform,
        double x, double y, double z, bool check_before_transforming,
        DLIList<RefEntity*> &entities_transformed,
        bool preview /*= false*/)
{
  CubitVector delta(x,y,z);

  //translate free, merged-away entities first
  DLIList<TopologyBridge*> free_ents; 
  get_merged_away_free_entities( entities_to_transform, free_ents );

  int i;
  if (!preview)
  {
    for( i=free_ents.size(); i--; )
    {
      TopologyBridge *bridge = free_ents.get_and_step();
      Curve *curve= CAST_TO( bridge, Curve );
      TBPoint *point = CAST_TO( bridge, TBPoint );

      if( curve || point )
      {
        GeometryEntity *geom = CAST_TO( bridge, GeometryEntity );
        GeometryQueryEngine* engine = geom->get_geometry_query_engine();
        CubitStatus result = engine->translate( geom, delta );
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
          return;
        }
      }
    }
  }

  RefFace *tmp_face;
  RefEdge *tmp_curve;
  RefVertex *tmp_vertex;
  CubitStatus result = CUBIT_SUCCESS;

  DLIList<Body*> bodies_to_translate;
  DLIList<BasicTopologyEntity*> ents_to_translate;
  CAST_LIST( entities_to_transform, bodies_to_translate, Body);

  if( bodies_to_translate.size() != entities_to_transform.size() )
  {
    for( int k=0; k<entities_to_transform.size(); k++ )
    {
      RefEntity *ent = entities_to_transform[k];

      if( ( tmp_face = CAST_TO( ent, RefFace ) ) != NULL )
        ents_to_translate.append( tmp_face );
      else if( (tmp_curve = CAST_TO( ent, RefEdge ) ) != NULL )
        ents_to_translate.append( tmp_curve );
      else if( (tmp_vertex = CAST_TO( ent, RefVertex ) ) != NULL )
        ents_to_translate.append( tmp_vertex );
    }
  }


  if( bodies_to_translate.size() )
  {
    DLIList<Body*> bodies_translated;
    result = translate( bodies_to_translate, 
      CubitVector(x,y,z), 
      &bodies_translated,
      check_before_transforming,
      preview );  

    if( result )
    {
      for( int k=0; k<bodies_translated.size(); k++ )
        entities_transformed.append( bodies_translated[k] );
    }
  }

  if( ents_to_translate.size() )
  {
    DLIList<BasicTopologyEntity*> btes_translated;
     result = translate( ents_to_translate, 
      CubitVector(x,y,z), 
      &btes_translated,
      check_before_transforming,
      preview );  

    if( result )
    {
      for( int k=0; k<btes_translated.size(); k++ )
        entities_transformed.append( btes_translated[k] );
    }
  }
}


CubitStatus GeometryQueryTool::translate( DLIList<Body*> &bodies,
                                          const CubitVector& delta,
                                          DLIList<Body*> *bodies_translated,
                                          bool check_to_transform,
                                          bool preview )
{
  CubitTransformMatrix xform;
  xform.translate( delta );    

  DLIList<RefEntity*> ents_transformed;

  for( int k=0; k<bodies.size(); k++ )
  {
    Body *body = bodies[k];

    if( check_to_transform )
      if (!okay_to_transform( body ))
        continue;

    if (preview)
    {
      DLIList<RefEdge*> edges;
      body->ref_edges(edges);
      if( edges.size() )
      {
        for (int i = 0; i < edges.size(); i++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
          {
            poly.transform(xform);
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
          {
            CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
            tmp_pt = xform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
      }
      else
      {
        //just draw the surfaces
        DLIList<RefFace*> faces;
        body->ref_faces( faces );
        for( int i=0; i<faces.size(); i-- )
        {
          GMem poly;
          faces.get_and_step()->get_graphics( poly );
          poly.transform(xform);

          int* facet_list = poly.facet_list();
          GPoint* plist = poly.point_list();

          GPoint p[3];
          for (i = 0; i < poly.fListCount; )
          {
            int sides = facet_list[i++];
            if (sides != 3)
            {
              i += sides;
              continue;
            }
            else
            {
              p[0] = plist[facet_list[i++]];
              p[1] = plist[facet_list[i++]];
              p[2] = plist[facet_list[i++]];
              GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
            }
          }
        }
      }
      GfxPreview::flush();
      continue;
    }

    BodySM* bodysm = body->get_body_sm_ptr();
    GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();
    CubitStatus result = engine->translate( bodysm, delta );
    if (result)
    {
      notify_intermediate_of_transform( body, xform );
      
      if( bodies_translated )
        bodies_translated->append( body );
      ents_transformed.append( body );
    }
    else
      PRINT_ERROR("Translate of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );
  }

  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &xform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;  
}

CubitStatus GeometryQueryTool::rotate( DLIList<RefEntity*> &entities_to_transform,  
                      const CubitVector& point,
                      const CubitVector& direction,
                      double angle,
                      bool check_before_transforming,
                      DLIList<RefEntity*> &entities_transformed,
                      bool preview /*= false*/)
{
  //rotate free, merged-away entities first
  DLIList<TopologyBridge*> free_ents; 
  get_merged_away_free_entities( entities_to_transform, free_ents );
  
  int i;
  if (!preview)
  {
    for( i=free_ents.size(); i--; )
    {
      TopologyBridge *bridge = free_ents.get_and_step();
      Curve *curve= CAST_TO( bridge, Curve );
      TBPoint *tmp_point = CAST_TO( bridge, TBPoint );

      if( curve || tmp_point )
      {
        GeometryEntity *geom = CAST_TO( bridge, GeometryEntity );
        GeometryQueryEngine* engine = geom->get_geometry_query_engine();
        CubitStatus result = engine->translate( geom, -point ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
          return result;
        }
        result = engine->rotate( geom, direction, angle ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::rotate failed.\n");
          return result;
        }
        result = engine->translate( geom, point ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
          return result;
        }
      }
    }
  }
  
  RefFace *tmp_face;
  RefEdge *tmp_curve;
  RefVertex *tmp_vertex;
  CubitStatus result = CUBIT_SUCCESS;

  DLIList<Body*> bodies_to_rotate;
  DLIList<BasicTopologyEntity*> ents_to_rotate;
  CAST_LIST( entities_to_transform, bodies_to_rotate, Body);

  if( bodies_to_rotate.size() != entities_to_transform.size() )
  {
    for( int k=0; k<entities_to_transform.size(); k++ )
    {
      RefEntity *ent = entities_to_transform[k];

      if( ( tmp_face = CAST_TO( ent, RefFace ) ) != NULL )
        ents_to_rotate.append( tmp_face );
      else if( (tmp_curve = CAST_TO( ent, RefEdge ) ) != NULL )
        ents_to_rotate.append( tmp_curve );
      else if( (tmp_vertex = CAST_TO( ent, RefVertex ) ) != NULL )
        ents_to_rotate.append( tmp_vertex );
    }
  }


  if( bodies_to_rotate.size() )
  {
    DLIList<Body*> bodies_rotated;
    
    result = rotate( bodies_to_rotate, 
      point, direction, angle,
      &bodies_rotated,
      check_before_transforming, preview); 

    if( result )
    {
      for( int k=0; k<bodies_rotated.size(); k++ )
        entities_transformed.append( bodies_rotated[k] );
    }
  }

  if( ents_to_rotate.size() )
  {
    DLIList<BasicTopologyEntity*> btes_rotated;
    
     result = rotate( ents_to_rotate, 
      point,
      direction,
      angle,
      &btes_rotated,
      check_before_transforming,
      preview );  

    if( result )
    {
      for( int k=0; k<btes_rotated.size(); k++ )
        entities_transformed.append( btes_rotated[k] );
    }
  }

  if( entities_transformed.size() )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryQueryTool::rotate( DLIList<Body*> &bodies,
                                       const CubitVector& axis,
                                       double angle,
                                       DLIList<Body*> *bodies_rotated,
                                       bool check_to_transform,
                                       bool preview )
{
  CubitTransformMatrix xform;
  xform.rotate( angle, axis );

  DLIList<RefEntity*> ents_transformed;

  for( int k=0; k<bodies.size(); k++ )
  {
    Body *body = bodies[k];

    if( check_to_transform )
      if (!okay_to_transform( body ))
        continue;

    if (preview)
    {
      DLIList<RefEdge*> edges;
      body->ref_edges(edges);
      if( edges.size() )
      {
        for (int i = 0; i < edges.size(); i++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
          {
            poly.transform(xform);
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
          {
            CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
            tmp_pt = xform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
      }
      else
      {
        //just draw the surfaces
        DLIList<RefFace*> faces;
        body->ref_faces( faces );
        for( int i=0; i<faces.size(); i-- )
        {
          GMem poly;
          faces.get_and_step()->get_graphics( poly );
          poly.transform(xform);

          int* facet_list = poly.facet_list();
          GPoint* plist = poly.point_list();

          GPoint p[3];
          for (i = 0; i < poly.fListCount; )
          {
            int sides = facet_list[i++];
            if (sides != 3)
            {
              i += sides;
              continue;
            }
            else
            {
              p[0] = plist[facet_list[i++]];
              p[1] = plist[facet_list[i++]];
              p[2] = plist[facet_list[i++]];
              GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
            }
          }
        }
      }
      GfxPreview::flush();
      continue;
    }

    BodySM* bodysm = body->get_body_sm_ptr();
    GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();
    CubitStatus result = engine->rotate( bodysm, axis, angle );

    if (result)
    {
      notify_intermediate_of_transform( body, xform );
      
      if( bodies_rotated )
        bodies_rotated->append( body );
      ents_transformed.append( body );     
    }
    else
      PRINT_ERROR("Rotate of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );
  }

  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &xform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryQueryTool::rotate( DLIList<Body*> &bodies,
                                       const CubitVector& point,
                                       const CubitVector& direction,
                                       double degrees,
                                       DLIList<Body*> *bodies_rotated,
                                       bool check_to_transform,
                                       bool preview /*=false*/)
{ 
  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix rot_mat;
  rot_mat.rotate( degrees, direction );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );
  
  CubitTransformMatrix total_transform;
  total_transform = mov_mat*rot_mat*prev_xform;    

  DLIList<RefEntity*> ents_transformed;

  for( int k=0; k<bodies.size(); k++ )
  {
    Body *body = bodies[k];

    if( check_to_transform )
      if (!okay_to_transform( body ))
        continue;

    if (preview)
    {
      DLIList<RefEdge*> edges;
      body->ref_edges(edges);
      if( edges.size() )
      {
        for (int i = 0; i < edges.size(); i++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
          {
            poly.transform( total_transform );
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
          {
            CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();          
            tmp_pt = total_transform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
      }
      else
      {
        //just draw the surfaces
        DLIList<RefFace*> faces;
        body->ref_faces( faces );
        for( int i=0; i<faces.size(); i-- )
        {
          GMem poly;
          faces.get_and_step()->get_graphics( poly );        
          poly.transform( total_transform );

          int* facet_list = poly.facet_list();
          GPoint* plist = poly.point_list();

          GPoint p[3];
          for (i = 0; i < poly.fListCount; )
          {
            int sides = facet_list[i++];
            if (sides != 3)
            {
              i += sides;
              continue;
            }
            else
            {
              p[0] = plist[facet_list[i++]];
              p[1] = plist[facet_list[i++]];
              p[2] = plist[facet_list[i++]];
              GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
            }
          }
        }
      }
      GfxPreview::flush();
      continue;
    }

    BodySM* bodysm = body->get_body_sm_ptr();
    GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();
    CubitStatus result = CUBIT_FAILURE;

    // Move to origin
    result = engine->translate( bodysm, -point );

    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( -point );
      notify_intermediate_of_transform( body, xform );
    }
    else
    {
      continue;
    }

    // Rotate about direction vector
    result = engine->rotate( bodysm, direction, degrees );
    if (result)
    {
      CubitTransformMatrix rot_mat;
      rot_mat.rotate( degrees, direction );
      notify_intermediate_of_transform( body, rot_mat );
    }
    else
    {
      continue;
    }

    result= engine->translate( bodysm, point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( point );
      notify_intermediate_of_transform( body, xform );

      if( bodies_rotated )
        bodies_rotated->append( body );
      ents_transformed.append( body );  
    }
    else
    {
      continue;
    }
  }

  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &total_transform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

void GeometryQueryTool::scale( DLIList<RefEntity*> &entities_to_transform, 
                              const CubitVector& point,
                              double scale_x, double scale_y, double scale_z, 
                              bool check_before_transforming, 
                              DLIList<RefEntity*> &entities_scaled,
                              bool preview /*= false*/)
{
  CubitVector factors(scale_x, scale_y, scale_z);

  //scale free, merged-away entities first
  DLIList<TopologyBridge*> free_ents; 
  get_merged_away_free_entities( entities_to_transform, free_ents );
  int i;
  if (!preview)
  {
      for( i=free_ents.size(); i--; )
      {
          TopologyBridge *bridge = free_ents.get_and_step();
          Curve *curve= CAST_TO( bridge, Curve );
          TBPoint *tbpoint = CAST_TO( bridge, TBPoint );

          if( curve || tbpoint )
          {
              GeometryEntity *geom = CAST_TO( bridge, GeometryEntity );
              GeometryQueryEngine* engine = geom->get_geometry_query_engine();
              CubitStatus result = engine->translate( geom, -point );
              if (CUBIT_SUCCESS != result) {
                PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
                return;
              }
              result = engine->scale( geom, factors );
              if (CUBIT_SUCCESS != result) {
                PRINT_ERROR("GeometryQueryEngine::scale failed.\n");
                return;
              }
              result = engine->translate( geom, point );
              if (CUBIT_SUCCESS != result) {
                PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
                return;
              }
          }
      }
  }
  

  CubitStatus result = CUBIT_SUCCESS;
  for(i=entities_to_transform.size(); i--; )
  {
    RefEntity *tmp_ent = entities_to_transform.get_and_step();
    Body *tmp_body;
    RefFace* tmp_face;
    RefEdge *tmp_curve;
    RefVertex *tmp_vertex;
    if( ( tmp_body = CAST_TO( tmp_ent, Body ) ) != NULL )
    {
      //non-uniform scaling
      if( scale_x != scale_y ||
          scale_x != scale_z ||
          scale_y != scale_z )
      {
        // use GMT version for non-uniform scaling b/c it updates topology if it changes
 
        result = GeometryModifyTool::instance()->scale(tmp_body,point, factors, 
                                                  check_before_transforming, preview,false);
        tmp_ent = tmp_body;
      }
      else
        result = scale( tmp_body,point, CubitVector(scale_x, scale_y, scale_z),
                 check_before_transforming, preview);
    }
    else if( ( tmp_face = CAST_TO( tmp_ent, RefFace ) ) != NULL )
    {
     // only allow scaling of RefFaces if preview is on
     if (!preview)
         continue;
     result = scale( tmp_face,point, CubitVector(scale_x, scale_y, scale_z),
                check_before_transforming, preview);
    }
    else if( ( tmp_curve = CAST_TO( tmp_ent, RefEdge ) ) != NULL )
    {
     result = scale( tmp_curve,point, CubitVector(scale_x, scale_y, scale_z),
                check_before_transforming, preview);
    }
    else if( ( tmp_vertex = CAST_TO( tmp_ent, RefVertex ) ) != NULL )
    {
     result = scale( tmp_vertex,point, CubitVector(scale_x, scale_y, scale_z),
                check_before_transforming, preview);
    }
    if(result)
      entities_scaled.append( tmp_ent );
  } 
  
}


CubitStatus GeometryQueryTool::scale( Body* body,const CubitVector& point, double factor, bool check_to_transform, bool preview )
{
  if( check_to_transform )
    if (!okay_to_transform( body ))
      return CUBIT_FAILURE;

  CubitTransformMatrix pre_form;
  pre_form.translate( -point );

  CubitTransformMatrix xform;
  xform.scale_about_origin( factor );

  CubitTransformMatrix post_form;
  post_form.translate( point );

 if (preview)
  {
    DLIList<RefEdge*> edges;
    body->ref_edges(edges);
    if( edges.size() )
    {
      for (int i = 0; i < edges.size(); i++)
      {
        GMem poly;
        if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
        {
          poly.transform(pre_form);
          poly.transform(xform);
          poly.transform(post_form);
          GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
        }
        else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
        {
          CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
          tmp_pt = pre_form*tmp_pt;
          tmp_pt = xform*tmp_pt;
          tmp_pt = post_form*tmp_pt;          
          GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
        }
      }
    }
    else
    {
      //just draw the surfaces
      DLIList<RefFace*> faces;
      body->ref_faces( faces );
      for( int i=0; i<faces.size(); i-- )
      {
        GMem poly;
        faces.get_and_step()->get_graphics( poly );
        poly.transform(pre_form);
        poly.transform(xform);
        poly.transform(post_form);

        int* facet_list = poly.facet_list();
        GPoint* plist = poly.point_list();
        
        GPoint p[3];
        for (i = 0; i < poly.fListCount; )
        {
          int sides = facet_list[i++];
          if (sides != 3)
          {
            i += sides;
            continue;
          }
          else
          {
            p[0] = plist[facet_list[i++]];
            p[1] = plist[facet_list[i++]];
            p[2] = plist[facet_list[i++]];
            GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
          }
        }
      }
    }
    GfxPreview::flush();
    return CUBIT_SUCCESS;
  }



  //if (preview)
  //{
  //  DLIList<RefEdge*> edges;
  //  body->ref_edges(edges);
  //  for (int i = 0; i < edges.size(); i++)
  //  {
  //    GMem poly;
  //    edges[i]->get_graphics(poly);
  //    poly.transform(pre_form);
  //    poly.transform(xform);
  //    poly.transform(post_form);
  //    GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
  //  }
  //  GfxPreview::flush();
  //  return CUBIT_SUCCESS;
  //}

  BodySM* bodysm = body->get_body_sm_ptr();
  GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();
  CubitStatus result = engine->translate( bodysm, -point );
  if (result)
  {
    notify_intermediate_of_transform( body, pre_form );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );


  result = engine->scale( bodysm, factor );
  if (result)
  {
    notify_intermediate_of_transform( body, xform );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );

  result = engine->translate( bodysm, point );
  if (result)
  {
    notify_intermediate_of_transform( body, post_form );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );
  
  notify_observers_of_transform( body );

  return result;
}

CubitStatus GeometryQueryTool::scale( Body *body,
                                     const CubitVector& point,
                                      const CubitVector& factors,
                                      bool check_to_transform,
                                      bool preview )
{
  if( check_to_transform )
    if (!okay_to_transform( body ))
      return CUBIT_FAILURE;

  CubitTransformMatrix pre_form;
  pre_form.translate( -point );

  CubitTransformMatrix xform;
  xform.scale_about_origin( factors );

  CubitTransformMatrix post_form;
  post_form.translate( point );

  if (preview)
  {
    DLIList<RefEdge*> edges;
    body->ref_edges(edges);
    if( edges.size() )
    {
      for (int i = 0; i < edges.size(); i++)
      {
        GMem poly;
        if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
        {
          poly.transform(pre_form);
          poly.transform(xform);
          poly.transform(post_form);
          GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
        }
        else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
        {
          CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
          tmp_pt = pre_form*tmp_pt;
          tmp_pt = xform*tmp_pt;
          tmp_pt = post_form*tmp_pt;          
          GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
        }
      }
    }
    else
    {
      //just draw the surfaces
      DLIList<RefFace*> faces;
      body->ref_faces( faces );
      for( int i=0; i<faces.size(); i-- )
      {
        GMem poly;
        faces.get_and_step()->get_graphics( poly );
        poly.transform(pre_form);
        poly.transform(xform);
        poly.transform(post_form);

        int* facet_list = poly.facet_list();
        GPoint* plist = poly.point_list();

        GPoint p[3];
        for (i = 0; i < poly.fListCount; )
        {
          int sides = facet_list[i++];
          if (sides != 3)
          {
            i += sides;
            continue;
          }
          else
          {
            p[0] = plist[facet_list[i++]];
            p[1] = plist[facet_list[i++]];
            p[2] = plist[facet_list[i++]];
            GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
          }
        }
      }
    }

    GfxPreview::flush();
    return CUBIT_SUCCESS;
  }

  BodySM* bodysm = body->get_body_sm_ptr();
  GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();

  CubitStatus result = engine->translate( bodysm, -point );
  if (result)
  {
    notify_intermediate_of_transform( body, pre_form );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );


  result = engine->scale( bodysm, factors );
  if (result)
  {
    notify_intermediate_of_transform( body, xform );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );

  result = engine->translate( bodysm, point );
  if (result)
  {
    notify_intermediate_of_transform( body, post_form );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );
  
  notify_observers_of_transform( body );

  return result;
}

void GeometryQueryTool::reflect( DLIList<RefEntity*> &entities_to_reflect,
                                const CubitVector& point,
                                const CubitVector& normal,
                            bool check_before_transforming,
                            DLIList<RefEntity*> &entities_transformed,
                            bool preview /*= false*/)
{
  //reflect free, merged-away entities
  DLIList<TopologyBridge*> free_ents; 
  get_merged_away_free_entities( entities_to_reflect, free_ents );
  
  int i;
  if (!preview)
  {
    for( i=free_ents.size(); i--; )
    {
      TopologyBridge *bridge = free_ents.get_and_step();
      Curve *curve= CAST_TO( bridge, Curve );
      TBPoint *tmp_point = CAST_TO( bridge, TBPoint );

      if( curve || tmp_point )
      {
        GeometryEntity *geom = CAST_TO( bridge, GeometryEntity );
        GeometryQueryEngine* engine = geom->get_geometry_query_engine();
        CubitStatus result = engine->translate( geom, -point ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
          return;
        }
        result = engine->reflect( geom, normal ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::reflect failed.\n");
          return;
        }
        result = engine->translate( geom, point ); 
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("GeometryQueryEngine::translate failed.\n");
          return;
        }
      }
    }
  }

  RefFace *tmp_face;
  RefEdge *tmp_curve;
  RefVertex *tmp_vertex;
  CubitStatus result;

  DLIList<Body*> bodies_to_reflect;
  DLIList<BasicTopologyEntity*> ents_to_reflect;
  CAST_LIST( entities_to_reflect, bodies_to_reflect, Body);

  if( bodies_to_reflect.size() != entities_to_reflect.size() )
  {
    for( int k=0; k<entities_to_reflect.size(); k++ )
    {
      RefEntity *ent = entities_to_reflect[k];

      if( ( tmp_face = CAST_TO( ent, RefFace ) ) != NULL )
        ents_to_reflect.append( tmp_face );
      else if( (tmp_curve = CAST_TO( ent, RefEdge ) ) != NULL )
        ents_to_reflect.append( tmp_curve );
      else if( (tmp_vertex = CAST_TO( ent, RefVertex ) ) != NULL )
        ents_to_reflect.append( tmp_vertex );
    }
  }

  if( bodies_to_reflect.size() )
  {
    DLIList<Body*> bodies_reflected;
    result = reflect( bodies_to_reflect,
      point, normal, 
      &bodies_reflected,      
      preview );  

    if( result )
    {
      for( int k=0; k<bodies_reflected.size(); k++ )
        entities_transformed.append( bodies_reflected[k] );
    }
  }

  if( ents_to_reflect.size() )
  {
    DLIList<BasicTopologyEntity*> btes_reflected;
     result = reflect( ents_to_reflect, 
      point,
      normal,
      &btes_reflected,
      check_before_transforming,
      preview );  

    if( result )
    {
      for( int k=0; k<btes_reflected.size(); k++ )
        entities_transformed.append( btes_reflected[k] );
    }
  }
}

CubitStatus GeometryQueryTool::reflect( DLIList<Body*> &bodies,
                                        const CubitVector& point,
                                        const CubitVector& normal,
                                        DLIList<Body*> *bodies_reflected,
                                        bool preview)
{
  Body *tmp_body;
  CubitStatus result = CUBIT_FAILURE;

  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix ref_mat;
  ref_mat.reflect(normal );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );

  CubitTransformMatrix total_transform;
  total_transform = mov_mat*ref_mat*prev_xform;
  
  DLIList<RefEntity*> ents_transformed;

  if (preview)
  {
    for (int i = 0; i < bodies.size(); i++)
    {
      DLIList<RefEdge*> edges;
      bodies[i]->ref_edges(edges);
      if( edges.size() )
      {
        for (int i = 0; i < edges.size(); i++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
          {
            poly.transform( total_transform );
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
          {
            CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();            
            tmp_pt = total_transform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
      }
      else
      {
        //just draw the surfaces
        DLIList<RefFace*> faces;
        bodies[i]->ref_faces( faces );
        for( int i=0; i<faces.size(); i-- )
        {
          GMem poly;
          faces.get_and_step()->get_graphics( poly );
          poly.transform( total_transform );

          int* facet_list = poly.facet_list();
          GPoint* plist = poly.point_list();

          GPoint p[3];
          for (i = 0; i < poly.fListCount; )
          {
            int sides = facet_list[i++];
            if (sides != 3)
            {
              i += sides;
              continue;
            }
            else
            {
              p[0] = plist[facet_list[i++]];
              p[1] = plist[facet_list[i++]];
              p[2] = plist[facet_list[i++]];
              GfxPreview::draw_polygon(p, 3, CUBIT_BLUE_INDEX, CUBIT_BLUE_INDEX, false);
            }
          }
        }
      }
      GfxPreview::flush();
    }
    return CUBIT_SUCCESS;
  }
  //first, reflect the underlying geometry
  int i;
  for( i=bodies.size(); i--;)
  {
    tmp_body = bodies.get_and_step();
    BodySM* bodysm = tmp_body->get_body_sm_ptr();
    GeometryQueryEngine* engine = bodysm->get_geometry_query_engine();

    // Move to origin
    result = engine->translate( bodysm, -point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( -point );
      notify_intermediate_of_transform( tmp_body, xform );
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        tmp_body->entity_name().c_str(), tmp_body->class_name(), tmp_body->id() );
      return result;
    }

    result = engine->reflect( bodysm, normal );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.reflect( normal);
      notify_intermediate_of_transform( tmp_body, xform );
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        tmp_body->entity_name().c_str(), tmp_body->class_name(), tmp_body->id() );
      return result;
    }
    
    result = engine->translate( bodysm, point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( point );
      notify_intermediate_of_transform( tmp_body, xform );

      if( bodies_reflected )
        bodies_reflected->append( tmp_body );
      ents_transformed.append( tmp_body );
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        tmp_body->entity_name().c_str(), tmp_body->class_name(), tmp_body->id() );
      return result;
    }
  }
    
  for( i=0; i<ents_transformed.size(); i++)
  {
    tmp_body = (Body*)ents_transformed[i];
    tmp_body = make_Body( tmp_body->get_body_sm_ptr() );
    if( tmp_body != ents_transformed[i] )
      ents_transformed[i] = tmp_body;
  }

  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &total_transform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;    
}

CubitStatus GeometryQueryTool::notify_intermediate_of_transform(
                                        TopologyEntity* entity,
                                        const CubitTransformMatrix& xform ) const
{
  for (IGESet::const_iterator iter = igeSet.begin(); iter != igeSet.end(); ++iter)
    (*iter)->notify_transform( entity->bridge_manager()->topology_bridge(), xform );
  return CUBIT_SUCCESS;
}

void GeometryQueryTool::notify_observers_of_transform(
                                        DLIList<RefEntity*> &ref_ents,
                                        const CubitTransformMatrix *xform ) const
{
  if( xform )
  {
    //gather all the children
    DLIList<RefEntity*> all_ref_ents;
    RefEntity::get_all_child_ref_entities( ref_ents, all_ref_ents );   
    all_ref_ents.uniquify_ordered();

    all_ref_ents += ref_ents;

    for(int i=0; i<all_ref_ents.size(); i++)
    {
      GeometryEvent evt(GeometryEvent::GEOMETRY_TRANSFORMED, all_ref_ents[i]);
      all_ref_ents[i]->notify_observers(&evt);
    }

    AppUtil::instance()->send_event(TransformEvent(*xform, all_ref_ents.as_vector() ));    
  }
  else
  {
    for( int i=0; i<ref_ents.size(); i++ )
    {
      RefEntity *ent = ref_ents[i];
      ent->notify_sub_all_observers( GeometryEvent::GEOMETRY_MODIFIED );    
    }
  }
}

void GeometryQueryTool::notify_observers_of_transform(
                                        RefEntity* ref_entity,
                                        const CubitTransformMatrix *xform ) const
{
  ref_entity->notify_sub_all_observers( GeometryEvent::GEOMETRY_MODIFIED );
  
  //CGMHistory::Event evt(CGMHistory::GEOMETRY_TRANSFORMED, ref_entity);
  //const_cast<CGMHistory&>(mHistory).add_event(evt);
}

CubitBoolean GeometryQueryTool::okay_to_transform( BasicTopologyEntity* bte ) const
{
  if (bte->get_parents() > 0)
  {
    PRINT_ERROR("Cannot transform %s (%s %d) : not a free %s\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id(), bte->class_name() );
    return CUBIT_FALSE;
  }

  DLIList<RefEntity*> list(1);
  list.append( bte );
  if (MergeTool::instance()->contains_merged_entities( list ))
  {
    PRINT_ERROR("Cannot tranmform %s (%s %d) : "
                "%s is merged or contains merged geomtry.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id(), bte->class_name() );
    return CUBIT_FALSE;
  }

  return CUBIT_TRUE;
}

CubitStatus GeometryQueryTool::translate( DLIList<BasicTopologyEntity*> &btes,
                                          const CubitVector& delta,
                                          DLIList<BasicTopologyEntity*> *btes_translated,
                                          bool check_to_transform,
                                          bool preview )
{
  CubitTransformMatrix xform;
  xform.translate(delta);

  DLIList<RefEntity*> transformed_ents;

  for( int k=0; k<btes.size(); k++ )
  {
    BasicTopologyEntity *bte = btes[k];

    if( check_to_transform )
      if (!okay_to_transform( bte ))
        continue;

    if (preview)
    {
      if(bte->dimension()==0)
      {
        DLIList<RefVertex*> points;
        bte->ref_vertices(points);
        for (int i = 0; i < points.size(); i++)
        {
          CubitVector temp(points[i]->center_point());
          temp = xform*temp;
          GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
        }
      }
      else
      {
        DLIList<RefEdge*> edges;
        bte->ref_edges(edges);
        for (int i = 0; i < edges.size(); i++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
          {
            poly.transform(xform);
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
          {
            CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
            tmp_pt = xform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
        if (edges.size() == 0)
        {
          DLIList<RefVertex*> points;
          bte->ref_vertices(points);
          for (int i = 0; i < points.size(); i++)
          {
            CubitVector temp(points[i]->center_point());
            temp=xform*temp;
            GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
          }
        }
      }
      GfxPreview::flush();
      continue;
    }

    GeometryEntity* geom = bte->get_geometry_entity_ptr();
    GeometryQueryEngine* engine = geom->get_geometry_query_engine();
    CubitStatus result = engine->translate( geom, delta );
    if (result)
    {
      notify_intermediate_of_transform( bte, xform );
      
      if( btes_translated )
        btes_translated->append( bte );
      transformed_ents.append( bte );
    }
    else
      PRINT_ERROR("Translate of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );
  }

  if( transformed_ents.size() )
  {
    notify_observers_of_transform( transformed_ents, &xform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryQueryTool::rotate( DLIList<BasicTopologyEntity*> &btes,
                                       const CubitVector& point,
                                       const CubitVector& direction,
                                       double angle,
                                       DLIList<BasicTopologyEntity*> *btes_rotated,
                                       bool check_to_transform,
                                       bool preview )
{
  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix rot_mat;
  rot_mat.rotate( angle, direction );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );
  
  CubitTransformMatrix total_transform;
  total_transform = mov_mat*rot_mat*prev_xform;  

  DLIList<RefEntity*> ents_transformed;

  for( int k=0; k<btes.size(); k++ )
  {
    BasicTopologyEntity *bte = btes[k];

    if( check_to_transform )
      if (!okay_to_transform( bte ))
        continue;

    if (preview)
    {
      if(bte->dimension()==0)
      {
        DLIList<RefVertex*> points;
        bte->ref_vertices(points);
        for (int i = 0; i < points.size(); i++)
        {
          CubitVector temp(points[i]->center_point());        
          temp = total_transform*temp;
          GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
        }
      }
      else
      {

        DLIList<RefEdge*> edges;
        bte->ref_edges(edges);

        for (int j = 0; j < edges.size(); j++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[j]->get_graphics(poly) )
          {          
            poly.transform(total_transform);
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[j]->start_vertex() == edges[j]->end_vertex() )
          {
            CubitVector tmp_pt = edges[j]->start_vertex()->coordinates(); 
            tmp_pt = total_transform*tmp_pt;          
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
        if (edges.size() == 0)
        {
          DLIList<RefVertex*> points;
          bte->ref_vertices(points);
          for (int i = 0; i < points.size(); i++)
          {
            CubitVector p(points[i]->center_point());          
            p = total_transform*p;
            GfxPreview::draw_point(p, CUBIT_BLUE_INDEX);
          }
        }
      }

      GfxPreview::flush();
      continue;
    }

    GeometryEntity* geom = bte->get_geometry_entity_ptr();
    GeometryQueryEngine* engine = geom->get_geometry_query_engine();
    CubitStatus result;
    // Move to origin
    result = engine->translate( geom, -point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( -point );
      notify_intermediate_of_transform( bte, xform );
    }
    else
    {
      PRINT_ERROR("Rotate of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      continue;
    }

    result = engine->rotate( geom, direction, angle );
    if (result)
    {
      CubitTransformMatrix rot_mat;
      rot_mat.rotate( angle, direction );
      notify_intermediate_of_transform( bte, rot_mat );
    }
    else
    {
      PRINT_ERROR("Rotate of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      continue;
    }
    result=engine->translate( geom, point );

    if (result)
    {
      CubitTransformMatrix mov_mat;
      mov_mat.translate( point );
      notify_intermediate_of_transform( bte, mov_mat );

      if( btes_rotated )
        btes_rotated->append( bte );
      ents_transformed.append( bte );
    }
    else
    {
      PRINT_ERROR("Rotate of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      continue;
    }
  }

  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &total_transform );
    return CUBIT_SUCCESS;
  }
   else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryQueryTool::rotate( DLIList<BasicTopologyEntity*> &btes,
                                       const CubitVector& axis,
                                       double angle,
                                       DLIList<BasicTopologyEntity*> *btes_rotated,
                                       bool check_to_transform,
                                       bool preview )
{ 
  CubitTransformMatrix xform;
  xform.rotate( angle, axis );

  DLIList<RefEntity*> ents_transformed;

  for( int k=0; k<btes.size(); k++ )
  {
    BasicTopologyEntity *bte = btes[k];

    if( check_to_transform )
      if (!okay_to_transform( bte ))
        continue;

    if (preview)
    {
      if(bte->dimension()==0)
      {
        DLIList<RefVertex*> points;
        bte->ref_vertices(points);
        for (int i = 0; i < points.size(); i++)
        {
          CubitVector temp(points[i]->center_point());
          temp = xform*temp;
          GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
        }
      }
      else
      {
        DLIList<RefEdge*> edges;
        bte->ref_edges(edges);

        for (int j = 0; j < edges.size(); j++)
        {
          GMem poly;
          if( CUBIT_SUCCESS == edges[j]->get_graphics(poly) )
          {
            poly.transform(xform);
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[j]->start_vertex() == edges[j]->end_vertex() )
          {
            CubitVector tmp_pt = edges[j]->start_vertex()->coordinates();
            tmp_pt = xform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
        if (edges.size() == 0)
        {
          DLIList<RefVertex*> points;
          bte->ref_vertices(points);
          for (int i = 0; i < points.size(); i++)
          {
            CubitVector p(points[i]->center_point());
            CubitVector q =xform*p;
            GfxPreview::draw_point(q, CUBIT_BLUE_INDEX);
          }
        }
      }

      GfxPreview::flush();
      continue;
    }

    GeometryEntity* geom = bte->get_geometry_entity_ptr();
    GeometryQueryEngine* engine = geom->get_geometry_query_engine();
    CubitStatus result = engine->rotate( geom, axis, angle );
    if (result)
    {
      notify_intermediate_of_transform( bte, xform );

      if( btes_rotated )
        btes_rotated->append( bte );
      ents_transformed.append( bte );    
    }
    else
      PRINT_ERROR("Rotate of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );
  }
  
  if( ents_transformed.size() )
  {
    notify_observers_of_transform( ents_transformed, &xform );
    return CUBIT_SUCCESS;
  }
   else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryQueryTool::scale( BasicTopologyEntity* bte,const CubitVector& point, double factor,
                                      bool check_to_transform, bool preview )
{
  if( check_to_transform )
    if (!okay_to_transform( bte ))
      return CUBIT_FAILURE;

  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix xform;
  xform.scale_about_origin( factor );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );


  if (preview)
  {

    if(bte->dimension()==0)
    {
      DLIList<RefVertex*> points;
      bte->ref_vertices(points);
      for (int i = 0; i < points.size(); i++)
      {
        CubitVector temp(points[i]->center_point());
        temp = prev_xform*temp;
        temp = xform*temp;
        temp = mov_mat*temp;
        GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
      }
    }
    else
    {

      DLIList<RefEdge*> edges;
      bte->ref_edges(edges);
      for (int i = 0; i < edges.size(); i++)
      {
        GMem poly;
        if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
        {
          poly.transform(prev_xform);
          poly.transform(xform);
          poly.transform(mov_mat);
          GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
        }
        else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
        {
          CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
          tmp_pt = prev_xform*tmp_pt;
          tmp_pt = xform*tmp_pt;
          tmp_pt = mov_mat*tmp_pt;          
          GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
        }
      }
    }
    GfxPreview::flush();
    return CUBIT_SUCCESS;
  }

  GeometryEntity* geom = bte->get_geometry_entity_ptr();
  GeometryQueryEngine* engine = geom->get_geometry_query_engine();
  CubitStatus result = engine->translate( geom, -point );
  if (result)
  {
    CubitTransformMatrix prev_xform;
    prev_xform.translate(-point);
    notify_intermediate_of_transform( bte, prev_xform );
  }
  else
  {
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );
    return result;
  }

  result = engine->scale( geom, factor );
  if (result)
  {
    CubitTransformMatrix xform;
    xform.scale_about_origin( factor );
    notify_intermediate_of_transform( bte, xform );
  }
  else
  {
    PRINT_ERROR("Scvale of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );
    return result;
  }
  result = engine->translate( geom, point );
  if (result)
  {
    CubitTransformMatrix mov_mat;
    mov_mat.translate( point );

    notify_intermediate_of_transform( bte, prev_xform );
    notify_observers_of_transform( bte );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
    bte->entity_name().c_str(), bte->class_name(), bte->id() );
  return result;
}


CubitStatus GeometryQueryTool::scale( BasicTopologyEntity* bte,
                                     const CubitVector& point,
                                      const CubitVector& factors,
                                      bool check_to_transform,
                                      bool preview )
{
  if( check_to_transform )
    if (!okay_to_transform( bte ))
      return CUBIT_FAILURE;

  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix xform;
  xform.scale_about_origin( factors );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );


  if (preview)
  {
    if(bte->dimension()==0)
    {
      DLIList<RefVertex*> points;
      bte->ref_vertices(points);
      for (int i = 0; i < points.size(); i++)
      {
        CubitVector temp(points[i]->center_point());
        temp = prev_xform*temp;
        temp = xform*temp;
        temp = mov_mat*temp;
        GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
      }
    }
    else
    {

      DLIList<RefEdge*> edges;
      bte->ref_edges(edges);
      for (int i = 0; i < edges.size(); i++)
      {
        GMem poly;  
        if( CUBIT_SUCCESS == edges[i]->get_graphics(poly) )
        {
          poly.transform(prev_xform);
          poly.transform(xform);
          poly.transform(mov_mat);
          GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
        }
        else if( edges[i]->start_vertex() == edges[i]->end_vertex() )
        {
          CubitVector tmp_pt = edges[i]->start_vertex()->coordinates();
          tmp_pt = prev_xform*tmp_pt;
          tmp_pt = xform*tmp_pt;
          tmp_pt = mov_mat*tmp_pt;          
          GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
        }
      }
    }
    GfxPreview::flush();
    return CUBIT_SUCCESS;
  }

  GeometryEntity* geom = bte->get_geometry_entity_ptr();
  GeometryQueryEngine* engine = geom->get_geometry_query_engine();
  CubitStatus result = engine->translate( geom, -point );
  if (result)
  {
    notify_intermediate_of_transform( bte, prev_xform );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );

  result = engine->scale( geom, factors );
  if (result)
  {
    notify_intermediate_of_transform( bte, xform );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );


  result = engine->translate( geom, point );
  if (result)
  {
    notify_intermediate_of_transform( bte, mov_mat );
    notify_observers_of_transform( bte );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      bte->entity_name().c_str(), bte->class_name(), bte->id() );



  return result;
}

CubitStatus GeometryQueryTool::reflect( DLIList<BasicTopologyEntity*> &btes,
                                        const CubitVector& point,
                                        const CubitVector& normal,
                                        DLIList<BasicTopologyEntity*> *btes_reflected,
                                        bool check_to_transform,
                                        bool preview )
{
  CubitTransformMatrix prev_xform;
  prev_xform.translate(-point);

  CubitTransformMatrix ref_mat;
  ref_mat.reflect(normal );

  CubitTransformMatrix mov_mat;
  mov_mat.translate( point );
  
  CubitTransformMatrix total_transform;
  total_transform = mov_mat*ref_mat*prev_xform;
  
  DLIList<RefEntity*> transformed_ents;

  for( int k=0; k<btes.size(); k++ )  
  {
    BasicTopologyEntity *bte = btes[k];

    if( check_to_transform )
      if (!okay_to_transform( bte ))
        return CUBIT_FAILURE;

    if (preview)
    {
      if(bte->dimension()==0)
      {
        DLIList<RefVertex*> points;
        bte->ref_vertices(points);
        for (int i = 0; i < points.size(); i++)
        {
          CubitVector temp(points[i]->center_point());        
          temp = total_transform*temp;
          GfxPreview::draw_point(temp, CUBIT_BLUE_INDEX);
        }
      }
      else
      {

        DLIList<RefEdge*> edges;
        bte->ref_edges(edges);
        for (int j = 0; j < edges.size(); j++)
        {
          GMem poly;        
          if( CUBIT_SUCCESS == edges[j]->get_graphics(poly) )
          {
            poly.transform( total_transform );
            GfxPreview::draw_polyline(poly.point_list(), poly.point_list_size(), CUBIT_BLUE_INDEX);
          }
          else if( edges[j]->start_vertex() == edges[j]->end_vertex() )
          {
            CubitVector tmp_pt = edges[j]->start_vertex()->coordinates();          
            tmp_pt = total_transform*tmp_pt;
            GfxPreview::draw_point( tmp_pt, CUBIT_BLUE_INDEX);
          }
        }
      }
      GfxPreview::flush();
      continue;
    }

    GeometryEntity* geom = bte->get_geometry_entity_ptr();
    GeometryQueryEngine* engine = geom->get_geometry_query_engine();
    CubitStatus result;

    // Move to origin
    result = engine->translate( geom, -point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( -point );
      notify_intermediate_of_transform( bte, xform );
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      return result;
    }

    result = engine->reflect( geom, normal );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.reflect( normal);
      notify_intermediate_of_transform( bte, xform );
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      return result;
    }

    result = engine->translate( geom, point );
    if (result)
    {
      CubitTransformMatrix xform;
      xform.translate( point );
      notify_intermediate_of_transform( bte, xform );

      if( btes_reflected )
        btes_reflected->append( bte );
      transformed_ents.append( bte );      
    }
    else
    {
      PRINT_ERROR("Reflect of %s (%s %d) failed.\n",
        bte->entity_name().c_str(), bte->class_name(), bte->id() );
      return result;
    }
  }

  if( transformed_ents.size() )
  {
    notify_observers_of_transform( transformed_ents, &total_transform );
    return CUBIT_SUCCESS;
  }
  else if( preview )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

void GeometryQueryTool::ige_remove_modified(DLIList<Surface*> &all_surfs,
                                            DLIList<Curve*> &all_curves,
                                            DLIList<TBPoint*> &all_points)
{
  IGESet::iterator itor;

  for (itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->remove_modified(all_surfs, all_curves, all_points);

  for (itor = igeSet.begin(); itor != igeSet.end(); ++itor)
    (*itor)->clean_out_deactivated_geometry();
}

CubitBoolean GeometryQueryTool::bodies_overlap( Body *body_ptr_1,
                                                Body *body_ptr_2 )
{
  BodySM *body1 = body_ptr_1->get_body_sm_ptr();
  BodySM *body2 = body_ptr_2->get_body_sm_ptr();

  if( is_intermediate_geometry( body_ptr_1 ) )
    return body_ptr_1->get_geometry_query_engine()->bodies_overlap( body1, body2 );
  else if( is_intermediate_geometry( body_ptr_2 ) )
    return body_ptr_2->get_geometry_query_engine()->bodies_overlap( body2, body1 );
  else if( body_ptr_1->get_geometry_query_engine() !=
           body_ptr_2->get_geometry_query_engine() )
  {
    PRINT_ERROR("Volumes must be of the same type (ACIS, SolidWorks, etc) to\n"
                "find if they overlap.\n");
    return CUBIT_FALSE;
  }
  else
    return body_ptr_1->get_geometry_query_engine()->bodies_overlap( body1, body2 );
}

CubitBoolean GeometryQueryTool::volumes_overlap( RefVolume *volume_1,
                                                 RefVolume *volume_2 )
{
  Lump *lump1 = volume_1->get_lump_ptr();
  Lump *lump2 = volume_2->get_lump_ptr();

  if( is_intermediate_geometry( volume_1 ) )
    return volume_1->get_geometry_query_engine()->volumes_overlap( lump1, lump2 );
  else if( is_intermediate_geometry( volume_2 ) )
    return volume_2->get_geometry_query_engine()->volumes_overlap( lump2, lump1 );
  else if( volume_1->get_geometry_query_engine() !=
           volume_2->get_geometry_query_engine() )
  {
    PRINT_ERROR("Volumes must be of the same type (ACIS, SolidWorks, etc) to\n"
                "find if they overlap.\n");
    return CUBIT_FALSE;
  }
  else
    return volume_1->get_geometry_query_engine()->volumes_overlap( lump1, lump2 );
}

void GeometryQueryTool::find_nonmanifold_curves(DLIList<RefVolume*> &vol_list,
                                                DLIList<RefEdge*> &curve_list)
{
  int i;
  for(i=vol_list.size(); i>0; i--)
  {
    RefVolume *vol = vol_list.get_and_step();
    DLIList<RefEdge*> vol_curves;
    vol->ref_edges(vol_curves);
    int j;
    for(j=vol_curves.size(); j>0; j--)
    {
      RefEdge *cur_curve = vol_curves.get_and_step();
      if(cur_curve->is_merged())
      {
        DLIList<RefFace*> curve_faces;
        cur_curve->ref_faces(curve_faces);
        bool merged_face_exists = false;
        while(curve_faces.size() && !merged_face_exists)
        {
          RefFace *cur_face = curve_faces.pop();
          if(cur_face->is_merged())
            merged_face_exists = true;
        }
        if(!merged_face_exists)
          curve_list.append(cur_curve);
      }
    }
  }
  curve_list.uniquify_ordered();
}

double GeometryQueryTool::estimate_merge_tolerance(DLIList<RefVolume*> &vol_list,
                                                    bool accurate_in,
                                                    bool report_in,
                                                    double lo_val_in, 
                                                    double hi_val_in,
                                                    int num_calculations_in,
                                                    bool return_calculations_in,
                                                    DLIList<double> *merge_tols,
                                                    DLIList<int> *num_proximities)
{
  double return_merge_tol = -1.0;  // return value of < 0.0 will mean failure
                                   // to find a merge tolerance

  DLIList<double> local_merge_tols;
  DLIList<int> local_num_proximities;
  bool report = report_in;
  double lo_val = 0.0;
  double hi_val = get_geometry_factor()*GEOMETRY_RESABS*10.0;  // 10 * merge tol -- arbitrary
  int num_calculations = num_calculations_in;

  if(lo_val_in != -1.0)
    lo_val = lo_val_in;
  if(hi_val_in != -1.0)
    hi_val = hi_val_in;

  if(hi_val > lo_val)
  {
    double cur_val = lo_val;
    double step = (hi_val - lo_val)/(double)num_calculations;
    int i;
    if(report)
    {
      PRINT_INFO("\n\nLooking for merge toleance...\n\n");
      PRINT_INFO("  Possible range: %f, %f\n", lo_val, hi_val);
      PRINT_INFO("  Number of steps: %d\n\n", num_calculations+1);
    }

    std::map <RefVertex*, DLIList<dist_vert_struct*>*> vert_dist_map;
    GeomMeasureTool::find_near_coincident_vertices_unique(vol_list, hi_val, vert_dist_map);

    int total_num_proximities = 0;
    for(i=0; i<=num_calculations; i++)
    {
      int cur_num_proximities = 0;
      local_merge_tols.append(cur_val);

      std::map <RefVertex*, DLIList<dist_vert_struct*>*>::iterator iter;
      for(iter=vert_dist_map.begin(); iter != vert_dist_map.end(); iter++ )
      {
        //RefVertex *vert = iter->first;
        DLIList<dist_vert_struct*> *struct_list = iter->second;
        int m;
        for(m=struct_list->size(); m>0; m--)
        {
          dist_vert_struct *dvs = struct_list->get();
          if(dvs->dist <= cur_val)
          {
  //          PRINT_INFO("Vertices %d and %d, distance: %.10lf\n", vert->id(), dvs->v2->id(), dvs->dist);
            struct_list->change_to(NULL);
            delete dvs;
            cur_num_proximities++;
          }
          struct_list->step();
        }
        struct_list->remove_all_with_value(NULL);
      }

      total_num_proximities += cur_num_proximities;
      local_num_proximities.append(total_num_proximities);

      if(report)
      {
        PRINT_INFO("  At merge tolerance = %f number of proximities = %d.\n", cur_val, total_num_proximities);
      }

      cur_val += step;
    }

    std::map <RefVertex*, DLIList<dist_vert_struct*>*>::iterator iter;
    for(iter=vert_dist_map.begin(); iter != vert_dist_map.end(); iter++ )
    {
      DLIList<dist_vert_struct*> *struct_list = iter->second;
      // I think all of the items in the lists should be gone
      // by now but just in case...
      while(struct_list->size())
        delete struct_list->pop();
      delete struct_list;
    }

    local_num_proximities.reset();
    local_merge_tols.reset();
    int num_total = local_merge_tols.size();
    if(num_total > 2)
    {
      int num_triplets = num_total - 2;
      int h, min_index, min_diff;
      DLIList<int> diffs;
      for(h=0; h<num_triplets; h++)
      {
        int num_begin = local_num_proximities.get();
        local_num_proximities.step(2);
        int num_end = local_num_proximities.get();
        local_num_proximities.back();
        int cur_diff = num_end - num_begin;
        if(h==0)
        {
          min_index = h;
          min_diff = cur_diff;
        }
        else
        {
          if(cur_diff < min_diff)
          {
            min_diff = cur_diff;
            min_index = h;
          }
        }
      }
      local_merge_tols.step(min_index+1);
      return_merge_tol = local_merge_tols.get();
    }
    else
      PRINT_ERROR("Unable to estimate merge tolerance.\n");

/*
    // Pick off a merge tolerance.
    local_num_proximities.reset();
    local_merge_tols.reset();
    DLIList<int> unique_num_proximities;
    DLIList<int> unique_counts;
    DLIList<double> tmp_merge_tol_list;
    int tmp_num_proximities;
    double tmp_merge_tol;

    double cur_merge_tol = local_merge_tols.get_and_step();
    int cur_num_prox = local_num_proximities.get_and_step();
    int cur_unique_counts = 1;
    // Loop over the whole size even though we have processed the 
    // first entry because we need to record the results after the
    // last entry.
    for(i=local_num_proximities.size(); i>0; i--)
    {
      if(i>1)
      {
        tmp_num_proximities = local_num_proximities.get_and_step();
        tmp_merge_tol = local_merge_tols.get_and_step();
      }
      else 
      {
        // On the last time in just give it a dummy value so we
        // can record the results from the last real entry.
        tmp_num_proximities = -1;
      }
      if(cur_num_prox == tmp_num_proximities)
      {
        cur_unique_counts++;
      }
      else
      {
        tmp_merge_tol_list.append(cur_merge_tol);
        unique_counts.append(cur_unique_counts);
        unique_num_proximities.append(cur_num_prox);
        cur_unique_counts = 1;
        cur_num_prox = tmp_num_proximities;
        cur_merge_tol = tmp_merge_tol;
      }
    }

    int max_index = -1;
    int cur_max_num_counts = 0;
    unique_counts.reset();
    unique_num_proximities.reset();
    tmp_merge_tol_list.reset();
    for(i=unique_counts.size(); i>0; i--)
    {
      int cur_num_counts = unique_counts.get();
      if(cur_num_counts > cur_max_num_counts)
      {
        cur_max_num_counts = cur_num_counts;
        max_index = unique_counts.get_index();
      }
      unique_counts.step();
    }

    if(max_index > -1)
    {
      tmp_merge_tol_list.step(max_index);
      return_merge_tol = tmp_merge_tol_list.get();
    }
    else
    {
      PRINT_ERROR("Unable to estimate merge tolerance.\n");
    }
    */

    if(report)
      PRINT_INFO("\nEstimated merge tolerance: %f\n", return_merge_tol);
  }
  else
    PRINT_ERROR("Range low value is larger than range high value.\n");

  return return_merge_tol;
}

void GeometryQueryTool::find_floating_volumes(DLIList<RefVolume*> &vol_list,
                                              DLIList<RefVolume*> &floating_list)
{
  int i;
  for(i=vol_list.size(); i>0; i--)
  {
    bool floating = true;
    RefVolume *vol = vol_list.get_and_step();
    DLIList<RefEdge*> vol_curves;
    DLIList<RefVertex*> vol_verts;
    DLIList<RefFace*> vol_surfs;
    vol->ref_edges(vol_curves);
    vol->ref_faces(vol_surfs);
    vol->ref_vertices(vol_verts);
    int j;
    for(j=vol_surfs.size(); j>0 && floating; j--)
    {
      RefFace *cur_surf = vol_surfs.get_and_step();
      if(cur_surf->is_merged())
        floating = false;
    }
    for(j=vol_curves.size(); j>0 && floating; j--)
    {
      RefEdge *cur_curve = vol_curves.get_and_step();
      if(cur_curve->is_merged())
        floating = false;
    }
    for(j=vol_verts.size(); j>0 && floating; j--)
    {
      RefVertex *cur_vert = vol_verts.get_and_step();
      if(cur_vert->is_merged())
        floating = false;
    }
    if(floating)
      floating_list.append(vol);
  }
  floating_list.uniquify_ordered();
}

void GeometryQueryTool::find_nonmanifold_vertices(DLIList<RefVolume*> &vol_list,
                                                DLIList<RefVertex*> &vertex_list)
{
  int i;
  for(i=vol_list.size(); i>0; i--)
  {
    RefVolume *vol = vol_list.get_and_step();
    DLIList<RefVertex*> vol_verts;
    vol->ref_vertices(vol_verts);
    int j;
    for(j=vol_verts.size(); j>0; j--)
    {
      RefVertex *cur_vert = vol_verts.get_and_step();
      if(cur_vert->is_merged())
      {
        DLIList<RefEdge*> vert_edges;
        cur_vert->ref_edges(vert_edges);
        bool merged_edge_exists = false;
        while(vert_edges.size() && !merged_edge_exists)
        {
          RefEdge *cur_edge = vert_edges.pop();
          if(cur_edge->is_merged())
            merged_edge_exists = true;
        }
        if(!merged_edge_exists)
          vertex_list.append(cur_vert);
      }
    }
  }
  vertex_list.uniquify_ordered();
}

CGMHistory& GeometryQueryTool::history()
{
  return mHistory;
}

void GeometryQueryTool::get_merged_away_free_entities( DLIList<RefEntity*> &ref_ents,
                                                       DLIList<TopologyBridge*> &free_ents )
{
  //determine if you have any free entities that were merged away 
  DLIList<RefEntity*> merged_ref_ents;
  MergeTool::instance()->contains_merged_entities( ref_ents, &merged_ref_ents );
  int i;
  for( i=merged_ref_ents.size(); i--; )
  {
    RefEntity *tmp_ent = merged_ref_ents.get_and_step();
    DLIList<TopologyBridge*> bridge_list;
    TopologyEntity *te = CAST_TO(tmp_ent, TopologyEntity );
    te->bridge_manager()->get_bridge_list( bridge_list );
    //bridge_list.reset();
    //bridge_list.step(); //don't want to first bridge...that's the real thing
    
    //check for free entities...
    int j;
    for( j=0; j<bridge_list.size(); j++ )
    {
      TopologyBridge *tmp_bridge = bridge_list.get_and_step();

      //ignore the representation bridge if it's a free vertex
      if( j==0 )
        if( tmp_ent->num_parent_ref_entities() == 0 )
          continue;

      DLIList<TopologyBridge*> parents;
      tmp_bridge->get_parents( parents );
      if( parents.size() == 0 ) 
        free_ents.append( tmp_bridge );
    }
  }
}


CubitStatus GeometryQueryTool::get_graphics( RefFace *ref_face,
                                             GMem *gmem,
                                             std::vector<RefEntity*> &facet_point_ownership_vector,
                                             std::vector<std::pair< RefEntity*, std::pair<int,int> > > &facetedges_on_refedges,
                                             unsigned short normal_tolerance, 
                                             double distance_tolerance, 
                                             double max_edge_length )
{
  Surface* surf_ptr = ref_face->get_surface_ptr();
  if (!surf_ptr)
  {
    PRINT_ERROR("RefFace %d is invalid -- no attached Surface.\n",ref_face->id());
    return CUBIT_FAILURE;
  }
   
  std::vector<TopologyBridge*> vertex_edge_to_point_vector;
  std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > facetedges_on_curve;

  CubitStatus rv = surf_ptr->get_geometry_query_engine()->
    get_graphics(surf_ptr, gmem, vertex_edge_to_point_vector, facetedges_on_curve,
    normal_tolerance, distance_tolerance, max_edge_length );
  if (CUBIT_SUCCESS != rv) {
    PRINT_ERROR("GeometryQueryEngine::get_graphics failed.\n");
    return rv;
  }

  facet_point_ownership_vector.resize( vertex_edge_to_point_vector.size(), NULL );
  RefEntity *ref_ent = NULL;
  for( unsigned int i = 0; i < vertex_edge_to_point_vector.size(); i++ )
  {    
    TopologyBridge *tb = vertex_edge_to_point_vector[i];     
    if( NULL == tb )
      continue;
    
    if( tb->topology_entity() )
    {
      ref_ent = CAST_TO( tb->topology_entity(), RefEntity ); 
      facet_point_ownership_vector[i] = ref_ent;
    }
    else
      assert(0);
  }

  for( unsigned int i = 0; i < facetedges_on_curve.size(); i++ )
  {    
    std::pair<TopologyBridge*, std::pair<int,int> > tmp_pair = facetedges_on_curve[i];
    TopologyBridge *tb = tmp_pair.first;
    
    if( tb->topology_entity() )
    {
      ref_ent = CAST_TO( tb->topology_entity(), RefEntity ); 
      facetedges_on_refedges.push_back( std::make_pair( ref_ent, tmp_pair.second ) ); 
    }
    else
    {
      assert(0);
    }
  }

  
  bool debug = false;

  if( debug )
  {    
    //GfxDebug::clear();
    //GPoint *pt_list = gmem->point_list();
    for( unsigned int i = 0; i < facet_point_ownership_vector.size(); i++ )
    {
      //RefEntity *ref_ent = facet_point_ownership_vector[i];

      //int color = 3;
      //if( ref_ent->dimension() == 1 )
        //color = 4;
      //else if( ref_ent->dimension() == 2 )
        //color = 5;
    }   

    GfxDebug::flush();
    GfxDebug::mouse_xforms();    
    

    GfxDebug::clear();
    //draw the curves
    for( unsigned int i=0; i<facetedges_on_refedges.size(); i++ )
    {
      std::pair<RefEntity*, std::pair<int,int> > tmp_pair;
      tmp_pair = facetedges_on_refedges[i];

      RefEntity* ref_ent = tmp_pair.first;
      std::pair<int,int> int_pair = tmp_pair.second;

      CubitVector pt0( gmem->point_list()[int_pair.first].x,
        gmem->point_list()[int_pair.first].y,
        gmem->point_list()[int_pair.first].z );

      CubitVector pt1( gmem->point_list()[int_pair.second].x,
        gmem->point_list()[int_pair.second].y,
        gmem->point_list()[int_pair.second].z );

      GfxDebug::draw_point(pt0, (ref_ent->id()%10) + 3  );
      GfxDebug::draw_point(pt1, (ref_ent->id()%10) + 3  );
      GfxDebug::draw_line( pt0, pt1, (ref_ent->id()%10) + 3 );      
    }
    GfxDebug::flush();      
    GfxDebug::mouse_xforms();

    int index = 0;
    GfxDebug::clear();
    for( int i=0; i<gmem->fListCount; i++ )
    {
      int step = gmem->facet_list()[index++];      
      //create pts      
      CubitVector pt0( gmem->point_list()[ gmem->facet_list()[index] ].x,
        gmem->point_list()[ gmem->facet_list()[index] ].y,
        gmem->point_list()[ gmem->facet_list()[index] ].z );       
      index++;

      CubitVector pt1( gmem->point_list()[ gmem->facet_list()[index] ].x,
        gmem->point_list()[ gmem->facet_list()[index] ].y,
        gmem->point_list()[ gmem->facet_list()[index] ].z );
      index++;

      CubitVector pt2( gmem->point_list()[ gmem->facet_list()[index] ].x,
        gmem->point_list()[ gmem->facet_list()[index] ].y,
        gmem->point_list()[ gmem->facet_list()[index] ].z );       
      index++;

       //draw lines
       GPoint three_pts[3];
       three_pts[0].x = pt0.x();
       three_pts[0].y = pt0.y();
       three_pts[0].z = pt0.z();    

       three_pts[1].x = pt1.x();
       three_pts[1].y = pt1.y();
       three_pts[1].z = pt1.z();    

       three_pts[2].x = pt2.x();
       three_pts[2].y = pt2.y();
       three_pts[2].z = pt2.z();                  

       GfxDebug::draw_polygon( three_pts, 3, i%20, i%20, true );
       GfxDebug::draw_line( pt0, pt1, i%20 );
       GfxDebug::draw_line( pt1, pt2, i%20 );
       GfxDebug::draw_line( pt0, pt2, i%20 );       
       i+=step;              
    }          
    GfxDebug::flush();
    GfxDebug::mouse_xforms();
  }




  return CUBIT_SUCCESS;  
}


CubitStatus GeometryQueryTool::get_graphics( Body *body, 
                                             GMem *g_mem,
                                             std::vector<RefFace*> &ref_face_to_facet_vector,
                                             std::vector<RefEntity*> &facet_point_ownership_vector,
                                             std::vector<std::pair< RefEntity*, std::pair<int,int> > > &facetedges_on_refedges,
                                             unsigned short normal_tolerance, 
                                             double distance_tolerance, 
                                             double max_edge_length )
{
  //get all the RefFaces from 
  DLIList<RefFace*> all_ref_faces; 
  body->ref_faces( all_ref_faces );    

  BodySM *body_sm = NULL;

  //get the underlying body in partition bodies
  if( GeometryQueryTool::instance()->is_intermediate_geometry( body ) )
  {      
    DLIList<TopologyBridge*> partition_body(1);
    body->get_body_sm_ptr()->get_geometry_query_engine()->get_underlying_bridges( body->get_body_sm_ptr(), partition_body );
    if( partition_body.size() == 1 )
      body_sm = static_cast<BodySM*>(partition_body.get());
  }
  else  
    body_sm = body->get_body_sm_ptr();        

  //make sure everything belongs to the same GQE  
  GeometryQueryEngine* gqe = body_sm->get_geometry_query_engine();  
  
  //get the graphics
  std::vector<Surface*> surfaces_to_facet_vector;
  std::vector<TopologyBridge*> tmp_facet_point_ownership_vector;
  std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > facetedges_on_curve;
  CubitStatus stat = gqe->get_graphics( body_sm, g_mem, surfaces_to_facet_vector, tmp_facet_point_ownership_vector,
    facetedges_on_curve, normal_tolerance, distance_tolerance, max_edge_length );

  if( CUBIT_FAILURE == stat )
    return CUBIT_FAILURE;

  //map Surfaces to MRefFaces
  Surface *current_surf = NULL;
  RefFace *current_ref_face = NULL;
  bool ignore_facets = false;

  unsigned int i;
  for( i=0; i<surfaces_to_facet_vector.size(); i++ )
  {    
    Surface *tmp_surf = surfaces_to_facet_vector[i];     
    
    if( tmp_surf != current_surf )
    {
      ignore_facets = false;
      current_surf = tmp_surf;

      if( tmp_surf->topology_entity() )
      {
        current_ref_face = CAST_TO( tmp_surf->topology_entity(), RefFace ); 
        if( current_ref_face )        
          ref_face_to_facet_vector.push_back( current_ref_face );
      }
      else
      {
        int times_through = 0;
        for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
        {          
          times_through++;
          DLIList<TopologyBridge*> bridge_list;
          (*itor)->get_tbs_with_bridge_manager_as_owner(tmp_surf, bridge_list);
          if( bridge_list.size() == 0 )
          {
            if( times_through == 2 )
              assert(0);
            continue;
          }

          if( bridge_list.size() == 1 )
          {
            current_ref_face = CAST_TO( bridge_list.get()->topology_entity(), RefFace );
            if( current_ref_face )
              ref_face_to_facet_vector.push_back( current_ref_face );              
            else 
            {
              PRINT_INFO("Having trouble finding the top-level RefFace\n");
              assert(0);
            }  
            break;
          }
          else if( bridge_list.size() > 1 ) // it is a partition...we'll ignore the facets
          {             
            ref_face_to_facet_vector.push_back( NULL );
            ignore_facets = true;
            break;
          }
        }     
      }
    }
    else
    {   
      if( ignore_facets )
        ref_face_to_facet_vector.push_back( NULL );
      else
        ref_face_to_facet_vector.push_back( current_ref_face );
    }
  }
  
  facet_point_ownership_vector.resize( tmp_facet_point_ownership_vector.size(), NULL );
  RefEntity *ref_ent = NULL;
  for( i=0; i<tmp_facet_point_ownership_vector.size(); i++ )
  {    
    TopologyBridge *tb = tmp_facet_point_ownership_vector[i];     
    if( NULL == tb )
      continue;
    
    if( tb->topology_entity() )
    {
      ref_ent = CAST_TO( tb->topology_entity(), RefEntity ); 
      facet_point_ownership_vector[i] = ref_ent;
    }
    else
    {
      int times_through = 0;
      for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
      {          
        times_through++;
        DLIList<TopologyBridge*> bridge_list;
        (*itor)->get_tbs_with_bridge_manager_as_owner(tb, bridge_list);
        if( bridge_list.size() == 0 )
        {
          if( times_through == 2 )
            assert(0);
          continue;
        }

        //composites
        if( bridge_list.size() == 1 )
        {
          ref_ent = CAST_TO( bridge_list.get()->topology_entity(), RefEntity );
          if( ref_ent )
            facet_point_ownership_vector[i] = ref_ent;
          else 
          {
            PRINT_INFO("Having trouble finding the top-level RefEdge or RefVertex\n");
            assert(0);
          }  
          break;
        }
        else if( bridge_list.size() > 1 ) // it is a partition...we'll ignore the facets
        {             
          facet_point_ownership_vector[i] = NULL;
          ignore_facets = true;
          break;
        }
      }     
    }
  }

  for( i=0; i<facetedges_on_curve.size(); i++ )
  {    
    std::pair<TopologyBridge*, std::pair<int,int> > tmp_pair = facetedges_on_curve[i];
    TopologyBridge *tb = tmp_pair.first;
    
    if( tb->topology_entity() )
    {
      ref_ent = CAST_TO( tb->topology_entity(), RefEntity ); 
      facetedges_on_refedges.push_back( std::make_pair( ref_ent, tmp_pair.second ) ); 
    }
    else
    {
      int times_through = 0;
      for (IGESet::iterator itor = igeSet.begin(); itor != igeSet.end(); ++itor)
      {          
        times_through++;
        DLIList<TopologyBridge*> bridge_list;
        (*itor)->get_tbs_with_bridge_manager_as_owner(tb, bridge_list);
        if( bridge_list.size() == 0 )
        {
          //if( times_through == 2 )
          //  assert(0);
          continue;
        }

        //composites
        if( bridge_list.size() == 1 )
        {
          ref_ent = CAST_TO( bridge_list.get()->topology_entity(), RefEntity );
          if( ref_ent )
            facetedges_on_refedges.push_back( std::make_pair( ref_ent, tmp_pair.second ) ); 
          else 
          {
            PRINT_INFO("Having trouble finding the top-level RefEdge or RefVertex\n");
            assert(0);
          }  
          break;
        }
        else if( bridge_list.size() > 1 ) // it is a partition...we'll ignore the facets
        {                  
          break;
        }
      }     
    }
  }


  bool debug = false;

  if( debug )
  {    
    //GfxDebug::clear();
    GPoint *pt_list = g_mem->point_list();
    for( unsigned int i = 0; i < facet_point_ownership_vector.size(); i++ )
    {
      RefEntity *ref_ent = facet_point_ownership_vector[i];      

      int color = 3;
      if( ref_ent->dimension() == 1 )       
        color = 4;
      else if( ref_ent->dimension() == 2 )       
        color = 5;   

      GfxDebug::draw_point( pt_list[i].x, pt_list[i].y, pt_list[i].z, color );
    }   

    GfxDebug::flush();
    GfxDebug::mouse_xforms();        
/*
    GfxDebug::clear();
    //draw the curves
    for( int i=0; i<facetedges_on_refedges.size(); i++ )
    {
      std::pair<RefEntity*, std::pair<int,int> > tmp_pair;
      tmp_pair = facetedges_on_refedges[i];

      RefEntity* ref_ent = tmp_pair.first;
      std::pair<int,int> int_pair = tmp_pair.second;

      CubitVector pt0( g_mem->point_list()[int_pair.first].x,
        g_mem->point_list()[int_pair.first].y,
        g_mem->point_list()[int_pair.first].z );

      CubitVector pt1( g_mem->point_list()[int_pair.second].x,
        g_mem->point_list()[int_pair.second].y,
        g_mem->point_list()[int_pair.second].z );

      GfxDebug::draw_point(pt0, (ref_ent->id()%10) + 3  );
      GfxDebug::draw_point(pt1, (ref_ent->id()%10) + 3  );
      GfxDebug::draw_line( pt0, pt1, (ref_ent->id()%10) + 3 );      
    }
    GfxDebug::flush();      
    GfxDebug::mouse_xforms();

    int index = 0;
    int color;

    index = 0;
    GfxDebug::clear();
    for( int i=0; i<g_mem->fListCount; i++ )
    {
      int step = g_mem->facet_list()[index++];      
      //create pts
      RefEntity *ref_ent = facet_point_ownership_vector[ g_mem->facet_list()[index] ];
      if( ref_ent->dimension() == 2 )
        color = ref_ent->id() + 3 ;
       CubitVector pt0( g_mem->point_list()[ g_mem->facet_list()[index] ].x,
         g_mem->point_list()[ g_mem->facet_list()[index] ].y,
         g_mem->point_list()[ g_mem->facet_list()[index] ].z );       
       index++;

       ref_ent = facet_point_ownership_vector[ g_mem->facet_list()[index] ];
       if( ref_ent->dimension() == 2 )
         color = ref_ent->id() + 3 ;       
       CubitVector pt1( g_mem->point_list()[ g_mem->facet_list()[index] ].x,
         g_mem->point_list()[ g_mem->facet_list()[index] ].y,
         g_mem->point_list()[ g_mem->facet_list()[index] ].z );
       index++;

       ref_ent = facet_point_ownership_vector[ g_mem->facet_list()[index] ];
       if( ref_ent->dimension() == 2 )
         color = ref_ent->id() + 3 ;       
       CubitVector pt2( g_mem->point_list()[ g_mem->facet_list()[index] ].x,
         g_mem->point_list()[ g_mem->facet_list()[index] ].y,
         g_mem->point_list()[ g_mem->facet_list()[index] ].z );       
       index++;

       //draw lines

       GPoint three_pts[3];
       three_pts[0].x = pt0.x();
       three_pts[0].y = pt0.y();
       three_pts[0].z = pt0.z();    

       three_pts[1].x = pt1.x();
       three_pts[1].y = pt1.y();
       three_pts[1].z = pt1.z();    

       three_pts[2].x = pt2.x();
       three_pts[2].y = pt2.y();
       three_pts[2].z = pt2.z();                  

       GfxDebug::draw_polygon( three_pts, 3, i%20, i%20, true );
       GfxDebug::draw_line( pt0, pt1, color );
       GfxDebug::draw_line( pt1, pt2, color );
       GfxDebug::draw_line( pt0, pt2, color );       
       i+=step;              
    }          
    GfxDebug::flush();
    GfxDebug::mouse_xforms(); 
    */
  }

  return CUBIT_SUCCESS;
}


CubitStatus GeometryQueryTool::get_point_containment( DLIList<Body*> &body_list,
                                                      DLIList<CubitVector> &point_list,
                                                      double tolerance,
                                                      bool allow_pts_in_multiple_bodies,
                                                      std::vector< std::pair<Body*, std::vector<int> > > &body_to_pt_indices )
{
  DLIList<TopologyEntity*> entity_list( body_list.size() );
  for( int i=body_list.size(); i--; )
    entity_list.append( body_list.get_and_step() );
    
  DLIList<TopologyBridge*> bridge_list;
  GeometryQueryEngine* gqe = common_query_engine( entity_list, bridge_list );

  //map for associating PartitionBody back to Body
  std::map<BodySM*, Body*> partition_body_map;
  
  //there is a slim chance that all bodies are ParititonBodies.  If so,
  //the gqe will be non-NULL...let's check if the first one is
  if( gqe && GeometryModifyTool::instance()->contains_partitions( body_list ) )
    gqe = NULL;
  
  //if the engine is still NULL, we could just have a partition body in there...
  if( gqe == NULL )
  {
    bridge_list.clean_out();

    for( int k=0; k<body_list.size(); k++ )
    {
      DLIList<TopologyBridge*> underlying_body;
      
      //if it is virtual
      if( body_list[k]->get_geometry_query_engine()->is_intermediate_engine() )
      {
        //get the underlying Body...in PartitionBody case, it should be a BodyACIS
        underlying_body.append( body_list[k]->bridge_manager()->topology_bridge() );
        IGESet::reverse_iterator itor;
        for (itor = igeSet.rbegin(); itor != igeSet.rend(); ++itor)
          (*itor)->export_geometry(underlying_body);
        
        //if it is different than the original Body TopologyBridge, we found it
        if( underlying_body.size() && body_list[k]->bridge_manager()->topology_bridge() != underlying_body.get() )
        {
          partition_body_map.insert( std::make_pair( (BodySM*)underlying_body.get(), body_list[k] ) );
          bridge_list.append( underlying_body.get() );
        }
      }
      else
        bridge_list.append( body_list[k]->bridge_manager()->topology_bridge() );
    }

    bridge_list.reset();
    gqe = bridge_list.get_and_step()->get_geometry_query_engine();

    for( int i=1; i<bridge_list.size(); i++ )
    {
      TopologyBridge *tmp_bridge = bridge_list.get_and_step();
      if( gqe != tmp_bridge->get_geometry_query_engine() )
      {
        gqe = NULL;
        break;
      }
    }
  }

  if( gqe == NULL )
  {
    PRINT_ERROR("Input geometry does not have the same geometry engine.\n");
    return CUBIT_FAILURE;
  }

  DLIList<BodySM*> body_sm_list;
  CAST_LIST(bridge_list, body_sm_list, BodySM);

  std::vector< std::pair<BodySM*, std::vector<int> > > bodysm_to_pt_indices;
  gqe->get_point_containment( body_sm_list, point_list, tolerance, allow_pts_in_multiple_bodies, bodysm_to_pt_indices );

  //convert the BodySMs back to BODYs  
  for( unsigned int i=0; i<bodysm_to_pt_indices.size(); i++ )
  {
    std::pair<BodySM*, std::vector<int> > tmp_pair = bodysm_to_pt_indices[i];    
    BodySM *body_sm = tmp_pair.first;    
    Body *tmp_body = dynamic_cast<Body*>( body_sm->topology_entity() );
    if( !tmp_body )
    {
      //PartitionBody case...look it up in the map
      std::map<BodySM*, Body*>::iterator tmp_iter = partition_body_map.find( body_sm );
      if( partition_body_map.end() != tmp_iter  )
      {
        body_to_pt_indices.push_back( std::make_pair( tmp_iter->second,  tmp_pair.second) );
      }
    }
    else
      body_to_pt_indices.push_back( std::make_pair( tmp_body, tmp_pair.second) );
  }

  return CUBIT_SUCCESS;
}

void GeometryQueryTool::validate_geometry_database()
{
  DLIList<GeometryQueryEngine*> gqe_list;
  this->get_gqe_list(gqe_list);
  for(int i=0; i<gqe_list.size(); i++)
  {
    GeometryQueryEngine *gqe = gqe_list[i];
    gqe->validate_geometry_database();
  }
}


