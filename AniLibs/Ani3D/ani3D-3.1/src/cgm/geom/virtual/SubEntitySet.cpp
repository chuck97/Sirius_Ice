#include "VGDefines.h"
#include "SubEntitySet.hpp"
#include "PartitionEntity.hpp"
#include "CoEdgeSM.hpp"
#include "DLIList.hpp"
#include "CubitSimpleAttrib.hpp"
#include "TDUniqueId.hpp"
#include "CubitVector.hpp"
#include "PartitionEngine.hpp"
#include "GeometryEntity.hpp"
#include "PartitionBody.hpp"
#include "PartitionSurface.hpp"
#include "CubitMessage.hpp"
#include <typeinfo>


const char* const PARTITION_DATA_ATTRIB_NAME = "PARTITION_ATTRIB";
const char* const PARTITION_GEOM_ATTRIB_NAME = "PARTITION_GEOM";
const char* const PARTITION_UID_ATTRIB_NAME  = "PARTITION_UID";

void SubEntitySet::strip_attributes( TopologyBridge* bridge )
{
  DLIList<CubitSimpleAttrib> list;

  bridge->get_simple_attribute( PARTITION_DATA_ATTRIB_NAME, list );
  while ( list.size() )
  {
    CubitSimpleAttrib csa = list.pop();
    bridge->remove_simple_attribute_virt( csa );
  }
 
  bridge->get_simple_attribute( PARTITION_GEOM_ATTRIB_NAME, list );
  while ( list.size() )
  {
    CubitSimpleAttrib csa = list.pop();
    bridge->remove_simple_attribute_virt( csa );
  }
 
  bridge->get_simple_attribute( PARTITION_UID_ATTRIB_NAME, list );
  while ( list.size() )
  {
    CubitSimpleAttrib csa = list.pop();
    bridge->remove_simple_attribute_virt( csa );
  }
}
  

SubEntitySet::SubEntitySet( TopologyBridge* real_entity, 
                            PartitionEntity* first_part )
  : bodyNext(0), bodyPtr(0), subEntityHead(0), lowerOrderHead(0), lastId(0), 
    layerNumber(SUBCOMP_PARTITION_LAYER), uniqueId(0)
{
//  if( real_entity->owner() )
//    real_entity->owner()->swap_bridge( real_entity, 
//        dynamic_cast<TopologyBridge*>(first_part) );
  
//  real_entity->owner( this );
  myEntity = real_entity;
  assert(!!myEntity);
  
    // Scan for existing partition geometry stored
    // in attributes on the real entity.  We need
    // to start new Ids after any existing ones
    // to avoid conflicts as the stored entities
    // are restored.
  DLIList<CubitSimpleAttrib> csa_list;
  
  real_entity->get_simple_attribute( PARTITION_GEOM_ATTRIB_NAME, csa_list );
  while( csa_list.size() )
  {
    CubitSimpleAttrib csa = csa_list.pop();
    if( csa.character_type() == PARTITION_GEOM_ATTRIB_NAME )
    {
      int id = csa.int_data_list()[0];
      if( id > lastId )
        lastId = id;
    }
  }
  
  real_entity->get_simple_attribute( PARTITION_UID_ATTRIB_NAME, csa_list );
  if ( csa_list.size() )
  {
    assert(csa_list.size() == 1);
    CubitSimpleAttrib csa = csa_list.pop();
    //real_entity->remove_simple_attribute_virt(csa);
    assert(csa.int_data_list().size() == 1);
    uniqueId = csa.int_data_list()[0];
    PartitionEngine::instance().add_to_id_map( this, uniqueId );
  }
  
  // add first partition into the SubEntitySet's linked list
  add_partition( first_part, 0 );
  
  csa_list.clean_out();
  real_entity->get_simple_attribute( PARTITION_DATA_ATTRIB_NAME, csa_list );
  while ( csa_list.size() )
  {
    CubitSimpleAttrib csa = csa_list.pop();
    real_entity->remove_simple_attribute_virt(csa);
    wrap_attribute(csa, first_part->entitySetId );
    real_entity->append_simple_attribute_virt(csa);
  }
}

void SubEntitySet::unwrap_attributes()
{
  DLIList<CubitSimpleAttrib> csa_list;
  if ( myEntity) 
    myEntity->get_simple_attribute( PARTITION_DATA_ATTRIB_NAME, csa_list );
  while( csa_list.size() )
  {
    CubitSimpleAttrib csa = csa_list.pop();
    myEntity->remove_simple_attribute_virt(csa);
    unwrap_attribute(csa);
    myEntity->append_simple_attribute_virt(csa);
  }
}
  

int SubEntitySet::get_unique_id()
{
  if( uniqueId == 0 )
  {
    uniqueId = TDUniqueId::generate_unique_id();
    CubitSimpleAttrib attrib( PARTITION_UID_ATTRIB_NAME, 0, 0, &uniqueId );
    myEntity->append_simple_attribute_virt( attrib );
    PartitionEngine::instance().add_to_id_map( this, uniqueId );
  }
  return uniqueId;
}

SubEntitySet::~SubEntitySet()
{
  if( uniqueId )
    PartitionEngine::instance().remove_from_id_map(this,uniqueId);
  assert( subEntityHead == 0 && lowerOrderHead == 0 );
  if( myEntity ) {
    strip_attributes(myEntity);
    if( myEntity->owner() == this )
      myEntity->owner(0);
  }
    
  if( bodyPtr )
  {
      // keep copy because calling bodyPtr->remove
      // means that bodyPtr will become NULL
    PartitionBody* body = bodyPtr;
    body->remove(*this);

    assert(body->entitySet != this);
    if( !body->has_children() )
    {
      BodySM* real_body = body->real_body();
      body->sub_entity_set().remove_bridge(real_body);
      if( body->owner() )
        body->owner()->swap_bridge( body, real_body, false );
      delete body;
    }
  }
}

void SubEntitySet::add_lower_order( PartitionEntity* partition )
{
  assert( partition->entitySet == 0 &&
          partition->entitySetNext == 0 );
  partition->entitySet = this;
  partition->entitySetNext = lowerOrderHead;
  lowerOrderHead = partition;
  partition->entitySetId = ++lastId;
}

void SubEntitySet::add_lower_order( PartitionEntity* partition,
                                    const CubitSimpleAttrib& attrib,
                                    int dimension,
                                    DLIList<CubitVector*>& points,
                                    DLIList<int>& facets,
                                    DLIList<int>& children,
                                    DLIList<int>& point_owners )
{
  assert( partition->entitySet == 0 &&
          partition->entitySetNext == 0 );
          
  int id, dim;
  CubitStatus result = read_geometry( id, dim, points, facets, children, 
                                      point_owners, attrib );
  assert( result && dimension == dim );
  if (CUBIT_SUCCESS != result || dimension != dim) {
    PRINT_ERROR("SubEntitySet::read_geometry failed or returned a mismatched dimension.\n");
    return;
  }
  
  partition->entitySet = this;
  partition->entitySetNext = lowerOrderHead;
  lowerOrderHead = partition;
  
  partition->entitySetId = id;
}
/*
void SubEntitySet::fix_duplicate_id( int id )
{
  for( int i = 0; i < 2; i++ )
  {
    PartitionEntit* list = i ? subEntityHead : lowerOrderHead;
    for( ; list; list = list->entitySetNext )
    {
      if( list->entitySetId == id )
        list->entitySetId = ++lastId;
    }
  }
}
*/
void SubEntitySet::add_partition( PartitionEntity* partition,
                                  PartitionEntity* after )
{
  assert( partition->entitySet == 0 && partition->entitySetNext == 0 );
  partition->entitySet = this;
  
  if( after )
  {
    assert( after->entitySet == this );
    partition->entitySetNext = after->entitySetNext;
    after->entitySetNext = partition;
  }
  else
  {
    partition->entitySetNext = subEntityHead;
    subEntityHead = partition;
  }
  partition->entitySetId = ++lastId;
}


void SubEntitySet::remove( PartitionEntity* partition )
{
  assert( partition->entitySet == this );
  partition->entitySet = 0;
  
  if( lowerOrderHead == partition )
  {
    lowerOrderHead = partition->entitySetNext;
  }
  else if( subEntityHead == partition )
  {
    subEntityHead = partition->entitySetNext;
  }
  else
  {
    PartitionEntity* curr = subEntityHead;
    while( curr && curr->entitySetNext != partition )
      curr = curr->entitySetNext;
    if( curr == 0 )
    {  
      curr = lowerOrderHead;
      while( curr && curr->entitySetNext != partition )
      	curr = curr->entitySetNext;
    }
    
    assert( curr && curr->entitySetNext == partition );
    curr->entitySetNext = partition->entitySetNext;
  }
  
  partition->entitySetNext = 0;

  if( !subEntityHead && !lowerOrderHead )
  {
    delete this;
  }
  
  partition->entitySetId = 0;
}

CubitStatus SubEntitySet::bridge_destroyed( TopologyBridge* bridge )
{
  if( myEntity && myEntity == bridge )
  {
    if( myEntity->owner() == this )
      myEntity->owner(0);
    myEntity = 0;
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

CubitStatus SubEntitySet::remove_bridge( TopologyBridge* bridge )
{
  if (myEntity && myEntity == bridge)
  {
    strip_attributes(myEntity);
    return bridge_destroyed(myEntity);
  }
  return CUBIT_FAILURE;
}

CubitStatus SubEntitySet::swap_bridge( TopologyBridge* remove,
                                       TopologyBridge* add,
                                       bool reversed )
{
  if( remove_bridge( remove ) )
  {
    if( add->owner() )
    {
      add->owner()->remove_bridge( add );
    }
    add->owner(this);
    myEntity = add;
    if (reversed)
      notify_reversed(add);
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}

void SubEntitySet::get_owners( DLIList<TopologyBridge*>& owner_list ) const
{
  for( PartitionEntity* ent = subEntityHead; ent; ent = ent->entitySetNext )
    owner_list.append( dynamic_cast<TopologyBridge*>(ent) );
    
  CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(get_entity());
  if( coedge && coedge->sense() == CUBIT_REVERSED )
    owner_list.reverse();
}

    
bool SubEntitySet::has_multiple_sub_entities() const
  { return subEntityHead && subEntityHead->entitySetNext; }

void SubEntitySet::print_debug_info( const char* prefix ) const
{
  if( !prefix ) prefix = "";  

  PRINT_INFO("%sSubEntitySet for %s %p\n", prefix, 
    myEntity ? fix_type_name(typeid(*myEntity).name()) : "TopologyBridge", 
    (void*)myEntity );

  PRINT_INFO("%s  SubEntities:\n", prefix );
  PartitionEntity* ent = subEntityHead;
  
  while( ent )
  {
    PartitionSurface *partition_surf = dynamic_cast<PartitionSurface*>(ent);
    if( partition_surf )
      PRINT_INFO("%p is a partition surface\n", (void*)partition_surf);

    PRINT_INFO("%s    %s %d (%p)\n", prefix, 
      fix_type_name(typeid(*ent).name()), ent->entitySetId, (void*)ent);
    ent = ent->entitySetNext;
  }

  PRINT_INFO("%s  Lower-order entities:\n", prefix );
  ent = lowerOrderHead;
  while( ent )
  {
    PRINT_INFO("%s    %s %d (%p)\n", prefix, 
      fix_type_name(typeid(*ent).name()), ent->entitySetId, (void*)ent);
    ent = ent->entitySetNext;
  }
  
  PRINT_INFO("%s  Body: %p\n", prefix, (void*)bodyPtr);
}

bool SubEntitySet::is_attribute( const CubitSimpleAttrib& csa, int id ) const
{
  if( csa.string_data_list().size() < 1 )
    return false;

  if( csa.character_type() != PARTITION_DATA_ATTRIB_NAME )
    return false;
  
  assert( csa.int_data_list().size() );
  if( id && csa.int_data_list()[0] != id )
    return false;
  
  return true;
}    

CubitStatus SubEntitySet::wrap_attribute( CubitSimpleAttrib& csa, int id ) const
{
  if( is_attribute( csa ) )
    return CUBIT_FAILURE;
    
  csa.string_data_list().insert(csa.string_data_list().begin(), CubitString(PARTITION_DATA_ATTRIB_NAME) );
  csa.int_data_list().insert(csa.int_data_list().begin(), id );
  return CUBIT_SUCCESS;
}

int SubEntitySet::unwrap_attribute( CubitSimpleAttrib& csa ) const
{
  if( ! is_attribute( csa ) )
    return 0;

  csa.string_data_list().erase(csa.string_data_list().begin());
  int id = csa.int_data_list()[0];
  csa.int_data_list().erase(csa.int_data_list().begin());
  
  return id;
}

void SubEntitySet::add_attribute( PartitionEntity* entity, const CubitSimpleAttrib& csa )
{
  assert( entity->entitySet == this );
  
  if( !myEntity )
    return;
  
  CubitSimpleAttrib copy = csa;
  if( wrap_attribute( copy, entity->entitySetId ) )
  {
    myEntity->append_simple_attribute_virt( copy );
  }
}

void SubEntitySet::rem_attribute( PartitionEntity* entity, const CubitSimpleAttrib& csa )
{
  assert( entity->entitySet == this );
  
  if( !myEntity )
    return;
  
  CubitSimpleAttrib copy = csa;
  if( wrap_attribute( copy, entity->entitySetId ) )
  {
    myEntity->remove_simple_attribute_virt( copy );
  }
}

void SubEntitySet::get_attributes( PartitionEntity* entity,
                                   DLIList<CubitSimpleAttrib>& list )
{
  assert( entity->entitySet == this );
  assert( !list.size() );
  list.clean_out();
  
  if( !myEntity )
    return;

  DLIList<CubitSimpleAttrib> tmp;
  myEntity->get_simple_attribute( PARTITION_DATA_ATTRIB_NAME, tmp );
  for(int i=0; i<tmp.size(); i++)
  {
    CubitSimpleAttrib& csa = tmp[i];
    if( is_attribute( csa, entity->entitySetId ) )
    {
      unwrap_attribute( csa );
      list.append(csa);
    }
  }
}
void SubEntitySet::get_attributes( PartitionEntity* entity, const char* name,
                                   DLIList<CubitSimpleAttrib>& list )
{
  list.clean_out();

  DLIList<CubitSimpleAttrib> tmp;
  get_attributes(entity,tmp);
  
  for ( int i=0; i<tmp.size(); i++)
  {
    if(tmp[i].character_type() == name)
    {
      list.append(tmp[i]);
    }
  }
}

void SubEntitySet::rem_all_attrib( PartitionEntity* entity )
{
  assert( entity->entitySet == this );
  
  if( !myEntity )
    return;

  DLIList<CubitSimpleAttrib> dead_list;
  myEntity->get_simple_attribute( PARTITION_DATA_ATTRIB_NAME, dead_list );
  for( int i = dead_list.size(); i--; )
  {
    const CubitSimpleAttrib& csa = dead_list.get_and_step();
    if( is_attribute( csa, entity->entitySetId ) )
      myEntity->remove_simple_attribute_virt( csa );
  }
}

int SubEntitySet::get_id( PartitionEntity* entity ) const
{
  if( entity->entitySet != this )
    return 0;
  return entity->entitySetId;
}

PartitionEntity* SubEntitySet::entity_from_id( int id ) const
{
  PartitionEntity* ent;
  for( ent = subEntityHead; ent; ent = ent->entitySetNext )
    if( ent->entitySetId == id )
      return ent;
  for( ent = lowerOrderHead; ent; ent = ent->entitySetNext )
    if( ent->entitySetId == id )
      return ent;
  return 0;
}

void SubEntitySet::get_sub_entities( DLIList<PartitionEntity*>& list ) const
{
  for( PartitionEntity* ent = subEntityHead; ent; ent = ent->entitySetNext )
    list.append(ent);
}


void SubEntitySet::get_lower_order( DLIList<PartitionEntity*>& list ) const
{
  for( PartitionEntity* ent = lowerOrderHead; ent; ent = ent->entitySetNext )
    list.append(ent);
}

void SubEntitySet::set_id( PartitionEntity* entity, int id )
{
  if( entity->entitySet != this )
    assert(0);
  else
  {
    entity->entitySetId = id;
    if( id > lastId )
      lastId = id;
  }
}

void SubEntitySet::renumerate( int lowest_value, bool higher_only )
{
  lastId = lowest_value - 1;
  
  PartitionEntity* ent;
  for( ent = subEntityHead; ent; ent = ent->entitySetNext )
    if( !higher_only || ent->entitySetId >= lowest_value )
      ent->entitySetId = ++lastId;
  for( ent = lowerOrderHead; ent; ent = ent->entitySetNext )
    if( !higher_only || ent->entitySetId >= lowest_value )
      ent->entitySetId = ++lastId;
}

void SubEntitySet::notify_reversed( TopologyBridge* bridge )
{
  assert( !bridge || bridge == myEntity );
  
    // reverse list order (only really necessary for curves)
  PartitionEntity* entity = subEntityHead;
  subEntityHead = 0;
  while( entity )
  {
    PartitionEntity* next = entity->entitySetNext;
    entity->entitySetNext = subEntityHead;
    subEntityHead = entity;
    entity = next;
  }
  
    // reverse the sense of all subentities
  for( entity = subEntityHead; entity; entity = entity->entitySetNext )
    entity->reverse_sense();  
}    


CubitStatus SubEntitySet::save_geometry( int id, int dimension, 
                                         DLIList<CubitVector*>* point_list,
                                         DLIList<int>* point_connectivity,
                                         DLIList<int>* topo_connectivity,
                                         DLIList<int>* point_owners,
                                         CubitSimpleAttrib& attrib )
{
  int i;

  attrib.int_data_list().clear();
  attrib.double_data_list().clear();
  
  attrib.int_data_list().push_back( id);
  attrib.int_data_list().push_back( dimension);
  int point_list_size = point_list ? point_list->size() : 0;
  attrib.int_data_list().push_back( point_list_size);
  int point_conn_size = point_connectivity ? point_connectivity->size() : 0;
  attrib.int_data_list().push_back( point_conn_size);
  int topo_conn_size = topo_connectivity ? topo_connectivity->size() : 0;
  attrib.int_data_list().push_back( topo_conn_size);
  
  if( point_list )
  {
    point_list->reset();
    for( i = point_list->size(); i--; )
    {
      CubitVector* p = point_list->get_and_step();
      attrib.double_data_list().push_back(p->x());
      attrib.double_data_list().push_back(p->y());
      attrib.double_data_list().push_back(p->z());
      delete p;
    }
    point_list->clean_out();
  }
    
  if( point_connectivity )
  {
    point_connectivity->reset();
    for( i = point_connectivity->size(); i--; )
      attrib.int_data_list().push_back(point_connectivity->get_and_step());
  }
    
  if( topo_connectivity )
  {
    topo_connectivity->reset();
    for( i = topo_connectivity->size(); i--; )
      attrib.int_data_list().push_back(topo_connectivity->get_and_step());
  }
  
  if ( point_owners )
  {
    point_owners->reset();
    for ( i = point_owners->size(); i--; )
      attrib.int_data_list().push_back(point_owners->get_and_step());
  }
 
  myEntity->append_simple_attribute_virt( attrib );

  return CUBIT_SUCCESS;
}


CubitStatus SubEntitySet::read_geometry( int& id, int& dimension, 
                                         DLIList<CubitVector*>& point_list,
                                         DLIList<int>& point_connectivity,
                                         DLIList<int>& topo_connectivity,
                                         DLIList<int>& point_owners,
                                         const CubitSimpleAttrib& attrib )
{
  int i;
  
    // read metadata
  int ioffset = 0;
  int doffset = 0;
  id              = attrib.int_data_list()[ioffset++];
  dimension       = attrib.int_data_list()[ioffset++];
  int point_count = attrib.int_data_list()[ioffset++];
  int conn_count  = attrib.int_data_list()[ioffset++];
  int topo_count  = attrib.int_data_list()[ioffset++];
  int owner_count = conn_count ? 2 * point_count : 0;
  
    // consistancy checks
  if( id < 0 || dimension < 0 || dimension > 3 ||
      point_count < 0 || conn_count < 0 || topo_count < 0 )
  {
    PRINT_ERROR("Invalid %s attribute read.  Corrupted file?\n"
                "\tId = %d\n\tDimension = %d\n\tPoint Count = %d\n"
                "\tPoint Connectivity Count = %d\n"
                "\tTopology Connectivity Count = %d\n",
                attrib.string_data_list()[0].c_str(),
                id, dimension, point_count, conn_count, topo_count );
    return CUBIT_FAILURE;
  }
  
  if( (size_t)(conn_count + topo_count + owner_count + 5) != attrib.int_data_list().size() ||
      attrib.int_data_list().size() < 5 )
  {
    PRINT_ERROR("Invalid %s attribute read.  Corrupted file?\n"
                "\tInsufficient/Excess integer data:\n"
                "\tAvailable integer data: %d\n"
                "\tMetadata              : 5\n"
                "\tPoint Connectivity    : %d\n"
                "\tTopology Connectivity : %d\n"
                "\tFacet Point Assoc.    : %d\n"
                "\tTotal Expected Ints   : %d\n",
                attrib.string_data_list()[0].c_str(),
                (int) attrib.int_data_list().size(),
                conn_count, topo_count, owner_count,
                conn_count + topo_count + owner_count + 5 );
    return CUBIT_FAILURE;
  }
                
  if( (size_t)(point_count * 3) != attrib.double_data_list().size() )
  {
    PRINT_ERROR("Invalid %s attribute read.  Corrupted file?\n"
                "\tInsufficient/Excess real data: %d for %d points.\n",
                attrib.string_data_list()[0].c_str(),
                (int) attrib.double_data_list().size(), point_count );
    return CUBIT_FAILURE;
  }
  
    // read point coordinates
  for( i = 0; i < point_count; i++ )
  {
    double x = attrib.double_data_list()[doffset++];
    double y = attrib.double_data_list()[doffset++];
    double z = attrib.double_data_list()[doffset++];
    point_list.append( new CubitVector(x,y,z) );
  }
  
    // read point connectivity (facets)
  for( i = 0; i < conn_count; i++ )
    point_connectivity.append( attrib.int_data_list()[ioffset++] );
  
    // read topology connectivity (children)
  for( i = 0; i < topo_count; i++ )
    topo_connectivity.append( attrib.int_data_list()[ioffset++] );
    
    // read facet point owners
  for( i = 0; i < owner_count; i++ )
    point_owners.append( attrib.int_data_list()[ioffset++] );
  
  return CUBIT_SUCCESS;
}


int SubEntitySet::get_geom_dimension( const CubitSimpleAttrib& attrib )
{
  assert(attrib.int_data_list().size() > 1);
  return attrib.int_data_list()[1];
}

int SubEntitySet::get_geom_id( const CubitSimpleAttrib& attrib )
{
  assert(attrib.int_data_list().size() > 1);
  return attrib.int_data_list()[0];
}

int SubEntitySet::get_segment_count( const CubitSimpleAttrib& attrib )
{
  assert(attrib.int_data_list().size() > 2);

  int point_count = attrib.int_data_list()[2];
  if ( point_count == 0 )
    return 0;

  assert( point_count > 1 );
  return point_count - 1;
}

CubitStatus SubEntitySet::save_geometry()
{
  if( !lowerOrderHead )
    return CUBIT_SUCCESS;
    
  DLIList<CubitSimpleAttrib> old_list;
  myEntity->get_simple_attribute(PARTITION_GEOM_ATTRIB_NAME, old_list);
  while (old_list.size())
  {
    CubitSimpleAttrib old_attrib = old_list.pop();
    myEntity->remove_simple_attribute_virt(old_attrib);
  }  
  
  CubitSimpleAttrib attrib(PARTITION_GEOM_ATTRIB_NAME);
  
  CubitStatus result = CUBIT_SUCCESS;
  for( int i = 0; i < 2; i++ )
  {
    PartitionEntity* list_node = i ? subEntityHead : lowerOrderHead;
    for( ; list_node; list_node = list_node->entitySetNext )
      if(dynamic_cast<GeometryEntity*>(list_node) && !list_node->save(attrib))
        result = CUBIT_FAILURE;
  }
  
  return result;
}
 

void SubEntitySet::strip_attributes() 
{
  strip_attributes( myEntity );
}
      
  
