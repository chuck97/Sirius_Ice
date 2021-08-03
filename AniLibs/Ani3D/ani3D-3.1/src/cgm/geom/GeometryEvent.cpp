#include "GeometryEvent.hpp"
#include <stddef.h>

GeometryEvent::GeometryEvent(Type type, RefEntity* entity)
  : mType(type), mRefEntity(entity)
{
}

GeometryEvent::~GeometryEvent()
{
}
    
RefEntity* GeometryEvent::get_entity() const
{
  return mRefEntity;
}

GeometryEvent::Type GeometryEvent::get_type() const
{
  return mType;
}

TopologyEvent::TopologyEvent(Type type, TopologyEntity* entity)
  : mType(type), mRefEntity(entity)
{
}

TopologyEvent::~TopologyEvent()
{
}

TopologyEntity* TopologyEvent::get_entity() const
{
  return mRefEntity;
}

TopologyEvent::Type TopologyEvent::get_type() const
{
  return mType;
}
  
GeometryIdSetEvent::GeometryIdSetEvent(RefEntity* ent, int old_id, int new_id)
  : GeometryEvent(GeometryEvent::ID_SET, ent), oldId(old_id), newId(new_id)
{
}
  
GeometryIdSetEvent::~GeometryIdSetEvent()
{
}

UnMergeEvent::UnMergeEvent( RefEntity *old_ptr, RefEntity *new_ptr )
  : GeometryEvent(GeometryEvent::NEW_ENTITY_UNMERGED, old_ptr),
    mOldEntity(old_ptr),
    mNewEntity(new_ptr)
{}

UnMergeEvent::~UnMergeEvent()
{
}

MergeEvent::MergeEvent( RefEntity *lost_entity,
                        RefEntity *kept_entity )
    : GeometryEvent(GeometryEvent::ENTITIES_MERGED, lost_entity)
{
  keptEntity = kept_entity;
}

MergeEvent::~MergeEvent()
{
}

RefEntity *MergeEvent::get_lost_entity() const
{
  return this->mRefEntity;
}

RefEntity *MergeEvent::get_kept_entity() const
{
  return keptEntity;
}

CompositeCombineEvent::CompositeCombineEvent(GeometryEvent::Type event_type)
  : GeometryEvent(event_type, NULL),
    mKeptEntity(NULL),
    mDeleteEntity(NULL)
{}

CompositeCombineEvent::CompositeCombineEvent(GeometryEvent::Type event_type,
                      RefEntity *keep_ptr,
                      RefEntity *delete_ptr )
  : GeometryEvent(event_type, NULL),
    mKeptEntity(keep_ptr),
    mDeleteEntity(delete_ptr)
{}

CompositeCombineEvent::~CompositeCombineEvent()
{
}

TransformEvent::TransformEvent( const CubitTransformMatrix &transform, std::vector<RefEntity*> ents )
  : GeometryEvent(GeometryEvent::GEOMETRY_TRANSFORMED, NULL ),
  transformation( transform )
{
  mRefEntity = NULL;
  refEnts = ents;  
}

const CubitTransformMatrix TransformEvent::get_transformation() const
{
  return transformation;
}

std::vector<RefEntity*> TransformEvent::get_ref_ents() const
{
  return refEnts;
}

TransformEvent::~TransformEvent()
{
}


