//- Class: CubitMeshEvent
//- Description: Class describing an mesh event such as mesh creation, deletion and modification.
//-
//- Owner:  

#ifndef GeometryEvent_hpp
#define GeometryEvent_hpp

#include "CubitEvent.hpp"
#include "CubitGeomConfigure.h"
#include "CubitTransformMatrix.hpp"
class RefEntity;
class TopologyEntity;

class CUBIT_GEOM_EXPORT GeometryEvent : public CubitEvent
{
  public:
    enum Type
    {
       TOPOLOGY_ENTITY_CONSTRUCTED,
         /* A TopologyEntity was created */
       TOPOLOGY_ENTITY_MODIFIED,
         /* A TopologyEntity was changed */
       TOPOLOGY_ENTITY_DESTRUCTED,
         /* A TopologyEntity was deleted */
       DEVELOPER_COMMAND_FLAG_MODIFIED,
        /* The 'set developer command on/off' was issued */
       MESH_SETTING_MODIFIED,
         /* a mesh setting on a geometry entity was modified */
       GEOMETRY_TOPOLOGY_MODIFIED,
         /* Both geometry and topology was
          modified. e.g. partitioned.
         */
       TOPOLOGY_MODIFIED,
         /* The topology of a RefVolume or one of its sub-entities
          * was changed.  Used for virtual topology. */
       GEOMETRY_MODIFIED,
         /* The geometry of a RefEntity was altered. */
       NEW_ENTITY_UNMERGED,
         /*A surface, curve, or vertex was unmerged. */
       FREE_REF_ENTITY_GENERATED,
         /* A vertex that is not a part of a curve,
          * or a curve that is not a part of a surface,
          * was just generated */
       TOP_LEVEL_ENTITY_DESTRUCTED,
         /* A body or free entity and all it's children
          * are about to be destructed.
          * The pointer is still valid at this point for all
          * calls except geometry queries. 
          * The pointer can be cast to a higher level entity
          * if desired. */

       ENTITY_NAME_CHANGED,
         /* the name for an entity changed */
       
       ENTITY_GEOMETRY_COLOR_CHANGED,
       ENTITY_VISIBILITY_CHANGED,
        /* the visibility of an entity changed */
       
       ENTITIES_MERGED,
         /* two entities are merged together. See MergeEvent. */
       ID_SET,
         /* an entity had its id set (changed). See IdSetEvent. */
       GROUP_MODIFIED,
         // Group was modified.  Sent only once after all
         // modifications are complete.
       GEOMETRY_TRANSFORMED,

       VIRTUAL_STATUS_CHANGED

    };

    GeometryEvent(Type type, RefEntity* entity);
    ~GeometryEvent();

    RefEntity* get_entity() const;
    Type get_type() const;

  protected:

    Type mType;
    RefEntity* mRefEntity;
};

class CUBIT_GEOM_EXPORT TopologyEvent : public CubitEvent
{
  public:
    enum Type
    {
       TOPOLOGY_ENTITY_CONSTRUCTED,
         /* A TopologyEntity was created */
       TOPOLOGY_ENTITY_MODIFIED,
         /* A TopologyEntity was changed */
       TOPOLOGY_ENTITY_DESTRUCTED
      /* A TopologyEntity was deleted */
    };

    TopologyEvent(Type type, TopologyEntity* entity);
    ~TopologyEvent();

    TopologyEntity* get_entity() const;
    Type get_type() const;

  private:

    Type mType;
    TopologyEntity* mRefEntity;
};

class CUBIT_GEOM_EXPORT GeometryIdSetEvent: public GeometryEvent
{
public:
  GeometryIdSetEvent(RefEntity* ent, int old_id, int new_id);
  ~GeometryIdSetEvent();

 int get_old_id() const { return oldId; }
 int get_new_id() const { return newId; }

private:
  int oldId;
  int newId;
};

class CUBIT_GEOM_EXPORT UnMergeEvent : public GeometryEvent
{
  public:

    UnMergeEvent( RefEntity *old_ptr, RefEntity *new_ptr );
    ~UnMergeEvent();

    RefEntity* const mOldEntity;

    RefEntity* const mNewEntity;
};

class CUBIT_GEOM_EXPORT TransformEvent: public GeometryEvent
{
public:
  TransformEvent( const CubitTransformMatrix &transformation, std::vector<RefEntity*> ents );
  ~TransformEvent();

 const CubitTransformMatrix get_transformation() const; 
 std::vector<RefEntity*> get_ref_ents() const;

private:
  std::vector<RefEntity*> refEnts;
  const CubitTransformMatrix transformation;
};

class CUBIT_GEOM_EXPORT MergeEvent: public GeometryEvent
{
public:
  MergeEvent( RefEntity *lost_entity, RefEntity *kept_entity );
  ~MergeEvent();

 RefEntity *get_lost_entity() const;
 RefEntity *get_kept_entity() const;

private:
  RefEntity *keptEntity;
};


class CUBIT_GEOM_EXPORT CompositeCombineEvent : public GeometryEvent
{
  public:
    CompositeCombineEvent(GeometryEvent::Type event_type);
    CompositeCombineEvent(GeometryEvent::Type event_type,
                          RefEntity *keep_ptr,
                          RefEntity *delete_ptr );
    ~CompositeCombineEvent();

    RefEntity* const mKeptEntity;
    RefEntity* const mDeleteEntity;
};



#endif
