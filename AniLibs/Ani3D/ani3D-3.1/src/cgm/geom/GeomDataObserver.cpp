//------------------------------------------------------------------------
// Class GeomDataObserver
// Description:  Observer class that stores/caches specific geometric
//               information (for example, the area for a surface.
//
// Author: David White
// Creation Date: 9/7/2003
//------------------------------------------------------------------------
#include "GeomDataObserver.hpp"
#include "RefEntity.hpp"
#include "GeometryEvent.hpp"

GeomDataObserver::GeomDataObserver(RefEntity* watched)
    : myRefEntity(watched)
{
  measureSet = CUBIT_FALSE;
    //initialize to something weird.
  myMeasure = -CUBIT_DBL_MAX;
}

GeomDataObserver::~GeomDataObserver()
{ 
  unregister_observable( myRefEntity ); 
}

GeomDataObserver* GeomDataObserver::get( RefEntity* on_this )
{
   DLIList<CubitObserver*> list;
   GeomDataObserver* eo = NULL;

   on_this->get_observer_list(list);
   for (int i = list.size(); i--; )
   {
     if ( (eo = dynamic_cast<GeomDataObserver*>(list.step_and_get()) ))
        break;
   }
   return eo;
}

GeomDataObserver* GeomDataObserver::create( RefEntity* on_this ) {
   GeomDataObserver* eo = get(on_this);
   if (eo)
     return eo;

   eo = new GeomDataObserver(on_this);
   eo->register_observable(on_this);
   return eo;
}

void GeomDataObserver::notify_observer(const CubitEvent* event)
{
   const GeometryEvent* geom_event = dynamic_cast<const GeometryEvent*>(event);

   if(geom_event)
   {
     assert(geom_event->get_entity() == myRefEntity);

     switch (geom_event->get_type())
     {
       case GeometryEvent::GEOMETRY_TOPOLOGY_MODIFIED:
       case GeometryEvent::TOPOLOGY_MODIFIED:
       case GeometryEvent::GEOMETRY_MODIFIED:
       case GeometryEvent::TOPOLOGY_ENTITY_DESTRUCTED:
         break;
       default:
         return;
     }

     /* don't call virtual functions or access
        class data after this! */
     delete this;
   }
}

