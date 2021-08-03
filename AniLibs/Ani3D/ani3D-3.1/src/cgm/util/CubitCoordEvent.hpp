//- Class: CubitCoorEvent
//- Description: Class describing Coordinate System events 
//- Owner:

#include "CubitEvent.hpp"
#include "CubitCoordinateSystem.hpp"

class CUBIT_UTIL_EXPORT CubitCoordEvent : public CubitEvent
{
  public:
    enum Type
    {
      COORDINATE_SYSTEM_CREATED,
      COORDINATE_SYSTEM_MODIFIED,
      COORDINATE_SYSTEM_DELETED
    };
    
    CubitCoordEvent(Type type, CubitCoordinateSystem* sys)
      : mType(type), mCubitCoordinateSystem(sys) {}
    ~CubitCoordEvent() {}

    CubitCoordinateSystem* get_entity() const {return mCubitCoordinateSystem; }
    Type get_type() const { return mType; }

  protected:

    Type mType;
    CubitCoordinateSystem* mCubitCoordinateSystem;
};
