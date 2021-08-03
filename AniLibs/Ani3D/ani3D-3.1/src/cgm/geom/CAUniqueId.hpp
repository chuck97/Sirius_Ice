//- Class:          CAUniqueId
//- Owner:          Greg Nielson
//- Description:    Cubit attribute for mesh interval.
//- Checked by:
//- Version:

#ifndef CA_UNIQUE_ID_HPP
#define CA_UNIQUE_ID_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

typedef std::map<long, long> UIDMap;

class CUBIT_GEOM_EXPORT CAUniqueId: public CubitAttrib
{
private:

  int uniqueId;
    //- the unique id tag for an entity

  static DLIList<CAUniqueId *> allCAUniqueIds;
    //- list of all CAUI's; used in actuate_all function

  static bool autoUniqueId;
    //- flag controlling whether uids are automatically created (even when no other 
    //- CA's request them)
    
  static UIDMap oldUIDToNewUID;

public:

  virtual ~CAUniqueId();

  CAUniqueId(RefEntity*, const CubitSimpleAttrib&);
    //- make a CAMI from a simple attribute

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset() {return CUBIT_SUCCESS;};
    //- don't need to do anything, as all the data gets assigned
    //- and not appended

  CubitSimpleAttrib cubit_simple_attrib();

  int unique_id() { return uniqueId;}

  void unique_id (int id) {uniqueId = id;};

  int int_attrib_type() {return CA_UNIQUE_ID;};

  static CubitStatus actuate_all();
    //- actuate all the CAUI's on the list, then empty the list

  static bool auto_unique_id();
  static void auto_unique_id(const bool flag);
    //- get/set autoUniqueId

  virtual void print();
    //- print the value of this attribute
 
  static UIDMap get_old_to_new_uid_map() { return oldUIDToNewUID; }  
  static void clear_out_old_to_new_map();

};

inline bool CAUniqueId::auto_unique_id()
{
  return autoUniqueId;
}

inline void CAUniqueId::auto_unique_id(const bool flag) 
{
  autoUniqueId = flag;
}

CubitAttrib* CAUniqueId_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

