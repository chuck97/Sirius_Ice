//-------------------------------------------------------------------------
// Filename      : OCCShell.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 11/16/07
//
// Owner         : 
//-------------------------------------------------------------------------

#ifndef OCC_SHELL_HPP
#define OCC_SHELL_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "ShellSM.hpp"

#include "TopoDS_Shell.hxx"

class OCCBody;
class OCCLump;
class OCCSurface;
class OCCLoop;
class OCCCurve;
class OCCPoint;
class BRepAlgoAPI_BooleanOperation;
class BRepBuilderAPI_ModifyShape;
class BRepBuilderAPI_MakeShape;
class LocOpe_SplitShape;
// ********** END CUBIT INCLUDES           **********

class OCCShell : public ShellSM
{
public:
  
  OCCShell(TopoDS_Shell *theShell, OCCSurface* surface = NULL);

  virtual ~OCCShell() ;
    //- Destructor.

  TopoDS_Shell* get_TopoDS_Shell() 
  {if(myTopoDSShell && !myTopoDSShell->IsNull()) assert (myTopoDSShell->ShapeType() == TopAbs_SHELL); return myTopoDSShell;}
  void set_TopoDS_Shell(TopoDS_Shell shell);

  OCCLump* my_lump() {return myLump;}
  OCCBody* my_body() {return myBody;}
  OCCSurface* my_surface() {return mySheetSurface;}
  void set_sheet_surface (OCCSurface* surface) {mySheetSurface = surface;}

  void set_lump(OCCLump* lump) {myLump = lump;}
  void set_body(OCCBody* body) {myBody = body;}

  bool is_sheet() const { return mySheetSurface != 0; }

  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual void append_simple_attribute_virt(const CubitSimpleAttrib&);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(const CubitSimpleAttrib&);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.

  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );

  //ModifyShape refers only to Transform and GTransform for now 1/10/11)
  CubitStatus update_OCC_entity( BRepBuilderAPI_ModifyShape *aBRepTrsf,
                                 BRepAlgoAPI_BooleanOperation *op = NULL);
  double measure(); //area of the shell
  static CubitStatus update_OCC_entity(TopoDS_Shell& old_shell,
                                       TopoDS_Shape& new_shell,
                                       BRepBuilderAPI_MakeShape *op,
                                       LocOpe_SplitShape* sp = NULL);

  DLIList<OCCSurface*> getMemberSurfaces() { return memberSurfaces; }
  void setMemberSurfaces(DLIList<OCCSurface*> surfaceList)
  {
    memberSurfaces = surfaceList;
  }
protected: 
  
private:
  TopoDS_Shell *myTopoDSShell;
  OCCSurface* mySheetSurface;
  OCCLump*    myLump;
  OCCBody*    myBody;
  DLIList<OCCSurface*> memberSurfaces;
};

#endif

