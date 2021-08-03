//-------------------------------------------------------------------------
// Filename      : OCCShell.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <stddef.h>
// ********** END STANDARD INCLUDES        **********
// ********** BEGIN CUBIT INCLUDES         **********
#include "CastTo.hpp"
#include "CubitUtil.hpp"

#include "OCCQueryEngine.hpp"
#include "OCCShell.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Surface.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCSurface.hpp"
#include "GfxDebug.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "TopoDS.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_ListOfShape.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
#include "BRepBuilderAPI_MakeShape.hxx"
#include "BRepBuilderAPI_ModifyShape.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepBuilderAPI_GTransform.hxx"
#include "LocOpe_SplitShape.hxx"
#include "TopoDS_Compound.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********


//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to an OCC SHELL.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::OCCShell(TopoDS_Shell *theShell, OCCSurface* surface)
{
  myTopoDSShell = theShell;
  mySheetSurface = surface;
  if(myTopoDSShell && !myTopoDSShell->IsNull())
    assert (myTopoDSShell->ShapeType() == TopAbs_SHELL);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::~OCCShell()
{
  if(myTopoDSShell)
  {
    myTopoDSShell->Nullify();
    delete (TopoDS_Shell*)myTopoDSShell;
    myTopoDSShell = NULL;
  }
}

void OCCShell::set_TopoDS_Shell(TopoDS_Shell shell)
{
  if(myTopoDSShell && shell.IsEqual(*myTopoDSShell))
    return;

  if (myTopoDSShell)
    myTopoDSShell->Nullify();

  *myTopoDSShell = shell;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCShell::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

void OCCShell::append_simple_attribute_virt(const CubitSimpleAttrib&)
{
}
void OCCShell::remove_simple_attribute_virt(const CubitSimpleAttrib& )
{
}
void OCCShell::remove_all_simple_attribute_virt()
{
}
CubitStatus OCCShell::get_simple_attribute(DLIList<CubitSimpleAttrib>&)
{
  return CUBIT_FAILURE;
}
CubitStatus OCCShell::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib>&)
  { return CUBIT_FAILURE; }

//-------------------------------------------------------------------------
// Purpose       : Query solid modeler topology
//
// Special Notes : 
//
// Author        : Jane Hu 
//
// Creation Date : 11/28/07
//-------------------------------------------------------------------------
void OCCShell::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
{ 
  if(mySheetSurface)
  {
    parents.append(mySheetSurface->my_lump());
    return;
  }

  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCBody * body = NULL;
  DLIList <OCCBody* > *bodies = oqe->BodyList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  bodies->size(); i++)
  {
     body = bodies->get_and_step();
     TopoDS_Shape* shape ; 
     body->get_TopoDS_Shape(shape);
     assert(shape != NULL); 
     TopExp::MapShapesAndAncestors(*shape,
				   TopAbs_SHELL, TopAbs_SOLID, M);
     TopoDS_Shell *shell = get_TopoDS_Shell(); 
     if (!M.Contains(*shell))
	continue;

     const TopTools_ListOfShape& ListOfShapes = 
				M.FindFromKey(*(get_TopoDS_Shell()));
     if (!ListOfShapes.IsEmpty()) 
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
	   TopoDS_Solid Solid = TopoDS::Solid(it.Value());
           int k = oqe->OCCMap->Find(Solid);
	   parents.append((OCCLump*)(oqe->OccToCGM->find(k))->second);
	 }
     } 
  }
}


void OCCShell::get_children_virt( DLIList<TopologyBridge*>& children )
{
  if(mySheetSurface)
  {
    children.append(mySheetSurface);
    return;
  }
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShell, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *surface = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
          if(surface)
	    children.append_unique(surface);
  }
}

//----------------------------------------------------------------
// Function: to update the core Shell
//           for any movement of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCShell::update_OCC_entity( BRepBuilderAPI_ModifyShape *aBRepTrsf,
                                        BRepAlgoAPI_BooleanOperation *op)
{
  if(mySheetSurface && op == NULL)
    return CUBIT_FAILURE;

  assert (aBRepTrsf != NULL || op != NULL);

  TopoDS_Shape shape;
  if(aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Shell());

  else if(!mySheetSurface)
  {
    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(*get_TopoDS_Shell()));
    if(shapes.Extent() == 0)
      shapes.Assign(op->Generated(*get_TopoDS_Shell()));
    if (shapes.Extent())
      shape = shapes.First();
    
    else if(op->IsDeleted(*get_TopoDS_Shell()))
      ;
    else
      return CUBIT_SUCCESS;
  } 
  TopoDS_Shell shell;
  if (!shape.IsNull())
    shell = TopoDS::Shell(shape);

  //set the surfaces
  DLIList<TopologyBridge *> surfaces;
  this->get_children_virt(surfaces);
  for (int i = 1; i <= surfaces.size(); i++)
  {
     OCCSurface *surface = CAST_TO(surfaces.get_and_step(), OCCSurface);
     surface->update_OCC_entity(aBRepTrsf, op);
  }
  OCCQueryEngine::instance()->update_OCC_map(*myTopoDSShell, shell);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Returns the area of the Shell
//
//-------------------------------------------------------------------------
double OCCShell::measure()
{
  GProp_GProps myProps;
  BRepGProp::SurfaceProperties(*myTopoDSShell, myProps);
  return myProps.Mass();
}

//----------------------------------------------------------------
// Function: TopoDS_Shape level function to update the core Shell
//           for any boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCShell::update_OCC_entity(TopoDS_Shell& old_shell,
                                        TopoDS_Shape& new_shell,
                                        BRepBuilderAPI_MakeShape *op,
                                        LocOpe_SplitShape* sp)
{
  //set the surfaces
  TopTools_IndexedMapOfShape M;
  TopoDS_Shape shape;
  TopExp::MapShapes(old_shell, TopAbs_FACE, M);
  TopTools_ListOfShape shapes;
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Face face = TopoDS::Face(M(ii));
    if (op)
    {
      shapes.Assign(op->Modified(face));
      if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(face));
    }
    else if(sp)
      shapes.Assign(sp->DescendantShapes(face));

    if(shapes.Extent() == 1)
      shape = shapes.First();
    else if(shapes.Extent() > 1)
    {
      shape = shapes.First() ;
    }
    else 
    {
      TopoDS_Face nullFace;
      shape = nullFace;
    }
    if(shapes.Extent() > 0 || (op && op->IsDeleted(face)))
      OCCSurface::update_OCC_entity(face,shape, op, NULL, sp);
  }
  if(!old_shell.IsSame(new_shell))
    OCCQueryEngine::instance()->update_OCC_map(old_shell, new_shell);
  return CUBIT_SUCCESS;
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

