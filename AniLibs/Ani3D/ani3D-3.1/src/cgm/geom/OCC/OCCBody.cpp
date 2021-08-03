//-------------------------------------------------------------------------
// Filename      : OCCBody.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/00
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CastTo.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "OCCBody.hpp"
#include "CubitSimpleAttrib.hpp"
#include "OCCQueryEngine.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "OCCSurface.hpp"
#include "CubitTransformMatrix.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCAttribSet.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "BRepBuilderAPI_ModifyShape.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepBuilderAPI_GTransform.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "gp_Ax1.hxx"
#include "gp_Ax2.hxx"
#include "gp_GTrsf.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Builder.hxx"
#include "TopoDS.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "Standard_Boolean.hxx"
#include "LocOpe_SplitShape.hxx"
#include "TopoDS_Compound.hxx"
//-------------------------------------------------------------------------
// Purpose       : A constructor with a list of lumps that are attached.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCBody::OCCBody(TopoDS_Compound *theShape, 
                 OCCSurface* surface, OCCShell* shell, Lump* lump)
{
  myTopoDSShape = theShape;
  if (surface != NULL)
    mySheetSurfaces.append(surface);
  if( shell != NULL)
    myShells.append(shell);
  if (lump != NULL)
    myLumps.append(lump);
  update_bounding_box();
  if (myTopoDSShape && !myTopoDSShape->IsNull())
    assert(myTopoDSShape->ShapeType() == TopAbs_COMPOUND);
}

TopoDS_Compound* OCCBody::get_TopoDS_Shape()
{
  if (myTopoDSShape && !myTopoDSShape->IsNull())
  {
    assert(myTopoDSShape->ShapeType() == TopAbs_COMPOUND);
    return myTopoDSShape;
  }
  return (TopoDS_Compound*)NULL;
}

void OCCBody::get_TopoDS_Shape(TopoDS_Shape *& shape)
{
  if (myTopoDSShape && !myTopoDSShape->IsNull())
  {
    assert(myTopoDSShape->ShapeType() == TopAbs_COMPOUND);
    shape =  myTopoDSShape;
  }
  else
  {
    DLIList<Lump*> lumps = this->lumps();
    DLIList<OCCShell*> shells = this->shells();
    DLIList<OCCSurface*> surfaces = this->my_sheet_surfaces();
    if(lumps.size() + shells.size() + surfaces.size() == 1)
    {
      if(lumps.size() == 1)
      {
        OCCLump* lump = CAST_TO(lumps.get(), OCCLump);
        assert(lump != NULL);
        shape = CAST_TO(lump, OCCLump)->get_TopoDS_Solid();
      }
      else if(shells.size() == 1)
        shape = shells.get()->get_TopoDS_Shell(); 
      else
        shape = surfaces.get()->get_TopoDS_Face();
    }
    else
      shape = NULL;
  }
}

void OCCBody::set_sheet_surfaces(DLIList<OCCSurface*> surfaces)
{
  mySheetSurfaces = surfaces;
}

void OCCBody::lumps(DLIList<Lump*>& my_lumps)
{
  myLumps = my_lumps;
}

CubitBoolean OCCBody::is_sheet_body()
{
  if(myLumps.size() == 0 && myShells.size() == 0 && mySheetSurfaces.size() == 1)
    return CUBIT_TRUE;
  return CUBIT_FALSE;
}

void OCCBody::set_TopoDS_Shape( TopoDS_Compound& theshape)
{
  if(!theshape.IsNull())
    assert(theshape.ShapeType() == TopAbs_COMPOUND);

  else
  {
    if(myTopoDSShape && !myTopoDSShape->IsNull())
      myTopoDSShape->Nullify();
    return;
  }
  
  if(myTopoDSShape && !myTopoDSShape->IsNull() && theshape.IsEqual(*myTopoDSShape))
    return;

  if (myTopoDSShape && !myTopoDSShape->IsNull() && theshape.IsPartner(*myTopoDSShape))
  {
    myTopoDSShape->Location(theshape.Location());
    myTopoDSShape->Orientation(theshape.Orientation());
  }

  else
  {
    if(myTopoDSShape)
      myTopoDSShape->Nullify();
    *myTopoDSShape = theshape;
  }
}

OCCBody::OCCBody(DLIList<Lump*>& my_lumps, 
                 DLIList<OCCShell*>& shells,
                 DLIList<OCCSurface*>& surfaces)
{
  myLumps = my_lumps;
  mySheetSurfaces = surfaces;
  myShells = shells;
  TopoDS_Compound* new_top = make_Compound(my_lumps, shells, surfaces);
  myTopoDSShape = new_top;
  assert (myTopoDSShape->ShapeType() == TopAbs_COMPOUND);
  update_bounding_box();
}

TopoDS_Compound* OCCBody::make_Compound(DLIList<Lump*>& my_lumps,
                                        DLIList<OCCShell*>& shells,
                                        DLIList<OCCSurface*>& surfaces)
{
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);
  for(int i = 0; i < my_lumps.size(); i ++)
  {
     OCCLump* lump = CAST_TO(my_lumps.get_and_step(), OCCLump);
     if(!lump)
     {
       PRINT_ERROR("Cannot create an OCC BodySM from the given lumps.\n"
                   "Possible incompatible geometry engines.\n");
       return (TopoDS_Compound *)NULL;
     }
     TopoDS_Solid * solid = CAST_TO(lump, OCCLump)->get_TopoDS_Solid();
     B.Add(Co, *solid);
  }
  for(int i = 0; i < shells.size(); i ++)
  {
     TopoDS_Shell * shell = shells.get_and_step()->get_TopoDS_Shell();
     B.Add(Co, *shell);
  }
  for(int i = 0; i < surfaces.size(); i ++)
  {
     TopoDS_Face * face = surfaces.get_and_step()->get_TopoDS_Face();
     B.Add(Co, *face);
  }

  TopoDS_Compound* new_top = new TopoDS_Compound(Co);
  return new_top; 
}

OCCBody::~OCCBody() 
{
  if (myTopoDSShape)
  {
    myTopoDSShape->Nullify();
    delete (TopoDS_Compound*)myTopoDSShape;
    myTopoDSShape = NULL;
  }
}

GeometryQueryEngine* OCCBody::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}

void OCCBody::append_simple_attribute_virt(const CubitSimpleAttrib &csa)
{ 
  if (myTopoDSShape != NULL)
  {
    OCCAttribSet::append_attribute(csa, *myTopoDSShape);
    return;
  }
  csa_list.append_unique(csa);
}
  
void OCCBody::remove_simple_attribute_virt(const CubitSimpleAttrib &csa)
{ 
  DLIList<Lump*> my_lumps;
  my_lumps = lumps();
  DLIList<OCCShell*> shells = this->shells();
  DLIList<OCCSurface*> surfaces = this->my_sheet_surfaces();

  if (myTopoDSShape != NULL)
  {
    OCCAttribSet::remove_attribute(csa, *myTopoDSShape);
    return;
  }
 
  else if (my_lumps.size() == 1)
  {
    OCCLump* lump = CAST_TO(my_lumps.get(), OCCLump);
    TopoDS_Solid* solid = lump->get_TopoDS_Solid();
    OCCAttribSet::remove_attribute(csa, *solid);
    if(!csa.isEmpty())
      csa_list.remove(csa);
    else
      csa_list.clean_out();
    return;
  }

  else if(shells.size() == 1)
  {
    TopoDS_Shell * shell = shells.get()->get_TopoDS_Shell();
    OCCAttribSet::remove_attribute(csa, *shell);
    if(!csa.isEmpty())
      csa_list.remove(csa);
    else
      csa_list.clean_out();
    return;
  }

  else if(surfaces.size() == 1)
  {
    TopoDS_Face* surf = surfaces.get()->get_TopoDS_Face();
    OCCAttribSet::remove_attribute(csa, *surf);
    if(!csa.isEmpty())
      csa_list.remove(csa);
    else
      csa_list.clean_out();
    return;
  } 
}


void OCCBody::remove_all_simple_attribute_virt()
{
  remove_simple_attribute_virt(CubitSimpleAttrib());
}
  
CubitStatus OCCBody::get_simple_attribute(DLIList<CubitSimpleAttrib>& csas)
{ 
  if (myTopoDSShape != NULL)
    return OCCAttribSet::get_attributes(*myTopoDSShape,csas);

  else
    csas = csa_list;
  return CUBIT_SUCCESS;
}

CubitStatus OCCBody::get_simple_attribute( const CubitString& name,
                                          DLIList<CubitSimpleAttrib>& csas )
{ 
  if (myTopoDSShape != NULL)
    return OCCAttribSet::get_attributes( name, *myTopoDSShape, csa_list );

  for(int i = 0 ; i < csa_list.size(); i ++)
  {
    const CubitSimpleAttrib& csa = csa_list.get_and_step();
    if(csa.string_data_list().size() > 0)
      if (csa.string_data_list()[0] == name)
       csas.append(csa);
  }
  return CUBIT_SUCCESS;
}

CubitStatus OCCBody::get_transforms( CubitTransformMatrix &tfm )
{
  return CUBIT_SUCCESS;
}


//----------------------------------------------------------------
// Function: copy
// Description: create a new copy of the body.
// Author: sjowen
//----------------------------------------------------------------
BodySM* OCCBody::copy()
{
  return (BodySM*)NULL;
}
    
//----------------------------------------------------------------
// Function: move
// Description: translate the body and its child entities
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::move(double dx, double dy, double dz)
{
  double tol = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
  if(fabs(dx) < tol && fabs(dy) < tol && fabs(dz) < tol)
    return CUBIT_SUCCESS;

  gp_Vec aVec(dx, dy, dz);
  gp_Trsf aTrsf;
  aTrsf.SetTranslation(aVec);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);

  return transform(aBRepTrsf);
}


//----------------------------------------------------------------
// Function: rotate
// Description: rotate the body and its child entities
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::rotate( double x, double y, double z, 
                             double angle )//in radians
{
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(x, y, z);
  gp_Ax1 anAxis(aOrigin, aDir);

  //a is angular value of rotation in radians
  gp_Trsf aTrsf;
  aTrsf.SetRotation(anAxis, angle);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);

  return transform(aBRepTrsf);
}

//----------------------------------------------------------------
// Function: transform
// Description: transform the body and its child entities
//              use a transform matrix
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::transform(BRepBuilderAPI_Transform& aBRepTrsf)
{
  TopoDS_Shape * shape;
  get_TopoDS_Shape(shape);
  aBRepTrsf.Perform(*shape);

  update_OCC_entity(&aBRepTrsf);
  // calculate for bounding box
  update_bounding_box();
  
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//              use a constant scale factor
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor )
{
  gp_Trsf aTrsf;
  aTrsf.SetScaleFactor(scale_factor);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  CubitStatus stat = transform(aBRepTrsf);
  return stat;
}

//----------------------------------------------------------------
// Function: scale
// Description: deforming transformation of the body and its child entities
// Author: Jane Hu 
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor_x,
                           double scale_factor_y,
                           double scale_factor_z )
{
  gp_GTrsf gTrsf; 
  gTrsf.SetValue(1,1, scale_factor_x);
  gTrsf.SetValue(2,2, scale_factor_y);
  gTrsf.SetValue(3,3, scale_factor_z);

  BRepBuilderAPI_GTransform gBRepTrsf(gTrsf);

  TopoDS_Shape * shape;
  get_TopoDS_Shape(shape);
  gBRepTrsf.Perform(*shape);

  update_OCC_entity(&gBRepTrsf);
  // calculate for bounding box
  update_bounding_box();

  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: reflect
// Description: reflect the body about an plane given the normal
//              vector.
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::reflect( double reflect_axis_x,
                              double reflect_axis_y,
                              double reflect_axis_z )
{
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(reflect_axis_x, reflect_axis_y,reflect_axis_z); 
  gp_Ax2 anAx2(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAx2);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  return transform(aBRepTrsf);
}

//----------------------------------------------------------------
// Function: private function to update the core compound and      
//           for any movement of the body.
// Note:     input shape must have the same number of Compound
//           as the body's lumps number.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::update_OCC_entity( BRepBuilderAPI_ModifyShape *aBRepTrsf,
                                       BRepAlgoAPI_BooleanOperation *op) 
{
  assert(aBRepTrsf != NULL || op != NULL);

  TopoDS_Compound compsolid;
  TopoDS_Shape shape;
  shape = aBRepTrsf->Shape();
  if(aBRepTrsf && myTopoDSShape)
  {
    compsolid = TopoDS::Compound(shape);
  
    if(OCCQueryEngine::instance()->OCCMap->IsBound(*myTopoDSShape) )
       OCCQueryEngine::instance()->update_OCC_map(*myTopoDSShape, shape);
    else if (!shape.IsEqual(*myTopoDSShape))
       set_TopoDS_Shape(compsolid);
  }

  //Boolean operation works only on one lump body
  //set the lumps
  DLIList<Lump *> lumps;
  lumps = this->lumps();
  for (int i = 1; i <= lumps.size(); i++)
  {
     OCCLump *lump = CAST_TO(lumps.get_and_step(), OCCLump);
     lump->update_OCC_entity(aBRepTrsf, op);
  }

  for(int i = 0; i < mySheetSurfaces.size(); i++)
  {
    OCCSurface* surface = mySheetSurfaces.get_and_step();
    surface->update_OCC_entity(aBRepTrsf, op);
  }
  for(int i = 0; i <myShells.size() ; i++)
  {
    OCCShell* occ_shell = myShells.get_and_step();
    occ_shell->update_OCC_entity(aBRepTrsf,op);
  }

  if (aBRepTrsf && !compsolid.IsNull())
    set_TopoDS_Shape(compsolid);

  update_bounding_box(); 

  //unset marks.
  DLIList<OCCCurve*> curves;
  DLIList<OCCPoint*> points;
  get_all_curves(curves);
  get_all_points(points);

  for(int i = 0; i < curves.size(); i++)
    curves.get_and_step()->set_myMarked(CUBIT_FALSE);

  for(int i = 0; i < points.size(); i++)
    points.get_and_step()->set_myMarked(CUBIT_FALSE);
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: TopoDS_Shape level function to update the core Body
//           for any Boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::update_OCC_entity(TopoDS_Shape& old_shape,
                                       TopoDS_Shape& new_shape,
                                       BRepBuilderAPI_MakeShape *op,
                                       LocOpe_SplitShape* sp)
{
  //set the Shells
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(old_shape, TopAbs_SOLID, M);
  TopTools_IndexedMapOfShape M_new;
  TopExp::MapShapes(new_shape, TopAbs_SOLID, M_new);
  TopTools_ListOfShape shapes;
  TopoDS_Shape shape;

  CubitBoolean updated = CUBIT_FALSE;	
  if(!old_shape.IsNull() && old_shape.ShapeType() == TopAbs_COMPOUND && 
     !new_shape.IsNull() && new_shape.ShapeType() == TopAbs_COMPOUND &&
     !old_shape.IsSame(new_shape))
  {
    //By updating underling solids, shells etc., the old_shape will get changed.
    //trying to make sure the the number of each entity in the old and new 
    //shapes are the same, which means that nothing is delete, that we can 
    //update the map here. Otherwise, when deleting solids, it'll delete the
    //the old body and create new body. This is Ok for general boolean operation    //except imprint when booleans are called, usually the original body are
    // supposed to be kept. 
    updated = CUBIT_TRUE;
    OCCQueryEngine::instance()->update_OCC_map(old_shape, new_shape);
  }
 
  DLIList<int> new_solid_nums;
  DLIList<int> unfound_nums;
  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Solid solid = TopoDS::Solid(M(ii));

    TopTools_ListOfShape shapes;
    if(op)
    {
      shapes.Assign(op->Modified(solid));
      if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(solid));
    }
    else if(sp)
      shapes.Assign(sp->DescendantShapes(solid));

    if (shapes.Extent() == 1)
      shape = shapes.First();

    else if(shapes.Extent() > 1)
    {
      //update all attributes first.
      TopTools_ListIteratorOfListOfShape it;
      it.Initialize(shapes);
      for(; it.More(); it.Next())
      {
        shape = it.Value();
        OCCQueryEngine::instance()->copy_attributes(solid, shape);
      } 
      shape = shapes.First();
    }

    else if(op->IsDeleted(solid))
    {
       if (M_new.Extent()== 1 && ii == 1)
         shape = M_new(1);
       else if(M_new.Extent()== 1 && ii > 1)
         shape.Nullify();
       else if(M_new.Extent() > 1)
       {
         GProp_GProps myProps;
         BRepGProp::VolumeProperties(solid, myProps);
         double bf_mass = myProps.Mass();
         gp_Pnt old_center = myProps.CentreOfMass();
         CubitBoolean found = CUBIT_FALSE;
         for(int l = 1; l <= M_new.Extent(); l++)
         {
           BRepGProp::VolumeProperties(M_new(l), myProps);
           double af_mass = myProps.Mass();
           double dTol = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
           if(fabs(bf_mass-af_mass) < dTol) //unchanged
           {
             gp_Pnt  new_center = myProps.CentreOfMass(); 
             if(new_center.IsEqual(old_center, dTol))
             {
               found = CUBIT_TRUE;
               shape = M_new(l);
               new_solid_nums.append(l);
               break;
             }
           }
         }
         if(!found)
         {
           unfound_nums.append(ii); 
           continue;
         }
       }
       else
         shape.Nullify();
    }
    else
    {
       shape = solid;
       continue;
    }

    if(shapes.Extent() > 0 || (op && op->IsDeleted(solid)))
      OCCLump::update_OCC_entity(solid, shape, op, sp);
  }

  if( unfound_nums.size() == 1 )
  {
    TopoDS_Solid solid = TopoDS::Solid(M(unfound_nums.get()));
    for(int kk = 1; kk <= M_new.Extent(); kk++)
    {
      if(!new_solid_nums.move_to(kk))
      {
        shape = M_new(kk);
        break;
      } 
    }
    OCCLump::update_OCC_entity(solid, shape, op, sp);
  }
  else if(unfound_nums.size() > 1)
  {
    shape.Nullify();
    for(int kk = 1; kk <=unfound_nums.size(); kk++)
    {
       TopoDS_Solid solid = TopoDS::Solid(M(unfound_nums.get_and_step()));
       OCCLump::update_OCC_entity(solid, shape, op, sp);
    }
  } 
  if(!old_shape.IsSame(new_shape) && !updated)
    OCCQueryEngine::instance()->update_OCC_map(old_shape, new_shape);
  return CUBIT_SUCCESS;
}
//----------------------------------------------------------------
// Function: update_bounding_box
// Description: calculate for bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
void OCCBody::update_bounding_box() 
{
  Bnd_Box box;
  TopoDS_Shape shape;
    
  for (int i = 0; i <  myLumps.size(); i++)
  {
     OCCLump *lump = CAST_TO(myLumps.get_and_step(), OCCLump);
     shape = *lump->get_TopoDS_Solid();
     BRepBndLib::Add(shape, box);
  }

  for(int i = 0; i < mySheetSurfaces.size(); i++)
  {
    OCCSurface* surface = mySheetSurfaces.get_and_step();
    shape = *surface->get_TopoDS_Face();
    BRepBndLib::Add(shape, box);
  }

  for(int i = 0; i <myShells.size() ; i++)
  {
    OCCShell* occ_shell = myShells.get_and_step();
    shape = *occ_shell->get_TopoDS_Shell();
    BRepBndLib::Add(shape, box);
  }
 
  //calculate the bounding box
  if(myLumps.size() + mySheetSurfaces.size() + myShells.size() == 0)
  {
    if(!myTopoDSShape)
      return;
    TopoDS_Shape shape = *myTopoDSShape;
    BRepBndLib::Add(shape, box);
  }
  
  double min[3], max[3];

  //get values
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

  //update boundingbox.
  boundingbox.reset(min, max);
}

//----------------------------------------------------------------
// Function: get_bounding_box
// Description: get the  bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
CubitBox OCCBody::get_bounding_box()
{
  return boundingbox ;
}

void OCCBody::get_parents_virt( DLIList<TopologyBridge*>& ) 
  {}
  
void OCCBody::get_children_virt( DLIList<TopologyBridge*>& lumps )
{
  for(int i = 0; i < mySheetSurfaces.size(); i++) 
  { 
    OCCSurface* surface = mySheetSurfaces.get_and_step();
    lumps.append(surface->my_lump());
  }

  for(int i = 0; i <myShells.size() ; i++)
  {
    OCCShell* occ_shell = myShells.get_and_step();
    lumps.append(occ_shell->my_lump());
  }

  for(int i = 0; i <myLumps.size(); i++)
    lumps.append( myLumps.get_and_step()); 
  return;
}

//-------------------------------------------------------------------------
// Purpose       : Find centroid and volume for lumps and shells
//
// Special Notes : 
//
// Author       : Jane Hu  
//
// Creation Date : 11/30/07
//-------------------------------------------------------------------------
CubitStatus OCCBody::mass_properties( CubitVector& centroid, 
                                      double& volume )
{
  if( myShells.size() == 0 && myLumps.size() == 0)
    return CUBIT_FAILURE;
  GProp_GProps myProps;
  TopoDS_Shape* pshape = myTopoDSShape;
  if(!pshape || pshape->IsNull())//single lump or shell or surface
  {
    DLIList<Lump*> lumps = this->lumps();
    if (lumps.size() > 0)
      pshape = CAST_TO(lumps.get(), OCCLump)->get_TopoDS_Solid();
  } 
  if(!pshape || pshape->IsNull())
    return CUBIT_FAILURE;
 
  BRepGProp::VolumeProperties(*pshape, myProps);
  volume = myProps.Mass();
  gp_Pnt pt = myProps.CentreOfMass(); 
  centroid.set(pt.X(), pt.Y(), pt.Z());
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Used to be OCCQueryEngine::is_point_in_body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment OCCBody::point_containment( const CubitVector &point, double /*tolerance*/ )
{
  CubitPointContainment pc_value;
  OCCLump *lump;

  int i;
  for( i=myLumps.size(); i--;)
  {
    lump = dynamic_cast<OCCLump*>(myLumps.get_and_step()); 
    pc_value = lump->point_containment( point );
    if( pc_value == CUBIT_PNT_INSIDE )
      return CUBIT_PNT_INSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  for(int i = 0; i < mySheetSurfaces.size(); i++)
  {
    OCCSurface* surface = mySheetSurfaces.get_and_step();
    pc_value = surface->point_containment( point );
    if( pc_value == CUBIT_PNT_INSIDE )
      return CUBIT_PNT_INSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  for(int i = 0; i <myShells.size() ; i++)
  {
    OCCShell* occ_shell = myShells.get_and_step();
    DLIList<TopologyBridge*> children;
    occ_shell->get_children_virt(children);
    for(int j = 0; j < children.size(); j++)
    {
      OCCSurface* surface = CAST_TO(children.get_and_step(), OCCSurface); 
      pc_value = surface->point_containment( point );
      if( pc_value == CUBIT_PNT_INSIDE )
        return CUBIT_PNT_INSIDE;
      else if( pc_value == CUBIT_PNT_BOUNDARY )
        return CUBIT_PNT_BOUNDARY;
    } 
  }
  return CUBIT_PNT_OUTSIDE;
}

//-------------------------------------------------------------------------
// Purpose       : return all surfaces in this body. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_surfaces(DLIList<OCCSurface*> &surfaces)
{
  TopoDS_Shape *shape;
  get_TopoDS_Shape(shape);

  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*shape, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *face = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCSurface* occ_face = CAST_TO(face, OCCSurface);
       if (occ_face)
         surfaces.append_unique(occ_face);
  }
}

//-------------------------------------------------------------------------
// Purpose       : return all curves in this body.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_curves(DLIList<OCCCurve*> &curves)
{
  TopoDS_Shape *shape;
  get_TopoDS_Shape(shape);

  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*shape, TopAbs_EDGE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *curve = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
       if (occ_curve)
         curves.append_unique(occ_curve);
  }
}

//-------------------------------------------------------------------------
// Purpose       : return all points in this body.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_points(DLIList<OCCPoint*> &points)
{
  TopoDS_Shape *shape;
  get_TopoDS_Shape(shape);

  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*shape, TopAbs_VERTEX, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *vertex = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCPoint* occ_point = CAST_TO(vertex, OCCPoint);
       if (occ_point)
         points.append_unique(occ_point);
  }

  DLIList<OCCSurface*> surfaces;
  this->get_all_surfaces(surfaces);
  for(int i = 0; i < surfaces.size(); i++)
  {
    OCCSurface* occ_surf = surfaces.get_and_step();
    points += occ_surf->get_hardpoints();
  }
}

