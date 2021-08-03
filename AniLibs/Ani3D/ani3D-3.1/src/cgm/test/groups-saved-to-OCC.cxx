#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "CubitCompat.hpp"
#include "CubitFacetData.hpp"
#include "CubitPointData.hpp"
#include "CubitQuadFacetData.hpp"
#include "DLIList.hpp"
#include "OCCModifyEngine.hpp"
#include "Surface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "BodySM.hpp"
#include "RefEntityFactory.hpp"
#include "RefFace.hpp"
#include "RefGroup.hpp"
#include "RefVertex.hpp"
#include "GroupingEntity.hpp"
#include "SenseEntity.hpp"

#include <iostream>

bool generateThingy(
  RefGroup* container,
  DLIList<RefEntity*>& model)
{
  GeometryModifyTool* gmt =
    GeometryModifyTool::instance();

  Body* body = gmt->pyramid(1., 4, 1., 1., 0.);
  if (!body) return false;

  container->add_ref_entity(body);
  // Add first face of first shell does not help:
  //container->add_ref_entity(body->get_first_sense_entity_ptr()->get_basic_topology_entity_ptr());
  model.append(body);
  // model.append(container);
  return true;
}

int main(int argc, char* argv[])
{
  const char *arg = "groups_saved.occ";
  // Initialize CGM
  CubitStatus s = InitCGMA::initialize_cgma( "OCC" );
  if (CUBIT_SUCCESS != s) return 1;

  CGMApp::instance()->attrib_manager()->auto_flag(CUBIT_TRUE);
  RefGroup* container =
    RefEntityFactory::instance()->construct_RefGroup("Foobly");
  DLIList<RefEntity*> model;

  //DLIList<TopologyBridge*> model;
  if (generateThingy(container, model))
    {
    DLIList<RefEntity*> containedEnts;
    container->get_sub_entities(containedEnts);
    std::cout << "Container has " << containedEnts.size() << " entries\n";
    int numExp;
    CubitString version = "13.1";
    CubitStatus s =
      CubitCompat_export_solid_model(model, arg, "OCC", numExp, version, NULL);
    std::cout << "Wrote \"" << arg << "\" (" << numExp << " entity, version \"" << version.c_str() << "\")\n";

    DLIList<RefGroup*> allGroups;
    GeometryQueryTool::instance()->ref_groups(allGroups);
    std::cout << "There are now " << allGroups.size() << " groups\n";
    int before = static_cast<int>(allGroups.size());

    // Now delete the entities and import the file to see
    // if the group was saved.
    std::cout << "Clearing geometry\n";
    GeometryQueryTool::instance()->delete_geometry();
    allGroups.clean_out();
    GeometryQueryTool::instance()->ref_groups(allGroups);
    std::cout << "There are now " << allGroups.size() << " groups\n";

    s = CubitCompat_import_solid_model(arg, "OCC",
      /* logfile_name = */     NULL,
      /* heal_step = */        CUBIT_TRUE,
      /* import_bodies = */    CUBIT_TRUE,
      /* import_surfaces = */  CUBIT_TRUE,
      /* import_curves = */    CUBIT_TRUE,
      /* import_vertices = */  CUBIT_TRUE,
      /* free_surfaces = */    CUBIT_TRUE,
      &model);
    if (s != CUBIT_SUCCESS)
      {
      std::cerr << "Could not import \"" << argv[1] << "\"\n";
      return 1;
      }
    std::cout << "Imported geometry\n";
    allGroups.clean_out();
    GeometryQueryTool::instance()->ref_groups(allGroups);
    std::cout << "There are now " << allGroups.size() << " groups";
    if (allGroups.size() != before)
      {
      std::cout << ", but there should have been " << before << "\n";
      return 1;
      }
    std::cout << ", as expected.\n";
    }

  return 0;
}
