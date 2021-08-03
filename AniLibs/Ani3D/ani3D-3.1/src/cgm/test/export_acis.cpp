#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"

int main()
{
  CubitStatus result = InitCGMA::initialize_cgma("ACIS");
  if (CUBIT_SUCCESS != result) return 1;

  // Create a cylinder
  double height = 10.0, major_rad = 3.0, tmp_minor = 3.0;
  Body* cylinder = GeometryModifyTool::instance()->cylinder(height, major_rad, tmp_minor, major_rad);
  if (NULL == cylinder) return 1;

  // Export it to an ACIS SAT file
  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported = 0;
  const char* filename = "export_acis.sat";
  const char* filetype = "ACIS_SAT";
  CubitString cubit_version("14.9");
  CubitStatus stat = CubitCompat_export_solid_model(ref_entity_list, filename, filetype, num_ents_exported, cubit_version);
  if (CUBIT_SUCCESS != stat) return 1;

  GeometryQueryTool::instance()->delete_geometry();

  return 0;
}
