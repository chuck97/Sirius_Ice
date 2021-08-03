#include "GeometryQueryEngine.hpp"
#include "CubitString.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "BodySM.hpp"



CubitStatus GeometryQueryEngine::set_export_version(const int,
                                                    const int) 
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::set_export_allint_version(const int) 
{
  return CUBIT_FAILURE;
}

int GeometryQueryEngine::get_allint_version() 
{
  return 100*get_major_version() + get_minor_version();
}

CubitStatus GeometryQueryEngine::list_engine_versions(CubitString &versions)
{
  versions = "(none specified)";
  return CUBIT_SUCCESS;
}

//Added functions for removing the geom classes dependecy on the virtual classes

 
CubitStatus GeometryQueryEngine::get_underlying_curves(Curve * curve_ptr,
                                 DLIList<TopologyBridge*>& curve_list)
{ return CUBIT_SUCCESS; }

CubitStatus GeometryQueryEngine::get_underlying_surfaces(Surface * surf_ptr,
                                 DLIList<TopologyBridge*>& surf_list)
{ return CUBIT_SUCCESS; }

CubitStatus GeometryQueryEngine::get_underlying_bridges(TopologyBridge* bridge_ptr,
                                 DLIList<TopologyBridge*>& bridge_list)
{ return CUBIT_SUCCESS; }


CubitStatus GeometryQueryEngine::get_graphics( BodySM *bodysm, GMem* g_mem,
                         std::vector<Surface*> &surface_to_facets_vector,
                         std::vector<TopologyBridge*> &vertex_edge_to_point_vector,
                         std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > &facet_edges_on_curves,
                         unsigned short normal_tolerance, 
                         double distance_tolerance, double max_edge_length ) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::get_graphics( Surface* surface_ptr,
                                               GMem *gmem,
                                               std::vector<TopologyBridge*> &vertex_edge_to_point_vector,
                                               std::vector<std::pair<TopologyBridge*, std::pair<int,int> > > &facet_edges_on_curves,
                                               unsigned short normal_tolerance, 
                                               double distance_tolerance, 
                                               double max_edge_length ) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::get_visible_entities( TopologyBridge *hidden_tb, 
                                                      DLIList<TopologyBridge*> &real_tbs )
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::delete_solid_model_entities(GeometryEntity* ref_entity_ptr) const
{
  return CUBIT_FAILURE;
}

CubitStatus GeometryQueryEngine::delete_topology_bridge(TopologyBridge* bridge) const
{
 
  // first, try the geometry entity

  GeometryEntity* ent = dynamic_cast<GeometryEntity*> (bridge);
  if (ent)
  {
    return delete_solid_model_entities(ent);
  }

  // now, try to cast to a BodySM and call that routine
  BodySM* body = dynamic_cast<BodySM*> (bridge);
  if (body)
  {
    delete_solid_model_entities(body);
    return CUBIT_SUCCESS;
  }

  PRINT_ERROR("GeometryQueryEngine::delete_solid_model_entities(TopologyBridge* bridge): no deletion...\n");
  return CUBIT_FAILURE;
}


CubitStatus GeometryQueryEngine::get_point_containment( DLIList<BodySM*> &body_list,
                                     DLIList<CubitVector> &point_list,
                                     double tolerance,
                                     bool allow_pts_in_multiple_bodies,
                                     std::vector< std::pair<BodySM*, std::vector<int> > > &bodysm_to_pt_indices )
{
  return CUBIT_FAILURE;
}
