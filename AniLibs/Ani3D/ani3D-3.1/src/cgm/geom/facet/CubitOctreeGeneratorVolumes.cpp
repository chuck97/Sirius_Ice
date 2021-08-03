#include "CubitOctreeGeneratorVolumes.hpp"
#include "CubitOctree.hpp"
#include "RefFace.hpp"
#include "TDOctreeRefFace.hpp"
#include "CubitFacet.hpp" 
#include "CubitOctreeNode.hpp" 
#include "CubitOctreeCell.hpp"
#include "CubitPoint.hpp"
#include "CubitFacetEdge.hpp"
#include "RefEdge.hpp"
#include "CpuTimer.hpp"
#include "FacetDataUtil.hpp"
#include "RefEdge.hpp"
#include "TDOctreeRefEdge.hpp"
#include "GeometryQueryEngine.hpp"
#include "FacetSurface.hpp"
#include "CubitOctreeConstants.hpp"
#include "ChollaEngine.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"
#include "OctreeFacetPointData.hpp"
#include "Body.hpp"

// Constructor for CubitOctree
CubitOctreeGeneratorVolumes::CubitOctreeGeneratorVolumes( DLIList<RefEntity*> &entity_list)
        :CubitOctreeGenerator( ), entityList(entity_list){
  
  entityList = entity_list;
  cubitOctree = new CubitOctree( this);
}

void CubitOctreeGeneratorVolumes::get_bounding_box( CubitVector &min, CubitVector &max ){
  RefEntity *volume;
  int index;
  CubitBox bound_box;

  min.x( CUBIT_DBL_MAX );
  min.y( CUBIT_DBL_MAX );
  min.z( CUBIT_DBL_MAX );
 
  max.x( -CUBIT_DBL_MAX );
  max.y( -CUBIT_DBL_MAX );
  max.z( -CUBIT_DBL_MAX );
  
  for( index = 0; index < entityList.size(); index++ ){
    volume =   entityList.get_and_step();

    bound_box =  volume->bounding_box();

    if( bound_box.minimum().x() < min.x() )
      min.x( bound_box.minimum().x());
    if( bound_box.minimum().y() < min.y() )
      min.y( bound_box.minimum().y());
    if( bound_box.minimum().z() < min.z() )
      min.z( bound_box.minimum().z());

    if( bound_box.maximum().x() > max.x() )
      max.x( bound_box.maximum().x());
    if( bound_box.maximum().y() > max.y() )
      max.y( bound_box.maximum().y());
    if( bound_box.maximum().z() > max.z() )
      max.z( bound_box.maximum().z());    
  }

}

// CubitOctree is generated based on facets of the FACE.
// Therefore small features are captured.
// First the CubitOctree grid is generated without any intersection
// calculation. The nodes are linked using 6 pointers to the adjacent nodes
// The cells are stored using CubitOctree data structure.  As the boundary points are need
// for MAT generation, intersection between the CubitOctree edges and facets are computed and 
// normals are also stored.
CubitBoolean CubitOctreeGeneratorVolumes::generate_lattice( void ){

    // initialize CubitOctree root
  cubitOctree->initialize_octree_generation();
  
    // Builds the tree fully till the level min_depthibcbtghs.lib
  cubitOctree->build_octree_till_min_depth( cubitOctree->get_root() );
  
    // Builds the tree based on facets' points, centroid, dimension ...
    // Note that the points on the curves belong to more than one face
    // Therefore the CubitOctree is subdivided till the maxDepth
  build_octree_till_skl_max_depth_based_on_facets();
 
    // Adjacents cells can be found by visiting the adjacent cells of every node
    // Cells are split till the difference in the levels of adjacent cells is 1.
    // At every node find the maximum and minimum cell.  Split the maximum cell if
    // it deffers by 2 or more in depth.

  
    //PRINT_DEBUG_157(" Number of nodes = %d", gridNodeVector.size() );
  return CUBIT_TRUE;
}
 

CubitBoolean CubitOctreeGeneratorVolumes::build_octree_till_skl_max_depth_based_on_facets( void ){
  
  DLIList< RefFace *> vol_ref_faces;
  DLIList<CubitFacet *> *facet_list;
  DLIList<CubitPoint*> *point_list;
  DLIList<CubitFacetEdge *> *ptr_facet_edge_list;
  DLIList<OctreeFacetPointData *> facet_point_data_list; //, vfacet_point_data_list;
  DLIList<CubitOctreeCell*> refine_cell_list;
  
    // DLIList<CubitFacet*> global_facet_list;
  

  int index;
  int i, j, k, l;
  CubitPoint *ptr_point;
  CubitFacet *ptr_facet, *ptr_adj_facet0, *ptr_adj_facet1;
  CubitFacetEdge *ptr_edge;
  RefFace *ptr_ref_face;
  OctreeFacetPointData *facet_point_data;
  CubitVector normal0, normal1;
  double internal_angle;
  Body *volume;
  TDOctreeRefFace *td_ref_face;

    // PR-CubitOctree has the property that the cell size depends on the density of the facet points
    // Higher the density the finer will be the facet cell size.
    // Step 1: The CubitOctree is built only till MaxDepth based on facet information
    //  1. Small features that have smaller curve length are detected using facet density
    //  2. Surface/Curve curvature are detected based on di-hedral angle between the adj facets.
    //     Also the aspect ratio of facets are used because it is found that when the surface has curvature on only one directio
    //     the facet aspect ratio will be very high.
  for( index = 0; index < entityList.size(); index++ ){
    volume =  CAST_TO( entityList.get_and_step(), Body );
    vol_ref_faces.clean_out();
    volume->ref_faces( vol_ref_faces );
    vol_ref_faces.reset();

      // get volume facets and insert points from those facets
  
    for( i = 0; i < vol_ref_faces.size(); i++ ){
      ptr_ref_face = vol_ref_faces.get_and_step();       
      td_ref_face = TDOctreeRefFace::get_td( ptr_ref_face );

      if( td_ref_face->get_visit() == CUBIT_FALSE ){
        td_ref_face->set_visit( CUBIT_TRUE );

          // grab facets here if stitching

          // add to global list for facet point insertion

          /*int num_parents = CAST_TO(ptr_ref_face,RefFace)->num_parent_ref_entities();
            if (num_parents == 1)
            {
            global_facet_list += *facet_list;
            }
        
            else
            {*/
          
          // Generated OctreeFacetPointData based on facet vertices    
        point_list = td_ref_face->get_ptr_cubit_point_list();    
        for( j = 0; j < point_list->size(); j++ ){
          ptr_point = point_list->get_and_step();
          facet_point_data = new OctreeFacetPointData( ptr_point->coordinates(), ptr_point );
          facet_point_data_list.push( facet_point_data );
        }
//        point_list.clean_out();
        
          //}
        

          // Generate OctreeFacetPointData based on facet centroid
        
        facet_list = td_ref_face->get_ptr_cubit_facet_list();
        for( j = 0; j < facet_list->size(); j++ ){
          ptr_facet = facet_list->get_and_step();
          facet_point_data = new OctreeFacetPointData( ptr_facet->center(), ptr_facet );
          facet_point_data_list.push( facet_point_data );
        }      
       
          /* Generate OctreeFacetPointData based on aspect ratio of facet 
             for( j = 0; j < facet_list->size(); j++ ){
             ptr_facet = facet_list->get_and_step();
             OctreeFacetPointData::generate_facet_point_data_at_slender_facet( ptr_facet, facet_point_data_list );
             }
          */

          // Generate OctreeFacetPointData based on surface curvature
        ptr_facet_edge_list = td_ref_face->get_ptr_cubit_facet_edge_list();
        for( j = 0; j < ptr_facet_edge_list->size(); j++ ){
          ptr_edge = ptr_facet_edge_list->get_and_step();
          if( ptr_edge->num_adj_facets() == 2 ){
              // mark adjacent facets          
            ptr_adj_facet0 = ptr_edge->adj_facet(0);
            ptr_adj_facet1 = ptr_edge->adj_facet(1);
            normal0 = ptr_adj_facet0->normal();
            normal1 = ptr_adj_facet1->normal();
                                            
            if (normal0.length_squared() > 0.0 && normal1.length_squared() > 0.0)
            {
              internal_angle = normal0.interior_angle(normal1);
            }
            else
            {
              internal_angle = 0.0;
            }
              //PRINT_INFO("%f ", internal_angle );
            if( internal_angle > ANGLE_BETWEEN_FACETS_NORMAL_FOR_OCTREE_BASED_ON_SURF_CURVATURE_FACTOR * ANG_FACET_EXTRACT[OCTREE_TIME_ACCURACY_LEVEL] ){
              OctreeFacetPointData::generate_facet_point_data_based_on_curvature( ptr_edge, /*internal_angle,*/ facet_point_data_list );
            }					
          }	
        }
          
          // now check for pairs of adjacents with dihedral angles along their shared edge
          // the CubitOctree should be refined in these regions so we have enough grid nodes to generate 3dmat points
          // disabling for now :( because there is a problem in FacetDataUtil's facet stitching code
        if (0)
        {
          DLIList<DLIList<RefEdge*> > edge_loops;
          ptr_ref_face->ref_edge_loops(edge_loops);
          for (i=0; i < edge_loops.size(); ++i)
          {
            DLIList<RefEdge*>& edge_loop = edge_loops.get_and_step();
            for (j=0; j < edge_loop.size(); ++j)
            {
              DLIList<RefEntity*> parents;
              RefEdge *edge = edge_loop.get_and_step();
              TDOctreeRefEdge *td_ref_edge = TDOctreeRefEdge::get_td(edge);
            
              if( td_ref_edge->get_visit() == CUBIT_FALSE )
              {
                td_ref_edge->set_visit( CUBIT_TRUE );
                edge->get_parent_ref_entities(parents);

                  // check Refface's for bad facets in case it hasn't been done
                  // already
                for (k=0; k < parents.size(); ++k)
                {
                  RefFace *parent = CAST_TO(parents.get_and_step(), RefFace);
                  if (parent == NULL) {continue;}
                  TDOctreeRefFace *td_parent = TDOctreeRefFace::get_td(parent);
                  if (td_parent == NULL) {continue;}
                  if (td_parent->get_create_2dmat())
                  {
                    td_parent->check_valid_facets(CUBIT_TRUE);
                  }
                }
              
                for (k=0; k < parents.size(); ++k)
                {
                  for (l=k+1; l < parents.size(); ++l)
                  {
                      // SVDrawTool::clear_non_retained();
                    DLIList<CubitFacet*> crease_facets;
                    DLIList<CubitFacetEdge*> crease_edges;
                    cubitOctree->find_cells_based_on_surface_angle(refine_cell_list, crease_edges, crease_facets,
                                                              CAST_TO(parents[k],RefFace), CAST_TO(parents[l],RefFace),
                                                              CUBIT_FALSE, 2.0*CUBIT_PI/3.0);
                      // will need to store crease_edges and facets

                      // SVDrawTool::mouse_xforms();
                  }
              
                }
              }
             
            }
          
          }
        }
        
      }
    }
  }
  
  
    /*
      DLIList<DLIList<CubitFacet*>*> shells;
      shells.append(&global_facet_list);
      CubitBoolean dummy;
      FacetDataUtil::stitch_facets(shells, OCTREE_EPSILON, dummy);
      if (shells.size() == 1)
      {
      FacetDataUtil::get_points(global_facet_list, point_list);
      for( j = 0; j < point_list.size(); j++ )
      {
      ptr_point = point_list.get_and_step();
      facet_point_data = new OctreeFacetPointData( ptr_point->coordinates(), ptr_point );
      facet_point_data_list.push( facet_point_data );
      }
      }
      else
      {
      PRINT_INFO("Stitch resulted in multiple shells!!!!!!!!!!\n");
    
      }
    */ 
  
 

    // Reset the visit to CUBIT_FALSE
  reset_td_octree_ref_face_visit( CUBIT_FALSE );
  reset_td_octree_ref_edge_visit( CUBIT_FALSE );

  facet_point_data_list.reset();
  for( i = 0; i < facet_point_data_list.size(); i++ ){
    facet_point_data =  facet_point_data_list.get_and_step();
      //SkeletonDebug::draw_point( facet_point_data->coordinates(), CUBIT_MAGENTA_INDEX );
    cubitOctree->subdivide_octree_based_on_facet_point( facet_point_data, cubitOctree->get_max_depth() );
  }

    // refine the cells on edges with large enough dihedral angles
    // disabling for now 
    //int target_depth = skeletonProxy->get_source_entity_max_depth();
    //cubitOctree->refine_cells_to_target_depth(refine_cell_list, target_depth);

  return CUBIT_TRUE;
}



CubitBoolean CubitOctreeGeneratorVolumes::find_intersection_between_octree_and_facets( DLIList<CubitOctreeNode *> &queue_for_mat_generation ){
  DLIList< RefFace *> vol_ref_faces;
  RefFace *ptr_ref_face;
  TDOctreeRefFace *td_ref_face;
  DLIList< CubitFacet *> *facet_list;
  CubitFacet *ptr_facet;
  DLIList<CubitOctreeCell *> CubitOctree_cell_list; 
  DLIList< CubitOctreeNode *> CubitOctree_grid_node_list;
  DLIList<CubitOctreeNode*> all_grid_nodes;
  
    //CubitVector facet_normal;
    //DLIList< CubitPoint *> *point_list;
  int i, j, k;
  int index;
  Body *volume;

//  CpuTimer total_idata_time;
  
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_back(), Body );  
    vol_ref_faces.clean_out();
    volume->ref_faces( vol_ref_faces );
      //vol_ref_faces.reset();
      //int cell_count = 0;
      //CubitBox facet_bbox;
      //PRINT_INFO(" \nRef Volume  = %d ", volume->id());
        
    for( i = 0; i < vol_ref_faces.size(); i++ ){
      ptr_ref_face = vol_ref_faces.get_and_step();
      td_ref_face = TDOctreeRefFace::get_td( ptr_ref_face );
      
      if(td_ref_face->get_visit() == CUBIT_FALSE){
        td_ref_face->set_visit( CUBIT_TRUE );

          //PRINT_INFO(" \nRef Face = %d ", ptr_ref_face->id());

        RefFace *current_ref_face = CAST_TO(ptr_ref_face, RefFace);
          /*if (current_ref_face->sense(CAST_TO(volume,Body)) == CUBIT_REVERSED)
            {
            temp_facet_normal *= -1.0;
            }*/
        
//        PRINT_INFO("Number of parent ref Entities for face %d: %d\n", current_ref_face->id(), current_ref_face->num_parent_ref_entities());

          // get the sense so we know which way to propagate the 3DMAT front inwards
        CubitSense ref_face_sense = CUBIT_FORWARD; //current_ref_face->sense(volume);
        
          // in the case of a merged face, figure out if two parent vols share the same
          // sizing function. If so, propagate in both directions
        if (current_ref_face->is_merged() /*&& current_ref_face->num_parent_ref_entities() == 2*/)
        {
          DLIList<RefEntity*> parent_ents;
          current_ref_face->get_parent_ref_entities(parent_ents);
          if (parent_ents.size() == 2)
          {
            Body *other_parent_vol = CAST_TO(parent_ents.get_and_step(),Body);
            if (other_parent_vol == volume) {other_parent_vol = CAST_TO(parent_ents.get(),Body);}
            if (other_parent_vol != NULL)
            {
              
                      // if the other parent volume has the same sizing function,
                      // pass CUBIT_UNKNOWN to CubitOctreeCell::set_color_and_idatas
                      // It will then set up the front normals to propagate in both directions
                      ref_face_sense = CUBIT_UNKNOWN;
              
            }
          }
        }
        
        facet_list = td_ref_face->get_ptr_cubit_facet_list();
          //PRINT_INFO("Surface %d: facets: %d\n", ptr_ref_face->id(), facet_list->size());
        
        for( j = 0; j < facet_list->size(); j++ ){
          ptr_facet = facet_list->get_and_step();

            //CubitVector temp_facet_normal = ptr_facet->normal();
        
          
            //SVDrawTool::draw_vector(ptr_facet->center(), ptr_facet->center()+temp_facet_normal, (i%8)+1);
          

//#ifndef NDEBUG
//          SkeletonDebug::draw_facet( ptr_facet, CUBIT_BLUE_INDEX );
//#endif
            //facet_bbox = ptr_facet->bounding_box();
            //draw_box( facet_bbox );
            //facet_normal = ptr_facet->normal();
            //if( ptr_facet->is_backwards() ){
            // facet_normal = facet_normal * -1;
            //}
            //GfxDebug::draw_line( ptr_facet->center().x(), ptr_facet->center().y(), ptr_facet->center().z(), ptr_facet->center().x()+facet_normal.x(), ptr_facet->center().y()+facet_normal.y(), ptr_facet->center().z()+facet_normal.z(), CUBIT_WHITE_INDEX );       
          CubitOctree_cell_list.clean_out();
          CubitOctree_grid_node_list.clean_out();
            //PRINT_INFO(" Before find_octree_cells_contatined \n");
          cubitOctree->find_octree_cells_contained_inside_bbox( ptr_facet, CubitOctree_cell_list );
            //PRINT_INFO(" Before mark_positive_and_negative_grid_nodes \n");
          cubitOctree->mark_positive_and_negative_octree_grid_nodes( ptr_facet, CubitOctree_cell_list, CubitOctree_grid_node_list  );
            //PRINT_INFO(" Before find_intersection_between_grid_edge_facet \n");

          if ( OCTREE_DEFAULT_INTERSECTION_METHOD[OCTREE_TIME_ACCURACY_LEVEL] == SAT_INTERSECT){
              //PRINT_INFO("%d cells in bbox\n", CubitOctree_cell_list.size());
            
            for (k=0; k < CubitOctree_cell_list.size(); ++k)
            {
              CubitOctreeCell *current_cell = CubitOctree_cell_list.get_and_step();
                // if (current_cell->get_mark() == CUBIT_FALSE) {continue;}
                //current_cell->set_mark(CUBIT_FALSE);
          
              if (current_cell->does_facet_intersect_octreecell(ptr_facet))
              {
                  /*  double corners[3];
                      double half_edge_length = current_cell->get_dimension()/2.0;
                      CubitVector center = current_cell->get_center();
                      center.get_xyz(corners);
                      float box[6] = {corners[0]-half_edge_length, corners[1]-half_edge_length, corners[2]-half_edge_length,
                      corners[0]+half_edge_length, corners[1]+half_edge_length, corners[2]+half_edge_length};
                      SVDrawTool::draw_cube(box, CUBIT_GREEN_INDEX, SVDrawTool::WIRE);*/

                current_cell->set_color_and_intersection_datas(ptr_facet, ptr_ref_face,
#ifdef USE_octree_BGMESH
                                                               cubitOctree->get_greycelllist(),
#endif
                                                               ref_face_sense);
                
                  //cubitOctree->get_greycelllist()->append(current_cell);
              }
                /*if (ref_face_sense) {
                  //  SVDrawTool::mouse_xforms();
                  }*/
              
            }

            for( k = 0; k < CubitOctree_grid_node_list.size(); k++ )
            {
              CubitOctreeNode *ptr_grid_node = CubitOctree_grid_node_list.get_and_step();
              if (ptr_grid_node->get_visit() == CUBIT_FALSE)
              {
                all_grid_nodes.append(ptr_grid_node);
                ptr_grid_node->set_visit(CUBIT_TRUE);
              }
          
              ptr_grid_node->set_mark( CUBIT_FALSE );	
              if( ptr_grid_node->get_halfspace_direction() == OCTREE_POSITIVE )
              {
                ptr_grid_node->set_halfspace_direction( OCTREE_NEGATIVE );
              }
            }
          }

          if (OCTREE_DEFAULT_INTERSECTION_METHOD[OCTREE_TIME_ACCURACY_LEVEL] == GRID_EDGE_INT)
          {
            cubitOctree->find_intersection_between_grid_edges_and_facet( CUBIT_OCTREE_VOLUME,  ptr_ref_face, ptr_facet, CubitOctree_grid_node_list );
          }
          

          
        }
          
      }
    }
  }
  
    // Reset the visit to CUBIT_FALSE
  reset_td_octree_ref_face_visit( CUBIT_FALSE );
    //SVDrawTool::mouse_xforms();
    // PRINT_INFO("profiling: time for intersecting and making idatas: %f\n", total_idata_time.cpu_secs());
  
  if (OCTREE_DEFAULT_INTERSECTION_METHOD[OCTREE_TIME_ACCURACY_LEVEL] == SAT_INTERSECT)
  {
      //CpuTimer front_color_time;
    for (i=0; i < all_grid_nodes.size(); ++i)
    {
      CubitOctreeNode *grid_node = all_grid_nodes.get_and_step();
      DLIList<OctreeIntersectionData*> *idata_list = grid_node->get_idata_list();
        //idata_list->uniquify_unordered();
      
      
      if (idata_list->size() > 0)
      {
        DLIList<OctreeIntersectionData*> closest_ones;
        //OctreeIntersectionData *closest = NULL;
        double min_list_dist = CUBIT_DBL_MAX;
        for (k=0; k < idata_list->size(); ++k)
        {
          if (idata_list->get()->is_merged()) {idata_list->step();}
          
          double temp_dist = idata_list->get()->get_length();
          if (temp_dist < min_list_dist - OCTREE_EPSILON)
          {
            closest_ones.clean_out();
//          closest = idata_list->get();
            closest_ones.append(idata_list->get());
            min_list_dist = temp_dist;
          }
          else if (fabs(temp_dist - min_list_dist) <= OCTREE_EPSILON)
          {
            closest_ones.append(idata_list->get());
          }
          idata_list->step();
        }
      
        bool white = false;

        for (j=0; j < closest_ones.size(); ++j)
        {
          if (closest_ones.get_and_step()->get_halfspace())
          {
            white = true;
            break;
          }
        }

        //closest = closest_ones.get();
//        if (!white) {SVDrawTool::draw_vector(grid_node->get_coord(), closest->get_int_point(), CUBIT_RED_INDEX);}
        
          // this is the closest facet encountered by the current grid node so far
          // update color of node based on this facet
        if (white)
        {
          grid_node->set_color(CUBIT_WHITE_INDEX);
          cubitOctree->get_whitenodelist()->append(grid_node);
          
        }
        else if (!white)
        {
//           if ( (grid_node->get_coord()-closest->get_facet_ptr()->center()) % (closest->get_normal()) < 0)
//           {
//             PRINT_INFO("Adding a white node to queue !!!!!!!!!!!!\n");
//           }

            // uncomment this
            //             for (j=0; j < grid_node->get_idata_list()->size(); ++j)
//           {
//             OctreeIntersectionData *idata = grid_node->get_idata_list()->get_and_step();
//             double proj_length = fabs((grid_node->get_coord() - idata->get_int_point())%idata->get_facet_normal());
//             idata->set_length(proj_length);
//             }
          
//           SVDrawTool::draw_vector(grid_node->get_coord(), closest->get_int_point(), CUBIT_RED_INDEX);
          
          grid_node->set_color(CUBIT_BLACK_INDEX);
          queue_for_mat_generation.append(grid_node);
        }
        else
        {
            //PRINT_INFO("Couldn't pick halfspace to color node!\n");
        }
      }
      grid_node->set_visit(CUBIT_FALSE);
    }
    
    for (i=0; i < cubitOctree->get_whitenodelist()->size(); ++i)
    {
      CubitOctreeNode *ptr_grid_node = cubitOctree->get_whitenodelist()->get_and_step();
      DLIList<OctreeIntersectionData*> *idata_list = ptr_grid_node->get_idata_list();
      for (k=0; k < idata_list->size(); ++k)
      {
        OctreeIntersectionData *idata = idata_list->get_and_step();
        if (idata != NULL) {delete idata;}
      }
      idata_list->clean_out();
      CubitOctreeNode *adj_node = NULL;
      for (k=0; k < 6; ++k)
      {
        adj_node = ptr_grid_node->get_adj_node(k);
        if (adj_node != NULL)
        {
          if (adj_node->get_color() != CUBIT_BLACK_INDEX && adj_node->get_color() != CUBIT_WHITE_INDEX)
          {
            adj_node->set_color(CUBIT_YELLOW_INDEX);
          }
            //SVDrawTool::draw_point(adj_node->get_coord(), adj_node->get_color());
        }
      }
    }
      // PRINT_INFO("profiling: time to color front nodes: %f\n", front_color_time.cpu_secs());
    
  }
  
 
  
    

    // PRINT_INFO("Testing: Cell Count = %d", cell_count );
  if (OCTREE_DEFAULT_INTERSECTION_METHOD[OCTREE_TIME_ACCURACY_LEVEL] == SAT_INTERSECT)
  {

   //     CpuTimer front_init;
    
    for( i = 0; i < queue_for_mat_generation.size(); i++ ){
      queue_for_mat_generation.get_and_step()->SAT_find_face_distance_average_normal( /*skeletonSizingFunction, queue_for_mat_generation*/ );    
    }
      //PRINT_INFO("profiling: time to init front: %f\n", front_init.cpu_secs());
       
  }

  
  return CUBIT_TRUE;
}



void CubitOctreeGeneratorVolumes::find_optimal_min_and_max_depth( RefEntity*entity, int &local_min_depth, int &local_max_depth )
{
  Body *volume = dynamic_cast<Body *> (entity);
  if( NULL == volume)
    return;
  
  CubitBox bbox = volume->bounding_box();
  
  // calculate min_depth by comparing min_range with max_range in log scale
  double max_range = std::max( std::max( bbox.x_range(), bbox.y_range()), bbox.z_range() );
  double min_range = std::min( std::min( bbox.x_range(), bbox.y_range()), bbox.z_range() );
  local_min_depth = (int)(log( max_range / min_range ) / log(2.0)) + 1;  // add plus 1 just to give enough one extra level of depth

  if( local_min_depth < MIN_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL] )
    local_min_depth = MIN_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL];
  if( local_min_depth > MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL] )
    local_min_depth = MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL];

  // calculate max_depth by comparing the min curve length with max_range in log scale
  DLIList< RefEdge *> list_edges;
  volume->ref_edges( list_edges );

  int i;
  RefEdge *ptr_edge;
  double min_edge_len = DBL_MAX;
  double len;
  CubitVector start_pt, mid_pt, end_pt;
  for( i = 0; i < list_edges.size(); i++ )
  {
    ptr_edge = list_edges.get_and_step();
    //len = ptr_edge->length_from_u( ptr_edge->start_param(), ptr_edge->end_param() );
    start_pt = ptr_edge->start_coordinates();
    end_pt = ptr_edge->end_coordinates();
    ptr_edge->mid_point(mid_pt);
    len = (mid_pt - start_pt).length() + (end_pt - mid_pt).length();
    if( len < min_edge_len )
    {
      min_edge_len = len;
    }
  }
  local_max_depth = (int)(log( max_range / min_edge_len) / log(2.0) ) + 1;

 /*
  // calculate max_depth by comparing the min mesh size with max_range in log scale
  double min_mesh_size = max_range;
  double mesh_size;
  for( i = 0; i < list_edges.size(); i++ )
  {
    ptr_edge = list_edges.get_and_step();
   
    mesh_size = max_range;
    MeshToolProxy* mesh_tool_ptr = MeshToolProxy::mesh_tool(ptr_edge);
    SizeIntervalType hardness;
    double tmp_size = mesh_tool_ptr->requested_interval_size(hardness);
    if (hardness > CALCULATED)
    {
      mesh_size = tmp_size/ 2.0;
    }
    if( mesh_size < min_mesh_size )
    {
      min_mesh_size = mesh_size;
    }
  }
  local_max_depth = std::max( local_max_depth, (int)(log( max_range / min_mesh_size) / log(2.0) ) + 1 );
*/
/*
   // calculate max_depth by comparing the min mesh size with max_range in log scale
  DLIList< RefFace *> list_faces;
  volume->ref_faces( list_faces );
  min_mesh_size = max_range;
  RefFace *ptr_face;
  for( i = 0; i < list_faces.size(); i++ )
  {
    ptr_face = list_faces.get_and_step();
   
    mesh_size = max_range;
    MeshToolProxy* mesh_tool_ptr = MeshToolProxy::mesh_tool(ptr_face);
    SizeIntervalType hardness;
    double tmp_size = mesh_tool_ptr->requested_interval_size(hardness);
    if (hardness > CALCULATED)
    {
      mesh_size = tmp_size/2.0;
    }
    if( mesh_size < min_mesh_size )
    {
      min_mesh_size = mesh_size;
    }
  }
  local_max_depth = std::max( local_max_depth, (int)(log( max_range / min_mesh_size) / log(2.0) ) + 1 );
*/
  
  if( local_max_depth > MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL] )
    local_max_depth = MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL];

  // make sure local_max_depth is atleast one depth higher than local_min_depth
  if( local_max_depth <= local_min_depth )
    local_max_depth = local_min_depth + 1;

  

}

// finds the optimal min and max depths of CubitOctree 
void CubitOctreeGeneratorVolumes::find_optimal_min_and_max_octree_depths( DLIList< RefEntity *> &entity_list, int &min_depth, int &max_depth )
{
  DLIList< int > list_min_depth;
  DLIList< int > list_max_depth;

  // traverse through the list of volumes and find optimal depth
  Body*volume;
  int index;
  int local_min_depth, local_max_depth;
  for( index = 0; index < entity_list.size(); index++ )
  {
    volume =  CAST_TO( entity_list.get_and_step(), Body );
    if( volume )
    {      
      CubitOctreeGeneratorVolumes::find_optimal_min_and_max_depth( volume, local_min_depth, local_max_depth );
      list_min_depth.append( local_min_depth );
      list_max_depth.append( local_max_depth );
    }
  }

  // find resultant min_depth and max_depth
  // for now let us take the floor of the average depth
  int avg_min_depth = 0;
  int avg_max_depth = 0;
  for( index = 0; index < list_min_depth.size(); index++ )
  {
    avg_min_depth += list_min_depth[index];
    avg_max_depth += list_max_depth[index];
  }
  avg_min_depth /= list_min_depth.size();
  avg_max_depth /= list_max_depth.size();

  // assign min of initial value and avg depths as final answer
  min_depth = std::min( min_depth, avg_min_depth );
  max_depth = std::min( max_depth, avg_max_depth );

}

void CubitOctreeGeneratorVolumes::reset_td_octree_ref_face_visit( const CubitBoolean type ){
  int index, i;
  Body *volume;
  DLIList<RefFace *> vol_ref_faces;
  
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_ref_faces.clean_out();
    volume->ref_faces( vol_ref_faces );
    for( i = 0; i < vol_ref_faces.size(); i++ ){
      TDOctreeRefFace::get_td( vol_ref_faces.get_and_step() )->set_visit( type );      
    }
  }
  
  return;
}


void CubitOctreeGeneratorVolumes::reset_td_octree_ref_edge_visit( const CubitBoolean type ){
  int index, i;
  Body *volume;
  DLIList<RefEdge *> vol_ref_edges;
  
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_ref_edges.clean_out();
    volume->ref_edges( vol_ref_edges );
    for( i = 0; i < vol_ref_edges.size(); i++ ){
      TDOctreeRefEdge::get_td( vol_ref_edges.get_and_step() )->set_visit( type );      
    }
  }
  
  return;
}

CubitBoolean CubitOctreeGeneratorVolumes::build_td_mref_faces( void ){
  DLIList< RefFace *> vol_ref_faces;
  RefFace *ptr_ref_face;
  TDOctreeRefFace *td_mref_face;
  DLIList<CubitFacet *> *facet_list;
  DLIList<CubitFacetEdge *> *facet_edge_list;
  DLIList<CubitPoint*> *point_list;
  int i;
  int index;
  Body *volume;
  
  CubitBoolean identified_engine = CUBIT_FALSE;
  CubitBoolean is_facet_geom = CUBIT_FALSE;
  
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_ref_faces.clean_out();
    volume->ref_faces( vol_ref_faces );
    identified_engine = CUBIT_FALSE;
    for( i = 0; i < vol_ref_faces.size(); i++ ){
      ptr_ref_face = vol_ref_faces.get_and_step();
      if( ptr_ref_face->marked() == CUBIT_FALSE ){
        if (identified_engine == CUBIT_FALSE)
        {
          if (ptr_ref_face->get_geometry_query_engine()->get_engine_version_string().find("Facet") != CubitString::npos)
          {
            is_facet_geom = CUBIT_TRUE;
          }
          identified_engine = CUBIT_TRUE;
        }
        
        
        ptr_ref_face->marked( CUBIT_TRUE );
        
        // Add TD for MRefFace
        TDOctreeRefFace::add_td(ptr_ref_face);
        facet_list = new DLIList<CubitFacet *>;
        facet_edge_list = new DLIList<CubitFacetEdge *>;		
        point_list	= new DLIList<CubitPoint*>;
    		
        //PRINT_DEBUG_157(" \n\nRef Face = %d \n ", ptr_ref_face->id());
        //facet_list->clean_out();
        //point_list->clean_out();
        //facet_edge_list->clean_out();
        
        if (is_facet_geom)
        {
          DLIList<CubitFacet*> facets;
          DLIList<CubitPoint*> points;
          
          FacetSurface *facet_surf_ptr = CAST_TO(ptr_ref_face->get_surface_ptr(),FacetSurface);
          if (facet_surf_ptr != NULL)
          {
            facet_surf_ptr->get_my_facets(facets,points);
          }
          FacetDataUtil::copy_facets(facets,*facet_list,*point_list,*facet_edge_list);
          /*PRINT_INFO("Grabbed %d facets from surface %d\n", facet_list->size(),ptr_ref_face);
           PRINT_INFO("Grabbed %d facet edges from surface %d\n", facet_edge_list->size(),ptr_ref_face);
           PRINT_INFO("Grabbed %d facet points from surface %d\n", point_list->size(),ptr_ref_face);*/
        }
        else
        {
          GMem gmem;
          //Surface* surf_ptr = ptr_ref_face->get_surface_ptr();
          //surf_ptr->get_geometry_query_engine()->get_graphics(surf_ptr, &gmem);
          int angle = 15;
          float distance = 0;
          //SVDrawTool::get_facet_tolerance(angle, distance);  
          angle = static_cast<int>(ANG_FACET_EXTRACT[OCTREE_TIME_ACCURACY_LEVEL]);
          //PRINT_INFO("\n Facets extracted with angle %d and distance %f \n", angle, distance );
          ptr_ref_face->get_graphics( gmem, angle, distance, 0 );
          ChollaEngine::get_facets(gmem, *facet_list, *point_list);
          
          /*
           // delete point_list;
           // point_list = new DLIList<CubitPoint*>;
           point_list->clean_out();
           
           DLIList<DLIList<CubitFacet *> *> shell_list;
           
           CubitBoolean dummy2;
           double tol = SKL_EPSILON;
           DLIList <CubitQuadFacet *> qfacet_list;
           FacetDataUtil::split_into_shells(*facet_list, qfacet_list,shell_list, dummy2);
           
           FacetDataUtil::stitch_facets(shell_list,tol,dummy2,false);
           //int nv,ne;
           //DLIList<CubitPoint*> un;
           
           //FacetDataUtil::merge_coincident_vertices(shell_list,tol,nv,ne,un);
           facet_list = shell_list.get();
           PRINT_INFO("skl: number of shells = %d\n", shell_list.size());
           
           FacetDataUtil::get_edges( *facet_list, *facet_edge_list );
           FacetDataUtil::get_points(*facet_list,*point_list);
           SVDrawTool::draw_facets(*facet_list,CUBIT_GREEN_INDEX);
           */  
          
          FacetDataUtil::get_edges( *facet_list, *facet_edge_list );
          /*
           // checks for problems with gfx facets
           // debug use only
           int v;
           for (v=0; v < facet_edge_list->size(); ++v)
           {
           CubitFacetEdge *temp_edge = facet_edge_list->get_and_step();
           if (temp_edge->num_adj_facets() != 2 && temp_edge->num_adj_facets() != 1)
           {
           PRINT_INFO("Skeleton sizing function, copying gfx facets: edge has %d adj facets!\n", temp_edge->num_adj_facets());
           SVDrawTool::draw_vector(temp_edge->point(0)->coordinates(), temp_edge->point(1)->coordinates(), CUBIT_RED_INDEX);
           SVDrawTool::mouse_xforms();
           }
           }
           int w;
           for (v=0; v < facet_list->size(); ++v)
           {
           CubitFacet *first = (*facet_list)[v];
           for (w=v+1; w < facet_list->size(); ++w)
           {
           CubitFacet *second = (*facet_list)[w];
           if (first == second)
           {
           PRINT_INFO("Skeleton Sizing Function: gfx facet list contains a duplicate, ids are %d,%d\n", first->id(), second->id());
           }
           }
           }
           */
        }
        
        td_mref_face = TDOctreeRefFace::get_td(ptr_ref_face);
        
        /*
         td_mref_face->set_ptr_cubit_facet_list(facet_list);
         td_mref_face->set_ptr_cubit_facet_edge_list(facet_edge_list);
         td_mref_face->set_ptr_cubit_point_list(point_list);
         */
        if( td_mref_face->get_ptr_cubit_facet_list() )
        {
          delete td_mref_face->get_ptr_cubit_facet_list(); 
          td_mref_face->set_ptr_cubit_facet_list(facet_list);
          
        }
        else
        {
          td_mref_face->set_ptr_cubit_facet_list(facet_list);
        }
        
        if( td_mref_face->get_ptr_cubit_facet_edge_list() )
        {
          delete td_mref_face->get_ptr_cubit_facet_edge_list();
          td_mref_face->set_ptr_cubit_facet_edge_list(facet_edge_list);
          
        }
        else
        {
          td_mref_face->set_ptr_cubit_facet_edge_list(facet_edge_list);
        }
        
        if( td_mref_face->get_ptr_cubit_point_list() )
        {
          delete td_mref_face->get_ptr_cubit_point_list();
          td_mref_face->set_ptr_cubit_point_list(point_list);
        }
        else
        {
          td_mref_face->set_ptr_cubit_point_list(point_list);
        }
        
        // check for bad facets and disable 2dmat if necessary
        td_mref_face->check_valid_facets(CUBIT_TRUE);
        //if (!create_2dmat) {td_mref_face->set_create_2dmat(CUBIT_FALSE);}
      }    
    }
  }
  
  // Reset the marker of faces
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_ref_faces.clean_out();  
    volume->ref_faces( vol_ref_faces );
    for( i = 0; i < vol_ref_faces.size(); i++ ){
      vol_ref_faces.get_and_step()->marked( CUBIT_FALSE );
    }
  }
  
  return CUBIT_TRUE;
}


CubitBoolean CubitOctreeGeneratorVolumes::build_td_mref_edges( void )
{
  DLIList< RefEdge *> vol_mref_edges;
  RefEdge *ptr_mref_edge;
  TDOctreeRefEdge *td_mref_edge;
  int i;
  int index;
  Body *volume;
  
  // Reset the marker of edges
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_mref_edges.clean_out();
    volume->ref_edges( vol_mref_edges );
    for( i = 0; i < vol_mref_edges.size(); i++ ){
      ptr_mref_edge = vol_mref_edges.get_and_step();
      if( ptr_mref_edge->marked() == CUBIT_FALSE ){
        ptr_mref_edge->marked( CUBIT_TRUE );
        // Add TD for MRefEdge
        TDOctreeRefEdge::add_td(ptr_mref_edge);
        td_mref_edge = TDOctreeRefEdge::get_td(ptr_mref_edge);
        td_mref_edge->generate_gpoint_list(ANG_FACET_EXTRACT[OCTREE_TIME_ACCURACY_LEVEL] * OCTREE_TOLERANCE_FOR_CURVE_DECIMATION_FACTOR );
      }
    }
  }
  
  // Reset the marker of edges
  for( index = 0; index < entityList.size(); index++ ){
    volume = CAST_TO( entityList.get_and_step(), Body );
    vol_mref_edges.clean_out();
    volume->ref_edges( vol_mref_edges );
    for( i = 0; i < vol_mref_edges.size();i++ ){
      vol_mref_edges.get_and_step()->marked( CUBIT_FALSE );
    }
  }
  
  return CUBIT_TRUE;
}



void CubitOctreeGeneratorVolumes::color_octreenode_via_grassfire( DLIList<CubitOctreeNode *> &queue_for_mat_generation )
{
  
  PriorityQueue<CubitOctreeNode *> heap( CubitOctreeNode::compare_function);
  CubitOctreeNode *ptr_grid_node;
  
  // ************************* For Testing Purpose ***************************
  // int i;
  //for( i = 0; i < queue_for_mat_generation.size(); i++ ){
  //ptr_grid_node = queue_for_mat_generation.get_and_step();
  
  //PRINT_DEBUG_157(" boundary node distance = %f \n",ptr_grid_node->get_distance() );
  //PRINT_DEBUG_157(" boundary node visit = %d \n",ptr_grid_node->get_visit() );
  
  //PRINT_DEBUG_157(" mat initial grid = %f %f %f %f \n", ptr_grid_node->get_distance(), ptr_grid_node->get_coord().x(), ptr_grid_node->get_coord().y(), ptr_grid_node->get_coord().z() );
  //GfxDebug::draw_point( ptr_grid_node->get_coord().x(), ptr_grid_node->get_coord().y(), ptr_grid_node->get_coord().z(), CUBIT_YELLOW_INDEX );
  //ptr_grid_node->display( this, DISTANCE );
  //}
  // ************************* For Testing Purpose ***************************
  
  //PRINT_DEBUG_157(" Num of initial points for mat generation = %d \n", queue_for_mat_generation.size() );
  
  // Initialize the heap
  while( queue_for_mat_generation.size() > 0 ){
    heap.push( queue_for_mat_generation.pop() );
  }
  
  // distance by default is CUBIT_DBL_MAX, and the white boundary nodes's distance is set to -1.
  //PRINT_DEBUG_157(" heap size = %d", heap.size() );
  //int draw_count = 1;
  
  
  while( heap.size() > 0 ){
    //        PRINT_INFO( " manhattan in process before opo() = %d \n",queue.size() );
    ptr_grid_node = heap.top();
    heap.pop();
    

    ptr_grid_node->find_distance_at_adj_node( &heap );
    
       //		PRINT_INFO( " manhattan in process after NumberAdjNode() = %d \n", queue.size() );
  }
  
}

CubitStatus CubitOctreeGeneratorVolumes::generate_full_octree( void )
{
  // all the functions of skeleton sizing function are called here
  CpuTimer octgentime;
  double time_octree_generation;
 	
  //ProgressTool *progressToolPtr = AppUtil::instance()->progress_tool();
  //assert(progressToolPtr != NULL );
  //char title[42]="Octree Generation in Progress ...";
  //progressToolPtr->start(0, 100, title, NULL, CUBIT_TRUE, CUBIT_TRUE);
	
  //PRINT_DEBUG_157("\nENTER: Octree Generation \n");

  
  CpuTimer preproctime;
  
    build_td_mref_faces();
    build_td_mref_edges();
  
  
  // this is used to store size and also used in 3D MAT generation
  CpuTimer gen_latt;

  int min_depth = MIN_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL];
  int max_depth = MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL];
  CubitOctreeGeneratorVolumes::find_optimal_min_and_max_octree_depths( entityList, min_depth, max_depth );
  get_octree_lattice()->set_max_depth(std::min( get_octree_lattice()->get_max_depth(), max_depth) );
  get_octree_lattice()->set_min_depth(std::max( get_octree_lattice()->get_min_depth(), min_depth) );
  //PRINT_INFO("\nOctree min_depth = %d, max_depth = %d\n", min_depth, max_depth);
  
  generate_lattice();
  //PRINT_INFO("profiling: time in generate_lattice: %f\n", gen_latt.cpu_secs());
  
 // PRINT_DEBUG_157("AFTER: octree lattice generation\n");
  
    
  //PRINT_DEBUG_157("BEFORE smooth transiton number of nodes = %d\n",gridNodeVector.size() );
  CpuTimer time_to_balance;
  get_octree_lattice()->establish_smooth_transition_of_cells( max_depth ); //OCTREE_MAX_DEPTH_OCTREE[OCTREE_TIME_ACCURACY_LEVEL] );
 // PRINT_INFO("profiling: time to balance octree: %f\n", time_to_balance.cpu_secs());
  
//  PRINT_DEBUG_157("AFTER smooth transiton \n" );
  
  // Finds the intersection between the facets and the octree
  CpuTimer intersection_time;
  DLIList<CubitOctreeNode *> queue_for_3Dmat_generation;
 find_intersection_between_octree_and_facets( queue_for_3Dmat_generation );
  //PRINT_INFO("profiling: intersection time: %f\n", intersection_time.cpu_secs());
  
  //PRINT_DEBUG_157("AFTER intersection between facets and octree \n" );
  
  
  //PRINT_DEBUG_157("\nDONE: Octree Generation \n" );
  //progressToolPtr->percent(1.0);

  if( AppUtil::instance()->interrupt() ){
    //CubitBoolean interruptStatus;
    PRINT_WARNING("\n Killed  octree generation...\n");
    //interruptStatus = CUBIT_TRUE;
    //progressToolPtr->end();
    return CUBIT_FAILURE;
  }	
  
  color_octreenode_via_grassfire( queue_for_3Dmat_generation );

  color_lattice_cell();
  
  
  time_octree_generation = octgentime.cpu_secs();
  
  
  
  
  PRINT_INFO("Total time for octree generation: %f\n", time_octree_generation);
  
  return CUBIT_SUCCESS;
}


// OPTIMIZE: priorityqueue -> dlist/queue
// EOF

