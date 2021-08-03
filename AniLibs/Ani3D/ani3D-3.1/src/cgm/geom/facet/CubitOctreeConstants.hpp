//-------------------------------------------------------------------------
// Filename      : CubitOctreeConstants.hpp
//
// Purpose       : Constants for the octree generation
//
// Creator       : William Roshan Quadros
//
// Creation Date : 9/2/2013
//
// Owner         : 
//-------------------------------------------------------------------------

#ifndef CUBIT_OCTREE_CONSTANTS_H 
#define CUBIT_OCTREE_CONSTANTS_H



#include "CubitDefines.h"

const int OCTREE_X = 0;
const int OCTREE_Y = 1;
const int OCTREE_Z = 2;
//- Cartician coordinate axis


//- during triliear interpolation if number of zero size nodes exceeds MAX_NUM_ZERO_SIZE_NODES
//- then average size is taken instead of zero.
const int MAX_NUM_ZERO_SIZE_NODES = 2;

const int CUTOFF_OCTREE_DEPTH = 3;
const int CUTOFF_BBOX_RANGE_RATIO = 17;
const int CUTOFF_NUM_NODES = 250;

enum CubitOctreeType{ CUBIT_OCTREE_VOLUME, CUBIT_OCTREE_FACE, CUBIT_OCTREE_EDGE };
// Depth settings for octree 
const int MIN_DEPTH_OCTREE[4] = { 0, 2, 5, 5};
const int MAX_DEPTH_OCTREE[4] = { 0, 7, 7, 7};
const int INCR_DEPTH_OCTREE_FOR_SKELETON = 0;  
const int OCTREE_MIN_DEPTH_OCTREE[4] = { 0, MIN_DEPTH_OCTREE[1] + INCR_DEPTH_OCTREE_FOR_SKELETON,MIN_DEPTH_OCTREE[2] + INCR_DEPTH_OCTREE_FOR_SKELETON, MIN_DEPTH_OCTREE[3] + INCR_DEPTH_OCTREE_FOR_SKELETON};
const int OCTREE_MAX_DEPTH_OCTREE[4] = { 0, MAX_DEPTH_OCTREE[1] + INCR_DEPTH_OCTREE_FOR_SKELETON,MAX_DEPTH_OCTREE[2] + INCR_DEPTH_OCTREE_FOR_SKELETON, MAX_DEPTH_OCTREE[3] + INCR_DEPTH_OCTREE_FOR_SKELETON};

enum OctreePosition { O_UNSET = -1, O_LEFT, O_RIGHT, O_BOTTOM, O_TOP, O_BACK, O_FRONT };


enum OctreeNodeConstant {NODE_SIZE, NODE_FACE_NUM, NODE_DISTANCE,NODE_NORMAL, NODE_TENSOR, NODE_ELEM_TENSOR };
// for display of size, face num, and distance at grid node

const double OCTREE_TOLERANCE_FOR_CURVE_DECIMATION_FACTOR = 0.5;

//- Default level of time-accuracy is set to 2 <1 to 3>
const int OCTREE_TIME_ACCURACY_LEVEL = 2;
// - This angle is used to insert virtual facetpointdata while generating PR-Octree
const double FACET_SPLITTING_INTERNAL_ANGLE = 10.0 * CUBIT_PI / 180.0;
//- graphics facets are splitted before generating the chordal axis.
//- after facets' classification, some class of facets are split only 
//- if the internal angle is large.  The new facets will contain internal angle
//-  approximately equal to FACET_SPLITTING_INTERNAL_ANGLE


enum OctreeSourceEntityType{ OCTREE_SIZE_DEFAULT=0 };

const double OCTREE_EPSILON = 0.000001;
const CubitBoolean OCTREE_POSITIVE = CUBIT_TRUE;
const CubitBoolean OCTREE_NEGATIVE = CUBIT_FALSE;

enum InterpolationMode{ INVERSE_SIZE, INVERSE_DISTANCE, INVERSE_SIZE_DISTANCE, MIN_DISTANCE, MIN_SIZE, AVERAGE };
//- Schemes for interpolation of sizing function 

enum InterpolationType{ INVERSE_LINEAR, INVERSE_QUADRATIC, INVERSE_CUBIC };
//- Options for inverse distance interpolation.
const InterpolationType DEFAULT_INVERSE_DISTANCE_SCHEME_INSIDE_CELL = INVERSE_LINEAR;

enum InterpolationModeInsideCell{ TRILINEAR_INTERPOLATION_INSIDE_CELL, INVERSE_DISTANCE_INTERPOLATION_INSIDE_CELL, MIN_DISTANCE_INTERPOLATION_INSIDE_CELL, MIN_SIZE_INTERPOLATION_INSIDE_CELL, INVERSE_SQUARE_DISTANCE_INTERPOLATION_INSIDE_CELL };
const InterpolationModeInsideCell DEFAULT_INTERPOLATION_INSIDE_CELL = INVERSE_SQUARE_DISTANCE_INTERPOLATION_INSIDE_CELL;



enum SizeType{ MESH_SIZE, SCALED_SIZE };
//- Size for meshing and scaled size for display (1-10)

enum DistanceMetric {PROJECTED_DIST, CAPSULE_DIST, MANHATTAN_DIST};


enum IntersectionMethod {GRID_EDGE_INT, SAT_INTERSECT};
//- These are the facet-octree cell intersection methods


const IntersectionMethod OCTREE_DEFAULT_INTERSECTION_METHOD[4] = { (IntersectionMethod)0, SAT_INTERSECT, SAT_INTERSECT, SAT_INTERSECT};   
const DistanceMetric defaultDistanceMetric = CAPSULE_DIST;
const double N_CLOSEST_FACETS_FACTOR_FOR_FRONT_NORMALS = 1.00;


//- Tolerance 
const double ANGLE_BETWEEN_FACETS_NORMAL_FOR_OCTREE_BASED_ON_SURF_CURVATURE_FACTOR = 0.5; 
const double SLENDER_FACET_INTERNAL_ANGLE =  10.0 * CUBIT_PI / 180.0;

const double ANG_FACET_EXTRACT[4]  =  { 15.0, 15.0, 10.0, 5.0 };

#endif

//EOF
