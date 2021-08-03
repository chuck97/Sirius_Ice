//- Class: CubitBox
//- Description: This file defines the CubitBox class.
//- Owner: Greg Sjaardema
//- Checked by:

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "CubitDefines.h"

CubitBox::~CubitBox()
{}

CubitBox::CubitBox() :
    minimum_(0.0, 0.0, 0.0),
    maximum_(0.0, 0.0, 0.0)
{}

CubitBox::CubitBox(const CubitVector &min, const CubitVector &max)
{
  minimum_.set (CUBIT_MIN (min.x(), max.x()),
                CUBIT_MIN (min.y(), max.y()),
                CUBIT_MIN (min.z(), max.z()));
  maximum_.set (CUBIT_MAX (min.x(), max.x()),
                CUBIT_MAX (min.y(), max.y()),
                CUBIT_MAX (min.z(), max.z()));
}

CubitBox::CubitBox(const CubitVector &min_max) :
    minimum_(min_max),
    maximum_(min_max)
{}

CubitBox::CubitBox(const CubitBox& copy_from):
    minimum_(copy_from.minimum_),
    maximum_(copy_from.maximum_)
{}

CubitBox::CubitBox(const std::vector<CubitVector> & pts)
{
  minimum_.set(0.0, 0.0, 0.0);
  maximum_.set(0.0, 0.0, 0.0);
  for ( size_t i=0; i<pts.size(); i++ )
  {
    double x = pts[i].x();
    double y = pts[i].y();
    double z = pts[i].z();
    minimum_.set(CUBIT_MIN(minimum_.x(), x),
                 CUBIT_MIN(minimum_.y(), y),
                 CUBIT_MIN(minimum_.z(), z));
    maximum_.set(CUBIT_MAX(maximum_.x(), x),
                 CUBIT_MAX(maximum_.y(), y),
                 CUBIT_MAX(maximum_.z(), z));
  }
}

void CubitBox::reset(const CubitVector& vector)
{
  minimum_ = vector;
  maximum_ = vector;
}

void CubitBox::reset(const CubitVector &min, const CubitVector &max)
{
  minimum_.set (CUBIT_MIN (min.x(), max.x()),
		CUBIT_MIN (min.y(), max.y()),
		CUBIT_MIN (min.z(), max.z()));
  maximum_.set (CUBIT_MAX (min.x(), max.x()),
		CUBIT_MAX (min.y(), max.y()),
		CUBIT_MAX (min.z(), max.z()));
}

void CubitBox::reset(const CubitBox &box)
{
  minimum_ = box.minimum_;
  maximum_ = box.maximum_;
}

CubitVector CubitBox::minimum() const
{
  return CubitVector(minimum_);
}

CubitVector CubitBox::maximum() const
{
  return CubitVector(maximum_);
}

double CubitBox::max_x() const
{
  return maximum_.x();
}

double CubitBox::max_y() const
{
  return maximum_.y();
}

double CubitBox::max_z() const
{
  return maximum_.z();
}

double CubitBox::min_x() const
{
  return minimum_.x();
}

double CubitBox::min_y() const
{
  return minimum_.y();
}

double CubitBox::min_z() const
{
  return minimum_.z();
}

CubitVector CubitBox::center() const
{
  return CubitVector(minimum_ + maximum_) / 2.0;
}

CubitVector CubitBox::diagonal() const
{
  return CubitVector(maximum_ - minimum_);
}

double CubitBox::minimum_range( void)
{
  return CUBIT_MIN(CUBIT_MIN( x_range(), y_range()), z_range());
}


double CubitBox::maximum_range( void)
{
  return CUBIT_MAX( CUBIT_MAX( x_range(), y_range()), z_range());
}

double CubitBox::x_range() const
{ return (maximum_.x() - minimum_.x()); }

double CubitBox::y_range() const
{ return (maximum_.y() - minimum_.y()); }

double CubitBox::z_range() const
{ return (maximum_.z() - minimum_.z()); }

CubitBox& CubitBox::operator=(const CubitBox &box)
{
  if (this != &box)
    {
      minimum_ = box.minimum_;
      maximum_ = box.maximum_;
    }
  return *this;
}

bool CubitBox::overlap( double tolerance , const CubitBox& box ) const
{
  //     | - note the '!'.  This condition checks that the boxes
  //     |   do NOT intersect and negates the result. 
  return ! ( (box.minimum_.x() - maximum_.x() > tolerance) ||
             (box.minimum_.y() - maximum_.y() > tolerance) ||
             (box.minimum_.z() - maximum_.z() > tolerance) ||
             (minimum_.x() - box.maximum_.x() > tolerance) ||
             (minimum_.y() - box.maximum_.y() > tolerance) ||
             (minimum_.z() - box.maximum_.z() > tolerance) );
}
             
  
  

// Union of this and box
CubitBox& CubitBox::operator|=(const CubitBox& box)
{
  minimum_.x(CUBIT_MIN(minimum_.x(), box.minimum_.x()));
  minimum_.y(CUBIT_MIN(minimum_.y(), box.minimum_.y()));
  minimum_.z(CUBIT_MIN(minimum_.z(), box.minimum_.z()));

  maximum_.x(CUBIT_MAX(maximum_.x(), box.maximum_.x()));
  maximum_.y(CUBIT_MAX(maximum_.y(), box.maximum_.y()));
  maximum_.z(CUBIT_MAX(maximum_.z(), box.maximum_.z()));
  return *this;
}

// Union of this and vector
CubitBox& CubitBox::operator|=(const CubitVector& vector)
{
  minimum_.x(CUBIT_MIN(minimum_.x(), vector.x()));
  minimum_.y(CUBIT_MIN(minimum_.y(), vector.y()));
  minimum_.z(CUBIT_MIN(minimum_.z(), vector.z()));

  maximum_.x(CUBIT_MAX(maximum_.x(), vector.x()));
  maximum_.y(CUBIT_MAX(maximum_.y(), vector.y()));
  maximum_.z(CUBIT_MAX(maximum_.z(), vector.z()));
  return *this;
}

// Intersection of this and box
CubitBox& CubitBox::operator&=(const CubitBox& box)
{
  minimum_.x(CUBIT_MAX(minimum_.x(), box.minimum_.x()));
  minimum_.y(CUBIT_MAX(minimum_.y(), box.minimum_.y()));
  minimum_.z(CUBIT_MAX(minimum_.z(), box.minimum_.z()));

  maximum_.x(CUBIT_MIN(maximum_.x(), box.maximum_.x()));
  maximum_.y(CUBIT_MIN(maximum_.y(), box.maximum_.y()));
  maximum_.z(CUBIT_MIN(maximum_.z(), box.maximum_.z()));

  if (minimum_.x() > maximum_.x() ||
      minimum_.y() > maximum_.y() ||
      minimum_.z() > maximum_.z())
  {
    minimum_.set(0.0, 0.0, 0.0);
    maximum_.set(0.0, 0.0, 0.0);
  }
  return *this;
}

CubitBox& CubitBox::operator*=(double scale)
{
  CubitVector center_vec = center();
  *this -= center_vec;
  minimum_ *= scale;
  maximum_ *= scale;
  *this += center_vec;
  return *this;
}

CubitBox& CubitBox::operator/=(double scale)
{
  if(scale == 0.0)
    throw ("Cannot Divide by Zero");
  //assert(scale != 0.0);
  *this *= 1/scale;
  return *this;
}

CubitBox& CubitBox::operator+=(const CubitVector& offset)
{
  minimum_ += offset;
  maximum_ += offset;
  return *this;
}

CubitBox& CubitBox::operator-=(const CubitVector& offset)
{
  minimum_ -= offset;
  maximum_ -= offset;
  return *this;
}

int CubitBox::operator<(const CubitBox& box) const
{
  return (box.minimum_.x() < minimum_.x() &&
          box.minimum_.y() < minimum_.y() &&
          box.minimum_.z() < minimum_.z() &&
          box.maximum_.x() > maximum_.x() &&
          box.maximum_.y() > maximum_.y() &&
          box.maximum_.z() > maximum_.z());
}

int CubitBox::operator<=(const CubitBox& box) const
{
  return (box.minimum_.x() <= minimum_.x() &&
          box.minimum_.y() <= minimum_.y() &&
          box.minimum_.z() <= minimum_.z() &&
          box.maximum_.x() >= maximum_.x() &&
          box.maximum_.y() >= maximum_.y() &&
          box.maximum_.z() >= maximum_.z());
}

int CubitBox::operator>(const CubitBox& box) const
{
  return (box.minimum_.x() > minimum_.x() &&
          box.minimum_.y() > minimum_.y() &&
          box.minimum_.z() > minimum_.z() &&
          box.maximum_.x() < maximum_.x() &&
          box.maximum_.y() < maximum_.y() &&
          box.maximum_.z() < maximum_.z());
}

int CubitBox::operator>=(const CubitBox& box) const
{
  return (box.minimum_.x() >= minimum_.x() &&
          box.minimum_.y() >= minimum_.y() &&
          box.minimum_.z() >= minimum_.z() &&
          box.maximum_.x() <= maximum_.x() &&
          box.maximum_.y() <= maximum_.y() &&
          box.maximum_.z() <= maximum_.z());
}

int CubitBox::operator>(const CubitVector& vect) const
{
  return (vect.x() > minimum_.x() && vect.x() < maximum_.x() &&
          vect.y() > minimum_.y() && vect.y() < maximum_.y() &&
          vect.z() > minimum_.z() && vect.z() < maximum_.z() );
}

int CubitBox::operator>=(const CubitVector& vect) const
{
  return (vect.x() >= minimum_.x() && 
          vect.x() <= maximum_.x() &&
          vect.y() >= minimum_.y() && 
          vect.y() <= maximum_.y() &&
          vect.z() >= minimum_.z() && 
          vect.z() <= maximum_.z() );
}

CubitBox operator|(const CubitBox& lhs, const CubitBox& rhs)
{
  return CubitBox(lhs) |= rhs;
}

CubitBox operator&(const CubitBox& lhs, const CubitBox& rhs)
{
  return CubitBox(lhs) &= rhs;
}

CubitBox operator*(double lhs, const CubitBox& rhs)
{
  return CubitBox(rhs) *= lhs;
}

CubitBox operator*(const CubitBox& lhs, double rhs)
{
  return CubitBox(lhs) *= rhs;
}

CubitBox operator/(const CubitBox& lhs, double rhs)
{
  return CubitBox(lhs) /= rhs;
}

CubitBox operator+(const CubitBox& lhs, const CubitVector& rhs)
{
  return CubitBox(lhs) += rhs;
}

CubitBox operator-(const CubitBox& lhs, const CubitVector& rhs)
{
  return CubitBox(lhs) -= rhs;
}

double CubitBox::distance_squared( const CubitVector& point ) const
{
  return (point - closest_point(point)).length_squared();
}



//-------------------------------------------------------------------------
// Purpose       : Find the closest point on the CubitBox
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/98
//-------------------------------------------------------------------------
CubitVector CubitBox::closest_point( const CubitVector& point ) const
{
  CubitVector result;
  
  if( point.x() < minimum_.x() )
    result.x( minimum_.x() );
  else if( point.x() > maximum_.x() )
    result.x( maximum_.x() );
  else
    result.x( point.x() );
  
  if( point.y() < minimum_.y() )
    result.y( minimum_.y() );
  else if( point.y() > maximum_.y() )
    result.y( maximum_.y() );
  else
    result.y( point.y() );
  
  if( point.z() < minimum_.z() )
    result.z( minimum_.z() );
  else if( point.z() > maximum_.z() )
    result.z( maximum_.z() );
  else
    result.z( point.z() );
  
  return result;
}

CubitBox::CubitBox(const double min[3], const double max[3])
{
  minimum_.set (CUBIT_MIN (min[0], max[0]),
                CUBIT_MIN (min[1], max[1]),
                CUBIT_MIN (min[2], max[2]));
  maximum_.set (CUBIT_MAX (min[0], max[0]),
                CUBIT_MAX (min[1], max[1]),
                CUBIT_MAX (min[2], max[2]));
}

void CubitBox::reset(const double min[3], const double max[3])
{
  minimum_.set (CUBIT_MIN (min[0], max[0]),
                CUBIT_MIN (min[1], max[1]),
                CUBIT_MIN (min[2], max[2]));
  maximum_.set (CUBIT_MAX (min[0], max[0]),
                CUBIT_MAX (min[1], max[1]),
                CUBIT_MAX (min[2], max[2]));
}

void CubitBox::get_corners( CubitVector vectors[8] ) const
{
  vectors[0] = minimum_;
  vectors[1] = CubitVector (maximum_.x(), minimum_.y(), minimum_.z());
  vectors[2] = CubitVector (maximum_.x(), maximum_.y(), minimum_.z());
  vectors[3] = CubitVector (minimum_.x(), maximum_.y(), minimum_.z());
  
  vectors[4] = CubitVector (minimum_.x(), minimum_.y(), maximum_.z());
  vectors[5] = CubitVector (maximum_.x(), minimum_.y(), maximum_.z());
  vectors[6] = maximum_;
  vectors[7] = CubitVector (minimum_.x(), maximum_.y(), maximum_.z());
}

int CubitBox::operator&&(const CubitBox& box) const 
{
    // Return false if there is no overlap
    // along at least one of the 3 axes.
  if (minimum_.x() > box.maximum_.x() ||
      maximum_.x() < box.minimum_.x() ||
      minimum_.y() > box.maximum_.y() ||
      maximum_.y() < box.minimum_.y() ||
      minimum_.z() > box.maximum_.z() ||
      maximum_.z() < box.minimum_.z() )
    return CUBIT_FALSE;
  
    // If you didn't return false...
  return CUBIT_TRUE;
}

int CubitBox::operator||(const CubitBox& box) const 
{
    // Return false if there is no overlap
    // along at least one of the 3 axes.
  if (minimum_.x() >= box.maximum_.x() ||
      maximum_.x() <= box.minimum_.x() ||
      minimum_.y() >= box.maximum_.y() ||
      maximum_.y() <= box.minimum_.y() ||
      minimum_.z() >= box.maximum_.z() ||
      maximum_.z() <= box.minimum_.z() )
    return CUBIT_FALSE;
  
    // If you didn't return false...
  return CUBIT_TRUE;
}

int CubitBox::operator<=(const CubitVector& vect) const
{
  return (vect.x() <= minimum_.x() ||
          vect.x() >= maximum_.x() ||
          vect.y() <= minimum_.y() ||
          vect.y() >= maximum_.y() ||
          vect.z() <= minimum_.z() ||
          vect.z() >= maximum_.z() );
}

/*
Ray-axis Aligned box intersection by
KAY, T. L., AND KAJIYA, J. T.
In Computer Graphics (SIGGRAPH '86 Proceedings) (Aug. 1986), D. C. Evans and R. J. Athay, Eds., vol. 20, pp. 269-278. 

Determines if there is an intersection between a ray and box but doesn't calculate the intersection.

*/
bool CubitBox::intersect(const CubitVector& ray_origin, const CubitVector& ray_direction)
{
    double tNear=CUBIT_DBL_MIN;
    double tFar=CUBIT_DBL_MAX;

    for(int i=0;i<3;i++)
    {
        //X plane
        if(fabs(ray_direction[i])<=CUBIT_DBL_MIN)
        {
            //Ray parallel to x plane
            if(ray_origin[i]<this->minimum_[i] || ray_origin[i]>this->maximum_[i])
            {
                //Not between
                return false;
            }
            else
            {
                return true;
            }
        }

        double t1 = (this->minimum_[i] - ray_origin[i]) / ray_direction[i];
        double t2 = (this->maximum_[i] - ray_origin[i]) / ray_direction[i];

        if(t1>t2)
        {
            double temp=t1;
            t1=t2;
            t2=temp;
        }
        if(t1>tNear)
        {
            tNear=t1;
        }

        if(t2<tFar)
        {
            tFar=t2;
        }

        if(tNear>tFar)
        {
            return false;
        }

        if(tFar<0)
        {
            return false;
        }

    }

    return true;
}

/* 
Fast Ray-Box Intersection
by Andrew Woo
from "Graphics Gems", Academic Press, 1990
*/

#define NUMDIM	3
#define RIGHT	0
#define LEFT	1
#define MIDDLE	2

//char HitBoundingBox(minB,maxB, origin, dir,coord)
//double minB[NUMDIM], maxB[NUMDIM];		/*box */
//double origin[NUMDIM], dir[NUMDIM];		/*ray */
//double coord[NUMDIM];				/* hit point */
bool CubitBox::intersect
(const CubitVector& ray_origin,
  const CubitVector& ray_direction,
  CubitVector& intersection_pt)
{
  bool inside = true;
  int quadrant[NUMDIM];
  register int i;
  int whichPlane;
  double maxT[NUMDIM];
  double candidatePlane[NUMDIM];

  /* Find candidate planes; this loop can be avoided if
  rays cast all from the eye(assume perpsective view) */
  for (i=0; i<NUMDIM; i++)
  {
    if((ray_origin)[i] < (this->minimum())[i])
    {
      quadrant[i] = LEFT;
      candidatePlane[i] = (this->minimum())[i];
      inside = false;
    }
    else if ( (ray_origin)[i] > (this->maximum())[i])
    {
      quadrant[i] = RIGHT;
      candidatePlane[i] = (this->maximum())[i];
      inside = false;
    }
    else	
    {
      quadrant[i] = MIDDLE;
    }
  }

  /* Ray origin inside bounding box */
  if(inside)
  {
    intersection_pt = ray_origin;
    return true;
  }

  /* Calculate T distances to candidate planes */
  for (i = 0; i < NUMDIM; i++)
  {
    if (quadrant[i] != MIDDLE && (ray_direction)[i] !=0.)
      maxT[i] = ( candidatePlane[i]-(ray_origin)[i]) / ((ray_direction)[i]);
    else
      maxT[i] = -1.;
  }

  /* Get largest of the maxT's for final choice of intersection */
  whichPlane = 0;
  for (i = 1; i < NUMDIM; i++)
  {
    if (maxT[whichPlane] < maxT[i])
      whichPlane = i;
  }

  /* Check final candidate actually inside box */
  double intersection_coords[3];

  if (maxT[whichPlane] < 0.)
    return false;

  for (i = 0; i < NUMDIM; i++)
  {
    if (whichPlane != i)
    {
      intersection_coords[i] = (ray_origin)[i] + maxT[whichPlane] * (ray_direction)[i];
      if (intersection_coords[i] < (this->minimum())[i] ||
          intersection_coords[i] > (this->maximum())[i])
        return false;
    }
    else
    {
      intersection_coords[i] = candidatePlane[i];
    }
  }

  return true;  /* ray hits box */
}
