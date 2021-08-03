//- Class: CubitPlane
//- Description: This file defines the CubitPlane class.
//- Owner: Greg Sjaardema
//- Checked by:


#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "CubitPlane.hpp"
#include "CubitMatrix.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"

CubitPlane::CubitPlane(const double A, const double B,
                       const double C, const double D) : normal_(A, B, C)
{
  double length = normal_.length();
  d_ = D / length;
  normal_ /= length;
}

CubitPlane::CubitPlane() : normal_(0.0, 0.0, 0.0)
{
  d_ = 0.0;
}

CubitPlane::CubitPlane(const CubitVector &Normal, const double D) :
  normal_(Normal)
{
  // make sure that Normal is not (0,0,0) causing a division by 0
  if(normal_.length()==0)
    throw std::invalid_argument ("Normal Length must not be 0");
  //assert(!normal_.length()==0);

  double length = normal_.length();
  d_ = D / length;
  normal_ /= length;
}
	
CubitPlane::CubitPlane(const CubitVector &Normal, const CubitVector &point) :
  normal_(Normal)
{
// plane_D = -(plane_normal.x()*X + plane_normal.y()*Y + plane_normal.z()*Z)
  normal_.normalize();
  d_ = -1.0 * (normal_.x()*point.x() + normal_.y()*point.y() +
	       normal_.z()*point.z());
}

void CubitPlane::set(const CubitVector &Normal, const CubitVector &point)
{
  normal_ = Normal;
  normal_.normalize();
  d_ = -1.0 * (normal_.x()*point.x() + normal_.y()*point.y() +
	       normal_.z()*point.z());
}

CubitPlane::CubitPlane(DLIList<CubitVector> &positions)
{
// Newells method for determining approximate plane equation.
// Plane equation is of the form:
// 0 = plane_D +plane_normal.x()*X + plane_normal.y()*Y + plane_normal.z()*Z
  if(positions.size() < 3)
    throw std::invalid_argument("Must have 3 or more Points");
  CubitVector vector_diff;
  CubitVector vector_sum;
  CubitVector ref_point = CubitVector (0.0, 0.0, 0.0);
  CubitVector tmp_vector;
  normal_ = CubitVector(0.0, 0.0, 0.0);
  CubitVector vector1, vector2;
  
  for (int i=positions.size(); i > 0; i--)
  {
    vector1 = positions.get_and_step();
    vector2 = positions.get();
    vector_diff = (vector2) - (vector1);
    ref_point += (vector1);
    vector_sum = (vector1) + (vector2);
    
    tmp_vector.set(vector_diff.y() * vector_sum.z(),
                   vector_diff.z() * vector_sum.x(),
                   vector_diff.x() * vector_sum.y());
    normal_ += tmp_vector;
  }
  double divisor = positions.size() * normal_.length();
//  if (divisor < .000000001)
//  {
//    PRINT_ERROR("Colinear basis vector in CubitPlane constructor");
//  }
  d_ = (ref_point % normal_) / divisor;
  normal_.normalize();
  normal_ *= -1.0;
}

bool CubitPlane::fit_points(const std::vector<CubitVector> &positions)
{
  //- use the method of linear regression
  //
  // by creating an Ax=b system where:
  //
  //            ___                                 ___
  //            |   x[i]*x[i],    x[i]*y[i],    x[i]  |
  //  A = sum_i |   x[i]*y[i],    y[i]*y[i],    y[i]  |
  //            |   x[i],         y[i],         n     |
  //            ---                                 ---
  //            ___         ___
  //            |  x[i]*z[i]  |
  //  b = sum_i |  y[i]*z[i]  |
  //            |  z[i]}      |
  //            ---         ---
  // But this assumes that z is a linear function of x and y (i.e. z component
  // of plane normal != 0.).  We check the rank of A to catch this case.  We
  // handle this case, by transposing the x or y coordinates for the z,
  // computing the normal, then transposing back (see itry below).

  for ( int itry=0; itry < 3; itry++ )
  {
    CubitMatrix A(3,3);
    std::vector<double> b(3);
    A.set(2,2,positions.size());
    for ( size_t ipt=0; ipt<positions.size(); ipt++ )
    {
      double x=0., y=0., z=0.;
      if ( itry == 0 )
      {
        x = positions[ipt].x();
        y = positions[ipt].y();
        z = positions[ipt].z();
      }
      else if ( itry == 1 )
      {
        x = positions[ipt].x();
        z = positions[ipt].y();
        y = positions[ipt].z();
      }
      else if ( itry == 2 )
      {
        z = positions[ipt].x();
        y = positions[ipt].y();
        x = positions[ipt].z();
      }
      double sum_xx = x*x + A.get(0,0);
      double sum_yy = y*y + A.get(1,1);
      double sum_xy = x*y + A.get(1,0);
      double sum_x  = x   + A.get(2,0);
      double sum_y  = y   + A.get(2,1);
      A.set(0,0,sum_xx);
      A.set(1,1,sum_yy);
      A.set(1,0,sum_xy); A.set(0,1,sum_xy);
      A.set(2,0,sum_x);  A.set(0,2,sum_x);
      A.set(1,2,sum_y);  A.set(2,1,sum_y);
      b[0] += x*z;
      b[1] += y*z;
      b[2] += z;
    }

    // Make sure we have full rank.  If not, it means:
    // if itry==0, plane is parallel to z-axis.
    // if itry==1, plane is parallel to z-axis & y-axis
    // if itry==2, points are colinear
    int rank = A.rank();
    if ( rank == 3 )
    {
      std::vector<double> coef(3);
      A.solveNxN( b, coef );
      if ( itry == 0 )
        normal_.set(coef[0], coef[1], -1);
      else if ( itry == 1 )
        normal_.set(coef[0], -1, coef[1]);
      else if ( itry == 2 )
        normal_.set(-1, coef[1], coef[0]);
      normal_.normalize();
      d_ = 0;
      for ( size_t ipt=0; ipt<positions.size(); ipt++ )
      {
        d_ += normal_.x()*positions[ipt].x() +
              normal_.y()*positions[ipt].y() +
              normal_.z()*positions[ipt].z();
      }
      d_ /= positions.size();
      return true;
    }
  }
  d_ = 0.0;
  normal_.set(0.,0.,0.);
  return false;
}

CubitPlane::CubitPlane(const CubitPlane& copy_from)
{
  normal_ = copy_from.normal_;
  d_ = copy_from.d_;
}

int CubitPlane::mk_plane_with_points( const CubitVector& V0, 
                                      const CubitVector& V1,
                                      const CubitVector& V2 )
{
    // First check to make sure that the points are not collinear
  
    // Vector going from Vertex 0 to Vertex 1	    
  CubitVector vector01 = V1 - V0;
  vector01.normalize();
  
    // Vector going from Vertex 0 to Vertex 2	  
  CubitVector vector02 = V2 - V0;
  vector02.normalize();
  
    // If the 3 points are collinear, then the cross product of these two
    // vectors will yield a null vector (one whose length is zero).
  normal_ = vector01 * vector02;
  
  if(normal_.length() <= CUBIT_RESABS)
  {
    PRINT_ERROR("Points are collinear.\n"
                "       Cannot create a CubitPlane object.\n");
    d_ = 0.0;
    return CUBIT_FAILURE;
  }
  
  normal_.normalize();
  
  double D0 =   -(normal_ % V0);
  double D1 =   -(normal_ % V1);
  double D2 =   -(normal_ % V2);
  d_ = (D0 + D1 + D2) / 3.0;
  
  return CUBIT_SUCCESS;
}

CubitVector CubitPlane::point_on_plane() const
{
  if ( fabs( normal_.x() ) > CUBIT_RESABS )
    return CubitVector(-d_ / normal_.x(), 0, 0);
  else if ( fabs( normal_.y() ) > CUBIT_RESABS )
    return CubitVector(0, -d_ / normal_.y(), 0);
  else if ( fabs( normal_.z() ) > CUBIT_RESABS )
    return CubitVector(0, 0, -d_ / normal_.z());
    // If A B and C are all zero, the plane is invalid,
    // Just return <0,0,0>
  return CubitVector(0,0,0);
}

//-  Plane assignment operator.
CubitPlane& CubitPlane::operator=(const CubitPlane &plane)
{
  if (this != &plane)
  {
    normal_ = plane.normal_;
    d_ = plane.d_;
  }
  return *this;
}

double CubitPlane::distance(const CubitVector &vector) const
{
  return normal_ % vector + d_;
}

CubitVector CubitPlane::intersect(const CubitVector &base,
                                  const CubitVector &direction) const
{
  double t = -(normal_ % base + d_) / (normal_ % direction);
  return (base + direction * t);
}

int CubitPlane::intersect(const CubitPlane &plane_2,
                          CubitVector &origin, CubitVector &vector) const
{
    // Code from Graphics Gems III, Georgiades
    // Calculate the line of intersection between two planes. 
    // Initialize the unit direction vector of the line of intersection in
    // xdir.
    // Pick the point on the line of intersection on the coordinate plane most
    // normal to xdir.
    // Return TRUE if successful, FALSE otherwise (indicating that the planes
    // don't intersect). The order in which the planes are given determines the
    // choice of direction of xdir.
    //
    // int GetXLine(vect4 *pl1, vect4 *plane_2, vect3 *vector, vect3 *xpt)
  double invdet;  // inverse of 2x2 matrix determinant
  vector = normal() * plane_2.normal();
  CubitVector dir2(vector.x()*vector.x(),
                   vector.y()*vector.y(),
                   vector.z()*vector.z());
  CubitVector plane1n = normal();
  CubitVector plane2n = plane_2.normal();
  
  if (dir2.z() > dir2.y() && dir2.z() > dir2.x() && dir2.z() > CUBIT_RESABS)
  {
      // then get a point on the XY plane
    invdet = 1.0 / vector.z();
    
      //solve < pl1.x * origin.x + pl1.y * origin.y = - pl1.w >
      //      < plane2n.x * origin.x + plane2n.y * origin.y = - plane2n.w >
    origin.set(plane1n.y() * plane_2.coefficient() -
               plane2n.y() * coefficient(),
               plane2n.x() * coefficient() -
               plane1n.x() * plane_2.coefficient(),
               0.0);
  }
  else if (dir2.y() > dir2.x() && dir2.y() > CUBIT_RESABS)
  {
      // then get a point on the XZ plane
    invdet = -1.0 / vector.y();
    
      // solve < pl1.x * origin.x + pl1.z * origin.z = -pl1.w >
      //       < plane2n.x * origin.x + plane2n.z * origin.z = -plane2n.w >
    origin.set(plane1n.z() * plane_2.coefficient() -
               plane2n.z() * coefficient(),
               0.0,
               plane2n.x() * coefficient() -
               plane1n.x() * plane_2.coefficient());
  }
  else if (dir2.x() > CUBIT_RESABS)
  {
      // then get a point on the YZ plane
    invdet = 1.0 / vector.x();
    
      // solve < pl1.y * origin.y + pl1.z * origin.z = - pl1.w >
      //             < plane2n.y * origin.y + plane2n.z * origin.z = - plane2n.w >
    origin.set(0.0,
               plane1n.z() * plane_2.coefficient() -
               plane2n.z() * coefficient(),
               plane2n.y() * coefficient() -
               plane1n.y() * plane_2.coefficient());
  }
  else // vector is zero, then no point of intersection exists
    return CUBIT_FALSE;
  
  origin *= invdet;
  invdet = 1.0 / (float)sqrt(dir2.x() + dir2.y() + dir2.z());
  vector *= invdet;
  return CUBIT_TRUE;
}

CubitVector CubitPlane::project( const CubitVector& point ) const
{
  return point - distance(point) * normal();
}
