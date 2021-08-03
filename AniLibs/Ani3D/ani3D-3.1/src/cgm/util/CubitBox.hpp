//- Class: CubitBox
//-
//- Description: This file defines the CubitBox class which represents
//- an axis-aligned rectangular box which can be uses as a bounding box.
//-
//- Owner: Greg Sjaardema
//- Checked by: 
//- Version: $Id: 

#ifndef CUBITBOX_HPP
#define CUBITBOX_HPP

#include "CubitVector.hpp"
#include "CubitBoxStruct.h"
#include "CubitUtilConfigure.h"
#include <vector>

class CUBIT_UTIL_EXPORT CubitBox
{
public:
    //- Heading: Constructors and Destructor
  CubitBox(); 
    //- Default constructor.
  
  CubitBox(const CubitVector &min, const CubitVector &max);
    //- Constructor: create box from two CubitVectors
  
  CubitBox(const double min[3], const double max[3] );
    //- Constructor: create box from two coordinates
  
  CubitBox(const CubitVector &min_max);
    //- Constructor: create box from one CubitVector
  
  CubitBox(const CubitBox& copy_from);  //- Copy Constructor
  
  CubitBox(const CubitBoxStruct& from);

  CubitBox(const std::vector<CubitVector> & pts);
  
  ~CubitBox(); 
    //- destructor
    
  CubitBox& bounding_box();
  
  
  void reset(const CubitVector &vector);
  void reset(const CubitVector &min, const CubitVector &max);
  void reset(const CubitBox &box);
  void reset(const double min[3], const double max[3]);
    //- reset ranges

  double max_x() const;
  double max_y() const;
  double max_z() const;

  double min_x() const;
  double min_y() const;
  double min_z() const;

  CubitVector minimum()  const;
  CubitVector maximum()  const;
  CubitVector center()   const;
  CubitVector diagonal() const;
    //- Return Box minimum/maximum/center
  
  void get_corners ( CubitVector corners[8] ) const;
    //- Fills 'corners' with the corners of this box.
    //- The order is:
    //-   0) minimum()
    //-   1-3) Front face (Constant minimum z-plane), normal out of box
    //-        using right hand rule, including corner[0].
    //-   4-7) Same as 0-3, but offset to back face
    //-        (constant maximum z-plane).  Normal of these last 4 points
    //-        is into box relative to back plane (same direction as
    //-        normal w/ first 4 points).  Maximum ends up at index 6.
  
  double x_range() const;
  double y_range() const;
  double z_range() const;
    //- x, y, and z range of the box (max - min)
  
  double minimum_range( void);
  double maximum_range( void);
  //- returns the maimum and maximum range 

  bool overlap( double tolerance, const CubitBox& other_box ) const;
    //- Check if boxes are within passed tolerance of each other.
    //- If tolerance is 0, use && or || operator.
  
  bool intersect(const CubitVector& ray_origin, const CubitVector& ray_direction,
      CubitVector& intersection_pt);
    //- Check if ray intersects box and returns an intersection point


  bool intersect(const CubitVector& ray_origin, const CubitVector& ray_direction);
    //- Check if ray intersects box but doesn't calculate an intersection point.




    //- Heading: Operators

    // Operators that modify {this}
  CubitBox& operator=(const CubitBox &box);
  CubitBox& operator|=(const CubitBox& box);
  CubitBox& operator|=(const CubitVector& vector);
  CubitBox& operator&=(const CubitBox& box);
  CubitBox& operator*=(double scale);
  CubitBox& operator/=(double scale);
  CubitBox& operator+=(const CubitVector& offset);
  CubitBox& operator-=(const CubitVector& offset);
    //- {=}  - Assignment
    //- {|=} - Union of {this} and {box}
    //- {&=} - Intersection (overlap) of {this} and {box}
    //- {*=} - Scale {this} about box center
    //- {/=} - Scale {this} about box center
    //- {+=} - Move {this} by {offset}  CubitVector
    //- {-=} - Move {this} by {-offset} CubitVector
  
    // Operators that check for containment
  int operator< (const CubitBox& box) const;
  int operator<=(const CubitBox& box) const;
  int operator> (const CubitBox& box) const;
  int operator>=(const CubitBox& box) const;
  int operator> (const CubitVector& vect) const;
  int operator>=(const CubitVector& vect) const;
  int operator<=(const CubitVector& vect) const;
  int operator&&(const CubitBox& box) const;
  int operator||(const CubitBox& box) const;
    //- {<}  - Is {this} completely surrounded by {box}?
    //- {>}  - Does {this} completely surround {box}?
    //- {<=,>=} - As above, but inner box may touch
    //-           boundary of outer box.
    //- {>}  - Is {vect} contained within {this}, but not on boundary?
    //- {>=} - Is {vect} contained within or on the boundary of {this}?
    //- {<=} - Is {vect} outside or on boundary of {this}?
    //- {&&} - Do {this} and {box} intersect?  Just butting against each
    //-        other also counts as an intersection. See {||}.
    //- {||} - Do {this} and {box} intersect?  Just butting against each
    //-        other does NOT count as an intersection.  See {&&}.
  
  CubitBox &operator=(const CubitBoxStruct &from);

  operator CubitBoxStruct() 
    {
      CubitBoxStruct to;
      to.minimum_ = minimum_;
      to.maximum_ = maximum_;
      return to;
    }

    // Operators that return a modification of {this}.
    // {this} itself is not modified.
  friend CUBIT_UTIL_EXPORT CubitBox operator|(const CubitBox& lhs, const CubitBox& rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator|(const CubitBox& lhs, const CubitVector& rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator&(const CubitBox& lhs, const CubitBox& rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator*(const CubitBox& lhs, double rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator*(double rhs, const CubitBox& lhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator/(const CubitBox& lhs, double rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator+(const CubitBox& lhs, const CubitVector& rhs);
  friend CUBIT_UTIL_EXPORT CubitBox operator-(const CubitBox& lhs, const CubitVector& rhs);
  
	double distance_squared( const CubitVector& position ) const;
  
  CubitVector closest_point( const CubitVector& position ) const;
    //R CubitVector
    //R- The closest point on the box to the passed position.
    //R- The passed position will be returned if it is within
    //R- the box.
    //I- A position from which to evaluate the closest point
    //I- on the box.
  
private:
  
  CubitVector minimum_; //- X, Y, and Z position of minimum corner
  CubitVector maximum_; //- X, Y, and Z position of maximum corner
};

inline CubitBox& CubitBox::operator=(const CubitBoxStruct &from)  
{
  minimum_ = from.minimum_;
  maximum_ = from.maximum_;
  return *this;
}

inline CubitBox::CubitBox(const CubitBoxStruct &from)  
{
  minimum_ = from.minimum_;
  maximum_ = from.maximum_;
}
#endif
