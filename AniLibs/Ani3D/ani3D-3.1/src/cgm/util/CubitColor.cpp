
#include "CubitColor.hpp"


CubitColor::CubitColor()
  : red(0), green(0), blue(0), alpha(0)
{
}

CubitColor::CubitColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
  : red(r), green(g), blue(b), alpha(a)
{
}

CubitColor::CubitColor(const CubitColor& other)
  : red(other.red), green(other.green), blue(other.blue), alpha(other.alpha)
{
}

bool CubitColor::operator==(const CubitColor& other) const
{
  return red == other.red && green == other.green && blue == other.blue && alpha == other.alpha;
}

bool CubitColor::operator!=(const CubitColor& other) const
{
  return !operator==(other);
}

//! create a CubitColor from double r,g,b. double value is between 0.0 and 1.0.
CubitColor CubitColor::rgb_from_double(double r, double g, double b, double a)
{
  return CubitColor(
      static_cast<unsigned char>(r*255.0),
      static_cast<unsigned char>(g*255.0),
      static_cast<unsigned char>(b*255.0),
      static_cast<unsigned char>(a*255.0)
      );
}

