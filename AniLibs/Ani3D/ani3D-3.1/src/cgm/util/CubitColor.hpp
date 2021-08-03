#ifndef CUBIT_COLOR_HPP
#define CUBIT_COLOR_HPP

#include "CubitUtilConfigure.h"

//! CubitColor class to represent RGBA colors.
//! Values are unsigned chars with a range of 0 to 255.
struct CUBIT_UTIL_EXPORT CubitColor
{
  CubitColor();
  CubitColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a=255);
  CubitColor(const CubitColor& other);

  bool operator==(const CubitColor& other) const;

  bool operator!=(const CubitColor& other) const;

  //! create a CubitColor from double r,g,b. double value is between 0.0 and 1.0.
  static CubitColor rgb_from_double(double r, double g, double b, double a=1.0);

  unsigned char red;
  unsigned char green;
  unsigned char blue;
  unsigned char alpha;
};


const CubitColor CUBIT_DEFAULT_COLOR = CubitColor();
const CubitColor CUBIT_BLACK         = CubitColor(0,0,0);
const CubitColor CUBIT_GREY          = CubitColor(127,127,127);
const CubitColor CUBIT_ORANGE        = CubitColor(255,165,0);
const CubitColor CUBIT_RED           = CubitColor(255,0,0);
const CubitColor CUBIT_GREEN         = CubitColor(0,255,0);
const CubitColor CUBIT_YELLOW        = CubitColor(255,255,0);
const CubitColor CUBIT_MAGENTA       = CubitColor(255,0,255);
const CubitColor CUBIT_CYAN          = CubitColor(0,255,255);
const CubitColor CUBIT_BLUE          = CubitColor(0,0,255);
const CubitColor CUBIT_WHITE         = CubitColor(255,255,255);
const CubitColor CUBIT_BROWN         = CubitColor(165,42,42);
const CubitColor CUBIT_GOLD          = CubitColor(255,215,0);
const CubitColor CUBIT_LIGHTBLUE     = CubitColor(173,216,230);
const CubitColor CUBIT_LIGHTGREEN    = CubitColor(0,204,0);
const CubitColor CUBIT_SALMON        = CubitColor(250,128,114);
const CubitColor CUBIT_CORAL         = CubitColor(255,127,80);
const CubitColor CUBIT_PINK          = CubitColor(255,192,203);

const int CUBIT_DEFAULT_COLOR_INDEX = -1;
const int CUBIT_BLACK_INDEX         =  0;
const int CUBIT_GREY_INDEX          =  1;
const int CUBIT_ORANGE_INDEX        =  2;
const int CUBIT_RED_INDEX           =  3;
const int CUBIT_GREEN_INDEX         =  4;
const int CUBIT_YELLOW_INDEX        =  5;
const int CUBIT_MAGENTA_INDEX       =  6;
const int CUBIT_CYAN_INDEX          =  7;
const int CUBIT_BLUE_INDEX          =  8;
const int CUBIT_WHITE_INDEX         =  9;
const int CUBIT_BROWN_INDEX         = 10;
const int CUBIT_GOLD_INDEX          = 11;
const int CUBIT_LIGHTBLUE_INDEX     = 12;
const int CUBIT_LIGHTGREEN_INDEX    = 13;
const int CUBIT_SALMON_INDEX        = 14;
const int CUBIT_CORAL_INDEX         = 15;
const int CUBIT_PINK_INDEX          = 16;

#endif

