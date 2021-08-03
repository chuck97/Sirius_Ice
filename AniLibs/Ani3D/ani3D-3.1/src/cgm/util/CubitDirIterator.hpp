//- Class:          CubitDirIterator
//- Description:    Class with functions to iterate over files in directories

#ifndef CubitDirIterator_HPP
#define CubitDirIterator_HPP

#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

// an iterator to iterate over files in a directory
class CUBIT_UTIL_EXPORT CubitDirIterator
{
public:
  CubitDirIterator(const CubitString& path,
                   const CubitString& pattern_match = "");
  virtual ~CubitDirIterator();
  void open(const CubitString& path,
            const CubitString& pattern_match = "");
  bool has_next();
  CubitString next();
protected:
  struct Helper;
  Helper* mHelper;
  CubitString mPattern;
  bool mDirs;
  bool atEnd;
  void cleanup();
};

#endif

