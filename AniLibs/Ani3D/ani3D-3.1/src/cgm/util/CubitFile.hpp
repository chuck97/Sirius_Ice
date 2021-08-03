//- Class:          CubitFile
//- Description:    Class with functions to get work with files

#ifndef CUBITFILE_HPP
#define CUBITFILE_HPP

#include <stdio.h>
#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

// wrapper for C FILE*
// Includes support for non-ascii filenames.
// Automatically opens/closes the file.
class CUBIT_UTIL_EXPORT CubitFile
{
public:
  // default constructor
  CubitFile();

  // constructor that also opens the file
  CubitFile(const CubitString& file, const char* mode);

  // closes the file, if open
  virtual ~CubitFile();

  // opens a file
  bool open(const CubitString& file, const char* mode);

  // closes the file
  void close();

  // Return the FILE* handle.
  // Returns NULL if status is not 'OK'.
  FILE* file() const;

  // returns whether FILE* is valid
  operator bool () const;

  // Return status of opening the file.
  // For example ENOENT for non-existant file.
  int error();

protected:
  FILE* mFile;
  int mError;
};


#endif

