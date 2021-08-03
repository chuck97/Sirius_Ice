//- Class:          CubitFile
//- Description:    Class with functions to work with files

#include "CubitFile.hpp"
#include <errno.h>

CubitFile::CubitFile()
  : mFile(NULL), mError(0)
{
}

CubitFile::CubitFile(const CubitString& file, const char* mode)
  : mFile(NULL)
{
  open(file, mode);
}

CubitFile::~CubitFile()
{
  close();
}

bool CubitFile::open(const CubitString& file, const char* mode)
{
  close();

  this->mError = 0;

#ifdef _WIN32
  this->mFile = _wfopen(CubitString::toUtf16(file).c_str(), CubitString::toUtf16(mode).c_str());
#else
  this->mFile = fopen(file.c_str(), mode);
#endif
  if(!this->mFile)
  {
    this->mError = errno;
  }

  return mError == 0;
}

void CubitFile::close()
{
  if(mFile)
  {
    fclose(mFile);
  }
  mFile = NULL;
  mError = 0;
}

FILE* CubitFile::file() const
{
  return this->mFile;
}

CubitFile::operator bool () const
{
  return this->mFile != NULL;
}

int CubitFile::error()
{
  return this->mError;
}

