//- Class:          CubitDirIterator
//- Description:    Class with functions to get files from a directory, etc.

#include "CubitDirIterator.hpp"
#include "CubitString.hpp"

#ifdef _WIN32
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <windows.h>
#else
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <dirent.h>
  #include <cstdlib>
  #include <sys/param.h>
  #include <unistd.h>
  #include <pwd.h>
  #include <fnmatch.h>
#endif
#include <errno.h>
#include <string.h>

struct CubitDirIterator::Helper
{
  Helper()
  {
#ifdef _WIN32
    mFileHandle = INVALID_HANDLE_VALUE;
#else
    mDirHandle = NULL;
    mFileHandle = NULL;
#endif
  }

#ifdef _WIN32
  WIN32_FIND_DATAW mDirHandle;
  HANDLE mFileHandle;
#else
  DIR* mDirHandle;
  dirent* mFileHandle;
#endif
};
  

CubitDirIterator::CubitDirIterator(const CubitString& path, const CubitString& pattern)
{
  this->mHelper = new Helper;
  open(path, pattern);
}

CubitDirIterator::~CubitDirIterator()
{
  cleanup();
  delete this->mHelper;
}

void CubitDirIterator::open(const CubitString& path, const CubitString& pattern)
{
  cleanup();

  if(pattern.length())
    mPattern = pattern;
  else
    mPattern = "";

  this->atEnd = false;
#ifdef _WIN32
  CubitString p = path;
  if (p.get_at(p.length()-1) != '\\')
    p += "\\";
  p += mPattern.length() == 0 ? "*" : mPattern;
  this->mHelper->mFileHandle = FindFirstFileW(CubitString::toUtf16(p).c_str(), &this->mHelper->mDirHandle);
  if(this->mHelper->mFileHandle == INVALID_HANDLE_VALUE)
    this->atEnd = true;
#else
  this->mHelper->mDirHandle = opendir(path.c_str());
  if(this->mHelper->mDirHandle)
  {
    do
    {
      this->mHelper->mFileHandle = readdir(this->mHelper->mDirHandle);
    } while (this->mHelper->mFileHandle && 
             (mPattern.length() != 0 && 
             fnmatch(mPattern.c_str(), this->mHelper->mFileHandle->d_name, 0))
             );
  }
  if(!this->mHelper->mDirHandle || !this->mHelper->mFileHandle)
    this->atEnd = true;
#endif
}

void CubitDirIterator::cleanup()
{
#ifdef _WIN32
  if(this->mHelper->mFileHandle != 0)
    FindClose(this->mHelper->mFileHandle);
  this->mHelper->mFileHandle = 0;
#else
  if(this->mHelper->mDirHandle)
    closedir(this->mHelper->mDirHandle);
  this->mHelper->mFileHandle = 0;
  this->mHelper->mDirHandle = 0;
#endif
}

bool CubitDirIterator::has_next()
{
  return !this->atEnd;
}

CubitString CubitDirIterator::next()
{
  CubitString file;
  if(this->atEnd)
    return file;

#ifdef _WIN32
  file = CubitString::toUtf8(this->mHelper->mDirHandle.cFileName);
  BOOL result = FindNextFileW(this->mHelper->mFileHandle, &this->mHelper->mDirHandle);
  if(result == 0)
    this->atEnd = true;
#else
  file = this->mHelper->mFileHandle->d_name;
  do
  {
    this->mHelper->mFileHandle = readdir(this->mHelper->mDirHandle);
  } while (this->mHelper->mFileHandle && 
           (mPattern.length() != 0 && 
           fnmatch(mPattern.c_str(), this->mHelper->mFileHandle->d_name, 0))
           );
  if(this->mHelper->mFileHandle == 0)
    this->atEnd = true;
#endif
  return file;
}

