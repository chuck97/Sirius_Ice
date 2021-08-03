
// We have 4 different implementations here
// Windows, Unix, , HP-UX

#include "CubitDynamicLoader.hpp"

#include <vector>
#include <sys/stat.h>

#include "CubitMessage.hpp"
#include "CubitString.hpp"

static std::vector<CubitString> gSearchPaths;

void CubitDynamicLoader::add_search_path(const char* path)
{
  gSearchPaths.push_back(CubitString(path));
}

static bool absolute_path(const char* path)
{
#ifdef _WIN32
  // either a '\' or a 'x:' drive letter means absolute path
  return path[0] == '\\' || (path[0] != '\0') ? (path[1] == ':') : 0;
#else
  // a '/' is always an absolute path
  return path[0] == '/';
#endif
}

CubitBoolean CubitDynamicLoader::library_exists(const char* library)
{
  struct stat buf;

  // if absolute path, test it directly
  if(absolute_path(library))
  {
    return stat(library, &buf) == 0 ? CUBIT_TRUE : CUBIT_FALSE;
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<(int)gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i] + CubitString("/") + library;
      if(stat(path.c_str(), &buf) == 0)
        return CUBIT_TRUE;
    }
  }
  
  // one more final attempt
  if (stat(library, &buf) == 0)
    return CUBIT_TRUE;

  return CUBIT_FALSE;
}

#ifdef _WIN32

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  // disable message boxes about failure to find dlls
  //UINT old_error_mode = SetErrorMode(SEM_NOOPENFILEERRORBOX | SEM_FAILCRITICALERRORS);
	UINT old_error_mode = SetErrorMode(SEM_FAILCRITICALERRORS);

  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    handle = LoadLibraryW(CubitString::toUtf16(lib).c_str());
  }
  else
  { 
    // try finding with our search paths
    for(unsigned int i=0; i<gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("\\")+lib;
      handle = LoadLibraryW(CubitString::toUtf16(path.c_str()).c_str());
      if(handle != InvalidLibraryHandle)
      {
        SetErrorMode(old_error_mode);
        return handle;
      }
    }
  }
  
  // one more final attempt
  handle = LoadLibraryW(CubitString::toUtf16(lib).c_str());

  SetErrorMode(old_error_mode);

  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  LPWSTR lpMsgBuf = NULL;
  DWORD dw = GetLastError(); 

  FormatMessageW(
    FORMAT_MESSAGE_ALLOCATE_BUFFER | 
    FORMAT_MESSAGE_FROM_SYSTEM |
    FORMAT_MESSAGE_IGNORE_INSERTS,
    NULL,
    dw,
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
    lpMsgBuf,
    0, NULL );

  CubitString error = CubitString::toUtf8(lpMsgBuf);
  LocalFree(lpMsgBuf);

  return error;
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return FreeLibrary(lib) == FALSE ? CUBIT_FAILURE : CUBIT_SUCCESS;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  return reinterpret_cast<void*>(GetProcAddress(lib, symbol));
}

const char* CubitDynamicLoader::library_prefix()
{
  return "";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".dll";
}

#elif defined(__hpux)

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  LibraryHandle handle = InvalidLibraryHandle;
  
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    handle = shl_load(lib, BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("/")+lib;
      handle = shl_load(path.c_str(), BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
      if(handle != InvalidLibraryHandle)
        return handle;
    }
  }
  
  // one more final attempt
  handle = shl_load(lib, BIND_DEFERRED | DYNAMIC_PATH | BIND_NONFATAL, 0L);
  return handle;
}

CubitString CubitDynamicLoader::get_error()
{
  return CubitString();
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return shl_unload(lib) == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  void* address;
  return shl_findsym(&lib, symbol, TYPE_PROCEDURE, &address) < 0 ? NULL : address;
}

const char* CubitDynamicLoader::library_prefix()
{
  return "lib";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".sl";
}

#else

#include <dlfcn.h>

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::InvalidLibraryHandle = NULL;

CubitDynamicLoader::LibraryHandle CubitDynamicLoader::load_library(const char* lib)
{
  // if absolute path, test it directly
  if(absolute_path(lib))
  {
    return dlopen(lib, RTLD_LAZY);
  }
  else
  { 
    // try finding with our search paths
    for(int i=0; i<(int)gSearchPaths.size(); i++)
    {
      CubitString path = gSearchPaths[i]+CubitString("/")+lib;
      LibraryHandle handle = dlopen(path.c_str(), RTLD_LAZY);
      if(handle != InvalidLibraryHandle)
        return handle;
    }
  }
  
  // one more final attempt
  return dlopen(lib, RTLD_LAZY );
}

CubitString CubitDynamicLoader::get_error()
{
  return dlerror();
}

CubitStatus CubitDynamicLoader::unload_library(CubitDynamicLoader::LibraryHandle lib)
{
  return dlclose(lib) == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void* CubitDynamicLoader::get_symbol_address(CubitDynamicLoader::LibraryHandle lib, const char* symbol)
{
  return dlsym(lib, symbol);
}

const char* CubitDynamicLoader::library_prefix()
{
  return "lib";
}

const char* CubitDynamicLoader::library_extension()
{
  return ".so";
}
#endif

