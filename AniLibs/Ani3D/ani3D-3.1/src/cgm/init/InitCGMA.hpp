#ifndef INIT_CGMA_HPP
#define INIT_CGMA_HPP

#include "CgmaInitConfigure.h"
#include "CubitDefines.h"

class CGMA_INIT_EXPORT InitCGMA
{
public:

  static CubitStatus initialize_cgma( const char* default_engine_name = 0 );

  static CubitStatus deinitialize_cgma();
};

#endif
