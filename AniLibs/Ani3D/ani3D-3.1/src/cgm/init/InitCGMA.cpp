#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "VirtualQueryEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "CubitUtil.hpp"
#include "CADefines.hpp"
#include "CubitAttrib.hpp"

#include <ctype.h>

// This is for fixing the segmentation fault with writing Attribute 43 (must be removed for future versions)
class CAAssemblyData : public CubitAttrib
{
public:
  CAAssemblyData(RefEntity*, const CubitSimpleAttrib&) : CubitAttrib(NULL) {}

  CubitSimpleAttrib cubit_simple_attrib() { return CubitSimpleAttrib(NULL, NULL, NULL); }

  CubitStatus actuate() { return CUBIT_SUCCESS; }

  CubitStatus update() { return CUBIT_SUCCESS; }

  int int_attrib_type() { return CA_ASSEMBLY_DATA; }

  CubitStatus reset() { return CUBIT_SUCCESS; }
};

CubitAttrib* CAAssemblyData_creator(RefEntity* entity, const CubitSimpleAttrib& p_csa)
{
  return new CAAssemblyData(entity, p_csa);
}

#ifdef CUBIT_CGM
#include <new>

class AcisModifyEngine
{
public:
  AcisModifyEngine();
  static AcisModifyEngine *instance_;
};

class AcisQueryEngine
{
public:
  AcisQueryEngine();
  static AcisQueryEngine *instance_;
};

class dummym
{
  char dummyvar[4096];
};

#ifdef ACIS_TWEAK_TOOL_CAT
class AcisTweakToolCAT
{
  public:
    static AcisTweakToolCAT* instance();
};
#endif

#elif defined(HAVE_ACIS) 
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#endif

#if defined(HAVE_OCC) 
#include "OCCQueryEngine.hpp"
#include "OCCModifyEngine.hpp"
#endif


static bool streq_nocase( const char* s, const char* t )
{
  for (; *s; ++s, ++t) 
    if (tolower(*s) != tolower(*t))
      return false;
  return !*t;
}

static bool has_been_initialized = false;
static char* first_engine_name = NULL;


CubitStatus InitCGMA::initialize_cgma( const char* default_engine_name )
{
  if( has_been_initialized ){
    // CGM is already initialized.  Return success if previous initialization had
    // the same parameter, failure otherwise.
    if( default_engine_name == first_engine_name ){
      return CUBIT_SUCCESS;
    }
    else if( default_engine_name && first_engine_name && 
	     streq_nocase( default_engine_name, first_engine_name )){
      return CUBIT_SUCCESS;
    }
    else{
      PRINT_ERROR( "initialize_cgma() called again, but default engines differ.\n" );
      return CUBIT_FAILURE;
    }
  }

  std::vector<CubitString> args;
  CGMApp::instance()->startup( args );
  GeometryModifyEngine* default_engine = 0;
  bool ignore_default = false;

  // Now register attribute 43 (must be removed as the bug in CUBIT is fixed in future versions)
  // Register CA_ASSEMBLY_DATA (a dummy attribute)
  CGMApp::instance()->attrib_manager()->register_attrib_type(CA_ASSEMBLY_DATA, "assembly data", "ASSEMBLY_DATA",
    CAAssemblyData_creator, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE, CUBIT_FALSE);

#ifdef CUBIT_CGM
  if (!AcisQueryEngine::instance_)
    AcisQueryEngine::instance_ = new (reinterpret_cast<AcisQueryEngine*>(new dummym)) AcisQueryEngine;
  if (!AcisModifyEngine::instance_)
    AcisModifyEngine::instance_ = new (reinterpret_cast<AcisModifyEngine*>(new dummym)) AcisModifyEngine;
  if (default_engine_name && streq_nocase("ACIS",default_engine_name))
    ignore_default = true; // We cannot set the default engine correctly
                       // because we cannot safely cast to GeometryModifyEngine
                       // without the complete definition of AcisModifyEngine.
                       // But it shouldn't matter, as ACIS is the default for
                       // Cubit anyway.
#ifdef ACIS_TWEAK_TOOL_CAT
  AcisTweakToolCAT::instance();
#endif

#elif defined(HAVE_ACIS)
  AcisQueryEngine::instance();
  AcisModifyEngine::instance();
  if (default_engine_name && streq_nocase("ACIS",default_engine_name))
    default_engine = AcisModifyEngine::instance();
#endif

#ifdef HAVE_OCC  
  OCCQueryEngine::instance();
  OCCModifyEngine::instance();
  if (default_engine_name && streq_nocase("OCC",default_engine_name))
    default_engine = OCCModifyEngine::instance();
#endif  

  FacetQueryEngine::instance();
  FacetModifyEngine::instance();
  VirtualQueryEngine::instance()->register_attributes();

  if (default_engine_name && streq_nocase("FACET",default_engine_name)) {
    default_engine = FacetModifyEngine::instance();
    FacetModifyEngine::instance()->set_modify_enabled(CUBIT_TRUE);
  }

  if(default_engine_name && !ignore_default) {
    if (!default_engine) {
      PRINT_ERROR("Invalid or unsupported engine: '%s'\n", default_engine_name);
      return CUBIT_FAILURE;
    }
    
    CubitStatus rval;
    rval = GeometryModifyTool::instance()->set_default_gme(default_engine);
    if (CUBIT_SUCCESS != rval)
      return rval;
    rval = GeometryQueryTool::instance()->set_default_gqe(default_engine->get_gqe());
    if (CUBIT_SUCCESS != rval)
      return rval;
  }
  
  // set has_been_initialized only if everything worked
  if( default_engine_name ){
    first_engine_name = CubitUtil::util_strdup(default_engine_name); 
  }
  has_been_initialized = true;

  return CUBIT_SUCCESS;
}

CubitStatus InitCGMA::deinitialize_cgma() 
{
  CGMApp::instance()->shutdown();
  has_been_initialized = false;
  return CUBIT_SUCCESS;
}

