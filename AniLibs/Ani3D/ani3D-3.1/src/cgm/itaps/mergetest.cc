/**
 * \file mergetest.cpp
 *
 * \brief mergetest, utility to merge files
 *
 */
#include "iGeom.h"
#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <string.h>
#include <math.h>


#define CHECK( STR ) if (err != iBase_SUCCESS) return print_error( STR, err, geom, __FILE__, __LINE__ )

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S

static bool print_error( const char* desc,
                         int err,
                         iGeom_Instance geom,
                         const char* file,
                         int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';

  std::cerr << "ERROR: " << desc << std::endl
            << "  Error code: " << err << std::endl
            << "  Error desc: " << buffer << std::endl
            << "  At        : " << file << ':' << line << std::endl
            ;

  return false; // must always return false or CHECK macro will break
}

typedef iBase_TagHandle TagHandle;
typedef iBase_EntityHandle GentityHandle;
typedef iBase_EntitySetHandle GentitysetHandle;

/* Frees allocated arrays for us */
template <typename T> class SimpleArray
{
  private:
    T* arr;
    int arrSize;
    int arrAllocated;

  public:
    SimpleArray() : arr(0) , arrSize(0), arrAllocated(0) {}
    SimpleArray( unsigned s ) :arrSize(s), arrAllocated(s) {
      arr = (T*)malloc(s*sizeof(T));
      for (unsigned i = 0; i < s; ++i)
        new (arr+i) T();
    }

    ~SimpleArray() {
      for (int i = 0; i < size(); ++i)
        arr[i].~T();
      free(arr);
    }

    T**  ptr()            { return &arr; }
    int& size()           { return arrSize; }
    int  size()     const { return arrSize; }
    int& capacity()       { return arrAllocated; }
    int  capacity() const { return arrAllocated; }

    typedef T* iterator;
    typedef const T* const_iterator;
    iterator       begin()       { return arr; }
    const_iterator begin() const { return arr; }
    iterator         end()       { return arr + arrSize; }
    const_iterator   end() const { return arr + arrSize; }


    T& operator[]( unsigned idx )       { return arr[idx]; }
    T  operator[]( unsigned idx ) const { return arr[idx]; }
};

#define ARRAY_INOUT( A ) A.ptr(), &A.capacity(), &A.size()
#define ARRAY_IN( A ) &A[0], A.size()

int main( int argc, char *argv[] )
{
    // Check command line arg
  #ifdef FORCE_OCC
    #ifndef HAVE_OCC
      #error "Cannot force use of OCC w/out OCC support"
    #endif
  std::string filename = STRINGIFY(SRCDIR) "/voltest.stp";
  std::string outfile = "merged.occ";
    std::string engine_opt = ";engine=OCC";
  #elif defined(HAVE_ACIS)
    std::string filename = STRINGIFY(SRCDIR) "/voltest.sat";
    std::string outfile = "merged.sat";
    std::string engine_opt = ";engine=ACIS";
  #elif defined(HAVE_OCC)
    std::string filename = STRINGIFY(SRCDIR) "/voltest.stp";
    std::string outfile = "merged.occ";
    std::string engine_opt = ";engine=OCC";
  #else
    std::string filename = STRINGIFY(SRCDIR) "/voltest.sat";
    std::string outfile = "merged.sat";
    std::string engine_opt;
  #endif
  if (argc == 1) {
    std::cout << "Using default input file: " << filename << std::endl;
    std::cout << "Using default output file: " << outfile << std::endl;
  }
  else if (argc >= 3) {
    filename = argv[1];
    outfile = argv[2];
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [geom_filename] [out file]" << std::endl;
    return 1;
  }

  // initialize  geom
  int err;
  iGeom_Instance geom;
  iGeom_newGeom( engine_opt.c_str(), &geom, &err, engine_opt.length() );

    // Print out Header information
  std::cout << "\n merge utility :\n\n";

  iGeom_load( geom, &filename[0], 0, &err, filename.length(), 0 );
  CHECK( "ERROR : can not load a geometry" );

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet( geom, &root_set, &err );
  CHECK( "ERROR : getRootSet failed!" );

    // print out the number of entities
  std::cout << "Model contents: " << std::endl;
  const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};
  for (int i = 0; i <= 3; ++i) {
    int count;
    iGeom_getNumOfType( geom, root_set, i, &count, &err );
    CHECK( "Error: problem getting entities after gLoad." );
    std::cout << gtype[i] << count << std::endl;
  }

  // get volumes, and call merge
  int count;
  iGeom_getNumOfType( geom, root_set, 3, &count, &err );
  CHECK( "Error: problem getting volume numbers after gLoad." );

  SimpleArray<iBase_EntityHandle> vols;

  iGeom_getEntities( geom, root_set, 3, ARRAY_INOUT(vols), &err );
  CHECK( "Error: problem getting volumes." );

  //  now imprint
  std::cout << "\n\nImprinting...." << std::endl;
  iGeom_imprintEnts(geom, ARRAY_IN(vols),&err);
  CHECK( "Error: problem imprinting volumes." );
  ////CHECK("Imprint failed.");

  std::cout << "\n\nMerging...." << std::endl;
  double dTol = 1e-4;
  iGeom_mergeEnts(geom, ARRAY_IN(vols), dTol, &err);
     ////CHECK("Merge failed.");
  CHECK( "Error: problem merging volumes." );

  iGeom_save(geom, outfile.c_str(), NULL, &err,
                strlen(outfile.c_str()), 0);
  CHECK( "Error: problem saving." );

  return 0;
}

