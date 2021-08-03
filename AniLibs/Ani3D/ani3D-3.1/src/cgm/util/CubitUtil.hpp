//-------------------------------------------------------------------------
// Filename      : CubitUtil.hpp
//
// Purpose       : This file contains utility functions that can be used
//                 throughout Cubit.
//
// Special Notes : This is a pure virtual class, to prevent instantiation.
//                 All functions are static, called like this:
//                 CubitUtil::function_name();
//
// Creator       : Darryl Melander
//
// Date          : 06/08/98
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------


#ifndef CUBIT_UTIL_HPP
#define CUBIT_UTIL_HPP

#include <cstring>
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CubitUtilConfigure.h"
#include "CubitString.hpp"

class CubitEntity;

class CUBIT_UTIL_EXPORT CubitUtil
{
public:
  
  static void set_digits(int new_digits);
  static int get_digits();

   static void convert_string_to_lowercase (char *string);
   static int  strcmp_case_insensitive     (const char *s1, const char *s2);
   static int  strncmp_case_insensitive    (const char *s1, const char *s2,
                                            int n);

  static void list_ids( const char *const heading,
                  const DLIList<CubitEntity*> &entity_list,
                  int should_sort = CUBIT_FALSE,
                  int report_once = CUBIT_FALSE,
                  int wrap = 80 );

  static void sort_and_print_ids( const char *const heading,
                           DLIList<int> &id_list,
                           int should_sort = CUBIT_FALSE,
                           int report_once = CUBIT_FALSE,
                           int wrap = 80);



  static void list_entity_ids( const char *pre_string, 
                               const DLIList<CubitEntity*> &entity_list, 
                               int width = 80, const char *post_string = "\n",
                               int sort = CUBIT_TRUE, int unique = CUBIT_TRUE,
                               int tab = 3, const char *sep_string = ",",
                               const char *post_string_none = "none\n");
  
  static void list_entity_ids( const char *pre_string, 
                               DLIList<int> &id_list,
                               int width = 80, 
                               const char *post_string = "\n",
                               int sort = CUBIT_TRUE, int unique = CUBIT_TRUE,
                               int tab_len = 3, const char *sep_string = ",",
                               const char *post_string_none = "none\n");
  
  static CubitString get_entity_ids_str( const char *pre_string, DLIList<int> &int_list, 
                                         int width = 80, const char *post_string = "\n",
                                         int sort = CUBIT_TRUE, int unique = CUBIT_TRUE,
                                         int left_tab = 3, const char *sep_string = ",",
                                         const char *post_string_none = "none\n");

  static void process_entity_ids( int method,
                                  CubitString &ret_str,
                                  const char *pre_string, 
                                  DLIList<int> &id_list,
                                  int max_len, 
                                  const char *post_string,
                                  int sort, int unique,
                                  int tab_len, const char *sep_string,
                                  const char* post_string_none = "none\n");
  
  static int int_len( int num );
    //- Finds the number of spaces required to print an integer number

  static void set_file_ptr( FILE* file_ptr );
  static void reset_file_ptr();
  static FILE* get_file_ptr();
  //- Used to optionally have list_entity_ids or get_entity_ids_str output dumped
  //- to a file as well.

  static CubitBoolean compare( const char* a, const char* b )
    { return ( a != NULL && b != NULL && 0 == strcmp(a, b) ); }

  static CubitSense opposite_sense(CubitSense sense);
    //- return the sense opposite from sense

  static CubitString get_temporary_filename();
    //- function to get a temporary filename, the file is created empty before return

  static int string_length( const char* string, int tabsize = 8 );
  static void print_columns( const std::vector<CubitString>& array,
                             const CubitString& indent = CubitString());

    //does the same thing as strdup... strdup is not supported by some
    // compilers
  static char* util_strdup(const char *s1);
    //a corresponding free for the above function... mainly calls free().
  static void util_strdup_free(char* s1){free(s1);}

  static void cubit_sleep(int duration_in_seconds);

  static CubitBoolean file_exist(const char* buffer);
  static CubitBoolean file_exist(const CubitString& buffer);

  // get an environment variable
  static CubitString getenv(const CubitString& str);
  static void setenv(const CubitString& var, const CubitString& value);

  static CubitString get_computer_name();
  static CubitString get_os();

  static int num_cpu();
  
private:
  CubitUtil(){}

  static FILE *fp;

};

inline void
CubitUtil::set_file_ptr( FILE* file_ptr )
{fp=file_ptr;}

inline void
CubitUtil::reset_file_ptr()
{fp=NULL;}

inline CubitSense CubitUtil::opposite_sense( CubitSense sense )
{
   assert( sense == CUBIT_UNKNOWN ||
           sense == CUBIT_FORWARD ||
           sense == CUBIT_REVERSED );
   if ( sense == CUBIT_UNKNOWN ) return CUBIT_UNKNOWN;
   else {
      CubitSense opp_sense = (CubitSense) (1 - sense);
      return opp_sense;
   }
}

inline FILE*
CubitUtil::get_file_ptr()
{return fp;}
#endif

