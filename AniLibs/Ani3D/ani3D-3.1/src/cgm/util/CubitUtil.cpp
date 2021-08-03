//-------------------------------------------------------------------------
// Filename      : CubitUtil.cc 
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

#include "CubitUtil.hpp"
#include "CubitString.hpp"
#include "CubitEntity.hpp"
#include "CubitMessage.hpp"

#include "AppUtil.hpp"
#include <ctype.h>
#include <time.h>
#include <sstream>


#ifdef _WIN32
#include "windows.h"
#else
#include <unistd.h>
#include <sys/utsname.h>
#endif

FILE *CubitUtil::fp = NULL;
static int displayDigits = -1;


int CubitUtil::get_digits()
{
  return displayDigits;
}

void CubitUtil::set_digits(int new_digits)
{
  displayDigits = new_digits;
}


void CubitUtil::convert_string_to_lowercase(char *string)
{
   register char *p = string;
   while (*p)
   {
     if (isascii(*p) && isupper(*p))
         *p = tolower(*p);
      p++;
   }
}

int CubitUtil::strcmp_case_insensitive (const char *s1, const char *s2)
{
   char c1, c2;
   
   do
   {
      c1 = *s1++;
      if(isupper(c1))
         c1 = tolower(c1);
      
      c2 = *s2++;
      if(isupper(c2))
         c2 = tolower(c2);
      
      if(c1 != c2)
         return c1 - c2;
      
   } while(c1 != '\0');
   
   return 0;
}

int CubitUtil::strncmp_case_insensitive (const char *s1, const char *s2,
                                         int n)
{
   char c1, c2;
   
   do
   {
      c1 = *s1++;
      if(isupper(c1))
         c1 = tolower(c1);
      
      c2 = *s2++;
      if(isupper(c2))
         c2 = tolower(c2);
      
      if(c1 != c2)
         return c1 - c2;
      
      n--;
   } while(c1 && n > 0);
   
   return 0;
}

void CubitUtil::list_ids( const char *const heading,
                            const DLIList<CubitEntity*> &entity_list,
                            int should_sort, int report_once,
                            int wrap )
{
  if ( entity_list.size() == 0 )
  {
    PRINT_INFO("  No %s.\n", heading );
    return;
  }
  DLIList <int> id_list( entity_list.size() );
  for ( int j = 0; j < entity_list.size(); j++ )
  {
    id_list.append( entity_list[j]->id() );
  }

  if ( id_list.size() == 1 )
  {
    PRINT_INFO("  The 1 %s id is %d.\n", heading, id_list[0]);
    return;
  }
  sort_and_print_ids( heading, id_list, should_sort, report_once, wrap );
}

void CubitUtil::sort_and_print_ids( const char *const heading,
                                      DLIList<int> &id_list,
                                      int should_sort, int report_once,
                                      int wrap )
{
    // sort, if desired
  if ( should_sort ) {
    id_list.sort();
  }

  if ( report_once ) {
    DLIList <int> id_list_2( id_list );
    id_list_2.reset();
    id_list.clean_out();
    id_list.append( id_list_2.get_and_step() );
    for ( int j = id_list_2.size()-1; j--; )
    {
      if ( id_list_2.get() != id_list_2.prev() )
        id_list.append( id_list_2.get() );
      id_list_2.step();
    }
  }

  if( wrap == -1 )
  {
    // print out ranges
    int begin = id_list.get_and_step();
    int end = begin;
    int current = -1;
    PRINT_INFO("  The %d %s ids are %d", id_list.size(), heading, begin);
    for (int i=id_list.size()-1; i > 0; i--) {
      current = id_list.get_and_step();
      if (current == end+1) {
        end++;
      }
      else {
        if (end == begin) {
          PRINT_INFO(", %d", current);
        }
        else if (end == begin+1) {
          PRINT_INFO(", %d, %d", end, current);
        }
        else {
          PRINT_INFO(" to %d, %d", end, current);
        }
        begin = current;
        end = begin;
      }
    }
    if (current == begin + 1)	{
      PRINT_INFO(", %d", current);
    }
    else if (current != begin) {
      PRINT_INFO(" to %d", current);
    }
    PRINT_INFO(".\n");
  }
  else
  {
    char pre_string[67];
    sprintf( pre_string, "  The %d %s ids are: ", id_list.size(),heading );
    CubitUtil::list_entity_ids( pre_string, id_list, wrap, ".\n", CUBIT_FALSE,
                                CUBIT_FALSE );
  }
}

void CubitUtil::list_entity_ids( const char *pre_string, 
                                 const DLIList<CubitEntity*> &entity_list, 
                                 int width, const char *post_string,
                                 int sort, int unique,
                                 int tab, const char *sep_string,
                                 const char *post_string_none )
{
  DLIList <int> id_list( entity_list.size() );
  for ( int i=0; i<entity_list.size(); i++ ) 
    id_list.append( entity_list.next(i)->id() );

  list_entity_ids( pre_string, id_list, width, post_string, sort,
                   unique, tab, sep_string, post_string_none );
}

void CubitUtil::list_entity_ids( const char *pre_string, 
                                 DLIList<int> &id_list,
                                 int width, 
                                 const char *post_string,
                                 int sort, int unique,
                                 int tab_len, const char *sep_string,
                                 const char *post_string_none )
{
  CubitString ret_str;
  process_entity_ids( 1, ret_str, pre_string, id_list, width, post_string,
                      sort, unique, tab_len, sep_string, post_string_none );
}

CubitString CubitUtil::get_entity_ids_str( const char *pre_string, 
                                           DLIList<int> &id_list, 
                                           int width, const char *post_string,
                                           int sort, int unique, int tab_len, 
                                           const char *sep_string,
                                           const char *post_string_none )
{
  CubitString ret_str;

  process_entity_ids( 0, ret_str, pre_string, id_list, width, post_string,
                      sort, unique, tab_len, sep_string, post_string_none );
  return ret_str;
}

void CubitUtil::process_entity_ids( int method,
                                    CubitString &ret_str,
                                    const char *pre_string, 
                                    DLIList<int> &id_list,
                                    int max_len, 
                                    const char *post_string,
                                    int sort, int unique,
                                    int tab_len, const char *sep_string,
                                    const char* post_string_none ) 
{
  // Method: 0 - to a string
  //         1 - to PRINT_INFO
  char temp[200];

  if ( id_list.size() == 0 ) {
    if( method )
      PRINT_INFO("%s%s", pre_string, post_string_none );
    else
    {
      sprintf( temp, "%s%s", pre_string, post_string_none );
      ret_str = temp;
    }
    if( fp )
      fprintf( fp, "%s%s", pre_string, post_string_none );
    return;
  }

  // sort
  if( sort )
  {
    id_list.sort();
    
    // make unique
    if( unique ) 
    {
      int i;
      DLIList <int> id_list_2( id_list );
      id_list_2.reset();
      id_list.clean_out();
      id_list.append( id_list_2.get_and_step() );
      for ( i=id_list_2.size()-1; i--; ) 
      {
        if ( id_list_2.get() != id_list_2.prev() )
          id_list.append( id_list_2.get() );
        id_list_2.step();
      }
    }
  }

  if( max_len < 0 )
    max_len = CUBIT_INT_MAX/2;
    
  // TODO: wrap prestring, if necessary
  if( method )
    PRINT_INFO( "%s", pre_string );
  else
    ret_str = pre_string;
  if( fp )
    fprintf( fp, "%s", pre_string );

  // Keep track of length printed
  int curr_len = strlen(pre_string);
  
  int num = 0;
  int begin = id_list.get();
  int previous = begin;
  int current;
  int comma = 0; // Is comma needed
  int beg_len, prev_len;
  int sep_len = strlen( sep_string );

  // Setup the tab
  char* tab = new char[tab_len+1];
  for( int i=0; i<tab_len; i++ )
     tab[i] = ' ';
  tab[tab_len] = '\0';

  // Loop until all the ids are printed.  Use ranges if possible.
  while( num < id_list.size()+1 )
  {
    current = id_list.get_and_step();
    num++;

    // Handle last entity
    if( num <= id_list.size() )
    {
      if( num==1 ) // Handle 1st time in loop
        continue;
      
      if( current==previous+1 )
      {
        previous = current;
        continue;
      }
    }

    // If we are here, we are no longer tracking a range and
    // need to print the range or a number.
    if( comma )
    {
      if( method )
        PRINT_INFO("%s", sep_string );
      else
        ret_str += sep_string;
      if( fp )
        fprintf( fp, "%s", sep_string );
      curr_len += sep_len;
    }

    if( begin==previous )
    {
      // a single number
      prev_len = int_len(previous);

      if( curr_len+1+prev_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d", tab, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d", tab, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d", tab, previous );
        curr_len = tab_len + prev_len;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d", previous );
        else
        {
          sprintf( temp, "%d", previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d", previous );
        curr_len = curr_len + prev_len;
      }
    }
    else if( previous==begin+1 )
    {
      // a range, but only 2 consecutive numbers
      prev_len = int_len(previous);
      beg_len = int_len(begin);

      // Print 1st
      if( curr_len+1+beg_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d%s", tab, begin, sep_string );
        }
        else
        {
          sprintf( temp, "\n%s%d%s", tab, begin, sep_string );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d%s", tab, begin, sep_string );
        curr_len = tab_len + beg_len + sep_len;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d%s", begin, sep_string );
        else
        {
          sprintf( temp, "%d%s", begin, sep_string );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d%s", begin, sep_string );
        curr_len = curr_len + beg_len + sep_len;
      }

      // Print 2nd
      if( curr_len+1+prev_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d", tab, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d", tab, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d", tab, previous );
        curr_len = tab_len + prev_len;
      }
      else
      {
        if( method )
          PRINT_INFO( " %d", previous );
        else
        {
          sprintf( temp, " %d", previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, " %d", previous );
        curr_len = curr_len + 1+prev_len;
      }
    }
    else
    {
      // a range of 3 or more consecutive numbers
      prev_len = int_len(previous);
      beg_len = int_len(begin);

      if( curr_len+beg_len+prev_len+5+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d to %d", tab, begin, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d to %d", tab, begin, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d to %d", tab, begin, previous );
        curr_len = tab_len + beg_len+prev_len+4;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d to %d", begin, previous );
        else
        {
          sprintf( temp, "%d to %d", begin, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d to %d", begin, previous );
        curr_len = curr_len + beg_len+4+prev_len;
      }
    }

    begin = current;
    previous = current;
    comma = 1;

  }

  //TODO: wrap poststring, if required
  if (post_string) {
    
    if( method )
      PRINT_INFO( "%s", post_string );
    else
      ret_str += post_string;
    if( fp )
      fprintf( fp, "%s", post_string );
  }
  
  delete [] tab;
}

#define INTABS(n) ((n) >= 0 ? (n) : (-(n)))
int CubitUtil::int_len( int num )
{
  int len = 0; // length of the string to hold the integer number
  unsigned long n; // absolute value of the integer value
  
  // If the number is negative, add 1 for the negative sign
  if (num < 0) len++;
  
  // Loop until the absolute value of the number reaches 0
  n = INTABS(num);
  do {
    // Increment the length and divide the number by 10
    len++;
    n /= 10;
  } while (n);

  return len;
}

namespace
{
  // Unix: Returns the TEMPDIR, TMP, or TEMP directory, whichever is set to a
  // writeable directory.  If neither is set, return /tmp.
  // Windows: Return path returned by GetTempPath() windows function.
  // If it doesn't return a writeable directory, use current directory.
  CubitString get_temp_directory()
  {
#ifdef _WIN32

    //get a place to put the temporary file
    const DWORD buf_size = MAX_PATH;
    wchar_t temp_path[buf_size];
    GetTempPathW(buf_size, temp_path);

    // If the path is not writeable, use the current directory instead
    DWORD atts = GetFileAttributesW(temp_path);
#if _MSC_VER > 1200 // after VC6.0
    if (atts == INVALID_FILE_ATTRIBUTES ||      // File doesn't exist
        (atts & FILE_ATTRIBUTE_DIRECTORY) == 0 || // File isn't a directory
        atts & FILE_ATTRIBUTE_READONLY)         // File is read only
#else
    if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0 || // File isn't a directory
        atts & FILE_ATTRIBUTE_READONLY)         // File is read only
#endif
    {
      if (DEBUG_FLAG(141))
      {
        PRINT_DEBUG_141("\nUsing cwd because ");
#if _MSC_VER > 1200
        if (atts == INVALID_FILE_ATTRIBUTES)
        {
          PRINT_DEBUG_141("directory doesn't exist: %s\n", CubitString::toUtf8(temp_path).c_str());
        }
        else if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0)
        {
          PRINT_DEBUG_141("file isn't a directory: %s\n", CubitString::toUtf8(temp_path).c_str());
        }
        else if (atts & FILE_ATTRIBUTE_READONLY)         
        {
          PRINT_DEBUG_141("directory is read only: %s\n", CubitString::toUtf8(temp_path).c_str());
        }
#else
        if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0)
        {
          PRINT_DEBUG_141("file isn't a directory: %s\n", CubitString::toUtf8(temp_path).c_str());
        }
        else if (atts & FILE_ATTRIBUTE_READONLY)         
        {
          PRINT_DEBUG_141("directory is read only: %s\n", CubitString::toUtf8(temp_path).c_str());
        }
#endif
      }
      temp_path[0] = '.';
      temp_path[1] = '\0';
    }
    else
    {
      PRINT_DEBUG_141("\nUsing GetTempPath: %s\n", CubitString::toUtf8(temp_path).c_str());
    }
    return CubitString::toUtf8(temp_path);

#else

    const char* tmpdir = "/tmp";
    const char* env_tmpdir = getenv("TMPDIR");
    if (!env_tmpdir)
      env_tmpdir = getenv("TMP");
    if (!env_tmpdir)
      env_tmpdir = getenv("TEMP");
    if(env_tmpdir)
      tmpdir = env_tmpdir;

    return tmpdir;

#endif
  }
}

CubitString CubitUtil::get_temporary_filename()
{

  CubitString ret_str;

#ifdef _WIN32
  
  //get a place to put the temporary file
  CubitString temp_path = get_temp_directory();

  // make an empty temporary and return the name for it
  wchar_t temp_file_name[MAX_PATH];
  if( GetTempFileNameW(CubitString::toUtf16(temp_path.c_str()).c_str(), L"CBT", 0, temp_file_name) != 0 )
    ret_str = CubitString::toUtf8(temp_file_name); 

#else

  CubitString tmpdir = get_temp_directory();
  const char* filepattern = "CBT.XXXXXX";
    //needs to be two longer because of the "/"?
  char *temp_file_name = new char[tmpdir.length() + strlen(filepattern) + 2];
  sprintf(temp_file_name, "%s/%s", tmpdir.c_str(), filepattern);

  // make an empty file and return the name for it
  int fd = mkstemp(temp_file_name);
  if( fd != -1 )
  {
    ret_str = temp_file_name; 
    // release the open done by mkstemp,
    // temporary file still exists
    close(fd);
  }
  delete [] temp_file_name;

#endif

  return ret_str;
}

void CubitUtil::print_columns( const std::vector<CubitString>& array,
                               const CubitString& indent)
{
  int term_height, term_width;
  
  if( ! AppUtil::instance()->get_terminal_size( term_height, term_width ) )
  {
    for( size_t i = 0; i < array.size(); i++ )
    {
      PRINT_INFO("%s%s\n", indent.c_str(), array[i].c_str() );
    }
    return;
  }
  
    // find lenth of longest string
  int maxlen = 0;
  for( size_t i = 0; i < array.size(); i++ )
  {
    int len = CubitUtil::string_length(array[i].c_str());
    if( len > maxlen )
      maxlen = len;
  }
  
  char* const line = new char[CUBIT_MAX(maxlen,term_width)+2];
  
    // calculate number of columns of output
  term_width -= string_length(indent.c_str());
  int width = maxlen + 1;
  int columns = term_width > width ? term_width / width : 1;
  
    // calculate number of rows of output
  int rows = array.size() / columns;
  if( array.size() % columns )
    rows++;
  
    // calculate the width of one column
  if (columns > 1)
    width = maxlen + (term_width - columns * maxlen) / (columns - 1);
  else
    width = term_width;
  
    // now write output
  for( int i = 0; i < rows; i++ )
  {
    size_t idx;
    const char* str;
    char* ptr = line + sprintf( line, "%s", indent.c_str() );
    for(int j = 0; j < columns - 1; j++ )
    {
      idx = j * rows + i;
      if (idx < array.size() )
        str = array[idx].c_str();
      else
        str = "";
        
      ptr += sprintf( ptr, "%-*s", width, str);
    }
    
    idx = (columns - 1) * rows + i;
    if (idx < array.size() )
      sprintf(ptr, "%s\n", array[idx].c_str());
    else
      sprintf(ptr, "\n");
    
    PRINT_INFO( "%s", line );
  }
  
  delete [] line;
}

int CubitUtil::string_length( const char* string, int tabsize )
{
  int result = 0;
  for( ; *string ; string++ )
  {
    if( *string == '\t' )
      result += tabsize;
    else if( *string >= ' ' )
      result ++;
  }
  return result;
}

//does the same thing as strdup... strdup is not supported by some
// compilers
char* CubitUtil::util_strdup(const char *s1)
{
#ifdef CUBIT_NO_STRDUP
  int len = strlen(s1)+1;
  char* ret_char = (char*) malloc ( (unsigned long) len * sizeof(char));
  strcpy(ret_char, s1);
  return ret_char;
#else
#ifdef _WIN32
  return _strdup(s1);
#else   
  return strdup(s1);
#endif  
#endif  
}

void CubitUtil::cubit_sleep(int duration_in_seconds)
{
#ifdef _WIN32
  ::Sleep(duration_in_seconds*1000);
#else
  sleep(duration_in_seconds);
#endif
}

CubitString CubitUtil::getenv(const CubitString& var)
{
#ifdef _WIN32
  return CubitString::toUtf8(_wgetenv(CubitString::toUtf16(var).c_str()));
#else
  return CubitString(::getenv(var.c_str()));
#endif
}

void CubitUtil::setenv(const CubitString& var, const CubitString& value)
{
#ifdef _WIN32
  CubitString tmp = var + "=" + value;
  _wputenv(CubitString::toUtf16(tmp).c_str());
#else
  ::setenv(var.c_str(), value.c_str(), 1);
#endif
}


CubitString CubitUtil::get_computer_name()
{
  CubitString name = "unknown";
#ifdef _WIN32
  wchar_t machine_buf[MAX_COMPUTERNAME_LENGTH + 1];
  DWORD machine_buf_length = MAX_COMPUTERNAME_LENGTH + 1;
  if (::GetComputerNameW(machine_buf, &machine_buf_length))
    name = CubitString::toUtf8(machine_buf);
#else
  struct utsname uname_data;
  if( uname( &uname_data ) >= 0 )
  {
    name = uname_data.nodename;
  }
#endif  
  return name;
}

CubitString CubitUtil::get_os()
{
#ifdef _WIN32
  CubitString os = "Microsoft Windows";
  OSVERSIONINFO osvi;
  
  ZeroMemory(&osvi, sizeof(OSVERSIONINFO));
  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
  
  if(GetVersionEx(&osvi))
  {
    std::stringstream str;
    str << " ";
    str << osvi.dwMajorVersion << "." << osvi.dwMinorVersion;
    str << " ";
    str << osvi.szCSDVersion;
    DWORD build = osvi.dwBuildNumber & 0xFFFF;
    str << " (Build " << build << ")";
    os += str.str().c_str();
  }
#else
  CubitString os = "unknown";
  struct utsname uname_data;
  if( uname( &uname_data ) >= 0 )
  {
    os = uname_data.sysname;
    os += " ";
    os += uname_data.release;
  }
#endif
  return os;
}

int CubitUtil::num_cpu()
{
  int num_cpu = 0;

  #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo( &sysinfo );
    num_cpu = sysinfo.dwNumberOfProcessors;  
  #else
    num_cpu = sysconf( _SC_NPROCESSORS_ONLN );
  #endif
    
  return num_cpu;
}
