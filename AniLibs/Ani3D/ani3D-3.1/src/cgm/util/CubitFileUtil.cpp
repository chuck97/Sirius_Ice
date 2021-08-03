//- Class:          CubitFileUtil
//- Description:    Class with functions to get files from a directory, etc.
//- Owner:          Steve Storm

#define NOMINMAX

#include "CubitFileUtil.hpp"
#include "CubitUtil.hpp"
#include "CubitString.hpp"
#include "CubitMessage.hpp"
#include "CubitDirIterator.hpp"

#ifdef _WIN32
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <windows.h>
  #ifndef PATH_MAX
    #define PATH_MAX _MAX_PATH
  #endif
  #include "shlwapi.h"
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

#include <algorithm>


#ifdef _WIN32
static const char* DIR_SEP_STR = "\\";
static const char DIR_SEP_CHAR = '\\';
#else
static const char* DIR_SEP_STR = "/";
static const char DIR_SEP_CHAR = '/';
#endif
   
const char* CubitFileUtil::separator()
{
  return DIR_SEP_STR;
}

CubitStatus
CubitFileUtil::get_current_working_directory( CubitString& wd )
{
#ifdef _WIN32
  wchar_t* buffer = _wgetcwd( NULL, 0 );
#else
  char* buffer = getcwd( NULL, 0 );
#endif
  if (!buffer)
  {
    PRINT_WARNING( "Unable to get new working directory\n" );
    return CUBIT_FAILURE;
  }
  else
  {
    // convert to string
#ifdef _WIN32
    wd = CubitString::toUtf8(buffer);
#else
    wd = buffer;
#endif

    // Add a slash at the end, if not already there
    int wd_len = wd.length();
    if( wd.c_str()[wd_len-1] != DIR_SEP_CHAR )
    {
      wd += DIR_SEP_STR;

      free(buffer);
    }

// TODO: need to figure out the right way to do this. This assumes
// variables of length PATH_MAX which is bad!
//#ifdef _WIN32
//    // Make sure format is compatible with full path format
//    CubitString full_path_str;
//    if( get_full_path_str( wd, full_path_str ) == CUBIT_SUCCESS )
//      strcpy( wd, full_path_str.c_str() );
//#endif

    return CUBIT_SUCCESS;
  }
}

CubitStatus CubitFileUtil::set_current_working_directory( const CubitString& wd )
{
#ifdef _WIN32
  int ret = _wchdir(CubitString::toUtf16(wd).c_str());
#else
  int ret = chdir(wd.c_str());
#endif
  return ret == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitString CubitFileUtil::add_name_to_path( const CubitString& path, const CubitString& name )
{
  CubitString result = path;
  // Add a slash at the end of the path, if not already there
  int path_len = result.length();
  if( result.c_str()[path_len-1] != DIR_SEP_CHAR )
  {
    result += DIR_SEP_STR;
  }
  // append the name to the end of the path
  result += name;
      
  return result;
}

CubitString CubitFileUtil::find_home_path(const CubitString& which_user)
{
  CubitString home_dir;

#ifdef _WIN32
  home_dir = CubitUtil::getenv("USERPROFILE");
#else
  if(which_user.length() == 0)
  {
    home_dir = CubitUtil::getenv("HOME");
    if( home_dir.length() == 0 )
    {
      struct passwd* userdata = getpwuid( getuid() );
      if( userdata )
        home_dir = userdata->pw_dir;
    }
  }
  else
  {
    struct passwd* userdata = getpwnam( which_user.c_str() );
    if(userdata)
      home_dir = userdata->pw_dir;
  }
#endif

  return home_dir;
}

CubitStatus
CubitFileUtil::create_directory( const CubitString& wd )
{
  // Create the directory
#ifdef _WIN32
  if (_wmkdir(CubitString::toUtf16(wd).c_str()) == -1)
  {  
    PRINT_WARNING( "Unable to create new directory\n" );
    return CUBIT_FAILURE;
  }
#else    
  if (mkdir(wd.c_str(), 0777) == -1)
  {  
    PRINT_WARNING( "Unable to create new directory\n" );
    return CUBIT_FAILURE;
  }
#endif    
  return CUBIT_SUCCESS;
}


CubitStatus CubitFileUtil::remove_file( const CubitString& file )
{
#ifdef _WIN32
  int status = _wremove(CubitString::toUtf16(file).c_str());
#else
  int status = remove(file.c_str());
#endif
  return status == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus CubitFileUtil::rename_file( const CubitString& old_file, const CubitString& new_file )
{
#ifdef _WIN32
  int status = _wrename(CubitString::toUtf16(old_file).c_str(), CubitString::toUtf16(new_file).c_str());
#else
  int status = rename(old_file.c_str(), new_file.c_str());
#endif
  return status == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus
CubitFileUtil::get_full_path_str( const CubitString& part,
                                  CubitString &full_path_str )
{
  CubitString my_part = CubitFileUtil::make_path_platform_compatible(part);

#ifdef _WIN32

  wchar_t* full = _wfullpath(NULL, CubitString::toUtf16(my_part).c_str(), 0);
  if(!full)
  {
    PRINT_ERROR( "problem getting full path to %s\n", part.c_str() );
    return CUBIT_FAILURE;
  }
  full_path_str = CubitString::toUtf8(full);
  free(full);
  
#else

  // we loop removing parts until realpath can resolve an existing path,
  // then add the non-existing parts back on.

  std::vector<CubitString> split_parts;
  CubitString trypart = part;

  if(!CubitFileUtil::is_absolute(trypart))
  {
    CubitString cwd;
    CubitFileUtil::get_current_working_directory(cwd);
    trypart = CubitFileUtil::add_name_to_path(cwd, trypart);
  }

  char full[PATH_MAX];
  while(trypart.length() && !realpath(trypart.c_str(), full))
  {
    CubitString split_part1, split_part2;
    CubitFileUtil::split_path(trypart, split_part1, split_part2);
    split_parts.push_back(split_part2);
    if(split_part1.length() == 0)
    {
      PRINT_ERROR( "problem getting full path to %s\n", part.c_str() );
      return CUBIT_FAILURE;
    }
    trypart = split_part1;
  }

  full_path_str = full;
  for(size_t i=0; i<split_parts.size(); i++)
  {
    full_path_str += CubitString("/") + split_parts[split_parts.size() - i - 1];
  }

#endif

  return CUBIT_SUCCESS;
}

CubitString
CubitFileUtil::make_path_platform_compatible( const CubitString& path)
{
  CubitString ret = path;
  for(size_t i=0; i<ret.length(); i++)
  {
#ifdef _WIN32
  // Replace '/' with '\\'
    if(ret.get_at(i) == '/' )
      ret.put_at(i, '\\');
#else
     // Replace '\\' with '/'
    if(ret.get_at(i) == '\\' )
      ret.put_at(i, '/');
#endif
  }
  return ret;
}

CubitString
CubitFileUtil::get_nice_filename( const CubitString& path )
{
  CubitString ret_str;

  CubitString dpart_in, fpart_in;
  split_path( path, dpart_in, fpart_in );

  CubitString wd;
  get_current_working_directory( wd );

  if( dpart_in == wd )
    ret_str = fpart_in;
  else
    ret_str = path;

  return ret_str;
}

void
CubitFileUtil::split_path( const CubitString& path, CubitString& dirpart, CubitString& filepart )
{
  CubitString mypath = path;
  while(mypath.length() && mypath.get_at(mypath.length()-1) == DIR_SEP_CHAR)
  {
    mypath = mypath.substr(0, mypath.length()-1);
  }
  size_t pos = mypath.find_last(DIR_SEP_CHAR);

  // No separator - could be filename or directory.  We assume
  // it's a directory.
  if(pos == CubitString::npos)
  {
    filepart = ".";
    dirpart = mypath;
  }
  else
  {
    filepart = mypath.substr(pos+1);
    dirpart = mypath.substr(0, pos);
  }

  // Add slash on end of dirpart if not already there
  if(dirpart.length() && dirpart.get_at(dirpart.length()-1) != DIR_SEP_CHAR)
    dirpart += DIR_SEP_STR;

  return;
}


CubitString
CubitFileUtil::get_file_extension(
  const CubitString& file,
  bool remove_version /* remove .1, .2, ...*/ )
{
  size_t dot_pos = file.find_last('.');
  size_t dot_pos2 = 0;

  if ( dot_pos == CubitString::npos )
    return "";

  if(remove_version)
  {
	  dot_pos2 = file.find_last('.',dot_pos);
    if ( dot_pos2 == CubitString::npos )
		  remove_version = false;
	  else if(!is_int_number( file.substr(dot_pos+1).c_str() ))
		  remove_version = false;
  }

  CubitString extension;
  if(!remove_version)
	  extension = file.substr(dot_pos);
  else
	  extension = file.substr(dot_pos2,dot_pos-dot_pos2);
  
  for ( size_t i = 0; i < extension.length(); i++ )
    extension.put_at(i, tolower( extension.get_at(i) ) );
    
  return extension;
  
}  //  get_file_extension()
                     
CubitBoolean 
CubitFileUtil::all_chars_are( char ch, const char *str )
{
  while( *str )
    if( *str++ != ch )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::is_int_number( const char *str )
{
  while( *str )
    if( !isdigit(*str++) )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::contains_char( char ch, const char *str )
{
  while( *str )
    if( *str++ == ch )
      return CUBIT_TRUE;
    return CUBIT_FALSE;
}

//CubitString
//CubitFileUtil::get_default_cubit_file_name()
//{
//  // Get a list of all files matching 'cubit*.cub', in current directory

//  CubitDirIterator dir_iter(".", "cubit*.cub");

//  int max_number = 0;

//  while(dir_iter.has_next())
//  {
//    CubitString f = dir_iter.next();

//    // cut off the "cubit" at the start and the ".cub" at the end
//    f = f.substr(5, f.length()-9);
//    if( is_int_number(f.c_str()) )
//    {
//      int num = atoi(f.c_str());
//      max_number = std::max(num, max_number);
//    }
//  }

//  max_number++;

//  CubitString num_str = CubitString::number(max_number);
//  if(max_number < 10)
//    num_str = CubitString("0") + num_str;

//  return CubitString("cubit") + num_str + CubitString(".cub");
//}

int
CubitFileUtil::get_next_filenumber( const CubitString& file_pattern,
                                    int &num_matches,
                                    int &first_id )
{
  // TODO: check this bug?
  // Continue if nothing after dot.  Also, I noticed the _WIN32 will match
  // file.cub against wildcard file.cub.*

  int max_number = 0;
  int min_number = CUBIT_INT_MAX;
  num_matches = 0;

  // get directory and file parts
  CubitString dir, file, full_file_pattern;
  CubitFileUtil::get_full_path_str(file_pattern, full_file_pattern);
  CubitFileUtil::split_path(full_file_pattern, dir, file);

  size_t wildcard_location = file.find("*");
  if(wildcard_location == CubitString::npos)
    return 1;

  CubitDirIterator dir_iter(dir, file);

  while(dir_iter.has_next())
  {
    CubitString f = dir_iter.next();

    // cut off the matched part
    f = f.substr(wildcard_location, 1 + f.length() - file.length());
    if( is_int_number(f.c_str()) )
    {
      int num = atoi(f.c_str());
      max_number = std::max(num, max_number);
      min_number = std::min(num, min_number);
    }
    num_matches++;
  }

  if( min_number != CUBIT_INT_MAX )
    first_id = min_number;

  return max_number+1;
}

bool CubitFileUtil::is_directory(const CubitString& path)
{
  off_t size;
  time_t time;
  int mode = 0;
  if(0 == file_info(path, size, time, mode))
  {
#ifdef _WIN32
  if( (_S_IFDIR & mode) )
#else
  if( S_ISDIR( mode ) )
#endif
    {
      return true;
    }
  }
  return false;
}

bool CubitFileUtil::path_exists(const CubitString& path)
{
  off_t size;
  time_t time;
  int mode;
  return 0 == file_info(path, size, time, mode);
}
   
bool CubitFileUtil::is_absolute(const CubitString& path)
{
#ifdef _WIN32
  return !PathIsRelativeW(CubitString::toUtf16(path).c_str());
#else
  return path.c_str()[0] == '/' ? true : false;
#endif
}
   
int  CubitFileUtil::file_info(const CubitString& path, off_t& size, time_t& time, int& mode)
{
#ifdef _WIN32
  // remove trailing separators
  CubitString mypath = path;
  if(mypath.length() >= 1 && mypath.get_at(mypath.length()-1) == DIR_SEP_CHAR)
    mypath = mypath.substr(0, mypath.length()-1);
  struct _stati64 file_info;
  int stat_result = _wstati64( CubitString::toUtf16(mypath).c_str(), &file_info );
#else
  struct stat file_info;
  int stat_result = lstat( path.c_str(), &file_info );
#endif

  if(stat_result == 0)
  {
    size = file_info.st_size;
    time = file_info.st_mtime;
    mode = file_info.st_mode;
    return 0;
  }
  return errno ? errno : ENOENT;
}

int CubitFileUtil::complete_filename(CubitString& line, int& num_additional_chars, bool& found_quote)
{
  // save the original length
  num_additional_chars = line.length();
  char* buffer = (char*) line.c_str();

  char *ptr_from;
  if ( (ptr_from = strrchr(buffer, '\'')) == NULL &&
       (ptr_from = strrchr(buffer, '\"')) == NULL)
  {
    return -1;
  }
  
  ptr_from++; // skip past the single/double quote
  
    // Separate the directory portion from the filename portion (if any)
  char *file;
  char *path;

#ifndef _WIN32
  if ((file = strrchr(ptr_from, '/')) == NULL) // No path
  {
    path = strdup(".");
    file = strdup(ptr_from);
  }
  else if (file == ptr_from)
  {
    path = strdup("/");
    file = strdup(ptr_from + 1);
  }
  else
  {
    path = strdup(ptr_from);
    char *end = strrchr(path, '/');
    *end = '\0';
    ++file;
    file = strdup(file);

    if ( *path == '~' )
    {
      const char *home = getenv("HOME");
      if( !home )
      {
        struct passwd* userdata = getpwuid( getuid() );
        if( userdata )
          home = userdata->pw_dir;
      }
      char *after_tilde = path;
      after_tilde++;
      std::string new_path = home;
      new_path += after_tilde;
      free(path);
      path = strdup( new_path.c_str() );
    }

  }
#else
  file = strrchr(ptr_from, '\\');
  bool slash = 0;
  if(file == NULL)
  {
    slash = 1;
    file = strrchr(ptr_from, '/');
  }
  if(file == NULL)
//   if ((file = strrchr(ptr_from, '\\')) == NULL
//       && file = strrchr(ptr_from, '/') == NULL) // No path
  {
    path = strdup(".");
    file = strdup(ptr_from);
  }
  else if (file == ptr_from)
  {
    path = strdup("\\");
    file = strdup(ptr_from + 1);
  }
  else
  {
    path = strdup(ptr_from);
    char *end;
    if(slash)
       end = strrchr(path, '/');
    else
       end = strrchr(path, '\\');
    
    *end = '\0';
    ++file;
    file = strdup(file);
  }
#endif  
    // Determine what type of files to match (suffix_strings)
    // Need to lowercase keyword and identifier.....
  char *tmp_str = strdup(buffer);
  char *keyword = strtok(tmp_str, " \t");
  CubitUtil::convert_string_to_lowercase(keyword);
  int lenkey = strlen(keyword);
  std::vector<std::string> suffix_strings;
  
  if ((strncmp("import", keyword, lenkey) == 0) ||
	  (strncmp("export", keyword, lenkey) == 0))
  {
    char *identifier = strtok(NULL, " \t");
    CubitUtil::convert_string_to_lowercase(identifier);
    int lenid  = strlen(identifier);

      //Let's set up a vector of strings
    std::vector<std::string> acis_strings, cubfile_strings,iges_strings,
       catia_strings,proe_strings,step_strings,fastq_strings, mesh_strings,
        facet_strings, avs_strings, stl_strings, ideas_strings, abaqus_strings,
        nastran_strings, presto_strings;

    acis_strings.push_back(".sat");
    acis_strings.push_back(".sab");
    acis_strings.push_back(".acis");
    acis_strings.push_back(".sat.");

    cubfile_strings.push_back(".cub");

    iges_strings.push_back(".igs");
    iges_strings.push_back(".iges");
    /*
#ifdef ACIS_CATIA_TRANSLATOR
    catia_strings.push_back(".model");
    catia_strings.push_back(".mod");
    catia_strings.push_back(".exp");
    catia_strings.push_back(".div");
#endif */

#ifdef GRANITE  
    proe_strings.push_back(".prt");
    proe_strings.push_back(".asm");
    proe_strings.push_back(".prt.");
    proe_strings.push_back(".asm.");
    proe_strings.push_back(".g");
#endif  
#ifdef CATIA
    catia_strings.push_back(".catpart");
    catia_strings.push_back(".catproduct");
    catia_strings.push_back(".ncgm");    
#endif  

    step_strings.push_back(".stp");
    step_strings.push_back(".step");

    fastq_strings.push_back(".fsq");

    mesh_strings.push_back(".g");
    mesh_strings.push_back(".e");
    mesh_strings.push_back(".exo");
    mesh_strings.push_back(".exoII");
    mesh_strings.push_back(".gen");

    facet_strings.push_back(".facets");
    facet_strings.push_back(".facet");
    facet_strings.push_back(".fac");
    facet_strings.push_back(".off");
    facet_strings.push_back(".OFF");
    
    avs_strings.push_back(".avs");

    stl_strings.push_back(".stl");

    ideas_strings.push_back(".unv");

    abaqus_strings.push_back(".inp");

    nastran_strings.push_back(".bdf");

    presto_strings.push_back(".i");
    
    if (strncmp("acis", identifier, lenid) == 0)
    {
      suffix_strings.assign(acis_strings.begin(),acis_strings.end());
    }
    else if (strncmp("cubit", identifier, lenid) == 0)
    {
      suffix_strings.assign(cubfile_strings.begin(), cubfile_strings.end());
    }
    else if (strncmp("iges", identifier, lenid) == 0)
    {
      suffix_strings.assign(iges_strings.begin(), iges_strings.end());
    }
    
    /*
#ifdef ACIS_CATIA_TRANSLATOR
    else if (strncmp("catia", identifier, lenid) == 0)
    {
      suffix_strings.assign(catia_strings.begin(), catia_strings.end());
    }
#endif */
#ifdef CATIA
    else if (strncmp("catia", identifier, lenid) == 0)
    {
      suffix_strings.assign(catia_strings.begin(), catia_strings.end());
    }
#endif

#ifdef GRANITE  
    else if (strncmp("proe", identifier, lenid) == 0)
    {
      suffix_strings.assign(proe_strings.begin(), proe_strings.end());
    }
#endif  
     else if (strncmp("step", identifier, lenid) == 0)
    {
      suffix_strings.assign(step_strings.begin(), step_strings.end());
    }
    else if (strncmp("fastq", identifier, lenid) == 0)
    {
      suffix_strings.assign(fastq_strings.begin(), fastq_strings.end());
    }
    else if ( strncmp("mesh", identifier, lenid) == 0 ||
              strncmp("free", identifier, lenid) == 0 )
    {
      suffix_strings.assign(mesh_strings.begin(), mesh_strings.end());
    }
    else if (strncmp("facets", identifier, lenid) == 0)
    {
      suffix_strings.assign(facet_strings.begin(), facet_strings.end());
    }
    else if (strncmp("avs", identifier, lenid) == 0)
    {
      suffix_strings.assign(avs_strings.begin(), avs_strings.end());
    }
    else if (strncmp("stl", identifier, lenid) == 0)
    {
      suffix_strings.assign(stl_strings.begin(), stl_strings.end());
    }
    else if (strncmp("ideas", identifier, lenid) == 0)
    {
      suffix_strings.assign(ideas_strings.begin(), ideas_strings.end());
    }
    else if (strncmp("abaqus", identifier, lenid) == 0)
    {
      suffix_strings.assign(abaqus_strings.begin(), abaqus_strings.end());
    }
    else if (strncmp("nastran", identifier, lenid) == 0)
    {
      suffix_strings.assign(nastran_strings.begin(), nastran_strings.end());
    }
    else if (strncmp("presto", identifier, lenid) == 0)
    {
      suffix_strings.assign(presto_strings.begin(), presto_strings.end());
    }
    else //user didn't specify anything, so let's add them all
    {
      suffix_strings.assign(acis_strings.begin(),acis_strings.end());
      suffix_strings.insert(suffix_strings.end(),cubfile_strings.begin(), cubfile_strings.end());
      suffix_strings.insert(suffix_strings.end(),iges_strings.begin(), iges_strings.end());
      /*
#ifdef   ACIS_CATIA_TRANSLATOR
      suffix_strings.insert(suffix_strings.end(),catia_strings.begin(), catia_strings.end());
#endif */

#ifdef CATIA
      suffix_strings.insert(suffix_strings.end(),catia_strings.begin(), catia_strings.end());
#endif 

#ifdef GRANITE  
      suffix_strings.insert(suffix_strings.end(),proe_strings.begin(), proe_strings.end());
#endif  
      suffix_strings.insert(suffix_strings.end(),step_strings.begin(), step_strings.end());
      suffix_strings.insert(suffix_strings.end(),fastq_strings.begin(), fastq_strings.end());
      suffix_strings.insert(suffix_strings.end(),mesh_strings.begin(), mesh_strings.end());
      suffix_strings.insert(suffix_strings.end(),facet_strings.begin(), facet_strings.end());
      suffix_strings.insert(suffix_strings.end(),avs_strings.begin(), avs_strings.end());
      suffix_strings.insert(suffix_strings.end(),stl_strings.begin(), stl_strings.end());
      suffix_strings.insert(suffix_strings.end(),ideas_strings.begin(), ideas_strings.end());
      suffix_strings.insert(suffix_strings.end(),abaqus_strings.begin(), abaqus_strings.end());
      suffix_strings.insert(suffix_strings.end(),nastran_strings.begin(), nastran_strings.end());
    }
  }
  else if (strncmp("open", keyword, lenkey) == 0)
  {
    suffix_strings.push_back(".cub");  // Not an actual suffix, just an indicator
  }
  else if (strncmp("playback", keyword, lenkey) == 0 ||
           strncmp("record", keyword, lenkey) == 0)
  {
    suffix_strings.push_back(".jou");
    suffix_strings.push_back(".test");
    suffix_strings.push_back(".cubit");
  }
  else if (strncmp("hardcopy", keyword, lenkey) == 0)
  {
    suffix_strings.push_back(".ps");
  }
  else if (strncmp("cd", keyword, lenkey) == 0)
  {
    suffix_strings.push_back("/");  // Not an actual suffix, just an indicator
  }
  else
  {
    suffix_strings.push_back("\0");
  }
    
  free(tmp_str);
  
  CubitString tmp_line(buffer);
  CubitString match = CubitFileUtil::list_matching_files(path, file, suffix_strings, tmp_line);
  
  if (match.length())
  {
    int file_length = strlen(file);
    match = match.c_str() + file_length;
    tmp_line += match;
    if (tmp_line.find("\'") != CubitString::npos)
	    tmp_line += "\'";
    else if (tmp_line.find("\"") != CubitString::npos)
      tmp_line += "\" ";
    else
    {
      PRINT_ERROR("INTERNAL ERROR: CubitFileUtil::complete_file -- Could not find "
                  "quote type.");
      tmp_line += "\' ";
    }
    //gl_in_quoted_string = !gl_in_quoted_string;
    //gl_fixup(gl_pos, gl_pos+additional_char);
    found_quote = !found_quote;
  }
  line = tmp_line;
  free(path);
  free(file);
  //gl_redraw();
  num_additional_chars = line.length() - num_additional_chars;
  
  return 1;
}

CubitString CubitFileUtil::list_matching_files(const char *path, const char *file, std::vector<std::string> suffixes, CubitString& line)
{
  CubitString filename;
  size_t      len;
  int         match = 0;
  int         width = 0;
  CubitString match_file;
  CubitString char_match;
  CubitBoolean list_files = CUBIT_FALSE;
  CubitBoolean dir_only = CUBIT_FALSE;
  CubitBoolean is_dir = CUBIT_FALSE;

  CubitDirIterator dirp(path);

  if (!dirp.has_next())
  {
    PRINT_INFO("\n");
    PRINT_ERROR("Invalid Directory: '%s'\n",
                path);
    return "";
  }

    // If we're looking for directories only, indicate it
  if (suffixes.size() > 0 && suffixes.front() ==  "/")
  {
    dir_only = CUBIT_TRUE;
  }
  
  len = strlen(file);
  while (dirp.has_next())
  {
    filename = dirp.next();

    // Skip the . and .. entries in the directory
    if (filename.c_str()[0] == '.' &&
        (filename.c_str()[1] == '\0'))
      continue;
      // Skip files if we only want directories

    is_dir = CubitFileUtil::is_directory(CubitFileUtil::add_name_to_path(path, filename));
        
    if (dir_only && !is_dir)
       continue;

      // Skip what doesn't match
    if (len && strncmp(filename.c_str(), file, len) != 0)
      continue;
    
      // We have a match, see if the suffix (if any) matches
    int sat_number_found = 0;

      // Only match ".prt.", ".sat." or ".asm." if they are followed by just digits.
    std::vector<std::string>::iterator suff;
    for(suff = suffixes.begin(); suff != suffixes.end(); suff++)
    {
      
      if(*suff == ".prt." || *suff == ".sat." || *suff == ".asm.")
      {        
        const char* c_ptr = strstr(filename.c_str(), suff->c_str());
        if(c_ptr)
        {
          c_ptr += 5;
            // Make sure the rest of the filename is just digits
          sat_number_found = 1;
          while(*c_ptr != '\0')
          {
            if (!isdigit(*c_ptr))
               sat_number_found = 0;
            c_ptr++;
          }
        }
      }
    }
    
    int found = 0;
    const char* dot = strrchr(filename.c_str(), '.');
    if(dot != NULL && !is_dir)
    {
      for(suff = suffixes.begin(); suff != suffixes.end(); suff++)
      {
#ifdef _WIN32
        if(stricmp(suff->c_str(), dot) == 0)
#else
        if(strcasecmp(suff->c_str(), dot) == 0)
#endif
        {
          found = 1;  //if the suffix matches
          break;
        }
      }
    }
    else  //If the filename has no suffix, then we'll add it
       found = 1;
            
      
    if(suffixes.begin()->size() == 0 || // If no suffix was sent in
        found == 1 || //or suffix matches or the filename has no suffix
        sat_number_found == 1)  //or if the suffix is something like ".sat."
    {
        // Save the filename of the first match
      if (++match == 1)
      {
        match_file = filename;
        if (list_files == CUBIT_FALSE)
          char_match = filename.c_str()+len;
      }
        // list_files is FALSE until we have found 2 or more matches
        // that have no initial characters in common
      else if (list_files == CUBIT_FALSE)
      {
          // If the initial character doesn't match, list instead.
        if (char_match.c_str()[0] != filename.c_str()[len])
        {
            // Indicate that we are listing
          list_files = CUBIT_TRUE;
            // Indicate that we are starting over
          dirp.open(path);
          match = 0;
            // Undo what has already been done
          continue;
        }
        
          // If there is a match, see how many characters match.
        for (unsigned int i = 0; i < char_match.length(); i++)
        {
          if (char_match.c_str()[i] != filename.c_str()[i+len])
          {
              // Chop off end of char_match
            char_match.put_at(i, '\0');
            break;
          }
        }
      }
        // If two or more matches, we want to 
        // output a list, so put out header and first two matches.
      else if (match == 2)
      {
        width = match_file.length() + filename.length() + 6;
          // See if we need to add a '/' to the first name
        char first_str[2];
        first_str[0] = first_str[1] = '\0';
        if(CubitFileUtil::is_directory(CubitFileUtil::add_name_to_path(path, match_file)))
          first_str[0] = '/';
        
          // Get ready to test second file
        bool tmp_is_dir = CubitFileUtil::is_directory(CubitFileUtil::add_name_to_path(path, filename));

        PRINT_INFO("\n\nPossible filename matches:\n  %s%s  %s%s",
                   match_file.c_str(), first_str, filename.c_str(),
                   tmp_is_dir ? "/" : "");
      }
      else  // This is at least the third match, and we
          // need to list the file names
      {
        bool tmp_is_dir = CubitFileUtil::is_directory(CubitFileUtil::add_name_to_path(path, filename));
        int dir_adjust = tmp_is_dir ? 1 : 0;
        width += filename.length() + 2 + dir_adjust;
        if (width > 79)
        {
          PRINT_INFO(" \n" );
          width = filename.length() + 2 + dir_adjust;
        }
        PRINT_INFO("  %s%s",  filename.c_str(),
                   dir_adjust ? "/" : "");
      }
    }
  }

  if (match == 1)
  {
      // See if match_file is a directory

    if(CubitFileUtil::is_directory(CubitFileUtil::add_name_to_path(path, match_file)))
    {
        // If our only match is a directory,
        // Add the directory name to the buffer, but
        // return NULL.
        // Otherwise, return match_file.
      
      // Look at 'line'. Find the first occurance of a '\'. Use that as
      // the suffix for the line if found. Otherwise use '/'.

      CubitString suffix = "/";
      if (line.find_first('\\') != CubitString::npos)
        suffix = "\\";

      line += match_file.substr(len) + suffix;

      match_file = "";
    }
  }
  else if (match > 1)
  {
      // If we didn't list files,
      // add the common chars to the buffer
    if (!list_files)
	  line += char_match;
    PRINT_INFO(" \n"  );
    match_file = "";
  }
  
  // Return the file name.
  // If it was a directory or there was more than one match, it returns an empty string.
  return match_file;
}

