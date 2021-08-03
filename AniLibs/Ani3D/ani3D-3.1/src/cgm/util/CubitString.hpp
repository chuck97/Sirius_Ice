//- Class: CubitString
//- 
//- Description: This file defines the CubitString class which is a 
//- simple implementation of a character string. Basic functionality 
//- is provided as well as equality/inequality tests, subscripting into 
//- a string, searching a string for another string.
//-

#if !defined(CUBITSTRING_HPP)
#define CUBITSTRING_HPP

#include <sstream>
#include <vector>
#include <string>

#include "CubitUtilConfigure.h"

// GCC compiler checking of format arguments
// Note:  For class member functions, the first argument is the
//        implicit 'this' pointer, so if the format string is the
//        first explicit argument A should be 2, not 1.
#ifdef __GNUC__
# define CUBIT_STRING_PRINTF_FORMAT(A,B) __attribute__ ((format (printf, A, B)))
#else
# define CUBIT_STRING_PRINTF_FORMAT(A,B)
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4251)  // hide warnings about template dll exports
#endif

//! String class that represents a UTF-8 string
class CUBIT_UTIL_EXPORT CubitString
{
public:
  //! max size of a string
  static const size_t npos;

  enum CaseSensitivity
  {
    CaseInsensitive,
    CaseSensitive
  };

  //! Default constructor
  CubitString();

  //! Default destructor
  ~CubitString();

  //! Copy Constructor
  CubitString(const CubitString& s);

  //! Create a string from a char*
  CubitString(const char* s);

  //- Create a string from a std::string
  CubitString(const std::string& s);

  CubitString(int i,char c);
  //- Create a string from a char

  //! Create a string from a integer or other stringstream recognized type
  template <typename T>
  static CubitString number(const T& i)
  {
    std::stringstream si;
    si << i;
    return CubitString(si.str().c_str());
  }

  //! Create a value from a CubitString.
  template <typename T>
  static T toNumber(const CubitString& str, bool *ok = NULL)
  {
    std::istringstream ss(str.c_str());
    T result;
    bool temp_ok = ss >> result;
    if(ok)
      (*ok) = temp_ok;
    return result;
  }

  //! Create a string from a double.
  //! Use either fixed point or scientific notation, whichever is shorter.
  //! s_length is the maximum string length: If s_length > 0, then
  //! string will contain no spaces and be close to s_length long
  //! without going over. Hence precision is variable.
  static CubitString number(double f, unsigned int s_length = 0, unsigned int sig_digits = 0);

  //! set a formated string with arguments
  static CubitString format(const char* fmt, ...) CUBIT_STRING_PRINTF_FORMAT(1,2);

  CubitString& operator=(const CubitString& s);
  
  CubitString& operator+=(const CubitString& s);

  CUBIT_UTIL_EXPORT friend CubitString operator+(const CubitString& s1, const CubitString& s2);

  CUBIT_UTIL_EXPORT friend bool operator<=(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator>=(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator<(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator>(const CubitString&, const CubitString&);
  
  bool operator==(const CubitString &s) const;
  bool operator!=(const CubitString &s) const;
  //- Predicates: Compare equality/inequality with other CubitStrings or
  //- with pointers to a list of characters.

  //! compare this string with another and return whether they are equal.
  //! one can also specify case sensitivity
  bool compare(const CubitString& s, CaseSensitivity cs = CaseSensitive);
  // note case insensitive comparison with unicode is not yet supported

  //! get a character at a position
  char get_at(size_t pos) const;
  //! set a character at a position
  void put_at(size_t pos, char c);

  //! get a substring starting at first and counting
  CubitString substr(size_t first, size_t count = CubitString::npos) const;

  //! remove len number of characters starting at a position pos.
  void erase(size_t pos, size_t len = npos);

  //! replace a portion of ths string with another string
  void replace(size_t pos, size_t len, const CubitString& str);

  //! replace every occurance of string with another string
  void replace(const CubitString& to_find, const CubitString& to_replace);

  //! convert to lower case (Note: this is not internationalized)
  void to_lower();
  static void to_lower(char *string);
  //! convert to upper case (Note: this is not internationalized)
  void to_upper();
  static void to_upper(char *string);

  //! trim off whitespace from the beginning and end of the string
  void trim();

  //! trim whitespace from the beginning and end of string, and
  //! each sequence of internal white space is replaced with a single space
  void simplify();

  //! split a string given a delimiter
  void tokenize( char delimiter, std::vector<CubitString> &strings ) const;

  
  //! find functions
  size_t find(const CubitString& s, size_t pos = 0) const;
  size_t find_first_of(const CubitString& s, size_t pos = 0) const;
  size_t find_first(char c, size_t pos = 0) const;
  size_t find_last (char c, size_t pos = CubitString::npos) const;


  //! return whether the string ends with str
  bool ends_with(const CubitString& str) const;
  //! return whether the string starts with str
  bool starts_with(const CubitString& str) const;

  //! returns whether the string is empty
  bool is_empty() const;

  //! returns the length of the string
  size_t length() const;

  //! returns a const char* for the string
  const char *c_str() const;
  
  //! returns a std::string
  const std::string& str() const;

  //! clear the string
  void clear();

  //! resize the string.  If growing, it is filled with null characters.
  void resize(size_t n);
  //! resize the string.  If growing, it is filled with the given character.
  void resize(size_t n, char c);

  // UTF-16/UTF-8 unicode conversions to support using Unicode aware Windows APIs.
  // Note: wchar_t on non-Windows platforms is UTF-32
#if defined(_WIN32)
  static std::wstring toUtf16(const char* str);
  static std::wstring toUtf16(const CubitString& str);
  static CubitString toUtf8(const wchar_t* str);
#endif

#if defined(_MSC_VER)
  static std::wstring toNative(const char* str);
  static std::wstring toNative(const CubitString& str);
#else
  static CubitString toNative(const char* str);
  static CubitString toNative(const CubitString& str);
#endif

  static std::wstring toWide(const CubitString& str);
  static CubitString toNarrow(const std::wstring& str);
  
private:
  std::string rep;
};

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

#endif

