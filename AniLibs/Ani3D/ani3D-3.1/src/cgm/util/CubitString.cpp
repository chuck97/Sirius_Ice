// CubitString.cpp


#include "CubitString.hpp"
#include "CubitDefines.h"

#ifdef _WIN32
#include <windows.h>
#define vsnprintf _vsnprintf
#endif

#ifdef _MSC_VER
#pragma warning ( 4 : 4291 4244 4305 4018 )
#define va_copy(dest, src) (dest = src)
#endif


#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <stdarg.h>

const size_t CubitString::npos = std::string::npos;

CubitString::CubitString()
{
}

CubitString::CubitString(const CubitString& s)
  : rep(s.rep)
{
}

CubitString::~CubitString()
{
}

CubitString::CubitString(const char *s)
{
  if(s)
  {
    rep = s;
  }
}

CubitString::CubitString(const std::string& s)
  : rep(s)
{
}

CubitString::CubitString(int i, char c)
  : rep(i,c)
{
}

CubitString CubitString::number(double f, unsigned int max_length, unsigned int sig_digits)
{
  std::stringstream str;
  if(max_length)
  {
    str.width(max_length);
  }
  if(sig_digits)
  {
    str.precision(sig_digits);
  }
  str << f;
  std::string ret = str.str();
  size_t i = ret.find_first_not_of(' ');
  if(i != std::string::npos)
    ret = ret.substr(i);

  // change precision to be short enough
  if (max_length) 
  {
    if(sig_digits == 0)
      sig_digits = max_length;

    while(ret.length() > max_length && sig_digits)
    {
      sig_digits--;
      str.precision(sig_digits);
      str.str(std::string());
      str << f;
      ret = str.str();
    }
  }
  return ret.c_str();
}

CubitString CubitString::format(const char* format, ...)
{
  va_list args, args2;
  va_start(args, format);
  // copy because first vsnprintf modifies args
  va_copy(args2, args);

  // how much room do we need for our string?
  int num = vsnprintf(NULL, 0, format, args);

  // string plus null terminator
  std::vector<char> str;
  str.resize(num+1);

  // print string
  num = vsnprintf(&str[0], num+1, format, args2);

  // remove extra null terminator
  str.resize(num);

  va_end(args);
  va_end(args2);

  return str.empty() ? CubitString() : CubitString(&str[0]);
}

CubitString& CubitString::operator=(const CubitString& s)
{
  rep = s.rep;
  return *this;
}

#ifdef _WIN32
# define strcasecmp _stricmp
#endif

bool CubitString::compare(const CubitString& s, CaseSensitivity cs)
{
  if(cs == CaseSensitive)
    return this->rep == s.rep;

  // someday we might need to fix this to be unicode aware
  if(strcasecmp(this->rep.c_str(), s.rep.c_str()) == 0)
    return true;
  return false;
}

char CubitString::get_at(size_t pos) const
{
  // some legacy code wants to look for the null terminator
  if(pos == rep.size())
    return 0;

  return rep.at(pos);
}

bool CubitString::is_empty() const
{
  return rep.empty();
}

size_t CubitString::length() const
{
  return rep.length();
}

const char * CubitString::c_str() const
{
  return rep.c_str();
}

const std::string& CubitString::str() const
{
  return rep;
}

void CubitString::clear()
{
  rep.clear();
}

void CubitString::resize(size_t n)
{
  rep.resize(n);
}

void CubitString::resize(size_t n, char c)
{
  rep.resize(n, c);
}


CubitString operator+(const CubitString& s1, const CubitString& s2)
{
  return CubitString(s1) += s2;
}

CubitString operator+(const CubitString& s1, const char *c2)
{
  return CubitString(s1) += c2;
}

CubitString operator+(const char *c1, const CubitString& s2)
{
  return CubitString(c1) += s2;
}

CubitString& CubitString::operator+=(const CubitString& s)
{
  rep += s.rep;
  return *this;
}

size_t CubitString::find(const CubitString& s, size_t pos) const
{
  if(pos != 0 && pos >= rep.size())
    throw std::out_of_range("find index out of range");
  return rep.find(s.rep, pos);
}

size_t CubitString::find_first_of(const CubitString& s, size_t pos) const
{
  return rep.find_first_of(s.rep, pos);
}

size_t CubitString::find_first(char c, size_t pos) const
{
  return rep.find(c, pos);
}

size_t CubitString::find_last(char c, size_t pos) const
{
  return rep.find_last_of(c, pos);
}


CubitString CubitString::substr(size_t first, size_t count) const
{
  if(first >= rep.size())
    return CubitString();

  if(count > rep.size() - first)
    count = rep.size() - first;

  return rep.substr(first, count).c_str();
}

void CubitString::erase(size_t pos, size_t len)
{
  rep.erase(pos, len);
}

void CubitString::replace(size_t pos, size_t len, const CubitString& str)
{
  rep.replace(pos, len, str.rep);
}

void CubitString::replace(const CubitString& to_find, const CubitString& to_replace)
{
  size_t len = to_find.length();
  size_t len2 = to_replace.length();
  for(size_t idx = rep.find(to_find.c_str()); idx != npos; idx = rep.find(to_find.c_str(), idx))
  {
    rep.replace(idx, len, to_replace.rep);
    idx += len2;
  }
}

bool CubitString::ends_with(const CubitString& str) const
{
  if(rep.size() >= str.rep.size())
  {
    size_t offset = rep.size() - str.rep.size();
    return std::equal(str.rep.begin(), str.rep.end(), rep.begin() + offset);
  }
  return false;
}

bool CubitString::starts_with(const CubitString& str) const
{
  if(rep.size() >= str.rep.size())
  {
    return std::equal(str.rep.begin(), str.rep.end(), rep.begin());
  }
  return false;
}

void CubitString::put_at(size_t pos, char c)
{
  rep.at(pos) = c;
}

void CubitString::to_lower() 
{
  for(size_t i=0; i<rep.size(); i++)
  {
    rep.at(i) = tolower(rep.at(i));
  }
}

void CubitString::to_upper() 
{
  for(size_t i=0; i<rep.size(); i++)
  {
    rep.at(i) = toupper(rep.at(i));
  }
}

void CubitString::to_lower(char *string)
{
    // convert this string to lower case
  char *p = string;
  while (*p != '\0')
  {
    *p = tolower (*p);
    p++;
  }
}

void CubitString::to_upper(char *string)
{
    // convert this string to upper case
  char *p = string;
  while (*p != '\0')
  {
    *p = toupper (*p);
    p++;
  }
}

void CubitString::trim()
{
  std::string::size_type start = rep.find_first_not_of(" \t\n\v\f\r");
  rep.erase(0, start);
  std::string::size_type end = rep.find_last_not_of(" \t\n\v\f\r");
  rep.erase(end+1, std::string::npos);
}

void CubitString::simplify()
{
  trim();
  for(std::string::size_type i = rep.find_first_of(" \t\n\v\f\r"); i != std::string::npos;
      i = rep.find_first_of(" \t\n\v\f\r", i+1))
  {
    std::string::size_type next_i = rep.find_first_not_of(" \t\n\v\f\r", i);
    if(next_i > i+1)
    {
      rep.erase(i+1, next_i-i-1);
      rep.at(i) = ' ';
    }
  }
}

void CubitString::tokenize( char delimiter, std::vector<CubitString> &strings ) const
{
  std::stringstream ss(rep);
  std::string s;
  while(std::getline(ss, s, delimiter))
    strings.push_back(CubitString(s.c_str()));
}

bool CubitString::operator==(const CubitString& s2) const
{
  return rep == s2.rep;
}

bool CubitString::operator!=(const CubitString& s2) const
{
  return rep != s2.rep;
}

bool operator<=(const CubitString& s1, const CubitString& s2)
{
  return s1.rep <= s2.rep;
}

bool operator>=(const CubitString& s1, const CubitString& s2)
{
  return s1.rep >= s2.rep;
}

bool operator<( const CubitString& s1, const CubitString& s2 )
{
  return s1.rep < s2.rep;
}

bool operator>( const CubitString& s1, const CubitString& s2 )
{
  return s1.rep > s2.rep;
}

#if defined(_WIN32)

std::wstring CubitString::toUtf16(const char* str)
{
  // there are other methods (e.g. C++11), but since its Windows specific to support Windows APIs, just use their API.
  std::wstring ret;
  int sz = MultiByteToWideChar(CP_UTF8, 0, str, -1, NULL, 0);
  if(sz)
  {
    std::vector<wchar_t> result(sz);
    MultiByteToWideChar(CP_UTF8, 0, str, -1, &result[0], sz);
    ret = &result[0];
  }
  return ret;
}

std::wstring CubitString::toUtf16(const CubitString& str)
{
  return toUtf16(str.c_str());
}

CubitString CubitString::toUtf8(const wchar_t* str)
{
  CubitString ret;
  int sz = WideCharToMultiByte(CP_UTF8, 0, str, -1, NULL, 0, NULL, NULL);
  if(sz)
  {
    std::vector<char> result(sz);
    WideCharToMultiByte(CP_UTF8, 0, str, -1, &result[0], sz, NULL, NULL);
    ret = &result[0];
  }
  return ret;
}

#endif

#if defined(_MSC_VER)

std::wstring CubitString::toNative(const char* str)
{
  return toUtf16(str);
}

std::wstring CubitString::toNative(const CubitString& str)
{
  return toUtf16(str);
}

#else

CubitString CubitString::toNative(const char* str)
{
  return str;
}

CubitString CubitString::toNative(const CubitString& str)
{
  return str;
}

#endif

std::wstring CubitString::toWide(const CubitString& str)
{
#ifdef _WIN32
  return toUtf16(str);
#else
  std::wstring wstr;
  size_t len = mbstowcs(NULL, str.c_str(), 0) + 1;
  if(len > 0)
  {
    std::vector<wchar_t> wchar(len);
    mbstowcs(&wchar[0], str.c_str(), len);
    wstr = &wchar[0];
  }
  return wstr;
#endif
}

CubitString CubitString::toNarrow(const std::wstring& wstr)
{
#ifdef _WIN32
  return toUtf8(wstr.c_str()).c_str();
#else
  CubitString str;
  size_t len = wcstombs(NULL, wstr.c_str(), 0) + 1;
  if(len > 0)
  {
    std::vector<char> chars(len);
    wcstombs(&chars[0], wstr.c_str(), len);
    str = &chars[0];
  }
  return str;
#endif
}



#ifdef TEST_STRING

void main() {
CubitString a = "Test ";
CubitString b = "String Class";
CubitString blank(' ');

CubitString c = a + b;
c+= b;
CubitString e = c;
CubitString f = c;
CubitString g = c;
g.put_at(1, 'Z');
cout << "G = " << g << endl;
cout << "C = " << c << endl;

c.put_at(1, 'X');
CubitString d = b;
d.put_at(1, 'Y');
}
#endif
