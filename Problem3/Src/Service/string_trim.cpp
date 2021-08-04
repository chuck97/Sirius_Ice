#include "string_trim.h"

using namespace std;

void ltrim(string &s) 
{
    s.erase(s.begin(), find_if(s.begin(), s.end(), 
    [](unsigned char ch) 
    {
        return !isspace(ch);
    }));
};

void rtrim(string &s) 
{
    s.erase(find_if(s.rbegin(), s.rend(), 
    [](unsigned char ch) 
    {
        return !isspace(ch);
    }).base(), s.end());
};

void trim(string &s) 
{
    ltrim(s);
    rtrim(s);
};