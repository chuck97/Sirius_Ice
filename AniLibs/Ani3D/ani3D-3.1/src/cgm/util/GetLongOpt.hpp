//- Class: GetLongOpt
//- Description: GetLongOpt manages the definition and parsing of 
//- long options. Long options can be abbreviated as long as there 
//- is no ambiguity. If an option requires a value, the value should
//- be separated from the option either by whitespace or an "=".
//- 
//- Other features: $
//- o	GetLongOpt can be used to parse options given through environments.$
//- o	GetLongOpt provides a usage function to print usage.$
//- o	Flags & options with optional or mandatory values are supported.$
//- o	The option marker ('-' in Unix) can be customized.$
//- o	Parsing of command line returns optind (see getopt(3)).$
//- o	Descriptive error messages.$
//-
//- Author: S Manoharan. Advanced Computer Research Institute. Lyon. France
//- Owner: Greg Sjaardema
//- Checked By:
//- Version $Id: 

#ifndef GETLONGOPT_HPP
#define GETLONGOPT_HPP

#include "CubitDefines.h"

#include <iostream>

#include <cstring>
#include "CubitUtilConfigure.h"
#include <vector>
#include <map>
#include "CubitString.hpp"

class CUBIT_UTIL_EXPORT GetLongOpt {
public:
  enum OptType { 
     Toggle, Valueless, OptionalValue, MandatoryValue
   };

  static const CubitString TOGGLE_ON;
  static const CubitString TOGGLE_OFF;

  GetLongOpt(const char optmark = '-'); //- Constructor
  //- Constructor for GetLongOpt takes an optional argument: the option
  //- marker. If unspecified, this defaults to '-', the standard (?)
  //- Unix option marker. 

  ~GetLongOpt();                        //- Destructor
  
  int parse(const std::vector<CubitString>& args, const CubitString& p = CubitString());
  //- A return
  //- value < 1 represents a parse error. Appropriate error messages
  //- are printed when errors are seen. {GetLongOpt::parse}, in its first
  //- form, takes two strings: the first one is the string to be
  //- parsed and the second one is a string to be prefixed to the
  //- parse errors. In ts second form, {GetLongOpt::parse} returns the
  //- the {optind} (see @i{getopt(3)}) if parsing is successful.

  int parse(const CubitString& opt_string, const CubitString &p);

  int enroll(const CubitString& opt, const OptType t,
       const CubitString& desc, const CubitString& val);
  //- Add an option to the list of valid command options.$
  //- {GetLongOpt::enroll} adds option specifications to its internal
  //- database. The first argument is the option sting. The second
  //- is an enum saying if the option is a flag ({GetLongOpt::Valueless}),
  //- if it requires a mandatory value ({GetLongOpt::MandatoryValue}) or
  //- if it takes an optional value ({GetLongOpt::OptionalValue}).
  //- The third argument is a string giving a brief description of
  //- the option. This description will be used by {GetLongOpt::usage}.
  //- GetLongOpt, for usage-printing, uses {$val} to represent values
  //- needed by the options. {<$val>} is a mandatory value and {[$val]}
  //- is an optional value. The final argument to {GetLongOpt::enroll}
  //- is the default string to be returned if the option is not
  //- specified. For flags (options with {Valueless}), use "" (empty
  //- string, or in fact any arbitrary string) for specifying {TRUE}
  //- and 0 (null pointer) to specify {FALSE}.
  
  CubitString retrieve(const CubitString& opt) const;
  //- Retrieve value of option$
  //- The values of the options that are enrolled in the database
  //- can be retrieved using {GetLongOpt::retrieve}. This returns a string
  //- and this string should be converted to whatever type you want.
  //- See @i{atoi}, @i{atof}, @i{atol} etc.
  //- If a "parse" is not done before retrieving all you will get
  //- are the default values you gave while enrolling!
  //- Ambiguities while retrieving (may happen when options are
  //- abbreviated) are resolved by taking the matching option that 
  //- was enrolled last. For example, -{v} will expand to {-verify}.$
  //- If you try to retrieve something you didn't enroll, you will
  //- get a warning message. 

  void options(std::ostream &outfile = std::cout) const;
  //- Print command line options only. Called by usage().

  void usage(std::ostream &outfile, const CubitString &p, const CubitString& u = CubitString()) const;
  //- Print usage information to {outfile}

private:
  struct Cell {
     CubitString option;	//- option name
     OptType type;		//- option type
     CubitString description;	//- a description of option
     CubitString value;	        //- value of option (string)
     int wasSet;                //- 0 if not set by user, 1 if set
     
     Cell() { option = description = value = 0; wasSet = 0; }
   };

  typedef std::map<CubitString, Cell> Table;
  Table mTable;

  char optmarker;		//- option marker
  
  int setcell(Cell &c, const CubitString &valtoken, const CubitString &nexttoken, const CubitString &name);
};

#endif // GETLONGOPT_HPP

