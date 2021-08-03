/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#include "GetLongOpt.hpp"
#include "CubitMessage.hpp"
#include "CubitUtil.hpp"
#include "CubitFileUtil.hpp"

const CubitString GetLongOpt::TOGGLE_ON = "True";
const CubitString GetLongOpt::TOGGLE_OFF = "False";

GetLongOpt::GetLongOpt(const char optmark)
{
   optmarker = optmark;
}

GetLongOpt::~GetLongOpt()
{
}

int GetLongOpt::enroll(const CubitString &opt, const OptType t,
const CubitString &desc, const CubitString &val)
{
   Cell c;
   c.option = opt;
   c.type = t;
   c.description = !desc.is_empty() ? desc : "no description available";
   c.value = val;
   c.wasSet = 0;

   this->mTable.insert(std::make_pair(c.option, c));

   return 1;
}

CubitString
GetLongOpt::retrieve(const CubitString &opt) const
{
  std::map<CubitString, Cell>::const_iterator iter;
  iter = mTable.find(CubitString(opt));
  if(iter == mTable.end())
  {
    PRINT_ERROR("GetLongOpt::retrieve - unenrolled option %c%s\n", optmarker, opt.c_str() );
    return 0;
  }
  return iter->second.value.c_str();
}

int GetLongOpt::parse(const std::vector<CubitString>& args, const CubitString &p)
{
  int this_optind = 1;
  int argc = args.size();
  
  CubitString pname = p;
  if(pname.is_empty())
  {
    CubitString dirpart;
    CubitFileUtil::split_path(args[0], dirpart, pname);
  }

  if ( argc-- <= 1 )
  {
    return this_optind;
  }

  size_t i = 1;
  
  while ( argc >= 1 )
  {
    CubitString token = args[i++].c_str();
    CubitString next_token = i<args.size() ? args[i] : "";
    --argc;
    
    if ( token.get_at(0) != optmarker || token.get_at(1) == optmarker )
      break;	/* end of options */
    
    ++this_optind;

    size_t token_index = token.find("=", 1);
    CubitString token_value = token.substr(token_index != CubitString::npos ? token_index : CubitString::npos);
    token = token.substr(1, token_index);
    size_t token_len = token.length();

    Table::iterator t;
    enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
    Table::iterator pc = mTable.end();	// pointer to the partially-matched cell
    bool toggle = true; // toggle state for partial match cell

    for ( t = mTable.begin(); t != mTable.end(); t++ )
    {
      bool match = strncmp(token.c_str(), t->second.option.c_str(), token_len) == 0;
      bool no = false;

      if(!match && t->second.type == Toggle &&
         token.length() > 2 && token.get_at(0) == 'n' && token.get_at(1) == 'o')
      {
        match = strncmp(token.c_str() + 2, t->second.option.c_str(), token_len-2) == 0;
        no = true;
      }

      if (match)
      {
        if ( t->second.option.length() == (no ? token_len-2 : token_len) )
        {
            /* an exact match found */
          if (t->second.type == Toggle && token_value.is_empty())
            token_value = no ? GetLongOpt::TOGGLE_OFF : GetLongOpt::TOGGLE_ON;
            
          int stat = setcell(t->second, token_value, next_token, pname);
          if ( stat == -1 ) return -1;
          else if ( stat == 1 )
          {
            i++; --argc; ++this_optind;
          }
          matchStatus = ExactMatch;
          break;
        }
        else
        {
            /* partial match found */
          matchStatus = PartialMatch;
          pc = t;
          toggle = !no;
        }
      } /* end if */
    } /* end for */
    
    if ( matchStatus == PartialMatch )
    {
      if (pc->second.type == Toggle && token_value.is_empty() )
        token_value = toggle ? GetLongOpt::TOGGLE_ON : GetLongOpt::TOGGLE_OFF;
        
      int stat = setcell(pc->second, token_value, next_token, pname);
      if ( stat == -1 ) return -1;
      else if ( stat == 1 ) 
      {
        ++i; --argc; ++this_optind;
      }
    }
    else if ( matchStatus == NoMatch )
    {
      PRINT_ERROR("%s: unrecognized option %c%s\n",
                  pname.c_str(), optmarker, token.c_str() );
      return -1;		/* no match */
    }
    
  } /* end while */
  
  return this_optind;
}

int GetLongOpt::parse(const CubitString& opt_string, const CubitString& p)
{
  std::vector<CubitString> args;
  opt_string.tokenize(' ', args);
  return parse(args, p);
}

/* ----------------------------------------------------------------
GetLongOpt::setcell returns
   -1	if there was an error
    0	if the nexttoken was not consumed
    1	if the nexttoken was consumed
------------------------------------------------------------------- */
int
GetLongOpt::setcell(Cell &c, const CubitString& valtoken, const CubitString& nexttoken,
                    const CubitString& name)
{
   switch ( c.type ) {
    case GetLongOpt::Toggle :
      if ( valtoken != GetLongOpt::TOGGLE_ON && valtoken != GetLongOpt::TOGGLE_OFF ) {
        PRINT_ERROR("%s: unsolicited value for flag %c[no]%s\n",
          name.c_str(), optmarker, c.option.c_str() );
        return -1;	/* unsolicited value specification */
      }
      c.value = valtoken;
      return 0;
      
    case GetLongOpt::Valueless :
      if ( valtoken.get_at(0) == '=' ) {
        PRINT_ERROR("%s: unsolicited value for flag %c%s\n",
          name.c_str(), optmarker,c.option.c_str() );
        return -1;	/* unsolicited value specification */
      }
      if (!c.wasSet)
      {
        c.value = GetLongOpt::TOGGLE_ON;
        c.wasSet = 1;
      }
      return 0;
    case GetLongOpt::OptionalValue :
      if ( valtoken.get_at(0) == '=' ) {
        c.value = valtoken.substr(1);
      }
      else {
        if ( !nexttoken.is_empty() && nexttoken.get_at(0) != optmarker ) {
          c.value = nexttoken;
          return 1;
    	  }
      }

      // explicit return here, just to make sure another if-case isn't 
      // put in which falls through to the next case (in the absence of
      // a break statement)
      return 0;
    case GetLongOpt::MandatoryValue :
      int return_val;
      if ( valtoken.get_at(0) == '=' ) {
        c.value = valtoken.substr(1);
        return_val = 0;
      }
      else {
        if ( !nexttoken.is_empty() && nexttoken.get_at(0) != optmarker ) {
          c.value = nexttoken;
          return_val = 1;
      }
      else {
        PRINT_ERROR("%s: mandatory value for %c%s\n",
          name.c_str(), optmarker, c.option.c_str() );
        return_val = -1;	/* mandatory value not specified */
        }
      }
      return return_val;
    default :
      break;
   }
   return -1;
}

void 
GetLongOpt::options(std::ostream &outfile) const
{
  Table::const_iterator iter;
  for(iter = mTable.begin(); iter != mTable.end(); ++iter)
  {
    const Cell& t = iter->second;
    outfile << "\t" << optmarker;
    if ( t.type == GetLongOpt::Toggle )
      outfile << "[no]";
    outfile << t.option.str();
    if ( t.type == GetLongOpt::MandatoryValue )
      outfile << " <$val>";
    else if ( t.type == GetLongOpt::OptionalValue )
      outfile << " [$val]";
    outfile << " (" << t.description.str() << ")\n";
  }
}

void
GetLongOpt::usage( std::ostream &outfile, const CubitString& p, const CubitString& ustring) const
{
   outfile << "usage: " << p.str() << " " << ustring.str() << "\n";
   options(outfile);
}

