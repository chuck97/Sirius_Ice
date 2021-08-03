
#include <fstream>
#include <iomanip>

#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
#include <cassert>
#include <cstring>
#include <vector>
#include <iostream>

#include "SettingHandler.hpp"

#ifdef _WIN32
#define vsnprintf _vsnprintf
//#define strdup _strdup
#endif

int CubitMessage::errorCount = 0;
int CubitMessage::warningCount = 0;
//
// Message type usage:
// PRINT_ERROR:    Message to tell user why the task did not complete.
// PRINT_WARNING:  Message to tell user why completed task may not be what
//                 was requested.
// PRINT_INFO:     Message to tell user about status and progress.
// PRINT_DEBUG:    Message to developer tied to a global debug flag.
// DIAGNOSTIC:     Message to developer.


struct DefaultCubitMessageHandler : public CubitMessageHandler
{
  void print_message(const char* message)
  {
    std::cout << message;
    std::cout.flush();
  }
  void print_error(const char* message)
  {
#ifdef XTERM
    char esc = 0x1B;
    // Turn on reverse video on VT102 (xterm also)
    // (0=normal, 1-bold, 4-underscore, 5-blink, 7-inverse)
    std::cout << esc << '[' << '7' << 'm';
#endif

    std::cout << message;
    std::cout.flush();

#ifdef XTERM
    std::cout << esc << '[' << '0' << 'm';
    std::cout.flush();
#endif

  }
};

static DefaultCubitMessageHandler mDefaultHandler;


CubitMessage* CubitMessage::instance_ = NULL;
CubitMessageHandler* CubitMessage::mHandler = &mDefaultHandler;
CubitMessageErrorHandler* CubitMessage::mErrorHandler = NULL;
int CubitMessage::infoFlag = CUBIT_TRUE;
int CubitMessage::diagnosticFlag = CUBIT_FALSE;
int CubitMessage::warningFlag = CUBIT_TRUE;
int CubitMessage::errorFlag = CUBIT_TRUE;
std::ofstream* CubitMessage::loggingStream = NULL;
CubitString* CubitMessage::loggingFile = NULL;
std::ofstream* CubitMessage::loggingErrorStream = NULL;
CubitString* CubitMessage::loggingErrorFile = NULL;
int CubitMessage::expectedStartErrorCount = -1;
int CubitMessage::expectedEndErrorCount = -1;
bool CubitMessage::expectedLessErrorCountAccepted = false;

CubitMessage* CubitMessage::instance()
{
  if (!instance_)
  {
    instance_ = new CubitMessage;
    if (!instance_)
    {
	  std::cerr << " *** Unable to instantiate message object ***" << std::endl;
      exit(1);
    }
  }
  return instance_;
}

void CubitMessage::free_instance()
{
  if (instance_)
    delete instance_;
  instance_ = NULL;
}

CubitMessage::CubitMessage()
{
  infoFlag      = CUBIT_TRUE;
  warningFlag   = CUBIT_TRUE;
  errorFlag     = CUBIT_TRUE;
  diagnosticFlag= CUBIT_FALSE;
  loggingStream = NULL;
  loggingFile   = NULL;
  loggingErrorStream = NULL;
  loggingErrorFile   = NULL;
  currentDebugFlag = CUBIT_DEBUG_1;

    // Initialize the debugFlag array
  static MessageFlag staticDebugFlag[] =
  {
    MessageFlag(  0, "UNUSED"),
    MessageFlag(  1, "Previously used; available for reuse."),
    MessageFlag(  2, "Whisker weaving information"),
    MessageFlag(  3, "Timing information for 3D Meshing routines."),
    MessageFlag(  4, "Graphics Debugging (DrawingTool)"),
    MessageFlag(  5, "FastQ debugging"),
    MessageFlag(  6, "Submapping graphics debugging"),
    MessageFlag(  7, "Knife progress whisker weaving information"),
    MessageFlag(  8, "Mapping Face debug / Linear Programing debug "),
    MessageFlag(  9, "Paver Debugging"),
    MessageFlag( 10, "WW: removed hex seam flag"),
    MessageFlag( 11, "Nodeset Associativity debugging"),
    MessageFlag( 12, "Fastq activity"),
    MessageFlag( 13, "Mesh entities"),
    MessageFlag( 14, "Previously used; available for reuse."),
    MessageFlag( 15, "Previously used; available for reuse."),
    MessageFlag( 16, "Previously used; available for reuse."),
    MessageFlag( 17, "Use Count debugging"),
    MessageFlag( 18, "Webcut debugging"),
    MessageFlag( 19, "Feature Merge / Unmerge debugging"),
    MessageFlag( 20, "Parallel meshing activity"),
    MessageFlag( 21, "Boundary Layer Tool Debugging"),
    MessageFlag( 22, "ExodusMesh sizing function debugging"),
    MessageFlag( 23, "Draw after joining chords in WW"),
    MessageFlag( 24, "SelfCrossingLoop (and derivatives) debug info"),
    MessageFlag( 25, "Extra invalidity checking in WW"),
    MessageFlag( 26, "Surface Smoothing"),
    MessageFlag( 27, "Primal Construction debugging, see also flag 70"),
    MessageFlag( 28, "Plastering debugging"),
    MessageFlag( 29, "Volume SubMapping"),
    MessageFlag( 30, "Previously used; available for reuse."),
    MessageFlag( 31, "CleanUp debugging"),
    MessageFlag( 32, "Previously used; available for reuse."),
    MessageFlag( 33, "Whisker Weaving inside chord list face drawing"),
    MessageFlag( 34, "If on Whisker Weaving doesn't merge sheets"),
    MessageFlag( 35, "If on WW query displays sheets before joining chords"),
    MessageFlag( 36, "Enable/Disable idr_keyword_debugger function"),
    MessageFlag( 37, "Previously used; available for reuse."),
    MessageFlag( 38, "WW hex formation messages"),
    MessageFlag( 39, "Doublet Pillower graphics output"),
    MessageFlag( 40, "Previously used; available for reuse."),
    MessageFlag( 41, "Previously used; available for reuse."),
    MessageFlag( 42, "Auto vertex type and sweep verification"),
    MessageFlag( 43, "Programmer Errors for SubMapping"),
    MessageFlag( 44, "Submapping Graphics Debugging"),
    MessageFlag( 45, "Pillow Sheet debugging"),
    MessageFlag( 46, "Paver breakout detection (expensive)"),
    MessageFlag( 47, "Extra LP debugging (see flag 8 also)"),
    MessageFlag( 48, "Previously used; available for reuse."),
    MessageFlag( 49, "Draws Face by Face Creation in Paving"),
    MessageFlag( 50, "Debugging for AutoSchemeSelect"),
    MessageFlag( 51, "Previously used; available for reuse."),
    MessageFlag( 52, "User Interface: If flag is enabled, filenames being\n"
                "\t\t\t\tused for input will be echoed and each input\n"
                "\t\t\t\tline will be echoed prior to being parsed."),
    MessageFlag( 53, "Surface Morpher debugging"),
    MessageFlag( 54, "Parser debugging"),
    MessageFlag( 55, "Previously used; available for reuse."),
    MessageFlag( 56, "Previously used; available for reuse."),
    MessageFlag( 57, "Relative Interval/Length setting"),
    MessageFlag( 58, "StcVertex debugging of Whisker Weaving"),
    MessageFlag( 59, "Previously used; available for reuse."),
    MessageFlag( 60, "StcVertex debugging of Looping"),
    MessageFlag( 61, "List number of points used in curve faceting"),
    MessageFlag( 62, "Print verbose information on group operations"),
    MessageFlag( 63, "Label Whisker Weaving diagrams tersely"),
    MessageFlag( 64, "No label on Whisker Weaving diagrams"),
    MessageFlag( 65, "Volume Morpher debugging"),
    MessageFlag( 66, "Print debug information on importing Pro/E geometry"),
    MessageFlag( 67, "List number of triangles used in surface faceting"),
    MessageFlag( 68, "Previously used; available for reuse."),
    MessageFlag( 69, "Previously used; available for reuse."),
    MessageFlag( 70, "STC Pillowing, see also flag 27"),
    MessageFlag( 71, "Previously used; available for reuse."),
    MessageFlag( 72, "DoubletPillower text messages"),
    MessageFlag( 73, "Auto Surface debugging (use new auto surf select)"),
    MessageFlag( 74, "Feature-based decomposition info"),
    MessageFlag( 75, "Many-to-many sweep imprint debugging"),
    MessageFlag( 76, "Virtual point and partition curve"),
    MessageFlag( 77, "Volume interval matching"),
    MessageFlag( 78, "Tipton Smoother jacobian modification enabler"),
    MessageFlag( 79, "Previously used; available for reuse."),
    MessageFlag( 80, "Previously used; available for reuse."),
    MessageFlag( 81, "Curve Morpher Debugging"),
    MessageFlag( 82, "Previously used; available for reuse."),
    MessageFlag( 83, "Previously used; available for reuse."),
    MessageFlag( 84, "Surface auto decomposition"),
    MessageFlag( 85, "U-SubMapping debugging"),
    MessageFlag( 86, "Virtual curve and partition surface"),
    MessageFlag( 87, "Composite curve and composite surface"),
    MessageFlag( 88, "Volume partitioning"),
    MessageFlag( 89, "Previously used; available for reuse."),
    MessageFlag( 90, "Geometry attributes"),
    MessageFlag( 91, "Smoothing Debug Output"),
    MessageFlag( 92, "Print name changed warnings"),
    MessageFlag( 93, "Hex Fix Up"),
    MessageFlag( 94, "Entity name attribute"),
    MessageFlag( 95, "Group imprint errors"),
    MessageFlag( 96, "GraftTool debugging"),
    MessageFlag( 97, "Previously used; available for reuse."),
    MessageFlag( 98, "Color code imported THEX meshes"),
    MessageFlag( 99, "Geometry creation"),
    MessageFlag(100, "Skew Control debugging"),
    MessageFlag(101, "Previously used; available for reuse."),
    MessageFlag(102, "CAEntityId debugging"),
    MessageFlag(103, "Print compact interval assignment constraints"),
    MessageFlag(104, "Report interval matching progress"),
    MessageFlag(105, "Previously used; available for reuse."),
    MessageFlag(106, "Mesh Cleaver debugging"),
    MessageFlag(107, "Midpoint_subdivision debugging"),
    MessageFlag(108, "Simulog tetmesher debugging"),
    MessageFlag(109, "Transition schemes debugging"),
    MessageFlag(110, "Mesh Defined Geometry"),
    MessageFlag(111, "TriAdvance mesher debugging"),
    MessageFlag(112, "Auto Detail Suppression"),
    MessageFlag(113, "Previously used; available for reuse."),
    MessageFlag(114, "Blend Finder Debugging"),
    MessageFlag(115, "Exporting Feature Debugging Files"),
    MessageFlag(116, "Sizing function tool data information"),
    MessageFlag(117, "Extra Information on Autoscheme Decision Making"),
    MessageFlag(118, "Blend finding optimization file"),
    MessageFlag(119, "Laminate Tool debugging"),
    MessageFlag(120, "Print unassociated node locations on import mesh"),
    MessageFlag(121, "Print verbose infeasible match interval messages"),
    MessageFlag(122, "Mesh-Based Geometry Debug Information"),
    MessageFlag(123, "Collect memory statistics from Tetmesher"),
    MessageFlag(124, "Print verbose Tetmesher debugging information"),
    MessageFlag(125, "Mesh refinement debugging"),
    MessageFlag(126, "Previously used; available for reuse."),
    MessageFlag(127, "SculptingTool debug flag"),
    MessageFlag(128, "Previously used; available for reuse."),
    MessageFlag(129, "Virtual Imprint Debugging"),
    MessageFlag(130, "Hexsheet Insertion Debugging"),
    MessageFlag(131, "Mesh Cutting Debugging"),
    MessageFlag(132, "Global Collection Smoothing"),
    MessageFlag(133, "Print verbose import mesh progress"),
    MessageFlag(134, "Previously used; available for reuse."),
    MessageFlag(135, "Keep WhiskerWeave data after meshing"),
    MessageFlag(136, "Previously used; available for reuse."),
    MessageFlag(137, "GJoin"),
    MessageFlag(138, "Parallel CGM timing"),
    MessageFlag(139, "RTree Debugging"),
    MessageFlag(140, "Previously used; available for reuse."),
    MessageFlag(141, "Settings save/restore"),
    MessageFlag(142, "Decompose Sweep Debugging"),
    MessageFlag(143, "Decomp Sweep Imprint Debugging"),
    MessageFlag(144, "Medial Axis/Chordal Axis Debugging"),
    MessageFlag(145, "Virtual Geometry Facet Operations"),
    MessageFlag(146, "Sector Tool Meshing Scheme"),
    MessageFlag(147, "Previously used; available for reuse."),
    MessageFlag(148, "Meshing Benchmarks"),
    MessageFlag(149, "MeshCutting Graphical debugging"),
    MessageFlag(150, "MBG to Acis conversion debugging"),
    MessageFlag(151, "Previously used; available for reuse."),
    MessageFlag(152, "Boundary Conditions Debugging"),
    MessageFlag(153, "Print Body information in Geometry operations"),
    MessageFlag(154, "Split Surface Debugging"),
    MessageFlag(155, "Meshing Benchmarks Summary"),
    MessageFlag(156, "CAMAL Paver CleanUp debuging"),
    MessageFlag(157, "Skeleton Sizing Function timing and counts"),
    MessageFlag(158, "Previously used; available for reuse."),
    MessageFlag(159, "Previously used; available for reuse."),
    MessageFlag(160, "Previously used; available for reuse."),
    MessageFlag(161, "Previously used; available for reuse."),
    MessageFlag(162, "Previously used; available for reuse."),
    MessageFlag(163, "Previously used; available for reuse."),
    MessageFlag(164, "Previously used; available for reuse."),
    MessageFlag(165, "Previously used; available for reuse."),
    MessageFlag(166, "Unconstrained Paving debugging"),
    MessageFlag(167, "Skeleton Sizing Function Debugging (messages)"),
    MessageFlag(168, "Tweak Target Multiple Debugging "),
    MessageFlag(169, "Enable Knupp affine transformation instead of Roca"),
    MessageFlag(170, "Previously used; available for reuse."),
    MessageFlag(171, "Previously used; available for reuse."),
    MessageFlag(172, "Previously used; available for reuse."),
    MessageFlag(173, "Previously used; available for reuse."),
    MessageFlag(174, "Previously used; available for reuse."),
    MessageFlag(175, "Previously used; available for reuse."),
    MessageFlag(176, "Enable UCP database checking"),
    MessageFlag(177, "Enable Unconstrained Plastering Debug Drawing"),
    MessageFlag(178, "Enable Harris instead of Parrish hex refinement"),
    MessageFlag(179, "Enable Camal Sweeper for UCP Front Advancements.\n"
                "\t\t\t\tIgnored if debug 189 is on"),
    MessageFlag(180, "DecompAide (decomposition helper) debugging"),
    MessageFlag(181, "MBG.  Draw curve paths."),
    MessageFlag(182, "UCP Detailed Debug Printing."),
    MessageFlag(183, "Previously used; available for reuse."),
    MessageFlag(184, "Enable old sheet refinement command."),
    MessageFlag(185, "Enable straddle elements on hardlines."),
    MessageFlag(186, "Disable parametric coordinates in TriAdvMesher"),
    MessageFlag(187, "Previously used; available for reuse."),
    MessageFlag(188, "Tolerant Triangle Meshing"),
    MessageFlag(189, "Previously used; available for reuse."),
    MessageFlag(190, "Count CAMAL calls to move_to and normal_at"),
    MessageFlag(191, "Auto clean messages"),
    MessageFlag(192, "Previously used; available for reuse."),
    MessageFlag(193, "Tetmesh with stand-alone INRIA execuable thru files"),
    MessageFlag(194, "Hex Mesh Matching Debug drawing"),
    MessageFlag(195, "Previously used; available for reuse."),
    MessageFlag(196, "Previously used; available for reuse."),
    MessageFlag(197, "CAMAL Paver debugging (no Cubit smoothing, etc.)"),
    MessageFlag(198, "Auto Midsurface debugging"),
    MessageFlag(199, "Angle smoothing debugging"),
    MessageFlag(200, "Paver quality data output"),
    MessageFlag(201, "Paver cleanup edge metrics"),
    MessageFlag(202, "Disable Paver cleanup 3-3 replace"),
    MessageFlag(203, "Disable Paver cleanup 3-offset-3/5 replace"),
    MessageFlag(204, "Enable Paver cleanup 3-valent quad cluster"),
    MessageFlag(205, "Enable Paver cleanup partial chord collapse"),
    MessageFlag(206, "Hex Mesh Matching, match chords one at a time."),
    MessageFlag(207, "Defeature and Geometry tolerant meshing"),
    MessageFlag(208, "Previously used; available for reuse."),
    MessageFlag(209, "Change sense of partial/full tet remesh (v = !v)"),
    MessageFlag(210, "Use tetgen tetmesher via files"),
    MessageFlag(211, "Use tetgen tetmesher via direct interface"),
    MessageFlag(212, "Create debugging groups when doing geometry/meshing association for parallel refinement"),
    MessageFlag(213, "Previously used; available for reuse."),
    MessageFlag(214, "Boundary Layers"),
    MessageFlag(215, "Command Parser"),
    MessageFlag(216, "Command Parser Detailed"),
    MessageFlag(217, "Graphics debugging for unite dissimilar mesh command"),
    MessageFlag(218, "Previously used; available for reuse."),
    MessageFlag(219, "Materials Interface"),
    MessageFlag(220, "Previously used; available for reuse."),
    MessageFlag(221, "Turn ON Boundary Layer Correction in New Sweeper; ignored if target evaluator is NULL."),
    MessageFlag(222, "Turn ON Scale Mesh option to maximize percentage of model which gets scaled."),
    MessageFlag(223, "Turn ON Scale Mesh Debug Drawing."),
    MessageFlag(224, "Turn ON Scale Mesh Extended Debug Drawing."),
    MessageFlag(225, "unassigned")
    
      // IMPORTANT!!!
      // If you add a new debug flag, make sure that you change
      // the result of CubitMessage::number_of_debug_flags().
      // In order to use this type of static initialization,
      // we can't use the sizeof operator, so change it manually.
  };
  debugFlag = staticDebugFlag;

  // Check initialization of debugFlag array.
  for (int i=number_of_debug_flags(); i > 0; i--)
  {
    debugFlag[i].setting = CUBIT_FALSE;
    assert(i == debugFlag[i].flagNumber);
    assert(debugFlag[i].description != NULL);
  }
}

CubitMessage::~CubitMessage()
{
  // Close all streams associated with debug flags.
  // If the same stream is being used for debug and logging, it
  // will get closed below.
  for (int i=number_of_debug_flags(); i > 0; i--)
    remove_debug_stream(i);

  // At this time, the only open streams possible are loggingStream and loggingErrorStream.
  if (loggingStream != NULL)
  {
    loggingStream->close();
    delete loggingStream;
    delete loggingFile;
  }
  if (loggingErrorStream != NULL)
  {
    loggingErrorStream->close();
    delete loggingErrorStream;
    delete loggingErrorFile;
  }

  // Set static instance_ to zero to indicated that we are dead.
  instance_ = 0;
}

void CubitMessage::delete_instance()
{
  delete instance_;
  instance_ = NULL;
}

int CubitMessage::number_of_debug_flags()
{
  return NUM_DEBUG_FLAGS;
//  return sizeof(debugFlag)/sizeof(debugFlag[0])-1;
}

void CubitMessage::internal_error ( const int message_type,
                                    std::ofstream *output_stream,
                                    const CubitString& msgbuf)
{
  int print_it = CUBIT_FALSE;

  CubitString prefix;

  switch (message_type)
  {
    case CUBIT_ERROR:
      if (errorFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "ERROR: ";
      }
      break;
    case CUBIT_ERROR_EXPECTED:
      if (errorFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "ERROR_EXPECTED: ";
      }
      break;
    case CUBIT_WARNING:
      if (warningFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "WARNING: ";
      }
      break;
    case CUBIT_INFO:
      if (infoFlag)
        print_it = CUBIT_TRUE;
      break;
    case CUBIT_DIAGNOSTIC:
      if (diagnosticFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "DIAGNOSTIC: ";
      }
      break;
    default:
      if (message_type >= CUBIT_DEBUG_1 && message_type <= number_of_debug_flags()+10)
      {
        if (debugFlag[message_type-10].setting) print_it = CUBIT_TRUE;
        break;
      }
  }

  if (print_it)
  {
      // loggingStream is used to journal error, warning, and info messages.
      // debug messages can also be journalled there by setting the
      // output stream for the debug flag to the same file.
    if (loggingStream != NULL && (message_type == CUBIT_ERROR ||
                                  message_type == CUBIT_WARNING ||
                                  message_type == CUBIT_INFO))
    {
      *loggingStream << prefix.c_str() << msgbuf.c_str();
      loggingStream->flush();
    }
      //loggingErrorStream is used to (if the user has requested it)
      // log only ERROR: messages
    if (loggingErrorStream != NULL && message_type == CUBIT_ERROR)
    {
      *loggingErrorStream << prefix.c_str() << msgbuf.c_str();
      loggingErrorStream->flush();
    }

    if (output_stream == NULL)
    {
      if(message_type == CUBIT_ERROR)
      {
        CubitString ctx;
        if(CubitMessage::mErrorHandler)
        {
          ctx = CubitMessage::mErrorHandler->error_context().c_str();
        }
        CubitString msg = prefix + ctx + msgbuf;
        CubitMessage::mHandler->print_error(msg.c_str());
      }
      else
      {
        CubitString msg = prefix + msgbuf;
        CubitMessage::mHandler->print_message(msg.c_str());
      }
    }
    else
    {
      *output_stream << prefix.c_str() << msgbuf.c_str();
      output_stream->flush();
    }
  }
}

int CubitMessage::print_error ( const CubitString& str )
{
  int error_type = CUBIT_ERROR;
  if(expectedStartErrorCount != -1 && expectedEndErrorCount > expectedStartErrorCount)
  {
    error_type = CUBIT_ERROR_EXPECTED;
  }
  else
  {
    add_to_error_count();
  }

  if(expectedStartErrorCount != -1)
  {
    expectedStartErrorCount++;
  }

  static char* cubit_ctest = getenv("CUBIT_CTEST");
  if(cubit_ctest && error_type == CUBIT_ERROR)
  {
    internal_error(CUBIT_INFO, NULL, "<DartMeasurement name=\"Error\" type=\"text/plain\">\n");
    internal_error(error_type, NULL, str);
    internal_error(CUBIT_INFO, NULL, "</DartMeasurement>\n");
  }

  internal_error(error_type, NULL, str);

  return CUBIT_FAILURE;
}

int CubitMessage::print_warning ( const CubitString& str )
{
  internal_error(CUBIT_WARNING, NULL, str);
  add_to_warning_count();
  return CUBIT_FAILURE;
}

int CubitMessage::print_info ( const CubitString& str )
{
  internal_error(CUBIT_INFO, NULL, str);
  return CUBIT_FAILURE;
}

int CubitMessage::is_debug_flag_set( int flag )
{
   if( DEBUG_FLAG( flag ))
   {
      currentDebugFlag = flag;
      return CUBIT_TRUE;
   }
   return CUBIT_FALSE;
}

int CubitMessage::print_debug( const CubitString& str )
{
  internal_error(currentDebugFlag+10,
                 debugFlag[currentDebugFlag].outputStream,
                 str);
  return CUBIT_FAILURE;
}

void CubitMessage::print_diagnostic ( const CubitString& str )
{
  internal_error(CUBIT_DIAGNOSTIC, NULL, str);
}

int CubitMessage::reset_error_count(int value)
{
  int current_value = errorCount;
  if (errorCount != value) {
    errorCount = value;
    PRINT_WARNING("Error count manually changed from %d to %d\n\n",
		  current_value, value);
  }
  return current_value;
}

int CubitMessage::error_count()
{
  return errorCount;
}

void CubitMessage::add_to_error_count()
{
  errorCount++;
}

int CubitMessage::reset_warning_count(int value)
{
  int current_value = warningCount;
  if (warningCount != value) {
    warningCount = value;
    PRINT_INFO("Warning count manually changed from %d to %d\n\n",
		  current_value, value);
  }
  return current_value;
}

int CubitMessage::warning_count()
{
  return warningCount;
}

void CubitMessage::add_to_warning_count()
{
  warningCount++;
}

void CubitMessage::output_debug_information(int from, int to, int step)
{
  if (to == -1)
    to = number_of_debug_flags();

  PRINT_INFO("Debug Flag Settings "
	     "(flag number, setting, output to, description):\n");
   for (int i=from; i <= to; i+=step) {
      debugFlag[i].output();
   }
  PRINT_INFO("\n");
}

void CubitMessage::output_debug_information(CubitString &match)
{
  int count = 0;
  for (int i=1; i <= number_of_debug_flags(); i++) {
    char *tmp = CubitUtil::util_strdup((char*)(debugFlag[i].description));
    if (tmp && strlen(tmp) > 0) {
      CubitString debug_description(tmp);
      debug_description.to_lower();
      if (debug_description.find(match, 0) < debug_description.length()) {
	if (count == 0) {
	  PRINT_INFO("Debug Flag Settings "
		     "(flag number, setting, output to, description):\n");
	}
	debugFlag[i].output();
	count++;
      }
    }
    CubitUtil::util_strdup_free(tmp);
  }
  if (count == 0) {
    PRINT_WARNING("No debug descriptions contain the "
		  "substring '%s'\n", match.c_str());
  }
  PRINT_INFO("\n");
}

void CubitMessage::output_logging_information()
{
  if (loggingStream != NULL)
     PRINT_INFO("logging           = On, log file = '%s'\n", loggingFile->c_str());
  else
     PRINT_INFO("logging           = Off\n");
  if (loggingErrorStream != NULL)
     PRINT_INFO("logging Errors    = On, log file = '%s'\n",loggingErrorFile->c_str());

}

void MessageFlag::output()
{
  CubitMessage::instance()->
    print_info(
        CubitString::format("%2d  %3s  %-16s   %s\n",
               flagNumber, (setting == 1 ? "ON " : "OFF"),
               (filename == NULL ? "terminal" : filename->c_str()),
               description)
        );
}

int CubitMessage::find_file_use(const CubitString &filename)
{
  if (filename == "terminal") {
    // remove_debug_stream has set the outputStream and filename to NULL.
    return -1;
  }

  // See if any of the other debug flags have this file open
  for (int i=number_of_debug_flags(); i > 0; i--) {
    if (debugFlag[i].filename && *(debugFlag[i].filename) == filename) {
      return i;
    }
  }
  if (loggingFile && *(loggingFile) == filename)
    return -2;

  if (loggingErrorFile && *(loggingErrorFile) == filename)
    return -3;

  return 0;
}

int CubitMessage::count_stream_users(const std::ofstream *stream)
{
  int match = 0;
  if (stream != NULL)
  {
    for (int i=number_of_debug_flags(); i > 0; i--)
    {
      if (debugFlag[i].outputStream == stream)
      {
        match++;
      }
    }

    if (loggingStream == stream)
      match++;
    if (loggingErrorStream == stream)
       match++;
  }
  return match;
}

void CubitMessage::set_logging_file_setting(const CubitString &filename, CubitBoolean resume_flag)
{
  // If logging is currently outputting to a file, close it if
  // it is the only thing using that file. (and the filenames don't match)
  if (loggingFile && *loggingFile == filename)
    return;

  if (loggingErrorFile && *loggingErrorFile == filename)
  {
    PRINT_ERROR("Can't set the logging file to be the same as the Error logging file.\n");
    return;
  }

  int users = count_stream_users(loggingStream);
  if (users == 1) { // Just us...
    loggingStream->close();
    delete loggingStream;
    loggingStream = NULL;
    delete loggingFile;
    loggingFile = NULL;
  }

  int match = find_file_use(filename);

  if (match == -1) // Filename is 'terminal'
    return;
  else if (match != 0)
  {
    loggingFile   = debugFlag[match].filename;
    loggingStream = debugFlag[match].outputStream;
  }
  else
  {
    loggingFile   = new CubitString(filename);
    if(resume_flag)
       loggingStream = new std::ofstream(CubitString::toNative(filename).c_str(), std::ios::out | std::ios::app);
    else
       loggingStream = new std::ofstream(CubitString::toNative(filename).c_str());
  }
}

void CubitMessage::set_debug_file_setting(const int index, const CubitString &filename)
{
  // If this flag is currently outputting to a file, close it if
  // this is the only flag using that file.
  remove_debug_stream(index);

  int match = find_file_use(filename);

  if (match == -1) // Filename is 'terminal'
    return;
  if (match == -2 || match == -3) {// Filename is same as loggingFile or loggingErrorFile;
    debugFlag[index].filename = loggingFile;
    debugFlag[index].outputStream = loggingStream;
  }
  else if (match == index)
    return;
  else if (match != 0) {
    debugFlag[index].filename = debugFlag[match].filename;
    debugFlag[index].outputStream = debugFlag[match].outputStream;
  }
  else {
    debugFlag[index].filename = new CubitString(filename);
    debugFlag[index].outputStream = new std::ofstream(CubitString::toNative(filename).c_str());
  }
}

void CubitMessage::remove_debug_stream(const int index)
{
  // NOTE: DO NOT USE PRINT_* CALLS, THIS IS CALLED FROM DESTRUCTOR.

  // Multiple debug flags may be using the same output stream,
  // Go through the list and count who is using this stream,
  // If only one use, close and delete the stream.
  if (debugFlag[index].outputStream == NULL)
    return;

  int match = count_stream_users(debugFlag[index].outputStream);

  if (match == 1) {
    debugFlag[index].outputStream->close();
    delete debugFlag[index].outputStream;
    delete debugFlag[index].filename;
  }
  debugFlag[index].filename = NULL;
  debugFlag[index].outputStream = NULL;
}

void CubitMessage::set_message_handler(CubitMessageHandler *handler)
{
  CubitMessage::mHandler = handler;
  if(CubitMessage::mHandler == NULL)
    CubitMessage::mHandler = &mDefaultHandler;
}

CubitMessageHandler* CubitMessage::get_message_handler()
{
  return CubitMessage::mHandler;
}

void CubitMessage::set_error_handler(CubitMessageErrorHandler *handler)
{
  CubitMessage::mErrorHandler = handler;
}

CubitMessageErrorHandler* CubitMessage::get_error_handler()
{
  return CubitMessage::mErrorHandler;
}

MessageFlag::MessageFlag(int flag_number, const char *desc)
    : flagNumber(flag_number), setting(CUBIT_FALSE),
      description(desc), filename(NULL), outputStream(NULL)
{
}

MessageFlag::MessageFlag()
{
  flagNumber   = 0;
  setting      = CUBIT_FALSE;
  description  = NULL;
  filename     = NULL;
  outputStream = NULL;
}


void CubitMessage::set_logging_file_setting(const char* filename)
{
  if (loggingFile && *loggingFile == filename)
     return;

  if (CubitUtil::compare(filename,"terminal")) { // Filename is 'terminal'
    if (loggingStream != NULL) {
      loggingStream->close();
      delete loggingStream;
      loggingStream = NULL;
      delete loggingFile;
      loggingFile = NULL;
    }
    return;
  }
  else {
    loggingFile   = new CubitString(filename);
    loggingStream = new std::ofstream(CubitString::toNative(filename).c_str(), std::ios::out | std::ios::app );
  }
}

void CubitMessage::set_error_logging_file_setting(const char* filename, CubitBoolean resume_flag)
{
  if (loggingErrorFile && *loggingErrorFile == filename)
     return;

  if(loggingFile && *loggingFile == filename)
  {
    PRINT_ERROR("Can't explicitly set the Error logging file to be the same as the logging file.\n");
    return;
  }
 
  if (CubitUtil::compare(filename,"terminal")) { // Filename is 'terminal'
    if (loggingErrorStream != NULL) {
      loggingErrorStream->close();
      delete loggingErrorStream;
      loggingErrorStream = NULL;
      delete loggingErrorFile;
      loggingErrorFile = NULL;
    }
    return;
  }
  else {
    loggingErrorFile   = new CubitString(filename);
    if(resume_flag)
        loggingErrorStream = new std::ofstream(CubitString::toNative(filename).c_str(), std::ios::out | std::ios::app );
    else
        loggingErrorStream = new std::ofstream(CubitString::toNative(filename).c_str());
  }
}

//Initialize all settings in this class
void CubitMessage::initialize_settings()
{

  SettingHandler::instance()->add_setting("Info",
                                          CubitMessage::set_info_flag,
					  CubitMessage::get_info_flag);

  /*SettingHandler::instance()->add_setting("Logging",
                                          CubitMessage::set_logging_file_setting,
					  CubitMessage::get_logging_file_setting);*/

  SettingHandler::instance()->add_setting("Diagnostic",
					 CubitMessage::set_diagnostic_flag,
					 CubitMessage::get_diagnostic_flag);

  SettingHandler::instance()->add_setting("Warning",
					  CubitMessage::set_warning_flag,
					  CubitMessage::get_warning_flag);
}

CubitString CubitMessage::logging_filename() const
{
  CubitString temp_string;
  if(loggingStream != NULL)
    temp_string = loggingFile->c_str();

  return temp_string;
}

CubitString CubitMessage::logging_errors_filename() const
{
  CubitString temp_string;
  if(loggingErrorStream != NULL)
     temp_string = loggingErrorFile->c_str();

  return temp_string;
}

void CubitMessage::start_expected_error_count(int error_count, bool less_than_accepted)
{
  this->expectedStartErrorCount = errorCount;
  this->expectedEndErrorCount = errorCount + error_count;
  this->expectedLessErrorCountAccepted = less_than_accepted;
}

void CubitMessage::stop_expected_error_count(const CubitString &message)
{
  bool has_error = false;
  if((expectedEndErrorCount != expectedStartErrorCount && expectedLessErrorCountAccepted == false) ||
     (expectedEndErrorCount < expectedStartErrorCount && expectedLessErrorCountAccepted == true))
  {
    has_error = true;
  }
  this->expectedStartErrorCount = -1;
  this->expectedEndErrorCount = -1;
  if(has_error)
  {
    CubitString msg = message;
    if(message.length() == 0)
    {
      msg = "unexpected errors";
    }
    print_error(msg + "\n");
  }
}

MessageFlag::~MessageFlag()
{
  // It is not safe to delete either the stream or the filename here
  // since multiple instances (debug flags) may refer to the same memory
}

#ifdef STANDALONE
void main() {
CubitMessage::instance()->output_debug_information(1, 10, 2);
CubitMessage::instance()->output_debug_information(12);
CubitMessage::instance()->set_debug_file(5, "Debug_Test.file");
DEBUG_FLAG(5, CUBIT_TRUE);
PRINT_DEBUG_5("This is a test\n");
CubitMessage::instance()->output_debug_information(5,5);
}
#endif
