//- Class: CubitMessage
//- Description: CubitMessage class - used for reporting messages to
//-              the user
//- Owner: Tim Tautges
//- Checked By:
//- Version:

#ifndef CUBITMESSAGE_HPP
#define CUBITMESSAGE_HPP

// *****************
// Class Declaration
// *****************

#include <fstream>

#include <cstdarg>
#include "CubitDefines.h"
#include "CubitUtilConfigure.h"
#include "CubitString.hpp"
#include "CubitMessageHandler.hpp"

// CUBIT debug flags; see CubitMessage.cc for a definition of each flag

#define CUBIT_ERROR                        0
#define CUBIT_WARNING                      1
#define CUBIT_INFO                         2
#define CUBIT_DIAGNOSTIC                   3
#define CUBIT_ERROR_EXPECTED               4
#define CUBIT_DEBUG_1                     11

#define PRINT_DEBUG_1(...) PRINT_DEBUG(1, __VA_ARGS__)
#define PRINT_DEBUG_2(...) PRINT_DEBUG(2, __VA_ARGS__)
#define PRINT_DEBUG_3(...) PRINT_DEBUG(3, __VA_ARGS__)
#define PRINT_DEBUG_4(...) PRINT_DEBUG(4, __VA_ARGS__)
#define PRINT_DEBUG_5(...) PRINT_DEBUG(5, __VA_ARGS__)
#define PRINT_DEBUG_6(...) PRINT_DEBUG(6, __VA_ARGS__)
#define PRINT_DEBUG_7(...) PRINT_DEBUG(7, __VA_ARGS__)
#define PRINT_DEBUG_8(...) PRINT_DEBUG(8, __VA_ARGS__)
#define PRINT_DEBUG_9(...) PRINT_DEBUG(9, __VA_ARGS__)
#define PRINT_DEBUG_10(...) PRINT_DEBUG(10, __VA_ARGS__)
#define PRINT_DEBUG_11(...) PRINT_DEBUG(11, __VA_ARGS__)
#define PRINT_DEBUG_12(...) PRINT_DEBUG(12, __VA_ARGS__)
#define PRINT_DEBUG_13(...) PRINT_DEBUG(13, __VA_ARGS__)
#define PRINT_DEBUG_14(...) PRINT_DEBUG(14, __VA_ARGS__)
#define PRINT_DEBUG_15(...) PRINT_DEBUG(15, __VA_ARGS__)
#define PRINT_DEBUG_16(...) PRINT_DEBUG(16, __VA_ARGS__)
#define PRINT_DEBUG_17(...) PRINT_DEBUG(17, __VA_ARGS__)
#define PRINT_DEBUG_18(...) PRINT_DEBUG(18, __VA_ARGS__)
#define PRINT_DEBUG_19(...) PRINT_DEBUG(19, __VA_ARGS__)
#define PRINT_DEBUG_20(...) PRINT_DEBUG(20, __VA_ARGS__)
#define PRINT_DEBUG_21(...) PRINT_DEBUG(21, __VA_ARGS__)
#define PRINT_DEBUG_22(...) PRINT_DEBUG(22, __VA_ARGS__)
#define PRINT_DEBUG_23(...) PRINT_DEBUG(23, __VA_ARGS__)
#define PRINT_DEBUG_24(...) PRINT_DEBUG(24, __VA_ARGS__)
#define PRINT_DEBUG_25(...) PRINT_DEBUG(25, __VA_ARGS__)
#define PRINT_DEBUG_26(...) PRINT_DEBUG(26, __VA_ARGS__)
#define PRINT_DEBUG_27(...) PRINT_DEBUG(27, __VA_ARGS__)
#define PRINT_DEBUG_28(...) PRINT_DEBUG(28, __VA_ARGS__)
#define PRINT_DEBUG_29(...) PRINT_DEBUG(29, __VA_ARGS__)
#define PRINT_DEBUG_30(...) PRINT_DEBUG(30, __VA_ARGS__)
#define PRINT_DEBUG_31(...) PRINT_DEBUG(31, __VA_ARGS__)
#define PRINT_DEBUG_32(...) PRINT_DEBUG(32, __VA_ARGS__)
#define PRINT_DEBUG_33(...) PRINT_DEBUG(33, __VA_ARGS__)
#define PRINT_DEBUG_34(...) PRINT_DEBUG(34, __VA_ARGS__)
#define PRINT_DEBUG_35(...) PRINT_DEBUG(35, __VA_ARGS__)
#define PRINT_DEBUG_36(...) PRINT_DEBUG(36, __VA_ARGS__)
#define PRINT_DEBUG_37(...) PRINT_DEBUG(37, __VA_ARGS__)
#define PRINT_DEBUG_38(...) PRINT_DEBUG(38, __VA_ARGS__)
#define PRINT_DEBUG_39(...) PRINT_DEBUG(39, __VA_ARGS__)
#define PRINT_DEBUG_40(...) PRINT_DEBUG(40, __VA_ARGS__)
#define PRINT_DEBUG_41(...) PRINT_DEBUG(41, __VA_ARGS__)
#define PRINT_DEBUG_42(...) PRINT_DEBUG(42, __VA_ARGS__)
#define PRINT_DEBUG_43(...) PRINT_DEBUG(43, __VA_ARGS__)
#define PRINT_DEBUG_44(...) PRINT_DEBUG(44, __VA_ARGS__)
#define PRINT_DEBUG_45(...) PRINT_DEBUG(45, __VA_ARGS__)
#define PRINT_DEBUG_46(...) PRINT_DEBUG(46, __VA_ARGS__)
#define PRINT_DEBUG_47(...) PRINT_DEBUG(47, __VA_ARGS__)
#define PRINT_DEBUG_48(...) PRINT_DEBUG(48, __VA_ARGS__)
#define PRINT_DEBUG_49(...) PRINT_DEBUG(49, __VA_ARGS__)
#define PRINT_DEBUG_50(...) PRINT_DEBUG(50, __VA_ARGS__)
#define PRINT_DEBUG_51(...) PRINT_DEBUG(51, __VA_ARGS__)
#define PRINT_DEBUG_52(...) PRINT_DEBUG(52, __VA_ARGS__)
#define PRINT_DEBUG_53(...) PRINT_DEBUG(53, __VA_ARGS__)
#define PRINT_DEBUG_54(...) PRINT_DEBUG(54, __VA_ARGS__)
#define PRINT_DEBUG_55(...) PRINT_DEBUG(55, __VA_ARGS__)
#define PRINT_DEBUG_56(...) PRINT_DEBUG(56, __VA_ARGS__)
#define PRINT_DEBUG_57(...) PRINT_DEBUG(57, __VA_ARGS__)
#define PRINT_DEBUG_58(...) PRINT_DEBUG(58, __VA_ARGS__)
#define PRINT_DEBUG_59(...) PRINT_DEBUG(59, __VA_ARGS__)
#define PRINT_DEBUG_60(...) PRINT_DEBUG(60, __VA_ARGS__)
#define PRINT_DEBUG_61(...) PRINT_DEBUG(61, __VA_ARGS__)
#define PRINT_DEBUG_62(...) PRINT_DEBUG(62, __VA_ARGS__)
#define PRINT_DEBUG_63(...) PRINT_DEBUG(63, __VA_ARGS__)
#define PRINT_DEBUG_64(...) PRINT_DEBUG(64, __VA_ARGS__)
#define PRINT_DEBUG_65(...) PRINT_DEBUG(65, __VA_ARGS__)
#define PRINT_DEBUG_66(...) PRINT_DEBUG(66, __VA_ARGS__)
#define PRINT_DEBUG_67(...) PRINT_DEBUG(67, __VA_ARGS__)
#define PRINT_DEBUG_68(...) PRINT_DEBUG(68, __VA_ARGS__)
#define PRINT_DEBUG_69(...) PRINT_DEBUG(69, __VA_ARGS__)
#define PRINT_DEBUG_70(...) PRINT_DEBUG(70, __VA_ARGS__)
#define PRINT_DEBUG_71(...) PRINT_DEBUG(71, __VA_ARGS__)
#define PRINT_DEBUG_72(...) PRINT_DEBUG(72, __VA_ARGS__)
#define PRINT_DEBUG_73(...) PRINT_DEBUG(73, __VA_ARGS__)
#define PRINT_DEBUG_74(...) PRINT_DEBUG(74, __VA_ARGS__)
#define PRINT_DEBUG_75(...) PRINT_DEBUG(75, __VA_ARGS__)
#define PRINT_DEBUG_76(...) PRINT_DEBUG(76, __VA_ARGS__)
#define PRINT_DEBUG_77(...) PRINT_DEBUG(77, __VA_ARGS__)
#define PRINT_DEBUG_78(...) PRINT_DEBUG(78, __VA_ARGS__)
#define PRINT_DEBUG_79(...) PRINT_DEBUG(79, __VA_ARGS__)
#define PRINT_DEBUG_80(...) PRINT_DEBUG(80, __VA_ARGS__)
#define PRINT_DEBUG_81(...) PRINT_DEBUG(81, __VA_ARGS__)
#define PRINT_DEBUG_82(...) PRINT_DEBUG(82, __VA_ARGS__)
#define PRINT_DEBUG_83(...) PRINT_DEBUG(83, __VA_ARGS__)
#define PRINT_DEBUG_84(...) PRINT_DEBUG(84, __VA_ARGS__)
#define PRINT_DEBUG_85(...) PRINT_DEBUG(85, __VA_ARGS__)
#define PRINT_DEBUG_86(...) PRINT_DEBUG(86, __VA_ARGS__)
#define PRINT_DEBUG_87(...) PRINT_DEBUG(87, __VA_ARGS__)
#define PRINT_DEBUG_88(...) PRINT_DEBUG(88, __VA_ARGS__)
#define PRINT_DEBUG_89(...) PRINT_DEBUG(89, __VA_ARGS__)
#define PRINT_DEBUG_90(...) PRINT_DEBUG(90, __VA_ARGS__)
#define PRINT_DEBUG_91(...) PRINT_DEBUG(91, __VA_ARGS__)
#define PRINT_DEBUG_92(...) PRINT_DEBUG(92, __VA_ARGS__)
#define PRINT_DEBUG_93(...) PRINT_DEBUG(93, __VA_ARGS__)
#define PRINT_DEBUG_94(...) PRINT_DEBUG(94, __VA_ARGS__)
#define PRINT_DEBUG_95(...) PRINT_DEBUG(95, __VA_ARGS__)
#define PRINT_DEBUG_96(...) PRINT_DEBUG(96, __VA_ARGS__)
#define PRINT_DEBUG_97(...) PRINT_DEBUG(97, __VA_ARGS__)
#define PRINT_DEBUG_98(...) PRINT_DEBUG(98, __VA_ARGS__)
#define PRINT_DEBUG_99(...) PRINT_DEBUG(99, __VA_ARGS__)
#define PRINT_DEBUG_100(...) PRINT_DEBUG(100, __VA_ARGS__)
#define PRINT_DEBUG_101(...) PRINT_DEBUG(101, __VA_ARGS__)
#define PRINT_DEBUG_102(...) PRINT_DEBUG(102, __VA_ARGS__)
#define PRINT_DEBUG_103(...) PRINT_DEBUG(103, __VA_ARGS__)
#define PRINT_DEBUG_104(...) PRINT_DEBUG(104, __VA_ARGS__)
#define PRINT_DEBUG_105(...) PRINT_DEBUG(105, __VA_ARGS__)
#define PRINT_DEBUG_106(...) PRINT_DEBUG(106, __VA_ARGS__)
#define PRINT_DEBUG_107(...) PRINT_DEBUG(107, __VA_ARGS__)
#define PRINT_DEBUG_108(...) PRINT_DEBUG(108, __VA_ARGS__)
#define PRINT_DEBUG_109(...) PRINT_DEBUG(109, __VA_ARGS__)
#define PRINT_DEBUG_110(...) PRINT_DEBUG(110, __VA_ARGS__)
#define PRINT_DEBUG_111(...) PRINT_DEBUG(111, __VA_ARGS__)
#define PRINT_DEBUG_112(...) PRINT_DEBUG(112, __VA_ARGS__)
#define PRINT_DEBUG_113(...) PRINT_DEBUG(113, __VA_ARGS__)
#define PRINT_DEBUG_114(...) PRINT_DEBUG(114, __VA_ARGS__)
#define PRINT_DEBUG_115(...) PRINT_DEBUG(115, __VA_ARGS__)
#define PRINT_DEBUG_116(...) PRINT_DEBUG(116, __VA_ARGS__)
#define PRINT_DEBUG_117(...) PRINT_DEBUG(117, __VA_ARGS__)
#define PRINT_DEBUG_118(...) PRINT_DEBUG(118, __VA_ARGS__)
#define PRINT_DEBUG_119(...) PRINT_DEBUG(119, __VA_ARGS__)
#define PRINT_DEBUG_120(...) PRINT_DEBUG(120, __VA_ARGS__)
#define PRINT_DEBUG_121(...) PRINT_DEBUG(121, __VA_ARGS__)
#define PRINT_DEBUG_122(...) PRINT_DEBUG(122, __VA_ARGS__)
#define PRINT_DEBUG_123(...) PRINT_DEBUG(123, __VA_ARGS__)
#define PRINT_DEBUG_124(...) PRINT_DEBUG(124, __VA_ARGS__)
#define PRINT_DEBUG_125(...) PRINT_DEBUG(125, __VA_ARGS__)
#define PRINT_DEBUG_126(...) PRINT_DEBUG(126, __VA_ARGS__)
#define PRINT_DEBUG_127(...) PRINT_DEBUG(127, __VA_ARGS__)
#define PRINT_DEBUG_128(...) PRINT_DEBUG(128, __VA_ARGS__)
#define PRINT_DEBUG_129(...) PRINT_DEBUG(129, __VA_ARGS__)
#define PRINT_DEBUG_130(...) PRINT_DEBUG(130, __VA_ARGS__)
#define PRINT_DEBUG_131(...) PRINT_DEBUG(131, __VA_ARGS__)
#define PRINT_DEBUG_132(...) PRINT_DEBUG(132, __VA_ARGS__)
#define PRINT_DEBUG_133(...) PRINT_DEBUG(133, __VA_ARGS__)
#define PRINT_DEBUG_134(...) PRINT_DEBUG(134, __VA_ARGS__)
#define PRINT_DEBUG_135(...) PRINT_DEBUG(135, __VA_ARGS__)
#define PRINT_DEBUG_136(...) PRINT_DEBUG(136, __VA_ARGS__)
#define PRINT_DEBUG_137(...) PRINT_DEBUG(137, __VA_ARGS__)
#define PRINT_DEBUG_138(...) PRINT_DEBUG(138, __VA_ARGS__)
#define PRINT_DEBUG_139(...) PRINT_DEBUG(139, __VA_ARGS__)
#define PRINT_DEBUG_140(...) PRINT_DEBUG(140, __VA_ARGS__)
#define PRINT_DEBUG_141(...) PRINT_DEBUG(141, __VA_ARGS__)
#define PRINT_DEBUG_142(...) PRINT_DEBUG(142, __VA_ARGS__)
#define PRINT_DEBUG_143(...) PRINT_DEBUG(143, __VA_ARGS__)
#define PRINT_DEBUG_144(...) PRINT_DEBUG(144, __VA_ARGS__)
#define PRINT_DEBUG_145(...) PRINT_DEBUG(145, __VA_ARGS__)
#define PRINT_DEBUG_146(...) PRINT_DEBUG(146, __VA_ARGS__)
#define PRINT_DEBUG_147(...) PRINT_DEBUG(147, __VA_ARGS__)
#define PRINT_DEBUG_148(...) PRINT_DEBUG(148, __VA_ARGS__)
#define PRINT_DEBUG_149(...) PRINT_DEBUG(149, __VA_ARGS__)
#define PRINT_DEBUG_150(...) PRINT_DEBUG(150, __VA_ARGS__)
#define PRINT_DEBUG_151(...) PRINT_DEBUG(151, __VA_ARGS__)
#define PRINT_DEBUG_152(...) PRINT_DEBUG(152, __VA_ARGS__)
#define PRINT_DEBUG_153(...) PRINT_DEBUG(153, __VA_ARGS__)
#define PRINT_DEBUG_154(...) PRINT_DEBUG(154, __VA_ARGS__)
#define PRINT_DEBUG_155(...) PRINT_DEBUG(155, __VA_ARGS__)
#define PRINT_DEBUG_156(...) PRINT_DEBUG(156, __VA_ARGS__)
#define PRINT_DEBUG_157(...) PRINT_DEBUG(157, __VA_ARGS__)
#define PRINT_DEBUG_158(...) PRINT_DEBUG(158, __VA_ARGS__)
#define PRINT_DEBUG_159(...) PRINT_DEBUG(159, __VA_ARGS__)
#define PRINT_DEBUG_160(...) PRINT_DEBUG(160, __VA_ARGS__)
#define PRINT_DEBUG_161(...) PRINT_DEBUG(161, __VA_ARGS__)
#define PRINT_DEBUG_162(...) PRINT_DEBUG(162, __VA_ARGS__)
#define PRINT_DEBUG_163(...) PRINT_DEBUG(163, __VA_ARGS__)
#define PRINT_DEBUG_164(...) PRINT_DEBUG(164, __VA_ARGS__)
#define PRINT_DEBUG_165(...) PRINT_DEBUG(165, __VA_ARGS__)
#define PRINT_DEBUG_166(...) PRINT_DEBUG(166, __VA_ARGS__)
#define PRINT_DEBUG_167(...) PRINT_DEBUG(167, __VA_ARGS__)
#define PRINT_DEBUG_168(...) PRINT_DEBUG(168, __VA_ARGS__)
#define PRINT_DEBUG_169(...) PRINT_DEBUG(169, __VA_ARGS__)
#define PRINT_DEBUG_170(...) PRINT_DEBUG(170, __VA_ARGS__)
#define PRINT_DEBUG_171(...) PRINT_DEBUG(171, __VA_ARGS__)
#define PRINT_DEBUG_172(...) PRINT_DEBUG(172, __VA_ARGS__)
#define PRINT_DEBUG_173(...) PRINT_DEBUG(173, __VA_ARGS__)
#define PRINT_DEBUG_174(...) PRINT_DEBUG(174, __VA_ARGS__)
#define PRINT_DEBUG_175(...) PRINT_DEBUG(175, __VA_ARGS__)
#define PRINT_DEBUG_176(...) PRINT_DEBUG(176, __VA_ARGS__)
#define PRINT_DEBUG_177(...) PRINT_DEBUG(177, __VA_ARGS__)
#define PRINT_DEBUG_178(...) PRINT_DEBUG(178, __VA_ARGS__)
#define PRINT_DEBUG_179(...) PRINT_DEBUG(179, __VA_ARGS__)
#define PRINT_DEBUG_180(...) PRINT_DEBUG(180, __VA_ARGS__)
#define PRINT_DEBUG_181(...) PRINT_DEBUG(181, __VA_ARGS__)
#define PRINT_DEBUG_182(...) PRINT_DEBUG(182, __VA_ARGS__)
#define PRINT_DEBUG_183(...) PRINT_DEBUG(183, __VA_ARGS__)
#define PRINT_DEBUG_184(...) PRINT_DEBUG(184, __VA_ARGS__)
#define PRINT_DEBUG_185(...) PRINT_DEBUG(185, __VA_ARGS__)
#define PRINT_DEBUG_186(...) PRINT_DEBUG(186, __VA_ARGS__)
#define PRINT_DEBUG_187(...) PRINT_DEBUG(187, __VA_ARGS__)
#define PRINT_DEBUG_188(...) PRINT_DEBUG(188, __VA_ARGS__)
#define PRINT_DEBUG_189(...) PRINT_DEBUG(189, __VA_ARGS__)
#define PRINT_DEBUG_190(...) PRINT_DEBUG(190, __VA_ARGS__)
#define PRINT_DEBUG_191(...) PRINT_DEBUG(191, __VA_ARGS__)
#define PRINT_DEBUG_192(...) PRINT_DEBUG(192, __VA_ARGS__)
#define PRINT_DEBUG_193(...) PRINT_DEBUG(193, __VA_ARGS__)
#define PRINT_DEBUG_194(...) PRINT_DEBUG(194, __VA_ARGS__)
#define PRINT_DEBUG_195(...) PRINT_DEBUG(195, __VA_ARGS__)
#define PRINT_DEBUG_196(...) PRINT_DEBUG(196, __VA_ARGS__)
#define PRINT_DEBUG_197(...) PRINT_DEBUG(197, __VA_ARGS__)
#define PRINT_DEBUG_198(...) PRINT_DEBUG(198, __VA_ARGS__)
#define PRINT_DEBUG_199(...) PRINT_DEBUG(199, __VA_ARGS__)
#define PRINT_DEBUG_200(...) PRINT_DEBUG(200, __VA_ARGS__)
#define PRINT_DEBUG_201(...) PRINT_DEBUG(201, __VA_ARGS__)
#define PRINT_DEBUG_202(...) PRINT_DEBUG(202, __VA_ARGS__)
#define PRINT_DEBUG_203(...) PRINT_DEBUG(203, __VA_ARGS__)
#define PRINT_DEBUG_204(...) PRINT_DEBUG(204, __VA_ARGS__)
#define PRINT_DEBUG_205(...) PRINT_DEBUG(205, __VA_ARGS__)
#define PRINT_DEBUG_206(...) PRINT_DEBUG(206, __VA_ARGS__)
#define PRINT_DEBUG_207(...) PRINT_DEBUG(207, __VA_ARGS__)
#define PRINT_DEBUG_208(...) PRINT_DEBUG(208, __VA_ARGS__)
#define PRINT_DEBUG_209(...) PRINT_DEBUG(209, __VA_ARGS__)
#define PRINT_DEBUG_210(...) PRINT_DEBUG(210, __VA_ARGS__)
#define PRINT_DEBUG_211(...) PRINT_DEBUG(211, __VA_ARGS__)
#define PRINT_DEBUG_212(...) PRINT_DEBUG(212, __VA_ARGS__)
#define PRINT_DEBUG_213(...) PRINT_DEBUG(213, __VA_ARGS__)
#define PRINT_DEBUG_214(...) PRINT_DEBUG(214, __VA_ARGS__)
#define PRINT_DEBUG_215(...) PRINT_DEBUG(215, __VA_ARGS__)
#define PRINT_DEBUG_216(...) PRINT_DEBUG(216, __VA_ARGS__)
#define PRINT_DEBUG_217(...) PRINT_DEBUG(217, __VA_ARGS__)
#define PRINT_DEBUG_218(...) PRINT_DEBUG(218, __VA_ARGS__)
#define PRINT_DEBUG_219(...) PRINT_DEBUG(219, __VA_ARGS__)
#define PRINT_DEBUG_220(...) PRINT_DEBUG(220, __VA_ARGS__)
#define PRINT_DEBUG_221(...) PRINT_DEBUG(221, __VA_ARGS__)
#define PRINT_DEBUG_222(...) PRINT_DEBUG(222, __VA_ARGS__)
#define NUM_DEBUG_FLAGS 224

#define PRINT_ERROR(...) CubitMessage::instance()->print_error(CubitString::format(__VA_ARGS__))
#define PRINT_WARNING(...) CubitMessage::instance()->print_warning(CubitString::format(__VA_ARGS__))
#define PRINT_INFO(...) CubitMessage::instance()->print_info(CubitString::format(__VA_ARGS__))
#define DIAGNOSTIC(...) CubitMessage::instance()->print_diagnostic(CubitString::format(__VA_ARGS__))
#define DIAGNOSTIC_FLAG CubitMessage::instance()->get_diagnostic_flag
#define DEBUG_FLAG CubitMessage::instance()->debug_flag
#define GET_INFO_FLAG CubitMessage::instance()->get_info_flag
#define SET_INFO_FLAG CubitMessage::instance()->set_info_flag
#define SET_WARNING_FLAG CubitMessage::instance()->set_warning_flag
#define GET_WARNING_FLAG CubitMessage::instance()->get_warning_flag
#define SET_ERROR_FLAG CubitMessage::instance()->set_error_flag
#define GET_ERROR_FLAG CubitMessage::instance()->get_error_flag
#define DEBUG_FLAG_SET CubitMessage::instance()->is_debug_flag_set
#define PRINT_DEBUG(x, ...) if(DEBUG_FLAG_SET(x)) CubitMessage::instance()->print_debug(CubitString::format(__VA_ARGS__))
#define PRINT_FREE CubitMessage::free_instance()

class CubitString;
class CubitMessage;

class CUBIT_UTIL_EXPORT MessageFlag
{
  friend class CubitMessage;
public:
  ~MessageFlag();
private:
  MessageFlag();
  MessageFlag(int flag_number, const char *desc);

  void output();

    // Member variables
  int flagNumber;
  int setting;
  const char *description;
  CubitString *filename;
  std::ofstream *outputStream;
};

class CUBIT_UTIL_EXPORT CubitMessage
{
protected:

  static CubitMessage* instance_;
  //- static pointer to unique instance of this class

  static CubitMessageHandler* mHandler;
  //- static pointer to the message output handler
  
  static CubitMessageErrorHandler* mErrorHandler;
  //- static pointer to the message context handler

  static MessageFlag staticDebugFlag[];

  MessageFlag *debugFlag;
  //- debug flag, used with internal_error

  static int infoFlag;
  //- info flag, used with internal_error

  static int warningFlag;
  //- warning flag, used with internal_error

  static int errorFlag;
  //- error flag, used with internal_error

  static int diagnosticFlag;
  //- diagnostic flag, used with internal_error

  int currentDebugFlag;

  static int errorCount;
  //- static variable to track the errors that occured in the
  //- a session.  Only gets set when PRINT_ERROR is called.

  static int expectedStartErrorCount;
  static int expectedEndErrorCount;
  static bool expectedLessErrorCountAccepted;

  static int warningCount;
  //- static variable to track the warnings that occured in the
  //- a session.  Only gets set when PRINT_WARNING is called.

  static std::ofstream *loggingStream;
  static CubitString *loggingFile;
  //- Stream pointer for logging of internal_error messages.
  //- If NULL, output goes to terminal only (except for debug)
  //- If Non-NULL, output goes to both terminal and stream.

  static std::ofstream *loggingErrorStream;
  static CubitString *loggingErrorFile;
  //- Stream pointer for logging of only ERROR messages.
  //- If NULL, ERROR output goes to normal places only
  //- If Non-NULL, output goes to both this stream, and all other places.

  void add_to_error_count();
  //- Increments the errorCount variable. Keep private (GDS).

  void add_to_warning_count();
  //- Increments the errorCount variable. Keep private (GDS).
  void set_debug_stream(const int index, std::ofstream *output_stream);
  //- Set the output stream for this debug flag to output_stream.

  void remove_debug_stream(const int index);
  //- Close and delete the stream if only one use.

  int find_file_use(const CubitString &filename);
  int count_stream_users(const std::ofstream *stream);

  CubitMessage ();
    //- Class Constructor. (Not callable by user code. Class is constructed
    //- by the {instance()} member function.

public:

  static CubitMessage* instance();
  static void free_instance();
  //- Controlled access and creation of the sole instance of this class.

  virtual ~CubitMessage();
  //- Class Destructor.

  static void delete_instance();

  void set_logging_file_setting(const CubitString &filename, CubitBoolean resume_flag = CUBIT_FALSE);
  void set_debug_file_setting(const int index, const CubitString &filename);

  CubitString logging_filename() const;
  CubitString logging_errors_filename() const;

  int is_debug_flag_set( int flag );
  //- for use with the PRINT_DEBUG macro only

  //static int get_debug_for_setting_handler(const int index)
  //{return staticDebugFlag[index].setting;};
  // static void set_debug_for_setting_handler(const int index, const int value)
  //{staticDebugFlag[index].setting = value;};
  //- for use with SettingHandler.cpp code only

  int  debug_flag(const int index);
  void debug_flag(const int index, const int flag);
  int  number_of_debug_flags();
  //- debug flag, used with internal_error

  virtual void set_debug_flag_gui(bool /*flag*/){};
  virtual int is_debug_flag_gui_set(){return 0;};
  virtual int print_debug_gui( const char* /* format*/ , ... ){return 0;};
  //- write out a debug message (from GUI only)
  //- used for GUI Debugging (CAT-only)

  static bool get_info_flag();
  static void set_info_flag(bool flag);
  //- info flag, used with internal_error

  static bool get_warning_flag();
  static void set_warning_flag(bool flag);
  //- warning flag, used with internal_error

  static bool get_error_flag();
  static void set_error_flag(bool flag);
  //- error flag, used with internal_error

  static bool get_diagnostic_flag();
  static void set_diagnostic_flag(bool flag);
  //- diagnostic flag, used with internal_error

  virtual void internal_error(const int message_type, std::ofstream *output_stream,
                              const CubitString& str);
  //- write out a debug/info/error/warning message

  int print_error(const CubitString& str);
  //- write out an error message

  int print_warning(const CubitString& str);
  //- write out a warning message

  int print_info(const CubitString& str);
  //- write out an info message

  int print_debug(const CubitString& str);
  //- write out a debug message

  void print_diagnostic(const CubitString& str);
  //- write out a diagnostic message

  int reset_error_count(int value = 0);
  //- Sets the errorCount variable to 0;
  //- Returns current value of variable.

  int error_count();
  //- Returns the value of the errorCount variable;
  //- This errorCount variable is incremented only if print_error is called;
  //- there is not a public  interface to only set the flag.
  //- My reasoning for that is that there should be
  //- some notification of an error so that the
  //- user can figure out why this function returns TRUE.

  int reset_warning_count(int value = 0);
  //- Sets the warningCount variable to 0;
  //- Returns current value of variable.

  int warning_count();
  //- Returns the value of the warningCount variable;
  //- This warningCount variable is incremented only if print_warning is called;
  //- there is not a public  interface to only set the flag.
  //- My reasoning for that is that there should be
  //- some notification of an error so that the
  //- user can figure out why this function returns TRUE.

  void output_debug_information(int from=1, int to=-1, int step=1);
  void output_debug_information(CubitString &match);
  void output_logging_information();

  static char* get_logging_file_setting();
  static void set_logging_file_setting(const char* file);
  static void set_error_logging_file_setting(const char* file, CubitBoolean resume_flag = CUBIT_FALSE);

  static void initialize_settings();

  static void set_message_handler(CubitMessageHandler *handler);
  static CubitMessageHandler* get_message_handler();
  
  static void set_error_handler(CubitMessageErrorHandler *handler);
  static CubitMessageErrorHandler* get_error_handler();

  //! start an expected error count
  //! also set whether an error count less than that is acceptable at stop time
  //! X number of errors will be suppressed, and will print out with a prefix of EXPECTED_ERROR:
  void start_expected_error_count(int error_count, bool less_than_accepted);
  //! stop an expected error count
  //! if the remaining number of expected errors is zero, we simply clean up.
  //! otherwise, an ERROR is reported
  //! if the CUBIT_CTEST environment variable is set, a cdash formatted message will be printed out
  void stop_expected_error_count(const CubitString& message);

}; // End of Class CubitMessage

inline int
CubitMessage::debug_flag(const int index)
{return debugFlag[index].setting;}

inline void
CubitMessage::debug_flag(const int index, const int flag)
{debugFlag[index].setting = flag;}

inline bool
CubitMessage::get_info_flag()
{return !!infoFlag;}

inline void
CubitMessage::set_info_flag(bool flag)
{infoFlag = flag;}

inline bool
CubitMessage::get_warning_flag()
{return !!warningFlag;}

inline void
CubitMessage::set_warning_flag(bool flag)
{warningFlag = flag;}

inline bool
CubitMessage::get_error_flag()
{return !!errorFlag;}

inline void
CubitMessage::set_error_flag(bool flag)
{errorFlag = flag;}

inline bool
CubitMessage::get_diagnostic_flag()
{return !!diagnosticFlag;}

inline void
CubitMessage::set_diagnostic_flag(bool flag)
{diagnosticFlag = flag;}

#endif

