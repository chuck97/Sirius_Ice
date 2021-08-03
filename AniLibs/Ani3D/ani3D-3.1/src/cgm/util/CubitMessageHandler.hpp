//- Class: CubitMessageHandler
//- Description: CubitMessageHandler class - used for reporting messages to
//-              the user

#ifndef CUBITMESSAGEHANDLER_HPP
#define CUBITMESSAGEHANDLER_HPP
#include "CubitUtilConfigure.h"
#include <string>

class CUBIT_UTIL_EXPORT CubitMessageHandler
{
public:
  CubitMessageHandler();
  virtual ~CubitMessageHandler();

  virtual void print_message(const char *message) = 0;
  virtual void print_error(const char *message) = 0;
};

class CUBIT_UTIL_EXPORT CubitMessageErrorHandler
{
public:
  CubitMessageErrorHandler();
  virtual ~CubitMessageErrorHandler();
  virtual std::string error_context() = 0;
};

#endif

