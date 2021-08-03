
#include "CubitProcess.hpp"
#include "CubitString.hpp"
#include "CubitFileUtil.hpp"
#include "CubitUtil.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <direct.h>

const char path_separator = '\\';
const char* path_separator_str = "\\";
#else
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <limits.h>
const char path_separator = '/';
const char* path_separator_str = "/";
#endif

#ifdef _WIN32
static CubitString quote_string(const CubitString& str)
{
  CubitString ret = "\"";  
  ret += str;
  ret += "\"";
  return ret;
}

static CubitString join_args(const std::vector<CubitString>& strings)
{
  CubitString joined_strings;
  for(unsigned int i=0; i<strings.size(); i++)
  {
    if(i!=0)
      joined_strings += " ";
    
    bool need_quote = strings[i].find(" ") != CubitString::npos;
    if(need_quote)
      joined_strings += "\"";
    joined_strings += strings[i];
    if(need_quote)
      joined_strings += "\"";
  }
  return joined_strings;
}
#endif

// given a relative or absolute path, make it absolute
static CubitString make_path_absolute(const CubitString& p)
{
  CubitString ret;
  if(!p.is_empty())
  {
    if ( CubitFileUtil::is_absolute(p.c_str()) )
    {
      ret = p;
    }
    else
    {
      // if any '/' character is in it, its relative to current directory
      CubitString wd;
      CubitFileUtil::get_current_working_directory(wd);
      ret = wd.c_str();
      ret += path_separator_str;
      ret += p;
    }
  }
  return ret;
}

#ifndef _WIN32
static sigset_t oldsig;
#endif

PidType CubitProcess::start(const CubitString& app, const std::vector<CubitString>& args, bool hide)
{
#ifndef _WIN32
  (void)hide;

  std::vector<const char*> c_args(args.size()+2);
  int idx = 0;
  c_args[idx++] = app.c_str();
  for(unsigned int i=0; i<args.size(); i++)
  {
    c_args[idx++] = args[i].c_str();
  }
  c_args[idx++] = NULL;
  
  CubitString app_real = find_executable(app);

  // temporarily block delivery of child signals
  // TODO: does this overwrite currently set signals for the process?
  // ANSWER: Yes it does, which is why we comment this out because it
  // can prevent the reaping of completed child processes.
  /*
  sigset_t newsig;
  sigemptyset(&newsig);
  sigaddset(&newsig, SIGCHLD);
  sigprocmask(SIG_BLOCK, &newsig, &oldsig);
  */

  pid_t pid = fork();
  if(pid < 0)
    return 0;

  if(pid == 0)
  {
    execv(app_real.c_str(), const_cast<char**>(&c_args[0]));
    perror(app_real.c_str());
    exit(EXIT_FAILURE);
  }

  return pid;
#else

  STARTUPINFOW si;
  PROCESS_INFORMATION pi;
  
  ZeroMemory( &si, sizeof(si) );
  si.cb = sizeof(si);
  ZeroMemory( &pi, sizeof(pi) );
  
  // hide child window
  if(hide)
  {
    si.dwFlags |= STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
  }

  CubitString real_app = find_executable(app);
  real_app.replace("/", "\\");
  
  CubitString call = quote_string(app);
  CubitString joined_args = join_args(args);
  call += " ";
  call += joined_args;

    // Start the child process. 
  if( CreateProcessW( CubitString::toUtf16(real_app).c_str(),  // path to cubit executable 
                      const_cast<wchar_t*>(CubitString::toUtf16(call).c_str()), // Command line. 
                      NULL,      // Process handle not inheritable. 
                      NULL,      // Thread handle not inheritable. 
                      TRUE,     // Set handle inheritance to TRUE.
                      0,         // No creation flags. 
                      NULL,      // Use parent's environment block. 
                      NULL,      // Use parent's starting directory. 
                      &si,       // Pointer to STARTUPINFO structure.
                      &pi )      // Pointer to PROCESS_INFORMATION structure.
      )
  {
    return pi;
  }
  pi.hProcess = 0;
  return pi;
#endif
}

int CubitProcess::wait(PidType pid)
{
#ifndef _WIN32
  int status;
  int result = -1;
  
  while(1)
  {
    int change_pid = waitpid(-1, &status, 0);
    if(change_pid == pid)
    {
      if(WIFEXITED(status))
      {
        result = WEXITSTATUS(status);
        break;
      }
      else if (WIFSIGNALED(status))
      {
        result = 128 + WTERMSIG(status);
        break;
      }
    }
  }
  
  sigprocmask(SIG_SETMASK, &oldsig, NULL);

  return result;
#else
  
  // Wait until child process exits.
  WaitForSingleObject( pid.hProcess, INFINITE );
    
  // Get its exit code
  DWORD exit_code = -1;
  GetExitCodeProcess(pid.hProcess, &exit_code);
    
  // Close process and thread handles. 
  CloseHandle( pid.hProcess );
  CloseHandle( pid.hThread );
  return exit_code;

#endif
}


int CubitProcess::execute(const CubitString& app, const std::vector<CubitString>& args, bool hide)
{
  PidType pid = CubitProcess::start(app, args, hide);
#if _WIN32
  if(pid.hProcess == 0)
    return -1;
#else
  if(pid == 0)
    return -1;
#endif
  return CubitProcess::wait(pid);
}
  
CubitString CubitProcess::find_executable(const CubitString& app)
{
  CubitString app_real = app;

#ifdef _WIN32
  if(!app_real.ends_with(".exe"))
	  app_real += ".exe";
  const char path_delimiter = ';';
#else
  const char path_delimiter = ':';
#endif

#ifdef _WIN32
  if(app_real.get_at(0) != path_separator && 
	  (app_real.length() > 1 && app_real.get_at(1) != ':'))
#else
  if(app_real.get_at(0) != path_separator)
#endif
  {
    if(app_real.find(path_separator_str) != CubitString::npos)
    {
      // if any '/' character is in it, its relative to current directory
      app_real = make_path_absolute(app_real);
    }
    else
    {
      // search PATH env. for this executable (note PATH may have relative directories)
      std::vector<CubitString> paths;
      std::stringstream ss(CubitUtil::getenv("PATH").c_str());
      std::string item;
      while(std::getline(ss, item, path_delimiter))
      {
        paths.push_back(item.c_str());
      }

      for(size_t i=0; i<paths.size(); i++)
      {
        CubitString p = paths[i];
        p += CubitFileUtil::separator();
        p += app_real;
        CubitString abs_p = make_path_absolute(p);
        if(CubitFileUtil::path_exists(abs_p.c_str()))
        {
          app_real = abs_p;
          break;
        }
      }
    }
  }
  return app_real;
}
