// stdafx.h : include file for standard include files
#pragma once

#ifdef _WIN32
  // additional declarations
  #define _CRT_SECURE_NO_DEPRECATE
  #define _CRT_SECURE_NO_WARNINGS
  #define _CRT_NON_CONFORMING_SWPRINTFS
  #pragma conform(forScope, off)
  #pragma warning(disable: 4288 4430 4311 4313 4365 4668) // ~CRT
  #pragma warning(disable: 4996) // deperacted POSIX
  #define _WIN32_WINNT 0x0502 
  #define NO_WARN_MBCS_MFC_DEPRECATION
  #define _AFX_NO_MFC_CONTROLS_IN_DIALOGS
  #ifndef _WIN32_WINNT        // "targetver.h" Specifies that the minimum required platform is Windows Vista.
  #define _WIN32_WINNT 0x0600 // Change this to the appropriate value to target other versions of Windows.
  #endif
  #include <shlwapi.h>
  #include <direct.h>
  #include <windows.h>
  #pragma comment(lib, "shlwapi.lib")
#else
  #define __min(a,b) (a<b?a:b)
  #define __max(a,b) (a<b?b:a)
  #include <fnmatch.h>
  #include <unistd.h>
  #include <time.h>
  #include <dlfcn.h>
  #pragma GCC diagnostic ignored "-Wmisleading-indentation"
  #pragma GCC diagnostic ignored "-Wunused-function"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>  
