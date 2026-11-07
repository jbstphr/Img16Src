// in might includes either
//  <along with stdafx.cpp> into project of VStudio-toolchain for precompiled headers 
//  <alone> into general project of other toolchains
#ifndef I16_INC_H_INCLUDED
#define I16_INC_H_INCLUDED
  #include <stdio.h>
  #include <memory.h>
  #include <stdlib.h>
  #include <sys/stat.h>
  #include <errno.h>  
  #ifdef _WIN32
    #define ATTRIB_EXPORT
    #ifndef WINVER
    #define WINVER 0x0600
    #endif
    #ifndef _WIN32_WINNT
    #define _WIN32_WINNT 0x0600
    #endif
    #define _CRT_SECURE_NO_DEPRECATE
    #define _CRT_SECURE_NO_WARNINGS
    #define _CRT_NON_CONFORMING_SWPRINTFS
    #pragma conform(forScope, off)
    #pragma warning(disable: 4288 4430 4311 4313 4996)
    #define NO_WARN_MBCS_MFC_DEPRECATION
    #define _AFX_NO_MFC_CONTROLS_IN_DIALOGS
    #define strcasecmp stricmp
    #include <windows.h>
  #else // ifLinux
    #define ATTRIB_EXPORT __attribute__((visibility("default"))) 
    #include <string.h>
    #include <thread>
    #pragma GCC diagnostic ignored "-Wmismatched-new-delete"
    #pragma GCC diagnostic ignored "-Wmisleading-indentation"
    #pragma GCC diagnostic ignored "-Wunused-function"
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wunused-variable"
  #endif
  #define INTERNAL_I16_LIBRARY // procedures' declarations 'as is'
  #include "stdafv.h"
  #include "../I16_Cmn/I16_Intf.h"

#endif // I16_INC_H_INCLUDED
