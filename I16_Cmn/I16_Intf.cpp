/** *************************************************************** **
 **                      I16 Process library                        **
 **  edition : Jan-2020, Mar-2026                       by S.B.     **
 This is DLL's loader for dynamic link a cpp-project against
 "I16Prc.dll". It is required when a process maps two or more
 different DLLs with the same SW-interface for different purposes (e.g.
 one for fast preview, and another for saving still images). Include
 this file "I16_Intf.cpp" into the project (along with "I16_Intf.h").
 The documentation is named "I16Prc.rtf".
 ** *************************************************************** **/

#include "stdafx.h"
  // Follow good MFC's pattern of unified app-wide precompiled header. Certain 
  // "toolchain+filesytems" might demand case-sensitive ~"stdafx", "StdAfx".

#include "I16_Intf.h"
#ifdef _WIN32 // Windows
  static void iErrMsg( const char *AAsset, int AErrCode ) {
    char qMsg[256]; sprintf_s( qMsg, sizeof(qMsg), "ERR: I16_load asset '%s'", AAsset );
    if (AErrCode) { char qErr[92], *qCh = qErr+sprintf_s( qErr, sizeof(qErr), "\n#%d ", AErrCode );
      if (FormatMessage( FORMAT_MESSAGE_FROM_SYSTEM, NULL, AErrCode, 0, qCh, (int)(sizeof(qErr)-(qCh-qErr)),
        NULL )<=0) strcat_s( qCh, (int)(sizeof(qErr)-(qCh-qErr)), "(error code)" );
      strcat_s( qMsg, sizeof(qMsg), qErr ); }; // if AErrCode
    // MessageBoxA( NULL, qMsg, "VPF-loader error", MB_OK|MB_ICONSTOP );
    printf( "%s\n", qMsg ); }

  static bool iGetPA( tI16_Intf *AIntf, const char *AFuncName, FARPROC *AFuncRef ) {
    *AFuncRef = GetProcAddress( AIntf->rHDLL, AFuncName ); if (*AFuncRef) return true;
    iErrMsg( AFuncName, 0 ); return false; }

  bool I16_MapDLL( const char *ADllPath, tI16_Intf *AIntf ) {
    memset( AIntf, 0, sizeof(tI16_Intf) ); const char
    #ifndef _WIN64
      *cDbgDmp = "?I16_DbgDmp@@YAXPBD@Z", *cInvoke = "?I16_Invoke@@YAHHHPBFHPAFH@Z", *cParStr = "?I16_ParStr@@YAPADH@Z";
    #else
      *cDbgDmp = "?I16_DbgDmp@@YAXPEBD@Z", *cInvoke = "?I16_Invoke@@YAHHHPEBFHPEAFH@Z", *cParStr = "?I16_ParStr@@YAPEADH@Z";
    #endif
    const char *cDllNam = strrchr( ADllPath, '\\' ); char qSavDir[192];
    if (!cDllNam) cDllNam = ADllPath; else { GetCurrentDirectoryA( 191, qSavDir ); char qDllDir[192];
        strcpy( qDllDir, ADllPath ); qDllDir[cDllNam++-ADllPath] = 0; SetCurrentDirectoryA( qDllDir ); };
    AIntf->rHDLL = LoadLibraryA( cDllNam ); if (cDllNam!=ADllPath) SetCurrentDirectoryA( qSavDir );
    if (!AIntf->rHDLL) { iErrMsg( ADllPath, GetLastError() ); return false; };
    if (!iGetPA( AIntf, "?I16_Init@@YAHHHH@Z", (FARPROC*)&AIntf->rInit )) goto labLeave;
    if (!iGetPA( AIntf, "?I16_Release@@YAXXZ", (FARPROC*)&AIntf->rRelease )) goto labLeave;
    if (!iGetPA( AIntf, "?I16_DbgDmpReset@@YAXXZ", (FARPROC*)&AIntf->rDbgDmpReset )) goto labLeave;
    if (!iGetPA( AIntf, cDbgDmp, (FARPROC*)&AIntf->rDbgDmp)) goto labLeave;
    if (!iGetPA( AIntf, cInvoke, (FARPROC*)&AIntf->rInvoke )) goto labLeave;
    if (!iGetPA( AIntf, "?I16_ParSet@@YAXHH@Z", (FARPROC*)&AIntf->rParSet )) goto labLeave;
    if (!iGetPA( AIntf, "?I16_ParGet@@YAHH@Z", (FARPROC*)&AIntf->rParGet )) goto labLeave;
    if (!iGetPA( AIntf, cParStr, (FARPROC*)&AIntf->rParStr)) goto labLeave;
    return true; labLeave: I16_ReleaseDLL( AIntf ); return false; }

  void I16_ReleaseDLL( tI16_Intf *AIntf ) {
    if (AIntf->rHDLL) { FreeLibrary( AIntf->rHDLL ); AIntf->rHDLL = NULL; }; }

#else // Linux (assume only 64-bit here)
  static void iErrMsg( const char *AAsset ) {
    char qMsg[256]; sprintf( qMsg, "I16_Load asset '%.99s' (%.99s)", AAsset, dlerror() );
    printf( "%s\n", qMsg ); }

  static bool iGetPA( tI16_Intf *AIntf, const char *AFuncName, void **AFuncRef ) {
    *AFuncRef = dlsym( AIntf->rHDLL, AFuncName ); if (*AFuncRef) return true;
    iErrMsg( AFuncName ); return false; }

  bool I16_MapDLL( const char *ADllPath, tI16_Intf *AIntf ) {
    memset( AIntf, 0, sizeof(tI16_Intf) );
    dlerror(); // clear existing errors
    AIntf->rHDLL = dlopen( ADllPath, RTLD_LAZY ); // Example: loading math lib
    if (!AIntf->rHDLL) { iErrMsg( ADllPath ); return false; };
    if (!iGetPA( AIntf, "_Z8I16_Initiii", (void**)&AIntf->rInit )) goto labLeave;
    if (!iGetPA( AIntf, "_Z11I16_Releasev", (void**)&AIntf->rRelease )) goto labLeave;
    if (!iGetPA( AIntf, "_Z15I16_DbgDmpResetv", (void**)&AIntf->rDbgDmpReset )) goto labLeave;
    if (!iGetPA( AIntf, "_Z10I16_DbgDmpPKc", (void**)&AIntf->rDbgDmp)) goto labLeave;
    if (!iGetPA( AIntf, "_Z10I16_InvokeiiPKsiPsi", (void**)&AIntf->rInvoke )) goto labLeave;
    if (!iGetPA( AIntf, "_Z10I16_ParSetii", (void**)&AIntf->rParSet )) goto labLeave;
    if (!iGetPA( AIntf, "_Z10I16_ParGeti", (void**)&AIntf->rParGet )) goto labLeave;
    if (!iGetPA( AIntf, "_Z10I16_ParStri", (void**)&AIntf->rParStr )) goto labLeave;
    return true; labLeave: I16_ReleaseDLL( AIntf ); return false; }

  void I16_ReleaseDLL( tI16_Intf *AIntf ) {
    if (AIntf->rHDLL) { dlclose( AIntf->rHDLL ); AIntf->rHDLL = NULL; }; }

#endif // Linux