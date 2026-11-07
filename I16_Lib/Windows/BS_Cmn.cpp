// BS_Cmn.cpp : general-purpose utilities (OS,EABI-specific)
// for Windows by S.B.
#include "stdafx.h"
#include <conio.h>
#include "math.h"
#include "windows.h"
#include "BS_Cmn.h"

//========== Threading service ============
void BS_Mutex( int AOp, CRITICAL_SECTION *ACS ) { // -- 2:create, 1:enter, -1:leave, -2:delete
  switch (AOp) { case 2: InitializeCriticalSection( ACS ); break; case -2: DeleteCriticalSection( ACS ); break;
                 case 1: EnterCriticalSection     ( ACS ); break; case -1: LeaveCriticalSection ( ACS ); break; }; }

const int cBS_SpawnThreads = 16; int xBS_SpawnThreads = 4;
typedef struct { tBS_SpawnThreadRec rR; tBS_SpawnThreadProc rProc; HANDLE rThHn,rEvent[2]; } tBS_ThreadParam_RT;
static DWORD _stdcall iSpawnedThreadPlain( void *AParam ) { tBS_ThreadParam_RT *qParam = (tBS_ThreadParam_RT*)AParam; qParam->rProc( &qParam->rR ); return 0; };
static DWORD _stdcall iBS_SpawnedThreadProc( void *AParam ) {
  tBS_ThreadParam_RT *qParam = (tBS_ThreadParam_RT*)AParam; do { 
    WaitForSingleObject( qParam->rEvent[0], INFINITE ); if (qParam->rProc) try { qParam->rProc( &qParam->rR ); } catch(...) {}; SetEvent( qParam->rEvent[1] );
  } while (qParam->rProc); return 0; } // do

static tBS_ThreadParam_RT  xBS_ThreadParam_RT[cBS_SpawnThreads];
static HANDLE              xBS_ThreadEvWait[cBS_SpawnThreads];
static CRITICAL_SECTION    xBS_ThreadCS;
static int                 xBS_ThreadPool = -1; void BS_SpawnThreads_State( bool AOn ) { // -- alloc,dealloc threads' pool
  if (AOn) { BS_Mutex( 2, &xBS_ThreadCS ); xBS_ThreadPool = 0; return; }; // initialize threads' pool
  if (xBS_ThreadPool<0) return; BS_Mutex( -2, &xBS_ThreadCS ); // deinitialize threads' pool 
  for (int qi=0;qi<xBS_ThreadPool;qi++) {
    tBS_ThreadParam_RT *qRT = xBS_ThreadParam_RT+qi; qRT->rProc = NULL; SetEvent( qRT->rEvent[0] );
    #if 1 // note that in context of DLLmain opened threads might never signal 
      Sleep( 1 );
    #else 
      if (AON<0 || WaitForSingleObject(qRT->rThHn,50)==WAIT_TIMEOUT) TerminateThread( qRT->rThHn, 2 );
    #endif
    CloseHandle( qRT->rThHn ); CloseHandle( qRT->rEvent[0] ); CloseHandle( qRT->rEvent[1] );
  }; xBS_ThreadPool = -1; } //for i

static void iSpawnThreadsParams( int ANThreads, int AFst, int ALst, tBS_SpawnThreadProc AThreadProc,  
  const void *AP0, const void *AP1, const void *AP2, const void *AP3, tBS_ThreadParam_RT *ART ) { ALst -= AFst; int qSpawn = 0;
  do { ART->rR.rP[0] = (void*)AP0; ART->rR.rP[1] = (void*)AP1; ART->rR.rP[2] = (void*)AP2; ART->rR.rP[3] = (void*)AP3; ART->rProc = AThreadProc; ART->rR.rSpawn = qSpawn++;
       ART->rR.rFst = AFst; int qSpan = ALst/ANThreads--; AFst += qSpan; ART->rR.rLst = AFst; ALst -= qSpan; ART++; } while (ALst); } // do
void BS_SpawnThreads( int ANThreads, int AFst, int ALst, tBS_SpawnThreadProc AThreadProc, const void *AP0, const void *AP1, const void *AP2, const void *AP3 ) {
  const bool cSpawnInStack = (ANThreads&cBS_SpawnPool)==0; ANThreads &= 0xFF;
  if (ANThreads>ALst-AFst)        ANThreads = ALst-AFst; else
  if (ANThreads<1)                ANThreads = xBS_SpawnThreads; else
  if (ANThreads>cBS_SpawnThreads) ANThreads = cBS_SpawnThreads;
  if (ANThreads==1) { tBS_SpawnThreadRec qR; qR.rP[0] = (void*)AP0; qR.rP[1] = (void*)AP1; qR.rP[2] = (void*)AP2; qR.rP[3] = (void*)AP3;
    qR.rSpawn = 0; qR.rFst = AFst; qR.rLst = ALst; AThreadProc( &qR ); return; }; // if
  if (cSpawnInStack) { // portable spawn in stack
    tBS_ThreadParam_RT qRTs[cBS_SpawnThreads]; iSpawnThreadsParams( ANThreads, AFst, ALst, AThreadProc, AP0, AP1, AP2, AP3, qRTs );
    HANDLE qThHn[cBS_SpawnThreads]; for (int qi1=ANThreads-1;qi1>=0;qi1--) qThHn[qi1] = CreateThread( NULL, 0, iSpawnedThreadPlain, &qRTs[qi1].rR, 0, NULL );
    WaitForMultipleObjects( ANThreads, qThHn, TRUE, INFINITE ); for (int qi2=0;qi2<ANThreads;qi2++) CloseHandle( qThHn[qi2] ); return;
  }; // if SpawnInStack
  BS_Mutex( 1, &xBS_ThreadCS ); while (xBS_ThreadPool<ANThreads) { // create threads' pool
    tBS_ThreadParam_RT *qRT = xBS_ThreadParam_RT+xBS_ThreadPool; qRT->rEvent[0] = CreateEvent( NULL, FALSE, FALSE, NULL );
    xBS_ThreadEvWait[xBS_ThreadPool] = qRT->rEvent[1] = CreateEvent( NULL, FALSE, FALSE, NULL );
    DWORD qThID; qRT->rThHn = CreateThread( NULL, 0, iBS_SpawnedThreadProc, qRT, 0, &qThID ); xBS_ThreadPool++; }; // while
  iSpawnThreadsParams( ANThreads, AFst, ALst, AThreadProc, AP0, AP1, AP2, AP3, xBS_ThreadParam_RT );
  for (int qi=0;qi<ANThreads;qi++) { tBS_ThreadParam_RT *qRT = xBS_ThreadParam_RT+qi; ResetEvent( qRT->rEvent[1] ); SetEvent( qRT->rEvent[0] ); };
  WaitForMultipleObjects( ANThreads, xBS_ThreadEvWait, TRUE, INFINITE ); BS_Mutex( -1, &xBS_ThreadCS ); }

static struct { tBS_Spawn1Proc rProc; HANDLE rThHn; } ixSPNs[3]; static DWORD WINAPI iSpawn1Proc( void *AParam ) { ixSPNs[mcVi(AParam)].rProc( 1 ); return 0; };
void BS_Spawn1( tBS_Spawn1Proc AProc, int ASlot ) { AProc( 0 ); ixSPNs[ASlot].rProc = AProc;
  DWORD qThID; ixSPNs[ASlot].rThHn = CreateThread( NULL, 0, iSpawn1Proc, mcIv(ASlot), 0, &qThID ); }
void BS_Spawn1Wait( int ASlot ) { int qTicks = BS_Ticks(); WaitForSingleObject( ixSPNs[ASlot].rThHn, INFINITE );
  char qS[88]; sprintf( qS, "inf Spawn1Wait(%d) %d<msec>", ASlot, BS_Ticks()-qTicks ); WriteLn( qS );
  CloseHandle( ixSPNs[ASlot].rThHn ); ixSPNs[ASlot].rProc( 2 ); }

//================= Commons ==================
char *BS_AppEnd( char *APath, const char *ASuffix ) {
  char *qPathEnd = strend( APath ); return qPathEnd+sprintf( qPathEnd, "\\%s", ASuffix ); }

bool BS_DirExists( const char *ADir ) { return (int)GetFileAttributes( ADir )!=-1; }
bool BS_DirCreate( const char *ADir ) { if (BS_DirExists( ADir )) return true; CreateDirectory( ADir, NULL ); return BS_DirExists( ADir ); }

char xBS_DirHome[192] = ".\\", xBS_Dir[40] = ""; bool BS_Dir( char *ADir ) {
  if (ADir==xBS_DirHome) { return GetCurrentDirectory( mcASz(xBS_DirHome), xBS_DirHome )>0; };
  if (!ADir) { xBS_Dir[0] = 0; return SetCurrentDirectory( xBS_DirHome )!=0; };
  char qS[256]; sprintf( qS, "%s\\%s", xBS_DirHome, ADir );
  BS_DirCreate( qS ); if (SetCurrentDirectory( qS )) { strncpy( xBS_Dir, ADir, mcASz(xBS_Dir) ); return true; };
  sprintf( qS, "BS_Dir(%s\\%s)", xBS_DirHome, ADir ); ErrMsg( qS ); return false; }

char xBS_LogPath[192], *xBS_LogName = NULL;
int BS_Write( const char *AFName, int AN, const tBS_Val *AData ) { if (!AFName || !*AFName) return 0; // -- dump frame
  if (xBS_LogName && AFName[0]!='\\' && AFName[0]!='/' && AFName[1]!=':') {
    strcpy( xBS_LogName, AFName ); AFName = xBS_LogPath; };
  FILE *qFO = fopen( AFName, "wb" );
  if (!qFO) { char qS[256]; sprintf( qS, "BS_Write(%s):=%d", AFName, errno ); ErrMsg( qS ); return -1; }; 
  fwrite( AData, cBS_Sz, AN, qFO ); fclose( qFO ); return 1; }

int BS_Read( const char *AFName, int ASizeBytes, tBS_Val *AData ) { if (!AFName || !*AFName) return 0; // -- read frame
  FILE *qFO = fopen( AFName, "rb" ); if (!qFO) { char qS[88];
    if (ASizeBytes<=0) { sprintf( qS, "WRN:BS_Read(%.44s):=%d", AFName, errno ); WriteLn( qS ); }
                  else { sprintf( qS, "BS_Read(%.44s):=%d", AFName, errno ); ErrMsg( qS ); };
    return -1; }
  int qRes = (int)fread( AData, 1, abs(ASizeBytes), qFO );  fclose( qFO ); return qRes; }

void BS_DatasWrt( const char *AFPrefix, int ANP, tBS_Val **ADatas, int AI1, int AI2 ) {
  char qS[256]; const char cC[25] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','I','J','K','L','M','N','O','P'};
  int qFlags = 0xFFFFFFFF; if (!AI2) { qFlags = AI1; AI1 = 0; AI2 = mcASz(cC); };
  for (int qi=AI1;qi<AI2;qi++) if (ADatas[qi] && (1<<qi)&qFlags) { sprintf( qS, "%s%c.raw", AFPrefix, cC[qi] ); BS_Write( qS, ANP, ADatas[qi] ); }; }
void BS_DatasWrt( const char *AFPrefix, int ANP, tBS_Val *A0, tBS_Val *A1, tBS_Val *A2, tBS_Val *A3, tBS_Val *A4, tBS_Val *A5 ) {
  int qFlg = 0; tBS_Val *qDs[6] = {A0,A1,A2,A3,A4,A5}; for (int qi=0;qi<mcASz(qDs);qi++) if (qDs[qi]) qFlg |= 1<<qi; BS_DatasWrt( AFPrefix, ANP, qDs, qFlg, 0 ); }

void BS_DatasAlloc( int ANDD, int AN, tBS_Val **ADatas ) { for (int qi=0;qi<ANDD;qi++) ADatas[qi] = new tBS_Val[AN]; } // -- alloc frames
void BS_DatasRelease( int ANDD, tBS_Val **ADatas ) { while (--ANDD>=0) if (ADatas[ANDD]) delete ADatas[ANDD]; } // -- dealloc frames

unsigned BS_Ticks() { return (unsigned)GetTickCount64(); }
void BS_Sleep( int AMsec ) { Sleep( AMsec ); }
int BS_GetCh() { return kbhit() ? getch() : 0; }
