// BS_Cmn.cpp : general-purpose utilities (OS,EABI-specific)
// for Linux by S.B.
#include "stdafx.h"
#include <sys/time.h>
#include <stdio.h>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h> // 'pthread' library shall be linked
#include "BS_Cmn.h"

//========== for Windows portability  ============
void BS_Mutex( int AOp, pthread_mutex_t *ACS ) { int qErr = 0;  switch (AOp) { // -- 2:create, 1:enter, -1:leave, -2:delete
  case 2: qErr = pthread_mutex_init( ACS, NULL ); break; case -2: qErr = pthread_mutex_destroy( ACS ); break;
  case 1: qErr = pthread_mutex_lock( ACS ); break; case -1: qErr = pthread_mutex_unlock( ACS ); break; };
  if (qErr) { char qS[40]; sprintf( qS, "BS_Mutex(%d):=%d", AOp, qErr ); ErrMsg( qS ); }; }
typedef struct { pthread_mutex_t m; pthread_cond_t c; bool s;} tBS_HEvent;
void CreateEvent( tBS_HEvent *e ) { BS_Mutex( 2, &e->m); pthread_cond_init( &e->c, NULL ); e->s = false; }
void SetEvent( tBS_HEvent *e ) { BS_Mutex( 1, &e->m ); e->s = true; pthread_cond_signal( &e->c ); BS_Mutex( -1, &e->m ); }
void WaitEvent( tBS_HEvent *e ) { BS_Mutex( 1, &e->m ); while (!e->s ) pthread_cond_wait( &e->c, &e->m ); e->s = false; BS_Mutex( -1, &e->m ); }
void DestroyEvent( tBS_HEvent *e ) { BS_Mutex( -2, &e->m ); pthread_cond_destroy( &e->c ); }

//========== Threading service ============
const int cBS_SpawnThreads = 16; int xBS_SpawnThreads = 4;
typedef struct { tBS_SpawnThreadRec rR; tBS_SpawnThreadProc rProc; pthread_t rThHn; tBS_HEvent rEvent[2]; } tBS_ThreadParam_RT;
void *iSpawnedThreadPlain( void *AParam ) { tBS_ThreadParam_RT *qP = (tBS_ThreadParam_RT*)AParam; qP->rProc( &qP->rR ); pthread_exit( NULL ); }
void *iBS_SpawnedThreadProc( void *AParam ) { tBS_ThreadParam_RT *qParam = (tBS_ThreadParam_RT*)AParam;
  do { WaitEvent( &qParam->rEvent[0] ); if (qParam->rProc) try { qParam->rProc( &qParam->rR ); } catch(...) {};
  SetEvent( &qParam->rEvent[1] ); } while (qParam->rProc); pthread_exit( NULL ); } // do

static tBS_ThreadParam_RT  xBS_ThreadParam_RT[cBS_SpawnThreads];
static pthread_mutex_t     xBS_ThreadCS; // critical section 
static int                 xBS_ThreadPool = -1; void BS_SpawnThreads_State( bool AOn ) { // -- alloc,dealloc threads' pool
  if (AOn) { BS_Mutex( 2, &xBS_ThreadCS ); xBS_ThreadPool = 0; return; }; // initialize threads' pool
  if (xBS_ThreadPool<0) return; BS_Mutex( -2, &xBS_ThreadCS ); // deinitialize threads' pool 
  for (int qi=0;qi<xBS_ThreadPool;qi++) {
    tBS_ThreadParam_RT *qRT = xBS_ThreadParam_RT+qi; qRT->rProc = NULL; SetEvent( &qRT->rEvent[0] ); pthread_join( qRT->rThHn, NULL );
    DestroyEvent( &qRT->rEvent[0] ); DestroyEvent( &qRT->rEvent[1] );
  }; xBS_ThreadPool = -1; } //for i
static void iSpawnThreadsParams( int ANThreads, int AFst, int ALst, tBS_SpawnThreadProc AThreadProc,
  const void *AP0, const void *AP1, const void *AP2, const void *AP3, tBS_ThreadParam_RT *ART ) { ALst -= AFst; int qSpawn = 0;
  do { ART->rR.rP[0] = (void*)AP0; ART->rR.rP[1] = (void*)AP1; ART->rR.rP[2] = (void*)AP2; ART->rR.rP[3] = (void*)AP3; ART->rProc = AThreadProc;
       ART->rR.rSpawn = qSpawn++; ART->rR.rFst = AFst; int qSpan = ALst/ANThreads--; AFst += qSpan;
       ART->rR.rLst = AFst; ALst -= qSpan; ART++; } while (ALst); } // do

void BS_SpawnThreads( int ANThreads, int AFst, int ALst, tBS_SpawnThreadProc AThreadProc, const void *AP0, const void *AP1, const void *AP2, const void *AP3 ) {
  const bool cSpawnInStack = (ANThreads&cBS_SpawnPool)==0; ANThreads &= 0xFF;
  if (ANThreads>ALst-AFst)        ANThreads = ALst-AFst; else
  if (ANThreads<1)                ANThreads = xBS_SpawnThreads; else
  if (ANThreads>cBS_SpawnThreads) ANThreads = cBS_SpawnThreads;
  if (ANThreads==1) { tBS_SpawnThreadRec qR; qR.rP[0] = (void*)AP0; qR.rP[1] = (void*)AP1; qR.rP[2] = (void*)AP2;
     qR.rP[3] = (void*)AP3; qR.rSpawn = 0; qR.rFst = AFst; qR.rLst = ALst; AThreadProc( &qR ); return; }; // if
  if (cSpawnInStack) { tBS_ThreadParam_RT qRTs[cBS_SpawnThreads]; // portable spawn in stack
    iSpawnThreadsParams( ANThreads, AFst, ALst, AThreadProc, AP0, AP1, AP2, AP3, qRTs );
    for (int qi1=ANThreads-1;qi1>=0;qi1--) { pthread_create( &qRTs[qi1].rThHn, NULL, iSpawnedThreadPlain, &qRTs[qi1] );
    pthread_join( qRTs[qi1].rThHn, NULL ); }; return; };
  BS_Mutex( 1, &xBS_ThreadCS );
  while (xBS_ThreadPool<ANThreads) { // create threads' pool
    tBS_ThreadParam_RT *qRT = xBS_ThreadParam_RT+xBS_ThreadPool;
    CreateEvent( &qRT->rEvent[0] ); CreateEvent( &qRT->rEvent[1] );
    pthread_create( &qRT->rThHn, NULL, iBS_SpawnedThreadProc, qRT );
    xBS_ThreadPool++; }; // while
  tBS_ThreadParam_RT *qRTs = xBS_ThreadParam_RT;
  iSpawnThreadsParams( ANThreads, AFst, ALst, AThreadProc, AP0, AP1, AP2, AP3, qRTs );
  for (int qi1=ANThreads-1;qi1>=0;qi1--) SetEvent( &qRTs[qi1].rEvent[0] ); sched_yield();
  for (int qi2=ANThreads-1;qi2>=0;qi2--) WaitEvent( &qRTs[qi2].rEvent[1] );
  BS_Mutex( -1, &xBS_ThreadCS ); }

static struct { tBS_Spawn1Proc rProc; pthread_t rThHn; } ixSPNs[3];
static void *iSpawn1Proc( void *AParam ) { ixSPNs[mcVi(AParam)].rProc( 1 ); return NULL; };
void BS_Spawn1( tBS_Spawn1Proc AProc, int ASlot ) { AProc( 0 ); ixSPNs[ASlot].rProc = AProc;
  pthread_create( &ixSPNs[ASlot].rThHn, NULL, iSpawn1Proc, mcIv(ASlot) ); }
void BS_Spawn1Wait( int ASlot ) { unsigned qTicks = BS_Ticks(); pthread_join( ixSPNs[ASlot].rThHn, NULL );
  char qS[88]; sprintf( qS, "inf Spawn1Wait(%d) %d<msec>", ASlot, BS_Ticks()-qTicks ); WriteLn( qS );
  ixSPNs[ASlot].rProc( 2 ); }

//================= Commons ==================
char *BS_AppEnd( char *APath, const char *ASuffix ) {
  char *qPathEnd = strend( APath ); return qPathEnd+sprintf( qPathEnd, "/%s", ASuffix ); }

bool BS_DirExists( const char *ADir ) { if (*ADir==0) return true;
  struct stat qStt; return stat(ADir,&qStt)!=-1 && S_ISDIR( qStt.st_mode ); }
bool BS_DirCreate( const char *ADir ) { if (BS_DirExists( ADir )) return true;
  mkdir( ADir, 0777 ); return BS_DirExists( ADir ); }

char xBS_DirHome[192] = "./", xBS_Dir[40] = ""; bool BS_Dir( char *ADir ) {
  if (ADir==xBS_DirHome) return getcwd( xBS_DirHome, mcASz(xBS_DirHome) )!=NULL;
  if (!ADir) { xBS_Dir[0] = 0; return !chdir( xBS_DirHome ); };
  char qS[256]; sprintf( qS, "%s/%s", xBS_DirHome, ADir );
  BS_DirCreate( qS ); if (!chdir( qS )) { strncpy( xBS_Dir, ADir, mcASz(xBS_Dir)-1 ); return true; };
  sprintf( qS, "BS_Dir(%s/%s)", xBS_DirHome, ADir ); ErrMsg( qS ); return false; }

char xBS_LogPath[192], *xBS_LogName = NULL;
int BS_Write( const char *AFName, int AN, const tBS_Val *AData ) { if (!AFName || !*AFName) return 0; // -- dump frame
  if (xBS_LogName && AFName[0]!='~' && AFName[0]!='/') {
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

unsigned BS_Ticks() { struct timespec qNow; clock_gettime( CLOCK_BOOTTIME, &qNow );
  return (int)(qNow.tv_sec*1000L+qNow.tv_nsec/1000000L); }

void BS_Sleep( int AMsec ) { usleep( AMsec*1000 ); }

int BS_GetCh() {
  struct termios qOldt, qNewt; tcgetattr( STDIN_FILENO, &qOldt ); qNewt = qOldt;
  qNewt.c_lflag &= ~(ICANON|ECHO); tcsetattr( STDIN_FILENO, TCSANOW, &qNewt );
  int qOldf = fcntl( STDIN_FILENO, F_GETFL, 0 ); fcntl( STDIN_FILENO, F_SETFL, qOldf|O_NONBLOCK );
  char qCh = getchar(); tcsetattr( STDIN_FILENO, TCSANOW, &qOldt );
  fcntl( STDIN_FILENO, F_SETFL, qOldf ); if (qCh==EOF) qCh=0; return qCh; }
