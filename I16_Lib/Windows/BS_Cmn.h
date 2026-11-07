// BS_Cmn.h : general-purpose utilities (OS,EABI-specific)
// for Windows
#ifndef BS_CMN_H_INCLUDED
#define BS_CMN_H_INCLUDED

//=========== Threading service ==============
extern int xBS_SpawnThreads; void BS_SpawnThreads_State( bool AOn ); // #of threads; alloc,dealloc threads' pool
typedef struct{void *rP[4]; int rSpawn,rFst,rLst;} tBS_SpawnThreadRec; typedef void(*tBS_SpawnThreadProc)(tBS_SpawnThreadRec*);
const int cBS_SpawnPool = 1<<16; // use narrow threads' pool (when not ported then this flag is simply igored)
void BS_SpawnThreads( int ANThreads, int AFst, int ALst, tBS_SpawnThreadProc AThreadProc, const void *AP0, const void *AP1=NULL, const void *AP2=NULL, const void *AP3=NULL );

typedef void (*tBS_Spawn1Proc)(int); // Spawn auxiliary operation; 0:init_copy, 1:operate, 2:final_copy
void BS_Spawn1( tBS_Spawn1Proc AProc, int ASlot=0 );
void BS_Spawn1Wait( int ASlot=0 );

//================= Commons ==================
#define cPSChr '\\' // path-separator character
#define cPSStr "\\" // path-separator string
char *BS_AppEnd( char *APath, const char *ASuffix );
bool BS_DirExists( const char *ADir ); // whether a directory exists
bool BS_DirCreate( const char *ADir ); // creates a directory
extern char xBS_DirHome[], xBS_Dir[];
bool BS_Dir( char *ADir ); // 'Dir=_DirHome':GetCurDir(_DirHome), '=nil':SetCurDir(_DirHome), '<>nil':SetCurDir(Dir)
extern char xBS_LogPath[], *xBS_LogName;
int BS_Write( const char *AFName, int AN, const tBS_Val *AData ); // dumps frame
int BS_Read( const char *AFName, int ASizeBytes, tBS_Val *AData ); // reads frame
void BS_DatasWrt( const char *AFPrefix, int ANP, tBS_Val **ADatas, int AI1=0, int AI2=8 ); // dump frames
void BS_DatasWrt( const char *AFPrefix, int ANP, tBS_Val *A0, tBS_Val *A1=NULL, tBS_Val *A2=NULL, tBS_Val *A3=NULL, tBS_Val *A4=NULL, tBS_Val *A5=NULL );
unsigned BS_Ticks(); // <msec> since machine starup
void BS_Sleep( int AMsec );
int BS_GetCh();

unsigned inline BS_Align( int AV ) { const int cBS_LnMaxLenA = 4; return (AV+cBS_LnMaxLenA-1)&~(cBS_LnMaxLenA-1); } // inc+align (for dir-iterator's buffer)
void inline BS_Cpy( int AN, const tBS_Val *ASrc, tBS_Val *ADst ) { if (AN>0 && ADst!=ASrc) memcpy( ADst, ASrc, AN*cBS_Sz ); }
void inline BS_RowsCpy( int AWSrc, const tBS_Val *ASrc, int AWDst, tBS_Val *ADst, int ARW, int ARH ) { for (int qi=0;qi<ARH;qi++) { BS_Cpy( ARW, ASrc, ADst ); ASrc += AWSrc; ADst += AWDst; }; }
void inline BS_DataDup( int ANDD, int AN, const tBS_Val *ASrc, tBS_Val **ADatasDst ) { for (int qi=0;qi<ANDD;qi++) BS_Cpy( AN, ASrc, ADatasDst[qi] ); }; // replica
void BS_DatasAlloc( int ANDD, int AN, tBS_Val **ADatas ); // allocs frames
void BS_DatasRelease( int ANDD, tBS_Val **ADatas ); // deallocs frames

#endif // BS_CMN_H_INCLUDED
