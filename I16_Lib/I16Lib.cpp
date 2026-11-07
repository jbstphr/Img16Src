// I16Lib.cpp : Vatech Fluoro Processing
// by S.B.
// source directory 'I16_GenC12' : general process (XML created for CUDA12.2)
// note1: for backward compatibility (~CUDA10) replace ~"12.2" in vcxproj to ~"10.2"
// note2: specific CUDA <might> needs specific toolchain; e.g. 
//        CUDA10 might ~never have targets|props.xml for MS-toolchain ver143,
//        CUDA12 might ~never install targets|props.xml for MS-toolchain ver140.

#include "stdafx.h"
#include "math.h"
#include "BS_GIni.h" // A wrapper around INI-files' (or registry's) API
#include "BS_Cmn.h"  // common routines (OS,EABI-specific)
#include "BS_Ut0.h"  // CPU basic-routines
#include "BS_Ut1.h"  // CUDA basic-routines (~interconnect with kernel.cu)
#include "BS_Ut2.h"  // CPU aplication-routines with standard dimension
#include "BS_Ut2a.h" // CPU aplication-routines addition
#include "BS_Ut3.h"  // CUDA mid-level routines
#include "BS_Ut4.h"  // beam-hardening, CLAHE (mostly CPU|CUDA is a parameter) 
#include "BS_Ut5.h"  // unsharpen's core for CPU|CUDA 
#include "BS_Ut5a.h" // air-estimation
#include "BS_Ut5c.h" // quad-laps
#include "BS_Ut5d.h" // extended mask-recount
#include "BS_Ut6.h"  // despeckle, thin_laps (CPU & CUDA with explicit dimension)
#include "BS_Ut7.h"  // ad-hoc service CPU-transforms
#include "BS_Ut7a.h" // collimator-estimation
#include "BS_Ut7b.h" // ribs-detection

//// Auxiliary parameters ////
int         xPars[2], // array of ints in the 'data-segment' 
            xN = 0, xModeId = 0; // reserved buffer size and processing mode

void I16_ParSet( int AParIdx, int AParVal ) { // set auxilary int parameters
  if (AParIdx<0 || AParIdx>=mcASz(xPars)) return; xPars[AParIdx] = AParVal; }
int I16_ParGet( int AParIdx ) { // get auxilary int parameters
  if (AParIdx<0 || AParIdx>=mcASz(xPars)) return 0x80000000; return xPars[AParIdx]; }
char *I16_ParStr( int AParIdx ) { return NULL; }; // stub

//// Debugs ////
void WriteLn( const char *AMsg ) { fprintf( stderr, "%s\n", AMsg ); } // -- WriteLn() to ~a console
void ErrMsg( const char *AMsg ) { // -- WriteLn( "ERR " Msg )
  char qMsg[256]; sprintf( qMsg, cTPS" ! %.222s", AMsg ); WriteLn( qMsg ); }

int         xDbgDmp = 0; // 0:minimalistic trace, 1:verbose trace, 2:dump images
void I16_DbgDmpReset() { // -- resets debug-dump
  BS_Dir( xBS_DirHome ); // remember home-directory
  try {
    fprintf( stderr, cTPS" DbgDmpReset purging older output :=%d]n", system( 
      #ifdef _WIN32
        "for /f %i in ('dir /a:d /b $*') do rd /s /q %i") );
      #else
        "find . -mindepth 1 -maxdepth 1 -type d -name '$*' -exec rm -rf {} +") );
      #endif
  } catch(...) { fprintf( stderr, cTPS" ! purging exception" ); }; }

void I16_DbgDmp( const char *ASubDir ) { // nil:minimalistic trace, "":verbose trace, "xyz":dump images
  // <>nil:interpret parameters creatively as a prefix of files' path; basically as 'separate subdirctory'
  if (ASubDir) { char qS[192]; strcpy( qS, ASubDir[0]=='$'?"":"$" );
    if (ASubDir[0]) { xDbgDmp = 2; strcat( qS, ASubDir ); 
      #if 0 // remove suffix ~"_0888" and dump debug ~"$00/*.raw" 
        char *qCh = strrchr( qS, '_' ); if (qCh && (size_t)(qCh-qS)>2) *qCh = 0;
      #endif
      BS_Dir( qS ); sprintf( qS, cTPS" DbgDmp to '%s'", ASubDir ); WriteLn( qS );
    } else
      xDbgDmp = 1;
  } else { xDbgDmp = 0; BS_Dir( NULL ); WriteLn( cTPS" DbgDmp off" ); }; }

//// application-functionality:statistical suppress of vertical laplacians ////
static void jTLLLoad( // load tThLoLap from initialization profile at current SectName() 
  GIniFile *AF, char *AKey, tThLoLap *ATLL, bool AComment ) { // ini-file, section, control record
  const char *cS = AF->DGetStr( AKey, NULL ); int qI[3]; float qF[1];
  int qN = cS ? sscanf(cS,"%d %d %d %f",qI+0,qI+1,qI+2,qF+0): 0, qMode = ATLL->n[0];
  if (qN>0) { ATLL->n[0] = qMode = qI[0]; if (qN>3) { ATLL->n[1] = qI[1]; ATLL->n[2] = qI[2]; ATLL->f[0] = qF[0]; }; };
  char qS[88] = {0}; if (!qMode) { if (!qN) strcpy( qS, "0" ); } else {
    if (qN<=3) sprintf( qS, "%d %d %d %.2f", qMode, ATLL->n[1], ATLL->n[2], ATLL->f[0] ); };
  if (qS[0]) {
    if (qMode==1) strcat( qS, "; SetFixScalar CGN0 CGN1 Att" ); else
    if (AComment)  strcat( qS, "; Mode RadX Y Att" );
    AF->SetStr( AKey, qS );
  }; }
static int jTLLLabel( char *AStr, tThLoLap *ATLL, int ALabelIdx ) { // debug-label's text about tThLoLap
  return ATLL->n[0] ? sprintf(AStr," %c=%d %d %d %.1f",'a'+ALabelIdx,ATLL->n[0],ATLL->n[1],ATLL->n[2],ATLL->f[0]) : 0; }

//// application-functionality: gray-scale-cut ////
float iiGammaByInp( float AGamma ) { return AGamma<1e-3f?1.0f:mcRng(AGamma,0.05f,9.99f); }
typedef struct {bool rEnabled,rInvert; int rGam[2]; float rCut[2],rGamma;} tGrayCut; // control record
static void jGrayCutLoad( GIniFile *AF, tGrayCut *AGC ) {
  int qI[16]; float qF[16]; char qS[256];
  if (Gini_sscanf( AF->DGetStr("nEnabled",NULL), qI )>0) AGC->rEnabled = qI[0]>0; else {
    sprintf( qS, "%d", (int)AGC->rEnabled ); AF->SetStr( "nEnabled", qS ); };
  if (!AGC->rEnabled) return;
  if (Gini_sscanf( AF->DGetStr("nInvert",NULL), qI )>0) AGC->rInvert = qI[0]>0; else {
    sprintf( qS, "%d", (int)AGC->rInvert ); AF->SetStr( "nInvert", qS ); };
  if (Gini_sscanf( AF->DGetStr("nGamInpOut",NULL), qI )>1) { AGC->rGam[0] = qI[0]; AGC->rGam[1] = qI[1]; } else {
    sprintf( qS, "%d %d", AGC->rGam[0], AGC->rGam[1] ); AF->SetStr( "nGamInpOut", qS ); };
  if (Gini_sscanf( AF->DGetStr("fCutDrkBrt",NULL), qF )>1) { AGC->rCut[0] = qF[0]; AGC->rCut[1] = qF[1]; } else {
    sprintf( qS, "%.4f %.4f", AGC->rCut[0], AGC->rCut[1] ); AF->SetStr( "fCutDrkBrt", qS ); };
  if (Gini_sscanf( AF->DGetStr("fGamma",NULL), qF )>0) { AGC->rGamma = iiGammaByInp(qF[0]); } else {
    sprintf( qS, "%.3f", AGC->rGamma ); AF->SetStr( "fGamma", qS ); }; }

static void jGrayCut( // GrayCut apply
  const tGrayCut *AGC, int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs, const char *ADbgPfx ) { 
  if (!AGC->rEnabled) return;
  char *lDbgSfx = NULL, qDbgPfx[88]; if (ADbgPfx) { strcpy( qDbgPfx, ADbgPfx ); lDbgSfx = strend(qDbgPfx); };
  const int cN = AW*AH; BS2_Dim( AW, AH );
  //int   qInpWhite = AGC->rGam[0], qInpShift = 0; // fit occasional unsigned inp gamut into rational diapazon
  //while (qInpWhite>1<<14) { qInpWhite >>= 1; qInpShift++; }; if (qInpShift) mcShiftRgt<tBS_Card>( cN, (tBS_Card*)AFrm, qInpShift );
  tBS2_HistCut qCut = { {0,AGC->rGam[0],AGC->rGam[1],0}, {AGC->rCut[0],AGC->rCut[1],0,0} };
  if (AGC->rInvert) mcSwap( int, qCut.n[cBS2_rhcBlack], qCut.n[cBS2_rhcWhite] );
  BS2_HistCut( &qCut, AFrm, NULL, ATs );
  if (1 && lDbgSfx) { strcpy( lDbgSfx, "bCut.raw" ); BS_Write( qDbgPfx, cN, AFrm ); };
  if (fabsf(AGC->rGamma-1)>1e-3f) {
    BS2_Gamma( AGC->rGamma, AGC->rGam[1], AFrm );
    if (1 && lDbgSfx) { strcpy( lDbgSfx, "cGamma.raw" ); BS_Write( qDbgPfx, cN, AFrm ); };
  }; }

static int jGrayCutLabel( char *AStr, const tGrayCut *AGC ) { // debug-label's text about gray-cut
  char *qCh = AStr+sprintf( AStr, "GC " );
  if (AGC->rEnabled) { qCh += sprintf( qCh, "%d=>%d", AGC->rGam[0], AGC->rGam[1] );
    if (AGC->rCut[0]>1e-6f || AGC->rCut[1]>1e-6f) qCh += sprintf( qCh, " %.3f,%.3f", AGC->rCut[0], AGC->rCut[1] );
    if (fabs(AGC->rGamma)-1>1e-3f) qCh += sprintf( qCh, " %.1f", AGC->rGamma );
    if (AGC->rInvert) qCh += sprintf( qCh, " inv" );
  } else qCh += sprintf( qCh, "disabled" );
  return (int)(qCh-AStr); }

//// optional application-functionality: collimator-estimation ////
static tBS_Val *xCollMsk;
static int jCollPreProcess( // collimator preprocess; returns # of occupied entries in Ts
  int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs, const char *ADbgPfx ) { 
  if (BS7a_CollimatorEnabled()<1) { xCollMsk = NULL; return 0; }; xCollMsk = ATs[0]; // keep entry Ts#0 for postprocess
  if (!BS7a_CollimatorPreProcess( AW, AH, AFrm, xCollMsk, ATs+1, ADbgPfx )) { xCollMsk = NULL; return 0; };
  return 1; } 
static void jCollPostProcess( // collimator postprocess
  int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs, const char *ADbgPfx ) { 
  if (xCollMsk) BS7a_CollimatorPostProcess( AW, AH, AFrm, xCollMsk, ATs+1, ADbgPfx ); };

//// optional application-functionality: ribs-detection ////
static int jRibsPreProcess( // ribs preprocess; returns # of occupied entries in Ts
  int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs, const char *ADbgPfx ) { 
  if (BS7b_RibsDetectionEnabled()<1) return 0; 
  BS7b_RibsDetection( AW, AH, AFrm, ATs, ADbgPfx );
  return 0; } 
static void jRibsPostProcess( // ribs postprocess
  int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs, const char *ADbgPfx ) { 
  if (BS7b_RibsDetectionEnabled()<1) return; };

//// auxiliary application-functionality: unsharpen ////
static void jjUSkey( GIniFile *AF, int AN, int *ADat, // - load integer key
  const char *AKey, const char *ACmt ) { 
  int qI[16]; memset( qI, 0, sizeof(qI) ); 
  if (xPars[cParIniSetDflt]<=0 && Gini_sscanf( AF->DGetStr(AKey,NULL), qI )>0) {
    for (int qn=0;qn<AN;qn++) ADat[qn] = qI[qn]; return; };
  char qS[256], *qCh = qS; for (int qn=0;qn<AN;qn++) { if (qn) *qCh++ = ' ';
    qCh += sprintf( qCh, "%d", ADat[qn] ); }; // for i
  if (ACmt) qCh += sprintf( qCh, "; %s", ACmt ); AF->SetStr( AKey, qS ); }
static void jjUSkey( GIniFile *AF, int AN, float *ADat, // - load float key
  const char *AKey, const char *ACmt, const int *APrecs=NULL ) {
  float qF[16]; memset( qF, 0, sizeof(qF) ); 
  if (xPars[cParIniSetDflt]<=0 && Gini_sscanf( AF->DGetStr(AKey,NULL), qF )>0) {
    for (int qn=0;qn<AN;qn++) ADat[qn] = qF[qn]; return; };
  char qS[256], *qCh = qS; for (int qn=0;qn<AN;qn++) { if (qn) *qCh++ = ' '; float qV = ADat[qn]; int qPrec = APrecs?APrecs[qn]:2;
    qCh += fabs(qV)<1e-9f || fabs(qV)>=(float)pow(10.0f,-qPrec) ? sprintf(qCh,"%.*f",qPrec,qV) : sprintf(qCh,"%g",qV); }; // for i
  if (ACmt) qCh += sprintf( qCh, "; %s", ACmt );  AF->SetStr( AKey, qS ); }
static void jjadj( tUnSharpAdjs *AUSAdjs, const void *AUS, // - set adjustment parameters
  int ABSCi, bool AIsFloat, float AKmn, float AKmx, const void *AAddrInUS ) {
  // BSCindex #0:"Brightness", #1:"Contrast", #2:"Sharpness"
  tUnSharpAEn *qE = AUSAdjs->rE+AUSAdjs->rN++;
  qE->rBCSi = mcRng( ABSCi, 0, 2 );
  qE->rIsFloat = AIsFloat; qE->rKmn = AKmn; qE->rKmx = AKmx;
  qE->rAddrOffs = (int)(size_t)((char*)AAddrInUS-(char*)AUS); }

typedef struct { float rn[3],rx[3]; } tAdjMnMx; // adjusments of few UnSharp's fields
  // (for Tag=1,=15 respectivelly) #0:Brightness #1:Contrast #2:for Sharpness 
#include "I16Lib.inl" // -- Frequently changing declarations for jUnsharp..., I16_Init.
                      // -- "jUnSharpLoad(...)" and "static const tUnSharp cUS??".

static int jUnSharpLabel( char *AStr, const tUnSharp *AUS ) { // debug-label's text about unsharp
  char *qCh = AStr+sprintf( AStr, "US " );
  if (AUS->rEnabled) {
    qCh += sprintf( qCh, "%d,%d %.1f,%.1f,%.1f", AUS->rA[0], AUS->rA[1], AUS->fA[0], AUS->fA[1], AUS->fA[2] );
    qCh += sprintf( qCh, " g%.1f,%.1f r%d,%d C%d", AUS->fB[0], AUS->fB[1], AUS->rC[0], AUS->rC[1], AUS->rD[0] );
    qCh += sprintf( qCh, " M%d|%d,%.1f", AUS->rJ[0], AUS->rJ[1], AUS->fJ[0] );
  } else qCh += sprintf( qCh, "disabled" ); // gamma rad_Lap cl_Blk
  return (int)(qCh-AStr); }

static int jUnSharpExtsLabel( char *AStr, const tUnSharpExts *AUSExts ) { // debug-label's text about unsharp
  char *qCh = AStr+sprintf( AStr, "EX " ); // the 'UnSharpExts' reports run-time statistics
  qCh += sprintf( qCh, "Big=%.2f|%.2f CLAblr=%d", AUSExts->rSttAveBig, AUSExts->rSttNrmDelBig, AUSExts->rSttCLABlr );
  if (AUSExts->rAirInfFlg) qCh += sprintf( qCh," Air=%X", AUSExts->rAirInfFlg );
  return (int)(qCh-AStr); }

static void jUnSharp( // UnSharp apply
  const tUnSharp *AUS, const tUnSharpAdjs *AUSAdjs, int AAdj, tUnSharpExts *AUSExts,
  int AW, int AH, tBS_Val *AFrm, tBS_Val **ATsH, tBS_Val **ATsD, const char *ADbgPfx ) {
  if (!AUS->rEnabled) return;
  const int     cN = AW*AH; char *lDbgSfx = NULL, qDbgPfx[88]; // local debug prefixx
  if (ADbgPfx) { strcpy( qDbgPfx, ADbgPfx ); ADbgPfx = qDbgPfx; lDbgSfx = strend(qDbgPfx); };
  CBS_Ticks     qTM; if (xDbgDmp==1) strcpy( qTM.rPfx, "jUS" ); qTM.Reset(); 
  tUnSharp      qUS; BS5_UnSharpRecAdjust( AW, AH, AUS, AUSAdjs, xModeId, AAdj, &qUS, AUSExts );
  // - here is the place for high-level adhoc alterations of tUnSharp,tUnSharpExts 
  #if 0 // ????
    const float cShear0 = -0.07f;
    // BS2A_ShearCPU( cBS2A_FOrg, cShear0, AW, AH, AFrm, AW, AH, ATsH[0] );
    BS2A_Shear( cHerLikFo, cShear0, AFrm, ATsH[0], ATsD );
    if (lDbgSfx) { strcpy( lDbgSfx, "3xShear.raw" ); BS_Write( ADbgPfx, cN, AFrm ); *lDbgSfx = 0; };
  #endif
  // - remove Poisson's noise
  BS6_DePoisson( AW, AH, &qUS, AUSExts, AFrm, ATsH, ATsD, ADbgPfx );
  #if 0 // ????
    BS2A_Shear( cHerLikFo|cHerLikFs, -cShear0, AFrm, ATsH[0], ATsD );
  #endif
  qTM.Mark( "BS6_DePoisson" ); 
  // - marginal bevels
  const int     mgn_Bev = qUS.rMgn[1], mgn_Cor = qUS.rMgn[2]; if (mgn_Bev) {
    BS_BevCorners( AW, AH, AFrm, mgn_Cor, mgn_Bev, -1 );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "9Bevels.raw" ); BS_Write( ADbgPfx, cN, AFrm ); };
    qTM.Mark( "mgn_Bev" ); }; 
  // - invoke the main presentation enhancement
  #if 0 // ????
    const float cShear1 = -cShear0*1.5f;
    BS2A_Shear( cHerLikFo|cHerLikFs, cShear1, AFrm, ATsH[0], ATsD );
  #endif
  if (ATsD) BS5_UnSharp_CUDA( &qUS, AUSExts, AW, AH, AFrm, ATsH, ATsD, ADbgPfx, lDbgSfx ); // CUDA unsharp
       else BS5_UnSharp_CPU( &qUS, AUSExts, AW, AH, AFrm, ATsH, ADbgPfx, lDbgSfx ); // CPU unsharp
  #if 0 // ????
    BS2A_Shear( cHerLikFo|cHerLikFs, -cShear1, AFrm, ATsH[0], ATsD );
  #endif
  qTM.Mark( "BS5_UnSharp" );
  if (mgn_Bev) { // - marginal bevels
    BS_BevCorners( AW, AH, AFrm, mgn_Cor, mgn_Bev+1, 0 );  
    //const tBS2_HistCut cHC = { {0,20000,1<<12,0}, {0.001f,0.001f,0,0} }; BS2_HistCut( &cHC, AFrm, NULL, ATsH );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Y9Bevels.raw" ); BS_Write( ADbgPfx, cN, AFrm ); };
    qTM.Mark( "mgn_Bev" ); }; 
} // end of jUnSharp

//// Commons ////
#ifdef _WIN32
BOOL APIENTRY DllMain( HMODULE AModule, DWORD AReason, LPVOID AReserved ) { switch (AReason) {
  case DLL_PROCESS_ATTACH:
    memset( xPars, 0, sizeof(xPars) ); // initialize auxiliary parameters
    Gini_StartUps( AModule ); break; // CWD etc
  }; return TRUE; } // case DLL_THREAD_DETACH: case DLL_PROCESS_DETACH: break;
#else // ifLinux
  __attribute__((constructor)) void DllMainLoad(void) { // Equivalent to DLL_PROCESS_ATTACH
    memset( xPars, 0, sizeof(xPars) ); // initialize auxiliary parameters
    Gini_StartUps(); } // CWD etc
#endif

//// Application APIs ////
tBS_Val     *xDaH[48], // host-buffers
            *xDaD[48]; // device-buffers
tThLoLap     xTLL[4]= {
  {0, 44, 0, 0.26f},
  {0, 1, 15, 1.00f},
  {0, 22, 2, 0.33f},
  {0, 2, 22, 0.33f} };

tGrayCut     xGC = { false, true, {(1<<14)-1,(1<<12)-1}, {1e-3f,1e-3f}, 0 }; // Gamma=0 means don't apply
tUnSharp     xUS; // set by _Init
tUnSharpAdjs xUSAdjs; // set by _Init
int          xPrintLabel = 0, // 0:no, -1:dark, 1:bright
             xPyrMxX = 0, xPyrMxY = 0, // Mx dnsample (arith_rgt_shift; affects dbg_dumps' dims);
             xShiftRgt = 0; // shift inp pixels right to fit into signed diapazon
// refer "writeln Invoke ... as WxH" in I16_Invoke about actual dimension.
// Pre_processes alone (~'ButtingAreaCorrection via MedianBands') typically need 'PyrMxX=1' or
// '=2', while 'Y=0' is OK. Post_processes (~'unsharpen') typically need 'PyrMxX=Y=1' or '=2'.

int I16_Init( int AW, int AH, int AModeId ) {
  I16_Release(); fprintf( stderr, cTPS" Init %dx%d,%d ...\n", AW, AH, AModeId );
  #ifdef _WIN32
    SYSTEM_INFO qSI; GetSystemInfo( &qSI ); xBS_SpawnThreads = qSI.dwNumberOfProcessors;
  #else // ifLinux
    xBS_SpawnThreads = std::thread::hardware_concurrency();
  #endif
  xBS_SpawnThreads = mcRng( xBS_SpawnThreads, 2, 8 );
  // it refers to "static tUnSharp cUS??" in VAPrcFlu.inl
  const tUnSharp *cUS = &cUS30; // fill default 'UnSharp-record' by 'ModeId'
  static const tAdjMnMx cAdjMnMxs[] {
    { {0.33f,0.60f,0.60f}, {1.8f,1.50f,1.60f} }, // default adjustment of BCS for Tag=1,15
    { {0.33f,0.88f,0.70f}, {1.8f,1.22f,1.50f} } }; // lesser contrast variations
  const tAdjMnMx  *cAdjMnMx = cAdjMnMxs+0;
  if (AModeId>=200){ cUS = &cUS200;} else
  if (AModeId>=102){ cUS = &cUS102;} else
  if (AModeId>=100){ cUS = &cUS100;} else
  if (AModeId>=90) { cUS = &cUS90; } else
  if (AModeId>=71) { cUS = &cUS71; } else
  if (AModeId>=70) { cUS = &cUS70; } else
  if (AModeId>=40) { cUS = &cUS40; } else
  if (AModeId>=32) { cUS = &cUS32; } else // cAdjMnMx = cAdjMnMxs+1;
  if (AModeId>=30) { cUS = &cUS30; } else
  if (AModeId>=22) { cUS = &cUS22; } else
  if (AModeId>=20) { cUS = &cUS20; } else
  if (AModeId>=10) { cUS = &cUS10; } else // cAdjMnMx = cAdjMnMxs+1; 
  if (AModeId>= 6) { cUS = &cUS06; } else
  if (AModeId>= 5) { cUS = &cUS05; } else
  if (AModeId>= 3) { cUS = &cUS03; } else
  if (AModeId>= 1) { cUS = &cUS01; };    // cAdjMnMx = cAdjMnMxs+1;
  memcpy( &xUS, cUS, sizeof(tUnSharp) );
  char     qS[256]; // ini-file name (in startup directory)
  char *qE = strrchr( xGIni_DLLName, '.' );
  sprintf( qE?qE:strend(qS), "%02d.ini", AModeId );
  GIniFile qIF; qIF.SetFName( xGIni_DLLPath ); fprintf( stderr, cTPS" \"%s\"\n", qIF.rFName );
  qIF.SetSName( "Settings" ); int qI[16];
  if (!xPars[cParCUDA]) if (Gini_sscanf( qIF.SGetStr(
    "nCUDA","0; -1:no, 0:auto, 1:mandatory"), qI )>0) xPars[cParCUDA] = qI[0];
  if (Gini_sscanf( qIF.SGetStr("nPrintLabel","0; -1:neg, 0:no, 1:pos, 2:inv"), qI )>0) xPrintLabel = qI[0];
  if (Gini_sscanf( qIF.SGetStr("nPyrMxXY","1 0; downshiftmx_X _Y"), qI )>0) { xPyrMxX = qI[0]; xPyrMxY = qI[1]; };
  if (Gini_sscanf( qIF.SGetStr("nShrShu", "0; ShiftRight"), qI )>0) xShiftRgt = qI[0];
  if (!xPars[cParIniSetDflt]) if (Gini_sscanf( qIF.SGetStr(
    "nIniSetDflt","0; 0:no, 1:yes"), qI )>0) xPars[cParIniSetDflt] = qI[0];
  qIF.SetSName( "TL" ); for (int qtTLL=0;qtTLL<mcASz(xTLL);qtTLL++) {
    sprintf( qS, "%c", 'a'+qtTLL ); jTLLLoad( &qIF, qS, xTLL+qtTLL, qtTLL==mcASz(xTLL)-1 ); };
  qIF.SetSName( "GC" ); jGrayCutLoad( &qIF, &xGC );
  #if 1 // balancing des_Gamma "UnSharp.fDes[0] in VaPrcFlu.inl" assumes "GC.Gamma=1"
    //if (fabs(xGC.rGamma-1)>1e-3f)  for (int qi=0;qi<mcASz(xUS.fB);qi++)
    //  xUS.fB[qi] = iiGammaByInp( xUS.fB[qi]/xGC.rGamma );
    if (fabs(xGC.rGamma-1)>1e-3f) { float des_Gamma = xUS.fDeG[0];
      if (des_Gamma>1e-3f) xUS.fDeG[0] = iiGammaByInp( des_Gamma/xGC.rGamma );
      if (xDbgDmp>0) { sprintf( qS, "ini des_Gamma by GC.Gamma %.1f/%.1f=%.2f",
        des_Gamma, xGC.rGamma, xUS.fDeG[0] ); WriteLn( qS ); }; };
  #endif
  qIF.SetSName( "CO" ); BS7a_CollimatorLoad( &qIF );
  qIF.SetSName( "RD" ); BS7b_RibsDetectionLoad( &qIF );
  qIF.SetSName( "AE" ); BS5a_AirEstimationLoad( &qIF );
  qIF.SetSName( "QL" ); BS5c_QuadLapsLoad( &qIF );
  qIF.SetSName( "ME" ); BS5d_MskLoad( &qIF );
  qIF.SetSName( "US" ); jUnSharpLoad( &qIF, &xUS, &xUSAdjs, cAdjMnMx ); // load 'US', init 'USAdjs'
  if (xPars[cParCUDA]>=0) { // Pars[ParCUDA] defines usage of CUDA 
    BC_DeviceFind(); if (xPars[cParCUDA]>0 && xBC_Err) return -1;
    if (!xBC_Err) { BC_DeviceReset(); if (!xBC_Err) xPars[cParCUDA] = 1; }; };
  xModeId = AModeId; xN = AW*AH;
  // xN <<= 2; // ????
  BS_SpawnThreads_State( true ); // threads' pool
  if (xPars[cParCUDA]<1) {
    BS_DatasAlloc( mcASz(xDaH), xN, xDaH ); // allocate host buffers
  } else {
    BS_DatasAlloc( mcASz(xDaH)/3, xN, xDaH ); // allocate host buffers
    BC_DatasAllocD( mcASz(xDaD), xN, xDaD ); // allocate device buffers
  };
  return AW|(AH<<16); }

void I16_Release() {
  if (xN<=0) return; fprintf( stderr, cTPS" Release\n" );
  if (xPars[cParCUDA]<1) {
    BS_DatasRelease( mcASz(xDaH), xDaH ); // deallocate host buffers
  } else {
    BS_DatasRelease( mcASz(xDaH)/3, xDaH ); // deallocate host buffers
    BC_DatasReleaseD( mcASz(xDaD), xDaD );
  };
  BS_SpawnThreads_State( false ); // threads' pool
  xN = 0; }

static int inline mcPyr( int APyr, int AL ) { return (AL>>APyr)<<APyr; }
int I16_Invoke( int AW, int AH, const tI16_Inp *ASrc, int AWSrc, tI16_Out *ADst, int AAdj ) {
  if (!ASrc) return 0; // nothing to do
  char       qS[256]; // adhoc string-buffer
  // align dimensions for <possible> downsampling only once to have uniform debug-frames
  const int  cAdjTag = (AAdj>>12)&15,
             cW = mcPyr( xPyrMxX, AW ), cH = mcPyr( xPyrMxY, AH ), cN = cW*cH;
  if (cN>xN) { ErrMsg( "Too large dimension" ); return -64; }; // fatal error
  if (xDbgDmp) fprintf( stderr, cTPS" Invoke (%dx%d %04X) as %dx%d\n", AW, AH, AAdj, cW, cH );
  int        qGx = (AW-cW)>>1, qGy = (AH-cH)>>1;
  const bool cHaveDn0 = (cW!=AW)||(cH!=AH);
  if (cHaveDn0) BS_RowsCpy( AW, ASrc+qGx+qGy*AW, cW, xDaH[0], cW, cH );
           else BS_Cpy( AW*AH, ASrc, xDaH[0] ); 
  if (xShiftRgt) mcShiftRgt<tBS_Card>( cN, (tBS_Card*)xDaH[0], xShiftRgt ); // shift inp pixels right
  if (xDbgDmp>1) BS_Write( "0inp.raw", cN, xDaH[0] );
  tBS_Val  **qDaHStack = xDaH+1, // stack of host-frames
           **qDaDStack = xPars[cParCUDA]>0 ? xDaD : NULL;
  jGrayCut( &xGC, cW, cH, xDaH[0], qDaHStack, xDbgDmp>1?"1":NULL ); // GrayCut apply
  qDaHStack += jCollPreProcess( cW, cH, xDaH[0], qDaHStack, xDbgDmp>1?"2A":NULL ); // collimator
  bool       qAntiScatter = false; // "antiscatter_grid's filter" for gen_medical by bitT#0 
  if ((xModeId<5 || (xModeId>9 && xModeId<90)) && cAdjTag&1) {
    BS6_ThLoLapReset(); BS_Cpy( cN, xDaH[0], *qDaHStack++ ); // (DaH++ keeps zero_bkg)
    for (int qtTLL=0;qtTLL<mcASz(xTLL);qtTLL++) { // statistical suppress of thin & long cartesian laplacians
      const char *cDbgPfx = NULL; if (xDbgDmp>1) { sprintf( qS, "2B%c", 'a'+qtTLL ); cDbgPfx = qS; };
      BS6_ThLoLapSuppress( xTLL+qtTLL, cW, cH, xDaH[0], qDaHStack, qDaDStack, cDbgPfx ); }; // of 'for tTLL'
    BS2_MimicDat( *--qDaHStack, xDaH[0], cN ); // (DaH-- restores zero_bkg as laplacian)
    qAntiScatter = true; }; // if AdjTag&1
  qDaHStack += jRibsPreProcess( cW, cH, xDaH[0], qDaHStack, xDbgDmp>1?"2C":NULL ); // ribs
  #if 0
    const int cWt = cW*8/8, cHt = cH*8/8;
    sprintf( qS, "TEMP %d x%d ", cWt, cHt ); WriteLn( qS );
    const float cShear = 0.07f;
    BS2A_RotateBySkew( 0, cShear, cW, cH, xDaH[0], cWt, cHt, xDaH[1] );
    BS_DatasWrt( "3x", cWt*cHt, xDaH+0, 0, 1 );
    //BS2A_RotateBySkew( -7777, -cShear, cWt, cHt, xDaH[0], cW, cH, xDaH[1] );
    //BS_DatasWrt( "3y", cN, xDaH+0, 0, 1 );
  #endif
  tUnSharpExts qUSExts; jUnSharp( &xUS, &xUSAdjs, AAdj, &qUSExts,
    cW, cH, xDaH[0], qDaHStack, qDaDStack, xDbgDmp>1?"":NULL ); // UnSharp apply
  jRibsPostProcess( cW, cH, xDaH[0], qDaHStack, xDbgDmp>1?"Ya":NULL );
  jCollPostProcess( cW, cH, xDaH[0], qDaHStack, xDbgDmp>1?"Yb":NULL );
  if (xPrintLabel) { // if print Label
    char *qCh = qS+sprintf( qS, "\"%04X\" %s stack=%d", AAdj, xPars[cParCUDA]>0?"CUDA":"CPU", qUSExts.rStackFrmBuf );
    if (qAntiScatter) { qCh += sprintf( qCh, "\nTL" ); // about antiscatter_grid's filter
      for (int qtTLL=0;qtTLL<mcASz(xTLL);qtTLL++) qCh += jTLLLabel( qCh, xTLL+qtTLL, qtTLL ); };
    qCh += sprintf( qCh, "\n" ); qCh += jGrayCutLabel( qCh, &xGC );
    qCh += sprintf( qCh, "\n" ); qCh += jUnSharpLabel( qCh, &xUS );
    qCh += sprintf( qCh, "\n" ); qCh += jUnSharpExtsLabel( qCh, &qUSExts );
    fprintf( stderr, cTPS" label:\n%s\n", qS );
    tBS_Val *qLabelDstImg = xDaH[0]; const tBS_Val *cLabelDstImg = qLabelDstImg;
    const int cLlab = __min(cW,cH); int qLabelGrayVal = 0;
    for (int qi=0;qi<cLlab;qi++) { qLabelGrayVal += cLabelDstImg[0]; cLabelDstImg += cW+1; }; qLabelGrayVal /= cLlab;
    const int cLabelGrayMx = (1<<(14+xShiftRgt))-1, cLabelGrayMn = xShiftRgt ? 0 : -cLabelGrayMx;
    if (abs(xPrintLabel)>1)
      qLabelGrayVal = qLabelGrayVal>cLabelGrayMx>>1 ? cLabelGrayMn : cLabelGrayMx;
    else { qLabelGrayVal = qLabelGrayVal/cLlab+(xPrintLabel<<14);
      qLabelGrayVal = mcRng( qLabelGrayVal, cLabelGrayMn, cLabelGrayMx ); };
    GIni_Lab8x8( cW, cH, qLabelDstImg, qS, qLabelGrayVal );
  }; // if print Label
  if (xDbgDmp>1) { BS_Write( "Zres.raw", cN, xDaH[0] );
    sprintf( qS, "%d %d; W H\n=%d %d; Wo Ho\n%d %d; pyrX pyrY\n%04X; ADJ",
      cW, cH, AW, AH, xPyrMxX, xPyrMxY, AAdj );
    int qL = (int)strlen( qS ); if (qL&1) qS[qL++] = 10; //LF
    BS_Write( "$RecPar.txt", qL>>1, (const tBS_Val*)qS ); };
  if (xShiftRgt) { BS_Rng( cN, xDaH[0], 0, (1<<(16-xShiftRgt))-1 ); // cardinal unharp's output here
    mcShiftLft<tBS_Card>( cN, (tBS_Card*)xDaH[0], xShiftRgt ); }; // shift out pixels right
  if (cHaveDn0) {
    BS_RowsCpy( cW, xDaH[0], AW, ADst+qGx+qGy*AW, cW, cH );
    qGx += (AW-cW*2)&1; qGy += (AH-cH*2)&1;
    BS_Bevels( AW, AH, ADst, qGx, qGx, qGy, qGy );
  } else
    BS_Cpy( cN, xDaH[0], ADst ); 
  return 1; }
