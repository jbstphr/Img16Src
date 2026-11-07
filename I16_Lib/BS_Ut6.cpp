// I16_Ut6.cpp : CPU ad-hoc utilities (with explicit dimension)
// by S.B.
#include "stdafx.h"
#include "BS_Cmn.h"
#include "math.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut3.h"
#include "BS_Ut5.h"
#include "BS_Ut6.h"

//// despeckle (for certain bands) ////

// CPU despeckle (similar code for ~general_denoise, ~horz_defects_correction).
// straightforward despeckle  org: ...............  res: ...............
// by median filter for            VVVVVVVVVVVVVVV       MMMMMMMMMMMMMMM
// enlisted bands of rows          ...............       ...............
static int ixDespeckleScalars[2]; // scalar parameters of cplb_Despeckle; #0:threshold; '#1=1':strong,'=0':general

void BS6_DespeckleScalars( int AThr, int AStrong ) { // -- despeckle scalar params
  ixDespeckleScalars[0] = AThr; ixDespeckleScalars[1] = AStrong; }

// general approach#0: for ~Unsharp (i.e. ~general); strong approach#1:for ~MBApply (i.e. for ~horz_defects_correction).
// general approach accounts median and isotropic_average with ~similar importance
// strong approach accounts less median but more anisotropic average (by vert_neighbors,ignores horz_neighbors)
static int __fastcall jAve3x3( const tBS_Val *AInp ) { return ((mcSum<tBS_Val,3>(AInp+0))+(mcSum<tBS_Val,3>(AInp+6)))/6; }
static int __fastcall jAve5x5( const tBS_Val *AInp ) { return ((mcSum<tBS_Val,10>(AInp+0))+(mcSum<tBS_Val,10>(AInp+15)))/20; }
static void jAve3x3( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { // average 3x3
  if (ixDespeckleScalars[1]) BS_CPK( 1, jAve3x3, AW, AH, ASrc, ADst ); else BS_Ave3x3( AW, AH, ASrc, ADst ); }
static void jAve5x5( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { // average 5x5
  if (ixDespeckleScalars[1]) BS_CPK( 2, jAve5x5, AW, AH, ASrc, ADst ); else BS_Ave5x5( AW, AH, ASrc, ADst ); }
static inline int mcDespeckleStrong( int AV, int AV1, int AV2, int AThr ) { // strong anisotropic approach
  int qD = ((AV1+AV2)>>1)-AV; if (abs(qD)>AThr) AV += qD; return AV; }
static int __fastcall cpb_DespeckleGen( tBS_Val *ADat ) { // general isotropic approach
  int qV = ADat[0], qV1 = ADat[1], qV2 = ADat[2]; qV1 -= qV; qV2 -= qV;
  int qD = abs(qV1)<abs(qV2)?qV1:qV2, cThr = ixDespeckleScalars[0], qA = abs(qD)-cThr; if (qA<=0) return qV;
  if (qA>cThr) qA = cThr; qA >>= 2; if (qD>0) qA = -qA; return qV+qD+qA; }
void BS6_Despeckle( int AW, int AH, int AH0, int AH1, // -- despeckle within a band
  tBS_Val *ADat, int ARadMed, tBS_Val *AMed, int ARadAve, tBS_Val *AAve ) { //
  // note: mcCPLSpawn in BS_Ut0.cpp doesn't '_SpawnPool' unless Dim1>99 (i.e. until #H>99)
  switch (ARadAve) { // isotropic or anisotropic average (Med is temp_buffer now)
    case 2: jAve3x3( AW, AH, ADat, AMed ); jAve3x3( AW, AH, AMed, AAve ); break; case  3: jAve5x5( AW, AH, ADat, AAve ); break;
    case 4: jAve5x5( AW, AH, ADat, AMed ); jAve5x5( AW, AH, AMed, AAve ); break; default: jAve3x3( AW, AH, ADat, AAve ); };
  BS_Med( mcRng(ARadMed,1,4), AW, AH, ADat, AMed ); 
  const int cN = AW*(AH-AH0-AH1), cOfs = AH0*AW; ADat += cOfs; AMed += cOfs; AAve += cOfs;
  if (ixDespeckleScalars[1]) { const int cThr = ixDespeckleScalars[0]; // strong anisotropic approach
    for (int qi=0;qi<cN;qi++) ADat[qi] = mcDespeckleStrong( ADat[qi], AMed[qi], AAve[qi], cThr );
  } else { // general isotropic approach
    tBS_Val *qDs[] = {ADat,AMed,AAve}; BS_CPB( cpb_DespeckleGen, mcASz(qDs)|cCPBNDD, cN, qDs );
  }; }

// CUDA general despeckle (as 'Strong=false' for CPU ).
void BS6_DespeckleCUDA( int AW, int AH, tBS_Val *ADat, // -- general despeckle
  int ARadMed, tBS_Val *AMed, int ARadAve, tBS_Val *AAve ) {
  BC_Dim( AW, AH ); switch (ARadAve) { // plain isotropic average (Med is temp_buffer now)
    case  2: BC_Ave3x3( ADat, AMed ); BC_Ave3x3( AMed, AAve ); break; case  3: BC_Ave5x5( ADat, AAve ); break;
    case  4: BC_Ave5x5( ADat, AMed ); BC_Ave5x5( AMed, AAve ); break; default: BC_Ave3x3( ADat, AAve ); };
  BC_Med( ADat, AMed, mcRng(ARadMed,1,4) );
  const int cThr = ixDespeckleScalars[0]; // standard isotropic approach
  extern void KDespeckle(int,tBS_Val*,const tBS_Val*,const tBS_Val*); KDespeckle( cThr, ADat, AMed, AAve ); }

//// statistical suppress of thin & long cartesian laplacians //// 
// Mode=1 "plain suppression" keeps ~fixed code (as the most robust statistically; though
// non-optimal in particular). Make another modes & procedure for experiments.
static const int ctll_NrmLg = 9,
                 ctll_Nrm = 1<<ctll_NrmLg; 
static tBS_Val  *xTLLMsk = NULL;

static void iThLoLapSet_CPU( int AW, int AH, // CPU suppression of TLL (set optional mask)
  tBS_Val *ADat, tBS_Val **ATs, int ACGN0, int ACGN1, float AAtt, const char *ADbgPfx ) {
  const int cN = AW*AH; BS2_Dim( AW, AH );
  BS_Cpy( cN, ADat, ATs[1] ); BS_DrA3x3( AW, AH, ATs[1] );
  const int cCorr = (int)(AAtt*ctll_Nrm); BS_Cpy( cN, ATs[1], ATs[0] ); 
  tBS2_Cut cCut = { {ACGN0,ACGN1,cCorr}, {0,0} }; BS2_Cut( &cCut, ATs[0] ); BS_VAdd( cN, ATs[0], ctll_Nrm-cCorr );
  if (ADbgPfx) { tBS_Val *qDs[] = {ATs[1],ATs[0]}; BS_DatasWrt( ADbgPfx, cN, qDs, 0, mcASz(qDs) ); };
  xTLLMsk = ATs[0]; }

static void iThLoLapSuppress_CPU( int AW, int AH, // CPU suppression of TLL (execute)
  tBS_Val *ADat, tBS_Val **ATs, int ARX, int ARY, float AAtt, const char *ADbgPfx ) {
  const int cN = AW*AH; tBS_Val *qT = ATs[0];
  BS_Cpy( cN, ADat, qT );
  if (ARX<ARY) { BS_MovAveX( ARX, AW, AH, qT ); BS_SubA( cN, ADat, qT ); BS_MovAveY( ARY, AW, AH, qT ); }  // vrt
          else { BS_MovAveY( ARY, AW, AH, qT ); BS_SubA( cN, ADat, qT ); BS_MovAveX( ARX, AW, AH, qT ); }; // hrz
  BS_VMul( cN, qT, AAtt ); if (xTLLMsk) BS_Mul( cN, xTLLMsk, qT, -ctll_NrmLg ); BS_Add( cN, qT, ADat );
  if (ADbgPfx) { tBS_Val *qDs[] = {ADat,qT}; BS_DatasWrt( ADbgPfx, AW*AH, qDs, 1, mcASz(qDs) ); }; }

static void iThLoLapSet_CUDA( int AW, int AH, // CUDA suppression of TLL (set optional mask)
  tBS_Val *ADat, tBS_Val **ATs, int ACGN0, int ACGN1, float AAtt, const char *ADbgPfx ) {
  const int cN = AW*AH; BC_Dim( AW, AH );
  BC_Move( cN*cBS_Sz, ADat, ATs[1] ); BC_DrA3x3( ATs[1], ATs+2 );
  const int cCorr = (int)(AAtt*ctll_Nrm); BC_Cpy( ATs[1], ATs[0] ); 
  tBS2_Cut cCut = { {ACGN0,ACGN1,cCorr}, {0,0} }; BS3_Cut( &cCut, ATs[0] ); BC_VAdd( ATs[0], ctll_Nrm-cCorr );
  if (ADbgPfx) { tBS_Val *qDs[] = {ATs[1],ATs[0]}; BC_DatasWrt( ADbgPfx, qDs, 0, mcASz(qDs) ); };
  xTLLMsk = ATs[0]; }

static void iThLoLapSuppress_CUDA( int AW, int AH, // plain CUDA suppression of TLL (limited radius)
  tBS_Val *ADat, tBS_Val **ATsD, int ARX, int ARY, float AAtt, const char *ADbgPfx ) {
  const int cN = AW*AH; BC_Dim( AW, AH ); tBS_Val *qDat = ATsD[0], *qT0 = ATsD[1], *qT1 = ATsD[2];
  BC_Move( cN*cBS_Sz, ADat, qDat );
  if (ARX<ARY) { BC_Av1( qDat, qT0, ARX ); BC_SubA( qDat, qT0 ); BC_Av5( qT0, qT1, ARY ); }  // vrt
          else { BC_Av5( qDat, qT0, ARY ); BC_SubA( qDat, qT0 ); BC_Av1( qT0, qT1, ARX ); }; // hrz
  BC_VMul( qT1, AAtt ); if (xTLLMsk) BC_Mul( xTLLMsk, qT1, -ctll_NrmLg ); BC_Add( qT1, qDat );
  if (ADbgPfx) { tBS_Val *qDs[] = {qDat,qT1}; BC_DatasWrt( ADbgPfx, qDs, 1, mcASz(qDs) ); };
  BC_Move( -cN*cBS_Sz, ADat, qDat ); }

void BS6_ThLoLapReset() { xTLLMsk = NULL; }

void BS6_ThLoLapSuppress( const tThLoLap *ATLL, // suppress thin & long laplacians
  int AW, int AH, tBS_Val *ADat, tBS_Val **ATsH, tBS_Val **ATsD, const char *ADbgPfx ) {
  switch (ATLL->n[cBS6_tllMode]) {
    case 0: return;
    case 1: { int qCGN0 = ATLL->n[cBS6_tllCGN0], qCGN1 = ATLL->n[cBS6_tllCGN1]; float qA = ATLL->f[cBS6_tlfAtt];
        if (ATsD) iThLoLapSet_CUDA( AW, AH, ADat, ATsD, qCGN0, qCGN1, qA, ADbgPfx );
             else iThLoLapSet_CPU ( AW, AH, ADat, ATsH, qCGN0, qCGN1, qA, ADbgPfx );
    }; break;
    case 2: { int qRX = ATLL->n[cBS6_tllRadX], qRY = ATLL->n[cBS6_tllRadY]; float qA = ATLL->f[cBS6_tlfAtt];
      if (!(qRX|qRY) || fabs(qA)<1e-3f) break; const int cRadMxCUDA = 15;
      if (ATsD && (qRX>cRadMxCUDA||qRY>cRadMxCUDA)) {
        fprintf( stderr, cTPS" WRN ThLoLapSuppress_CUDA(%d,%d) radius>%d\n", qRX, qRY, cRadMxCUDA );
        qRX = __min(qRX,cRadMxCUDA); qRY = __min(qRY,cRadMxCUDA); // ATsD = NULL;
      };
      if (ATsD) iThLoLapSuppress_CUDA( AW, AH, ADat, ATsD+1, qRX, qRY, qA, ADbgPfx );
           else iThLoLapSuppress_CPU ( AW, AH, ADat, ATsH+1, qRX, qRY, qA, ADbgPfx );
    }; break;
    default: ErrMsg( "BS6_ThLoLapSuppress params" );
  }; } // switch

//// Radon's despeckle ////
static float ix_RDLimK; static void iRadonDespeckleParams( int *AR, int *AFlg, float ALimK ) {
  int qR = *AR, qFlg = *AFlg; if (qR<1) qR = 1;
  if (qFlg) { if (qR<2) qFlg &= ~4; }
       else { if (qR<2) qFlg = 3; else if (qR<3) qFlg = 6; else qFlg = 7; };
  *AR = qR; *AFlg = qFlg; ix_RDLimK = ALimK*1.1f; } // !!! balance most important threshold against CUDA

template<int cNDD> static int __fastcall mccpb_RadonDespeckle( tBS_Val *ADat ) { // -- Radon's despeckle CPU
  int qMn = ADat[1], qMx = qMn, qAv = qMn; for (int qi=2;qi<=cNDD;qi++) { 
    int qv = ADat[qi]; if (qMn>qv) qMn = qv; else if (qMx<qv) qMx = qv; qAv += qv; };
  qAv /= cNDD; const float cLimK = ix_RDLimK; int qV0 = ADat[0];
  qMn = qAv+(int)(cLimK*(qMn-qAv)); if (qV0<=qMn) { qV0 = qMn; goto labLeave; };
  qMx = qAv+(int)(cLimK*(qMx-qAv)); if (qV0>qMx) qV0 = qMx; labLeave:
  ADat[0] = qV0-qAv; ADat[1] = qAv; return 2; } 

static void iRadonDespeckle_CPU( // -- Radon's despeckle CPU
  int AW, int AH, tBS_Val *ADat, tBS_Val *AAve, tBS_Val **ATs, int AR, int AFlg, float ALimK ) {
  if (!AFlg) return; // inp Dat:org_frame; out Dat:laplacian, Ave:average
  iRadonDespeckleParams( &AR, &AFlg, ALimK ); int qNDD = 0, qDirs[8];
  if (AFlg&1) { qDirs[qNDD++] = 1; qDirs[qNDD++] = 5; };
  if (AFlg&2) { qDirs[qNDD++] = 3; qDirs[qNDD++] = 7; };
  if (AFlg&4) { qDirs[qNDD++] = 2; qDirs[qNDD++] = 4; qDirs[qNDD++] = 6; qDirs[qNDD++] = 8; };
  tBS_Val *qTs[9]; qTs[0] = ADat; qTs[1] = AAve; memcpy( qTs+2, ATs, (qNDD-1)*sizeof(void*) );
  const int cN = AW*AH; BS_DataDup( qNDD, cN, qTs[0], qTs+1 );
  BS_MovAve( AR, AW, AH, qTs+1, qDirs, qNDD ); // BS_DatasWrt( "qq", cN, ATs, 0, qNDD );
  tBS_CPBProc qProc; switch (qNDD) { //
    case  2: qProc = mccpb_RadonDespeckle<2>; break;
    case  4: qProc = mccpb_RadonDespeckle<4>; break;
    case  6: qProc = mccpb_RadonDespeckle<6>; break;
    default: qProc = mccpb_RadonDespeckle<8>; };
  BS_CPB( qProc, qNDD+1, cN, qTs ); }

static void iRadonDespeckle_CUDA(
  tBS_Val *ADat, tBS_Val *AAve, tBS_Val **ATs, int AR, int AFlg, float ALimK ) { // -- Radon's despeckle CUDA
  // inp Dat:org_frame; out Dat:laplacian, Ave:average
  iRadonDespeckleParams( &AR, &AFlg, ALimK ); tBS_Val *qA[8];
  memcpy( qA, ATs, sizeof(qA) ); int qNDD = 0;
  if (AFlg&1) { BC_Av1( ADat, qA[qNDD++], AR ); BC_Av5( ADat, qA[qNDD++], AR ); };
  if (AFlg&2) { BC_Av3( ADat, qA[qNDD++], AR ); BC_Av7( ADat, qA[qNDD++], AR ); };
  if (AFlg&4) { BC_Av2( ADat, qA[qNDD++], AR ); BC_Av4( ADat, qA[qNDD++], AR );
                BC_Av6( ADat, qA[qNDD++], AR ); BC_Av8( ADat, qA[qNDD++], AR ); };
  // BC_DatasWrt( "qq", qA, 0, qNDD );
  extern void KRadonDespeckle(float,tBS_Val*,tBS_Val*,tBS_Val**); 
  memset( qA+qNDD, 0, (8-qNDD)*sizeof(void*) ); KRadonDespeckle( ALimK, ADat, AAve, qA ); }

//// Emphasize scrape of high-freq grid (test-pattern with lines) ////
// Term "scrape" here means symm high-freq laps ~"--+--", "-0+0-". Don't emphasize non-symm laps
// ~"+++--", "--+++". 
// CPU-version (rather slow,acccurate): deals with "directional_hypothetical_grids". CUDA-version
// (rather fast,rough) deals with "summary_hypothetical_grid".
static inline int mcScrapeGridAdd( int AR, int A0, int A1 ) { A0 += A1; A0 -= A0>>(3-AR); return A0; }
static inline int iScrapeGrid( const tBS_Val *ADat ) { // ~symmetric_high_freqs in a pixel
  int qv = ADat[0]; if (!qv) return 0; int qu = 0;
  for (int qr=1;qr<3;qr++) { int q0 = ADat[-qr], q1 = ADat[qr];
    if (qv<0) { if (q0<=0 || q1<=0) continue; q0 = mcScrapeGridAdd(qr,q0,q1); qu += __max(qv,-q0); }
         else { if (q0>=0 || q1>=0) continue; q0 = mcScrapeGridAdd(qr,q0,q1); qu += __min(qv,-q0); };
  }; return qu; };
static int icpl_ScrapeGrid( int AL, tBS_Val *ADat ) { // trace_shader reveals ~symmetric_high_freqs
  const int cR = 2; if (AL<cR*5) { BS_ZFil( AL, ADat ); return AL; }
  const int cL = BS_Align( AL ); tBS_Val *qT1 = ADat+cL, *qT2 = qT1+cL;
  mcMovAve<tBS_Val>( AL, 1, ADat, qT1 ); mcMovAve<tBS_Val>( AL, 2, ADat, qT2 );
  BS_SubS( AL, ADat, qT1 ); BS_SubS( AL, ADat, qT2 );
  // dir_hypo_grid estimates as ~"{scrape(R=1)*3+scrape(R=2)*2}/5"
  BS_ZFil( cR, ADat ); BS_ZFil( cR, ADat+AL-cR ); for (int qi=cR;qi<AL-cR;qi++) {
    ADat[qi] = (iScrapeGrid(qT1+qi)*5+(iScrapeGrid(qT2+qi)*3))>>3; }; return AL; }
static int __fastcall icpb_ScrapeGrid( tBS_Val *ADat ) { int M = ADat[0], Z = M; // result ~ave_limb stat
  for (int qi=1;qi<4;qi++) { int u = ADat[qi]; if (abs(M)>abs(u)) M = u; Z += u; }; return (Z-M)>>1; }
static void iScrapeGrid_CPU( int AW, int AH, const char *ADbgPfx, char *ADbgSfx, // emphasize high_freqs's grid in CPU
  tBS_Val *AFrm, tBS_Val **ATs, float AKgr ) { if (AKgr<1e-3f) return;            // with detailed diagnostic
  const int cN = AW*AH, cNDD = 4, cDirs[cNDD] = {1,3,5,7}; BS_DataDup( cNDD, cN, AFrm, ATs );
  BS_CPLA( icpl_ScrapeGrid, cNDD, cDirs, AW, AH, ATs ); // estimation dir_hypo_grids
  BS_CPB( icpb_ScrapeGrid, cCPBNDD|4, cN, ATs ); // BS_CPB( mccpb_Ave<cNDD>, cCPBNDD|cNDD, cN, ATs );  
  BS_VMul( cN, ATs[0], AKgr ); BS_Add( cN, ATs[0], AFrm ); // Zfrm+=des_KG*SrapeGrid 
  if (1 && ADbgPfx) { strcpy( ADbgSfx, "6Axgrd.raw" ); BS_Write( ADbgPfx, cN, ATs[0] ); };
  if (0 && ADbgPfx) { strcpy( ADbgSfx, "6Ayinp.raw" ); BS_Write( ADbgPfx, cN, AFrm ); }; }

static void iScrapeGrid_CUDA( const char *ADbgPfx, char *ADbgSfx,                // emphasize high_freqs's grid in CUDA
  tBS_Val *AFrm, tBS_Val **ATs, float AKgr ) { if (AKgr<1e-3f) return;             // with detailed diagnostic
  BC_Ave4( AFrm, ATs, 1 ); BC_Ave4( AFrm, ATs+4, 2 ); for (int qi=0;qi<8;qi++) BC_SubS( AFrm, ATs[qi] );
  extern void KScrapeGrid(tBS_Val**); // inp #0..3:lap(R=1), #4..7:lap(R=2); out #0:sum_hypo_grid
  KScrapeGrid( ATs ); // estimate sum_hypo_grid
  BC_VMul( ATs[0], AKgr ); BC_Add( ATs[0], AFrm ); // Zfrm+=des_KG*SrapeGrid 
  if (1 && ADbgPfx) { strcpy( ADbgSfx, "6Axgrd.raw" ); BC_Write( ADbgPfx, 0, ATs[0] ); };
  if (0 && ADbgPfx) { strcpy( ADbgSfx, "6Ayinp.raw" ); BC_Write( ADbgPfx, 0, AFrm ); }; }

//// Dispatch high-freqs's scrapes ////
// "empirical reshuffle of high_freqs" assumes that
// a) Ave has ~clean mid_freqs; 
// b) ~residual Aves' high_freqs are more regular than Lap's ~original ones;
// c) Lap shall be blurred before ~'<Radon>_Noise(Lap)'.
static void iScrapeDbgTrc_CPU( const char *ADbgPfx, const tBS_Val *ADat, tBS_Val **ATs, int ARadHF ) {
  const int cW = BS2_W(), cH = BS2_H(), cN = cW*cH;
#if 0 // !!!!
  BS_Write( ADbgPfx, cN, ADat ); } 
#else
  BS2_MovAve2D( ARadHF, ADat, ATs[0], cW, cH ); BS_SubS( cN, ADat, ATs[0] ); BS_Write( ADbgPfx, cN, ATs[0] ); }
#endif

static void iScrapenPFetch_CPU( int AW, int AH, const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF, // fetch non-proc-scrape in CPU
  tBS_Val *AFrm, tBS_Val *ASnP, tBS_Val **ATs, float AKnP, int ALimnP ) {                            // with detailed diagnostic
  if (1 && ADbgPfx) { strcpy( ADbgSfx, "6Ba_Frm.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, AFrm, ATs, Adiag_RadHF ); }; // trace by diag_RadHF
  if (!ASnP && fabs(AKnP)<1e-3f) return; const int cN = AW*AH; if (!ASnP) ASnP = *ATs++; BS_Ave3x3( AW, AH, AFrm, ASnP );
  BS_SubS( cN, AFrm, ASnP ); BS_RngA( cN, ASnP, ALimnP ); // fixed_scrape_radius=1
  if (1 && ADbgPfx && Adiag_RadHF!=1) { strcpy( ADbgSfx, "6Bb_SnP.raw" ); BS_Write( ADbgPfx, cN, ASnP ); };
  if (fabs(AKnP)<1e-3f) return; BS_Cpy( cN, ASnP, ATs[0] ); BS_VMul( cN, ATs[0], AKnP ); BS_Add( cN, ATs[0], AFrm ); 
  if (0 && ADbgPfx) { strcpy( ADbgSfx, "6Bc_Frm.raw" ); BS_Write( ADbgPfx, cN, AFrm ); }; } // Frm+=des_KnP*ScrapenP

static void iScrapeStdDispatch_CPU( int AW, int AH, const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF, // dispatch standard scrape in CPU
  tBS_Val *AAve, tBS_Val *ALap, tBS_Val **ATs, const float AKdSS[4] ) { const int cN = AW*AH;             // with detailed diagnostic
  if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Ca_Ave.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, AAve, ATs, Adiag_RadHF );
    strcpy( ADbgSfx, "6Ca_Lap.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, ALap, ATs, Adiag_RadHF ); };
  if ((fabs(AKdSS[0])>1e-3f) || (fabs(AKdSS[2])>1e-3f)) { BS_Ave3x3( AW, AH, AAve, ATs[0] ); BS_SubS( cN, AAve, ATs[0] ); };
  if ((fabs(AKdSS[1])>1e-3f) || (fabs(AKdSS[3])>1e-3f)) { BS_Ave3x3( AW, AH, ALap, ATs[1] ); BS_SubS( cN, ALap, ATs[1] ); };
  if (BS2_Sum( 1.0f, AAve, AKdSS[0], ATs[0], AKdSS[1], ATs[1] )) { if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Cc_Ave.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, AAve, ATs+2, Adiag_RadHF ); }; };
  if (BS2_Sum( 1.0f, ALap, AKdSS[2], ATs[0], AKdSS[3], ATs[1] )) { if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Cc_Lap.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, ALap, ATs+2, Adiag_RadHF ); }; };
} // of ScrapeDispatch_CPU

static void iScrapenPAdd_CPU( int AW, int AH, const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF,               // add non_processed scrape
  tBS_Val *AFrm, tBS_Val *ASnP, tBS_Val *ALoF, tBS_Val **ATs, int ARadBlr, int Arcu0, int Arcu1, float AKSnP ) { // with detailed diagnostic
  if (AKSnP<1e-3f) return; const int cN = AW*AH; const bool cAbs = ARadBlr<0; ARadBlr = abs(ARadBlr);
  if (ARadBlr && (Arcu0!=Arcu1)) { tBS_Val *qMsk = ATs[0]; BS_Cpy( cN, ALoF, qMsk );
    BS2_MovAve2Dp( ARadBlr%100, qMsk, ATs+1, ARadBlr/100, AW, AH ); BS_SubS( cN, ALoF, qMsk ); if (cAbs) BS_Abs( cN, qMsk );
    if (1 && ADbgPfx) { strcpy( ADbgSfx, "8Blol.raw" ); BS_Write( ADbgPfx, cN, qMsk ); }; // low_freq laplacian
    tBS2_Cut cCut = {{Arcu0,Arcu1,512},{0,0}}; BS2_Cut( &cCut, qMsk );
    if (1 && ADbgPfx) {strcpy( ADbgSfx, "8Cmsk.raw" ); BS_Write( ADbgPfx, cN, qMsk ); };
    BS_Mul( cN, qMsk, ASnP, -9 );
    if (1 && ADbgPfx) {strcpy( ADbgSfx, "8Dsnp.raw" ); BS_Write( ADbgPfx, cN, ASnP ); }; }; // if RadBlr
  BS_VMul( cN, ASnP, AKSnP ); BS_Add( cN, ASnP, AFrm ); }

static void iScrapeDbgTrc_CUDA( const char *ADbgPfx, const tBS_Val *ADat, tBS_Val **ATs, int ARadHF ) {
  BS3_MovAve2D( ARadHF, ADat, ATs[0], ATs[1] ); BC_SubS( ADat, ATs[0] ); BC_Write( ADbgPfx, 0, ATs[0] ); }

static void iScrapenPFetch_CUDA( const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF, // fetch non-proc-scrape in CUDA
  tBS_Val *AFrm, tBS_Val *ASnP, tBS_Val **ATs, float AKnP, int ALimnP ) {             // with detailed diagnostic
  if (1 && ADbgPfx) { strcpy( ADbgSfx, "6Ba_Frm.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, AFrm, ATs, Adiag_RadHF ); }; // trace by diag_RadHF
  if (!ASnP && fabs(AKnP)<1e-3f) return; if (!ASnP) ASnP = *ATs++; BC_Ave3x3( AFrm, ASnP );
  BC_SubS( AFrm, ASnP ); BC_RngA( ASnP, ALimnP ); // fixed_scrape_radius=1
  if (1 && ADbgPfx && Adiag_RadHF!=1) { strcpy( ADbgSfx, "6Bb_SnP.raw" ); BC_Write( ADbgPfx, 0, ASnP ); };
  if (fabs(AKnP)<1e-3f) return; BC_Cpy( ASnP, ATs[0] ); BC_VMul( ATs[0], AKnP ); BC_Add( ATs[0], AFrm );
  if (0 && ADbgPfx) { strcpy( ADbgSfx, "6Bc_Frm.raw" ); BC_Write( ADbgPfx, 0, AFrm ); }; } // Frm+=des_KnP*ScrapenP

static void iScrapeStdDispatch_CUDA( const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF, // dispatch standard scrape in CUDA
  tBS_Val *AAve, tBS_Val *ALap, tBS_Val **ATs, const float AKdSS[4] ) {                    // with detailed diagnostic
  if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Ca_Ave.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, AAve, ATs, Adiag_RadHF );
    strcpy( ADbgSfx, "6Ca_Lap.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, ALap, ATs, Adiag_RadHF ); };
  if ((fabs(AKdSS[0])>1e-3f) || (fabs(AKdSS[2])>1e-3f)) { BC_Ave3x3( AAve, ATs[0] ); BC_SubS( AAve, ATs[0] ); };
  if ((fabs(AKdSS[1])>1e-3f) || (fabs(AKdSS[3])>1e-3f)) { BC_Ave3x3( ALap, ATs[1] ); BC_SubS( ALap, ATs[1] ); };
  if (BS3_Sum( 1.0f, AAve, AKdSS[0], ATs[0], AKdSS[1], ATs[1] )) { if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Cc_Ave.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, AAve, ATs+2, Adiag_RadHF ); }; };
  if (BS3_Sum( 1.0f, ALap, AKdSS[2], ATs[0], AKdSS[3], ATs[1] )) { if (1 && ADbgPfx) { // trace by diag_RadHF
    strcpy( ADbgSfx, "6Cc_Lap.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, ALap, ATs+2, Adiag_RadHF ); }; };
} // of ScrapeDispatch_CUDA

static void iScrapenPAdd_CUDA( const char *ADbgPfx, char *ADbgSfx, int Adiag_RadHF,                              // add non_processed scrape
  tBS_Val *AFrm, tBS_Val *ASnP, tBS_Val *ALoF, tBS_Val **ATs, int ARadBlr, int Ar0SnP, int Ar1SnP, float AKSnP ) { // with detailed diagnostic
  if (AKSnP<1e-3f) return; const bool cAbs = ARadBlr<0; ARadBlr = abs(ARadBlr);
  if (ARadBlr && (Ar0SnP!=Ar1SnP)) { tBS_Val *qMsk = ATs[0]; BC_Cpy( ALoF, qMsk );
    BS3_MovAve2Dp( ARadBlr%100, qMsk, ATs+1, ARadBlr/100 ); BC_SubS( ALoF, qMsk ); if (cAbs) BC_Abs( qMsk );
    if (1 && ADbgPfx) {strcpy( ADbgSfx, "8Blol.raw" ); BC_Write( ADbgPfx, 0, qMsk ); }; // low_freq laplacian
    tBS2_Cut cCut = {{Ar0SnP,Ar1SnP,512},{0,0}}; BS3_Cut( &cCut, qMsk );
    if (1 && ADbgPfx) { strcpy( ADbgSfx, "8Cmsk.raw" ); BC_Write( ADbgPfx, 0, qMsk ); };
    BC_Mul( qMsk, ASnP, -9 );
    if (1 && ADbgPfx) { strcpy( ADbgSfx, "8Dsnp.raw" ); BC_Write( ADbgPfx, 0, ASnP ); }; }; // if RadBlr
  BC_VMul( ASnP, AKSnP ); BC_Add( ASnP, AFrm ); }

//// Combined Poisson's denoise ////
//  if Thr then                  it traces "1*.raw"  
//    Despeckle                  frm:="the_org"
//  if RadR then                 it traces "6*.raw"
//    if KSnP then
//      ScrapeGrid               frm+=scrape_grid
//      Scrape_nProc_Fetch       qfrm-=scrape_nproc (separate non-processing scrape)
//    Radon_Radial_Despeckle     lap ave; lap contains mid_freqs with ~aliased high_freqs
//    if exists(Scrape) then     lap ave; reshuffle high_freqs between lap & ave
//      ScrapeStd_Dispatch                ave contains low_freqs with ~smooth high_freqs
//  if RadN then                 it traces "3*.raw"; "Radon|Shannon_noise for mid_freq laplacians" 
//                               noise:=Noise(RadN,lap); 
//    if RadN2 then              noise-=KN2*MaxLap(RadN2,noise); "refine" the noise
//                               lap-=noise
//  if RadR then                 frm:=Kmf*lap+1.0*ave
//    Scrape_nProc_Add           frm+=scrape_nproc
extern float iiGammaByInp(float);
void BS6_DePoisson( int AW, int AH, const tUnSharp *AUS, tUnSharpExts *AUSExts,
  tBS_Val *AFrm, tBS_Val **ATsH, tBS_Val **ATsD, const char *ADbgPfx ) { // Poisson_noise removal
  CBS_Ticks     qTM; if (xDbgDmp==1) strcpy( qTM.rPfx, "jUS" ); qTM.Reset(); 
  const int     cN = AW*AH; char *lDbgSfx = NULL, qDbgPfx[88]; // local debug prefix
  if (ADbgPfx) { strcpy( qDbgPfx, ADbgPfx ); ADbgPfx = qDbgPfx; lDbgSfx = strend(qDbgPfx); };
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "4AInp.raw" ); BS_Write( ADbgPfx, cN, AFrm ); };
  BS2_Dim( AW, AH );
  // - despeckle_stage's params
  int     des_Thr   = AUS->rDeG[0]; // initial despeckle by threshold
  float   des_Gamma = iiGammaByInp(AUS->fDeG[0]),
          des_Kgr   = AUS->fDeG[1], // ScrapeGrid...(...,Kgr); ~0
          des_KnP   = AUS->fDeG[2], // ScrapeFetch...(...,KnP); ~0.6
          des_KSnP  = AUS->fDeG[3]; // ~0.2
  int     des_RadR  = AUS->rDeR[0]; // Radon_Radial_Despeckle
  float   des_KR    = AUS->fDeR[0], // RadonDespeckle(...,KR); ~0.6
          des_KU    = AUS->fDeR[1], // ~0.3; "_MxLapXLimbo(des_KU,...,-1)"
          des_Kmf   = AUS->fDeR[2]; // Denoise_Res:=Ave+Kmf*Mid_Freq; ~1.0
  float   des_KdSS[4]; memcpy( des_KdSS, AUS->fDeD, sizeof(des_KdSS) ); // A+a* A+l* L+a* L+l*
  int     des_RadN  = AUS->rDeN[0], // Radon_Shannon_denoise
          des_LgN   = AUS->rDeN[1], // log2 of shannon_noise's threshold
          des_RadN2 = AUS->rDeN[2], // morpho_residuals' radius
          des_ZN2   = AUS->rDeN[3]; // morpho_residuals' zero-range?
  float   des_KN    = AUS->fDeN[0], // -=KN*shannon_noise
          des_KN2   = AUS->fDeN[1]; // shannon_noise-=KN2*morpho_residuals
  int     des_RadSnP= AUS->rDeS[0], // ~203
          des_r0SnP = AUS->rDeS[1], // ~111
          des_r1SnP = AUS->rDeS[2]; // ~-3456
  const int cFlgR[] = {3,6,7}, // RadR=1... (few sentinels are below...)
            des_FlgR = cFlgR[mcRng(des_RadR,1,mcASz(cFlgR))-1], // Radon's despeckle
            diag_RadHF = __max(AUS->rA[0],1); // diagnostic output of high-freqs
  int     des_LimnP = (int)(abs(des_Thr)*2.4f); // ScrapeFetch...(...,LimT); ~Thr
  if (des_KnP>1.0f) des_LimnP = (int)(des_LimnP/sqrt(des_KnP)); 
  if (!des_RadR || des_KN<1e-3f) des_RadN = 0; if (des_KN2<1e-3f) des_RadN2 = 0;
  if (!ATsD) { // == CPU
    tBS_Val *pFrm = ATsH[0], **qTsH = ATsH+1; BS_Cpy( cN, AFrm, pFrm );
    tBS_Val *pOrg = NULL; if (des_RadR) { pOrg = *qTsH++; BS_Cpy( cN, pFrm, pOrg ); }; // to care about zero_bacground 
    tBS_Val *pSnP = NULL; // non_processed part of fetched scrape
    tBS_Val *pLoF = NULL; // low freqs (moving average)
    // strcpy( lDbgSfx, "0aq.raw" ); BS_Write( ADbgPfx, cN, qTsH[-1] ); // rarely trace
    // strcpy( lDbgSfx, "0aq.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, pFrm, qTsH, 1 ); // rarely trace
    if (fabs(des_Gamma-1)>1e-3f) { BS2_Gamma( des_Gamma, 1<<14, pFrm ); if (1 && ADbgPfx) { // balance_gamma
      strcpy( lDbgSfx, "4BInpGamma.raw" ); BS_Write( ADbgPfx, cN, pFrm ); }; qTM.Mark( "des_Gamma" ); };
    if (des_Thr>0) { // initial despeckle by threshold
      if (1 && ADbgPfx) { // trace by diag_RadHF
        strcpy( lDbgSfx, "5AInp.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, pFrm, qTsH+0, diag_RadHF ); };
      BS6_DespeckleScalars( des_Thr, 0 ); // strong=0: means general_despeckle
      if (1 && ADbgPfx) BS_Cpy( cN, pFrm, qTsH[2] );
      BS6_Despeckle( AW, AH, 0, 0, pFrm, 1, qTsH[0], 1, qTsH[1] ); 
      if (1 && ADbgPfx) { // trace 'Inp-Med', 'Inp-Ave'
        BS_SubS( cN, qTsH[2], qTsH[0] ); BS_SubS( cN, qTsH[2], qTsH[1] ); BS_SubA( cN, pFrm, qTsH[2] );
        strcpy( lDbgSfx, "5BThr" ); BS_DatasWrt( ADbgPfx, cN, qTsH, 0, 3 ); };
      qTM.Mark( "des_Thr" ); };
    if (des_RadR) { // Radon's despeckle (it traces "6*.raw")
      if (des_KSnP>1e-3f) pSnP = *qTsH++; pLoF = *qTsH++; 
      iScrapeGrid_CPU( AW, AH, ADbgPfx, lDbgSfx, pFrm,
        qTsH, des_Kgr ); // emphasize high_freqs's grid
      iScrapenPFetch_CPU( AW, AH, ADbgPfx, lDbgSfx, diag_RadHF, pFrm, pSnP,
        qTsH, des_KnP, des_LimnP ); // fetch non-proc scrape
      iRadonDespeckle_CPU( AW, AH, pFrm, pLoF,
        qTsH, des_RadR, des_FlgR, des_KR ); // Dat=>Lap,Ave
      iScrapeStdDispatch_CPU( AW, AH, ADbgPfx, lDbgSfx, diag_RadHF, pLoF, pFrm,
        qTsH, des_KdSS ); // dispatch standard scrape
      qTM.Mark( "des_Radon" ); };
    if (des_RadN) { // Radon|Shannon-denoise for mid_freq laplacians
      BS_Noise( des_RadN, AW, AH, pFrm, qTsH[0], qTsH+1, des_LgN );
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Anoi.raw" ); BS_Write( ADbgPfx, cN, qTsH[0] ); };
      //if (1) { strcpy( lDbgSfx, "7B" ); BS_DatasWrt( ADbgPfx, cN, qTsH+1, 0, 8 ); }; // rarely trace ????
      if (des_RadN2) { BS_Cpy( cN, qTsH[0], qTsH[1] );
        BS_MxLapStd( des_RadN2, AW, AH, qTsH[1], qTsH+2, 0.0f, des_ZN2, 0 );
        //if (1) { strcpy( lDbgSfx, "7D0" ); BS_DatasWrt( ADbgPfx, cN, qTsH+2, 0, 8 ); }; // rarely trace
        BS_VMul( cN, qTsH[1], des_KN2 );
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Cml.raw" ); BS_Write( ADbgPfx, cN, qTsH[1] ); };
        BS_SubA( cN, qTsH[1], qTsH[0] );
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Dnoi.raw" ); BS_Write( ADbgPfx, cN, qTsH[0] ); };
      }; // if des_RadN2
      BS_VMul( cN, qTsH[0], des_KN ); BS_Add( cN, qTsH[0], pFrm ); // ~"(Org+MxLap)/2"
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "7EMfq.raw" ); BS_Write( ADbgPfx, cN, pFrm ); };
      if (des_KU>1e-3f) { BS2_MxLapXLimbo( 1, des_KU, 0, 0, 0.0f, pFrm, qTsH+0, AW, AH ); // maxlap with extended limbo
        if (1 && ADbgPfx) {
          strcpy( lDbgSfx, "7FMfq.raw" ); BS_Write( ADbgPfx, cN, pFrm ); }; }; 
      qTM.Mark( "des_Noi" ); }; // if des_RadN
    if (des_RadR) {
      if (pLoF) BS2_Sum( des_Kmf, pFrm, 1.0f, pLoF );
      iScrapenPAdd_CPU( AW, AH, ADbgPfx, lDbgSfx, diag_RadHF, // add non_processed scrape
                        pFrm, pSnP, pLoF, qTsH, des_RadSnP, des_r0SnP, des_r1SnP, des_KSnP ); 
      if (pOrg) BS2_MimicDat( pOrg, pFrm, cN ); // to care about zero_bacground            
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "8XlapA0.raw" ); iScrapeDbgTrc_CPU( ADbgPfx, pFrm, qTsH+0, diag_RadHF ); }; };
    BS_Cpy( cN, pFrm, AFrm );
  } else { BS3_Dim( AW, AH ); // == CUDA ==
    tBS_Val *pFrm = ATsD[0], **qTsH = ATsD+1; BC_Move( cN*cBS_Sz, AFrm, pFrm );
    tBS_Val *pOrg = NULL; if (des_RadR) { pOrg = *qTsH++; BC_Cpy( pFrm, pOrg ); }; // to care about zero_bacground 
    tBS_Val *pSnP = NULL; // non_processed part of fetched scrape
    tBS_Val *pLoF = NULL; // low freqs (moving average)
    if (fabs(des_Gamma-1)>1e-3f) { BS3_Gamma( des_Gamma, 1<<14, pFrm ); if (1 && ADbgPfx) { // balance_gamma
      strcpy( lDbgSfx, "4BInpGamma.raw" ); BC_Write( ADbgPfx, cN, pFrm ); }; qTM.Mark( "des_Gamma" ); };
    if (des_Thr>0) { // initial despeckle by threshold
      if (1 && ADbgPfx) { // trace by diag_RadHF
        strcpy( lDbgSfx, "5AlapA0.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, pFrm, qTsH+0, diag_RadHF ); };
      BS6_DespeckleScalars( des_Thr, 0 ); // strong=0: means general_despeckle
      if (1 && ADbgPfx) BC_Cpy( pFrm, qTsH[2] );
      BS6_DespeckleCUDA( AW, AH, pFrm, 1, qTsH[0], 1, qTsH[1] );
      if (1 && ADbgPfx) { // trace 'Inp-Med', 'Inp-Ave'
        BC_SubS( qTsH[2], qTsH[0] ); BC_SubS( qTsH[2], qTsH[1] ); BC_SubA( pFrm, qTsH[2] );
        strcpy( lDbgSfx, "5BThr" ); BC_DatasWrt( ADbgPfx, qTsH, 0, 3 ); };
      qTM.Mark( "des_Thr" ); }; 
    if (des_RadR) { // Radon's despeckle (it traces "6*.raw")
      if (des_KSnP>1e-3f) pSnP = *qTsH++; pLoF = *qTsH++; 
      iScrapeGrid_CUDA( ADbgPfx, lDbgSfx, pFrm,
        qTsH, des_Kgr ); // emphasize high_freqs's grid
      iScrapenPFetch_CUDA( ADbgPfx, lDbgSfx, diag_RadHF, pFrm, pSnP,
        qTsH, des_KnP, des_LimnP ); // fetch non-proc scrape
      iRadonDespeckle_CUDA( pFrm, pLoF,
        qTsH, des_RadR, des_FlgR, des_KR ); // Dat=>Lap,Ave
      iScrapeStdDispatch_CUDA( ADbgPfx, lDbgSfx, diag_RadHF, pLoF, pFrm,
        qTsH, des_KdSS ); // dispatch standard scrape
      qTM.Mark( "des_Radon" ); }; 
    if (des_RadN) { // Radon|Shannon-denoise for mid_freq laplacians
      BC_Noise( pFrm, qTsH[0], qTsH+1, des_RadN, des_LgN );
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Anoi.raw" ); BC_Write( ADbgPfx, cN, qTsH[0] ); };
      if (des_RadN2) { BC_Cpy( qTsH[0], qTsH[1] );
        BC_MxLapStd( des_RadN2, qTsH[1], qTsH+2, 0.0f, des_ZN2, 0 ); 
        BC_VMul( qTsH[1], des_KN2 );
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Cml.raw" ); BC_Write( ADbgPfx, cN, qTsH[1] ); };
        //if (0 && ADbgPfx) { strcpy( lDbgSfx, "3Dq" ); BC_DatasWrt( ADbgPfx, qTsH+2, 0, 4 ); }; // rarely trace
        BC_SubA( qTsH[1], qTsH[0] );
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "7Dnoi.raw" ); BC_Write( ADbgPfx, cN, qTsH[0] ); };
      }; // if des_RadN2
      BC_VMul( qTsH[0], des_KN ); BC_Add( qTsH[0], pFrm ); // ~"(Org+MxLap)/2"
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "7EMfq.raw" ); BC_Write( ADbgPfx, cN, pFrm ); };
      if (des_KU>1e-3f) { BS3_MxLapXLimbo( 1, des_KU, 0, 0, 0.0f, pFrm, qTsH+0, ATsH[0] ); // maxlap with extended limbo
        if (1 && ADbgPfx) { // strcpy( lDbgSfx, "3Gatmp" ); BC_DatasWrt( ADbgPfx, qTsH+0, 0, 2 );
          strcpy( lDbgSfx, "7FMfq.raw" ); BC_Write( ADbgPfx, cN, pFrm ); }; };
      qTM.Mark( "des_Noi" ); }; // if des_RadN
    if (des_RadR) {
      if (pLoF) BS3_Sum( des_Kmf, pFrm, 1.0f, pLoF );
      iScrapenPAdd_CUDA( ADbgPfx, lDbgSfx, diag_RadHF, // add non_processed scrape
                         pFrm, pSnP, pLoF, qTsH, des_RadSnP, des_r0SnP, des_r1SnP, des_KSnP );
      if (pOrg) BS3_MimicDat( pOrg, pFrm ); // to care about zero_bacground 
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "8XlapA0.raw" ); iScrapeDbgTrc_CUDA( ADbgPfx, pFrm, qTsH+0, diag_RadHF ); }; };
    BC_Move( -cN*cBS_Sz, AFrm, pFrm );
  }; // of CUDA
}
