// I16_Ut5.cpp : unsharpen's core for CPU|CUDA 
// frequent change needs 'synchronization|equivalence' between CPU & CUDA
// by S.B.
// cooperative debug with another data-series might require (in BS5_UnSharp_CPU) code
// ~"int qOk = BS_Read( "D:\\SimC\\VPF\\I16_Tst\\q\\$a11_0000\\GC0.raw", cN*2, pH[2] );"

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h" // common routines (OS,EABI-specific)
#include "BS_Ut0.h" // CPU basic-routines
#include "BS_Ut1.h" // CUDA basic-routines (~interconnect with kernel.cu)
#include "BS_Ut2.h" // CPU aplication-routines with standard dimension
#include "BS_Ut2a.h" // CPU aplication-routines addition
#include "BS_Ut3.h" // CUDA mid-level routines
#include "BS_Ut4.h" // beam-hardening, CLAHE (mostly CPU|CUDA is a parameter) 
#include "BS_Ut5.h"
#include "BS_Ut6.h" // despeckle, thin_laps (CPU & CUDA with explicit dimension)
#include "BS_Ut7.h" // ad-hoc service CPU-transforms
#include "BS_Ut5a.h" // "air-estimation"
#include "BS_Ut5c.h" // "quad-laps"
#include "BS_Ut5d.h" // "extended mask-recount"

static CBS_Ticks ixTM; // ticks' mark
static void iW( const char *AMsg ) { // error-tags here are ~informal (don't affect behavior)
  char qS[256]; sprintf( qS, cTPS" %.250s", AMsg ); WriteLn( qS ); }

//// adjusts control-record, fills transient parameters ////
void BS5_UnSharpRecAdjust( int AW, int AH, const tUnSharp *AUSOrg, const tUnSharpAdjs *AUSAdjs,
  int AModeId, int AAdj, tUnSharp *AUS, tUnSharpExts *AUSExts ) { // modify tUnSharp,tUnSharpExts
  char        qS[256]; memcpy( AUS, AUSOrg, sizeof(tUnSharp) );
  bool        qHaveBCS = false; int qBCS[3];
  for (int qi0=0;qi0<3;qi0++) { int qT = (AAdj>>((2-qi0)*4))&15;
    if (!qT) qT = 8; qT -= 8; qBCS[qi0] = qT; if (qT) qHaveBCS = true; };
  if (!qHaveBCS) {
    sprintf( qS, "USAdj(%Xh): none", AAdj ); iW (qS );
  } else { int qHaveChg = 0; // if HaveBCS
    const int cN = AUSAdjs->rN; for (int qn=0;qn<cN;qn++) {
      const tUnSharpAEn *cE = AUSAdjs->rE+qn;
      int     qBCSi = cE->rBCSi; if (qBCSi<0 || qBCSi>=3) {
        sprintf( qS, "USAdj ERR: #%d.BCSi=%d", qn, qBCSi ); iW( qS ); continue; };
      int     qChg = qBCS[qBCSi]; if (!qChg) continue;
      float   qMul = qChg<0 ? mcIpl(-qChg, 7, 1.0f, cE->rKmn )
                            : mcIpl( qChg, 7, 1.0f, cE->rKmx );
      char   *qAddr = (char*)AUS+cE->rAddrOffs;
      if (cE->rIsFloat) { // interpolate a float
        float *qPtr = (float*)qAddr; float qO = *qPtr, qV = qO*qMul; *qPtr = qV;
        sprintf( qS, "us#%02d %g=>%g", qn, qO, qV ); iW( qS );
      } else { // interpolate an int
        int *qPtr = (int*)qAddr; int qO = *qPtr, qV = (int)(qO*qMul); *qPtr = qV;
        sprintf( qS, "us#%02d %d=>%d", qn, qO, qV ); iW( qS );
      }; 
      qHaveChg++;
    };
    sprintf( qS, "USAdj(%Xh): %d changes", AAdj, qHaveChg ); iW( qS );
  }; // if HaveBCS
  memset( AUSExts, 0, sizeof(tUnSharpExts) ); // initialize tUnSharpExts
  AUSExts->rAirEn = BS5a_AirEstimationEnabled(); // con-coursing "air-estimation"
  #ifdef I16_UT5C_H_INCLUDED
  AUSExts->rQuadEn = BS5c_QuadLapsEnabled(); // con-coursing "quad-laps"
  #endif
  AUSExts->rMskEn = BS5d_MskEnabled(); // con-coursing "extended mask-recount"
  AUSExts->rAdj = AAdj; AUSExts->rModeId = AModeId;
  // assign introral bevels,corners (AUS.mgn[?])=0 if those are zeros
  switch (AModeId) {
    case 5: case 6:
      if (!AUS->rMgn[1]) { AUS->rMgn[1] = (AW*345+1619)/1620; // 282, 312 345 !!!
        sprintf( qS, "USAdj: mgnBev:=%d", AUS->rMgn[1] ); iW( qS );
      }; break;
  }; } // end of BS5_UnSharpRecAdjust

//// reconsider tUnSharp in CPU|CUDA-wrapper before call main unsharp ////
static void iiUnsharpRecDnSampleRad( int *AUSParRad, float AKs, bool ADimDnEn ) {
  int qRad = *AUSParRad; // optional 'BS_DimDn' could be quite non-liner ...
  int qRadMod = qRad/100, qRadFraq = qRad-qRadMod*100; 
  qRadFraq = (int)(qRadFraq*AKs+0.4999f); // Radius *= K_spatial_dilation
  if (ADimDnEn && !qRadMod && qRadFraq>1) { qRadMod = 1; qRadFraq >>= 1; }; 
  *AUSParRad = (qRadMod*100)|qRadFraq; }
static void jUnSharpRecReConsider( int AW, int AH, tUnSharp *AUS, tUnSharpExts *AUSExts ) {
  char qS[256]; // consider and change US^ for a) "spatial downsample" b)"binning"
  if (AUSExts->rModeId<5 || AUSExts->rModeId>9) return;
  const float a_DnSpDi = AUS->fA[4]; if (a_DnSpDi<1e-3f) return; // 
  const bool  cAccuFlg = (AUSExts->rAdj & (1<<12))!=0, // accuracy flag "1BCS"
              cFrm1x1 = AW*AH>999*1799; // no binning
  if (cAccuFlg==cFrm1x1) return;
  if (cFrm1x1) AUSExts->rDnSample = true; // physical spatial downsample
  // downsample-related factors
  AUS->rMgn[0] >>= 1; // mgn_Msk
  // note: 'mgn_Bev' Mgn#1 doesn't participate in unsharp
  float     cKs = a_DnSpDi, // spatial dilation
            cKm = 1.0f/(float)pow(a_DnSpDi,0.25f);  // magnitudes amplification
  sprintf( qS, "US DnSmpl=%d fact_spa=%.2f mag=%.2f by accu:%d %s %04dx%04d",
    (int)AUSExts->rDnSample, cKs, cKm, (int)cAccuFlg, cFrm1x1?"1x1":"2x2", AW, AH ); iW( qS );
  AUS->rA[0] = (AUS->rA[0]>>1)+(AUS->rA[0]&1);
  AUS->fA[0] *= cKm; // a_KHF
  iiUnsharpRecDnSampleRad( AUS->rC+0, cKs, false ); // c_RadBas
  iiUnsharpRecDnSampleRad( AUS->rC+1, cKs, true  ); // c_RadBig
  AUS->fC[0] *= cKm; // c_RefBig
  AUS->rD[0] -= 1; // d_Blk d_Blr  //AUS->rD[1] >>= 1;
  AUS->rE[0] -= 1; // e_Blk e_Blr  //AUS->rE[1] >>= 1; 
  AUS->rJ[0] = (int)(AUS->rJ[0]*cKs); // j_RadMM
  AUS->fL[0] *= (float)pow(cKs,0.95f); // l_KEdg reduce magnitude of adjust ~non-robust hessian-radius
  AUS->fSA[2] *= cKm*1.25f; // s_Za#2
  // AUS->fSB[0] *= cKm; // s_Zb#0
  // AUS->fSB[1] *= cKm*1.0f; // s_Zb#1
  } // end of jUnSharpRecReConsider

//// miscelanneous commons ////
static const int
  cNrmGamLg = 14, cNrmGam = 1<<cNrmGamLg, // normalized 'gamut=1" 
  cNrmAveLg =  9, cNrmAve = 1<<cNrmAveLg; // normalized average="some_standard_to compare"

static float iiNrmAve( int AN, tBS_Val *ADat, tBS_Val *ATmpAbs=NULL, int ANrmAve=cNrmAve ) {
  if (ATmpAbs) { BS_Cpy( AN, ADat, ATmpAbs ); BS_Abs( AN, ATmpAbs ); } else ATmpAbs = ADat; // 'absolute average' or 'plain signed avegare'
  float qA = mcAveLL( AN, ATmpAbs ); BS_VMul( AN, ADat, ANrmAve/(qA+0.5f) ); return qA; }

static const int cQRad = 4, cQRadO = 1; // anisotropic blur of "quad-laps" in _RadonMM

static int iiBodyOne_Rblr( float ABodyOne, int AW, int AH ) {
  int qRblr = (int)((tLL)(20*ABodyOne)*AW*AH>>20); return mcRng( qRblr, 1, 8 ); }

//// run-time hard-coding changes of control parameters by initially known stats ////
//  Init: E.SttNrmDelBig := MetricNrmDel(c_RefBig,cNrmDelLapBig,E.SttAveBig);
//        where "c_RefBig" is from INI-file, "E.SttAveBig" calculates at "stage C".
//  Usage: "XYZ_eff := iiByNrmDel(XYZ_ini,Exts.SttNrmDelBig,-0.3,0.3)".
static float iiMetricNrmDel( float ARefNrm, float ANrmDel, float AMetric ) { // relative metric's delta (against a 'norm=1')
  if (ARefNrm<1e-3f) return 0; float qDel = AMetric-ARefNrm; if (fabs(qDel)<1e-3f) return 0;
  qDel /= ARefNrm*ANrmDel; return mcRng( qDel, -1, 1 ); }
template<typename tN> static tN iiByNrmDel( float ADel, tN AV, float ADel1Neg, float ADel1Pos ) { // metric_eff from metric_ini
  if (ADel<-1e-3f) return AV+(tN)(AV*ADel1Neg*(-ADel)); if (ADel>1e-3f) return AV+(tN)(AV*ADel1Pos*ADel); return AV; }

//// helpers for _UnSharp_... ////
static int iiHistCutOut( tBS2_HistCut *AHistCut, // fills control-record for ~final "output histogram's cut"
  float AEstimCut0, float AEstimCut1, const tUnSharp *AUS, bool ADbgTrc ) {
  const float // pixels ('~<threshold(rhfCut0)' become downgraded negatives (smaller magnitude)
    u_BlOuDG = AUS->fU[0], // BlackOutDnGray ~0.2f ) downgrades output negative gray
    u_Air0Mn = AUS->fU[1], // DetectedAir0Mn ~3E-3f
    u_Air1Mn = AUS->fU[2], // DetectedAir1Mn ~1E-5f
    cBlackInpGapNeg = 0, // this gamut's centinel should activate quite infrequently;
    cBlackInpGapPos = 0.5f; // it limits histogram's length (to process a bit faster)
  const int cBlackOutGapNeg = (int)(cNrmGam*u_BlOuDG*0.15f), // *'~0.1*BlOuDG'
            cWhiteOutGapPos = (int)(cNrmGam*0.04f), // for ~'add HiF in unsharp wrapper'
            mgn_Msk = AUS->rMgn[0]; // ignore margins when build histogram
  tBS2_HistCut cHistCut = {
  {-(int)(cNrmGam*cBlackInpGapNeg), // Black
    (int)(cNrmGam*(1+cBlackInpGapPos)), // White
    cNrmGam-cBlackOutGapNeg-cWhiteOutGapPos, // GaO
    abs(mgn_Msk) }, // |mgn_Msk| => HistCut.n[Mgn]
  { 0, 0, // fCut0 fCut1
    u_BlOuDG, 0 } }; // fTD0 fTD1 
  memcpy( AHistCut, &cHistCut, sizeof(tBS2_HistCut) );
  AHistCut->f[cBS2_rhfCut0] = __max(AEstimCut0,u_Air0Mn);
  AHistCut->f[cBS2_rhfCut1] = __max(AEstimCut1,u_Air1Mn);
  if (ADbgTrc) { char qS[256]; sprintf( qS, "HisCutOut:=%d %d..%d=>%d %.3f,%.3f",
    cBlackOutGapNeg, AHistCut->n[cBS2_rhcBlack], AHistCut->n[cBS2_rhcWhite],
    AHistCut->n[cBS2_rhcGaO], AHistCut->f[cBS2_rhfCut0], AHistCut->f[cBS2_rhfCut1] ); iW( qS ); };
  return cBlackOutGapNeg; } // returns cBlackOutGapNeg to increment output of 'HistCut(...)'

static int iiHaveNZ( int AL, const float *ACoeffs ) { // -1:have_neg, 0:all_zeros, 1:have_pos
  bool qHavePos = false; for (int qi=0;qi<AL;qi++) { float qV = ACoeffs[qi]; if (qV<-1e-3f) return -1;
   if (qV>1e-3f) qHavePos = true; }; return qHavePos?1:0; }

//// UnSharp CPU ////
static void iUnSharp_CPU( const tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH,
  tBS_Val *AHiF, tBS_Val *AHiFblr, tBS_Val **ATs,
  const char *ADbgPfx, char *lDbgSfx ) { const int cN = AW*AH; char qS[256];
  tBS_Val     **qStack; // ~last used frame
  tBS_Val      *qLapAcc = *ATs++; BS_ZFil( cN, qLapAcc ); float qKLapAcc = 0; // ad-hoc laplacians' accumulator
  if (AUSExts->rAirEn) { // start concourse "air-estimation"
    if (ADbgPfx) sprintf( lDbgSfx, "QA" ); BS5a_AirEstimation( AW, AH, AHiFblr, ADbgPfx, AUSExts->rAdj ); };
  // - stage B "gamma-correction"; #0:cGammaSpa #1:cGammaCLA #2:cGammaAtt #3:|Der3x3|
  tBS_Val     **pB = ATs,
               *cGammaSpa = pB[0],   // for spatial metrics (a gamma-corrected frame)
               *cGammaCLA = pB[1],   // for CLAHE
               *cGammaAtt = pB[2],   // for common attenuation
               *cGSDer3x3 = pB[3]; { // global shape |Der3x3|
    const float b_GammaSpa = AUS->fB[0], b_GammaCLA = AUS->fB[1], b_GammaAtt = AUS->fB[2];
    BS_Cpy( cN, AHiF,    cGammaSpa ); if (fabsf(b_GammaSpa-1)>1e-3f) BS2_Gamma( b_GammaSpa, cNrmGam, cGammaSpa );
    BS_Cpy( cN, AHiFblr, cGammaCLA ); if (fabsf(b_GammaCLA-1)>1e-3f) BS2_Gamma( b_GammaCLA, cNrmGam, cGammaCLA );
    BS_Cpy( cN, AHiFblr, cGammaAtt ); if (fabsf(b_GammaAtt-1)>1e-3f) BS2_Gamma( b_GammaAtt, cNrmGam, cGammaAtt );
    BS_Cpy( cN, cGammaCLA, cGSDer3x3 ); BS_DrA3x3( AW, AH, cGSDer3x3, pB+4 ); iiNrmAve( cN, cGSDer3x3 );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "B" ); BS_DatasWrt( ADbgPfx, cN, pB, 0, 4 ); };
    ixTM.Mark( "B gamma-correction" ); };
  // - stage C "basic spatial metrics"; #0:BlrBas #1:BlrBig #2:|LapBig|_nrm #3:LapBas #4:LapBig
  tBS_Val     **pC = pB+4,
               *cGSLapBig = pC[2]; { // global shape |LapBig|_nrm
    const int   c_RadBas = AUS->rC[0], c_RadBig = AUS->rC[1], c_RadMM = AUS->rC[2]; // ~4, 205 or 109
    const float c_RefBig = AUS->fC[0], c_DelBig = AUS->fC[1], c_KMM = AUS->fC[2]; // reference norm=1 of |LapBig|
    BS_Cpy( cN, cGammaSpa, pC[0] ); BS2_MovAve2Dp( c_RadBas%100, pC[0], pC+5, c_RadBas/100 ); // #0:BlrBas
    BS_Cpy( cN, cGammaSpa, pC[3] ); BS_SubA( cN, pC[0], pC[3] ); // #3:LapBas
    if ( 0 && c_RadMM>0 && c_KMM>-1e-3f) {
      BS_MxLapStd( c_RadMM, AW, AH, pC[3], pC+5, c_KMM ); // #3:(1-KMM)*LapBas+KMM*MaxMag(LapBas)
      BS_Cpy( cN, pC[0], cGammaSpa ); BS_Add( cN, pC[3], cGammaSpa ); }; // update GammaSpa
    BS_Cpy( cN, cGammaSpa, pC[1] ); BS2_MovAve2Dp( c_RadBig%100, pC[1], pC+5, c_RadBig/100 ); // #1:BlrBig
    BS_Cpy( cN, cGammaSpa, pC[4] ); BS_SubA( cN, pC[1], pC[4] ); // #4:LapBig
    BS_Cpy( cN, pC[4], cGSLapBig ); BS_Abs( cN, cGSLapBig );
    AUSExts->rSttAveBig = iiNrmAve( cN, cGSLapBig ); // #2:|LapBig|_nrm
    AUSExts->rSttNrmDelBig = iiMetricNrmDel( c_RefBig, c_DelBig, AUSExts->rSttAveBig );
    if (1 && ADbgPfx) { // report UnSharpExts.AveBig, relative <against reference> |LapBig|
      sprintf( qS, "USstC ave|Big|=%.2f NrmDelBig=%.2f", AUSExts->rSttAveBig, AUSExts->rSttNrmDelBig ); iW( qS );
      strcpy( lDbgSfx, "C" ); BS_DatasWrt( ADbgPfx, cN, pC, 0, 5 ); };
    ixTM.Mark( "C basic spatial metrics" ); };
  // - stage D "CLAHE"; #0:plain_CLAHE #1:mskBodyOne #2:CLAInc #3:|Der3x3CLA|
  tBS_Val     **pD = pC+5,
               *cGSDer3x3CLA = pD[3], // global shape |Der3x3CLA|
               *qmskBodyOne = NULL; {
    const int   d_Blk = AUS->rD[0], d_Blr = AUS->rD[1], // ~6,9
                d_IncRadBlr = 1;
    const float d_BodyOne = AUS->fD[0], d_ExtMn = AUS->fD[1], d_ExtMx = AUS->fD[2];
    AUSExts->rSttCLABlr = d_Blr<=0 ? d_Blr : iiByNrmDel( AUSExts->rSttNrmDelBig, d_Blr, -0.1f, 0.25f ); // actual d_Blr
    tBS4_CLAHE_Exts qCLAExts = {d_ExtMn,d_ExtMx}; if (ADbgPfx) strcpy( lDbgSfx, "Da" );
    BS4_CLAHE_Transform( // #0:plain_CLAHE #1:BodyOne
      d_Blk, cNrmGam, AUSExts->rSttCLABlr, cGammaCLA, pD, false, ADbgPfx, &qCLAExts );
    if (0 && ADbgPfx) { strcpy( lDbgSfx, "Db" ); BS_DatasWrt( ADbgPfx, cN, pD, 0, 2 ); };
    BS_SubA( cN, cGammaCLA, pD[0] ); // #0:CLAInc:=CLAHE(Gamma1)-Gamma1
    BS2_MovAve2D( d_IncRadBlr, pD[0], pD[2] ); // #2:CLAInc blurred
    BS_Cpy( cN, cGammaCLA, pD[0] ); BS_Add( cN, pD[2], pD[0] ); // #0:plain_CLAHE
    if (d_BodyOne>1e-3f) { qmskBodyOne = pD[1]; // #1:BodyOne 
      BS_VMul( cN, qmskBodyOne, cNrmAve*1.2f/((1<<cBS4_XLATGamLg)*d_BodyOne) );
      BS_Rng( cN, qmskBodyOne, -999999, cNrmAve-1 );
      BS2_MovAve2Dp( iiBodyOne_Rblr(d_BodyOne,AW,AH), qmskBodyOne, pD+4, 2 );
      BS_Mul( cN, qmskBodyOne, pD[2], -cNrmAveLg ); }; // attenuate CLAInc
    BS_Cpy( cN, pD[0], cGSDer3x3CLA ); BS_DrA3x3( AW, AH, cGSDer3x3CLA, pD+4 ); iiNrmAve( cN, cGSDer3x3CLA );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Dc" ); BS_DatasWrt( ADbgPfx, cN, pD, 0, 4 ); };
    ixTM.Mark( "D CLAHE" ); qStack = pD+4; };
  // - stage EA "to attenuate mskGDi,mskEdg; #0:mskAtt
  tBS_Val      *qmskAtt = NULL; const int e_Blk = AUS->rE[0]; if (e_Blk>1) {
    const int   e_Blr = AUS->rE[1], e_Lim = AUS->rE[2], e_Zero = AUS->rE[3]; // Lim is against e_Blk=6
    tBS_Val   **pEA = qStack, *qEt = pEA[3]; BS_Cpy( cN, cGammaAtt, qEt );
    BS4_FrmNrm_CPU( cNrmGam, qEt, pEA ); qmskAtt = pEA[0]; // #0:normalized GammaAtt
    BS4_CLAHE_TraAveRef( e_Blk, e_Blr, qEt, pEA, false, (tBS_Val*)-1 ); // #0:hist_counts of #t 
    if (0 && ADbgPfx) { strcpy( lDbgSfx, "EAb_His.raw" ); BS_Write( ADbgPfx, cN, qmskAtt ); };
    const int csh = (e_Blk-6)<<1, cLim = csh>0 ? e_Lim<<csh : e_Lim>>-csh; BS_Neg( cN, qmskAtt );
    if (0 && csh) { sprintf( qS, "USstEA Blk=%d Lim %d=>%d", e_Blk, e_Lim, cLim ); iW( qS ); };
    BS_VAdd( cN, qmskAtt, cLim ); BS_Rng( cN, qmskAtt, 0 ); BS_VMul( cN, qmskAtt, (float)cNrmAve/cLim );
    if (e_Zero) BS_VAdd( cN, qmskAtt, -e_Zero );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "EAc_Att.raw" ); BS_Write( ADbgPfx, cN, qmskAtt ); };
    ixTM.Mark( "EA MskAtt" ); qStack = pEA+1;
  }; // if eBlk
  // - stage F "the masks"; #0:mskThn #1:mskThk #2:mskGDi #3:.. 
  tBS_Val     **pF = qStack, // the reference: _ReadMe.txt, "Generation of normalized masks at stage F"  
               *qmskThn = pF[0], // ~'Thin'; merge LapAda,LapMM to LapAdaEff
               *qmskThk = pF[1], // ~'Thick'; merge LapBas,LapBig to LapAda
               *qmskGDi = pF[2], // ~'General dilation'; dilate HiFrem,LapAdaEff
               *qmskEdg = pF[3]; // ~ a bump at lap's zero-izoline
  bool          qBigGThn_blr; {  // true when ~blur LapMM should go rather to big GThn
    float f_Att = AUS->fF[0], f_EmLthn = AUS->fF[1], f_EmLedg = AUS->fF[2], // ~ 6.0 0.1 0.05
          f_MxLbias = AUS->fF[3]; // ~ 1.15
    float f_K[3][3]; memcpy( f_K, AUS->fFK, sizeof(AUS->fFK) ); // general shapes from basic stats
    int   f_C[4][5]; memcpy( f_C, AUS->rFC, sizeof(f_C) ); // normalized masks by ~shapes
    qBigGThn_blr = f_C[0][0]>f_C[0][1]; if (1) { sprintf( qS, "USstF BigGThn_blr=%d rcuBlack,White={%d,%d}",
      (int)qBigGThn_blr, f_C[0][0], f_C[0][1] ); iW( qS ); };
    int   f_RadThk = f_C[1][4], // if (f_C[1][0]>f_C[1][1]) AUSExts->rSttFlg |= 4; // f_Cb[0]>f_Cb[1]
          qf_Khn = 0; // f_K[k] has non-zero coefs
    tBS_Val **qFd = pF, **qFb = qFd+4, **qFc = qFb+mcASz(f_K), **qT = qFc+3;
    { // <shapes> @b0:~GSh @b1:~GDi @b2:~GThn
      #if 0 // ???? temp
        BS_Cpy( cN, cGSDer3x3, qT[0] ); BS2_MovAve2Dp( 3, qT[0], qT+1, 1 ); BS_SubA( cN, qT[0], cGSDer3x3 );
        BS_Cpy( cN, cGSDer3x3CLA, qT[0] ); BS2_MovAve2Dp( 4, qT[0], qT+1, 1 ); BS_SubA( cN, qT[0], cGSDer3x3CLA );
      #endif
      if (1 && ADbgPfx) { tBS_Val *cTs[] = {cGSDer3x3,cGSDer3x3CLA,cGSLapBig};      // optional summary
        strcpy( lDbgSfx, "Fa" ); BS_DatasWrt( ADbgPfx, cN, cTs, 0, mcASz(cTs) ); }; // trace of inp_shapes
      for (int qk=0;qk<mcASz(f_K);qk++) { // K0*|Der3x3|+K1*|Der3x3CLA|+K2*|LapBig| (all normalized)
        int qHaveNZ = iiHaveNZ( 3, f_K[qk] ); if (!qHaveNZ) continue;
        qf_Khn |= 1<<qk; BS_Cpy( cN, cGSDer3x3, qFb[qk] ); // B#3=|Der3x3|
        BS2_Sum( f_K[qk][0], qFb[qk], f_K[qk][1], cGSDer3x3CLA, f_K[qk][2], cGSLapBig );
        iiNrmAve( cN, qFb[qk],qHaveNZ<0?qT[0]:NULL ); };
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Fb" ); BS_DatasWrt( ADbgPfx, cN, qFb, qf_Khn, 0 ); };
      ixTM.Mark( "Fb shapes" ); };
    { // <hard-coded logic of spatial filters> @c0:MskThk @c1:MskGDi @c2:MskEdg
      BS_Cpy( cN, qFb[0], qFc[0] ); BS2_MovAve2Dp( f_RadThk%100, qFc[0], qT, f_RadThk/100 ); // @c0:MskThk
      BS_Cpy( cN, qFc[0], qFc[2] ); BS_SubS( cN, qFb[0], qFc[2] ); // @c2:MskEdg
      if (!(qf_Khn&(1<<1))) { BS_Cpy( cN, qFc[2], qFb[1] ); BS_Abs( cN, qFb[1] ); // if no Kb then GDi:=abs(MskEdg)
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "FB1.raw" ); BS_Write( ADbgPfx, cN, qFb[1] ); }; };
      if (!(qf_Khn&(1<<2))) { BS_Cpy( cN, qFc[2], qFb[2] ); BS_Neg( cN, qFb[2] ); // if no Kc then GThn:=-MskEdg
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "FB2.raw" ); BS_Write( ADbgPfx, cN, qFb[2] ); }; };
      BS_Cpy( cN, qFb[1], qFc[1] ); // @c1:MskGdi:=GDi+MskAtt*KGDi
      if (qmskAtt) { BS_VMul( cN, qmskAtt, f_Att ); BS_Add( cN, qmskAtt, qFc[1] ); };
      if (1 && ADbgPfx) { // index+semantic in this trace_dump
        strcpy( lDbgSfx, "Fc0Thk.raw" ); BS_Write( ADbgPfx, cN, qFc[0] );
        strcpy( lDbgSfx, "Fc1GDi.raw" ); BS_Write( ADbgPfx, cN, qFc[1] );
        strcpy( lDbgSfx, "Fc2Edg.raw" ); BS_Write( ADbgPfx, cN, qFc[2] ); };
      ixTM.Mark( "Fc Masks" ); };
    { // cut gamut's of denormalized MskXYZ to normalized mskXYZ (and perhaps blur)
      BS_Cpy( cN, qFb[2], qmskThn ); BS_Cpy( cN, qFc[0], qmskThk );
      BS_Cpy( cN, qFc[1], qmskGDi ); BS_Cpy( cN, qFc[2], qmskEdg );
      int   qf_Chnzr = 0; for (int qc=0;qc<4;qc++) { // f_C[c] cut has nonzero range
        if (f_C[qc][0]==f_C[qc][1]) continue; qf_Chnzr |= 1<<qc; // masks within 0..NrmAve
        int cTg0 = f_C[qc][2], cGaO = cNrmAve-cTg0-f_C[qc][3];
        if (cGaO>2) { tBS2_Cut cCut = { {f_C[qc][0],f_C[qc][1],cGaO}, {0,0} };
          BS2_Cut( &cCut, qFd[qc] ); if (cTg0) BS_VAdd( cN, qFd[qc], cTg0 ); }; // if GaO
        int qRad = f_C[qc][4]; BS2_MovAve2Dp( qRad%100, qFd[qc], qT, qRad/100 ); }; // for c
      if (!(qf_Chnzr&(1<<2))) BS_Cpy( cN, qFd[1], qFd[2] ); // if no Cc then mskGDi := mskThk
      if (f_MxLbias>1e-3f) { // mxlap for mskGdi, mskEdg
        int qRadThk = (f_RadThk%100)<<(f_RadThk/100); if (qRadThk>8) qRadThk = 8;
        BS_Med( 2, AW, AH, qmskGDi, qT[0] ); BS_Add( cN, qT[0], qmskGDi ); mcShiftRgt( cN, qmskGDi, 1 );
        BS_MxLapStd( qRadThk+2, AW, AH, qmskGDi, qT, f_MxLbias, 0, -1, 0 );
        if (f_MxLbias>1) BS_Rng( cN, qmskGDi, -99999, cNrmAve-1 );
        BS_Mul( cN, qmskEdg, qmskEdg, -cNrmAveLg ); // <possibly> square edge_mask
        BS_MxLapStd( qRadThk, AW, AH, qmskEdg, qT, f_MxLbias, 0, -1, 0 );
        if (f_MxLbias>1) BS_Rng( cN, qmskEdg, -99999, cNrmAve-1 ); };
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Fd" ); BS_DatasWrt( ADbgPfx, cN, qFd, qf_Chnzr, 0 ); };
      ixTM.Mark( "Fd norm masks" ); };
    if (BS2_AccMim( qLapAcc, qT[0], (float)fabs(f_EmLthn), qmskThn, f_EmLthn<0?pC[3]:pC[4],
                                    (float)fabs(f_EmLedg), qmskEdg, f_EmLedg<0?pC[3]:pC[4] )) { // Edge mimic Laplacian
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "FeEmL.raw" ); BS_Write( ADbgPfx, cN, qLapAcc ); };
      ixTM.Mark( "Ff acc EmL" ); qKLapAcc = 1; };
    if (AUSExts->rMskEn) { if (ADbgPfx) strcpy( lDbgSfx, "QC" ); // start concourse "extended mask-recount"
      BS5d_Msk( AW, AH, false, ADbgPfx, pC[1], qmskThk, qmskGDi ); };
    qStack = qFd+4; }; // last used frame
  // - stage H "LapAdaptive"; #0:LapMid #1:BlrAda #2:LapAda
  // Merge LapBas,LapBig to LapAda by mskThk; i.e. attenuate (LapBas-LapBig).
  tBS_Val     **pH = qStack; {
    const int   h_Lim = AUS->rH[0], h_LimP = AUS->rH[1]; // truncate LapAda
    BS_Cpy( cN, pC[0], pH[0] ); BS_SubA( cN, pC[1], pH[0] ); // #0:LapMid:=BlrBas-BlrBig (also '=LapBig-LapBas')
    BS_Cpy( cN, pH[0], pH[1] ); BS_Mul( cN, qmskThk, pH[1], -cNrmAveLg ); // #1:LapMid_dil:=mskThk*LapMid
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Ha" ); BS_DatasWrt( ADbgPfx, cN, pH, 0, 2 ); };
    BS_Add( cN, pC[1], pH[1] ); // #1:BlrAda:=BlrBig+LapMid_dil
    BS_Cpy( cN, cGammaSpa, pH[2] ); BS_SubA( cN, pH[1], pH[2] ); // #2:LapAda:=GammaSpa-BlrAda
    BS_RngA( cN, pH[2], h_Lim, h_LimP );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Hb" ); BS_DatasWrt( ADbgPfx, cN, pH, 1, 3 ); };
    ixTM.Mark( "H LapAdaptive" ); qStack = pH+3; };
  // - stage J "Radon & LapAdaEff"; #0:LapMM #1:LapAdaEff
  // Get LapMM, then merge LapAda,LapMM to LapAdaEff by mskThn; i.e. attenuate (LapAda-LapMM).
  tBS_Val     **pJ = qStack; {
    const int   j_RadMM = AUS->rJ[0], j_RadMMO = AUS->rJ[1], j_RadENm = AUS->rJ[2], j_Lim = AUS->rJ[3];
    const float j_KMM = AUS->fJ[0], j_KENm = AUS->fJ[1];
    tBS_Val   **qqQuadLaps = NULL;
    #ifdef I16_UT5C_H_INCLUDED
    tBS_Val    *qqbufs[cBS5c_QuadLaps]; if (AUSExts->rQuadEn)
      qqQuadLaps = BS5c_QuadLapsAlloc(AW,AH,qqbufs,j_RadMMO,cQRadO); // CPU-memory for "quad-laps"
    #endif
    const tBS2_RadonMM cBS2_RadonMM = {j_RadMM,j_RadMMO,cQRad,cQRadO,j_RadENm};
    tBS_Val    *qENm = j_KENm>1e-3f ? pJ[3] : NULL;
    BS2_RadonMM( &cBS2_RadonMM, pH[2], // Radon's denoise of ~LapAdaptive
      pJ[2], pJ[0], pJ[1], qENm, pJ+4, qqQuadLaps ); // Dev Blr Srp ENm
    BS_VMul( cN, pJ[0], j_KMM ); // #0:LapMM*KMM
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Ja" ); BS_DatasWrt( ADbgPfx, cN, pJ, 0, qENm?4:3 ); };
    BS_SubA( cN, pJ[0], pJ[1] ); BS_Mul( cN, qmskThn, pJ[1], -cNrmAveLg ); // #1:=(Srp-Blr)*=mskThn
    BS_Add( cN, pJ[0], pJ[1] ); // #1:LapAdaEff (+=Blr)
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jb" ); BS_DatasWrt( ADbgPfx, cN, pJ, 1, 2 ); };
    if (qENm) { // add ~limit_ENm_by_Dev
      BS2_RngLapByConductor( pJ[2], qENm ); BS_VMul( cN, qENm, j_KENm );
      BS_Add( cN, qENm, pJ[1] ); // #1:LapAdaEff += ENm*KENm
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jc" ); BS_DatasWrt( ADbgPfx, cN, pJ, (1<<1)|(1<<3), 0 ); }; };
    if (BS_RngA( cN, pJ[1], j_Lim ))
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jd_LapAdaEff_J1.raw" ); BS_Write( ADbgPfx, cN, pJ[1] ); };
    #ifdef I16_UT5C_H_INCLUDED
    if (AUSExts->rQuadEn) { // start concourse "quad-laps"; the 'LapOrg' could vary as ~...
      if (ADbgPfx) strcpy( lDbgSfx, "QB" ); // ... H#0:LapMid H#2:LapAda J#0:LapMM
      BS5c_QuadLaps( AW, AH, qqQuadLaps, false, ADbgPfx, pJ[1], qmskThk ); };
    #endif
    ixTM.Mark( "J RadonDenoise" ); qStack = pJ+2; };
  // - stage LA "edge"; #0:mskEdg_blr, LapAcc+=Edg
  { const int   la_RadHess = AUS->rL[0]; const float la_KEdg = AUS->fL[0]; if (la_RadHess && la_KEdg>1e-3f) {
      tBS_Val **pLA = qStack, *qEdg = pLA[0]; BS_Cpy( cN, qmskEdg, qEdg ); // ~hessian of 1st derivation from mskEdg
      BS2_MimicDat( pJ[0], qEdg ); // #0:lapedg=LapMM
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAa" ); BS_DatasWrt( ADbgPfx, cN, pLA, 0, 1 ); };
      BS2_Hess( abs(la_RadHess), pLA ); // #0:Hessian #1:LapEdgBlr
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAb" ); BS_DatasWrt( ADbgPfx, cN, pLA, 0, 3 ); };
      BS2_MimicDat( pLA[1], qEdg ); BS_VMul( cN, qEdg, la_KEdg ); // #0:Edg
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAcEdg.raw" ); BS_Write( ADbgPfx, cN, qEdg ); };
      BS_Add( cN, qEdg, qLapAcc ); qKLapAcc = 1;
      if (true) BS_Cpy( cN, pLA[2], qmskEdg ); // reassign by mskEdgBlr
      ixTM.Mark( "LA Edge" ); qStack = pLA+1; }; };
  // - stage LB "ridge"; LapAcc+=Rdg
  { const float lb_KRdg = AUS->fL[1]; if (lb_KRdg>1e-3f) {
      tBS_Val **pLB = qStack; BS_Cpy( cN, pJ[1], pLB[0] ); // could be either J#0,#1 (either LapMM,LapAdaEff) 
      BS2_RidgeByConductors( pLB[0], cGSLapBig, cGSDer3x3, cGSDer3x3CLA ); // ... |LapBig|_nrm |Der3x3|_nrm |Der3x3CLA|_nrm
      BS_VMul( cN, pLB[0], lb_KRdg ); BS_Add( cN, pLB[0], qLapAcc ); qKLapAcc = 1; 
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LBaRdg.raw" ); BS_Write( ADbgPfx, cN, pLB[0] ); };
      ixTM.Mark( "LB Ridge" ); }; };
  // - stage LC "maxmag"
  { bool        qlc_chg = false;
    const int   lc_RadMM = AUS->rL[1]; const float lc_KMM = AUS->fL[2]; if (lc_RadMM>0 && lc_KMM>-1e-3f) {
      BS_MxLapStd( lc_RadMM, AW, AH, qLapAcc, qStack, lc_KMM ); qlc_chg = true; };
    const int   lc_Lim = AUS->rL[2]; qlc_chg |= BS_RngA( cN, qLapAcc, lc_Lim );
    if (qlc_chg) { if (0 && ADbgPfx) { strcpy( lDbgSfx, "LCa_lapacc.raw" ); BS_Write( ADbgPfx, cN, qLapAcc ); };
      ixTM.Mark( "LC LapAcc refine" ); qKLapAcc = 1; }; };
  // - stage QA "wait for air-estimation"
  tAirEstimationWait qAEW; memset( &qAEW, 0, sizeof(tAirEstimationWait) ); if (AUSExts->rAirEn) {
    BS5a_AirEstimationWait( &qAEW ); tBS_Val **pQ = NULL;
    if (qAEW.rAiF[0]) {
      pQ = qStack; qStack += 1; BS_Cpy( cN, qAEW.rAiF[0], pQ[0] );
      if (1 && ADbgPfx ) { strcpy( lDbgSfx, "QAzBkgPtr.raw" ); BS_Write( ADbgPfx, cN, pQ[0] ); }; };
    for (int qt=0;qt<mcASz(qAEW.rAiF);qt++) if (qAEW.rAiF[qt]) delete qAEW.rAiF[qt]; // release ~'BkrPtr from air-detection'
    if (pQ) qAEW.rAiF[0] = pQ[0];
    AUSExts->rAirInfFlg = qAEW.rDbgInfFlg; };
  #ifdef I16_UT5C_H_INCLUDED
  // - stage QB "wait for quad-laps"
  if (AUSExts->rQuadEn) { tBS_Val *qT = BS5c_QuadLapsWait(); if (qT) {
    BS_Add( cN, qT, qLapAcc ); qKLapAcc = 1; delete qT; }; };
  #endif
  // - stage QC "wait for extended mask-recount"
  if (AUSExts->rMskEn) BS5d_MskWait();
  // - stage N "global dilation"; LapAdaEff*=mskGDiRes, HiF:=mskGDiRes
  { const float n_Thk = AUS->fN[0], n_GDi = AUS->fN[1], n_Edg = AUS->fN[2]; // ~0 0.4 0.2
    tBS_Val   **pN = qStack, *mskGDiRes = AHiF;
    BS2_MskSum( mskGDiRes, pN[0], cNrmAve, n_Thk, qmskThk, n_GDi, qmskGDi, n_Edg, qmskEdg );
    if (qAEW.rAiF[0]) BS_Mul( cN, qAEW.rAiF[0], mskGDiRes, -8 ); // mskGDiRes*=Air
    BS_Mul( cN, mskGDiRes, pJ[1], -cNrmAveLg ); // J#1:LapAdaEff_dil (by mskGDiRes)
    if (1 && ADbgPfx) { 
      strcpy( lDbgSfx, "Na_GDiRes_X1.raw" ); BS_Write( ADbgPfx, cN, mskGDiRes );
      strcpy( lDbgSfx, "Nb_LapAdaEff_J1.raw" ); BS_Write( ADbgPfx, cN, pJ[1] ); };
    ixTM.Mark( "N Dilation" ); };
  // - stage S "summa summarum"; #0:BodySum #1:LapSum #2:ResSum
  tBS_Val     **pS = qStack; { // s_Za2*HiFRem adds later
    const float s_Za0 = AUS->fSA[0], s_Za1 = AUS->fSA[1], s_Za3 = AUS->fSA[3], // #0:BodySum
                s_Zb0 = AUS->fSB[0], s_Zb1 = AUS->fSB[1], s_Zb2 = AUS->fSB[2]; // #1:LapSum
    BS_Cpy( cN, pH[1], pS[0] ); BS2_Sum( s_Za0, pS[0], s_Za1, cGammaCLA,
      1-s_Za0-s_Za1, AHiFblr, (float)fabs(s_Za3), pH[s_Za3<-1e-3f?1:0] );
    BS_Cpy( cN, pD[2], pS[1] ); BS2_Sum( s_Zb0, pS[1], s_Zb1, pH[2], s_Zb2, pJ[1], qKLapAcc, qLapAcc );
    // don't multiply ~"large_low_freqs of CLAInc in LapSum" on ~"GDiRes"
    if (qmskBodyOne) BS_Mul( cN, qmskBodyOne, pS[1], -cNrmAveLg ); // attenuate LapSum
    if (qAEW.rAiF[0]) BS_Mul( cN, qAEW.rAiF[0], pS[1], -8 ); // #1:LapSum*=Air
    BS_Cpy( cN, pS[0], pS[2] ); BS_AddRg( cN, pS[1], pS[2] ); // #2:ResSum
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "S" ); BS_DatasWrt( ADbgPfx, cN, pS, 0, 3 );
      if (qKLapAcc>1e-3f) { strcpy( lDbgSfx, "S_lapacc.raw" ); BS_Write( ADbgPfx, cN, qLapAcc ); } };
    ixTM.Mark( "S summa summarum" ); };
  // - stage U "histogram cut out"; #0:Res
  tBS_Val      **pU = pS+3; { tBS2_HistCut qHistCut;
    int          qGapNeg = iiHistCutOut( &qHistCut, qAEW.rAir[0], qAEW.rAir[1], AUS, ADbgPfx!=NULL );
    BS_Cpy( cN, pS[2], pU[0] ); BS2_HistCut( &qHistCut, pU[0], NULL, pU+1 );
    BS_VAdd( cN, pU[0], qGapNeg ); // compensate the most of expecting negatives
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "U" ); BS_DatasWrt( ADbgPfx, cN, pU, 0, 1 );
      sprintf( qS, "USstU CatchMM=%d,%d", xBS2_rhCatchMM[0], xBS2_rhCatchMM[1] ); iW( qS ); };
    ixTM.Mark( "U histogram cut out" ); };
  // - stage V "finalize"; #0:ResBkGM
  tBS_Val     **pV = pU+1; BS_Cpy( cN, pU[0], pV[0] );
  const int     mgn_Msk = AUS->rMgn[0]; if (mgn_Msk>0) { // blur surrounding frame's margins when '>0'
    BS_Cpy( cN, AHiFblr, pV[1] ); BS2_VMimicNZ( pV[1], 1<<8 );
    BS2_MovAve2Dp( 2, pV[1], pV+2, 1 ); BS_Mul( cN, pV[0], pV[1], -8 ); // BS_Beth( cN, NULL, pV[0], pV[1], -8 );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "V" ); BS_DatasWrt( ADbgPfx, cN, pV, 0, 2 ); };
    ixTM.Mark( "V finalize" ); };
  BS_Rng( cN, pV[0], 0 ); BS_Cpy( cN, pV[0], AHiFblr ); // non negative output data
  AUSExts->rStackFrmBuf = (int)(pV-ATs); } // end of iUnSharp_CPU

//// unsharp wrapper CPU ////
void BS5_UnSharp_CPU( tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs,
  const char *ADbgPfx, char *lDbgSfx ) { // -- UnSharp CPU (declaration for ~"no downsample")
  const int     cN = AW*AH; BS2_Dim( AW, AH ); char qS[256], qSCWD[40];
  if (xDbgDmp==1) strcpy( ixTM.rPfx, "US" ); ixTM.Reset(); // tick-marks
  jUnSharpRecReConsider( AW, AH, AUS, AUSExts );
  // - stage A "high freqs preliminaries"; #0:HiFrem #1:HiF #2:HiFblr
  tBS_Val     **pA = ATs, **qTs = pA+3; {
    int         a_RadHF = AUS->rA[0], a_RadMM = AUS->rA[1], a_TGMM = AUS->rA[2],
                a_ZMM = AUS->rA[3], a_Lim = AUS->rA[4], a_LimP = AUS->rA[5];
    float       a_KHF = AUS->fA[0], a_KMM = AUS->fA[1], a_KAO = (float)fabs(AUS->fA[2]), a_KCPU = AUS->fA[3];
    // here and below #0:HiFrem (laplacian), #1:HiF (main unsharp's inp), #2:HifBlr (blurred org)
    BS_Cpy( cN, AFrm, pA[0] ); BS2_MovAve2D( a_RadHF, pA[0], pA[1] );
    BS_Cpy( cN, pA[1], pA[2] ); BS_SubA( cN, pA[1], pA[0] );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Aa" ); BS_DatasWrt( ADbgPfx, cN, pA, 0, 1 ); }; // #0:HiFrem
    if (a_RadMM) { if (ADbgPfx) a_KAO *= -1; // KAO<0:copies dbgtrc's frames to Ts^
      BS2_MxLapXLimbo( a_RadMM, a_KMM, a_TGMM, a_ZMM, a_KAO, pA[0], qTs, AW, AH ); // maxlap with extended limbo
      if (1 && ADbgPfx && a_KAO<-1e-3f) { // optional trace; ~auto-oscillator <if activated>
        strcpy( lDbgSfx, "Aat" ); BS_DatasWrt( ADbgPfx, cN, qTs, (1<<0)|(1<<2), 0 ); }; }; // if RadMM<>0
    if (a_KCPU<-1e-3f && fabs(a_KCPU+1)>1e-3f ) BS_VMul( cN, pA[0], -a_KCPU );
    BS_RngA( cN, pA[0], a_Lim, a_LimP ); BS2_Sum( 1, pA[1], a_KHF, pA[0] ); // limit HiFrem; HiF += HiFrem*KHF
    if (1 && ADbgPfx) { sprintf( qS, "USstA RadHF=%d KHF=%.2f KCPU=%.2f", a_RadHF, a_KHF, a_KCPU );
      iW( qS ); strcpy( lDbgSfx, "Ab" ); BS_DatasWrt( ADbgPfx, cN, pA, 0, 3 ); }; // #0:HiFrem #1:HiF #2:HiFBlr
    ixTM.Mark( "A enter" ); };
  // #0:HiFrem stays unchanged; unsharp's inp is #1:HiF #2:HifBlr
  int           qw = AW, qh = AH;
  if (AUSExts->rDnSample) { const int cn = cN>>2; qw >>= 1; qh >>= 1; // if DnSample
    BS2_Dim( qw, qh ); // default frames' dimension
    BS_DimDn( qw, qh, pA[1], qTs[0] ); BS_Cpy( cn, qTs[0], pA[1] );
    BS_DimDn( qw, qh, pA[2], qTs[0] ); BS_Cpy( cn, qTs[0], pA[2] );
    if (ADbgPfx) { strcpy( qSCWD, xBS_Dir ); sprintf( qS, "$%04dx%04d", qw, qh ); BS_Dir( qS );
      if (0) { strcpy( lDbgSfx, "A" ); BS_DatasWrt( ADbgPfx, cn, pA, 1, 3 ); }; };
  }; // if DnSample
  // principal UnSharp's output: normalized mskGDi in A#1, the Frame in A#2
  iUnSharp_CPU( AUS, AUSExts, qw, qh, pA[1], pA[2], qTs, ADbgPfx, lDbgSfx ); // unsharp without "step A"
  if (AUSExts->rDnSample) { const int cn = cN>>2; // if DnSample
    if (ADbgPfx) { if (0) { strcpy( lDbgSfx, "X" ); BS_DatasWrt( ADbgPfx, cn, pA, 1, 3 ); }; BS_Dir( qSCWD ); };
    BS_DimUp( qw, qh, pA[1], qTs[0] ); BS_DimUp( qw, qh, pA[2], qTs[1] );
    BS2_Dim( AW, AH ); // default frames' dimension
    BS_Cpy( cN, qTs[0], pA[1] ); BS_DimUpC( AW, AH, qTs[1], pA[2] );
  }; // if DnSample
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "X" ); BS_DatasWrt( ADbgPfx, cN, pA, 0, 3 ); };
  BS_Mul( cN, pA[1], pA[0], -cNrmAveLg ); // HiFrem *= mskGDi
  const float   s_Za2 = AUS->fSA[2]; BS2_Sum( 1, pA[2], s_Za2, pA[0] ); // Frm += HiFrem
  if (0 && ADbgPfx) { strcpy( lDbgSfx, "Xa.raw" ); BS_Write( ADbgPfx, cN, pA[2] ); };
  BS_Cpy( cN, pA[2], AFrm );
  ixTM.Mark( "X Leave" ); ixTM.Mark( NULL ); } // end of BS5_UnSharp_CPU

//// CUDA miscelanneous commons ////
static float iiNrmAveCUDA( tBS_Val *ADat, tBS_Val *ATD, tBS_Val *ATH, tBS_Val *ATmpAbs=NULL, int ANrmAve=cNrmAve ) {
  if (ATmpAbs) { BC_Cpy( ADat, ATmpAbs ); BC_Abs( ATmpAbs ); } else ATmpAbs = ADat; // 'absolute average' or 'plain signed avegare'
  float qA = (float)BC_Mean( ATmpAbs, ATD, ATH ); BC_VMul( ADat, ANrmAve/(qA+0.5f) );
  return qA*0.98f; } // 2024-12-25: adjustment point for "Lap(R)|_CPU ~= |Lap(R)|_CUDA*K" statistically

//// CUDA UnSharp ////
static void iUnSharp_CUDA( const tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH,
  tBS_Val *AHiF, tBS_Val *AHiFblr, tBS_Val **ATsH, tBS_Val **ATsD,
  const char *ADbgPfx, char *lDbgSfx ) { const int cN = AW*AH; char qS[256]; 
  tBS_Val     **qStack; // last used frame
  tBS_Val      *qLapAcc = *ATsD++; BC_ZFil( cN, qLapAcc ); float qKLapAcc = 0; // ad-hoc laplacians' accumulator
  if (AUSExts->rAirEn) { // start concourse "air-estimation"
    tBS_Val *qFrmH = *ATsH++; BC_Move( -cN*cBS_Sz, qFrmH, AHiFblr ); 
    if (ADbgPfx) sprintf( lDbgSfx, "QA" ); BS5a_AirEstimation( AW, AH, qFrmH, ADbgPfx, AUSExts->rAdj ); };
  // - stage B "gamma-correction"; #0:cGammaSpa #1:cGammaCLA #2:cGammaAtt #3:|Der3x3|
  tBS_Val     **pB = ATsD,
               *cGammaSpa = pB[0],   // for spatial metrics (a gamma-corrected frame)
               *cGammaCLA = pB[1],   // for CLAHE
               *cGammaAtt = pB[2],   // for common attenuation
               *cGSDer3x3 = pB[3]; { // global shape |Der3x3|
    const float b_GammaSpa = AUS->fB[0], b_GammaCLA = AUS->fB[1], b_GammaAtt = AUS->fB[2];
    BC_Cpy( AHiF,    cGammaSpa ); if (fabsf(b_GammaSpa-1)>1e-3f) BS3_Gamma( b_GammaSpa, cNrmGam, cGammaSpa );
    BC_Cpy( AHiFblr, cGammaCLA ); if (fabsf(b_GammaCLA-1)>1e-3f) BS3_Gamma( b_GammaCLA, cNrmGam, cGammaCLA );
    BC_Cpy( AHiFblr, cGammaAtt ); if (fabsf(b_GammaAtt-1)>1e-3f) BS3_Gamma( b_GammaAtt, cNrmGam, cGammaAtt );
    BC_Cpy( cGammaCLA, cGSDer3x3 ); BC_DrA3x3( cGSDer3x3, pB+4 ); iiNrmAveCUDA( cGSDer3x3, pB[4], ATsH[0] );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "B" ); BC_DatasWrt( ADbgPfx, pB, 0, 4 ); };
    ixTM.Mark( "B gamma-correction" ); };   
  // - stage C "basic spatial metrics"; #0:BlrBas #1:BlrBig #2:|LapBig|_nrm #3:LapBas #4:LapBig
  tBS_Val     **pC = pB+4,
               *cGSLapBig = pC[2]; { // global shape |LapBig|_nrm
    const int   c_RadBas = AUS->rC[0], c_RadBig = AUS->rC[1], c_RadMM = AUS->rC[2]; // ~4, 205 or 109
    const float c_RefBig = AUS->fC[0], c_DelBig = AUS->fC[1], c_KMM = AUS->fC[2]; // reference norm=1 of |LapBig|
    BC_Cpy( cGammaSpa, pC[0] ); BS3_MovAve2Dp( c_RadBas%100, pC[0], pC+5, c_RadBas/100 ); // #0:BlrBas
    BC_Cpy( cGammaSpa, pC[3] ); BC_SubA( pC[0], pC[3] ); // #3:LapBas
    if (c_RadMM>0 && c_KMM>-1e-3f) {
      BC_MxLapStd( c_RadMM, pC[3], pC+5, c_KMM ); // #3:(1-KMM)*LapBas+KMM*MaxMag(LapBas)
      BC_Cpy( pC[0], cGammaSpa ); BC_Add( pC[3], cGammaSpa ); }; // update GammaSpa
    BC_Cpy( cGammaSpa, pC[1] ); BS3_MovAve2Dp( c_RadBig%100, pC[1], pC+5, c_RadBig/100 ); // #1:BlrBig
    BC_Cpy( cGammaSpa, pC[4] ); BC_SubA( pC[1], pC[4] ); // #4:LapBig
    BC_Cpy( pC[4], cGSLapBig ); BC_Abs( cGSLapBig );
    AUSExts->rSttAveBig = iiNrmAveCUDA( cGSLapBig, pC[5], ATsH[0] ); // #2:|LapBig|_nrm
    AUSExts->rSttNrmDelBig = iiMetricNrmDel( c_RefBig, c_DelBig, AUSExts->rSttAveBig );
    if (1 && ADbgPfx) { // report UnSharpExts.AveBig, relative <against reference> |LapBig|
      sprintf( qS, "USstC ave|Big|=%.2f NrmDelBig=%.2f", AUSExts->rSttAveBig, AUSExts->rSttNrmDelBig ); iW( qS );
      strcpy( lDbgSfx, "C" ); BC_DatasWrt( ADbgPfx, pC, 0, 5 ); };
    ixTM.Mark( "C basic spatial metrics" ); };
  // - stage D "CLAHE"; #0:plain_CLAHE #1:mskBodyOne #2:CLAInc #3:|Der3x3CLA|
  tBS_Val     **pD = pC+5,
               *cGSDer3x3CLA = pD[3], // global shape |Der3x3CLA|
               *qmskBodyOne = NULL; { 
    const int   d_Blk = AUS->rD[0], d_Blr = AUS->rD[1], // ~6,9
                d_IncRadBlr = 1;
    const float d_BodyOne = AUS->fD[0], d_ExtMn = AUS->fD[1], d_ExtMx = AUS->fD[2];
    AUSExts->rSttCLABlr = d_Blr<=0 ? d_Blr : iiByNrmDel( AUSExts->rSttNrmDelBig, d_Blr, -0.1f, 0.25f ); // actual d_Blr
    tBS4_CLAHE_Exts qCLAExts = {d_ExtMn,d_ExtMx}; if (ADbgPfx) strcpy( lDbgSfx, "Da" );
    BS4_CLAHE_Transform( // #0:plain_CLAHE #1:BodyOne
      d_Blk, cNrmGam, AUSExts->rSttCLABlr, cGammaCLA, pD, true, ADbgPfx, &qCLAExts );
    if (0 && ADbgPfx) { strcpy( lDbgSfx, "Db" ); BC_DatasWrt( ADbgPfx, pD, 0, 2 ); };
    BC_SubA( cGammaCLA, pD[0] ); // #0:CLAInc:=CLAHE(Gamma1)-Gamma1
    BS3_MovAve2D( d_IncRadBlr, pD[0], pD[2], pD[4] ); // #2:CLAInc blurred
    BC_Cpy( cGammaCLA, pD[0] ); BC_Add( pD[2], pD[0] ); // #0:plain_CLAHE
    if (d_BodyOne>1e-3f) { qmskBodyOne = pD[1]; // #1:BodyOne 
      BC_VMul( qmskBodyOne, cNrmAve*1.2f/((1<<cBS4_XLATGamLg)*d_BodyOne) );
      BC_Rng( qmskBodyOne, -999999, cNrmAve-1 );
      BS3_MovAve2Dp( iiBodyOne_Rblr(d_BodyOne,AW,AH), qmskBodyOne, pD+4, 2 );
      BC_Mul( qmskBodyOne, pD[2], -cNrmAveLg ); }; // attenuate CLAInc
    BC_Cpy( pD[0], cGSDer3x3CLA ); BC_DrA3x3( cGSDer3x3CLA, pD+4 ); iiNrmAveCUDA( cGSDer3x3CLA, pD[4], ATsH[0] );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Dc" ); BC_DatasWrt( ADbgPfx, pD, 0, 4 ); };
    ixTM.Mark( "D CLAHE" ); qStack = pD+4; };
  // - stage EA "to attenuate mskGDi,mskEdg; #0:mskAtt
  tBS_Val      *qmskAtt = NULL; const int e_Blk = AUS->rE[0]; if (e_Blk>1) {
    const int   e_Blr = AUS->rE[1], e_Lim = AUS->rE[2], e_Zero = AUS->rE[3]; // Lim is against eBlk=6
    tBS_Val   **pEA = qStack, *qEt = pEA[3]; BC_Cpy( cGammaAtt, qEt );
    BS4_FrmNrm_CUDA( cNrmGam, qEt, pEA ); qmskAtt = pEA[0]; // #0:normalized GammaAtt
    BS4_CLAHE_TraAveRef( e_Blk, e_Blr, qEt, pEA, true, (tBS_Val*)-1 ); // #0:hist_counts of #t
    if (0 && ADbgPfx) { strcpy( lDbgSfx, "EAb_His.raw" ); BC_Write( ADbgPfx, cN, qmskAtt ); };
    const int csh = (e_Blk-6)<<1, cLim = csh>0 ? e_Lim<<csh : e_Lim>>-csh; BC_Neg( qmskAtt );
    if (0 && csh) { sprintf( qS, "USstEA Blk=%d Lim %d=>%d", e_Blk, e_Lim, cLim ); iW( qS ); };
    BC_VAdd( qmskAtt, cLim ); BC_Rng( qmskAtt, 0 ); BC_VMul( qmskAtt, (float)cNrmAve/cLim );
    if (e_Zero) BC_VAdd( qmskAtt, -e_Zero );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "EAc_Att.raw" ); BC_Write( ADbgPfx, cN, qmskAtt ); };
    ixTM.Mark( "EA MskAtt" ); qStack = pEA+1;
  }; // if eBlk
  // - stage F "the masks"; #0:mskThn #1:mskThk #2:mskGDi #3:.. 
  tBS_Val     **pF = qStack, // the reference: _ReadMe.txt, "Generation of normalized masks at stage F"  
               *qmskThn = pF[0], // ~'Thin'; merge LapAda,LapMM to LapAdaEff
               *qmskThk = pF[1], // ~'Thick'; merge LapBas,LapBig to LapAda
               *qmskGDi = pF[2], // ~ dilate HiFrem,LapAdaEff
               *qmskEdg = pF[3]; // ~ a bump at lap's zero-izoline
  bool          qBigGThn_blr; {  // true when ~blur LapMM should go rather to big GThn
    float f_Att = AUS->fF[0], f_EmLthn = AUS->fF[1], f_EmLedg = AUS->fF[2], // ~ 6.0 0.1 0.05
          f_MxLbias = AUS->fF[3]; // ~ 1.15
    float f_K[3][3]; memcpy( f_K, AUS->fFK, sizeof(AUS->fFK) ); // general shapes from basic stats
    int   f_C[4][5]; memcpy( f_C, AUS->rFC, sizeof(f_C) ); // normalized masks by ~shapes
    qBigGThn_blr = f_C[0][0]>f_C[0][1]; if (1) { sprintf( qS, "USstF BigGThn_blr=%d rcuBlack,White={%d,%d}",
      (int)qBigGThn_blr, f_C[0][0], f_C[0][1] ); iW( qS ); };
    int   f_RadThk = f_C[1][4], // if (f_C[1][0]>f_C[1][1]) AUSExts->rSttFlg |= 4; // f_Cb[0]>f_Cb[1]
          qf_Khn = 0; // f_K[k] has non-zero coefs
    tBS_Val **qFd = pF, **qFb = qFd+4, **qFc = qFb+mcASz(f_K), **qT = qFc+2;
    { // <shapes> @b0:~GSh @b1:~GDi @b2:~GThn
      if (1 && ADbgPfx) { tBS_Val *cTs[] = {cGSDer3x3,cGSDer3x3CLA,cGSLapBig};  // optional summary
        strcpy( lDbgSfx, "Fa" ); BC_DatasWrt( ADbgPfx, cTs, 0, mcASz(cTs) ); }; // trace of inp_shapes
      for (int qk=0;qk<mcASz(f_K);qk++) { // K0*|Der3x3|+K1*|Der3x3CLA|+K2*|LapBig| (all normalized)
        int qHaveNZ = iiHaveNZ( 3, f_K[qk] ); if (!qHaveNZ) continue;
        qf_Khn |= 1<<qk; BC_Cpy( cGSDer3x3, qFb[qk] );
        BS3_Sum( f_K[qk][0], qFb[qk], f_K[qk][1], cGSDer3x3CLA, f_K[qk][2], cGSLapBig );
        iiNrmAveCUDA( qFb[qk], qT[0], ATsH[0], qHaveNZ<0?qT[1]:NULL ); };
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Fb" ); BC_DatasWrt( ADbgPfx, qFb, qf_Khn, 0 ); };
      ixTM.Mark( "Fb shapes" ); };
    { // <hard-coded logic of spatial filters> @c0:MskThk @c1:MskGDi @c2:MskEdg
      BC_Cpy( qFb[0], qFc[0] ); BS3_MovAve2Dp( f_RadThk%100, qFc[0], qT, f_RadThk/100 ); // @c0:MskThk
      BC_Cpy( qFc[0], qFc[2] ); BC_SubS( qFb[0], qFc[2] ); // @c2:MskEdg
      if (!(qf_Khn&(1<<1))) { BC_Cpy( qFc[2], qFb[1] ); BC_Abs( qFb[1] ); // if no Kb then GDi:=abs(MskEdg)
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "FB1.raw" ); BC_Write( ADbgPfx, cN, qFb[1] ); }; };
      if (!(qf_Khn&(1<<2))) { BC_Cpy( qFc[2], qFb[2] ); BC_Neg( qFb[2] ); // if no Kc then GThn:=-MskEdg
        if (1 && ADbgPfx) { strcpy( lDbgSfx, "FB2.raw" ); BC_Write( ADbgPfx, cN, qFb[2] ); }; };
      BC_Cpy( qFb[1], qFc[1] ); // @c1:MskGdi:=GDi+MskAtt*KGDi
      if (qmskAtt) { BC_VMul( qmskAtt, f_Att ); BC_Add( qmskAtt, qFc[1] ); };
      if (1 && ADbgPfx) { // index+semantic in this trace_dump
        strcpy( lDbgSfx, "Fc0Thk.raw" ); BC_Write( ADbgPfx, cN, qFc[0] );
        strcpy( lDbgSfx, "Fc1GDi.raw" ); BC_Write( ADbgPfx, cN, qFc[1] );
        strcpy( lDbgSfx, "Fc2Edg.raw" ); BC_Write( ADbgPfx, cN, qFc[2] ); };
      ixTM.Mark( "Fc Masks" ); };
    { // cut gamut's of denormalized MskXYZ to mskXYZ (and perhaps blur)
      BC_Cpy( qFb[2], qmskThn ); BC_Cpy( qFc[0], qmskThk );
      BC_Cpy( qFc[1], qmskGDi ); BC_Cpy( qFc[2], qmskEdg );
      int   qf_Chnzr = 0; for (int qc=0;qc<4;qc++) { // f_C[c] cut has nonzero range
        if (f_C[qc][0]==f_C[qc][1]) continue; qf_Chnzr |= 1<<qc; // masks within 0..NrmAve
        int cTg0 = f_C[qc][2], cGaO = cNrmAve-cTg0-f_C[qc][3];
        if (cGaO>2) { tBS2_Cut cCut = { {f_C[qc][0],f_C[qc][1],cGaO}, {0,0} };
          BS3_Cut( &cCut, qFd[qc] ); if (cTg0) BC_VAdd( qFd[qc], cTg0 ); }; // if GaO
        int qRad = f_C[qc][4]; BS3_MovAve2Dp( qRad%100, qFd[qc], qT, qRad/100 ); }; // for c
      if (!(qf_Chnzr&(1<<2))) BC_Cpy( qFd[1], qFd[2] ); // if no Cc then mskGDi := mskThk
      if (f_MxLbias>1e-3f) { // mxlap for mskGdi, mskEdg
        int qRadThk = (f_RadThk%100)<<(f_RadThk/100); if (qRadThk>8) qRadThk = 8;
        BC_Med( qmskGDi, qT[0], 2 ); BC_Add( qT[0], qmskGDi ); BC_VMul( qmskGDi, 0.5f );
        BC_MxLapStd( qRadThk+2, qmskGDi, qT, f_MxLbias, -1, 0, 0 );
        if (f_MxLbias>1) BC_Rng( qmskGDi, -99999, cNrmAve-1 );
        BC_Mul( qmskEdg, qmskEdg, -cNrmAveLg ); // <possibly> square edge_mask
        BC_MxLapStd( qRadThk, qmskEdg, qT, f_MxLbias, -1, 0, 0 );
        if (f_MxLbias>1) BC_Rng( qmskEdg, -99999, cNrmAve-1 ); };
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Fd" ); BC_DatasWrt( ADbgPfx, qFd, qf_Chnzr, 0 ); };
      ixTM.Mark( "Fd norm masks" ); };
    if (BS3_AccMim( qLapAcc, qT[0], (float)fabs(f_EmLthn), qmskThn, f_EmLthn<0?pC[3]:pC[4],
                                    (float)fabs(f_EmLedg), qmskEdg, f_EmLedg<0?pC[3]:pC[4] )) { // Edge mimic Laplacian
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "FeEmL.raw" ); BC_Write( ADbgPfx, cN, qLapAcc ); };
      ixTM.Mark( "Ff acc EmL" ); qKLapAcc = 1; };
    if (AUSExts->rMskEn) { if (ADbgPfx) strcpy( lDbgSfx, "QC" ); // start concourse "extended mask-recount"
      BS5d_Msk( AW, AH, true, ADbgPfx, pC[1], qmskThk, qmskGDi ); };
    qStack = qFd+4; }; // last used frame
  // - stage H "LapAdaptive"; #0:LapMid #1:BlrAda #2:LapAda
  // Merge LapBas,LapBig to LapAda by mskThk; i.e. attenuate (LapBas-LapBig).
  tBS_Val     **pH = qStack; {
    const int   h_Lim = AUS->rH[0], h_LimP = AUS->rH[1]; // truncate LapAda
    BC_Cpy( pC[0], pH[0] ); BC_SubA( pC[1], pH[0] ); // #0:LapMid:=BlrBas-BlrBig (also '=LapBig-LapBas')
    BC_Cpy( pH[0], pH[1] ); BC_Mul( qmskThk, pH[1], -cNrmAveLg ); // #1:LapMid_dil:=mskThk*LapMid
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Ha" ); BC_DatasWrt( ADbgPfx, pH, 0, 2 ); };
    BC_Add( pC[1], pH[1] ); // #1:BlrAda:=BlrBig+LapMid_dil
    BC_Cpy( cGammaSpa, pH[2] ); BC_SubA( pH[1], pH[2] ); // #2:LapAda:=GammaSpa-BlrAda
    BC_RngA( pH[2], h_Lim, h_LimP );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Hb" ); BC_DatasWrt( ADbgPfx, pH, 1, 3 ); };
    ixTM.Mark( "H LapAdaptive" ); qStack = pH+3; };
  // - stage J "Radon & LapAdaEff"; #0:LapMM #1:LapAdaEff
  // Get LapMM, then merge LapAda,LapMM to LapAdaEff by mskThn; i.e. attenuate (LapAda-LapMM).
  tBS_Val     **pJ = qStack; {
    const int   j_RadMM = AUS->rJ[0], j_RadMMO = AUS->rJ[1], j_RadENm = AUS->rJ[2], j_Lim = AUS->rJ[3];
    const float j_KMM = AUS->fJ[0], j_KENm = AUS->fJ[1];
    tBS_Val    *qqbufs[cBS5c_QuadLaps], **qqQuadLaps = AUSExts->rQuadEn ? 
                BS5c_QuadLapsAlloc(AW,AH,qqbufs,j_RadMMO,cQRadO) : NULL; // CPU-memory for "quad-laps"
    const tBS2_RadonMM cBS2_RadonMM = {j_RadMM,j_RadMMO,cQRad,cQRadO,j_RadENm};
    tBS_Val    *qENm = j_KENm>1e-3f ? pJ[3] : NULL;
    BS3_RadonMM( &cBS2_RadonMM, pH[2], // Radon's denoise of ~LapAdaptive
      pJ[2], pJ[0], pJ[1], qENm, pJ+4, qqQuadLaps ); // Dev Blr Srp ENm
    BC_VMul( pJ[0], j_KMM ); // #0:LapMM*KMM
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Ja" ); BC_DatasWrt( ADbgPfx, pJ, 0, qENm?4:3 ); };
    BC_SubA( pJ[0], pJ[1] ); BC_Mul( qmskThn, pJ[1], -cNrmAveLg ); // #1:=(Srp-Blr)*=mskThn
    BC_Add( pJ[0], pJ[1] ); // #1:LapAdaEff (+=Blr)
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jb" ); BC_DatasWrt( ADbgPfx, pJ, 1, 2 ); };
    if (qENm) { // add ~limit_ENm_by_Dev
      BS3_RngLapByConductor( pJ[2], qENm ); BC_VMul( qENm, j_KENm );
      BC_Add( qENm, pJ[1] ); // #1:LapAdaEff += ENm*KENm
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jc" ); BC_DatasWrt( ADbgPfx, pJ, (1<<1)|(1<<3), 0 ); }; };
    if (BC_RngA( pJ[1], j_Lim ))
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "Jd_LapAdaEff_J1.raw" ); BC_Write( ADbgPfx, cN, pJ[1] ); };
    if (AUSExts->rQuadEn) { // start concourse "quad-laps"; the 'LapOrg' could vary as ~...
      if (ADbgPfx) strcpy( lDbgSfx, "QB" ); // ... H#0:LapMid H#2:LapAda J#0:LapMM
      BS5c_QuadLaps( AW, AH, qqQuadLaps, true, ADbgPfx, pJ[1], qmskThk ); };
    ixTM.Mark( "J RadonDenoise" ); qStack = pJ+2; };
  // - stage LA "edge"; #0:mskEdg_blr, LapAcc+=Edg
  { const int   la_RadHess = AUS->rL[0]; const float la_KEdg = AUS->fL[0]; if (la_RadHess && la_KEdg>1e-3f) {
      tBS_Val **pLA = qStack, *qEdg = pLA[0]; BC_Cpy( qmskEdg, qEdg ); // ~hessian of 1st derivation from mskEdg
      BS3_MimicDat( pJ[0], qEdg ); // #0:lapedg=LapMM
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAa" ); BC_DatasWrt( ADbgPfx, pLA, 0, 1 ); };
      BS3_Hess( abs(la_RadHess), pLA ); // #0:Hessian #1:LapEdgBlr
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAb" ); BC_DatasWrt( ADbgPfx, pLA, 0, 3 ); };
      BS3_MimicDat( pLA[1], qEdg ); BC_VMul( qEdg, la_KEdg ); // #0:Edg
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LAcEdg.raw" ); BC_Write( ADbgPfx, cN, qEdg ); };
      BC_Add( qEdg, qLapAcc ); qKLapAcc = 1;
      if (true) BC_Cpy( pLA[2], qmskEdg ); // optionally reassign by mskEdgBlr
      ixTM.Mark( "LA Edge" ); qStack = pLA+1; }; };
  // - stage LB "ridge"; LapAcc+=Rdg
  { const float lb_KRdg = AUS->fL[1]; if (lb_KRdg>1e-3f) {
      tBS_Val **pLB = qStack; BC_Cpy( pJ[1], pLB[0] ); // could be either J#0,#1 (either LapMM,LapAdaEff) 
      BS3_RidgeByConductors( pLB[0], cGSLapBig, cGSDer3x3, cGSDer3x3CLA ); // ... |LapBig|_nrm |Der3x3|_nrm |Der3x3CLA|_nrm
      BC_VMul( pLB[0], lb_KRdg ); BC_Add( pLB[0], qLapAcc ); qKLapAcc = 1;
      if (1 && ADbgPfx) { strcpy( lDbgSfx, "LBaRdg.raw" ); BC_Write( ADbgPfx, cN, pLB[0] ); };
      ixTM.Mark( "LB Ridge" ); }; };
  // - stage LC "maxmag"
  { bool        qlc_chg = false;
    const int   lc_RadMM = AUS->rL[1]; const float lc_KMM = AUS->fL[2]; if (lc_RadMM>0 && lc_KMM>-1e-3f) {
      BC_MxLapStd( lc_RadMM, qLapAcc, qStack, lc_KMM ); qlc_chg = true; };
    const int   lc_Lim = AUS->rL[2]; qlc_chg |= BC_RngA( qLapAcc, lc_Lim );
    if (qlc_chg) { if (0 && ADbgPfx) { strcpy( lDbgSfx, "LCa_lapacc.raw" ); BC_Write( ADbgPfx, cN, qLapAcc ); };
      ixTM.Mark( "LC LapAcc refine" ); qKLapAcc = 1; }; };
  // - stage QA "wait for air-estimation"
  tAirEstimationWait qAEW; memset( &qAEW, 0, sizeof(tAirEstimationWait) ); if (AUSExts->rAirEn) {
    BS5a_AirEstimationWait( &qAEW ); tBS_Val **pQ = NULL;
    if (qAEW.rAiF[0]) { pQ = qStack; qStack += 1; BC_Move( cN*cBS_Sz, qAEW.rAiF[0], pQ[0] );
      if (1 && ADbgPfx ) { strcpy( lDbgSfx, "QAzBkgPtr.raw" ); BC_Write( ADbgPfx, cN, pQ[0] ); }; };
    for (int qt=0;qt<mcASz(qAEW.rAiF);qt++) if (qAEW.rAiF[qt]) delete qAEW.rAiF[qt]; // release ~'BkrPtr from air-detection'
    if (pQ) qAEW.rAiF[0] = pQ[0];
    AUSExts->rAirInfFlg = qAEW.rDbgInfFlg; };
  // - stage QB "wait for quad-laps"
  if (AUSExts->rQuadEn) { tBS_Val *qT = BS5c_QuadLapsWait(); if (qT) {
    BC_Move( cN*cBS_Sz, qT, qStack[0] ); BC_Add( qStack[0], qLapAcc ); qKLapAcc = 1; delete qT; }; };
  // - stage QC "wait for extended mask-recount"
  if (AUSExts->rMskEn) BS5d_MskWait();
  // - stage N "global dilation"; LapAdaEff*=mskGDiRes, HiF:=mskGDiRes
  { const float n_Thk = AUS->fN[0], n_GDi = AUS->fN[1], n_Edg = AUS->fN[2]; // ~0 0.4 0.2
    tBS_Val   **pN = qStack, *mskGDiRes = AHiF;
    BS3_MskSum( mskGDiRes, pN[0], cNrmAve, n_Thk, qmskThk, n_GDi, qmskGDi, n_Edg, qmskEdg );
    if (qAEW.rAiF[0]) BC_Mul( qAEW.rAiF[0], mskGDiRes, -8 ); // mskGDiRes*=Air
    BC_Mul( mskGDiRes, pJ[1], -cNrmAveLg ); // J#1:LapAdaEff_dil (by mskGDiRes)
    if (1 && ADbgPfx) { 
      strcpy( lDbgSfx, "Na_GDiRes_X1.raw" ); BC_Write( ADbgPfx, cN, mskGDiRes );
      strcpy( lDbgSfx, "Nb_LapAdaEff_J1.raw" ); BC_Write( ADbgPfx, cN, pJ[1] ); };
    ixTM.Mark( "N Dilation" ); };
  // - stage S "summa summarum"; #0:BodySum #1:LapSum #2:ResSum
  tBS_Val     **pS = qStack; { // s_Za2*HiFRem adds later
    const float s_Za0 = AUS->fSA[0], s_Za1 = AUS->fSA[1], s_Za3 = AUS->fSA[3], // #0:BodySum
                s_Zb0 = AUS->fSB[0], s_Zb1 = AUS->fSB[1], s_Zb2 = AUS->fSB[2]; // #1:LapSum
    BC_Cpy( pH[1], pS[0] ); BS3_Sum( s_Za0, pS[0], s_Za1, cGammaCLA,
      1-s_Za0-s_Za1, AHiFblr, (float)fabs(s_Za3), pH[s_Za3<-1e-3f?1:0] );
    BC_Cpy( pD[2], pS[1] ); BS3_Sum( s_Zb0, pS[1], s_Zb1, pH[2], s_Zb2, pJ[1], qKLapAcc, qLapAcc );
    // don't multiply ~"large_low_freqs of CLAInc in LapSum" on ~"GDiRes"
    if (qmskBodyOne) BC_Mul( qmskBodyOne, pS[1], -cNrmAveLg ); // attenuate LapSum
    if (qAEW.rAiF[0]) BC_Mul( qAEW.rAiF[0], pS[1], -8 ); // #1:LapSum *= Air
    BC_Cpy( pS[0], pS[2] ); BC_AddRg( pS[1], pS[2] ); // #2:ResSum
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "S" ); BC_DatasWrt( ADbgPfx, pS, 0, 3 );
      if (qKLapAcc>1e-3f) { strcpy( lDbgSfx, "S_lapacc.raw" ); BC_Write( ADbgPfx, cN, qLapAcc ); }; };
    ixTM.Mark( "S summa summarum" ); };
  // - stage U "histogram cut out"; #0:Res
  tBS_Val     **pU = pS+3; { tBS2_HistCut qHistCut;
    int         qGapNeg = iiHistCutOut( &qHistCut, qAEW.rAir[0], qAEW.rAir[1], AUS, ADbgPfx!=NULL );
    BC_Cpy( pS[2], pU[0] ); BS3_HistCut( &qHistCut, pU[0], NULL, pU+1 );
    BC_VAdd( pU[0], qGapNeg ); // compensate the most of expecting negatives
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "U" ); BC_DatasWrt( ADbgPfx, pU, 0, 1 ); 
      sprintf( qS, "USstU CatchMM=%d,%d", xBS2_rhCatchMM[0], xBS2_rhCatchMM[1] ); iW( qS ); };
    ixTM.Mark( "U histogram cut out" ); };
  // - stage V "finalize"; #0:ResBkGM
  tBS_Val     **pV = pU+1; BC_Cpy( pU[0], pV[0] );
  const int     mgn_Msk = AUS->rMgn[0]; if (mgn_Msk>0) { // blur surrounding frame's margins when '>0'
    BC_Move( cN*sizeof(tBS_Val), AHiFblr, pV[1] ); BS3_VMimicNZ( pV[1], 1<<8 );
    BS3_MovAve2Dp( 2, pV[1], pV+2, 1 ); BC_Mul( pV[0], pV[1], -8 );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "V" ); BC_DatasWrt( ADbgPfx, pV, 0, 2 ); };
    ixTM.Mark( "V finalize" ); };
  BC_Rng( pV[0], 0 ); BC_Cpy( pV[0], AHiFblr ); // non negative output data
  AUSExts->rStackFrmBuf = (int)(pV-ATsD); } // end of iUnSharp_CUDA

//// CUDA unsharp wrapper ////
void BS5_UnSharp_CUDA( tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH, tBS_Val *AFrm, tBS_Val **ATsH, tBS_Val **ATsD,
  const char *ADbgPfx, char *lDbgSfx ) { // -- UnSharp CUDA (declaration for ~"no downsample")
  const int     cN = AW*AH; char qS[256], qSCWD[40];
  if (xDbgDmp==1) strcpy( ixTM.rPfx, "US" ); ixTM.Reset(); // tick-marks
  BS2_Dim( AW, AH ); BS3_Dim( AW, AH ); // default frames' dimension
  jUnSharpRecReConsider( AW, AH, AUS, AUSExts );
  // - stage A "high freqs preliminaries"; #0:HiFrem #1:HiF #2:HiFblr
  tBS_Val     **pA = ATsD, **qTsD = ATsD+3; {
    int         a_RadHF = AUS->rA[0], a_RadMM = AUS->rA[1], a_TGMM = AUS->rA[2],
                a_ZMM = AUS->rA[3], a_Lim = AUS->rA[4], a_LimP = AUS->rA[5];
    float       a_KHF = AUS->fA[0], a_KMM = AUS->fA[1], a_KAO = (float)fabs(AUS->fA[2]), a_KCPU = AUS->fA[3];
    // here and below #0:HiFrem (laplacian), #1:HiF (main unsharp's inp), #2:HifBlr (blurred org)
    BC_Move( cN*cBS_Sz, AFrm, pA[0] ); BS3_MovAve2D( a_RadHF, pA[0], pA[1], qTsD[0] );
    BC_Cpy( pA[1], pA[2] ); BC_SubA( pA[1], pA[0] );
    if (1 && ADbgPfx) { strcpy( lDbgSfx, "Aa" ); BC_DatasWrt( ADbgPfx, pA, 0, 1 ); }; // #0:HiFrem
    if (a_RadMM) { if (ADbgPfx) a_KAO *= -1; // KAO<0:copies dbgtrc's frames to Ts^
      BS3_MxLapXLimbo( a_RadMM, a_KMM, a_TGMM, a_ZMM, a_KAO, pA[0], qTsD, ATsH[0] ); // maxlap with extended limbo
      if (1 && ADbgPfx && a_KAO<-1e-3f) { // optional trace; ~auto-oscillator <if activated>
        strcpy( lDbgSfx, "Aat" ); BC_DatasWrt( ADbgPfx, qTsD, (1<<0)|(1<<2), 0 ); }; }; // if RadMM<>0
    if (a_KCPU>1e-3f && fabs(a_KCPU-1)>1e-3f ) BC_VMul( pA[0], a_KCPU );
    BC_RngA( pA[0], a_Lim, a_LimP ); BS3_Sum( 1, pA[1], a_KHF, pA[0] ); // limit HiFrem; HiF += HiFrem*KHF
    if (1 && ADbgPfx) { sprintf( qS, "USstA RadHF=%d KHF=%.2f KCPU=%.2f", a_RadHF, a_KHF, a_KCPU );
      iW( qS ); strcpy( lDbgSfx, "Ab" ); BC_DatasWrt( ADbgPfx, pA, 0, 3 ); }; // #0:HiFrem #1:HiF #2:HiFBlr
    ixTM.Mark( "A enter" ); };
  // #0:HiFrem stays unchanged; unsharp's inp is #1:HiF, #2:HifBlr
  int           qw = AW, qh = AH;
  if (AUSExts->rDnSample) { const int cn = cN>>2; qw >>= 1; qh >>= 1;  // if DnSample
    BS2_Dim( qw, qh ); BS3_Dim( qw, qh ); // default frames' dimension
    BC_DimDn( pA[1], qTsD[0] ); BC_DimDn( pA[2], qTsD[1] );
    BC_Cpy( qTsD[0], pA[1] ); BC_Cpy( qTsD[1], pA[2] );
    if (ADbgPfx) { strcpy( qSCWD, xBS_Dir ); sprintf( qS, "$%04dx%04d", qw, qh ); BS_Dir( qS );
      if (0) { strcpy( lDbgSfx, "A" ); BC_DatasWrt( ADbgPfx, pA, 1, 3 ); }; };
  }; // if DnSample
  // principal UnSharp's output: normalized mskGDi in A#1, the Frame in A#2
  iUnSharp_CUDA( AUS, AUSExts, qw, qh, pA[1], pA[2], ATsH, qTsD, ADbgPfx, lDbgSfx ); // unsharp without "step A"
  if (AUSExts->rDnSample) { // if DnSample
    if (ADbgPfx) { if (0) { strcpy( lDbgSfx, "X" ); BC_DatasWrt( ADbgPfx, pA, 1, 3 ); }; BS_Dir( qSCWD ); };
    BC_DimUp( pA[1], qTsD[0] ); BC_DimUp( pA[2], qTsD[1] );
    BS2_Dim( AW, AH ); BS3_Dim( AW, AH ); // default frames' dimension
    BC_Cpy( qTsD[0], pA[1] ); BC_DimUpC( qTsD[1], pA[2] ); 
  }; // if DnSample
  BC_Mul( pA[1], pA[0], -cNrmAveLg ); // HiFrem *= mskGDi
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "X" ); BC_DatasWrt( ADbgPfx, pA, 0, 3 ); };
  const float   s_Za2 = AUS->fSA[2]; BS3_Sum( 1, pA[2], s_Za2, pA[0] ); // Frm += HiFrem
  if (0 && ADbgPfx) { strcpy( lDbgSfx, "Xa.raw" ); BC_Write( ADbgPfx, cN, pA[2] ); };
  BC_Move( -cN*cBS_Sz, AFrm, pA[2] );
  ixTM.Mark( "X Leave" ); ixTM.Mark( NULL ); } // end of BS5_UnSharp_CUDA
