
// I16_Ut5a.cpp : air-estimation in CPU (within downsampled frame)
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
//#include "BS_Ut1.h"
#include "BS_Ut2.h"
//#include "BS_Ut3.h"
#include "BS_Ut4.h"
#include "BS_Ut5.h"
#include "BS_Ut5a.h"

//// miscelanneous commons ////  
static int ixAirEstimationEnabled = 0; // 0:disabled, 1:Lims, 2:Teeth
static void iW( const char *AMsg ) { // printed error_tags here are ~informal (don't affect behavior) 
  char qS[256]; sprintf( qS, "ae %.252s", AMsg ); WriteLn( qS ); }
static inline int iDiv32_1000( int N, int D ) { return (int)(((tLL)N*1000+(D>>1)-1)/D); } // N/D <in 0.001>

static struct { // Lims' air-estimation 
  float // rK,    // Org+ref_average(counts)*K
        rGamma;   // Gamma (Org=>Org_gamma,Org_8bit)
  int   rDerDec,  // cDerDec (|Der3x3|=>Der_lap)
        rCntLim,  // cLim (Sgnd_ivls=>Sgnd_ivls_lim)
        rBlk,     // CLAHE's block (Sgnd_ivls_lim=>Ref_ave)
        rcuBlack, // rcuBlack (A_metric_nrm=target0)
        rcuWhite, // rcuWhite (A_metric_nrm=255)
        rcuTg0,   // "target0"
        rRad;     // blur(A_metric_nrm)
} ixAL = {
  0.9f, // Gamma (Dp#0=>Dp#0,Dp#1)
  555,  // cDerDec (Dp#2=>Dp#2)
  333,  // cLim (Dp#3=>Dp#3)
  5,    // CLAHE's block (Dp#1,Dp#3=>Dp#4)
  2222, 4567, 11, // rcuBlack White Tg0 (Dp#4=>Dp#5)
  203 // blur_rad (Dp#5=>Dp#6)
};
// the main attention is "Dp#4=>Dp#5" (i.e. A_metric_den=>A_metric_nrm)

static void iAirEstimationLoad_Lims( GIniFile *AF ) {
  #if 0
  const char *cS; char qS[256];
  cS = AF->GetStr( "1a", NULL ); if (!cS || !sscanf( cS, "%d %d %d %f", &ixAL.rDerDec, &ixAL.rCntLim, &ixAL.rBlk, &ixAL.rGamma )) {
    sprintf( qS, "%d %d %d %.2f; Dec Lim Blk Gamma", ixAL.rDerDec, ixAL.rCntLim, ixAL.rBlk, ixAL.rGamma ); AF->SetStr( "1a", qS ); };
  cS = AF->GetStr( "1b", NULL ); if (!cS || !sscanf( cS, "%d %d %d %d", &ixAL.rcuBlack, &ixAL.rcuWhite, &ixAL.rcuTg0, &ixAL.rRad )) {
    sprintf( qS, "%d %d %d %d; Blk Wht Tg0 Rad", ixAL.rcuBlack, ixAL.rcuWhite, ixAL.rcuTg0, ixAL.rRad ); AF->SetStr( "1b", qS ); };
  #endif
}

static void iAirEstimationLoad_Teeth( GIniFile *AF ) {}

static const int
  cNrmGaLg = 14, cNrmGa = 1<<cNrmGaLg;

//// air-esimation's core ////
struct tAiWr { public: // - envelope of the 'AirEstimation'
  int rW,rH; const tBS_Val *cOrgPtr; tBS_Val *rBkgPtr; // frames 
  int rAiA[2]; float rAir[2]; // output scalars
  int rcTag,rcB,rcC; // inp adjustments (for optional reference)
  char rDbgPfx[256], *rDbgSfx; // debug-dump's strings,
  bool rDbgTrc; // trace at least min messages
  int  rDbgInfFlg; // short dbg-info flags
  void DbgSet( const char *ADbgPfx ) { if (ADbgPfx) rDbgSfx = rDbgPfx+sprintf( rDbgPfx, "%s_", ADbgPfx ); rDbgTrc = true; }
  void AirEstimation() { switch (ixAirEstimationEnabled) {
    case 1: iAirEstimation_Lims(); break;
    case 2: iAirEstimation_Teeth(); break; }; }

  protected:
  //                     Lims' air-estimation in short:
  // Distinguish ~large air's areas by spatial counts of Laplacian(|Derivation|)".
  // Input image Org "0..Org_gamut" has approx relation "Air < Edges < Objects".
  //   step a) Get Laplacian(|Derivation|) as Der_lap:=|Der3x3|-cDerDec.
  //   step b) Sgnd_ivls:=directional_counts(Der_lap). Expect here statistically
  // "Air_average < -cLim < Objects_mn < 0 < Edges", where "100 < cLim < 300"
  // (for subjective Air|Objects of our downsampled images ~500x500). We also
  // constraint the interval as ~"Sgnd_ivls_lim:=Rng(-cLim,cLim)".
  //   step c) Ref_ave:=_CLAHE_TraRef(Org_8bit,Sgnd_ivls_lim); i.e. spatially
  // attenuate signed intervals by ~CLAHE.
  //   step d) A_metric_den:=Org+K*Ref_ave. From 
  // <in Org>          "Air < Edges < Objects" and
  // <in Ref_ave>      "Air < Objects < Edges" expect, that
  // <in A_metric_den> "Air < Objects,Edges".
  //   step e) Cut "A_metric_den=>A_metric_nrm" (so that "0=Air,255=Objects").
  //                    Expecting profile-shapes:
  // ~Obj:    _____                  :                            :
  //     :   /     \               ~0:  /\     /\                 :  /\_____/\'
  // ~Air:__/       \__              :_/  \___/  \_             ~0:_/         \_
  //          "Org"              "Sgnd_ivls","Ref_ave"              "A_metric"
  //                    Expecting 2D-views:
  // "Der": +++0    "Der_lap": +++-  "Ref_ave": -0+-   "A_metric_den: 00 80 FF 00
  //      : 00+0  "Sgnd_ivls": --+-           : -0+0   "A_metric_nrm: 00 80 FF 80
  //      : 00+0             : --+-           : -0+-                : 00 80 FF 00
  //                    Notes about ~ "the constants":
  //   - cDerDec serves as <instead_of> local average of |Der3x3|. Choose it so
  // that "Der_lap>0" only in Edges (i.e. "Der_lap<0" in Air,Objects).
  //   - cLim constraints Sgnd_ivls. An accounting of "|Sgnd_ivls|>|cLim|" risk
  // CLAHE's overflow. Choose cLim so that (statistically approx):
  // "Sgnd_ivls(Air_average) < -cLim < Sgnd_ivls(Air_mx) < Sgnd_ivls(Objects_mn)".
  //   - cK is attenuator "K:=cK*Org_gamut/cLim"; cK~=0.2..0.3.
  void iAirEstimation_Lims() { char qS[256]; // 'OrgPtr' shall not be modifed until the end
    const int   cWo = rW, cHo = rH, cNo = cWo*cHo,
                cPyr = 2, cWp = cWo>>cPyr, cHp = cHo>>cPyr, cNp = cWp*cHp;
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "A.raw" ); BS_Write( rDbgPfx, cNo, cOrgPtr ); }; // trace 'Org'
    { int qZo = 0, qZn = 0; // sentinels ~"here shall be ~<25% of zeros ..."
      for (int qi=0;qi<cNo;qi++) { int qV = cOrgPtr[qi]; if (!qV) qZo++; else if (qV<0) qZn++; }
      qZo = iDiv32_1000( qZo, cNo ); qZn = iDiv32_1000( qZn, cNo ); // <in 0.1%>
      if (rDbgTrc) { sprintf( qS, "Zo=%.2f%% Zn=%.2f%%", qZo*0.1f, qZn*0.1f ); iW( qS ); };
      if (qZo>250 && qZn<50) { // adhoc sentinel_A 
        tBS_Val *qB = rBkgPtr = new tBS_Val[cNo]; BS_ZFil( cNo, qB );
        for (int qi=0;qi<cNo;qi++) { int qV = cOrgPtr[qi]; if (!qV) continue; qB[qi] = qV>0?256:-256; }; 
        tBS_Val *qT = new tBS_Val[cNo]; BS2_MovAve2Dp( ixAL.rRad%100, qB, &qT, ixAL.rRad/100, cWo, cHo ); delete qT;
        if (1 && rDbgSfx) { strcpy( rDbgSfx, "B.raw" ); BS_Write( rDbgPfx, cNo, qB ); };
        if (rDbgTrc) iW( "WRN estimation by adhoc sentinel_A" ); return; }; // if adhoc sentinelA
    }; // of sentinel
    tBS_Val    *qD0[2]; BS_DatasAlloc( mcASz(qD0), cNo, qD0 );
    tBS_Val    *qO = qD0[0], *qT = qD0[1]; BS_Cpy( cNo, cOrgPtr, qO );
    BS_ScaleDn( cPyr, cWo, cHo, qO, qT ); // #0:downsampled OrgPtr
    tBS_Val    *qDp[6]; BS_DatasAlloc( mcASz(qDp), cNp, qDp ); // frames' buffer
    { // -- preparations
      BS_DataDup( 3, cNp, qO, qDp );
      BS2_Gamma( ixAL.rGamma, cNrmGa, qDp[0], 0, cNp ); // Dp#0:Org_gamma
      BS_VMul( cNp, qDp[1], 255.0f/(cNrmGa-1) ); }; // Dp#1:Org_8bit
    { // -- step a) Der_lap:=|Der3x3|-cDerDec 
      BS_DrA3x3( cWp, cHp, qDp[2] ); // Dp#2:Der_lap
      BS_VAdd( cNp, qDp[2], -abs(ixAL.rDerDec) ); };
    int         qCntLim;
    { // -- step b) Sgnd_ivls:=directional_counts(Der_lap). 
      tBS_Val **qLimbo = qDp+2; // Dp#3:Sgnd_ivls
      const int cNDD = 4, cDirs[cNDD] = {1,3,5,7}; //const int cNDD = 4, cDirs[cNDD] = {2,4,6,8};
      BS_DataDup( cNDD-1, cNp, qLimbo[0], qLimbo+1 ); 
      BS2_Met( qLimbo, cBS_MetC|cBS_MetByPitch, NULL, cDirs, cNDD, cWp, cHp );
      BS_CPB( mccpb_Ave<cNDD>, cNDD|cCPBNDD, cNp, qLimbo );
      float      qadj = 1.0f; // adjust resulting estimation of Air
      if (rcB<8) qadj = mcIpl(7-rcB,6,0.9f,0.5f); else // increase Air
      if (rcB>8) qadj = mcIpl(rcB-9,6,1.1f,2.0f); // decrease Air
      qCntLim = (int)(ixAL.rCntLim*qadj+0.4999f);
      BS_Rng( cNp, qDp[2], -qCntLim, qCntLim ); // Dp#2:Sgnd_ivls_lim
      if (rDbgTrc) { sprintf( qS, "adj(t%Xcs):=%.2f CntLim(%d):=%d",
        rcB, qadj, ixAL.rCntLim, qCntLim ); iW( qS ); }; };
    { // -- step c) Ref_ave:=_CLAHE_TraRef(Org_8bit,Sgnd_ivls).
      const int cBlr = 3; // Dp#3:Ref_ave
      BS4_CLAHE_TraRefAir( ixAL.rBlk, cBlr, qDp[1], qDp[2], qDp+3, cWp, cHp ); };
    { // -- step d) A_metric_den:=Org+K*Ref_ave.
      const float cK = 0.44f*cNrmGa/qCntLim; // hardcoded attenuation of linear combination 
      BS_VMul( cNp, qDp[3], cK ); // Dp#3: *=K
      BS_Cpy( cNp, qDp[3], qDp[4] ); BS_Add( cNp, qDp[0], qDp[4] ); // Dp#4:A_Metric_den
      if (rDbgTrc) { sprintf( qS, "K(%d):=%.2f", qCntLim, cK ); iW( qS ); }; };
    const int cTg0 = ixAL.rcuTg0; // Dp#5:A_Metric
    { // -- step e) Normalize A_metric_den=>A_metric_nrm.
      const tBS2_Cut cCut = {{ixAL.rcuBlack,ixAL.rcuWhite,257-cTg0},{0,0}}; 
      BS_Cpy( cNp, qDp[4], qDp[5] ); BS2_Cut( &cCut, qDp[5], cNp ); if (cTg0) BS_VAdd( cNp, qDp[5], cTg0 );
      int       qBlrRad = (ixAL.rRad%100)<<(ixAL.rRad/100-cPyr); // Dp#5:blur(A_Metric)
      BS2_MovAve2D( __max(qBlrRad,1), qDp[5], qT, cWp, cHp ); BS_Cpy( cNp, qT, qDp[5] ); };
    if (1 && rDbgSfx) for (int qi=0;qi<mcASz(qDp);qi++) { // 'if debug_dump'
      // #0:Org_gamma #1:Org_8bit #2:Der_lap_lim #3:Ref_ave*K #4:A_Metric_den #5:A_Metric_nrm
      BS_Cpy( cNp, qDp[qi], qO ); BS_ScaleUp( cPyr, cWp, cHp, qO, qT, false );
      sprintf( rDbgSfx, "B%d.raw", qi ); BS_Write( rDbgPfx, cNo, qO ); }; // for i
    if (1) { // -- assignment of BkgPtr
      rBkgPtr = new tBS_Val[cNo]; BS_Cpy( cNp, qDp[5], rBkgPtr ); 
      BS_ScaleUp( cPyr, cWp, cHp, rBkgPtr, qT, false ); }; // -- end of assignment of BkgPtr
    if (1) { // -- numeric outputs
      const tBS_Val *cMnrm = qDp[5];
      int       qC = 0; for (int qi=0;qi<cNp;qi++) if (cMnrm[qi]<cTg0*2) qC++;
      rAir[0] = 0.25f*cNp/(qC+1); 
      rAiA[0] = (int)(rAir[0]*cNrmGaLg);
      if (rDbgTrc) { sprintf( qS, "est %d..%d (%.3f..%.3f)",
        rAiA[0], rAiA[1], rAir[0], rAir[1] ); iW( qS ); }; }; // -- end of numeric outputs
    BS_DatasRelease( mcASz(qD0), qD0 ); // -- release buffers
    BS_DatasRelease( mcASz(qDp), qDp );
  } // end of AirEstimation_Lims

  //                     Teeth's air-estimation in short
  // Soft-tissues (background)  |========== occurs =======||=== occurs_integrated ===|
  // are dark but have ~random              /^^^^^^^\              ________-------^^^  
  // derivations. The metric     _^_/^^^^^^^         \_^__  _------ 
  // is # of occurences.         Air     Otsu1       Impla  Air     Otsu1       Impla
  // Otsu1 denotes ~'soft cortical bones'. Necessary note: commonplace 'unsigned Otsu
  // criterion' randomly confuses ~'Otsu1 vs. Impla'. This code uses signed numbers
  // instead of ~'silly multiple Otsu-thresholds'.
  
  static int iBOtsuXmx( const int *AY, int AXI, int AY0, int AY1 ) { int qM = AY0-AY[0], qXm = 0; // robust signed Otsu
    for (int qx=0;qx<AXI;qx++) { int qm = AY[qx]-mcIpl(qx,AXI,AY0,AY1); if (qM<qm) { qM = qm; qXm = qx; }; }; return qXm; }
  static int iBOtsuXmn( const int *AY, int AXI, int AY0, int AY1 ) { int qM = AY0-AY[0], qXm = 0; // robust signed Otsu
    for (int qx=0;qx<AXI;qx++) { int qm = AY[qx]-mcIpl(qx,AXI,AY0,AY1); if (qM>qm) { qM = qm; qXm = qx; }; }; return qXm; }
  static int iBAveSub( int AN, int *AHist ) { if (AN<=0) return 0; int qAve = mcAve<int>( AN, AHist ); // calc & subtract segment's average
    for (int qi=0;qi<AN;qi++) AHist[qi] -= qAve; return qAve; }
  static inline bool mcBIvl( int L, int LenMn, int x0, int *I0, int *I1 ) {
    if (L<LenMn) return false; *I0 = x0; *I1 = x0+L; return true; }
  static int iBLftPosIvl( const int *AY, int AX1, int ALenMn, int *AI0, int *AI1 ) { // 1st positive interval on the left
    int qL = 0; for (int qx=0;qx<AX1;qx++) if (AY[qx]>0) qL++; else if (qL) {
      if (mcBIvl( qL, ALenMn, qx-qL+1, AI0, AI1 )) return 1; qL = 0; };
    if (mcBIvl( qL, ALenMn, AX1-qL+1, AI0, AI1 )) return 2; return 0; }
  static int iBRgtPosIvl( const int *AY, int AX0, int AX1, int ALenMn, int *AI0, int *AI1 ) { // 1st positive interval on the right
    int qL = 0; for (int qx=AX1;qx>AX0;qx--) if (AY[qx]>0) qL++; else if (qL) {
      if (mcBIvl( qL, ALenMn, qx, AI0, AI1 )) return 1; qL = 0; };
    if (mcBIvl( qL, ALenMn, AX0, AI0, AI1 )) return 2; return 0; }

  static void iFillPTCIvl( int *APTC, int AX0, int AX1, int AY0, int AY1 ) {
    AX1 -= AX0; for (int qx=0;qx<AX1;qx++) APTC[qx+AX0] = mcIpl( qx, AX1, AY0, AY1 ); }
  void iiDbgRowMark( tBS_Val *ARow, int AOGa, int AMark ) { // for debug-trace only
    const int cWo1 = rW-1; int qCen = AOGa*cWo1>>12; qCen = mcRng( qCen, 2, cWo1-2 );
    for (int qi=qCen-2;qi<=qCen+2;qi++) ARow[qi] = AMark; }
  
  void iAirEstimation_Teeth() { char qS[256]; // 'OrgPtr' shall not be modifed until the end
    const int   cWo = rW, cHo = rH, cNo = cWo*cHo,
                cPyr = 1, cWp = cWo>>cPyr, cHp = cHo>>cPyr, cNp = cWp*cHp,
                cOGa = 1<<12, cOGaA = BS_Align( cOGa+1 ); // aligned size of a histogram's buffer
    tBS_Val    *qD0[2]; BS_DatasAlloc( mcASz(qD0), cNo, qD0 );
    tBS_Val    *qOp = qD0[0], *qT = qD0[1]; BS_Cpy( cNo, cOrgPtr, qOp ); // #0:downsampled OrgPtr
    BS_ScaleDn( cPyr, cWo, cHo, qOp, qT ); mcShiftRgt( cNp, qOp, cNrmGaLg-12 ); 
    const int   cMgnX = 3, cMgnY = 24, // the margins
                cWm = cWp-cMgnX*2, cHm = cHp-cMgnY*2, cNm = cWm*cHm; // marginal region of histogram
    int        *qH0 = new int[cOGaA*3], // allocate histograms' buffers; occurs (i.e. "the counts") in H#0
               *qH1 = qH0+cOGaA, *qH2 = qH1+cOGaA; // occurs_integrated in H#1; occursDiff in H2
    { // -- loop through the 'basic' histograms in H#0 ("occurs"), H#1 ("occursIg")
      BS_GamHist( cWm, cHm, qOp+cMgnY*cWp+cMgnX, cWp, (unsigned*)qH0, cOGa );
      const float cScaleOccu = 256.0f*cOGa/cNm, cScaleOccuIg = 4096.0f/cNm;
      tLL      qOccuIg = 0; for (int qi2=0;qi2<cOGa;qi2++) { int qOccu = qH0[qi2]; qOccuIg += qOccu;
        qH0[qi2] = (int)(qOccu*cScaleOccu+0.4999f); qH1[qi2] = (int)(qOccuIg*cScaleOccuIg+0.4999f); };
      mcMovAve<unsigned>( cOGa, 3, (unsigned*)qH0, (unsigned*)qH2 ); memcpy( qH0, qH2, cOGa*sizeof(int) );
    }; // -- end of loop through the 'basic' histograms
    int         qAir, qOtsu0, qOtsu1, qOtsu2, qImpla; // the marks (within gamut OGa=0..4096)
    { // -- calculate marks
      // marks step a) robust "modified Otsu criterion" in occuIg
      qOtsu1 = iBOtsuXmn( qH1, cOGa, 0, 4096 ); // ~cortical bones
      memcpy( qH2, qH0, cOGa*sizeof(int) ); // divide <the gamut on two sub-intervals> & conquer... 
      int       qAveLft = iBAveSub( qOtsu1, qH2 ),
                qAveRgt = iBAveSub( cOGa-qOtsu1, qH2+qOtsu1 );
      // marks step b) ~1st left peak of occurences (determine Air,Otsu0)
      switch (iBLftPosIvl( qH2, qOtsu1, 33, &qAir, &qOtsu0 )) { // on the left
        case 0: { rDbgInfFlg |= 0x10; // no long enough positive interval
          qAir = iBOtsuXmn( qH1, qOtsu1, 0, qH1[qOtsu1] );
          qOtsu0 = iBOtsuXmx( qH1, qOtsu1, 0, qH1[qOtsu1] );
          if (qAir>qOtsu0) mcSwap( int, qAir, qOtsu0 );
        }; break;
        case 2: { rDbgInfFlg |= 0x20; // 1st positive interval is at the right end
          qOtsu0 = (int)(qAir*0.2f);
          qAir = (int)(qAir*0.1f);
        }; break;
        default: if (qAir>678) { rDbgInfFlg |= 0x30; // too far from the left
          qOtsu0 = qAir;
          qAir = iBOtsuXmn( qH1, qOtsu0, 0, qH1[qOtsu0] );
        } else { // most normal
        }; }; // switch
      // marks step c) ~1st right peak of occurences
      switch (iBRgtPosIvl( qH2, qOtsu1, cOGa-1, 11, &qOtsu2, &qImpla )) { // on the right
        case 0: { rDbgInfFlg |= 1; // no long enough positive interval
          qImpla = 4095;
          qOtsu2 = (qOtsu1+qImpla)/2;
        }; break;
        case 2: { rDbgInfFlg |= 2; // 1st positive interval is at the left end
          qOtsu2 = (qOtsu1*1+qImpla*4)/5;
          qImpla = (qOtsu2*1+4095*1)/2;
        }; break;
        default: if (qImpla<3210) { rDbgInfFlg |= 3; // too far from the right
          qImpla = 4095;
          qOtsu2 = (qOtsu1+qImpla)/2;
        } else { // most normal 
          qOtsu2 = (qImpla*3+4095*1)/4;
          qImpla = (qOtsu2*1+4095*1)/2;
        }; }; // switch
      if (rDbgTrc) { // 'if debug_dump' (trace text details here)
        sprintf( qS, "%04d..%04d:  air ots0 ots1 ots2 impl DbgInfFlg=%X", qAveLft, qAveRgt, rDbgInfFlg ); iW( qS );
        float qK = (float)(cWo-1)/(cOGa-1);
        sprintf( qS, "row0..%04d: %04d %04d %04d %04d %04d", cWo,
          (int)(qAir*qK), (int)(qOtsu0*qK), (int)(qOtsu1*qK), (int)(qOtsu2*qK), (int)(qImpla*qK) ); iW( qS );
        sprintf( qS, "pix0..%04d: %04d %04d %04d %04d %04d", cOGa, 
          qAir, qOtsu0, qOtsu1, qOtsu2, qImpla ); iW( qS ); };
    }; // -- end of calculate marks
    if (rDbgSfx) { // 'if debug_dump' (trace images here)
      //               _            row#1        _____/               |  |  | 
      //  row#0       / \          occurs       /         row#2    |__|__|__|_|
      // occurs: ^___/   \__^    integral: ____/          marks: Air   Otsu?  Impla
      tBS_Val  *qQ[1]; BS_DatasAlloc( mcASz(qQ), cNo, qQ );
      BS_Cpy( cNp, qOp, qQ[0] ); BS_ScaleUp( cPyr, cWp, cHp, qQ[0], qT, false ); // p#0 in T#1
      #if 1 // gray-rescale 0..4096=>0..cWp (match "as is occurs col#" in row#0 aginst pixel values)
        const tBS2_Cut cCut = {{0,4096,cWo},0,0}; BS2_Cut( &cCut, qQ[0], cNo );
      #endif
      tBS_Val  *qDbg = qQ[0]; // map shapes H#0..H#3 to beginning rows of Dbg^ (map range 0..gamut to range 0..width)
      for (int qio=0;qio<cWo;qio++) { // occurences, occurences' integral
        int qgi = qio*(cOGa-1)/(cWo-1); // grayidx:=i*gamut/width
        qDbg[qio+cWo*0] = qH0[qgi]; qDbg[qio+cWo*1] = qH1[qgi]; // qDbg[qio+cWo*2] = qgi;
        qDbg[qio+cWo*2] = qH2[qgi];
      }; // end of of 'for i2'
      // marks of Air etc
      tBS_Val *qDbgRowOtsu = qDbg+cWo*3; BS_ZFil( cWo, qDbgRowOtsu ); iiDbgRowMark( qDbgRowOtsu, qAir, 1800 );
      iiDbgRowMark( qDbgRowOtsu, qOtsu0, 2000 ); iiDbgRowMark( qDbgRowOtsu, qOtsu1, 2100 );
      iiDbgRowMark( qDbgRowOtsu, qOtsu2, 2200 ); iiDbgRowMark( qDbgRowOtsu, qImpla, 2300 );
      strcpy( rDbgSfx, "B" ); BS_DatasWrt( rDbgPfx, cNo, qQ, 0, 1 );
      BS_DatasRelease( mcASz(qQ), qQ );
    }; // -- end of 'if debug_dump'
    { // -- assignment of BkgPtr
      int      *qPTC = qH0; // use large pre-allocated buffer 
      iFillPTCIvl( qPTC, 0,      qAir,     0,   8 );
      iFillPTCIvl( qPTC, qAir,   qOtsu0,   8,  99 );
      iFillPTCIvl( qPTC, qOtsu0, qOtsu1,  99, 256 );
      iFillPTCIvl( qPTC, qOtsu1, qOtsu2, 256, 244 );
      iFillPTCIvl( qPTC, qOtsu2, qImpla, 244, 222 );
      iFillPTCIvl( qPTC, qImpla, 4096,   222, 128 );
      CBS_XLat  qXLat; qXLat.Init( qPTC, 0, 4096 ); BS_Cpy( cNp, qOp, qT ); qXLat.Apply( cNp, qOp );
      rBkgPtr = new tBS_Val[cNo]; BS_Cpy( cNp, qOp, rBkgPtr );
      int       qBlrRad = (ixAL.rRad%100)<<(ixAL.rRad/100-cPyr); // ~const '3'
      BS2_MovAve2D( __max(qBlrRad,1), rBkgPtr, qT, cWp, cHp ); BS_Cpy( cNp, qT, rBkgPtr );
      BS_ScaleUp( cPyr, cWp, cHp, rBkgPtr, qT, false ); 
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "C.raw" ); BS_Write( rDbgPfx, cNo, rBkgPtr ); };
    }; // -- end of assignment of BkgPtr
    if (1) { // -- numeric outputs
      rAiA[0] = (qAir+qOtsu0)/2; rAiA[1] = qImpla; 
      rAir[0] = qH1[rAiA[0]]/4096.0f; rAir[1] = 1.0f-qH1[rAiA[1]]/4096.0f;
      const float cAir0Mn = 1.0e-3f; if (rAir[0]<cAir0Mn) rAir[1] = cAir0Mn;
      const float cAir1Mn = 7.5e-4f; if (rAir[1]<cAir1Mn) rAir[1] = cAir1Mn;
      rAiA[0] = rAiA[0]<<(cNrmGaLg-12); rAiA[1] = (4096-rAiA[1])<<(cNrmGaLg-12); // ~scale from ~12 to 14bits
      if (rDbgTrc) { sprintf( qS, "Air %d..%d (%.3f..%.3f)",
        rAiA[0], rAiA[1], rAir[0], rAir[1] ); iW( qS ); };
    }; // -- end of numeric outputs
    delete qH0;
    BS_DatasRelease( mcASz(qD0), qD0 );
  }; // end of AirEstimation_Teeth

}; // end of tAiWr

static tAiWr *jxAiWrSpawn1; static void jAiWrSpawn1( int ATag ) { // could be inline or spawned (via the only exemplar of tAiWr)
  if (ATag!=1) return; jxAiWrSpawn1->AirEstimation(); } 

//// exports ////
void BS5a_AirEstimationLoad( GIniFile *AF ) { switch (atoi( AF->SGetStr("nEnabled","0") )) { // -- load estimation parameters
  case  1: ixAirEstimationEnabled = 1; iAirEstimationLoad_Lims( AF ); break;
  case  2: ixAirEstimationEnabled = 2; iAirEstimationLoad_Teeth( AF ); break;
  default: ixAirEstimationEnabled = 0; }; }

int BS5a_AirEstimationEnabled() { return ixAirEstimationEnabled; } // -- ">0: enabled some or another mode"

void BS5a_AirEstimation( int AW, int AH, const tBS_Val *AFrm, const char *ADbgPfx, int AAdj ) { // -- forks air-detection thread
  static tAiWr qAiWr; memset( &qAiWr, 0, sizeof(tAiWr) ); qAiWr.DbgSet( ADbgPfx ); // this static locates in DSeg
  qAiWr.rW = AW; qAiWr.rH = AH; qAiWr.cOrgPtr = AFrm;
  qAiWr.rcTag = (AAdj>>12)&15; qAiWr.rcB = (AAdj>>8)&15; qAiWr.rcC = (AAdj>>4)&15;
  jxAiWrSpawn1 = &qAiWr;
  #if 1
    BS_Spawn1( jAiWrSpawn1, 0 ); // spawn
  #else
    jAiWrSpawn1( 1 ); jxAiWrSpawn1 = NULL; // inline execution
  #endif
}

void BS5a_AirEstimationWait( tAirEstimationWait *AAEW ) { // -- wait for forked air-detection thread
  if (!jxAiWrSpawn1) return; BS_Spawn1Wait( 0 );
  AAEW->rAiA[0] = jxAiWrSpawn1->rAiA[0]; AAEW->rAiA[1] = jxAiWrSpawn1->rAiA[1];
  AAEW->rAir[0] = jxAiWrSpawn1->rAir[0]; AAEW->rAir[1] = jxAiWrSpawn1->rAir[1];
  AAEW->rAiF[0] = jxAiWrSpawn1->rBkgPtr;
  AAEW->rDbgInfFlg = (ixAirEstimationEnabled<<8) | (jxAiWrSpawn1->rDbgInfFlg&255); }
