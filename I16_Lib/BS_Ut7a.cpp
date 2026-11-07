// I16_Ut7a.cpp : collimator-estimation (CPU-only)
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut5.h"
#include "BS_Ut7.h" // ad-hoc service CPU-transforms
#include "BS_Ut7a.h"

//// commons with the scope of local 'DSeg' ////  
static int ixCollimatorEnabled = 0; int BS7a_CollimatorEnabled() { return ixCollimatorEnabled; } // -- >0:is enabled

// Input & output for "process mskGDi by ixGDi:
void BS7a_CollimatorLoad( GIniFile *AF ) { // load estimation parameters (called by I16_Init)
  ixCollimatorEnabled = atoi( AF->SGetStr( "nEnabled", "0" ) ); if (ixCollimatorEnabled<=0) return;
}

//// collimator-esimation's core ////
// Call Preprocess.
//  = Optional downsample
//  = Determine "slightly inflated nocollimator area"; expect relevant |Lap| inside.
// stA: entry stats: fast_descend, |Lap_neg|;
//      assign LapThr (threshold of ~valuable |Lap| within the frame).
// stB: MskEst<imation> by fast_descend ("=0:coll","1:uncertain","2:nocoll")
//      MskEstBin ("=0:coll","1:nocoll");
// stC: reproject MskEstBin (set peripheral spots to zero).
// stD: restrict MskEst by reprojected MskEstBin.
// stE: fast_ascend of |Lap| within restricted MskEst.
// stF: refill ~<the margings> of fast_ascend.
// stG: calc MaskNrm by fast_ascend ("=0:coll","511:nocoll").
// stH: inflate MskNrm to smallest rectangle.
// stI: blur & crimp MskNrm.
// stJ: extrapolate Blr within "coll";
//      Exp := Msk>T ? Blr : Blr_extrapolated.
// = Optional upscale of MskNrm, Exp.
// = Replace "coll" by the Exp in Frm:
//   Frm := MskNrm>T ? Frm : Exp (synthetic org frame)
//   Lap := Frm-Blr; Frm := Blr+Lap*MskNrm
// Call core process (grid suppression, unsharp, ...).
// Call PostProcess.
//   Frm := MskNrm>T? Frm : const
struct tCollimatorWr { // wrapper-record
  public: // -- publics
  static const int cOrgGaLg = 14, cOrgGa = 1<<cOrgGaLg, cOrgGaA = (int)(cOrgGa*0.5f), // gamut of original data
                   cMskGaLg =  9, cMskGa = 1<<cMskGaLg; // gamut of the masks
  // = initialize wrapper-record
  void Init( int AWo, int AHo, tBS_Val *AFrm, tBS_Val *AMsk, tBS_Val **ATs, bool ADbg ) {
    rWo = AWo; rHo = AHo; rFrm = AFrm; rMsk = AMsk; rTs = ATs;
    rPyr = 2; rWp = rWo>>rPyr; rHp = rHo>>rPyr;
    rDbgDir = rDbgSfx = NULL; if (ADbg) { //  debug dumps
      strcpy( rDbgPfx, xBS_DirHome ); rDbgDir = strend( rDbgPfx ); *rDbgDir++ = cPSChr;
      rDbgSfx = rDbgDir+sprintf( rDbgDir, "$c%04dx%04d", rWp, rHp );
      if (BS_DirCreate( rDbgPfx )) fprintf( stderr, "rd %s\n", rDbgPfx ); *rDbgSfx++ = cPSChr;
    }; } // of Init

  // = process downsampled frame
  bool Process() {
    const int cWo = rWo, cHo = rHo, cNo = cWo*cHo,
              cPyr = rPyr, cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    tBS_Val **qD0 = rTs+2;
    tBS_Val  *qOp = qD0[0]; BS_Cpy( cNo, rFrm, qOp ); // #0:downsampled OrgPtr
    BS_ScaleDn( cPyr, cWo, cHo, qOp, qD0[1] ); BS_Rng( cNp, qOp, 1 ); // "a pix=0" might be a "helper tag"
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "AAorg.raw" ); BS_Write( rDbgPfx, cNp, qOp ); };
    tBS_Val  *qTs[20]; BS_DatasAlloc( mcASz(qTs), cNp, qTs ); // allocate temp_limbo
    char      qS[88];
    CBS_Ticks qTM; if (xDbgDmp==1) strcpy( qTM.rPfx, "Collim" ); qTM.Reset(); 
    // stA: entry stats
    tBS_Val **pLap = qTs+0; // #0:OrgBlr #1:Lap #2:|Lap_neg|
    int       qLapThr = iTheLap( qOp, pLap ); // threshold of ~valuable |Lap_neg|
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "ABlap" ); BS_DatasWrt( rDbgPfx, cNp, pLap, 0, 3 ); };
    tBS_Val **pDes = pLap+3; // #0:fast_descend iterate twice to avoid occasional rays
    iAlgoProc( iiDescend, qOp, pDes+0 ); iAlgoProc( iiDescend, pDes[0], pDes+1, cOrgGa-1 );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "ADdes" ); BS_DatasWrt( rDbgPfx, cNp, pDes, 0, 2 ); };
    BS_Cpy( cNp, pDes[1], pDes[0] ); // ignore 1st iteration of fast_descend (as we have the 2nd one)
    // stB: MskEst<imation> by fast_descend
    tBS_Val **pMsk = pDes+1; tHisStt qHisStt; // stats
    iMskByDescend( pDes[0], pLap[2], qLapThr, pLap[0], pMsk, &qHisStt ); // #0:MskEst
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "BAmsk" ); BS_DatasWrt( rDbgPfx, cNp, pMsk, 0, 1 ); };
    if (qHisStt.rCollArea100<5 || qHisStt.rCollArea100>93) {
      sprintf( qS, cTPS" WRN CO Process CollArea100=%d. No collimator.", qHisStt.rCollArea100 );
      WriteLn( qS ); return false; };
    if (qHisStt.rCollArea100>50) {
      // stC: reproject MskEstBin
      tBS_Val   **pMsR = pMsk+1;
      { tBS_Val  *qMsBin = pMsR[0]; BS_Cpy( cNp, pMsk[0], qMsBin ); // #0:MskEstBin
        for (int qi=0;qi<cNp;qi++) qMsBin[qi] = qMsBin[qi]>1?1:0; };
      iAlgoProc( iiReproHull, pMsR[0], pMsR+1, -1 ); // #1:MskEstBin_reprojected
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "Cmsr" ); BS_DatasWrt( rDbgPfx, cNp, pMsR, 0, 2 ); };
      // stD: restrict MskEst by MskEstBin_reprojected
      { tBS_Val  *qMre = pMsR[1], *qMsk = pMsk[0];
        for (int qi=0;qi<cNp;qi++) if (qMre[qi]<1) qMsk[qi] = 0; };
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "Dmsk" ); BS_DatasWrt( rDbgPfx, cNp, pMsk, 0, 1 ); };
    }; // if CollArea>...
    // ensure that mask correlates with valuable |Lap_neg|
    if (!iMskCorrelatesToLap( pMsk[0], pLap[2], qLapThr ))
      return false;
    // stE: fast_ascend
    tBS_Val **pAsc = pMsk+1;
    BS_Cpy( cNp, pLap[2], pAsc[0] ); BS_Mul( cNp, pMsk[0], pAsc[0], -1 ); // asc_inp:=Lap*Msk shr 1
    iAlgoProc( iiAscend, pAsc[0], pAsc+1 );
    iAlgoProc( iiAscend, pAsc[1], pAsc+2 );
    iAlgoProc( iiAscend, pAsc[2], pAsc+3 );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Easc" ); BS_DatasWrt( rDbgPfx, cNp, pAsc, 0, 4 ); };
    BS_Cpy( cNp, pAsc[3], pAsc[0] ); // forget about transient iterations of fast_asecnd
    if (qHisStt.rCollArea100<=33) {
      // stF: refill ~<the margings> of fast_ascend
      iRefillAscend( pMsk[0], pAsc );
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "Frfa" ); BS_DatasWrt( rDbgPfx, cNp, pAsc, 0, 3 ); };
    }; // if CollArea<...
    // stG: calc MaskNrm by fast_ascend
    tBS_Val **pMbA = pAsc+1;
    iMskByAscend( pAsc[0], pMsk[0], pLap[2], qLapThr, pLap[0], pMbA ); // #0:MskNrm
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Gmba" ); BS_DatasWrt( rDbgPfx, cNp, pMbA, 0, 1 ); };
    if (qHisStt.rCollArea100>50) {
      // stH: inflate MskNrm to smallest rectangle
      tBS_Val  **pRcA = pMbA+1;
      { tBS_Val *qM = pMbA[0]; for (int qi=0;qi<cNp;qi++) if (qM[qi]) qM[qi] = 1; };
      iTiltRect( pMbA[0], pRcA+0 ); // #0:rectangular mask
      { tBS_Val *qM = pRcA[0]; for (int qi=0;qi<cNp;qi++) if (qM[qi]) qM[qi] = cMskGa; };
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "Htrc" ); BS_DatasWrt( rDbgPfx, cNp, pRcA, 0, 1 ); };
      BS_Cpy( cNp, pRcA[0], pMbA[0] ); // replace MskNrm by the rectangle
    }; // if CollArea>...
    // stI: blur & crimp MskNrm
    tBS_Val **pCri = pMbA+1;
    BS2_MovAve2D( 5, pMbA[0], pCri[0], cWp, cHp ); // #0:MskNrm_blr
    { // threshold of MskNrm_blr for extrapolation
      tBS_Val  *qBlr = pLap[0], *qMsk = pCri[0], *qExp = pCri[1]; // #1:Exp_inp
      const int cT = (int)(cMskGa*0.99f); for (int qi=0;qi<cNp;qi++) qExp[qi] = qMsk[qi]>cT ? qBlr[qi] : 0;
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "Icri" ); BS_DatasWrt( rDbgPfx, cNp, pCri, 0, 2 ); }; };
    // stJ: extrapolate Blr within "coll"; "Exp := MskNrm_blr>T ? Blr : Blr_extrapolated"
    tBS_Val **pExp = pCri+1;
    iExtrapolate( pExp ); // #0:Exp_out
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Jexp" ); BS_DatasWrt( rDbgPfx, cNp, pExp, 0, 1 ); };
    // leave
    BS_Cpy( cNp, pCri[0], rMsk );   BS_ScaleUp( cPyr, cWp, cHp, rMsk, qOp, false ); // upscale MskNrm_blr
    BS_Cpy( cNp, pExp[0], rTs[0] ); BS_ScaleUp( cPyr, cWp, cHp, rTs[0], qOp, false ); // upscale Exp
    sprintf( qS, cTPS" inf CO Process used ~%d of %d bufs", (int)(pExp+3-qTs), mcASz(qTs) ); WriteLn( qS );
    BS_DatasRelease( mcASz(qTs), qTs ); // release temp_limbo
    qTM.Mark( "CollimatorWr_Process" );
    return true; } // of Process

  protected: // -- intrinsics
  tBS_Val *rFrm, *rMsk, **rTs; // org frame, nocoll-mask, temp buffers from caller
  int      rWo, rHo; // org frame's dimension
  int      rPyr, rWp, rHp; // the "pyramid"
  char     rDbgPfx[256], *rDbgDir, *rDbgSfx; // debug-dump's strings

  // = helper: invoke AlgoProc; optional skew-iterations
  typedef void(tAlgoProc)(tBS_Val**,int,int);
  void iAlgoProc( tAlgoProc AAlgoProc, const tBS_Val *ASrc, tBS_Val **ARes, int ASkewFillMgn=-9999, int AW=0, int AH=0 ) {
    if (!AW) { AW = rWp; AH = rHp; }; const int cN = AW*AH;
    BS_DataDup( 8, cN, ASrc, ARes ); AAlgoProc( ARes, AW, AH );
    if (ASkewFillMgn<-1) return;
    tBS_Val **qTs = ARes+1;
    const int cShear[][2] = {{-12,0},{-5,0},{0,-12},{0,-5}};
    for (int qskw=0;qskw<mcASz(cShear);qskw++) { BS_Cpy( cN, ASrc, qTs[0] ); // skew-buffer in "temps"
      BS7_SkewH( cShear[qskw][0], true, qTs[0], AW, AH ); BS7_SkewV( cShear[qskw][1], true, qTs[0], AW, AH );
      if (ASkewFillMgn>=0) { BS7_SkewFillMgn( false, qTs[0], ASkewFillMgn, AW, AH );
                             BS7_SkewFillMgn( true, qTs[0], ASkewFillMgn, AW, AH ); };
      BS_DataDup( 7, cN, qTs[0], qTs+1 ); AAlgoProc( qTs, AW, AH );
      BS7_SkewV( -cShear[qskw][1], true, qTs[0], AW, AH ); BS7_SkewH( -cShear[qskw][0], true, qTs[0], AW, AH ); 
      BS_CPB( mccpb_Mn<2>, 2|cCPBNDD, cN, ARes );
    }; } // for skw, of iAlgoProc

  // = ~fast_descend
  static int iicpl_Descend( int AL, tBS_Val *ADat ) { mcWSMAd16( AL, ADat ); return AL; }
  static void iiDescend( tBS_Val **ADat, int AW, int AH ) { const int cN = AW*AH;
    BS_CPLA( iicpl_Descend, 8, cBS_AllDirsT, AW, AH, ADat );
    BS_CPB( mccpb_Mn<8>, 8|cCPBNDD, cN, ADat ); } // mccpb_Mx<8> icpb_Descend

  // = returns threshold of ~valuable |Lap_neg|
  static void iiTheLapCombine( // combine |Lap_neg| with Lap_pos at margins
    const tBS_Val *A0, tBS_Val *A1, int AL, int AStride, int AO ) { A0 += AO; A1 += AO;
    for (int qi=0;qi<AL;qi++) { *A1 = __max( *A1, -*A0 ); A0 += AStride; A1 += AStride; }; }
  int iTheLap( const tBS_Val *AOrg, tBS_Val **ARes ) { // #0:Blr #1:Lap #2:|Lap_neg|
    char       qS[88];
    const int  cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    BS2_MovAve2D( 5, AOrg, ARes[0], cWp, cHp ); BS_Cpy( cNp, AOrg, ARes[1] ); BS_SubA( cNp, ARes[0], ARes[1] );
    BS_Cpy( cNp, ARes[1], ARes[2] ); BS_Neg( cNp, ARes[2] ); BS_Rng( cNp, ARes[2], 0 );
    for (int qi=0;qi<7;qi++) { // combine |Lap_neg| with Lap_pos at margins
      iiTheLapCombine( ARes[1], ARes[2], cWp, 1,   qi*cWp ); iiTheLapCombine( ARes[1], ARes[2], cWp, 1,   cWp*(cHp-1-qi) );
      iiTheLapCombine( ARes[1], ARes[2], cHp, cWp, qi );     iiTheLapCombine( ARes[1], ARes[2], cHp, cWp, cWp-1-qi ); };
    int        qA = 1; // ~AveEffective
    tBS_Val   *qLapNeg = ARes[2]; // |Lap_neg|
    for (int qpass=0;qpass<2;qpass++) { tLL qZ = 0; int qC = 0;
      for (int qi=0;qi<cNp;qi++) { int qV = qLapNeg[qi]; if (qV>=qA) { qZ += qV; qC++; }; };
      qA = (int)(qZ/qC); }; // for pass
    sprintf( qS, cTPS" inf CO TheLap LapThr=%d", qA ); WriteLn( qS );
    return qA; }

  // = determine initial nocoll-mask by fast_descend (involve also |Lap|)
  typedef struct { int rNoCollThr, rCollThr, rCollArea100; } tHisStt;
  void iMskByDescend( const tBS_Val *ADes, // out #0:MskEst "0:coll","1:uncertain","2:nocoll"
    const tBS_Val *ALap, int ALapThr, const tBS_Val *AOrgBlr, tBS_Val **ARes, tHisStt *AR ) { 
    // this code examines descend's "plain histogram" and "reference histogram by |Lap|".
    const int  cWp = rWp, cHp = rHp, cNp = cWp*cHp; char       qS[88];
    unsigned  *qHis = new unsigned[cOrgGa]; // the histograms of fast_descend
    // plain histogram of occurrences
    BS_GamHist( cWp, cHp, ADes, cWp, qHis, cOrgGa );
    int        qOtsuI1 = 0, qOccuZ1 = 0; { // Otsu1 for ~max_derivative
      const int cotsuI0 = (int)(cOrgGa*0.5f), cotsuIC = cOrgGa-cotsuI0; // start search OtsuI1 from OtsuI0
      int      qoccuZ0 = mcSum<unsigned>( cotsuI0, qHis ), qnp0 = cNp-qoccuZ0;
        int    qM = 0x7FFFFFFF, qZ = 0; for (int qi=0;qi<cotsuIC;qi++) { qZ += qHis[qi+cotsuI0];
        int    qoccuA = (int)((tLL)qi*qnp0/cotsuIC), qotsuD = qZ-qoccuA; // "mn(occuZ[i]-occuAscend[i])"
        if (qM>qotsuD) { qM = qotsuD; qOtsuI1 = qi; qOccuZ1 = qZ; }; }; // for i
      qOtsuI1 += cotsuI0; qOccuZ1 += qoccuZ0; }; // of Otsu1
    AR->rCollArea100 = (cNp-qOccuZ1)*100/cNp;
    int        qOtsuI3 = qOtsuI1; { // Otsu3 for ~warranted_collimator
      int      qM = -0x7FFFFFFF, qZ = qOccuZ1; for (int qi=qOtsuI1;qi<cOrgGa;qi++) { qZ += qHis[qi];
        int    qoccuA = (int)(((tLL)qi*cNp)>>cOrgGaLg), qotsuD = qZ-qoccuA; // "mx(occuZ[i]-occuAscend[i])"
        if (qM<qotsuD) { qM = qotsuD; qOtsuI3 = qi; }; }; }; // of Otsu3
    // histogram of reference
    int        qLapThr = (int)(ALapThr*0.1f);
    BS_Cpy( cNp, ALap, ARes[1] ); BS_VAdd( cNp, ARes[1], -qLapThr ); BS_Rng( cNp, ARes[1], 0, 63 );
    BS_GamHist( cWp, cHp, ADes, cWp, qHis, cOrgGa, ARes[1] );
    int        qLapI1 = cOrgGa-1+mcHistNP( -mcSum(cOrgGa,qHis)/9999-1, qHis+cOrgGa-1 );
    delete qHis; 
    sprintf( qS, cTPS" inf CO MskByDescend CollArea=%d%% Otsu=%d..%d Lap=%d",
      AR->rCollArea100, qOtsuI1, qOtsuI3, qLapI1 ); WriteLn( qS );
    { // #0:MskEst
      float cKNoCollThr, cKCollThr;
      if (AR->rCollArea100>50) { cKNoCollThr = 0.22f; cKCollThr = 0.50f; }
                          else { cKNoCollThr = 0.44f; cKCollThr = 0.14f; };
      int     qNoCollThr = AR->rNoCollThr = qOtsuI1-(int)((qOtsuI3-qOtsuI1)*cKNoCollThr),
              qCollThr = AR->rCollThr = qOtsuI1+(int)((qOtsuI3-qOtsuI1)*cKCollThr);
      tBS_Val *qMsk = ARes[0];
      for (int qi=0;qi<cNp;qi++) { int qV = ADes[qi];
        if (qV<qNoCollThr) qV = 2; else
        if (qV<qCollThr)   qV = 1; else
                           qV = 0;
        qMsk[qi] = qV; }; }; // for i
    sprintf( qS, cTPS" inf CO MskByDescend NoColl=%d Coll=%d",
      AR->rNoCollThr, AR->rCollThr ); WriteLn( qS ); 
  } // of iMskByDescend

  // = reproject(binary MskImg)
  static int iicpl_ReproHull( int AL, tBS_Val *ADat ) { ADat[0] = mcSum(AL,ADat); return 1; } // the 'tCPLProcInp'
  static bool iiprj_ReproHull( int ADirI, int AL, tBS_Val *APrj ) { // the 'tProjectionProc' for 'tangled projection'
    tBS_Val *qT = new tBS_Val[AL]; BS_MovAve1D( AL, 9, APrj, qT ); BS_Cpy( AL, qT, APrj ); delete qT; // blur Prj
    mcWSMPd16( AL, APrj ); return true; }
  static int __fastcall iicpb_ReproHull( tBS_Val *ADat ) { mcSortBAs<tBS_Val,8,2>( ADat ); return (ADat[0]+ADat[1])>>1; }
  static void iiReproHull( tBS_Val **ADat, int AW, int AH ) { // The inp is ~solid binary mask:
    const int cNp = AW*AH, cWH[2] = {AW,AH}; // For non-exotic situations, the Dirs here are _AllDirsT.
    BS7_Projection( // invoke tCPLProcInp and tProjectionProc
      iicpl_ReproHull, iiprj_ReproHull, 8, cBS_AllDirsT, ADat, cWH );
    BS_CPB( iicpb_ReproHull, cCPBNDD|8, cNp, ADat ); } // mccpb_Mn<cNDD>

  static int iicpl_Ascend( int AL, tBS_Val *ADat ) { mcWSMPa16( AL, ADat ); return AL; }
  static void iiAscend( tBS_Val **ADat, int AW, int AH ) {
    BS_CPLA( iicpl_Ascend, 8, cBS_AllDirsT, AW, AH, ADat ); BS_CPB( mccpb_Mx<8>, 8|cCPBNDD, AW*AH, ADat ); }
   // static int __fastcall icpb_Ascend( tBS_Val *ADat ) { mcSortBDe<tBS_Val,8,3>( ADat ); return (ADat[0]+ADat[1]+ADat[2])/3; }

  // = refill peripeherel ares of fast_ascend
  static int iiRefill( tBS_Val *ADat, int AI0, int AI1 ) {
    if (abs(AI0-AI1)<3) return AI0;
    int qV0 = ADat[AI0]; if (!qV0) return AI0;
    int qInc = AI0?-1:1, qI0 = AI0, qM = qV0;
    for (int qi=AI0+qInc;qi!=AI1;qi+=qInc) { int qV = ADat[qi];
      if (qM>qV) { qI0 = qi; break; }; // 1s descend
      qM = qV; }; // for i
    if (qM>qV0) { if (AI0) BS_VFil( AI0-qI0+1, ADat+qI0, qM );
                      else BS_VFil( qI0+1, ADat, qM ); };
    return qI0; }
  static int iicpl_Refill( int AL, tBS_Val *ADat ) {
    if (AL<10) return 0;
    int qI0 = iiRefill( ADat, 0, AL ); if (qI0>=AL-1) return AL;
    int qI1 = iiRefill( ADat, AL-1, qI0-1 );
    return (qI0>0) || (qI1<AL-1) ? AL : 0; }
  void iRefillAscend( const tBS_Val *AMask, tBS_Val **ADs ) { 
    const int  cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    tBS_Val   *qAsc = ADs[0]; // inp "raw fast_ascend"
    BS_Cpy( cNp, qAsc, ADs[1] ); // prepare inp of refill in #1 ("0:don't touch", ">0:try to refill") 
    { tBS_Val *qP = ADs[1]; for (int qi=0;qi<cNp;qi++) { int qV = AMask[qi];
        if (qV<1) qP[qi] = 0; else if (qV>1) qP[qi] = __max(qP[qi],1); }; }; // for i
    const int cDirs[] = {1,3,5,7}, cNDD = mcASz(cDirs);
    BS_DataDup( cNDD-1, cNp, ADs[1], ADs+2 );
    BS_CPLA( iicpl_Refill, cNDD, cDirs, cWp, cHp, ADs+1 );
    BS_CPB( mccpb_Mx<cNDD+1>, cCPBNDD|(cNDD+1), cNp, ADs ); } // #0:fast_ascend with refilled margins

  // = Src:=Src>Thr?1:0, where "Thr=Tdrk when OrgBlr=0", "=Tbrt when OrgBlr=1". 
  void iMskNrm( const tBS_Val *AOrgBlr, const tBS_Val *ASrc, int ATdrk, int ATbrt, tBS_Val **ARes ) { 
    const int  cNp = rWp*rHp; tBS_Val *qMskNrm = ARes[0];
    for (int qi=0;qi<cNp;qi++) { int qT = mcIpl( AOrgBlr[qi], cOrgGa-1, ATdrk, ATbrt ); qMskNrm[qi] = ASrc[qi]<qT?0:cMskGa-1; }; }

  // = determine nocoll-mask by fast_ascend (involve also |Lap| etc)
  void iMskByAscend( const tBS_Val *AAsc, // out #0:MskNrm ("0:coll","MskGa:nocoll")
    const tBS_Val *AMskEst, const tBS_Val *ALap, int ALapThr, const tBS_Val *AOrgBlr, tBS_Val **ARes ) { 
    char        qS[88];
    const int   cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    tBS_Val    *qLap = ARes[1];
    #if 0 // don't use Lap_ref
      const int  cThrBrt = (int)(ALapThr*1.1f),
                 cThrDrk = (int)(ALapThr*1.2f);
      if (rDbgSfx) BS_ZFil( cNp, qLap ); // zero Lap_ref for dbg trace
    #else // Lap_ref:="Lap*Msk shr 1"-"LapThr*K"
      BS_Cpy( cNp, ALap, qLap ); BS_Mul( cNp, AMskEst, qLap, -1 ); // #1:Lap_ref
      BS_VAdd( cNp, qLap, -(int)(ALapThr*0.5f) ); BS_Rng( cNp, qLap, 0 );
      unsigned *qHis = new unsigned[cOrgGa]; // the histograms of fast_descend
      BS_GamHist( cWp, cHp, AAsc, cWp, qHis, cOrgGa, qLap );
      int       qLapThr = mcHistNP( mcSum(cOrgGa,qHis)/333+1, qHis );
      delete qHis;
      const int cThrBrt = (int)(qLapThr*0.9f),
                cThrDrk = qLapThr-1;
    #endif
    sprintf( qS, cTPS" inf CO MskByAscend LapThr=%d,%d Thr0=%d,%d", ALapThr, qLapThr, cThrBrt, cThrDrk ); WriteLn( qS );
    iMskNrm( AOrgBlr, AAsc, cThrDrk, cThrBrt, ARes ); } // #0:MskNrm

  // = extend mask to rectangle
  static int iicpl_TiltRect( int AL, tBS_Val *ADat ) { BS_VFil( AL, ADat, mcSum(AL,ADat) ); return AL; }
  static int __fastcall iicpb_TiltRect( tBS_Val *ADat ) { return __min(ADat[0],ADat[1])>5 ? 1 : 0; }
  void iTiltRect( const tBS_Val *AMsk, tBS_Val **ARes ) {
    const int cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    BS_DataDup( 8, cNp, AMsk, ARes ); BS_CPLA( iicpl_TiltRect, 8, cBS_AllDirsT, cWp, cHp, ARes );
    char      qS[88]; strcpy( qS, cTPS" inf CO Rect"); char *qCh = strend(qS);;
    int       qZ0, qDir0 = 0; { // process full limbo of 8 directions to determine Z0 & Dir0 within 
      for (int qi=0;qi<4;qi++) {
        tBS_Val *qDs[] = {ARes[qi],ARes[qi+4]}; BS_CPB( iicpb_TiltRect, 2|cCPBNDD, cNp, qDs );
        int     qz = mcSum( cNp, ARes[qi] ); qCh = qCh+sprintf( qCh, " %d", qz );
        if (!qi) qZ0 = qz; else if (qZ0>qz) { qZ0 = qz; qDir0 = qi; }; }; // for i
      BS_Cpy( cNp, ARes[qDir0], ARes[0] ); WriteLn( qS ); }; // of process full limbo
    // skew around Dir0
    int       qDir2[2] = {qDir0+1,qDir0+5};
    tBS_Val **qTs = ARes+1; qCh = qS+sprintf( qS, "     skew around Dir%d", qDir0 );
    const int cShear[][2] = {{-6,-6},{-11,-11},{-16,-16},{16,16},{11,11},{6,6}};
    for (int qskw=0;qskw<mcASz(cShear);qskw++) { BS_Cpy( cNp, AMsk, qTs[0] ); // skew-buffer in "temps"
      BS7_SkewH( cShear[qskw][0], true, qTs[0], cWp, cHp ); BS7_SkewV( cShear[qskw][1], true, qTs[0], cWp, cHp );
      //sprintf( rDbgSfx, "qskw%d.raw", qskw ); BS_Write( rDbgPfx, cNp, qTs[0] ); //!!!!!
      BS_Cpy( cNp, qTs[0], qTs[1] );
      BS_CPLA( iicpl_TiltRect, 2, qDir2, cWp, cHp, qTs );
      BS_CPB( iicpb_TiltRect, 2|cCPBNDD, cNp, qTs );
      int qz = mcSum( cNp, qTs[0] ); qCh = qCh+sprintf( qCh, " %d", qz ); if (qZ0>qz) { 
        BS7_SkewV( -cShear[qskw][1], true, qTs[0], cWp, cHp ); BS7_SkewH( -cShear[qskw][0], true, qTs[0], cWp, cHp );
        qZ0 = qz; BS_Cpy( cNp, qTs[0], ARes[0] ); }; // if Z0>z
    }; // for skw
    WriteLn( qS ); } // of iRect

  // = extrapolate non-zero data from center to peripheral
  static const int cExtRad = 33;
  static bool iiExtrapolate( tBS_Val *ADat, int AI0, int AI1 ) {
    if ((abs(AI0-AI1)<cExtRad) || ADat[AI0]) return false; int qInc = AI0?-1:1;
    for (int qi=AI0+qInc;qi!=AI1;qi+=qInc) if (ADat[qi]) {
      ADat += qi; int qV0 = *ADat, qL1 = AI0?AI0-qi:qi, qL0 = __min(cExtRad,qL1);
      for (int qj0=0;  qj0<qL0;qj0++) { ADat -= qInc; *ADat = mcIpl(qj0,cExtRad,qV0,cOrgGaA); };
      for (int qj1=qL0;qj1<qL1;qj1++) { ADat -= qInc; *ADat = cOrgGaA; };
      return true; }; return false; } // if Dat[i]
  static int iicpl_Extrapolate( int AL, tBS_Val *ADat ) { if (AL<cExtRad*1) return 0;
    bool q0 = iiExtrapolate( ADat, 0, AL ), q1 = iiExtrapolate( ADat, AL-1, -1 );
    return q0 || q1 ? AL : 0; }
  static int __fastcall iicpb_Extrapolate( tBS_Val *ADat ) { // ~sum of non-zero entries
    tBS_Val qDat[8]; int qC = 0, qZ = 0;
    for (int qi2=0;qi2<8;qi2++) { int qV = ADat[qi2]; if (qV<1) continue; qZ += qV; qDat[qC++] = qV; };
    if (qC<1) return cOrgGaA; if (qC<4) return qZ/qC;
    int qMn = qDat[0], qMx = qMn;
    for (int qi4=0;qi4<qC;qi4++) { int qV = qDat[qi4]; if (qMn>qV) qMn = qV; else if (qMx<qV) qMx = qV; };
    return (qZ-qMn-qMx)/(qC-2); }
  void iExtrapolate( tBS_Val **ADs ) {
    const int cWp = rWp, cHp = rHp, cNp = cWp*cHp;
    const int *cDirs = cBS_AllDirsT, cNDD = 8;
    BS_DataDup( cNDD-1, cNp, ADs[0], ADs+1 );
    BS_CPLA( iicpl_Extrapolate, cNDD, cDirs, cWp, cHp, ADs );
    BS_CPB( iicpb_Extrapolate, cNDD|cCPBNDD, cNp, ADs ); } // #0:extrapolations' average

  // = ensure that mask correlates with valuable |Lap_neg|
  static void itppr_MskCorrelatesToLap( tBS_SpawnThreadRec *AP ) { // invokes LnSpawnProc for each line
    tCollimatorWr *qSelf = (tCollimatorWr*)AP->rP[0]; const int cWp = qSelf->rWp, cHp = qSelf->rHp;
    unsigned *qHiss = (unsigned*)AP->rP[1]; for (int qh=AP->rFst;qh<AP->rLst;qh++)
      BS_GamHist( cWp, cHp, (tBS_Val*)AP->rP[2], cWp, qHiss+qh*4, 3, qh?(tBS_Val*)AP->rP[3]:NULL ); }
  bool iMskCorrelatesToLap( const tBS_Val *AMsk, const tBS_Val *ALap, int ALapThr ) {
    unsigned qHiss[8] = {0}; BS_SpawnThreads( // occurences in #0, reference sums in #1
      2, 0, 2, itppr_MskCorrelatesToLap, this, qHiss, AMsk, ALap );
    const unsigned *cHis = qHiss+0, *cRef = qHiss+4; char qS[88];
    int qMet = (int)((tLL)cRef[2]*cHis[0]/((tLL)cRef[0]*cHis[2]+1));
    if (qMet<14) { sprintf( qS, cTPS" WRN CO Process MskCorrelatesToLap:=%d. No collimator.", qMet );
      WriteLn( qS ); return false; };
    sprintf( qS, cTPS" inf CO Process MskCorrelatesToLap:=%d.", qMet );
    WriteLn( qS ); return true; }

}; // end of tCollimatorWr

bool BS7a_CollimatorPreProcess( int AW, int AH, tBS_Val *AFrm, // collimator-estimation
  tBS_Val *AMsk, tBS_Val **ATs, const char *ADbgPfx ) { const int cN = AW*AH;
  tCollimatorWr qCoWr; memset( &qCoWr, 0, sizeof(tCollimatorWr) );
  qCoWr.Init( AW, AH, AFrm, AMsk, ATs, ADbgPfx!=NULL );
  if (!qCoWr.Process()) return false;
  const int cT = (int)(tCollimatorWr::cMskGa*0.8f);
  tBS_Val *qExp = ATs[0]; for (int qi=0;qi<cN;qi++) if (AMsk[qi]<cT) AFrm[qi] = qExp[qi];
  tBS_Val *qBlr = ATs[1]; BS_Cpy( cN, AFrm, qBlr );
  BS2_MovAve2Dp( 3, qBlr, ATs+2, 1, AW, AH ); BS_SubA( cN, qBlr, AFrm );
  BS_Mul( cN, qExp, AFrm, -tCollimatorWr::cOrgGaLg );
  BS_Add( cN, qBlr, AFrm );
  if (1 && ADbgPfx) { tBS_Val *qTs[] = {AFrm,AMsk,ATs[0]};
    BS_DatasWrt( ADbgPfx, cN, qTs, 0, mcASz(qTs) ); }; 
  return true; } 

void BS7a_CollimatorPostProcess( int AW, int AH, tBS_Val *AFrm, // collimator-fixation
  tBS_Val *AMsk, tBS_Val **ATs, const char *ADbgPfx ) { const int cN = AW*AH;
  const int cT = (int)(tCollimatorWr::cMskGa*0.9f);
  for (int qi=0;qi<cN;qi++) if (AMsk[qi]<cT) AFrm[qi] = tCollimatorWr::cOrgGaA;
  if (1 && ADbgPfx) { tBS_Val *qTs[] = {AFrm};
    BS_DatasWrt( ADbgPfx, cN, qTs, 0, mcASz(qTs) ); }; } 
