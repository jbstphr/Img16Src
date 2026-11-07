// I16_Ut7b.cpp : ribs-detection 
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_GIni.h"
#include "BS_Ut0.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut5.h"
#include "BS_Ut6.h" // CPU|CUDA ad-hoc utilities (with explicit dimension)
#include "BS_Ut7.h" // ad-hoc service CPU-transforms
#include "BS_Ut7b.h" // rib-cage core
#include "BS_Ut7c.h" // rib-cage core's services

//// commons with the scope of local 'DSeg' ////  
static int ixRDEnabled = 0; int BS7b_RibsDetectionEnabled() { return ixRDEnabled; } // -- >0:is enabled
static struct { //
  float fLRL[3]; // 
} ixRD = {
  {1.45f,0.2f,0.33f}
};

void BS7b_RibsDetectionLoad( GIniFile *AF ) { // load estimation parameters (called by I16_Init)
  ixRDEnabled = atoi( AF->SGetStr( "nEnabled", "0" ) ); if (ixRDEnabled<=0) return;
}

static const int
  cLftRgtLungsInp = 3, // #inp Dp[] of ProcessLftRgtLung
  cLftRgtLungsOut = 5, // #out Dp[] of ProcessLftRgtLung
  cLftRgtLungsAll = __max(cLftRgtLungsInp,cLftRgtLungsOut), // io buffers of ProcessLftRgtLung
  cLftRgtLungsAve = 8, // inp Ap[] of ProcessLftRgtLung
  cLftRgtFROI     = cLftRgtLungsOut+1, // iobuffers of ProcessLftRgtFROI
  cMaskLg         = 12, // mask "1 shl..." in this code
  cLargeBlrRad    = 20; // large blur radius (within downsampled lung area)

// --  transient buffers summary
// -- ProcessLftRgtLung_inp: #0:Org #1:AveAve #2:Mask
// --                  _out: #0:removed_ribs #1:ribs_splt #2:Mask #3:fast_descend #4:add_aux
// -- ProcessLftRgtFROI_inp: #0:Org ProcessLftRgtLung_out
// --                  _out: #0:Org-ribs

// -- the flow of split masks
// -- LRLung_MaxMagRefine_inp: maxmag
// --                    _out: maxmag_rf mm_sbfdet mm_sbfmul
// -- LRLung_MMAmbiToShrp_inp: maxmag_rf mm_sbfdet mm_sbfmul
// --                    _out: maxmag_shrp
// -- LRLung_RRAmbiToShrp_inp: removed_ribs_ambi mm_sbfdet; RRAmbiToShrp(InvCas=0)
// --                    _out: removed_ribs ribs_splt
// -- ProcessLftRgtFROI
// --   RRAmbiToShrp(InvCas=1)_inp: removed_ribs ribs_splt
// --                         _out: removed_ribs
// --   RRAmbiToShrp(InvCas=2)_inp: removed_ribs ribs_splt_as_sbfdet
//   "mm_sbfdet" by maxmag; to determine split in MMAmbiToShrp,RRAmbiToShrp(InvCas=0).
//   "mm_sbfmul" by maxmag; to mul split in MMAmbiToShrp.
//   "ribs_splt" by removed_ribs; direct split in RRAmbiToShrp(InvCas=1),
//                                to re-determine split in RRAmbiToShrp(InvCas=2).

typedef struct {
  int rRibFrmPix;
} tLungStt;
 
/////////////////////////////////////////
//// detect ribs within certain lung ////
/////////////////////////////////////////
struct tLRLWr { // - LeftRightLung process' Wrapper
  public: // -- publics of LRLWr
  tLungStt rStt; // it fills the stats

  // = initialize wrapper-record
  void Init( bool ARgt, int AWh, int AHh, tBS_Val **ADh, tBS_Val **AAh, bool ADbg ) {
    memset( &rStt, 0, sizeof(tLungStt) );
    rRgt = ARgt; rWh = AWh; rHh = AHh;
    memcpy( rDh, ADh, sizeof(rDh) ); memcpy( rAh, AAh, sizeof(rAh) );
    rDbgPfx = rDbgSfx = NULL; if (ADbg) { // debug-dumps
      strcpy( rDbgBuf, xBS_DirHome ); rDbgPfx = rDbgBuf; char *qDbgDir = strend( rDbgPfx ); *qDbgDir++ = cPSChr;
      rDbgSfx = qDbgDir+sprintf( qDbgDir, "$%s%04dx%04d", rRgt?"r":"l", AWh, AHh );
      if (BS_DirCreate( rDbgPfx )) fprintf( stderr, "rd %s\n", rDbgPfx ); *rDbgSfx++ = cPSChr;
    }; }

  // = process a lung in "SplitLftRgtLungs"
  void ProcessLftRgtLung() {
    // -- ProcessLftRgtLung_inp ('Dh[LftRgtLungsInp] of Wh*Hh').
    // -- ProcessLftRgtLung_out ('Dh[LftRgtLungsOut]').
    char      qS[256];
    const int cW = rWh, cH = rHh, cN = cW*cH; tBS_Val // the identifiers of passed 'Dh[i]:
      *cOrg = rDh[0],    // inp Dh#0: the origin
      *cAveAve = rDh[1], // inp Dh#1: heavy blurred origin
      *cMask = rDh[2];   // inp Dh#2: Mask
    tBS_Val  *qP[40]; BS_DatasAlloc( mcASz(qP), cN, qP ); // local application buffers
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "AAinp" ); BS_DatasWrt( rDbgPfx, cN, rDh, 0, cLftRgtLungsInp ); };
    { int qZ = 0, qC = 0; for (int qi=0;qi<cN;qi++) { // determine ribcage's ~saturation
        int qV = cMask[qi]; if (qV>0 && qV<(1<<12)-1) { qZ += cOrg[qi]; qC++; }; }; // for
      rStt.rRibFrmPix = qC>cN/200 ? (int)(qZ*1.45f/qC) : 9999; };
    // = estimate primary ribs (posterior arcs /-\)
    tBS_Val **pP = qP+0; BS_Cpy( cN, cOrg, pP[0] );
    jLRLung_InitialConvex( rStt.rRibFrmPix, cMask, pP[0] ); // P#0:the_convex
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "ABcvx.raw" ); BS_Write( rDbgPfx, cN, pP[0] ); };
    tBS_Val **pXcsSlo = qP+1; tBS7_AscArcFrmPf qAFP; // #0:excess_slow
    jLRLung_ExcessSlow( cMask, pP[0], pXcsSlo, &qAFP ); // "C.raw"
    if (0 && rDbgSfx) { tBS_Val **pTmp = pXcsSlo+1; BS_DataDup( 2, cN, cOrg, pTmp ); // print summary comparative results
      BS_SubA( cN, pXcsSlo[0], pTmp[1] );
      strcpy( rDbgSfx, "CY" ); BS_DatasWrt( rDbgPfx, cN, pTmp, 0, 2 ); };
    tBS_Val **pPriRib = pXcsSlo+1; // primary_rib's lap #0:pri_tomm #1:pri_torfn #2:mgncon #3:crd
    jLRLung_PrimaryRibs( cOrg, cMask, pXcsSlo[0], pPriRib ); // "F.raw"
    if (1 && rDbgSfx) { tBS_Val **pTmp = pPriRib+4; BS_DataDup( 2, cN, cOrg, pTmp ); // print summary comparative results
      BS_SubA( cN, pPriRib[1], pTmp[1] );
      strcpy( rDbgSfx, "FY" ); BS_DatasWrt( rDbgPfx, cN, pTmp, 0, 2 ); };
    tBS_Val **pMxMag = pPriRib+4; // #0:maxmag (~blurred primary_ribs)
    jLRLung_MaxMag( cMask, pPriRib[0], pMxMag, &qAFP ); // "G.raw"
    tBS_Val **pMMrfn = pMxMag+1; // out#0:ribs_rfn #1:maxmag_rfn #2:maxmag_sms #3:mgncon_rfn
    // <as we're downsampled now> then edges' refinement_approaches cost less (being imperfecrt)
    jLRLung_MaxMagRefine( pPriRib[1], pPriRib[2], pPriRib[3], pMxMag[0], pMMrfn, &qAFP ); // "H.raw"
    if (1 && rDbgSfx) { tBS_Val **pTmp = pMMrfn+4; BS_DataDup( 3, cN, cOrg, pTmp ); // print summary comparative results
      BS_SubA( cN, pMMrfn[0], pTmp[1] ); BS_SubA( cN, pMMrfn[1], pTmp[2] );
      strcpy( rDbgSfx, "HY" ); BS_DatasWrt( rDbgPfx, cN, pTmp, 0, 3 ); };
    tBS_Val **pMMshr = pMMrfn+4; // #0:maxmag_shrp
    jLRLung_MMAmbiToShrp( pMMrfn[0], pMMrfn[1], pMMrfn[2], pMMshr ); // "I.raw"
    if (1 && rDbgSfx) { tBS_Val **pTmp = pMMshr+1; BS_DataDup( 3, cN, cOrg, pTmp ); // print summary comparative results
      BS_SubA( cN, pMMshr[0], pTmp[1] );
      strcpy( rDbgSfx, "IY" ); BS_DatasWrt( rDbgPfx, cN, pTmp, 0, 2 ); };
    // = estimate clavicle (also anterior down-arcs)
    tBS_Val **pCinp = pMMshr+1; // #0:initial clavicle's shape #1:it's convex 
    BS_Cpy( cN, cOrg, pCinp[0] ); BS_SubA( cN, pMMshr[0], pCinp[0] );
    BS_Cpy( cN, pCinp[0], pCinp[1] ); jLRLung_InitialConvex( rStt.rRibFrmPix, cMask, pCinp[1] );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Ocinp" ); BS_DatasWrt( rDbgPfx, cN, pCinp, 0, 2 ); };
    tBS_Val **pClav = pCinp+2; // #0:fast_descend #1:clavicle_xcs #2:xs_mgncon
    jLRLung_EstimateClavicle( cMask, pCinp[1], pClav ); // "P.raw"
    if (0 && rDbgSfx) { tBS_Val **pTmp = pClav+3; BS_DataDup( 1, cN, pCinp[0], pTmp ); // print clavicle's estimation
      BS_SubA( cN, pClav[1], pTmp[0] ); 
      strcpy( rDbgSfx, "OY" ); BS_DatasWrt( rDbgPfx, cN, pTmp, 0, 1 ); };
    tBS_Val **pRrib = pClav+3; // #0:removed_ribs_ambi (ribs_ambi := "maxmag_shrp+clavicle")
    jLRLung_RemoveClavicle( pPriRib[3], pCinp[0], pClav[0], pClav[1], pRrib ); // "R.raw"; OCinp0.raw=>RW0.raw
    // compare changes of 
    tBS_Val **pRshp = pRrib+1; // #0:removed_ribs #1:ribs_splt
    jLRLung_RRAmbiToShrp( cOrg, pRrib[0], pMMrfn[2], pRshp ); // "S.raw"; RW0.raw=>SV0.raw
    // pass results to the caller (point to change)
    sprintf( qS, cTPS" inf RibCage ProcessLftRgtLung(%s) used ~%d of %d bufs",
      rRgt?"rgt":"lft", (int)(pRshp+3-qP), mcASz(qP) ); WriteLn( qS );    
    BS_Cpy( cN, pRshp[0], rDh[0] ); // out Dh#0:removed_ribs
    BS_Cpy( cN, pRshp[1], rDh[1] ); // out Dh#1:ribs_splt
    BS_Cpy( cN, cMask, rDh[2] );    // out Dh#2:Mask
    BS_Cpy( cN, pClav[0], rDh[3] ); // out Dh#3:fast_descend
    BS_Cpy( cN, pXcsSlo[0], rDh[4] ); BS_VMul( cN, rDh[4], 0.2f ); // out Dh#4:add_aux
    BS_VMul( cN, pMMrfn[3], 0.33f ); BS_Add( cN, pMMrfn[3], rDh[4] );
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "Wout" ); BS_DatasWrt( rDbgPfx, cN, rDh, 0, cLftRgtLungsOut ); };
    // release buffers
    BS_DatasRelease( mcASz(qP), qP ); } // of ProcessLftRgtLung

  protected: // -- intrinsics of LRLWr
  bool     rRgt;
  int      rWh, rHh; // "half" of the "pyramid"
  tBS_Val *rDh[cLftRgtLungsAll], *rAh[cLftRgtLungsAve];
  char     rDbgBuf[256], *rDbgPfx, *rDbgSfx; // debug-dump's strings

  // = inp:Org, out:inital_convex
  void jLRLung_InitialConvex( int ARibFrmPix, const tBS_Val *AMask, tBS_Val *ADat ) {
    const int cN = rWh*rHh; BS_VAdd( cN, ADat, -ARibFrmPix ); BS_Mul( cN, AMask, ADat, -cMaskLg ); }

  // = inp:the_convex, out#0:excess_slow
  void jLRLung_ExcessSlow( const tBS_Val *AMask, const tBS_Val *AConvex, tBS_Val **ARes, tBS7_AscArcFrmPf *AAFP ) {
    char      qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH;
    tBS_Val **pPxcs = ARes+1, **qArc = pPxcs+4; // local application buffers
    const int cArc = mcASz(cBS7_AscendDirs);
    BS7c_DescendXcs( // #0:xcs, Arc:descends
      cArc, cBS7_AscendDirs, AConvex, pPxcs, qArc, cW, cH );
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "CAdes" ); BS_DatasWrt( rDbgPfx, cN, qArc, 0, cArc ); };
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "CAxcs" ); BS_DatasWrt( rDbgPfx, cN, pPxcs, 0, 1 ); };
    tBS_Val  **pPpof = pPxcs+1; 
    BS7c_ArcToMet( AMask, qArc, pPpof[0], cW, cH ); // #0:inp_of_AscendFrmPf
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "CBinp.raw" ); BS_Write( rDbgPfx, cN, pPpof[0] ); };
    if (rDbgSfx) strcpy( rDbgSfx, "CBPf" ); BS7c_AscArcFrmPf(
      pPpof[0], AAFP, 1?rDbgPfx:NULL, cW, cH ); // frame polyfit's stats
    sprintf( qS, cTPS" inf RibCage LRL ExcessSlow(%s) FrmPf:=*%.1f+%d",
      rRgt?"rgt":"lft", AAFP->rA[1], (int)AAFP->rA[0] ); WriteLn( qS );
    // roughly adjust peripheral's inclines by "A1" (SkewDescendXcs depends mostly on this)
    tBS_Val  **pPxce = pPpof;
    float      qA0 = AAFP->rA[0], qA1 = AAFP->rA[1], qAm = qA1, qInc =-0.33f; tLL qZm;
    for (int qpie=0;qpie<5;qpie++) { if (rDbgSfx) sprintf( rDbgSfx, "CD%d", qpie );
      BS7c_PfSkewDescendXcs( AConvex, qA0, qA1, pPxce+0, pPxce+1, 0?rDbgPfx:NULL, cW, cH );
      tBS_Val *qTs2[2] = { pPxcs[0], pPxce[0] }; BS_CPB( mccpb_Mn<2>, 2|cCPBNDD, cN, qTs2 );
      tLL      qz = mcSumLL( cN, pPxce[0] ); if (!qpie) { qZm = qz; qA1 += qInc; } else {
        if (qZm>qz) { qZm = qz; qAm = qA1; qInc *= 0.66f; qA1 += qInc; }
               else { qInc *= -1; qA1 = qAm+qInc; }; };
    }; // for pie
    // leave
    AAFP->rA[1] = qAm; // keep adjusted tAscArcFrmPf.A1
    BS_Cpy( cN, pPxcs[0], ARes[0] ); BS_Mul( cN, AMask, ARes[0], -cMaskLg ); // #0:excess_slow accurate (by all BS9_AscendDirs)
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "CWxcs" ); BS_DatasWrt( rDbgPfx, cN, ARes, 0, 1 ); };
    sprintf( qS, cTPS" inf RibCage LRL ExcessSlow(%s) complete. AFP.A1=%.2f",
      rRgt?"rgt":"lft", AAFP->rA[1] );
   WriteLn( qS ); }
  
  // = inp:excess_slow, out prim_rib_laps #0:pri_tomm #1:pri_torfn #2:mgncon #3:crd
  void jLRLung_PrimaryRibs( const tBS_Val *AOrg, const tBS_Val *AMask, const tBS_Val *AXcs, tBS_Val **ARes ) {
    char qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH;
    BS_Cpy( cN, AOrg, ARes[0] ); BS_SubA( cN, AXcs, ARes[0] );
    BS_Cpy( cN, ARes[0], ARes[1] ); BS_MovAveY( cLargeBlrRad, cW, cH, ARes[1] );
    BS_SubA( cN, ARes[1], ARes[0] ); // Res#0:lap
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "FBprib.raw" ); BS_Write( rDbgPfx, cN, ARes[0] ); };
    BS_Cpy( cN, ARes[0], ARes[1] ); // Res#0:pri_tomm #1:pri_torfn
    #if 0 // optionally pri_tomm += a_combination(quadratures '8-4' '1-5' '2-6')
      tBS_Val **pQua = ARes+4; for (int qi=0;qi<4;qi++) { 
        BS_Cpy( cN, rAh[(qi-1)&7], pQua[qi] ); BS_SubA( cN, rAh[(qi+3)&7], pQua[qi] ); };
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "FCAquad" ); BS_DatasWrt( rDbgPfx, cN, pQua, 0, 3 ); };
      BS_CPB( mccpb_MxA<3>, 3|cCPBNDD, cN, pQua ); // #1:weighed sum of the quadratures
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "FCBquad.raw" ); BS_Write( rDbgPfx, cN, pQua[0] ); };
      BS_VMul( cN, pQua[0], 0.5f ); BS_Add( cN, pQua[0], ARes[0] ); // pri_tomm += bias_addition
    #endif
    BS7c_MgnConSimple( ARes[0], true, ARes+2, ARes+3, cW, cH ); // #2:mgncon
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "FDmgncon.raw" ); BS_Write( rDbgPfx, cN, ARes[2] ); };
    BS_SubA( cN, ARes[2], ARes[0] ); // Res#0:pri_tomm
    // leave
    BS_Mul( cN, AMask, ARes[0], -cMaskLg ); // #0:pri_tomm
    BS_Mul( cN, AMask, ARes[1], -cMaskLg ); // #1:pri_torfn
    BS_Mul( cN, AMask, ARes[2], -cMaskLg ); // #2:mgncon
    BS_Cpy( cN, ARes[0], ARes[3] ); BS7c_CrdRgtTopToLftBtm( ARes[3], cW, cH ); // #3:crd
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "FWprib" ); BS_DatasWrt( rDbgPfx, cN, ARes, 0, 4 ); };
    sprintf( qS, cTPS" inf RibCage LRL PrimaryRibs(%s) complete", rRgt?"rgt":"lft" ); WriteLn( qS ); }

  // = inp:Lap,tAscArcFrmPf.A; out:maxLap(against variations of A[0])
  void jLRLung_MaxMag( const tBS_Val *AMask, const tBS_Val *ALap, tBS_Val **ARes, tBS7_AscArcFrmPf *AAFP ) {
    char       qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH;
    tBS_Val  **pP = ARes+2; // local application buffers
    // adjust global incline by "A0" (SkewEstimate depends mostly on this)
    float      qA0 = AAFP->rA[0], qA1 = AAFP->rA[1], qAm = qA0, qInc = 12; tLL qZm;
    for (int qpie=0;qpie<5;qpie++) {
      if (rDbgSfx) sprintf( rDbgSfx, "GB%d", qpie );
      int      ct = __min(qpie,1); // inside icpl_PfSkewMaxMag Rad~=20
      BS7c_PfSkewMaxMag( ALap, qA0, qA1, qpie==4, pP+ct, pP+2, 0?rDbgPfx:NULL, cW, cH );
      tBS_Val *qBuf = pP[2]; // to calc the absolute
      BS_Cpy( cN, pP[ct], qBuf ); BS_Abs( cN, qBuf ); tLL qz = mcSumLL( cN, qBuf ); 
      if (!qpie) { qZm = qz; qA0 += qInc; } else {
        BS_CPB( mccpb_MxA<2>, 2|cCPBNDD, cN, pP );
        if (qZm<qz) { qZm = qz; qAm = qA0; qInc *= 0.66f; qA0 += qInc; }
               else { qInc *= -1; qA0 = qAm+qInc; }; };
    }; // for pie; P#0 contains initial MaxMag
    AAFP->rA[0] = qA0; AAFP->rA[1] = qA1; // <possibly> keep adjusted tAscArcFrmPf.A
    BS_Mul( cN, AMask, pP[0], -cMaskLg ); // apply on the mask
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "GWmm.raw" ); BS_Write( rDbgPfx, cN, pP[0] ); };
    BS_Cpy( cN, pP[0], ARes[0] ); // Res#0:MaxMag_initial
    sprintf( qS, cTPS" inf RibCage LRL MaxMag(%s) complete", rRgt?"rgt":"lft" ); WriteLn( qS ); }

  // = empirically refine whole maxmag; out#0:ribs_rfn #1:maxmag_rfn #2:maxmag_sms #3:mgncon_rfn 
  void jLRLung_MaxMagRefine( const tBS_Val *ARibs, const tBS_Val *AMgnCon, const tBS_Val *ACrd, const tBS_Val *AMaxMag,
    tBS_Val **ARes, tBS7_AscArcFrmPf *AAFP ) {
    char      qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH, cNr = 4;
    tBS_Val **pP = ARes+cNr; tBS_Val const *cInps[] = {ARibs,AMgnCon,ACrd,AMaxMag};
    for (int qi2=0;qi2<mcASz(cInps);qi2++) { BS_Cpy( cN, cInps[qi2], pP[qi2] );
      BS7c_PfSkewFrm( false, pP[qi2], AAFP->rA[0], AAFP->rA[1], 0, NULL, cW, cH ); };
    if (1 && rDbgPfx) { strcpy( rDbgSfx, "HAinp" ); BS_DatasWrt( rDbgPfx, cN, pP, 0, mcASz(cInps) ); };
    tBS_Val **pR1 = pP+mcASz(cInps); // #0:ribs_rfn #1:mgncon_rfn 
    if (rDbgSfx) strcpy( rDbgSfx, "HB" ); BS7c_RibsRefine( pP[0], pP[1], pP[2], pR1, 1?rDbgPfx:NULL, cW, cH ); // "Ribs -= MgnCon"
    tBS_Val **pR2 = pR1+2; // #0:maxmag_rfn #1:maxmag_sms
    if (rDbgSfx) strcpy( rDbgSfx, "HC" ); BS7c_MaxMagRefine( pR1[0], pP[3], pR2, 1?rDbgPfx:NULL, cW, cH ); // empirically refine MaxMag
    tBS_Val **pR3 = pR2+2; // #0:maxmag_rfn 
    if (rDbgSfx) strcpy( rDbgSfx, "HD" ); BS7c_CombineRefine( pR1[0], pR2[0], pR3, 1?rDbgPfx:NULL, cW, cH ); // combine refinenements
    tBS_Val *cRess[cNr] = {pR1[0],pR3[0],pR2[1],pR1[1]}; // Res#0:#ribs_rfn #1:maxmag_rfn #2:maxmag_sms #3:mgncon_rfn 
    for (int qi8=0;qi8<cNr;qi8++) { BS_Cpy( cN, cRess[qi8], ARes[qi8] ); 
      BS7c_PfSkewFrm( true, ARes[qi8], AAFP->rA[0], AAFP->rA[1], 0, NULL, cW, cH ); 
      BS_MovAveX( 0, cW, cH, ARes[qi8] ); };
    if (1 && rDbgPfx) { strcpy( rDbgSfx, "HWmmrf" ); BS_DatasWrt( rDbgPfx, cN, ARes, 0, cNr ); };
    sprintf( qS, cTPS" inf RibCage LRL MaxMagRefine(%s) complete", rRgt?"rgt":"lft" ); WriteLn( qS ); }

  // = refine_approach of the eadges
  void jLRLung_MMAmbiToShrp( const tBS_Val *AShrp, const tBS_Val *AMaxMag,
    const tBS_Val *ASMSdet, tBS_Val **ARes ) {
    // inp:primary_ribs1,maxmag_rf,mm_sbfdet,mm_sbfmul; out:maxmag_shrp
    // Laplacians' approach "MaxMag=>DescLap" is ~easy (as the target is "known exatly").
    if (1 && rDbgSfx) strcpy( rDbgSfx, "I" ); BS7c_MMAmbiToShrp( // approach ambient shape to exact sharp one
      AShrp, AMaxMag, ASMSdet, ARes+0, ARes+1, rDbgPfx, rWh, rHh ); }

  // = inp:Convex(~"OrgConvex-maxLap");
  // out#0:fast_descend #1:clavicle_xcs #2:xs_mgncon
  void jLRLung_EstimateClavicle( const tBS_Val *AMask, const tBS_Val *AConvex, tBS_Val **ARes ) {
    // This code cares about peripehral skew-aliases from ~BS7_DirDescends.
    // Note: instead of "AMask" one could apply here ~"if AConvex[i]>=0 then dat[i] := 0"
    char      qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH;
    tBS_Val **pP = ARes+3; // local application buffers
    if (1 && rDbgSfx) strcpy( rDbgSfx, "PAdes" ); BS7c_DirDescends(
      AConvex, pP, pP+8, 0 && rDbgPfx?rDbgPfx:NULL, cW, cH ); // eight dir-descends
    if (1 && rDbgSfx) strcpy( rDbgSfx, "PBcla" ); BS7c_Clavicle(
      pP, ARes, pP+8, 1 && rDbgPfx?rDbgPfx:NULL, cW, cH ); // #0:fast_descend #1clavicle_xcs
    BS_Mul( cN, AMask, ARes[1], -cMaskLg ); // remove peripheral aliases
    BS7c_MgnConSimple( ARes[1], true, ARes+2, pP, cW, cH ); // #2:xcs_mgncon
    BS_Mul( cN, AMask, ARes[2], -cMaskLg ); // remove peripheral aliases
    BS_VMul( cN, ARes[2], 1.25f ); // adhoc multiplier
    BS_SubA( cN, ARes[2], ARes[1] ); BS_Rng( cN, ARes[1], 0 ); // excess -= excess_mgncon
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "PW"); BS_DatasWrt( rDbgPfx, cN, ARes, 0, 3 ); };
    sprintf( qS, cTPS" inf RibCage LRL EstimateClavicle(%s) complete", rRgt?"rgt":"lft" ); WriteLn( qS ); }

  // = inp:#initial clavicle's shape, fast_descend, clavicle_excess; out#0:removed_ribs_ambi
  void jLRLung_RemoveClavicle( const tBS_Val *ACrd, const tBS_Val *AClaInp, const tBS_Val *AFstDes, const tBS_Val *AXcs, tBS_Val **ARes ) {
    char      qS[88]; const int cW = rWh, cH = rHh, cN = cW*cH;
    tBS_Val **pP = ARes+1; BS_Cpy( cN, AFstDes, pP[1] );
    BS7c_CrdOfFastDescend( 0.4f, 0.01f, true, pP[1], cW, cH ); // #1:convex_mask (by fast_descend)
    BS2_MovAve2D( 3, pP[1], pP[0], cW, cH ); // #0:blur convex_mask
    BS_Cpy( cN, ACrd, pP[1] ); BS_VMul( cN, pP[1], 0.9f );
    BS_Cpy( cN, pP[0], pP[2] ); BS_SubA( cN, pP[1], pP[2] );
    BS_Rng( cN, pP[2], 0 ); BS2_MovAve2Dp( 3, pP[2], pP+3, 1, cW, cH ); // #2:coordinate_map (summary)
    BS_Cpy( cN, AXcs, pP[3] ); BS_Mul( cN, pP[2], pP[3], -cMaskLg ); // #3:clavicle
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "RB" ); BS_DatasWrt( rDbgPfx, cN, pP, 0, 4 ); };
    BS_Cpy( cN, AClaInp, ARes[0] ); BS_SubA( cN, pP[3], ARes[0] ); // Res#0:removed_ribs_ambi
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "RW" ); BS_DatasWrt( rDbgPfx, cN, ARes, 0, 1 ); };
    sprintf( qS, cTPS" inf RibCage LRL RemoveClavicle(%s) complete", rRgt?"rgt":"lft" ); WriteLn( qS ); }

  // = refine_approach of the eadges; inp:removed_ribs_ambi,mm_sbfdet; out#0:ribs_shrp,#1:ribs_splt
  void jLRLung_RRAmbiToShrp( const tBS_Val *AOrg, const tBS_Val *ARRibAmbi, const tBS_Val *AMM_SbFdet, tBS_Val **ARes ) {
    // inp:primary_ribs_lap,maxmag_rf,mm_sbfdet; out#0:removed_ribs_shrp, #1:ribs_splt
    if (rDbgSfx) strcpy( rDbgSfx, "S" ); BS7c_RRAmbiToShrp( 0, // approach ambient removed_ribs to ribs_sharp
      AOrg, ARRibAmbi, AMM_SbFdet, ARes, ARes+2, rDbgPfx, rWh, rHh ); } // invocation case #0

}; // end of LRLWr

////////////////////////////////////////////////
//// tune frame's ROI by ribs' substraction ////
////////////////////////////////////////////////
struct tLRFROIWr { // - FullROI tune' Wrapper
  public: // -- publics of LRFROIWr
  // = initialize wrapper-record
  void Init( bool ARgt, int AWr, int AHr, tBS_Val **ADr, const tLungStt *AStt, bool ADbg ) {
    rRgt = ARgt; rWr = AWr; rHr = AHr;
    memcpy( rDr, ADr, sizeof(rDr) ); memcpy( &rStt, AStt, sizeof(tLungStt) );
    rDbgPfx = rDbgSfx = NULL; if (ADbg) { // debug-dumps
      strcpy( rDbgBuf, xBS_DirHome ); rDbgPfx = rDbgBuf; char *qDbgDir = strend( rDbgPfx ); *qDbgDir++ = cPSChr;
      rDbgSfx = qDbgDir+sprintf( qDbgDir, "$%s%04dx%04d", rRgt?"xr":"xl", AWr, AHr );
      if (BS_DirCreate( rDbgPfx )) fprintf( stderr, "rd %s\n", rDbgPfx ); *rDbgSfx++ = cPSChr;
    }; }

  // = process full-scale ROI in "TuneLftRgtFROIs"
  void TuneLftRgtFROI() {
    // ProcessLftRgtFROI_inp: #0:Org ProcessLftRgtLung_out ('Dr[LftRgtFROI] of Wr*Hr')
    //                  _out: #0:Org-ribs
    char      qS[256];
    const int cW = rWr, cH = rHr, cN = cW*cH; tBS_Val // the identifiers of passed 'Dh[i]:
      *cOrg = rDr[0],    // inp Dr#0:Org
      *cRRibs = rDr[1],  // inp Dr#1:removed_ribs
      *cSplit = rDr[2],  // inp Dr#2:ribs_splt
      *cMask = rDr[3],   // inp Dr#3:Mask
      *cConvex = rDr[4], // inp Dr#4:fast_descend (the "convex")
      *cAddAux = rDr[5]; // inp Dr#5:add_aux
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Ainp" ); BS_DatasWrt( rDbgPfx, cN, rDr, 0, mcASz(rDr) ); };
    tBS_Val  *qP[16]; BS_DatasAlloc( mcASz(qP), cN, qP ); // local application buffers
    if (rDbgSfx) strcpy( rDbgSfx, "S" ); BS7c_RRAmbiToShrp( 1,
      cOrg, cRRibs, cSplit, qP+0, qP+1, 1?rDbgPfx:NULL, cW, cH ); // invocation case #1
    BS_Add( cN, cAddAux, qP[0] );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "Wout.raw" ); BS_Write( rDbgPfx, cN, qP[0] ); };
    BS_Cpy( cN, qP[0], rDr[0] ); // out Dr#0:removed_ribs
    sprintf( qS, cTPS" inf RibCage ProcessLftRgtFROI(%s) used ~%d bufs",
      rRgt?"rgt":"lft", mcASz(qP) ); WriteLn( qS );
    // release buffers
    BS_DatasRelease( mcASz(qP), qP ); } // of ProcessLftRgtFROI

  private: // -- intrinsics of LRFROIWr
  bool     rRgt;
  int      rWr, rHr; // "half_rois" of full frame
  tBS_Val *rDr[cLftRgtFROI];
  tLungStt rStt;
  char     rDbgBuf[256], *rDbgPfx, *rDbgSfx; // debug-dump's strings
}; // end of LRFROIWr

///////////////////////////////
//// ribs-detection's core ////
///////////////////////////////
struct tRiWr { // - Ribs-Detection Wrapper
  // identifiers' suffixes Wp,Hp - "pyramid's downsample" of full frame
  //                       Wh,Hh - "half" of the "pyramid"
  public: // -- publics of RiWr
  int rW, rH; // original dimension
  tBS_Val *cInps[1]; tBS_Val **rTs; // frames W*H
  char rDbgPfx[256], *rDbgDir, *rDbgSfx; // debug-dump's strings
  void DbgSet() { strcpy( rDbgPfx, xBS_DirHome ); rDbgDir = strend( rDbgPfx ); *rDbgDir++ = cPSChr; }
  void DbgDirSet( const char *APfx, int AW, int AH ) { if (!rDbgDir) return;
    rDbgSfx = rDbgDir+sprintf( rDbgDir, "$%s%04dx%04d", APfx, AW, AH );
    if (BS_DirCreate( rDbgPfx )) fprintf( stderr, "rd %s\n", rDbgPfx ); *rDbgSfx++ = cPSChr; }

  void RibsDetection() { // the procedure of ribs-detection
    // global_shape=>xcs_zerocenter=>xcscnt_sig            =>the_mask=>SplitLftRgtLungs
    //                             =>xcscnt_nrm=>repro_hull=>
    const int cW = rW, cH = rH, cN = cW*cH, cWp = cW>>cPyr, cHp = cH>>cPyr, cNp = cWp*cHp;
    DbgDirSet( "c", cWp, cHp ); 
    // allocate buffers & pyramid's downsample
    tBS_Val    *qDp[cLftRgtLungsAll]; BS_DatasAlloc( mcASz(qDp), cNp, qDp ); // sharing frames' buffer 
    tBS_Val    *qAp[cLftRgtLungsAve]; BS_DatasAlloc( mcASz(qAp), cNp, qAp ); // dir_averages of Dp#0
    tBS_Val    *qEp[20]; BS_DatasAlloc( mcASz(qEp), cNp, qEp ); // private frames' buffer 
    char        qS[256]; // adhoc string-buffer; ticks' marks are ~meaningless when one dumps debug frames
    CBS_Ticks   qTM; if (xDbgDmp==1 && !rDbgSfx) strcpy( qTM.rPfx, "RibCage" ); qTM.Reset(); 
    // downsample "Inps W*H" into "Dp Wp*Hp"
    BS_Cpy( cN, cInps[0], rTs[0] );
    BS_Cpy( cNp, BS_ScaleDnA(cPyr,cW,cH,rTs[0],rTs[1],false), qDp[0] );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "AaDp" ); BS_DatasWrt( rDbgPfx, cNp, qDp, 0, 1 ); };
    tBS_Val   **qD = qDp+1; // local application buffers
    // dir_averages (in standalone buffer for simplicity of indexing only)
    BS_DataDup( 8, cNp, qDp[0], qAp ); BS_MovAve( cLargeBlrRad, cWp, cHp, qAp );
    BS_CPB( mccpb_Ave<8>, 8|cCPBNDD, cNp, qAp, qD[0] ); // D#0:AveAve
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "AcAp" ); BS_DatasWrt( rDbgPfx, cNp, qAp, 0, mcASz(qAp) );
      strcpy( rDbgSfx, "AdD0.raw" ); BS_Write( rDbgPfx, cNp, qD[0] ); };
    // excess of "aliased ascend" against global shape
    tBS_Val   **qExcs = qEp; BS_Cpy( cNp, qDp[0], qExcs[0] );
    iRDXcsMetric( cWp, cHp, qExcs[0], qExcs+1 );// #0:global_shape=>xcs_metric
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "BbXcs" ); BS_DatasWrt( rDbgPfx, cNp, qExcs, 0, 1 ); };
    iRDXcsZeroCenter( cWp, cHp, qExcs[0], qExcs+1 ); // #0:xcs_metric=>xcs_zerocenter
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "BcXcs" ); BS_DatasWrt( rDbgPfx, cNp, qExcs, 0, 1 ); };
    tBS_Val   **qEcnt = qExcs+1; iRDXcsCnt( cWp, cHp, qExcs[0], qEcnt ); // out#1:xcscnt_sig, #2:xcscnt_nrm
      // #0:xcs_zerocenter_counts for debug; the outputs here are spatial metrics ~"-W..+W".
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "CCnt" ); BS_DatasWrt( rDbgPfx, cNp, qEcnt, 0, 3 ); };
    tBS_Val   **qEhul = qEcnt+3; iRDReproHull( cWp, cHp, qEcnt[2], qEhul ); // #0:repro_hull
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "DHul" ); BS_DatasWrt( rDbgPfx, cNp, qEhul, 0, 2 ); };
    // empirical consolidation xcscnt_sig,repro_hull => the_mask
    // CutRog.rcuWhite normalizes xcscnt_sig. CutRep.rcuWhite normalizes repro_hull.
    tBS_Val   **qEcnd = qEhul+2; iRDEmpConso( cWp, cHp, qEcnt[1], qEhul[1], qEcnd );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "ECnd" ); BS_DatasWrt( rDbgPfx, cNp, qEcnd, 0, 3 ); };
    // separate lungs (and then process them)
    // passing Dp#0:BlrAda #1:AveAve #2:mask_blur
    BS_Cpy( cNp, qEcnd[2], qD[1] );
    BS_VMul( cNp, qD[1], (1<<cMaskLg)/6 ); BS2_MovAve2Dp( 6, qD[1], qEp, 0, cWp, cHp );
    BS_DatasRelease( mcASz(qEp), qEp );
    qTM.Mark( "RibCage RibsDet_begin" );
    sprintf( qS, cTPS" inf RibCage RibsDetection_begin used ~%d of %d bufs",
      (int)(qEcnd+2-qEp), mcASz(qEp) ); WriteLn( qS );
    jSplitLftRgtLungs( cWp, cHp, qDp, qAp ); // calls ppr_ProcessLftRgtLung
    qTM.Mark( "RibCage RibsDet SplitLftRgtLungs" );
    // "SplitLftRgtLungs" released Dp, Ap and outputs RChLft, RChRgt, DhLft, DhRgt, LungSttLft, LungSttRgt
    jTuneLftRgtFROIs(); // calls ppr_TuneLftRgtROI
    qTM.Mark( "RibCage RibsDet TuneLftRgtFROIs" ); }

  protected: // -- intrinsics of RiWr
  static const int cPyr = 2,
                   cRougLg = 7,          // The 'rough' assumed to be '-Mg..Mg-1' when compared
                   cRougMg = 1<<cRougLg, // but '0..2*Mg-1' while combined.
                   cReprLg = 7,          // The 'reprojection' assumed to be '0..Mg-1' while compared
                   cReprMg = 1<<cReprLg; // & combined.

  // = for RibsDetection: global_shape=>xcs_metric (aliased shape of both lungs)
  static int iicpl_RDXcsMetric( int AL, tBS_Val *ADat ) { // trace_shader
    tBS_Val *qDat = ADat+BS_Align( AL ); BS_Cpy( AL, ADat, qDat ); // PeriToCenter ascend's excesses 
    mcWSMPa16( AL, ADat ); BS_SubA( AL, qDat, ADat ); return AL; } // ascend's excess of "combined halfs"
  static int __fastcall icpb_RDXcsMetric( tBS_Val *ADat ) { // pixel_shader
    // simplest implementation could be just "min(dir_exess)", or "sum of few smallest ones"
    // this encounts "dir_excess>emprical_threshold"; output plato is "1 shl 10"
    int cT = ADat[0], qN = 0; for (int qi=1;qi<9;qi++) if (ADat[qi]>cT) qN++; return qN<<7; }
  void iRDXcsMetric( int AWp, int AHp, tBS_Val *ADat, tBS_Val **ATs ) {
    const int cNp = AWp*AHp; BS_DataDup( 9, cNp, ADat, ATs );
    { // simplest implementation could be single "const=777..2222"
      // this calculates "emprical_threshold=func(global_shape)"
      tBS_Val *qThr = ATs[0]; const int cR = 5555; for (int qi=0;qi<cNp;qi++) {
      int qV = qThr[qi]-cR; qV = mcRng( qV, 0, cR ); qThr[qi] = mcIpl( qV, cR, 3000, 1234 ); }; // for
      if (0 && rDbgSfx) { strcpy( rDbgSfx, "BAB" ); BS_DatasWrt( rDbgPfx, cNp, ATs, 0, 1 ); }; }; 
    BS_CPLA( iicpl_RDXcsMetric, 8, cBS_AllDirsT, AWp, AHp, ATs+1 ); // trace_shader
    BS_CPB( icpb_RDXcsMetric, cCPBNDD|9, cNp, ATs ); // pixel_shader 
    if (0 && rDbgSfx) { strcpy( rDbgSfx, "BAC" ); BS_DatasWrt( rDbgPfx, cNp, ATs, 0, 1 ); };
    BS2_MovAve2D( 4, ATs[0], ADat, AWp, AHp ); } // xcs_blur

  // = for RibsDetection: xcs_metric=>xcs_zerocenter (aliased shape of both lungs)
  static int iRDXcsZeroOfs( int AL, tBS_Val *ADat ) { // calc zero_ofs; xcs_metric's plato is "1 shl 10"
    const int cXcsThr = 555, cM = 256, cC0 = (int)(cM*0.05f), cC1 = (int)(cM*0.5f)-cC0, // count Pix>XcsThr
              cZ0 = 999, cZ1 = 333; // if Cnt<C0 then zero_ofs:=cZ0; if Cnt>C1 then zero_ofs:=cZ1
    int qC = 0; for (int qi=0;qi<AL;qi++) if (ADat[qi]>cXcsThr) qC++; // count pixels>XcsThr
    qC = qC*cM/AL-cC0; qC = mcRng( qC, 0, cC1 ); ADat[0] = mcIpl( qC, cC1, cZ0, cZ1 ); return 1; }
  static int iRDXcsZeroSub( int AL, tBS_Val *ADat ) { BS_VAdd( AL-1, ADat+1, -ADat[0] ); ADat[0] = ADat[1]; return AL; } // -zero_ofs
  static void iRDXcsZeroCenter( int AWp, int AHp, tBS_Val *ADat, tBS_Val **ATs ) { // xcs_zerocenter
    BS_CPLCCols( iRDXcsZeroOfs, AWp, AHp, ADat );
    BS_Rng( AWp, ADat, 0 ); mcIplIvls<tBS_Val,0>( AWp, ADat );
    BS_MovAve1D( AWp, 9, ADat, ATs[0] ); BS_Cpy( AWp, ATs[0], ADat );
    BS_CPLCCols( iRDXcsZeroSub, AWp, AHp, ADat ); }

  // = for RibsDetection: the counts (i.e. spatial metricss) of xcs_zerocenter
  // Calculates xcscnt_sig. Then treats lungs' interims at neck & backbone:
  //   for consoludation: negate uncertainties,
  //   for reprojection: refill neg horz interims as positive (to have ~solid shape).
  static int iicpl_RDXcsCnt1( int AL, tBS_Val *ADat ) { // negate uncertain interim at neck,backbone
    const int cw = AL/6, cTcnt = 99, cAcnt = -(int)(AL*0.95f)-1; tBS_Val *qRow = ADat+(AL>>1)-(cw>>1); 
    for (int qc=0;qc<cw;qc++) if (qRow[qc]<cTcnt) qRow[qc] = cAcnt; return AL; }
  static int iicpl_RDXcsCnt2( int AL, tBS_Val *ADat ) { // refill neg interim at neck,backbone
    const int cC = AL>>1, cWb = AL/6, cWa = cWb/4; int qI0, qI1;
    const int cC0 = cC-cWb; for (int qc0=cC-1-cWa;qc0>cC0;qc0--) if (ADat[qc0]>0) {
      qI0 = qc0; goto lab0; }; return 0; lab0:
    const int cC1 = cC+cWb; for (int qc1=cC+cWa;  qc1<cC1;qc1++) if (ADat[qc1]>0) {
      qI1 = qc1; goto lab1; }; return 0; lab1:
    mcIplIvl<tBS_Val>( ADat+qI0, qI1-qI0 ); return qI1; }
  void iRDXcsCnt( int AWp, int AHp, const tBS_Val *AExcessZC, tBS_Val **ARes ) { //
    const int cNp = AWp*AHp; BS_Cpy( cNp, AExcessZC, ARes[0] );
    BS2_Me8( ARes[0], cBS_MetC|cBS_Me8MxA, NULL, NULL, AWp, AHp ); // #0:xcs_zerocenter_counts
    BS_Cpy( cNp, ARes[0], ARes[1] ); // #1:xcscnt_sig for consoludation
    BS_CPLCRows( iicpl_RDXcsCnt1, AWp, AHp, ARes[1] );
    const int cHgap0 = (int)(AHp*0.1f), cHgap1 = (int)(AHp*0.2f);
    BS_Cpy( cNp, ARes[1], ARes[2] ); // #2:xcscnt_refill for reprojection
    BS_CPLCRows( iicpl_RDXcsCnt2, AWp, AHp-cHgap0-cHgap1, ARes[2]+cHgap0*AWp ); } 

  // = for RibsDetection: reproject xcscnt_refill to <~polygon-shaped> repro_hull
  static int iicpl_RDReproHull( int AL, tBS_Val *ADat ) { ADat[0] = mcSum(AL,ADat); return 1; } // the 'tCPLProcInp'
  static bool iiprj_RDReproHull( int ADirI, int AL, tBS_Val *APrj ) { // the 'tProjectionProc' for 'tangled projection'
    tBS_Val *qT = new tBS_Val[AL]; BS_MovAve1D( AL, 11, APrj, qT ); BS_Cpy( AL, qT, APrj ); delete qT; // blur Prj
    mcWSMPd16( AL, APrj ); return true; }
  static int __fastcall iicpb_RDReproHull( tBS_Val *ADat ) { mcSortBAs<tBS_Val,8,2>( ADat ); return (ADat[0]+ADat[1])>>1; }
  void iRDReproHull( int AWp, int AHp, const tBS_Val *AMet, tBS_Val **ARes ) {
    // The inp is denormalized spatial metric as global object's shape. The inp of _Projection should:
    //   be ~solid to produce ~solid hull_like projections,
    //   have ~few LSBs to prevent overflow of BS_Val by summation in trace_shader.
    const int cNp = AWp*AHp; BS_Cpy( cNp, AMet, ARes[0] ); BS_Rng( cNp, ARes[0], 0, (AWp+AHp)>>2 );
    const float cK = 555555.5f/(mcSumLL<tBS_Val>(cNp,ARes[0])+1); // here is just one of possible ways to normalize
    BS_VMul( cNp, ARes[0], cK ); BS_Rng( cNp, ARes[0], 0, 9 ); // to prevent overflow  
    const int cWH[2] = {AWp,AHp}; // For non-exotic situations, the Dirs here are _AllDirsT.
    BS_DataDup( 8, cNp, ARes[0], ARes+1 ); BS7_Projection( // invoke tCPLProcInp and tProjectionProc
      iicpl_RDReproHull, iiprj_RDReproHull, 8, cBS_AllDirsT, ARes+1, cWH );
    BS_CPB( iicpb_RDReproHull, cCPBNDD|8, cNp, ARes+1 ); } // mccpb_Mn<cNDD>

  // = for RibsDetection: consolidate rough_combined (it combines rough_signed & rough_reproj)
  // rough_signed=xcscnt_sig is '0..2*RogMg-1' while combined, '-RogMg..RogMg-1' when compared.
  // rough_reproj=repro_hull is '0..RepMg-1' when compared & combined.
  static inline int ii_RDEmpConsoComb( int ARoug, int ARepr ) { return (ARoug<<cReprLg)|ARepr; }
  static inline void ii_RDEmpConsoUnComb( unsigned APack, int &ARoug, int &ARepr ) { 
    ARoug = (APack>>cReprLg)-cRougMg; ARepr = APack&((1<<cReprLg)-1); } // the 'rough' here becomes 'signed again'
  static int iicpl_RDEmpConsoHrz( int AL, tBS_Val *ADat ) { // horz segment of rough_combined
    int qi0; for (qi0=0;qi0<AL-(AL>>1);qi0++) { int qRough, qRepr;
      ii_RDEmpConsoUnComb( ADat[qi0], qRough, qRepr ); if (qRough>10 && qRepr>6) goto lab1; };
    BS_ZFil( AL, ADat ); return AL; lab1:
    int qi1; for (qi1=AL-1;qi1>(AL>>1);qi1--) { int qRough, qRepr;
      ii_RDEmpConsoUnComb( ADat[qi1], qRough, qRepr ); if (qRough>10 && qRepr>6) goto lab2; };
    BS_ZFil( AL, ADat ); return AL; lab2:
    BS_ZFil( AL, ADat ); BS_VFil( qi1-qi0, ADat+qi0, 1 ); return AL; } 
  static int iicpl_RDEmpConsoVrt( int AL, tBS_Val *ADat ) { // vert segment of rough_combined
    int qi0; for (qi0=0;qi0<AL-(AL>>1);qi0++) { int qRough, qRepr;
      ii_RDEmpConsoUnComb( ADat[qi0], qRough, qRepr ); if (qRough>10 && qRepr>6) goto lab1; };
    BS_ZFil( AL, ADat ); return AL; lab1:
    int qi1; for (qi1=qi0+99;qi1<AL-11;qi1++) { int qRough, qRepr;
      ii_RDEmpConsoUnComb( ADat[qi1], qRough, qRepr ); if (qRough<-20 && qRepr<40) goto lab2; };
    BS_ZFil( AL, ADat ); return AL; lab2:
    BS_ZFil( AL, ADat ); BS_VFil( qi1-qi0, ADat+qi0, 1 ); return AL; } 
  void iRDEmpConso( int AWp, int AHp, const tBS_Val *ARougSig, const tBS_Val *ARougRep, tBS_Val **ARes ) {
    // Logical consolidation by the rough_combined: specific collocation of rough_reproj and rough_signed
    // constraints the begin and the end of a cross-chord within lung-area.
    // Out MaskCons should be bright in lungs; but dark in backbone,peripheral.
    const int cNp = AWp*AHp; // Masks#0: rough_signed, #1:rough_reproj, #2:rough_combined
    tBS_Val  *qRog = ARes[0]; BS_Cpy( cNp, ARougSig, qRog ); // signed metric
    const int cRogInpMx = (AWp+AHp)>>2; // to assign cCutRep.rcuBlack,White
    // cCutRep.rcuGaO should be cRougMg*2 (hence Rog^ is ready to combine)
    tBS2_Cut  cCutRog = {{-cRogInpMx,cRogInpMx,cRougMg*2},{0,0}}; BS2_Cut( &cCutRog, qRog, cNp );
    tBS_Val  *qRep = ARes[1]; BS_Cpy( cNp, ARougRep, qRep ); // reprojected metric
    #if 0 // use histogram's cut
      tBS2_HistCut cCutRepr = {{0,1024,cReprMg,__min(AWp,AHp)/9},{0.44f,0.01f,0,0}}; BS2_HistCut( &cCutRep, qqRep, NULL, ARes+2, AWp, AHp );
    #else // cCutRep.rcuWhite has to be ~= above average rough_reproj within object
      tBS2_Cut cCutRep = {{0,777,cReprMg},{0,0}}; BS2_Cut( &cCutRep, qRep, cNp );
    #endif
    tBS_Val **qTs = ARes+2;
    tBS_Val  *qBuf = qTs[0]; for (int qi2=0;qi2<cNp;qi2++) // combined metric (input of consolidation)
      qBuf[qi2] = ii_RDEmpConsoComb( qRog[qi2], qRep[qi2] );
    const int cDirsHrz[] = {8,1,2}, cNDDh = mcASz(cDirsHrz),
              cDirsVrt[] = {4,5,6}, cNDDv = mcASz(cDirsVrt), cNDD = cNDDh+cNDDv;
    BS_DataDup( cNDD-1, cNp, qTs[0], qTs+1 ); // consolidated directional masks (the 'CDMs')
    BS_CPLA( iicpl_RDEmpConsoHrz, cNDDh, cDirsHrz, AWp, AHp, qTs ); // horizontal CDMs
    BS_CPLA( iicpl_RDEmpConsoVrt, cNDDv, cDirsVrt, AWp, AHp, qTs+cNDDh ); // vertical CDMs
    for (int qi4=1;qi4<cNDD;qi4++) BS_Add( cNp, qTs[qi4], qBuf ); // sum CDMs
    if (1) for (int qi6=0;qi6<cNp;qi6++) { int qV = ARougSig[qi6];
      if (qV<11) qBuf[qi6] = 0; else
      if (qV>cRougMg*2-11 && qBuf[qi6]>0) qBuf[qi6] = cNDD; };
  } // end of iRDEmpConso

  // = separate lungs "SplitLftRgtLungs":
  // the output of SplitLftRgtLungs (in Pyr downscaled)
  int      rRChLft[4], rRChRgt[4]; // cutting_half_rect Lft,Rgt
  tBS_Val *rDhLft[cLftRgtLungsAll], *rDhRgt[cLftRgtLungsAll]; // info in leading #LftRgtLungsOut entries
  tLungStt rSttLft, rSttRgt; // the stats

  // Receives frames 'Dp[0..LftRgtLungsInp-1] of Wp*Hp'. SLRLCutRects determines 'cutting_half_rect Lft,Rgt'.
  // SLRLCpyBuf splits 'Dp#i' to 'DhLft#i', 'DhRgt#i'. Then it calls ProcessLftRgtLung twice and combines
  // #LftRgtLungsOut returned entries of 'DhLft','DhRgt' back to 'Dp'.
  void jSplitLftRgtLungs( int AWp, int AHp, tBS_Val **ADp, tBS_Val **AAp ) { 
    char qS[256]; // adhoc string-buf
    int  cRChLft[4], cRChRgt[4]; iSLRLCutRects( AWp, AHp, ADp[2], cRChLft, cRChRgt ); // cutting_half_rect Lft,Rgt
    // allocate buffers for ProcessLftRgtLung
    int cWhLft = cRChLft[2]-cRChLft[0], cHhLft = cRChLft[3]-cRChLft[1], cNhLft = cWhLft*cHhLft, 
        cWhRgt = cRChRgt[2]-cRChRgt[0], cHhRgt = cRChRgt[3]-cRChRgt[1], cNhRgt = cWhRgt*cHhRgt;
    tBS_Val *qDhLft[cLftRgtLungsAll], *qDhRgt[cLftRgtLungsAll],
            *qAhLft[cLftRgtLungsAve], *qAhRgt[cLftRgtLungsAve];
    BS_DatasAlloc( mcASz(qDhLft), cNhLft, qDhLft ); BS_DatasAlloc( mcASz(qDhRgt), cNhRgt, qDhRgt ); 
    BS_DatasAlloc( mcASz(qAhLft), cNhLft, qAhLft ); BS_DatasAlloc( mcASz(qAhRgt), cNhRgt, qAhRgt );
    // split input entries Dp[0..LftRgtLungsInp-1] in continuous order
    for (int qi1=0;qi1<cLftRgtLungsInp;qi1++) { // pass Dp[0..LftRgtLungsInp-1]
      iSLRLCpyBuf( AWp, AHp, ADp[qi1], cRChLft, qDhLft[qi1],false, false );
      iSLRLCpyBuf( AWp, AHp, ADp[qi1], cRChRgt, qDhRgt[qi1], true, false ); };
    for (int qi2=0;qi2<cLftRgtLungsAve;qi2++) { // pass Ap[0..cLftRgtLungsAve-1]
      const int cFlipReIdx[8] = {0,7,6,5,4,3,2,1}; // assume cBS_AllDirsT from UT0.h
      iSLRLCpyBuf( AWp, AHp, AAp[qi2], cRChLft, qAhLft[cFlipReIdx[qi2]],false, false );
      iSLRLCpyBuf( AWp, AHp, AAp[qi2], cRChRgt, qAhRgt[qi2], true, false ); };
    BS_DatasRelease( cLftRgtLungsAve, AAp ); // note that these were allocated
    BS_DatasRelease( cLftRgtLungsAll, ADp ); // by the caller RibsDetection
    // call "ProcessLftRgtLung" 
    tLRLWr qLRLWr[2];
    qLRLWr[0].Init( false, cWhLft, cHhLft, qDhLft, qAhLft, rDbgDir!=NULL );
    qLRLWr[1].Init( true,  cWhRgt, cHhRgt, qDhRgt, qAhRgt, rDbgDir!=NULL );
    BS_SpawnThreads( 0, 0, 2, ppr_ProcessLftRgtLung, qLRLWr ); // !!!!
    memcpy( rRChLft, cRChLft, sizeof(rRChLft) ); memcpy( rRChRgt, cRChRgt, sizeof(rRChRgt) );
    memcpy( rDhLft, qDhLft, sizeof(rDhLft) ); memcpy( rDhRgt, qDhRgt, sizeof(rDhRgt) );
    memcpy( &rSttLft, &qLRLWr[0].rStt, sizeof(tLungStt) ); memcpy( &rSttRgt, &qLRLWr[1].rStt, sizeof(tLungStt) );
    sprintf( qS, cTPS" inf RibCage SplitLftRgtLungs done %dx%d %d:%d||%d:%d %dx%d",
      cWhLft, cHhLft, cRChLft[2], cRChLft[1], cRChRgt[0], cRChRgt[1], cWhRgt, cHhRgt ); WriteLn( qS );
    BS_DatasRelease( mcASz(qAhLft), qAhLft ); BS_DatasRelease( mcASz(qAhRgt), qAhRgt ); } // of jSplitLftRgtLungs

  // = for SplitLftRgtLungs: lft&rgt cutting rectangles
  static int iicpl_SLRLPrj( int AL, tBS_Val *ADat ) { unsigned qV = mcSum( AL,ADat )>>7; ADat[0] = __min(qV,0x7FFF); return 1; }
  static void iSLRLPrj( int AWp, int AHp, const tBS_Val *AMaskCons, tBS_Val *APrj, bool AByCols ) {
    const int cNp = AWp*AHp; tBS_Val *qT = new tBS_Val[cNp]; BS_Cpy( cNp, AMaskCons, qT );
    if (AByCols) { BS_CPLCCols( iicpl_SLRLPrj, AWp, AHp, qT ); BS_Cpy( AWp, qT, APrj );  }
            else { BS_CPLCRows( iicpl_SLRLPrj, AWp, AHp, qT, true ); mcGath<tBS_Val>( AHp, qT, APrj, AWp ); }; delete qT; }
  static void iSLRLPrjCrimpMgn( int AN, const tBS_Val *APrj, unsigned AZCrimp, int *AMgn, int AGap0, int AGap1 ) { // vrt or hrz limits of both lungs
    // treat the projection as 'as a histogram'; get and crimp horizontal or vertical 'unsigned margings'
    int qMgn0 = mcHistNP<tBS_Val>(AZCrimp,APrj)-AGap0, qMgn1 = -mcHistNP<tBS_Val>(-(int)(AZCrimp),APrj+AN-1)-AGap1;
    AMgn[0] = __max(qMgn0,0); AMgn[1] = __min(qMgn1,AN/2); }
  static void iSLRLInnerHrzIvl( int AWp, tBS_Val *APrj, int *AIvl ) { // horz limits of backbone's interval in the center
    int qW = AWp>>1; // the output: RCutLft[rgt]:=Ivl[0], RCutRgt[lft]:=Ivl[1].
    unsigned qZCrimp0 = (unsigned)(mcSum( qW, APrj )*0.005f), qZCrimp1 = (unsigned)(mcSum( AWp-qW, APrj+qW )*0.005f);
    int qBnd0 = mcHistNP<tBS_Val>(-(int)(qZCrimp0),APrj+qW-1)+14, qBnd1 = mcHistNP<tBS_Val>(qZCrimp1,APrj+qW)-14;
    if (qBnd0>qBnd1) qBnd0 = qBnd1 = (qBnd0+qBnd1)>>1;
    AIvl[0] = qW+qBnd0; AIvl[1] = qW+qBnd1; }
  void iSLRLCutRects( int AWp, int AHp, const tBS_Val *AMaskCons, int *ARChLft, int *ARChRgt ) { 
    // Inp MaskCons<olidated> should be bright in lungs; but dark in backbone,peripheral.
    // treat the projection as 'as a histogram'; get and crimp lft,rgt- or top,btm-bounds
    tBS_Val *qPrj = new tBS_Val[__max(AWp,AHp)];
    iSLRLPrj( AWp, AHp, AMaskCons, qPrj, true );
    int qZCrimp = ((unsigned)(mcSum(AWp,qPrj)))>>8; if (qZCrimp<9) qZCrimp = 9; //
    int qMgnHrz[2]; iSLRLPrjCrimpMgn( AWp, qPrj, qZCrimp, qMgnHrz, 16, 16 );
    iSLRLPrj( AWp, AHp, AMaskCons, qPrj, false );
    int qMgnVrt[2]; iSLRLPrjCrimpMgn( AHp, qPrj, qZCrimp, qMgnVrt, 12, 38 );
    // get projection within inner_height to determine backbone location in the center
    int qTopInn = qMgnVrt[0], qBtmInn = AHp-qMgnVrt[1], qHInn = qBtmInn-qTopInn;
    qTopInn += (int)(qHInn*0.1f); qBtmInn -= (int)(qHInn*0.1f);
    iSLRLPrj( AWp, qBtmInn-qTopInn, AMaskCons+qTopInn*AWp, qPrj, true );
    int qIvlInn[2]; iSLRLInnerHrzIvl( AWp, qPrj, qIvlInn );
    delete qPrj;
    // assign cutting_half_rects
    ARChLft[1] = ARChRgt[1] = qMgnVrt[0]; 
    ARChLft[3] = ARChRgt[3] = qMgnVrt[0]+((AHp-qMgnVrt[1]-qMgnVrt[0])&~3);
    ARChLft[0] = qMgnHrz[0]; ARChRgt[2] = AWp-qMgnHrz[1];    
    ARChLft[2] = qIvlInn[0]; ARChRgt[0] = qIvlInn[1];
    if (1) { // synchronize widths
      int qWl = ARChLft[2]-ARChLft[0], qWr = ARChRgt[2]-ARChRgt[0], qWm = __max(qWl,qWr);
      const int cWm[] = {ARChLft[2],AWp-ARChRgt[0]};
      for (int qi=0;qi<mcASz(cWm);qi++) if (qWm>=cWm[qi]) qWm = cWm[qi]-1;
      qWm &= ~3; // optional alignment
      ARChLft[0] = ARChLft[2]-qWm; ARChRgt[2] = ARChRgt[0]+qWm; };
    if (0) { // biggest rects 
      const int  cCenX = AWp>>1, cWh = cCenX, cHh = AHp*3/4, cNh = cWh*cHh;
      ARChLft[0] = cCenX-cWh+1; ARChLft[1] = 0; ARChLft[2] = cCenX; ARChLft[3] = 0+cHh;
      ARChRgt[0] = cCenX; ARChRgt[1] = 0; ARChRgt[2] = cCenX+cWh; ARChRgt[3] = 0+cHh; };
    if (1 && rDbgSfx) { // debug temporary 
      DbgDirSet( "c", AWp, AHp ); 
      const int cNp = AWp*AHp; tBS_Val *qT = new tBS_Val[cNp]; BS_Cpy( cNp, AMaskCons, qT );
      BS_DrawRect( AWp, AHp, qT, ARChLft, 9999 ); 
      BS_DrawRect( AWp, AHp, qT, ARChLft, 9999 ); BS_DrawRect( AWp, AHp, qT, ARChRgt, 9998 );
      char qS[256]; sprintf( qS, "%04d,%04d %04dx%04d\n%04d,%04d %04dx%04d",
        ARChLft[0], ARChLft[1], ARChLft[2]-ARChLft[0], ARChLft[3]-ARChLft[1],
        ARChRgt[0], ARChRgt[1], ARChRgt[2]-ARChRgt[0], ARChRgt[3]-ARChRgt[1] );
      const int qRInn[4] = {qIvlInn[0],qTopInn,qIvlInn[1],qBtmInn}; BS_DrawRect( AWp, AHp, qT, qRInn, -1 ); 
      GIni_Lab8x8( AWp, AHp, qT, qS, 9997 );
      strcpy( rDbgSfx, "ICutRects.raw" ); BS_Write( rDbgPfx, cNp, qT  );
      delete qT; };
  } // end of iSLRLCutRects

  // = for SplitLftRgtLungs: 'Dp[i]<=>Dh[i]' within 'RCh'
  static void iSLRLCpyBuf( int AWp, int AHp, tBS_Val *ADp, const int *ARCh, tBS_Val *ADh,
    bool ARgt, bool ACpyBck ) { const int cWh = ARCh[2]-ARCh[0], cHh = ARCh[3]-ARCh[1];
    if (ACpyBck) { if (!ARgt) mcFlipHrz( cWh, cHh, ADh, cWh ); BS_RowsCpy( cWh, ADh, AWp, ADp+ARCh[1]*AWp+ARCh[0], cWh, cHh ); }
            else { BS_RowsCpy( AWp, ADp+ARCh[1]*AWp+ARCh[0], cWh, ADh, cWh, cHh ); if (!ARgt) mcFlipHrz( cWh, cHh, ADh, cWh ); }; }

  // = for SplitLftRgtLungs: invoke tLRLWr.ProcessLftRgtLung
  static void ppr_ProcessLftRgtLung( tBS_SpawnThreadRec *AP ) {
    for (int qi=AP->rFst;qi<AP->rLst;qi++) ((tLRLWr*)AP->rP[0])[qi].ProcessLftRgtLung(); }

  // = tune full-frame ROIs "TuneLftRgtFROIs":
  // Forms input of ProcessLftRgtFROI as ~"Org and output of SplitLftRgtLungs".
  // Then calls ProcessLftRgtFROI twice and overwrites tuned FROIs in Org.
  void jTuneLftRgtFROIs() { // creates and executes two FROIWrappers
    char qS[256]; // adhoc string-buf
    const int cW = rW, cH = rH, cN = cW*cH, cWp = cW>>cPyr, cHp = cH>>cPyr, cNp = cWp*cHp;
    // "SplitLftRgtLungs" outputs RChLft[4], RChRgt[4],
    // DhLft[cLftRgtLungsAll], DhRgt[cLftRgtLungsAll] (leading LftRgtLungsOut entries contain info)
    int cRChLft[4]; iSLRLOrgRect( rRChLft, cRChLft );
    int cRChRgt[4]; iSLRLOrgRect( rRChRgt, cRChRgt );
    // allocate buffers for ProcessLftRgtFROI
    int pWhLft = rRChLft[2]-rRChLft[0], pHhLft = rRChLft[3]-rRChLft[1], pNhLft = pWhLft*pHhLft, 
        pWhRgt = rRChRgt[2]-rRChRgt[0], pHhRgt = rRChRgt[3]-rRChRgt[1], pNhRgt = pWhRgt*pHhRgt;
    int cWhLft = cRChLft[2]-cRChLft[0], cHhLft = cRChLft[3]-cRChLft[1], cNhLft = cWhLft*cHhLft, 
        cWhRgt = cRChRgt[2]-cRChRgt[0], cHhRgt = cRChRgt[3]-cRChRgt[1], cNhRgt = cWhRgt*cHhRgt;
    tBS_Val *qDhLft[cLftRgtFROI]; BS_DatasAlloc( mcASz(qDhLft), cNhLft, qDhLft ); 
    tBS_Val *qDhRgt[cLftRgtFROI]; BS_DatasAlloc( mcASz(qDhRgt), cNhRgt, qDhRgt ); 
    // pass Org entry<ies> from cInps
    iSLRLCpyBuf( cW, cH, cInps[0], cRChLft, qDhLft[0], false, false );
    iSLRLCpyBuf( cW, cH, cInps[0], cRChRgt, qDhRgt[0], true,  false );
    // pass output of SplitLftRgtLungs (point to change)
    for (int qi=0;qi<cLftRgtLungsOut;qi++) { 
      BS_Cpy( pNhLft, rDhLft[qi], qDhLft[qi+1] );
      BS_Cpy( cNhLft, BS_ScaleUpA(cPyr,pWhLft,pHhLft,qDhLft[qi+1],rTs[0],true), qDhLft[qi+1] );
      BS_Cpy( pNhRgt, rDhRgt[qi], qDhRgt[qi+1] );
      BS_Cpy( cNhRgt, BS_ScaleUpA(cPyr,pWhRgt,pHhRgt,qDhRgt[qi+1],rTs[0],true), qDhRgt[qi+1] ); };
    // call "ProcessLftRgtFROI" 
    tLRFROIWr qLRFROIWr[2]; memset( qLRFROIWr, 0, sizeof(qLRFROIWr) );
    qLRFROIWr[0].Init( false, cWhLft, cHhLft, qDhLft, &rSttLft, rDbgDir!=NULL );
    qLRFROIWr[1].Init( true,  cWhRgt, cHhRgt, qDhRgt, &rSttRgt, rDbgDir!=NULL );
    BS_SpawnThreads( 0, 0, 2, ppr_TuneLftRgtROI, qLRFROIWr ); // !!!!
    // pass tuned FROIs to cInps
    iSLRLCpyBuf( cW, cH, cInps[0], cRChLft, qDhLft[0], false, true );
    iSLRLCpyBuf( cW, cH, cInps[0], cRChRgt, qDhRgt[0], true,  true );
    // leave
    sprintf( qS, cTPS" inf RibCage TuneLftRgtFROIs done %dx%d %d:%d||%d:%d %dx%d",
      cWhLft, cHhLft, cRChLft[2], cRChLft[1], cRChRgt[0], cRChRgt[1], cWhRgt, cHhRgt ); WriteLn( qS );
    BS_DatasRelease( mcASz(qDhLft), qDhLft ); BS_DatasRelease( mcASz(qDhRgt), qDhRgt ); } // of jTuneLftRgtFROIs

  // = for TuneLftRgtFROIs: cutting_half_rect Lft,Rgt in org frame
  void iSLRLOrgRect( const int *AP, int *AO ) { for (int qi=0;qi<4;qi++) AO[qi] = AP[qi]<<cPyr; }

  // = for TuneLftRgtFROIs: invoke tLRFROIWr.ProcessLftRgtFROI
  static void ppr_TuneLftRgtROI( tBS_SpawnThreadRec *AP ) {
    for (int qi=AP->rFst;qi<AP->rLst;qi++) ((tLRFROIWr*)AP->rP[0])[qi].TuneLftRgtFROI(); }

}; // end of RiWr

//// exports ////  
void BS7b_RibsDetection( int AW, int AH, tBS_Val *AFrm, // -- ribs-detection preprocess
  tBS_Val **ATs, const char *ADbgPfx ) { const int cN = AW*AH;
  tRiWr qRiWr; memset( &qRiWr, 0, sizeof(tRiWr) );
  qRiWr.rW = AW; qRiWr.rH = AH; qRiWr.cInps[0] = AFrm; qRiWr.rTs = ATs;
  if (ADbgPfx!=NULL) qRiWr.DbgSet();
  qRiWr.RibsDetection(); 
  if (1 && ADbgPfx) {
    //tBS_Val *qTs[] = {AFrm,qRiWr.rTs[0],qRiWr.rTs[1],qRiWr.rTs[2]};
    tBS_Val *qTs[] = {AFrm};
    BS_DatasWrt( ADbgPfx, cN, qTs, 0, mcASz(qTs) ); };
}