// I16_Ut4.cpp : the routines with (~mostly) CPU|CUDA as parameter
// the routines here use fixed frames; CPU:BS2_Dim=>BS2_W,BS2_H, CUDA:BS3_Dim=>BS3_W,BS3_H
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut3.h"
#include "BS_Ut4.h"

//// class "CLAHED as self" incapsulates control-parameters and work-buffers ////
class tBS4_CLAHED { protected: // - intrinsics
  static const int
    cStt     = 256*256, // statistics's blocks
    cXLATLg  = cBS4_XLATGamLg, // XLATs' gamut binary logarithm
    cXLATGam = 1<<cXLATLg, cXLATGa1 = cXLATGam-1, cXLATGa2 = cXLATGa1-1;
  char      rDbgPfx[256], *rDbgSfx; // debug-dump's strings,
  int       rNx, rNy, rSz,  // #blocks by x&y, #pixels_in_a_block
            rBlr; // XLAT's blur radius (when XLATSht=8)
  float     rSLCurvExtMn, rSLCurvExtMx;
  bool      rRef;
  unsigned *rHists; // histograms interpret first as 'unsigned' then as 'int'
  unsigned *rHistC; // counts (used by ~TraAveRef; unused by main 'Transform')
  tBS_Val   rSttMn[cStt], rSttGa[cStt], // S-like curve inp:Mn0,Mx0 (by "iHistStt")
            rSttMn1[cStt], rSttMx1[cStt], // S-like curve out:Mn1,Mx1  
            rStt[cStt]; // adhoc temp buffer Nx*Ny for a statistic (mainly for debug-print)
  int       rSttNP[cStt]; // sum of pixels' metric

  static int iiRadBlr( int ABlr ) { return ABlr>0 ? ABlr<<(cXLATLg-8) : ABlr; } // normalize blur against 'XLATLg=8'

  // helpers for debug-print tiles Nx*Ny (using 'adhoc temp buffer')
  inline int iiDbgTileX( int ATgtCol ) { return (rNx+2)*ATgtCol+2; }
  inline int iiDbgTileY( int ATgtRow ) { return (rNy+2)*ATgtRow+2; }
  static const int cDbgTileCols = 11, cDbgTileRows = 9; // # tile-rulers by X & Y
  void iiDbgPrintRulers( tBS_Val *ATgt ) { // org frame's col- & row-indexes
    const int cW = BS2_W(), cBx = cW/rNx, cDx = (cW-cBx*rNx)>>1;
    tBS_Val *qTgtCol0 = ATgt+iiDbgTileX(0); for (int qx=0;qx<rNx;qx++) qTgtCol0[qx] = cBx*qx+cDx;
    for (int qj=0;qj<cDbgTileCols;qj++) BS_Cpy( rNx, qTgtCol0, ATgt+iiDbgTileX(qj) );
    const int cH = BS2_H(), cBy = cH/rNy, cDy = (cH-cBy*rNy)>>1;
    tBS_Val *qTgtRow0 = ATgt+iiDbgTileY(0)*cW; for (int qy=0;qy<rNy;qy++) qTgtRow0[qy*cW] = cBy*qy+cDy;
    for (int qj=0;qj<cDbgTileRows;qj++) BS_RowsCpy( cW, qTgtRow0, cW, ATgt+iiDbgTileY(qj)*cW, 1, rNy ); }
  void iiDbgPrintTile( const tBS_Val *ASrc, int ATgtCol, int ATgtRow, tBS_Val *ATgt ) { // copy a tile
    const int cW = BS2_W(); BS_RowsCpy( rNx, ASrc, cW, ATgt+iiDbgTileY(ATgtRow)*cW+iiDbgTileX(ATgtCol), rNx, rNy ); }
  void iDbgPrintHists( int ATgtRow, tBS_Val *ATgt ) { // print contents of 'Hists'
    const int cW = BS2_W(), cNxy = rNx*rNy; tBS_Val *qStt = rStt; // via 'adhoc temp buffer'
    int qjd = cW/rNx-1; qjd = __min(qjd,9); for (int qj=0;qj<=qjd;qj++) { int *qHists = (int*)rHists+qj*cXLATGa1/qjd;
    for (int qi=0;qi<cNxy;qi++) qStt[qi] = mcRng( qHists[qi<<cXLATLg], -0x7FFF, 0x7FFF );
    iiDbgPrintTile( qStt, qj, ATgtRow, ATgt ); }; }
  void iiDbgPrintSttGa( int ATgtCol, int ATgtRow, tBS_Val *ATgt, bool ATag ) { // print SttGa with indication
    const int cW = BS2_W(), cNxy = rNx*rNy; tBS_Val *qStt = rStt; // indicate small or big values as ~gradient
    for (int qi=0;qi<cNxy;qi++) { const int cSttGaDbgThr = 11; // indicate '<','>' vs. adhoc threshold
    int qV = rSttGa[qi]; qStt[qi] = (qV<cSttGaDbgThr)^ATag?qi:qV; };
    iiDbgPrintTile( qStt, ATgtCol, ATgtRow, ATgt ); }

  void iHistStt( int An ) { // get statistics of single histogram (block #n)
    unsigned *qHist = rHists+(An<<cXLATLg);
    const int cSz = rRef ? mcSum<unsigned,cXLATGam>(qHist) : rSz, // if there was a reference
              cNP = cSz-qHist[cXLATGa1], // sum of metric_pixels without reserved entry #XLATGa1
    #if 1
              cNPmn = cNP>>7, // to ~rough on the btm
              cNPmx = cNP>>7; // to ~rough on the top
    int       qMn0 = mcHistNP<unsigned>( cNPmn, qHist ),
              qMx0 = cXLATGa1+mcHistNP<unsigned>( -cNPmx, qHist+cXLATGa2 );
    #else // !!!!
              cNPmn = 0, cNPmx = 0, qMn0 = 0, qMx0 = cXLATGa1;
    #endif          
    rSttNP[An] = cNP-cNPmn-cNPmx; rSttMn[An] = qMn0; rSttGa[An] = qMx0-qMn0; } // approx objective gamut is ~Mn0..Mx0~
  static void ppr_CLAHEStt( tBS_SpawnThreadRec *AP ) { // invoke HistStt...
    tBS4_CLAHED *qSelf = (tBS4_CLAHED*)AP->rP[0]; int cLst = AP->rLst; for (int qn=AP->rFst;qn<cLst;qn++) qSelf->iHistStt( qn ); }
  
  void iHistSttSLkCrv( int An ) { // contrast_limitation of single histogram by S-like curve (block #n)
    int qMn0 = rSttMn[An], qMx0 = qMn0+rSttGa[An], qGa0 = qMx0-qMn0+1; // inp:Mn0,Mx0
    const int cextLim = (qGa0<<3)+cXLATGa1/3 ; // limit of gamut's extension
    int qextMn = int(qMn0*rSLCurvExtMn); if (qextMn>cextLim) qextMn = cextLim;
    int qMn1 = qMn0-qextMn; if (qMn1<0) qMn1 = 0; // decrement Mn0
    int qextMx = (int)((cXLATGa1-qMx0)*rSLCurvExtMx); if (qextMx>cextLim ) qextMx = cextLim; // increment Mx0
    int qMx1 = qMx0+qextMx; if (qMx1>cXLATGa1) qMx1 = cXLATGa1;
    rSttMn1[An] = qMn1; rSttMx1[An] = qMx1; } // out:Mn1,Mx1

  void iiHistEquFill1to1( int An, int ADbgMrk ) { unsigned *qHist = rHists+(An<<cXLATLg);
    for (int qi=0;qi<cXLATGam;qi++) qHist[qi] = qi; rStt[An] = ADbgMrk; }
  void iHistEqu( int An ) { // equalize single histogram (block #n)
    // this procedure operates rather 'unsigneds' (e.g. for ~auto-subst 'div=>shr' when used)
    const unsigned cNP = rSttNP[An]; // '#pixes within Mn0,Mx0' should be ~> 'BlockSize/4'
    if (cNP<(unsigned)rSz>>2) { iiHistEquFill1to1( An, -1 ); return; }; // initial check '-1'
    unsigned *qHist = rHists+(An<<cXLATLg); // non-histogram pixels '=0' in BeamHardMap, '#255' in Hist
    const int cNPbh = qHist[cXLATGa1], cNPbhmx = (int)(cNP*0.33f); // max allowed #NPbh
    if (1 && cNPbh>=cNPbhmx) { iiHistEquFill1to1( An, -2 ); return; };  ; // initial check '-2'
    unsigned  qMn0 = rSttMn[An], qMx0 = qMn0+rSttGa[An], // S-like curve; inp:Mn0,Mx0,out:Mn1,Mx1 
              qMn1 = rSttMn1[An], qMx1 = rSttMx1[An]; // signed vars to unsigned
    if (rDbgSfx) { rSttMn[An] = qMn1; rStt[An] = qMx1; }; // for debug's dump
    const unsigned c2 = qMn0; // 'i2' fills [0.. Mn0[
    for (unsigned qi2=0;qi2<c2;qi2++) qHist[qi2] = mcIpl(qi2,c2,0,qMn1);
    const unsigned c4 = qMx0-qMn0; // 'i4' fills [Mn0..Mx0[
    #if 0 // linear mid-segment of S-like-curve
      for (unsigned qi4=0;qi4<c4;qi4++) qHist[qMn0+qi4] = mcIpl(qi4,c4,qMn1,qMx1 );
    #else // equalized mid-segment of S-like-curve
      const unsigned cD1 = qMx1-qMn1-1, cNPa = (cNP>>1)-1; tLLU qU = 0; 
      for (unsigned qi4=0;qi4<c4;qi4++) { unsigned qV = qHist[qi4+qMn0];
        qHist[qi4+qMn0] = (unsigned)((qU*cD1+cNPa)/cNP)+qMn1; qU += qV; };
    #endif
    const unsigned c6 = cXLATGam-qMx0; // 'i6' fills [Mx0..Gamut[
    for (unsigned qi6=0;qi6<c6;qi6++) qHist[qMx0+qi6] = mcIpl(qi6,c6,qMx1,cXLATGam);
    if (1 && cNPbh>cNPbhmx>>7) // ints here (as 'i-Hist[i]' is signed)
      for (int qi=0;qi<cXLATGam;qi++) qHist[qi] = mcIpl(cNPbh,cNPbhmx,(int)qHist[qi],qi); }
  static void ppr_CLAHEDEqu( tBS_SpawnThreadRec *AP ) { // invoke iHistEqu...
    tBS4_CLAHED *qSelf = (tBS4_CLAHED*)AP->rP[0]; int cLst = AP->rLst;
    for (int qn=AP->rFst;qn<cLst;qn++) qSelf->iHistEqu( qn ); }

  static void ppr_CLAHEDCrB( tBS_SpawnThreadRec *AP ) { // cross-blur XLATs at the same gray-levels
    tBS4_CLAHED *qSelf = (tBS4_CLAHED*)AP->rP[0]; //const tBS_Val *cSttGa = qSelf->rSttGa;
    int cNx = qSelf->rNx, cNy = qSelf->rNy, cNxy = cNx*cNy, cLst = AP->rLst;
    tBS_Val *qBuf0 = new tBS_Val[cStt], *qBuf1 = new tBS_Val[cStt];
    for (int qi=AP->rFst;qi<cLst;qi++) {
      int *qHist = ((int*)qSelf->rHists)+qi;
      for (int qn2=0;qn2<cNxy;qn2++) qBuf0[qn2] = qHist[qn2<<cXLATLg]; // gather 'int to BS_Val'
      BS_Ave3x3( cNx, cNy, qBuf0, qBuf1 ); // BS_Med1( cNx, cNy, qBuf0, qBuf1 );
      for (int qn4=0;qn4<cNxy;qn4++) qHist[qn4<<cXLATLg] = qBuf1[qn4]; // scatter 'BS_Val to int'
    }; // for i
    delete qBuf0; delete qBuf1; }

  void iHistBlr( int An, int *ATmp ) { // blur single histogram (block #n)
    int qRad = rBlr;
    if (qRad>0) { if (qRad>cXLATGam>>1) qRad = cXLATGam>>1; }
    else switch (qRad) { // blur wide hists more, blur narrow hists less
      case -1: { // for skull ~ necessary to avoid over-equalize in ~brain
        const int cGaCe = (int)(cXLATGam*20.0f/256); // the most object-like Ga0
        int       qGaDi = cGaCe-abs(cGaCe-rSttGa[An]); if (qGaDi<0) qGaDi = 0;
        qRad = mcIpl( qGaDi, cGaCe, (int)(22.0f/256*cXLATGam), (int)(12.0f/256*cXLATGam) );
      }; break;
      default: if (!An) ErrMsg( "CLAHEHistBlr wrong RadBlr" ); return; }; // switch
    int *qHist = (int*)rHists+(An<<cXLATLg);
    mcMovAve<int>( cXLATGam, qRad, (int*)qHist, ATmp ); memcpy( qHist, ATmp, cXLATGam*sizeof(int) ); }

  static void ppr_CLAHEDBlr( tBS_SpawnThreadRec *AP ) { // invoke HistBlr...
    tBS4_CLAHED *qSelf = (tBS4_CLAHED*)AP->rP[0]; int *qTmp = new int[qSelf->rSz], cLst = AP->rLst;
    for (int qn=AP->rFst;qn<cLst;qn++) qSelf->iHistBlr( qn, qTmp ); delete qTmp; }

  void iHistRefAve( int An ) { // average single refernce-histogram by (block #n)
    // this procedure operates rather 'unsigneds' (e.g. for ~auto-subst 'div=>shr' when used)
    const unsigned cNP = rSttNP[An]; // '#pixes within Mn0,Mx0' should be ~> 'BlockSize/4'
    int *qHist = (int*)rHists+(An<<cXLATLg), *qHisC = (int*)rHistC+(An<<cXLATLg);
    for (int qi2=0;qi2<cXLATGam;qi2++) { int qC = qHisC[qi2]; if (qC) qHist[qi2] /= qC; }; } 
    // for (int qi2=0;qi2<cXLATGam;qi2++) { int qC = qHisC[qi2]; if (qC) { qC = qHist[qi2]/qC; qHist[qi2] = mcRg(qC); }; }; }

  static void ppr_CLAHEDRefAve( tBS_SpawnThreadRec *AP ) { // invoke iHistRefAve...
    tBS4_CLAHED *qSelf = (tBS4_CLAHED*)AP->rP[0]; int cLst = AP->rLst;
    for (int qn=AP->rFst;qn<cLst;qn++) qSelf->iHistRefAve( qn ); }

  public: // - interface

  void DefineExts( tBS4_CLAHE_Exts *AExts ) { // beyond basic BS_GamCLAHE,BC_GamCLAHE
    memset( this, 0, sizeof(tBS4_CLAHED) ); if (!AExts) return;
    rSLCurvExtMn = AExts->rSLCurvExtMn; rSLCurvExtMx = AExts->rSLCurvExtMx; }

  void DbgSet( const char *ADbgPfx ) { if (ADbgPfx) rDbgSfx = rDbgPfx+sprintf( rDbgPfx, "%s_", ADbgPfx ); }

  int Transform( int AB, int AGam, int ABlr, const tBS_Val *AFrmOrg, // -- plain CLAHE
    const tBS_Val *ARef, tBS_Val **ADs, bool ACUDA ) { // out #0:CLAHE #1:~8bit_FrmOrg #2:CLAInc
    #ifdef _DEBUG
      const bool cDbgTrc = true;
    #else
      const bool cDbgTrc = rDbgSfx!=NULL;
    #endif
    const int  cW = ACUDA?BS3_W():BS2_W(), cH = ACUDA?BS3_H():BS2_H(), cN = cW*cH; char qS[256];
    rNx = cW>>AB; rNy = cH>>AB; rSz = 1<<(AB*2); rRef = ARef!=NULL; // == set fields
    if (ACUDA) BC_Cpy( AFrmOrg, ADs[0] ); else BS_Cpy( cN, AFrmOrg, ADs[0] ); 
    const int  cCuPxNu = 99; tBS2_HistCut cParCut = { {0,AGam,cXLATGa1,0},{cCuPxNu,cCuPxNu} };
    if (ACUDA) BS3_HistCut( &cParCut, ADs[0], NULL, ADs+1 );
          else BS2_HistCut( &cParCut, ADs[0], NULL, ADs+1 ); // 'full input' histogram
    int        qCatchMM[2]; memcpy( qCatchMM, xBS2_rhCatchMM, sizeof(qCatchMM) );
    if (cDbgTrc) { sprintf( qS, cTPS" CLA%s(%d,%d,%d) HistCut:=%d..%d=>%d %.32s", ACUDA?"CUDA":"CPU",
      AB, AGam, ABlr, qCatchMM[0], qCatchMM[1], cXLATGa1, rDbgPfx ); WriteLn( qS ); };
    // out#1:=the_org in CLAHE_internal_scale 
    if (ACUDA) BC_Cpy( ADs[0], ADs[1] ); else BS_Cpy( cN, ADs[0], ADs[1] );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "A.raw" ); // dbg_dump of the_org in CLAHE_internal_scale
      if (ACUDA) BC_Write( rDbgPfx, cN, ADs[1] ); else BS_Write( rDbgPfx, cN, ADs[1] ); };
    // == build and treat histograms
    if (ACUDA) rHists = BC_GamCLAHEHists( ADs[0], AB, cXLATGam, ARef );
          else rHists = BS_GamCLAHEHists( cW, cH, ADs[0], AB, cXLATGam, ARef );
    if (!rHists) { ErrMsg( "CLAHETra Hist=nil" ); return -1; }; // histogram of 'counts' or of 'refs_sum'
    // - local buffers
    const int  cNxy = rNx*rNy, cPrintHistRow = 3; tBS_Val *qStt = rStt; // local buffer of 'Nxy'
    tBS_Val   *qDbg = NULL; if (rDbgSfx) { // schematic_charts (tile_row #0)
      qDbg = new tBS_Val[cN]; BS_ZFil( cN, qDbg ); iiDbgPrintRulers( qDbg );
      iDbgPrintHists( cPrintHistRow+0, qDbg ); }; // 'if debug_dump'
    // - equalize histograms with contrast limitation
    BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEStt, this ); // - get hists' stats SttNP[],Mn[],Ga[]
    if (qDbg) { // schematic_charts (tile_row #0)
      iiDbgPrintTile( rSttMn, 0, 0, qDbg ); iiDbgPrintTile( rSttGa, 1, 0, qDbg );
      for (int qi=0;qi<cNxy;qi++) { int qV = rSttNP[qi]; if (rRef) qV /= 10;
        qStt[qi] = mcRng(qV,-0x7FFF,0x7FFF); }; iiDbgPrintTile( qStt, 2, 0, qDbg ); }; // 'if debug_buf'
    BS_Med( 1, rNx, rNy, rSttMn, qStt ); BS_Cpy( cNxy, qStt, rSttMn );
    BS_Med( 1, rNx, rNy, rSttGa, qStt ); BS_Cpy( cNxy, qStt, rSttGa );
    if (qDbg) { // schematic_charts (tile_row #1)
      iiDbgPrintTile( rSttMn, 0, 1, qDbg ); iiDbgPrintTile( rSttGa, 1, 1, qDbg );
      iiDbgPrintSttGa( 2, 1, qDbg, false ); }; // // 'if debug_buf'
    for (int qn=0;qn<cNxy;qn++) iHistSttSLkCrv( qn ); // - contrast_limitation by S-like curve
    BS_Ave3x3( rNx, rNy, rSttMn1, qStt ); BS_Cpy( cNxy, qStt, rSttMn1 );
    BS_Ave5x5( rNx, rNy, rSttMx1, qStt ); BS_Cpy( cNxy, qStt, rSttMx1 );
    BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDEqu, this ); // - equalize
    if (qDbg) { // schematic_charts (tile_row #2)
      iiDbgPrintTile( rSttMn1, 0, 2, qDbg ); iiDbgPrintTile( rSttMx1, 1, 2, qDbg );
      iiDbgPrintSttGa( 2, 2, qDbg, true ); iDbgPrintHists( cPrintHistRow+2, qDbg ); }; // 'if debug_buf'
    BS_SpawnThreads( 0, 0, cXLATGam, ppr_CLAHEDCrB, this ); if (qDbg) iDbgPrintHists( cPrintHistRow+3, qDbg ); // - crossblur
    rBlr = iiRadBlr( ABlr ); if (rBlr) { // - blur if Blr<>0
      BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDBlr, this ); if (qDbg) iDbgPrintHists( cPrintHistRow+4, qDbg ); };
    if (qDbg) { // dbg_dump schematic_charts
      #if 1 // (for convenient observation only) upscale schematic_charts 
        int qdW = iiDbgTileX( cDbgTileCols ), qdH = iiDbgTileY( cDbgTileRows );
        tBS_Val *qT = new tBS_Val[cN], *qBuf0 = qDbg, *qBuf1 = qT;
        BS_RowsCpy( cW, qBuf0, qdW, qBuf1, qdW, qdH );
        for (int qUpScale=0;qUpScale<3;qUpScale++) { if ((qdW>cW>>1)||(qdH>cH>>1)) break;
          mcSwap( tBS_Val*, qBuf0, qBuf1 ); BS_DimUp( qdW, qdH, qBuf0, qBuf1 ); qdW <<= 1; qdH <<= 1; }; // for UpScale
        BS_ZFil( cN, qBuf0 ); BS_RowsCpy( qdW, qBuf1, cW, qBuf0, qdW, qdH ); BS_Cpy( cN, qBuf0, qDbg ); delete qT;
      #endif // (for convenient observation only)
      strcpy( rDbgSfx, "C.raw" ); BS_Write( rDbgPfx, cN, qDbg ); }; // 'if debug_buf'
    // == apply CLAHE
    if (ACUDA) BC_GamCLAHEApply( ADs[0], AB, cXLATGam, (int*)rHists );
          else BS_GamCLAHEApply( cW, cH, ADs[0], AB, cXLATGam, (int*)rHists );
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "D.raw" ); // dbg_dump of out in CLAHE_internal_scale
      if (ACUDA) BC_Write( rDbgPfx, cN, ADs[0] ); else BS_Write( rDbgPfx, cN, ADs[0] ); }; // 'if debug_dump'
    // - get 'narrow output' histogram
    if (ACUDA) BC_GamHist( ADs[0], cW, rHists, cXLATGa1, NULL ); 
          else BS_GamHist( cW, cH, ADs[0], cW, rHists, cXLATGa1, NULL ); 
    int        qTransMM[2] = { mcHistNP<unsigned>(cCuPxNu,rHists), cXLATGa1-abs(mcHistNP<unsigned>(-cCuPxNu,rHists+cXLATGa2)) };
    if (ACUDA) { // - restore org histogram
      BC_VAdd( ADs[0], -qTransMM[0] );
      BC_VMul( ADs[0], (float)(qCatchMM[1]-qCatchMM[0]+1)/(qTransMM[1]-qTransMM[0]+1) );
      //if (cXLATLg!=8) BC_VMul( ADs[1], 1.0f/(1<<(cXLATLg-8)) );
    } else {
      BS_VAdd( cN, ADs[0], -qTransMM[0] );
      BS_VMul( cN, ADs[0], (float)(qCatchMM[1]-qCatchMM[0]+1)/(qTransMM[1]-qTransMM[0]+1) );
      //if (cXLATLg!=8) BS_VMul( cN, ADs[1], 1.0f/(1<<(cXLATLg-8)) );
    };
    //if (AbhMap) { tBS_Val *qT = ADs[0]; for (int qi=0;qi<cN;qi++) if (!cBS2_bhMap[qi]) qT[qi] += 1; };
    if (qDbg) delete qDbg; delete rHists; return 2; } // end of Transform

  static const tBS2_HistCut* iFrmNrmHistCutRec( int AGam ) { // -- normalized frame 
    static const int cCuPxNu = 11; static tBS2_HistCut cParCut = {{0,AGam,cXLATGam,0},{cCuPxNu,cCuPxNu}}; return &cParCut; }
  static void FrmNrm_CPU( int AGam, tBS_Val *AFrm, tBS_Val **ATs2, int AW=0, int AH=0 ) { 
    if (AGam>cXLATGam+1) BS2_HistCut( iFrmNrmHistCutRec(AGam), AFrm, NULL, ATs2, AW, AH ); }
  static void FrmNrm_CUDA( int AGam, tBS_Val *AFrm, tBS_Val **ATs2 ) { 
    if (AGam>cXLATGam+1) BS3_HistCut( iFrmNrmHistCutRec(AGam), AFrm, NULL, ATs2 ); }

  int TraAveRef( int AB, int ABlr, const tBS_Val *AFrmNrm, // -- ref_averages
    int ANRefs, tBS_Val **ARefs, tBS_Val **ADs, bool ACUDA ) {    
    const int  cW = BS2_W(), cH = BS2_H(), cN = cW*cH; char qS[256];
    rNx = cW>>AB; rNy = cH>>AB; rSz = 1<<(AB*2); rBlr = iiRadBlr( ABlr );  // == set fields
    if (ACUDA) BC_Cpy( AFrmNrm, ADs[ANRefs] );
          else BS_Cpy( cN, AFrmNrm, ADs[ANRefs] );
    // BS_Write( "qqA.raw", cN, ADs[ANRefs] );
    if (ACUDA) rHistC = BC_GamCLAHEHists( ADs[ANRefs], AB, cXLATGam, NULL );
          else rHistC = BS_GamCLAHEHists( cW, cH, ADs[ANRefs], AB, cXLATGam, NULL );
    if (!rHistC) { ErrMsg( "CLAHETraAveRef Hists=nil" ); return -1; }; // histogram of 'counts' or of 'refs_sum'
    const int  cNxy = rNx*rNy; for (int qi2=0;qi2<ANRefs;qi2++) { // loop by references
      if (ACUDA) BC_Cpy( ADs[ANRefs], ADs[qi2] );
            else BS_Cpy( cN, ADs[ANRefs], ADs[qi2] );
      const bool cHistCountsOnly = ARefs[qi2]==(tBS_Val*)-1;
      if (cHistCountsOnly) { const int cHistsSz = cNxy*cXLATGam; // == histogram counts only
        rHists = new unsigned[cHistsSz]; memcpy( rHists, rHistC, cHistsSz*sizeof(int) );
      } else { // == ref_averages
        if (ACUDA) rHists = BC_GamCLAHEHists( ADs[qi2], AB, cXLATGam, ARefs[qi2] ); 
              else rHists = BS_GamCLAHEHists( cW, cH, ADs[qi2], AB, cXLATGam, ARefs[qi2] ); };
      if (!rHists) { if (ACUDA) BC_ZFil( cN, ADs[qi2] ); else BS_ZFil( cN, ADs[qi2] );
        sprintf( qS, "CLAHETraAveRef Histz[%d]=nil", qi2 ); ErrMsg( qS ); continue; };
      #if 0 // zero histogram counts for "pix=0"
        BS_VScat( cNxy, rHists, cXLATGam, 0 );
      #endif
      if (!cHistCountsOnly) BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDRefAve, this ); // - invoke iHistRefAve...
      BS_SpawnThreads( 0, 0, cXLATGam, ppr_CLAHEDCrB, this ); //if (qDbg) iDbgPrintHists( cPrintHistRow+5, qDbg ); // - crossblur
      if (rBlr) { // - blur if Blr<>0
        BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDBlr, this ); }; // if (qDbg) iDbgPrintHists( cPrintHistRow+6, qDbg ); 
      if (ACUDA) BC_GamCLAHEApply( ADs[qi2], AB, cXLATGam, (int*)rHists );
            else BS_GamCLAHEApply( cW, cH, ADs[qi2], AB, cXLATGam, (int*)rHists ); // == apply CLAHE
      delete rHists; }; // for i2
    if (rHistC) delete rHistC; return ANRefs; } // end of TraAveRef

  int TraRefAir( int AB, int ABlr, const tBS_Val *AFrmNrm, // -- ref_average for voluntary dimension
    const tBS_Val *ARef, tBS_Val **ADs, int AW, int AH ) {
    if (!AW) { AW = BS2_W(); AH = BS2_H(); }; 
    rNx = AW>>AB; rNy = AH>>AB; rSz = 1<<(AB*2); rBlr = iiRadBlr( ABlr ); // == set fields
    const int cN = AW*AH, cNxy = rNx*rNy;
    BS_DataDup( 2, cN, AFrmNrm, ADs+0 ); // BS_Write( "qqA.raw", cN, ARef ); //???
    rHistC = BS_GamCLAHEHists( AW, AH, ADs[0], AB, cXLATGam, NULL ); 
    if (!rHistC) { ErrMsg( "CLAHERef HistC=nil" ); return -1; }; // histogram of 'counts' or of 'refs_sum'
    rHists = BS_GamCLAHEHists( AW, AH, ADs[0], AB, cXLATGam, ARef ); // == build and treat histograms
    if (!rHists) { ErrMsg( "CLAHERef Hists=nil" ); return -1; }; // histogram of 'counts' or of 'refs_sum'
    BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDRefAve, this ); // - invoke iHistRefAve...
    BS_SpawnThreads( 0, 0, cXLATGam, ppr_CLAHEDCrB, this ); // if (qDbg) iDbgPrintHists( 3, qDbg ); // - crossblur
    if (rBlr) { // - blur if Blr<>0
      BS_SpawnThreads( 0, 0, cNxy, ppr_CLAHEDBlr, this ); }; // if (qDbg) iDbgPrintHists( 4, qDbg );
    int qRet = 1; BS_GamCLAHEApply( AW, AH, ADs[0], AB, cXLATGam, (int*)rHists ); // == apply CLAHE
    if (0) BS_GamCLAHEApply( AW, AH, ADs[qRet++], AB, cXLATGam, (int*)rHistC ); // == apply CLAHE !!!!
    delete rHistC; delete rHists;
    return qRet; } // end of TraRef
}; // end of tBS2_CLAHED

//// exports ////
int BS4_CLAHE_Transform( int AB, int AGam, int ABlr, const tBS_Val *AFrmOrg, // -- plain CLAHE 
  tBS_Val **ADs, bool ACUDA, const char *ADbgPfx, tBS4_CLAHE_Exts *AExts, const tBS_Val *ARef ) {
  tBS4_CLAHED *qCL = new tBS4_CLAHED(); qCL->DefineExts( AExts ); 
  // qCL->DbgSet( ADbgPfx ); // enable internal debug dumps ~"DbgPfx_A.raw"
  int qRes = qCL->Transform( AB, AGam, ABlr, AFrmOrg, ARef, ADs, ACUDA ); delete qCL;
  return qRes; }

// -- normalized frame for misc CLAHE-procedures
void BS4_FrmNrm_CPU( int AGam, tBS_Val *AFrm, tBS_Val **ATs2, int AW, int AH ) { tBS4_CLAHED::FrmNrm_CPU( AGam, AFrm, ATs2, AW, AH ); };
void BS4_FrmNrm_CUDA( int AGam, tBS_Val *AFrm, tBS_Val **ATs2 ) { tBS4_CLAHED::FrmNrm_CUDA( AGam, AFrm, ATs2 ); };

int BS4_CLAHE_TraAveRef( int AB, int ABlr, const tBS_Val *AFrmNrm, // -- ref_averages
  tBS_Val **ADs, bool ACUDA, tBS_Val *ARef0, tBS_Val *ARef1, tBS_Val *ARef2, tBS_Val *ARef3 ) {
  tBS_Val *qRefs[4] = {ARef0,ARef1,ARef2,ARef3}; int qNRefs = 1; for (;qNRefs<mcASz(qRefs);qNRefs++) if (!qRefs[qNRefs]) break;
  tBS4_CLAHED *qCL = new tBS4_CLAHED; qCL->DefineExts( NULL );
  int qRes = qCL->TraAveRef( AB, ABlr, AFrmNrm, qNRefs, qRefs, ADs, ACUDA ); delete qCL;
  return qRes; }

int BS4_CLAHE_TraRefAir( int AB, int ABlr, const tBS_Val *AFrmNrm, // -- ref_average for voluntary dimension
  const tBS_Val *ARef, tBS_Val **ADs, int AW, int AH ) { 
  tBS4_CLAHED *qCL = new tBS4_CLAHED; qCL->DefineExts( NULL );
  int qRes = qCL->TraRefAir( AB, ABlr, AFrmNrm, ARef, ADs, AW, AH ); delete qCL;
  return qRes; }
