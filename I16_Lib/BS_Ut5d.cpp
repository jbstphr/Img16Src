// I16_Ut5d.cpp : mask-recount (CPU-only)
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut5.h"
#include "BS_Ut7.h"
#include "BS_Ut5d.h"

//// commons with the scope of local 'DSeg' ////  

static int ixMskEnabled = 0; int BS5d_MskEnabled() { return ixMskEnabled; } // -- >0:is enabled

static struct { // used by 'Wrapper.mskGDi_ReCount'
  int rT[3], // trinarization thresholds, optional CLAHE's block-size
      rC[5]; // Cuts
} ixGDi = {
  {299, 399, 0}, // Thr0 Thr1 Blk
  {-123, 123, 156, 0, 102} // rg0 rg1 tg0 tg1 Rad
};

void BS5d_MskLoad( GIniFile *AF ) { // load estimation parameters (called by I16_Init)
  ixMskEnabled = atoi( AF->SGetStr( "nEnabled", "0" ) ); if (ixMskEnabled<=0) return;
  if (xPars[cParIniSetDflt]>0) { char qS[256]; // only override ini-profile
    sprintf( qS, "%d %d; t0 t1 Blk", ixGDi.rT[0], ixGDi.rT[1] );
    AF->SetStr( "nGDi3", qS );
    sprintf( qS, "%d %d %d %d %d; rg0 rg1 tg0 tg1 rad", ixGDi.rC[0], ixGDi.rC[1], ixGDi.rC[2], ixGDi.rC[3], ixGDi.rC[4] );
    AF->SetStr( "nGDc5", qS );
    return; }; int qI[16];
  if (Gini_sscanf(AF->DGetStr("nGDi3",NULL),qI)>0) memcpy( ixGDi.rT, qI, sizeof(ixGDi.rT) );
  if (Gini_sscanf(AF->DGetStr("nGDc5",NULL),qI)>0) memcpy( ixGDi.rC, qI, sizeof(ixGDi.rC) );
}

// Input & output for "process mskGDi by ixGDi":
//   inp: normalized mskGDi "eFd2.raw"
//   out: mskGDi redefined by spatial counts
// Elementary process-flow:
//   stage a) trinarization by thresholds: "-1 in flat areas", "1 in edges",
// "0 otherwize".
//   stage b) TraceShader counts intervals of const sign by 8 directions;
// PixelShader determines the "metric of flatness" (e.g. ~max of signed counts).
//   stage c) Normalize & blur the "metric of flatness".
//   ca) Optional CLAHE.
//   cb) Cut.
//   cd) Blur (if "ca)" is ON then "BlrRad~<5" else "BlrRad~>5".

struct tMskWr { // wrapper-class
  // -- debug dumps
  char     rDbgPfx[256], *rDbgSfx; // debug-dump's strings
  void DbgSet( const char *ADbgPfx ) { // debug-trace
    if (ADbgPfx) rDbgSfx = rDbgPfx+sprintf( rDbgPfx, "%s_", ADbgPfx ); }

  bool     rCUDA; int rW, rH; // frame's dimension
  // -- the frame<s> to update (and/or to consider)
  tBS_Val *rImg, *rmskThk, *rmskGDi; 

  // -- recount mskGDi by parameters in "ixGDi"
  
  #if 0 // interval:=length
  static int icplGDi_TraceShader( int AL, tBS_Val *ADat ) { // TraceShader
    mcIvlC( AL, ADat ); return AL; }
  #else // pos_interval:=length, neg_interval:=dist_from_center
  static int icplGDi_TraceShader( int AL, tBS_Val *ADat ) { // TraceShader
    int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; // simplify loop's condition
    if (qV>0) { if (qL<0) { qV = qL; qL = 1; goto lab1; }; qL++; continue; };
    if (qV<0) { if (qL>0) { qV = qL; qL =-1; goto lab1; }; qL--; continue; };
    if (!qL) continue; qV = qL; qL = 0; lab1: if (qV>0) BS_VFil( qV, ADat+qi-qV, qV ); else {
    qV = -qV; tBS_Val *qDat1 = ADat+qi, *qDat0 = qDat1-qV; qV >>= 1; qDat0[qV] = -(qV<<1);
    for (int qj=0;qj<qV;qj++) qDat0[qj] = qDat1[-qj] = -(qV+qj+1); }; }; return AL; } // for i
  #endif

  static int __fastcall icpbGDi_PixelShader( tBS_Val *ADat ) { // PixelShader
    int qV = 0; for (int qi=0;qi<4;qi++) qV += __max(ADat[qi],ADat[qi+4]); return qV; }

  void imskGDi_ReCount( const tBS_Val *AImg, const tBS_Val *AmskThk, tBS_Val *AmskGDi, tBS_Val **ATs ) {
    const int cW = rW, cH = rH, cN = cW*cH; // recount mskGDi by ixGDi in local DSeg
    tBS_Val     *qT = ATs[0]; BS_Cpy( cN, AmskGDi, qT );
    { const int  cT0 = ixGDi.rT[0], cT1 = ixGDi.rT[1]; // ~trinarization by threshold
      for (int qi=0;qi<cN;qi++) { int qV = qT[qi]; if (qV>cT1) qV = 1; else if (qV>cT0) qV = 0;
        else qV = -1; qT[qi] = qV; } // "-1 in flat areas", "1 in known edges", "0 otherwise"
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "GDiA_tt.raw" ); BS_Write( rDbgPfx, cN, qT ); }; };
    { // get the "metric of flatness"
      BS_DataDup( 7, cN, qT, ATs+1 );
      BS_CPLA( icplGDi_TraceShader, 8, cBS_AllDirsT, cW, cH, ATs ); // TraceShader counts flat areas
      BS_CPB( icpbGDi_PixelShader, 8|cCPBNDD, cN, ATs ); // PixelShader determines ~mincount
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "GDiB_mc.raw" ); BS_Write( rDbgPfx, cN, qT ); }; };
    { // optionally blur the "metric of flatness" by CLAHE 
      const int  cB = ixGDi.rT[2]; if (cB>4 && cB<9) { // blocks' size
      BS4_CLAHE_TraAveRef( cB, 0, AImg, ATs+1, false, qT ); BS_Cpy( cN, ATs[1], qT );
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "GDiCa_clab.raw" ); BS_Write( rDbgPfx, cN, qT ); }; }; };
    { // normalize & finally blur the "metric of flatness"
      const int *cC = ixGDi.rC, cTg0 = cC[2], cRad = cC[4]; // masks within 0..NrmAve
      tBS2_Cut   cCut = { {cC[0], cC[1], 512-cTg0-cC[3]}, {0,0} };
      BS2_Cut( &cCut, qT ); if (cTg0) BS_VAdd( cN, qT, cTg0 );
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "GDiCb_msk.raw" ); BS_Write( rDbgPfx, cN, qT ); };
      BS2_MovAve2Dp( cRad%100, qT, ATs+1, cRad/100, cW, cH );
      if (1 && rDbgSfx) { strcpy( rDbgSfx, "GDiCc_mskb.raw" ); BS_Write( rDbgPfx, cN, qT ); }; };
    // leave
    BS_Cpy( cN, qT, AmskGDi ); }

  // -- process
  void Process() {
    const int cW = rW, cH = rH, cN = cW*cH;
    tBS_Val  *qPs[3]; 
    if (rCUDA) {
      BS_DatasAlloc( mcASz(qPs), cN, qPs );
      BC_Move( -cN*cBS_Sz, qPs[0], rImg );
      BC_Move( -cN*cBS_Sz, qPs[1], rmskThk );
      BC_Move( -cN*cBS_Sz, qPs[2], rmskGDi );
    } else {
      qPs[0] = rImg; qPs[1] = rmskThk; qPs[2] = rmskGDi;
    };
    tBS_Val  *qTs[8]; BS_DatasAlloc( mcASz(qTs), cN, qTs ); // allocate temp_limbo
    if (1 && rDbgSfx) { strcpy( rDbgSfx, "AinpGDi.raw" ); BS_Write( rDbgPfx, cN, qPs[2] ); };
    imskGDi_ReCount( qPs[0], qPs[1], qPs[2], qTs );
    if (rCUDA) {
      BC_Move( cN*cBS_Sz, qPs[2], rmskGDi );
      BS_DatasRelease( mcASz(qPs), qPs );
    };
    BS_DatasRelease( mcASz(qTs), qTs ); } // release temp_limbo
}; // end of tMskWr

static tMskWr *jxMkWrSpawn1; static void jMkWrSpawn1( int ATag ) { if (ATag==1) jxMkWrSpawn1->Process(); }
  // it could be inline or spawned (via the only exemplar of tAiWr)

void BS5d_Msk( int AW, int AH, bool ACUDA, const char *ADbgPfx, // -- forks mask-recount thread
  const tBS_Val *AImg, const tBS_Val *AmskThk, tBS_Val *AmskGDi ) { // these parameters may change by implementation
  // all implementations update mskGDi
  static tMskWr qMkWr; memset( &qMkWr, 0, sizeof(tMskWr) ); // this static locates in DSeg
  qMkWr.DbgSet( ADbgPfx );
  qMkWr.rW = AW; qMkWr.rH = AH; qMkWr.rCUDA = ACUDA;
  qMkWr.rImg = (tBS_Val*)AImg; qMkWr.rmskThk = (tBS_Val*)AmskThk; qMkWr.rmskGDi = AmskGDi; 
  jxMkWrSpawn1 = &qMkWr; BS_Spawn1( jMkWrSpawn1, 3 ); } //  inline execution "jMkWrSpawn1( 1 )"

void BS5d_MskWait() { BS_Spawn1Wait( 3 ); } // -- wait for forked mask-recount thread

