// by S.B.
#include "stdafx.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut2a.h" // Skew/Shear
#include "math.h"

//== Ticks' mark ==
//=================
CBS_Ticks::CBS_Ticks() { rPfx[0] = 0; Reset(); }
void CBS_Ticks::Reset() { rTicks0 = rTicks = BS_Ticks(); }
void CBS_Ticks::Mark( const char *AMsg ) { int qTicks = BS_Ticks(); if (rPfx[0]) { char qS[128];
  if (AMsg) sprintf( qS, cTPS" %s %.2f \"%s\"", rPfx, abs(qTicks-rTicks)*0.001f, AMsg );
       else { sprintf( qS, cTPS" %s %.2f", rPfx, abs(qTicks-rTicks0)*0.001f ); rTicks0 = qTicks; };
  WriteLn( qS ); }; rTicks = qTicks; } // if Pfx[0]

//== 2D-iterators-service (array-accessors) ==
//============================================
static void iLnSetDim( tBS_LIRec &AR, int AW, int AH ) { AR.rW = AW; AR.rH = AH; AR.rMaxLen = BS_Align( mcsqrti(AW*AW+AH*AH)+2 ); AR.rPeriod = 0; } // set array's dimension
static void inline imcLnSetInc( tBS_LIRec &AR, int AXInc, int AYInc, int APeriod ) { AR.rXInc = AXInc; AR.rYInc = AYInc; AR.rPeriod = APeriod; }
static void inline imcLnSetISt( tBS_LIRec &AR, int AStX, int AStY, int AStXInc, int AStYInc ) { AR.rStX = AStX<0?AR.rW+AStX:AStX; AR.rStY = AStY<0?AR.rH+AStY:AStY; AR.rStXInc = AStXInc; AR.rStYInc = AStYInc; }
static void iLnSetDir( tBS_LIRec &AR, int ADir, int APass ) { AR.rDir = ADir; AR.rPass = APass; switch (ADir) { // initialize line-iterator-record
  case  1: imcLnSetInc( AR, 1, 0, 0 ); imcLnSetISt( AR, 0, 0, 0, 1 ); break; // horz '-' lft->rgt
  case  5: imcLnSetInc( AR, 0, 1, 0 ); imcLnSetISt( AR, 0, 0, 1, 0 ); break; // vert '|' top->btm
  case  3: imcLnSetInc( AR, 1, 1, 0 ); if (!APass) imcLnSetISt( AR, 0, 0, 1, 0 ); else imcLnSetISt( AR, 0, 1, 0, 1 ); break; // diag '\' lft,top->rgt,btm
  case  7: imcLnSetInc( AR,-1, 1, 0 ); if (!APass) imcLnSetISt( AR, 0, 0, 1, 0 ); else imcLnSetISt( AR,-1, 1, 0, 1 ); break; // diag '/' rgt,top->lft,btm
  case    2: imcLnSetInc( AR, 1, 0, 3 ); if (!APass) imcLnSetISt( AR, 0, 0, 0, 1 ); else imcLnSetISt( AR, 3, 0, 3, 0 ); break; // diag '\-' lft,top->rgt,btm
  case  4: imcLnSetInc( AR, 0, 1, 3 ); if (!APass) imcLnSetISt( AR, 0, 0, 1, 0 ); else imcLnSetISt( AR, 0, 3, 0, 3 ); break; // diag '|\' lft,top->rgt,btm
  case  6: imcLnSetInc( AR, 0, 1,-3 ); if (!APass) imcLnSetISt( AR,-1, 0,-1, 0 ); else imcLnSetISt( AR,-1, 3, 0, 3 ); break; // diag '/|' rgt,top->lft,btm
  default: imcLnSetInc( AR,-1, 0, 3 ); if (!APass) imcLnSetISt( AR,-1, 0, 0, 1 ); else imcLnSetISt( AR,-4, 0,-3, 0 ); break; // diag '-/' rgt,top->lft,btm
  }; } // switch

static int iLnIndexes( tBS_LIRec &AR, int *AIndexes ) { // retursns straight-line's reindex
  unsigned  qX = (unsigned)AR.rStX, qW = (unsigned)AR.rW; if (qX>=qW) return 0;
  unsigned  qY = (unsigned)AR.rStY, qH = (unsigned)AR.rH; if (qY>=qH) return 0;
  const int cXInc = AR.rXInc, cYInc = AR.rYInc; int qL = 0, qPeriod = AR.rPeriod; if (!qPeriod) {
    do { AIndexes[qL++] = qY*qW+qX; qX += cXInc; qY += cYInc; } while (qX<qW && qY<qH);
  } else { int qPInc = 1; if (qPeriod<0) { qPeriod = -qPeriod; qPInc = -1; };
    if (cXInc) do { for (int qi=0;qi<qPeriod;qi++) { AIndexes[qL++] = qX+qW*qY; qX += cXInc; if (qX>=qW) goto labDone; }; qY += qPInc; } while (qY<qH);
          else do { for (int qi=0;qi<qPeriod;qi++) { AIndexes[qL++] = qX+qW*qY; qY += cYInc; if (qY>=qH) goto labDone; }; qX += qPInc; } while (qX<qW); };
  labDone: AR.rStX += AR.rStXInc; AR.rStY += AR.rStYInc; return qL; } // advance start-coordinates of straight-line
static void itppr_LnIterator( tBS_SpawnThreadRec *AP ) { // invokes LnSpawnProc for each line
  tBS_SpawnLIProc    qProc = (tBS_SpawnLIProc)AP->rP[0];
  const int         *cDirs =(const int*)AP->rP[1];
  int                qW = mcVi(AP->rP[2]), qH = (qW>>16)&0xFFFF; qW &= 0xFFFF;
  tBS_SpawnLIRec     qSPR; iLnSetDim( qSPR.rLIR, qW, qH );
  int               *qBufI = new int[qSPR.rLIR.rMaxLen];
  qProc( -1, 0, &qSPR, qBufI, (tBS_Val**)AP->rP[3] ); // LnSpawnProc(-1)
  const int cLst = AP->rLst; for (int qDiri=AP->rFst;qDiri<cLst;qDiri++) {
    qSPR.rDir = cDirs[qDiri]&0xF; qSPR.rLC = 0;
    qProc( -3, qDiri, &qSPR, qBufI, (tBS_Val**)AP->rP[3] ); // LnSpawnProc(-3)
    int              qPass = (cDirs[qDiri]&0xF0)>>4, qPassE;
    if (qPass) { qPassE = qPass; qPass--; } // whether a dir has two passes
          else qPassE = cBS_LnDirPass[qSPR.rDir-1];
    do { iLnSetDir( qSPR.rLIR, qSPR.rDir, qPass ); 
         lab1: int qLL = iLnIndexes( qSPR.rLIR, qBufI ); // LnSpawnProc(>0)
         if (qLL) { qBufI[qLL] = 0; qProc( qLL, qDiri, &qSPR, qBufI, (tBS_Val**)(AP->rP[3]) ); qSPR.rLC++; goto lab1; };
    } while (++qPass<qPassE);
    qProc( -4, qDiri, &qSPR, qBufI, (tBS_Val**)(AP->rP[3]) ); // LnSpawnProc(-4)
  }; // for
  qProc( -2, 0, &qSPR, qBufI, (tBS_Val**)(AP->rP[3]) ); // LnSpawnProc(-2)
  delete qBufI; }

void BS_SpawnLI( int ANDD, const int *ADirs, tBS_SpawnLIProc APrcc, int AW, int AH, tBS_Val **ADatas ) {
  // wrapper of LnIterator for 2D-array 'WxH' (with versatile ad-hoc parameters)
  BS_SpawnThreads( cBS_SpawnPool, 0, ANDD, itppr_LnIterator, (const void*)APrcc, ADirs, mcIv(AW|(AH<<16)), ADatas ); }

int BS_GetProjectionBufI( int ADir, int AW, int AH, int *ABufI ) {
  int qL = 0; tBS_LIRec qLIR; iLnSetDim( qLIR, AW, AH );
  for (int qPass=0;qPass<cBS_LnDirPass[ADir-1];qPass++) { iLnSetDir( qLIR, ADir, qPass );
    labRep: ABufI[qL++] = qLIR.rStX+qLIR.rStY*AW;
    qLIR.rStX += qLIR.rStXInc; if (qLIR.rStX<0 || qLIR.rStX>=AW) goto labCont;
    qLIR.rStY += qLIR.rStYInc; if (qLIR.rStY<0 || qLIR.rStY>=AH) goto labCont; goto labRep;
    labCont: if (!qPass) switch (ADir) { // implementation-specific reshuffle
      case 2: case 3: case 4: case 6: case 8: mcFlip<int>( qL, ABufI ); break;
    }; };
  switch (ADir) { // implementation-specific reshuffle
    case 2: case 8: mcFlip<int>( qL, ABufI ); break; };
  return qL; }

//====== areas' consolidation pipeline =======
//============================================
// 2D-array accessors for areas' consolidation pipeline (for ad-hoc usage as well)
static inline tBS_Val *mcCPLBuf( int AL ) { return (tBS_Val*)malloc( BS_Align(AL)*cBS_Sz*4 ); } // reserve ad-hoc buffer's 
static void ilpr_CPLA( int ALL, int ADiri, tBS_SpawnLIRec *ASPR, int *ABufI, tBS_Val **ADatas ) { // CPL step#A
  if (ALL<0) { switch (ALL) {
    case -1: ASPR->rDat[0] = mcCPLBuf( ASPR->rLIR.rMaxLen ); break;
    case -2: free(ASPR->rDat[0]); break; }; // switch
    return; }; // if LL<0
  tBS_Val *qDataBuf = (tBS_Val*)ASPR->rDat[0], *qData = ADatas[1+ADiri]; BS_Gath( ALL, qData, qDataBuf, ABufI );
  BS_Scat( ((tBS_CPLProc)ADatas[0])(ALL,qDataBuf), qDataBuf, qData, ABufI ); }
void BS_CPLA( tBS_CPLProc AProc, int ANDD, const int *ADirs, int AW, int AH, tBS_Val **ADatas ) {
  if (ANDD<3) { int qOnlyRows = -1, qOnlyCols = -1; for (int qj=0;qj<ANDD;qj++) switch (ADirs[qj]) {
      case 1: qOnlyRows = qj; break; case 5: qOnlyCols = qj; break; default: goto lab1; }; // switch
    if (qOnlyRows>=0) BS_CPLCRows( AProc, AW, AH, ADatas[qOnlyRows] );
    if (qOnlyCols>=0) BS_CPLCCols( AProc, AW, AH, ADatas[qOnlyCols] ); return; }; lab1: // if NDD<3
  tBS_Val *qDs[17]; qDs[0] = (tBS_Val*)AProc; memcpy( qDs+1, ADatas, ANDD*sizeof(void*) ); BS_SpawnLI( ANDD, ADirs, ilpr_CPLA, AW, AH, qDs ); }

#define mctppr_CPB1( NDD ) { for (int qi=AP->rFst;qi<cLst;qi++) { for (int qj=0;qj<NDD;qj++) qBuf[qj] = qDatas[qj][qi]; qDatas[0][qi] = qProc( qBuf ); }; }
static void itppr_CPB( tBS_SpawnThreadRec *AP ) { const int cLst = AP->rLst, cNDD = mcVi(AP->rP[1]); int qNDD = cNDD&(cCPBNDD-1); // pipeline step #B (by a pixel)
  tBS_CPBProc qProc = (tBS_CPBProc)AP->rP[0]; tBS_Val *qDatas[cCPBNDD]; memcpy( qDatas, AP->rP[2], qNDD*sizeof(void*) ); tBS_Val qBuf[cCPBNDD];
  if (!(cNDD&cCPBNDD)) { for (int qi=AP->rFst;qi<cLst;qi++) {
    for (int qj1=0;qj1<qNDD;qj1++) qBuf[qj1] = qDatas[qj1][qi]; int qRes = qProc( qBuf ); for (int qj2=0;qj2<qRes;qj2++) qDatas[qj2][qi] = qBuf[qj2]; };
  } else switch (qNDD) { 
    case 1: mctppr_CPB1( 1 ); break; case 2: mctppr_CPB1( 2 ); break; case 3: mctppr_CPB1( 3 ); break;
    case 4: mctppr_CPB1( 4 ); break; case 5: mctppr_CPB1( 5 ); break; case 6: mctppr_CPB1( 6 ); break;
    case 7: mctppr_CPB1( 7 ); break; case 8: mctppr_CPB1( 8 ); break; default: mctppr_CPB1( qNDD );
  }; } // switch NDD
void BS_CPB( tBS_CPBProc AProc, int ANDD, int AN, tBS_Val **ADatas, void *AOut ) { tBS_Val *qBuf[cCPBNDD]; memcpy( qBuf, ADatas, (ANDD&(cCPBNDD-1))*sizeof(void*) );
  if (AOut) { int qNOut = (ANDD>>8)&(cCPBNDD-1); if (qNOut) for (int qo=0;qo<qNOut;qo++) qBuf[qo] = ((tBS_Val**)AOut)[qo]; else { qBuf[0] = (tBS_Val*)AOut; qNOut = 1; };
  for (int qo=0;qo<qNOut;qo++) BS_Cpy( AN, ADatas[qo], qBuf[qo] ); }; BS_SpawnThreads( cBS_SpawnPool, 0, AN, itppr_CPB, (void*)AProc, mcIv(ANDD&255), qBuf ); } // if Out

static inline void mcCPLSpawn( int AD1, int AD2, tBS_SpawnThreadProc AProc0, void *AP0, void *AP1, void *AP2, void *AP3 ) {
  BS_SpawnThreads( AD1<=99||AD1*AD2<2048?1:cBS_SpawnPool, 0, AD1, AProc0, AP0, AP1, AP2, AP3 ); };
static void itppr_CPLCRows( tBS_SpawnThreadRec *AP ) { // pipeline step #C (along rows)
  int qW = mcVi(AP->rP[1]); const int cFst = AP->rFst; tBS_Val *qData = ((tBS_Val*)AP->rP[2])+qW*cFst, *qDataBuf = mcCPLBuf( qW ); const int cLst = AP->rLst;
  for (int qi=cFst;qi<cLst;qi++) { BS_Cpy( qW, qData, qDataBuf ); BS_Cpy( ((tBS_CPLProc)AP->rP[0])(qW,qDataBuf), qDataBuf, qData ); qData += qW; }; delete qDataBuf; }
static void itppr_CPLCRowsNonBuf( tBS_SpawnThreadRec *AP ) { // pipeline step #C (along rows, non buffered)
  int qW = mcVi(AP->rP[1]); const int cFst = AP->rFst; tBS_Val *qData = ((tBS_Val*)AP->rP[2])+qW*cFst; const int cLst = AP->rLst; for (int qi=cFst;qi<cLst;qi++) { ((tBS_CPLProc)AP->rP[0])( qW, qData ); qData += qW; }; }
void BS_CPLCRows( tBS_CPLProc AProc, int AW, int AH, tBS_Val *AData, bool ANonBuf ) { mcCPLSpawn( AH, AW, ANonBuf?itppr_CPLCRowsNonBuf:itppr_CPLCRows, (void*)AProc, mcIv(AW), AData, NULL ); }
static void itppr_CPLCCols( tBS_SpawnThreadRec *AP ) { // pipeline step #C (along columns)
  int qW = mcVi(AP->rP[1]), qH = qW>>16; qW &= 0xFFFF; const int cFst = AP->rFst; tBS_Val *qData = ((tBS_Val*)AP->rP[2])+cFst, *qDataBuf = mcCPLBuf( qH ); const int cLst = AP->rLst;
  for (int qi=cFst;qi<cLst;qi++) { mcGath<tBS_Val>( qH, qData, qDataBuf, qW ); mcScat<tBS_Val>( ((tBS_CPLProc)AP->rP[0])(qH,qDataBuf), qDataBuf, qData, qW ); qData++; }; delete qDataBuf; }
void BS_CPLCCols( tBS_CPLProc AProc, int AW, int AH, tBS_Val *AData ) { mcCPLSpawn( AW, AH, itppr_CPLCCols, (void*)AProc, mcIv(AW|(AH<<16)), AData, NULL ); }

template <int cR,typename tBuf,typename tSrc,typename tDst,typename tProc>static void mcBS_CPKThread( tBS_SpawnThreadRec *AP ) { // only inside BS_CPK
  int qT = mcVi(AP->rP[0]); const int cW = qT&0x7FFF, cH = (qT>>16)&0x7FFF, cD = cR*2+1, cFst = AP->rFst, cLst = AP->rLst;
  const tSrc *cSrc = (const tSrc*)AP->rP[1], *cSrcRows[cD]; // source array and columns' buffer
  for (int qy=0;qy<cD;qy++) { int qrow = cFst-cR+qy; if (qrow<0) qrow = 0; cSrcRows[qy] = cSrc+qrow*cW; };
  tProc qProc = (tProc)AP->rP[3]; // kernel procedure
  for (int qrow=cFst;qrow<cLst;qrow++) { tBuf qBuf[cD*cD]; // pixels' buffer [2R+1][2R+1]
    for (int qy=0;qy<cD;qy++) { tBuf *qBp = qBuf+cD*qy; const tSrc *cSp = cSrcRows[qy], csp0 = cSp[0]; // fill pixels' buffer for col#0
      for (int qx0=1;qx0<cR+2;qx0++) qBp[qx0] = csp0; for (int qx1=cR+2;qx1<cD;qx1++) qBp[qx1] = cSp[qx1-cR-1];
    }; // for y
    tDst *qDst = ((tDst*)AP->rP[2])+qrow*cW; // the destination's row
    for (int qcol=0;qcol<cW;qcol++) { memmove( qBuf, qBuf+1, (cD*cD-1)*sizeof(tBuf) ); // advance pixels' buffer
      int qx = qcol+cR; if (qx>=cW) qx = cW-1; for (int qy=0;qy<cD;qy++) qBuf[cD-1+cD*qy] = cSrcRows[qy][qx]; qDst[qcol] = (tDst)qProc( qBuf );
    }; // for col
    memmove( cSrcRows, cSrcRows+1, (cD-1)*sizeof(void*) ); // advance columns' buffer
    int qy = qrow+cR+1; if (qy>=cH) qy = cH-1; cSrcRows[cD-1] = cSrc+qy*cW;
  }; } // for row
void BS_CPK( int AR, tBS_CPKProc AProc, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { tBS_SpawnThreadProc qProc; switch (AR) {
  case 1: qProc = mcBS_CPKThread<1,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break; case 2: qProc = mcBS_CPKThread<2,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break;
  case 3: qProc = mcBS_CPKThread<3,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break; case 4: qProc = mcBS_CPKThread<4,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break;
  case 5: qProc = mcBS_CPKThread<5,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break; case 6: qProc = mcBS_CPKThread<6,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break;
  case 7: qProc = mcBS_CPKThread<7,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break; case 8: qProc = mcBS_CPKThread<8,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break;
  case 9: qProc = mcBS_CPKThread<9,tBS_Val,tBS_Val,tBS_Val,tBS_CPKProc>; break; default: { ErrMsg( "BS_CPK unsupp R" ); return; }; };
  mcCPLSpawn( AH, AW, qProc, mcIv((AH<<16)|AW), (void*)ASrc, ADst, (void*)AProc ); }

void BS_CPKInt( int AR, tBS_CPKProcInt AProc, int AW, int AH, const int *ASrc, tBS_Val *ADst ) { //
  tBS_SpawnThreadProc qProc = mcBS_CPKThread<1,int,int,tBS_Val,tBS_CPKProcInt>; 
  mcCPLSpawn( AH, AW, qProc, mcIv((AH<<16)|AW), (void*)ASrc, ADst, (void*)AProc ); }

//================ 1D arrays =================
//============================================
void BS_VMul( int AL, tBS_Val *ADat, float AV ) { // Dat *= V
  if (fabs(AV-1)<1e-4f) return; if (fabs(AV)<1e-4f) { BS_ZFil( AL, ADat ); return; };
  if (fabs(AV+1)<1e-4f) { BS_Neg( AL, ADat ); return; }; // Dat *= -1
  if (fabs(AV)>1.01f) { int qM = (int)(AV*256+0.4999f);  for (int qi=0;qi<AL;qi++) { int qv =( ADat[qi]*qM+127)>>8; ADat[qi] = mcRng(qv,-0x7FFF,0x7FFF); }; return; };
  if (fabs(AV)>1e-3f) { int qM = (int)(AV*1024+0.4999f); for (int qi=0;qi<AL;qi++) ADat[qi] = (ADat[qi]*qM+511)>>10; return; };
  for (int qi=0;qi<AL;qi++) ADat[qi] = (int)(ADat[qi]*AV+0.4999f); }
void BS_AddMul( int AL, const tBS_Val *AAdd, int *AAcc, float AV ) {
  if (fabs(AV-1)<1e-3f) { for (int qi=0;qi<AL;qi++) AAcc[qi] += AAdd[qi]; return; };
  if (fabs(AV+1)<1e-3f) { for (int qi=0;qi<AL;qi++) AAcc[qi] -= AAdd[qi]; return; };
  const int cSh = 8, cK = (int)(AV*(1<<cSh)+0.4999f); // Acc += V*Add
  if (!cK) return; for (int qi=0;qi<AL;qi++) AAcc[qi] += (AAdd[qi]*cK+(1<<(cSh-1))-1)>>cSh; }
bool BS_Rng( int AL, tBS_Val *ADat, int AMn, int AMx ) {
  if ((AMn<=-0x7FFF && AMx>=0x7FFF) || (AMn>=AMx)) return false;
  if (AMn<=-0x7FFF) { if (AMx<0x7FFF) for (int qi=0;qi<AL;qi++) { int qV = ADat[qi]; if (qV>AMx) ADat[qi] = AMx; }; return true; };
  if (AMx>= 0x7FFF) { for (int qi=0;qi<AL;qi++) { int qV = ADat[qi]; if (qV<AMn) ADat[qi] = AMn; }; return true; };
  for (int qi=0;qi<AL;qi++) { int qV = ADat[qi]; if (qV<AMn) ADat[qi] = AMn; else if (qV>AMx) ADat[qi] = AMx; }; return true; }

bool BS_RngA( int AL, tBS_Val *ALap, int AMx, int AMxP ) { if (AMx<=0) return false; 
  if (AMxP<1 || AMxP>=0x7FFF) for (int qi=0;qi<AL;qi++) { int qV = ALap[qi], qA = abs(qV)-AMx; // -Mx<V<Mx
    if (qA<=0) continue; if (qA>AMx) qA = AMx; qA = AMx+(qA>>2); ALap[qi] = qV<0?-qA:qA; }
  else for (int qi=0;qi<AL;qi++) { int qV = ALap[qi]; // -Mx<V<MxP
    if (-qV>AMx) { qV = -qV-AMx; if (qV>AMx) qV = AMx; ALap[qi] = -AMx-(qV>>2); } else
    if (qV>AMxP) { qV -= AMxP; if (qV>AMxP) qV = AMxP; ALap[qi] = AMxP+(qV>>2); }; }; return true; };

template <int cNrm> static void iMul( int AL, const tBS_Val *AMul, tBS_Val *ADat ) {
  const int cA = (1<<(cNrm-1))-1; for (int qi=0;qi<AL;qi++) ADat[qi] = (AMul[qi]*ADat[qi]+cA)>>cNrm; }
void BS_Mul( int AL, const tBS_Val *AMul, tBS_Val *ADat, int ANrm ) { // Mul:=Org*Mul/Nrm; Nrm<0 denotes "shr |Nrm|"
  if (ANrm>0) { for (int qi=0;qi<AL;qi++) ADat[qi] = AMul[qi]*ADat[qi]/ANrm; return; };
  switch (ANrm) {
  case  -1: iMul< 1>( AL, AMul, ADat ); break; case  -2: iMul< 2>( AL, AMul, ADat ); break; case  -3: iMul< 3>( AL, AMul, ADat ); break;  
  case  -4: iMul< 4>( AL, AMul, ADat ); break; case  -5: iMul< 5>( AL, AMul, ADat ); break; case  -6: iMul< 6>( AL, AMul, ADat ); break;  
  case  -7: iMul< 7>( AL, AMul, ADat ); break; case  -8: iMul< 8>( AL, AMul, ADat ); break; case  -9: iMul< 9>( AL, AMul, ADat ); break;
  case -10: iMul<10>( AL, AMul, ADat ); break; case -11: iMul<11>( AL, AMul, ADat ); break; case -12: iMul<12>( AL, AMul, ADat ); break;
  case -13: iMul<13>( AL, AMul, ADat ); break; case -14: iMul<14>( AL, AMul, ADat ); break; case -15: iMul<15>( AL, AMul, ADat ); break;
  default : for (int qi=0;qi<AL;qi++) { int qv = AMul[qi]*ADat[qi]; ADat[qi] = mcRg(qv); }; }; }
    
//void BS_Beth( int AL, const tBS_Val *ADat0, tBS_Val *ADat, const tBS_Val *AMul, int ANrm ) { // Dat:=Dat0+(Dat-Dat0)*Mul/Nrm
//  if (ADat0) BS_SubA( AL, ADat0, ADat ); BS_Mul( AL, AMul, ADat, ANrm ); if (ADat0) BS_Add( AL, ADat0, ADat ); }

//======== Gray-scale gamut's XLAT ==========
//===========================================
void BS_GamHist( int AW, int AH, const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam, const tBS_Val *ARef ) { memset( AHist, 0, AGam*sizeof(AHist[0]) ); AGam--;
  if (ARef) for (int qr=0;qr<AH;qr++) { for (int qc=0;qc<AW;qc++) { int qV = ADat[qc]; AHist[mcRng(qV,0,AGam)] += ARef[qc]; }; ADat += AStride; ARef += AStride; }
       else for (int qr=0;qr<AH;qr++) { for (int qc=0;qc<AW;qc++) { int qV = ADat[qc]; AHist[mcRng(qV,0,AGam)]++; }; ADat += AStride; }; }

// fixed code needs: 'K0=V[1]'<=0, 'K1=V[3]'>=0
static void iTvk( int AIdxDir, const int *AX, int AGam, int ATD, int *AV ) { // determine terminal_derivation
  const int *cX = AIdxDir>0?AX+0:AX+AGam-1; AV[0] = cX[0]; // 'V','K' in XLATApply
  if (ATD>0) AV[1] = (cX[0]-cX[__min(ATD,AGam>>1)*AIdxDir])*cBS_XLatTD/ATD; else {
  if (AX[0]>AX[AGam-1]) AIdxDir *= -1; AV[1] = __max(-cBS_XLatTD,ATD)*AIdxDir; }; } // if descend then inverse
void CBS_XLat::Init( const int *AXLat, int AMn, int AGam, int ATD0, int ATD1 ) { cX = AXLat; rM = AMn; rG = AGam;
  iTvk( 1, AXLat, AGam, ATD0, rV+0 ); iTvk( -1, AXLat, AGam, ATD1, rV+2 ); }

static void itppr_XLATApply( tBS_SpawnThreadRec *AP ) { CBS_XLat *qThis = (CBS_XLat*)AP->rP[0]; tBS_Val *qDat = (tBS_Val*)AP->rP[1]; 
  const int cV0 = qThis->rV[0], cK0 = qThis->rV[1], cV1 = qThis->rV[2], cK1 = qThis->rV[3], cMn = qThis->rM, cGam1 = qThis->rG-1, *cXLAT = qThis->cX, cLst = AP->rLst;
  for (int qi=AP->rFst;qi<cLst;qi++) { int qv = qDat[qi]-cMn; if (qv<0) { qv = cV0-(qv*cK0>>cBS_XLatTDLg); qDat[qi] = qv<-0x7FFF?-0x7FFF:qv; continue; };
    if (qv>cGam1) { qv = cV1+((qv-cGam1)*cK1>>cBS_XLatTDLg); qDat[qi] = qv>0x7FFF?0x7FFF:qv; continue; }; qDat[qi] = cXLAT[qv]; }; } // for i
void CBS_XLat::Apply( int AN, tBS_Val *ADat ) { if (cX) BS_SpawnThreads( cBS_SpawnPool|2, 0, AN, itppr_XLATApply, this, ADat ); }

void BS_GamXLATFit( int AGam, int *AXLAT ) { AGam--; int qMx = AXLAT[AGam]; if (!qMx || qMx==AGam) return; // scale so that XLAT[Gam-1]:=Gam-1
  const int cSh = 8, cScale = (AGam<<cSh)/qMx; for (int qj=1;qj<=AGam;qj++) AXLAT[qj] = (AXLAT[qj]*cScale)>>cSh; }
       
static const int cBmn = 5; // min block's dim (as power of two)
#define mcCLAHEbl(B,Dim,N,D0) { N = Dim>>B; D0 = (Dim-(N<<B))>>1; } // block's macro
typedef struct { int rNx,rNy,rX0,rY0,rSz,rStride; unsigned rGam; } tCLAHE; // hists[Ny][Nx]; skip X0 from lft, Y0 from top
static bool iCLAHE( int AW, int AH, int AB, int AGam, tCLAHE *AP ) { if (AB<cBmn || AB>cBmn+5) return false; // get blocks
  mcCLAHEbl( AB, AW, AP->rNx, AP->rX0 ); mcCLAHEbl( AB, AH, AP->rNy, AP->rY0 ); AP->rSz = AGam*AP->rNx*AP->rNy; AP->rStride = AW; AP->rGam = AGam; return AP->rNx>1 && AP->rNy>1; }
template <int cB> static void ippr_CLAHEHists( tBS_SpawnThreadRec *AP ) { const tCLAHE *cP = (tCLAHE*)AP->rP[0];
  const tBS_Val *cDat = (tBS_Val*)AP->rP[1]; unsigned *qHists = (unsigned*)AP->rP[2]; const unsigned cGam = cP->rGam;
  const int      cNx = cP->rNx, cStride = cP->rStride, cFst = AP->rFst, cLst = AP->rLst; cDat += cStride*cFst; 
  const tBS_Val *cRef = (const tBS_Val*)AP->rP[3]; if (cRef) { cRef += cStride*cFst; // here is a reference
    for (int qrow=cFst;qrow<cLst;qrow++) { const int cNy0 = (qrow>>cB); // blk by Y
      for (int qcol=0;qcol<(cNx<<cB);qcol++) { const int cNx0 = (qcol>>cB); unsigned qV = cDat[qcol]; if (qV>=cGam) qV = cGam-1; // blk by X
        #ifdef WIN32
          InterlockedExchangeAdd( ((LONG*)qHists)+(cNy0*cNx+cNx0)*cGam+qV, cRef[qcol] );
        #else // GCC
          __sync_fetch_and_add( qHists+(cNy0*cNx+cNx0)*cGam+qV, cRef[qcol] ); 
        #endif
      }; cDat += cStride; cRef += cStride; }; // for col,row
  } else { // no reference; commonplace histogram of itself
    for (int qrow=cFst;qrow<cLst;qrow++) { const int cNy0 = (qrow>>cB); // blk by Y
      for (int qcol=0;qcol<(cNx<<cB);qcol++) { const int cNx0 = (qcol>>cB); unsigned qV = cDat[qcol]; if (qV>=cGam) qV = cGam-1; // blk by X
        #ifdef WIN32
          InterlockedIncrement( ((LONG*)qHists)+(cNy0*cNx+cNx0)*cGam+qV );
        #else // GCC
          __sync_fetch_and_add( qHists+(cNy0*cNx+cNx0)*cGam+qV, 1 );
        #endif
        }; cDat += cStride; }; // for col,row
  }; }
unsigned *BS_GamCLAHEHists( int AW, int AH, const tBS_Val *ADat, int AB, int AGam, const tBS_Val *ARef ) { // -- get Hists[Ny][Nx]
  if (AB<=0) { unsigned *qHist = new unsigned[AGam]; AB *= -1; int qo = AB+AB*AW; // whole frame (allo margin)
    BS_GamHist( AW-AB*2, AH-AB*2, ADat+qo, AW, qHist, AGam, ARef?ARef+qo:NULL ); return qHist; };
  tCLAHE qP; if (!iCLAHE( AW, AH, AB, AGam, &qP )) return NULL; // get blocks
  unsigned *qHists = new unsigned[qP.rSz]; memset( qHists, 0, qP.rSz*sizeof(unsigned) ); tBS_SpawnThreadProc qProc; switch (AB) {
    case cBmn+0:qProc = ippr_CLAHEHists<(cBmn+0)>; break; case cBmn+1:qProc = ippr_CLAHEHists<(cBmn+1)>; break;
    case cBmn+2:qProc = ippr_CLAHEHists<(cBmn+2)>; break; case cBmn+3:qProc = ippr_CLAHEHists<(cBmn+3)>; break;
    case cBmn+4:qProc = ippr_CLAHEHists<(cBmn+4)>; break; case cBmn+5:qProc = ippr_CLAHEHists<(cBmn+5)>; break;
    default: return NULL; }; const int co = qP.rY0*AW+qP.rX0; ADat += co; if (ARef) ARef += co;
  BS_SpawnThreads( cBS_SpawnPool, 0, qP.rNy<<AB, qProc, &qP, ADat, qHists, ARef );
  #if 0
    const int cNxy = qP.rNx*qP.rNy, cNp = 1<<(AB*2); for (int qn=0;qn<cNxy;qn++) {
      int qZ = mcSum<unsigned>( AGam, qHists+AGam*qn ); if (qZ!=cNp) ErrMsg( "GamCLAHEHists" ); };
  #endif
  return qHists; }

template <int cB> static inline void mcCLAHEio( int N, int D, int P, int *IO ) { const int M = (1<<cB)+(1<<(cB-1))+D; P -= M; // index,offset-macro
  if (P<=0) { IO[0] = 0; IO[1] = ((P+M)<<cB)/M; return; }; P += 1<<cB; D = P>>cB;
  if (D>=N-2) { IO[0] = N-2; IO[1] = ((P-((N-2)<<cB))<<cB)/M; return; }; IO[0] = D; IO[1] = P-(D<<cB); }
template <int cB> static void ippr_CLAHEApply( tBS_SpawnThreadRec *AP ) { const tCLAHE *cP = (tCLAHE*)AP->rP[0];
  tBS_Val  *qDat = (tBS_Val*)AP->rP[1]; const int *cHists = (const int*)AP->rP[2]; const unsigned cGam = cP->rGam;
  const int cNx = cP->rNx, cNy = cP->rNy, cX0 = cP->rX0, cY0 = cP->rY0, cStride = cP->rStride, cFst = AP->rFst, cLst = AP->rLst;
  qDat += cStride*cFst; for (int qrow=cFst;qrow<cLst;qrow++) { // loop by rows
    int     qIOy[2]; mcCLAHEio<cB>( cNy, cY0, qrow, qIOy ); qIOy[0] *= cNx;
    for (int qcol=0;qcol<cStride;qcol++) { int qIOx[2]; mcCLAHEio<cB>( cNx, cX0, qcol, qIOx ); // loop by cols
      //int   qI = qIOy[0]+qIOx[0]; unsigned qV = qDat[qcol]; qV = mcRng( qV, 1, cGam-2 ); const int *cHistV = cHists+qV;
      //int   qV00 = abs(cHistV[qI*(int)cGam+1]-cHistV[qI*(int)cGam-1]),
      //      qV01 = abs(cHistV[(qI+1)*(int)cGam+1]-cHistV[(qI+1)*(int)cGam-1]),
      //      qV10 = abs(cHistV[(qI+cNx)*(int)cGam+1]-cHistV[(qI+cNx)*(int)cGam-1]),
      //      qV11 = abs(cHistV[(qI+cNx+1)*(int)cGam+1]-cHistV[(qI+cNx+1)*(int)cGam-1]);
      int   qI = qIOy[0]+qIOx[0]; unsigned qV = qDat[qcol]; if (qV>=cGam) qV = cGam-1; const int *cHistV = cHists+qV;
      int   qV00 = cHistV[qI*cGam], qV01 = cHistV[(qI+1)*cGam], qV10 = cHistV[(qI+cNx)*cGam], qV11 = cHistV[(qI+cNx+1)*cGam];
      const int cc = (1<<(cB-1))-1;
      qV00 += ((qV01-qV00)*qIOx[1]+cc)>>cB; qV10 += ((qV11-qV10)*qIOx[1]+cc)>>cB; qDat[qcol] = qV00+(((qV10-qV00)*qIOy[1]+cc)>>cB);
  }; qDat += cStride; }; } // for col,row
void BS_GamCLAHEApply( int AW, int AH, tBS_Val *ADat, int AB, int AGam, int *AHists ) { // -- apply Hists[Ny][Nx]
  if (AB<=0) { CBS_XLat qXLat; qXLat.Init( AHists, 0, AGam, 0, 0 ); qXLat.Apply( AW*AH, ADat ); return; }; // whole frame
  tCLAHE qP; if (!AHists || !iCLAHE( AW, AH, AB, AGam, &qP )) return; // get blocks
  tBS_SpawnThreadProc qProc; switch (AB) {
    case cBmn+0:qProc = ippr_CLAHEApply<(cBmn+0)>; break; case cBmn+1:qProc = ippr_CLAHEApply<(cBmn+1)>; break;
    case cBmn+2:qProc = ippr_CLAHEApply<(cBmn+2)>; break; case cBmn+3:qProc = ippr_CLAHEApply<(cBmn+3)>; break;
    case cBmn+4:qProc = ippr_CLAHEApply<(cBmn+4)>; break; case cBmn+5:qProc = ippr_CLAHEApply<(cBmn+5)>; break;
    default: return; }; BS_SpawnThreads( cBS_SpawnPool, 0, AH, qProc, &qP, ADat, AHists ); }

//============= moving averages =============
//===========================================
void BS_MovAve1D( int AL, int AR, const tBS_Val *AInp, tBS_Val *AOut ) { // ~mcMovAve<BS_Val,...> with sentinels
  if (AR<=-AL/2) { int qL = 0, qA = 0; for (int qi=0;qi<AL;qi++) { int qN, qV = AInp[qi]; //R=-9999: ave within lap's bounds
    if (qV>0) { if (qL<0) { qN = 1; goto lab1; }; qL++; qA += qV; continue; };
    if (qV<0) { if (qL>0) { qN =-1; goto lab1; }; qL--; qA += qV; continue; };
    AOut[qi] = 0; if (!qL) continue; qN = 0; lab1: qL = abs(qL);
    BS_VFil( qL, AOut+qi-qL, qA/qL ); qL = qN; qA = qV; }; // for i
    if (qL) { qL = abs(qL); BS_VFil( qL, AOut+AL-qL, qA/qL ); }; return; }; // 
  if (AR<0) { int qL = 0; for (int qi=0;qi<AL;qi++) { int qV = AInp[qi]; // R<0: ave within lap by fixed radius
    if (qV>0) { if (qL<0) { qV =-qL; qL = 1; goto lab2; }; qL++; continue; };
    if (qV<0) { if (qL>0) { qV = qL; qL =-1; goto lab2; }; qL--; continue; };
    AOut[qi] = 0; if (!qL) continue; qV = abs(qL); qL = 0; lab2:
    mcMovAve<tBS_Val>( qV, -AR, AInp+qi-qV, AOut+qi-qV ); }; // for i
    if (qL) { qL = abs(qL); mcMovAve<tBS_Val>( qL, -AR, AInp+AL-qL, AOut+AL-qL ); }; return; };
  if (AR<2 || AR>15 || AR>=AL/2) {
    mcMovAve<tBS_Val>( AL, AR, AInp, AOut );
    return; };
  switch (AR) { // normal moving average
    case  2: mcMovAve<tBS_Val, 2>( AL, AInp, AOut ); break; case  3: mcMovAve<tBS_Val, 3>( AL, AInp, AOut ); break;
    case  4: mcMovAve<tBS_Val, 4>( AL, AInp, AOut ); break; case  5: mcMovAve<tBS_Val, 5>( AL, AInp, AOut ); break;
    case  6: mcMovAve<tBS_Val, 6>( AL, AInp, AOut ); break; case  7: mcMovAve<tBS_Val, 7>( AL, AInp, AOut ); break;
    case  8: mcMovAve<tBS_Val, 8>( AL, AInp, AOut ); break; case  9: mcMovAve<tBS_Val, 9>( AL, AInp, AOut ); break;
    case 10: mcMovAve<tBS_Val,10>( AL, AInp, AOut ); break; case 11: mcMovAve<tBS_Val,11>( AL, AInp, AOut ); break;
    case 12: mcMovAve<tBS_Val,12>( AL, AInp, AOut ); break; case 13: mcMovAve<tBS_Val,13>( AL, AInp, AOut ); break;
    case 14: mcMovAve<tBS_Val,14>( AL, AInp, AOut ); break; case 15: mcMovAve<tBS_Val,15>( AL, AInp, AOut ); break;
  }; } // switch R
static void mcMovAve1D( int AMode, int AL, int AR, const tBS_Val *AInp, tBS_Val *AOut ) { BS_MovAve1D( AL, AR, AInp, AOut ); }
void iilpr_MovAvePostFix( int ALL, tBS_SpawnLIRec *ASPR, int *ABufI, tBS_Val *AData ) {
  switch (ALL) {
    case -1: { const int csz = ASPR->rLIR.rMaxLen*cBS_Sz; ASPR->rDat[0] = malloc(csz); ASPR->rDat[1] = malloc(csz); }; break;
    case -2: { free( ASPR->rDat[0] ); free( ASPR->rDat[1] ); }; break;
    case -4: if (ASPR->rInt[1]>0) { int qW = ASPR->rLIR.rW, qH = ASPR->rLIR.rH; // can't correct mix plus|minus when R<0
      switch (ASPR->rLIR.rDir) {
        case 2: case 8: { // nearly horz-dir needs horz-smooth to fix uneven vertical spacing
          tBS_Val   *qBuf0 = (tBS_Val*)ASPR->rDat[0];
          for (int qr=0;qr<qH;qr++) { memcpy( qBuf0, AData, qW*cBS_Sz );
            for (int qc=2;qc<qW-2;qc++) AData[qc] = (qBuf0[qc-1]+qBuf0[qc]+qBuf0[qc+1])/3; AData += qW; };
        }; break;
        case 4: case 6: { // nearly vert-dir needs vert-smooth to fix uneven horz spacing
          tBS_Val   *qDat3[3] = { (tBS_Val*)ASPR->rDat[0], (tBS_Val*)ASPR->rDat[1], (tBS_Val*)ABufI };
          memcpy( qDat3[0], AData, qW*cBS_Sz ); AData += qW;
          memcpy( qDat3[1], AData, qW*cBS_Sz );
          int        qr0 = 0; // internal index of row #0 (by modulo 3)
          for (int qr=2;qr<qH-1;qr++) {
            tBS_Val *qBuf0 = (tBS_Val*)qDat3[qr0],       // row #0 'current-1'
                    *qBuf1 = (tBS_Val*)qDat3[(qr0+1)%3], // row #1 (matches to current data-row)
                    *qBuf2 = (tBS_Val*)qDat3[(qr0+2)%3]; // row #2 'current+1'
            memcpy( qBuf2, AData+qW, qW*cBS_Sz );
            for (int qc=0;qc<qW;qc++) AData[qc] = (qBuf0[qc]+qBuf1[qc]+qBuf2[qc])/3;
            qr0 = (qr0+1)%3; AData += qW; }; // for r
        }; break; // if Int[1]>0
      }; // switch Dir
    }; break;
  }; } // switch LL

static void ilpr_MovAve( int ALL, int ADiri, tBS_SpawnLIRec *ASPR, int *ABufI, tBS_Val **AData ) {
  if (ALL<=0) {
    iilpr_MovAvePostFix( ALL, ASPR, ABufI, AData[2+ADiri] );
    if (ALL==-3) {
      int      qNP = ASPR->rLIR.rW*ASPR->rLIR.rH;
      tBS_Val *qFrm = AData[2+ADiri]; if (!qFrm) AData[2+ADiri] = qFrm = (tBS_Val*)malloc(qNP*cBS_Sz);
      if (AData[1]) BS_Cpy( qNP, AData[1], qFrm );
      int      qRad = mcVi(AData[0]), qMode = ASPR->rInt[0] = (qRad>>0)&0xFF; // mode
      ASPR->rInt[1] = BS_SDist( ASPR->rDir, qRad>>8 ); // radius
    }; return;
  }; // if LL<0
  tBS_Val     *qBuf1 = (tBS_Val*)ASPR->rDat[0], *qBuf2 = (tBS_Val*)ASPR->rDat[1], *qFrm = (tBS_Val*)AData[2+ADiri];
  BS_Gath( ALL, qFrm, qBuf1, ABufI );
  mcMovAve1D( ASPR->rInt[0], ALL, ASPR->rInt[1], qBuf1, qBuf2 ); 
  BS_Scat( ALL, qBuf2, qFrm, ABufI ); }
void BS_MovAve( int AR, int AW, int AH, tBS_Val **ADatas, const int *ADirs, int ANDD ) { // threads of 2D-movave 
  // if (AR<1) return; 
  const int AMode = 0; // reserved parameter
  tBS_Val *qDs[40]; qDs[0] = (tBS_Val*)mcIv((AR<<8)|AMode); qDs[1] = NULL;
  if (ANDD>mcASz(qDs)-2) ANDD = mcASz(qDs)-2; memcpy( qDs+2, ADatas, ANDD*sizeof(void*) );
  BS_SpawnLI( ANDD, ADirs, ilpr_MovAve, AW, AH, qDs ); }

static int ixMovAveXYRad; static int icplc_MovAveXY( int AL, tBS_Val *AData ) { tBS_Val *qBlr = AData+BS_Align(AL); mcMovAve1D( 0, AL, ixMovAveXYRad, AData, qBlr ); BS_Cpy( AL, qBlr, AData ); return AL; }
void BS_MovAveX( int AR, int AW, int AH, tBS_Val *AData ) { if (AR) { ixMovAveXYRad = AR; BS_CPLCRows( icplc_MovAveXY, AW, AH, AData ); }; } // horz average
void BS_MovAveY( int AR, int AW, int AH, tBS_Val *AData ) { if (AR) { ixMovAveXYRad = AR; BS_CPLCCols( icplc_MovAveXY, AW, AH, AData ); }; } // vert average

static int __fastcall icplb_MovAve2D2( tBS_Val *AData ) { return (AData[0]+AData[1])/2; }
static int __fastcall icplb_MovAve2D4( tBS_Val *AData ) { return (AData[0]+AData[1]+AData[2]+AData[3])/4; }
void BS_MovAve2D( int AR, int AW, int AH, tBS_Val *AData, unsigned AFlg ) { // engineering 2D-movave 
  const int   cN = AW*AH, AMode = 0; // reserved parameter
  tBS_Val    *qData[4]; qData[0] = AData; BS_DatasAlloc( 3, cN, qData+1 );
  if ((AFlg&4) && (AR>2)) { BS_DataDup( 3, cN, qData[0], qData+1 ); // 22 deg
    const int cDirs[4] = {2,4,6,8}; BS_MovAve( AR, AW, AH, qData, cDirs, 4 );
    BS_CPB( icplb_MovAve2D4, 4|cCPBNDD, cN, qData ); };
  if (AFlg&2) { BS_Cpy( cN, qData[0], qData[1] ); // 45 deg
    tBS_Val  *qDataT[6] = { (tBS_Val*)mcIv((AR<<8)|AMode), NULL, qData[0], qData[1], qData[0], qData[1] };
    const int cDirs[4] = {3|(1<<4),7|(1<<4),3|(2<<4),7|(2<<4)}; BS_SpawnLI( 4, cDirs, ilpr_MovAve, AW, AH, qDataT );
    BS_CPB( icplb_MovAve2D2, 2|cCPBNDD, cN, qData ); };
  if (AFlg&1) { BS_MovAveX( AR, AW, AH, AData ); BS_MovAveY( AR, AW, AH, AData ); }; // horz/vert
  BS_DatasRelease( 3, qData+1 ); }

//======= application spatial filters =======
//===========================================
template<int cR> int __fastcall mcAveNxN( const tBS_Val *AInp ) { const int cN = cR*2+1; return mcAve<tBS_Val,(cN*cN)>( AInp ); }
void BS_Ave3x3( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { BS_CPK( 1, mcAveNxN<1>, AW, AH, ASrc, ADst ); }
void BS_Ave5x5( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { BS_CPK( 2, mcAveNxN<2>, AW, AH, ASrc, ADst ); }

template<int cR> int __fastcall mcMedNxN( const tBS_Val *AInp ) { const int cN = (cR<<1)+1; tBS_Val qB[cN*cN];
  memcpy( qB, AInp, sizeof(qB) ); mcSortBAs<tBS_Val,cN*cN,(cN*cN>>1)+1>( qB ); return qB[cN*cN>>1]; }
template<int cR> int __fastcall mcXMedNxN( const tBS_Val *AInp ) { const int cN = (cR<<1)+1; tBS_Val qB[cN*cN]; // approx sort
  BS_Cpy( cN*cN, AInp, qB ); for (int qi1=0;qi1<cN;qi1++) mcSortBAs<tBS_Val,cN,(cN>>1)+1>( qB+cN*qi1 );
  for (int qi2=0;qi2<cN;qi2++) qB[qi2] = qB[qi2*cN+(cN>>1)]; mcSortBAs<tBS_Val,cN,(cN>>1)+1>( qB ); return qB[cN>>1]; }
void BS_Med( int AR, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { switch (AR) { // mcMedNxN<2>? mcMedNxN<3>?
  case 1: BS_CPK( 1, mcMedNxN <1>, AW, AH, ASrc, ADst ); break; case 2: BS_CPK( 2, mcMedNxN <2>, AW, AH, ASrc, ADst ); break;
  case 3: BS_CPK( 3, mcXMedNxN<3>, AW, AH, ASrc, ADst ); break; case 4: BS_CPK( 4, mcXMedNxN<4>, AW, AH, ASrc, ADst ); break;
  case 5: BS_CPK( 5, mcXMedNxN<5>, AW, AH, ASrc, ADst ); break; case 6: BS_CPK( 6, mcXMedNxN<6>, AW, AH, ASrc, ADst ); break;
  case 7: BS_CPK( 7, mcXMedNxN<7>, AW, AH, ASrc, ADst ); break; case 8: BS_CPK( 8, mcXMedNxN<8>, AW, AH, ASrc, ADst ); break;
  case 9: BS_CPK( 9, mcXMedNxN<9>, AW, AH, ASrc, ADst ); break; default: ErrMsg( "BS_Med params" ); }; }

// temp
template<int cR> int __fastcall mcMedBNxN( const tBS_Val *AInp ) { const int cN = (cR<<1)+1; tBS_Val qB[cN*cN];
  memcpy( qB, AInp, sizeof(qB) ); mcSortBAs<tBS_Val,cN*cN,cN*cN-1>( qB ); return qB[cN*cN-1]-qB[0]; }

void BS_MedB( int AR, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) {
  BS_CPK( 1, mcMedBNxN<1>, AW, AH, ASrc, ADst );
}

// absolute derivation Rad=1
static int icpl_DrA1( int AL, tBS_Val *ADat ) { int qV0 = ADat[0], qV1 = ADat[1]; ADat[0] = 0; AL--; // absolute derivation
  for (int qi=1;qi<AL;qi++) { int qV2 = ADat[qi+1]; ADat[qi] = abs(qV2-qV0); qV0 = qV1; qV1 = qV2; }; ADat[AL] = 0; return AL+1; }
static int icpl_Dr1( int AL, tBS_Val *ADat ) { int qV0 = ADat[0], qV1 = ADat[1]; ADat[0] = 0; AL--; // signed derivation
  for (int qi=1;qi<AL;qi++) { int qV2 = ADat[qi+1]; ADat[qi] = qV2-qV0; qV0 = qV1; qV1 = qV2; }; ADat[AL] = 0; return AL+1; }
void BS_Dr1( int AW, int AH, const tBS_Val *ADat, tBS_Val *ADrX, tBS_Val *ADrY, bool Asgn ) { const int cN = AW*AH; // -- XYders (R=1)
  if (ADrY) { BS_Cpy( cN, ADat, ADrY ); BS_CPLCCols( Asgn?icpl_Dr1:icpl_DrA1, AW, AH, ADrY ); }; 
  if (ADrX) { BS_Cpy( cN, ADat, ADrX ); BS_CPLCRows( Asgn?icpl_Dr1:icpl_DrA1, AW, AH, ADrX, true ); }; }
void BS_DrA3x3( int AW, int AH, tBS_Val *ADat ) { const int cN = AW*AH; tBS_Val *qDrY = new tBS_Val[cN]; // -- sum(abs(XYders))
  BS_Dr1( AW, AH, ADat, ADat, qDrY, false ); BS_Add( cN, qDrY, ADat ); delete qDrY; }
void BS_DrA3x3( int AW, int AH, tBS_Val *ADat, tBS_Val **ADrXY ) { const int cN = AW*AH;  // -- sum(abs(XYders)); keep the Ders
  BS_Dr1( AW, AH, ADat, ADrXY[0], ADrXY[1], false ); BS_Cpy( cN, ADrXY[0], ADat ); BS_Add( cN, ADrXY[1], ADat ); }

//======= Radon-Shannon-denoise =======
//=====================================
#if 1
template <int cThrLg, int cReg> static int __fastcall imccpb_Noise( tBS_Val *ADat8 ) { 
  // inp:eight values "Blr(Org,rad,dir]-Org"; out:"basic RadonShannon_noise".
  // "RadonDeviation~=0":max_noise (i.e subtract thick lap), ">=1 shl ThrLg":min_noise (i.e. subtract thin lap)
  int qD = 0; for (int qi2=0;qi2<4;qi2++) qD += abs(ADat8[qi2]-ADat8[qi2+4]);
  if (!qD) return ADat8[0]; int qN = 8-(qD>>cThrLg); if (qN==8) return mcSum<tBS_Val,8>( ADat8 )/(8+cReg);
  if (qN<1) qN = 1; int qZ = 0; for (int qi4=0;qi4<qN;qi4++) { int qMnJ = 0, qMn = ADat8[qi4];
    for (int qj=qi4+1;qj<8;qj++) { int qV = ADat8[qj]; if (abs(qMn)>abs(qV)) { qMn = qV; qMnJ = qj; }; }; // for j
    if (qMnJ) ADat8[qMnJ] = ADat8[qi4]; qZ += qMn; }; return qZ/(qN+cReg); } // for i4
#else
template <int cThrLg, int cReg> static int __fastcall imccpb_Noise( tBS_Val *ADat8 ) { 
  int qD = 0; for (int qi2=0;qi2<4;qi2++) qD += abs(ADat8[qi2]-ADat8[qi2+4]);
  if (!qD) return ADat8[0]; int qN = 8-(qD>>cThrLg); if (qN==8) return 8;
  if (qN<1) qN = 1; int qZ = 0; for (int qi4=0;qi4<qN;qi4++) { int qMnJ = 0, qMn = ADat8[qi4];
    for (int qj=qi4+1;qj<8;qj++) { int qV = ADat8[qj]; if (abs(qMn)>abs(qV)) { qMn = qV; qMnJ = qj; }; }; // for j
    if (qMnJ) ADat8[qMnJ] = ADat8[qi4]; qZ += qMn; };
  if (qN<1)
    qN = qN;
  return qN; } // for i4
#endif
#if 1 // subjective operation
  static int iNoiseLimboMn( tBS_Val *ADat ) {
    tBS_Val qDat[8]; BS_Cpy( 8, ADat, qDat );
    for (int qi=0;qi<8;qi++) {
      int qU = qDat[(qi-1)&7], qV = qDat[(qi+1)&7];
      if (abs(qU)>abs(qV)) qU = qV; qV = qDat[qi]; 
      if ((qV<0 && qU<qV-1) || (qV>0 && qU>qV+1)) qV = (qV+qU)>>1;
      ADat[qi] = qV; };
    return 8;
  }
#endif
void BS_Noise( int AR, int AW, int AH, // -- basic RadonShannon_noise
  const tBS_Val *AOrg, tBS_Val *ANoise, tBS_Val **ATs8, int AThrLg ) {
  static const tBS_CPBProc cProcs[4][12] = {
    {imccpb_Noise< 5,0>,imccpb_Noise< 6,0>,imccpb_Noise< 7,0>,imccpb_Noise< 8,0>,imccpb_Noise< 9,0>,imccpb_Noise<10,0>,
     imccpb_Noise<11,0>,imccpb_Noise<12,0>,imccpb_Noise<13,0>,imccpb_Noise<14,0>,imccpb_Noise<15,0>,imccpb_Noise<16,0>},
    {imccpb_Noise< 5,1>,imccpb_Noise< 6,1>,imccpb_Noise< 7,1>,imccpb_Noise< 8,1>,imccpb_Noise< 9,1>,imccpb_Noise<10,1>,
     imccpb_Noise<11,1>,imccpb_Noise<12,1>,imccpb_Noise<13,1>,imccpb_Noise<14,1>,imccpb_Noise<15,1>,imccpb_Noise<16,1>},
    {imccpb_Noise< 5,2>,imccpb_Noise< 6,2>,imccpb_Noise< 7,2>,imccpb_Noise< 8,2>,imccpb_Noise< 9,2>,imccpb_Noise<10,2>,
     imccpb_Noise<11,2>,imccpb_Noise<12,2>,imccpb_Noise<13,2>,imccpb_Noise<14,2>,imccpb_Noise<15,2>,imccpb_Noise<16,2>},
    {imccpb_Noise< 5,3>,imccpb_Noise< 6,3>,imccpb_Noise< 7,3>,imccpb_Noise< 8,3>,imccpb_Noise< 9,3>,imccpb_Noise<10,3>,
     imccpb_Noise<11,3>,imccpb_Noise<12,3>,imccpb_Noise<13,3>,imccpb_Noise<14,3>,imccpb_Noise<15,3>,imccpb_Noise<16,3>} };
  int cReg = AThrLg/100; AThrLg %= 100; AThrLg -= 5; if (cReg>=mcASz(cProcs) || AThrLg<0 || AThrLg>11) {
    ErrMsg( "BS_Noise (threshold's log2)" ); return; }; 
  const int cN = AW*AH; BS_DataDup( 8, cN, AOrg, ATs8 ); BS_MovAve( AR, AW, AH, ATs8 );
  #if 1 // subjective operation
    BS_CPB( iNoiseLimboMn, 8, cN, ATs8 );
  #endif
  for (int qi=0;qi<8;qi++) BS_SubA( cN, AOrg, ATs8[qi] );
  BS_CPB( cProcs[cReg][AThrLg], 8|cCPBNDD, cN, ATs8, ANoise ); }

void BS_MxLapRngZ( int AZ, int AN, tBS_Val *ALap ) { if (AZ<=0) return; // -- "if |Lap[i]|<Z then Lap[i]>>=2"
  for (int qi=0;qi<AN;qi++) { int qV = ALap[qi], qA = abs(qV); if (qA<AZ) ALap[qi] = (qV<0)?-(qA>>2):qA>>2; }; }

static inline int iMx4( const tBS_Val *A ) { int M = A[0]; for (int i=1;i<4;i++) { int v = A[i<<1]; if (abs(M)<abs(v)) M = v; }; return M; }
static int __fastcall mccpb_MxLapStd8( tBS_Val *ADat ) { return (iMx4(ADat)+iMx4(ADat+1))>>1; }
int BS_MxLapStd( // -- reduce shot_noise after RadonShannon_noise
  int AR, int AW, int AH, tBS_Val *ALap, tBS_Val **ATs, float AK, int AZ0, int AZ1, float AShear ) { 
  // call-sequence: "_MxLapRngZ(Z0); MaxLap; _MxLapRngZ(Z1)"
  if (AR<1) return 0; const int cN = AW*AH; BS_MxLapRngZ( AZ0, cN, ALap );
  const float cEps = 1e-3f; const bool cJoin = (AK>cEps)&&(AK<1-cEps);
  if (cJoin) BS_Cpy( cN, ALap, ATs[9] ); int cNDD; switch (AR) {
  case  1: { const int cDir4[4] = {1,3,5,7}; BS_DataDup( 4, cN, ALap, ATs ); BS_MovAve( 1, AW, AH, ATs, cDir4, 4 );
    BS_CPB( mccpb_MxA<4>, 4|cCPBNDD, cN, ATs, ALap ); BS_MxLapRngZ( BS_MxLapZ1dflt(AZ0,AZ1), cN, ALap );
    if (AK<cEps) { BS_Ave3x3( AW, AH, ALap, ATs[4] ); BS_Add( cN, ATs[4], ALap ); mcShiftRgt( cN, ALap, 1 ); };
    cNDD = 4; }; break;
  default: {
    BS2A_ShearCPU( 0, AShear, ALap, ATs[0] );
    BS_DataDup( 8, cN, ALap, ATs ); BS_MovAve( AR, AW, AH, ATs );
    BS_CPB( AR>2?mccpb_MxLapStd8:mccpb_MxA<8>, 8|cCPBNDD, cN, ATs, ALap );
    BS2A_ShearCPU( 0, -AShear, ALap, ATs[0] );
    BS_MxLapRngZ( BS_MxLapZ1dflt(AZ0,AZ1), cN, ALap );
    if (AK<cEps) { BS_Ave3x3( AW, AH, ALap, ATs[8] ); BS_Cpy( cN, ATs[8], ALap ); };
    cNDD = 8; }; }; // switch R
  if (AK<cEps) return cNDD; BS_VMul( cN, ALap, AK );
  if (cJoin) { BS_VMul( cN, ATs[9], 1-AK ); BS_Add( cN, ATs[9], ALap ); };
  return cNDD; } // returns #NDD of dir averages in TS^

//=========== down-/up-sample 2x2 ============
//============================================
static void itppr_DimDn( tBS_SpawnThreadRec *AP ) { // downsample: 2W*2H => W*H
  const int cW = mcVi(AP->rP[0]), cLst = AP->rLst; const tBS_Val *cSrc = (const tBS_Val*)AP->rP[1];
  tBS_Val *qDst = (tBS_Val*)AP->rP[2]; for (int qr=AP->rFst;qr<cLst;qr++) { const int co = qr*cW; for (int qc=0;qc<cW;qc++) {
    const int c2 = co*4+qc*2; qDst[co+qc] = (cSrc[c2]+cSrc[c2+1]+cSrc[c2+cW*2]+cSrc[c2+cW*2+1])>>2; }; }; }
static void itppr_DimDnM( tBS_SpawnThreadRec *AP ) { // downsample: 2W*2H => W*H
  const int cW = mcVi(AP->rP[0]), cLst = AP->rLst; const tBS_Val *cSrc = (const tBS_Val*)AP->rP[1];
  tBS_Val *qDst = (tBS_Val*)AP->rP[2]; for (int qr=AP->rFst;qr<cLst;qr++) { const int co = qr*cW; for (int qc=0;qc<cW;qc++) {
    const int c2 = co*4+qc*2; int qM = cSrc[c2], qV = cSrc[c2+1]; if (abs(qM)<abs(qV)) qM = qV; qV = cSrc[c2+cW*2];
    if (abs(qM)<abs(qV)) qM = qV; qV = cSrc[c2+cW*2+1]; if (abs(qM)<abs(qV)) qM = qV; qDst[co+qc] = qM; }; }; }
void BS_DimDn( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst, bool AM ) { BS_SpawnThreads( 2, 0, AH, AM?itppr_DimDnM:itppr_DimDn, mcIv(AW), ASrc, ADst ); }

static void itppr_DimUp( tBS_SpawnThreadRec *AP ) { // upsample: W*H => 2W*2H
  const int cW = mcVi(AP->rP[0]), cLst = AP->rLst; const tBS_Val *cSrc = (const tBS_Val*)AP->rP[1];
  tBS_Val *qDst = (tBS_Val*)AP->rP[2]; for (int qr=AP->rFst;qr<cLst;qr++) { const int co = qr*cW; for (int qc=0;qc<cW;qc++) {
    const int c2 = co*4+qc*2; qDst[c2] = qDst[c2+1] = qDst[c2+cW*2] = qDst[c2+cW*2+1] = cSrc[co+qc]; }; }; }
void BS_DimUp( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { BS_SpawnThreads( cBS_SpawnPool|2, 0, AH, itppr_DimUp, mcIv(AW), ASrc, ADst ); }

static void itppr_DimUpC( tBS_SpawnThreadRec *AP ) { // correction after upsample: W*H => W*H
  const int cW = mcVi(AP->rP[0]), cFst = AP->rFst, cBeg = cFst*cW, cLst = AP->rLst, cW1 = cW-1;
  const tBS_Val *cSrc = ((const tBS_Val*)AP->rP[1])+cBeg; tBS_Val *qDst = ((tBS_Val*)AP->rP[2])+cBeg;
  for (int qr=cFst;qr<cLst;qr++) { int qy = (((qr&1)<<1)-1)*cW; *qDst++ = *cSrc++;
  for (int qc=1;qc<cW1;qc++) { int qx = ((qc&1)<<1)-1, qV = cSrc[qx]-cSrc[-qx]+cSrc[qy]-cSrc[-qy]; *qDst++ = (qV>>2)+*cSrc++; }; *qDst++ = *cSrc++; }; }
void BS_DimUpC( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ) { 
  BS_SpawnThreads( cBS_SpawnPool|2, 1, AH-1, itppr_DimUpC, mcIv(AW), ASrc, ADst );
  BS_Cpy( AW, ASrc, ADst ); BS_Cpy( AW, ASrc+AW*(AH-1), ADst+AW*(AH-1) ); }

tBS_Val *BS_ScaleDnA( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM ) { for (int qst=1;qst<=Apyr;qst++) { // pyramid's dn-scale
  Aw0 >>= 1; Ah0 >>= 1; BS_DimDn( Aw0, Ah0, ADat, ATmp, AM ); mcSwap( tBS_Val*, ADat, ATmp ); }; return ADat; }
void BS_ScaleDn( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM ) { BS_Cpy( (Aw0>>Apyr)*(Ah0>>Apyr), BS_ScaleDnA(Apyr,Aw0,Ah0,ADat,ATmp,AM), ADat ); }
tBS_Val *BS_ScaleUpA( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC ) { for (int qst=Apyr;qst>0;qst--) { // pyramid's up-scale
  BS_DimUp( Awp, Ahp, ADat, ATmp ); Awp <<= 1; Ahp <<= 1; if (AC) BS_DimUpC( Awp, Ahp, ATmp, ADat ); else mcSwap( tBS_Val*, ADat, ATmp ); }; return ADat; }
void BS_ScaleUp( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC ) { BS_Cpy( (Awp<<Apyr)*(Ahp<<Apyr), BS_ScaleUpA(Apyr,Awp,Ahp,ADat,ATmp,AC), ADat ); }

//============== for margins ================
//===========================================
static void itppr_Bevels( tBS_SpawnThreadRec *AP ) {
  int      qW = mcVi(AP->rP[0]), qH = qW>>16; qW &= 0x7FFF;
  tBS_Val *qData = (tBS_Val*)AP->rP[1];
  int      qM0x = mcVi(AP->rP[2]), qM1x = qM0x>>16; qM0x &= 0x7FFF;
  int      qM0y = mcVi(AP->rP[3]), qM1y = qM0y>>16; qM0y &= 0x7FFF;
  int      qStd = 0; switch (AP->rSpawn) {
    case 0: qStd =-qW; qData += (qM1y-1)*qW;     break; // top
    case 1: qStd = qW; qData += (qH-qM1y)*qW;    break; // btm
    case 2: qStd = -1; qData += qM1y*qW+qM1x-1;  break; // lft
    case 3: qStd =  1; qData += qM1y*qW+qW-qM1x; break; // rgt
  }; // switch
  int      qL = abs(qStd)>1?qW:qH-qM1y*2; tBS_Val *qB0 = new tBS_Val[qL], *qB1 = new tBS_Val[qL];
  if (abs(qStd)>1) for (int qr=1;qr<=qM1y;qr++) { // top, btm
      tBS_Val *qDat = qData; if (qr>qM1y-qM0y) qDat -= qStd; BS_Cpy( qL, qDat, qB0 );
      if (qM1x>0) { BS_VFil( qM1x, qB0, qB0[qM1x] ); BS_VFil( qM1x, qB0+qW-qM1x, qB0[qW-qM1x-1] ); };
      mcMovAve<tBS_Val>( qL, qr, qB0, qB1 ); BS_Cpy( qL, qB1, qData ); qData += qStd;
    } // for r
  else for (int qr=1;qr<=qM1x;qr++) { // lft, rgt
      tBS_Val *qDat = qData; if (qr>qM1x-qM0x) qDat -= qStd; mcGath<tBS_Val>( qL, qDat, qB0, qW );
      mcMovAve<tBS_Val>( qL, qr, qB0, qB1 ); mcScat<tBS_Val>( qL, qB1, qData, qW ); qData += qStd;
    }; // for r
  delete qB0; delete qB1; }
void BS_Bevels( int AW, int AH, tBS_Val *ADat, int AM0x, int AM1x, int AM0y, int AM1y ) { if (AM0y<0) AM0y = AM0x; if (AM1y<0) AM1y = AM1x;
  BS_SpawnThreads( cBS_SpawnPool|4, 0, 4, itppr_Bevels, mcIv((AH<<16)|AW), ADat, mcIv((AM1x<<16)|AM0x), mcIv((AM1y<<16)|AM0y) ); }

static void iVFilLft( int AL, tBS_Val *ADst, int AV ) { if (AV<0) AV = ADst[ AL]; while (--AL>=0) *ADst++ = AV; } // Dat:=V lft-to-rgt
static void iVFilRgt( int AL, tBS_Val *ADst, int AV ) { if (AV<0) AV = ADst[-AL]; while (--AL>=0) *ADst-- = AV; } // Dat:=V rgt-to-lft
static void iVAddTop( int AL, int AW, tBS_Val *ADst ) { int qV = ADst[AL*AW]; while (--AL>=0) { *ADst = (*ADst+qV)>>1; ADst += AW; }; } // Dat:=.../2 top-to-btm
static void iVAddBtm( int AL, int AW, tBS_Val *ADst ) { int qV = ADst[-AL*AW]; while (--AL>=0) { *ADst = (*ADst+qV)>>1; ADst -= AW; }; } // Dat:=.../2 btm-to-top
static void ippr_FillBev( tBS_SpawnThreadRec* AP ) {
  int qW = mcVi(AP->rP[0]), qH = qW>>16; qW &= 0xFFFF; tBS_Val *qFrame = (tBS_Val*)AP->rP[1];
  int qBevel = mcVi(AP->rP[2]), qV = (short)(qBevel>>16); qBevel &= 0xFFFF;
  for (int qCorner=AP->rFst;qCorner<AP->rLst;qCorner++) switch (qCorner) {
    case 0: { for (int qi=0;qi<qBevel;qi++) iVFilLft( qBevel-qi, qFrame+(qH-1-qi)*qW, qV ); // lft,btm
      if (qV<0) for (int qi=0;qi<qBevel;qi++) iVAddBtm( qBevel-qi, qW, qFrame+(qH-1)*qW+qi ); }; break;
    case 1: { for (int qi=0;qi<qBevel;qi++) iVFilRgt( qBevel-qi, qFrame+(qH-1-qi)*qW+qW-1, qV ); // rgt,btm
      if (qV<0) for (int qi=0;qi<qBevel;qi++) iVAddBtm( qBevel-qi, qW, qFrame+(qH-1)*qW+qW-1-qi ); }; break;
    case 2: { for (int qi=0;qi<qBevel;qi++) iVFilLft( qBevel-qi, qFrame+qi*qW, qV ); // lft,top
      if (qV<0) for (int qi=0;qi<qBevel;qi++) iVAddTop( qBevel-qi, qW, qFrame+qi ); }; break;
    case 3: { for (int qi=0;qi<qBevel;qi++) iVFilRgt( qBevel-qi, qFrame+qi*qW+qW-1, qV ); // rgt,top
      if (qV<0) for (int qi=0;qi<qBevel;qi++) iVAddTop( qBevel-qi, qW, qFrame+(qW-1)-qi ); }; break;
  }; } // switch
void BS_BevCorners( int AW, int AH, tBS_Val *AFrame, int ABevCorners, int ABevel, int AV ) {
  //const int cBevelCorners = ABevel>>16; if (!cBevelCorners) return; ABevel &= 0xFFFF;
  if (!ABevCorners || ABevel<3) return;
  BS_SpawnThreads( cBS_SpawnPool|ABevCorners, 0, ABevCorners, ippr_FillBev, mcIv((AH<<16)|AW), AFrame, mcIv((AV<<16)|ABevel) ); }

void BS_DrawRect( int AW, int AH, tBS_Val *AData, const int ARect[4], int AMark ) { // -- draw rectangle
  const int cX0 = mcRng(ARect[0],0,AW-1), cY0 = mcRng(ARect[1],0,AH-1), cX1 = mcRng(ARect[2],cX0+1,AW),
    cY1 = mcRng(ARect[3],cY0+1,AH), cW = cX1-cX0, cH = cY1-cY0; AData += cX0+cY0*AW; 
  BS_VScat( cH, AData, AW, AMark ); BS_VScat( cH, AData+cW-1, AW, AMark );
  BS_VFil( cW, AData, AMark ); BS_VFil( cW, AData+AW*(cH-1), AMark ); }
