// I16_Ut7.cpp : ad-hoc service CPU-transforms (for collimator, ribcage etc)
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut7.h"

//================= commons ==================
//============================================
static int ixW, ixH, ixDbgDmp; // standard frame-size for local operations
void BS7_Dim( int AW, int AH, int ADbgDmp ) { ixW = AW; ixH = AH; ixDbgDmp = ADbgDmp; } // -- assigns standard frame's dimension

//======= alter dirlayer's projection ========
//============================================
// fw_ & bk_projection for parallel beams; centinel dbg_prints in callbacks against scramble by multithreading 
static void ppr_Projection( tBS_SpawnThreadRec *AP ) { // gather & scatter line_start pixels of 2D-array to/from 1D-buffer
  // inp: dirlayers with a metric in line_start pixels; out:altered metrics in line_start pixels
  int  qW = mcVi(AP->rP[1]), qH = qW>>16; qW &= 0x7FFF;  const int cL = qW+qH, *cDirs = (const int*)AP->rP[2];
  int *qBufI = new int[cL]; tBS_Val *qPrj = new tBS_Val[cL], **qDatas = (tBS_Val**)AP->rP[3]; 
  for (int qdiri=AP->rFst;qdiri<AP->rLst;qdiri++) { int qL = BS_GetProjectionBufI( cDirs[qdiri], qW, qH, qBufI ); BS_Gath( qL, qDatas[qdiri], qPrj, qBufI );
    if (((tBS7_ProjectionProc)AP->rP[0])( qdiri, qL, qPrj )) BS_Scat( qL, qPrj, qDatas[qdiri], qBufI ); }; delete qPrj; delete qBufI; } // for i
static int cpl_Projection( int AL, tBS_Val *ADat ) { BS_VFil( AL-1, ADat+1, ADat[0] ); return AL; }; // apply filtered ortho-metric
void BS7_Projection( tBS_CPLProc ACPLProcInp, tBS7_ProjectionProc AProjectionProc, int ANDD, const int *ADirs, tBS_Val **ADatas, const int *AWH ) { // -- wrapper for '_ProjectionProc'
  const int cW = AWH?AWH[0]:ixW, cH = AWH?AWH[1]:ixH; if (ACPLProcInp) BS_CPLA( ACPLProcInp, ANDD, ADirs, cW, cH, ADatas ); // ~calc ortho-metric
  BS_SpawnThreads( ANDD, 0, ANDD, ppr_Projection, (const void*)AProjectionProc, mcIv((cH<<16)|cW), ADirs, ADatas ); // ~filter projected ortho-metric
  BS_CPLA( cpl_Projection, ANDD, ADirs, cW, cH, ADatas ); } // ~apply filtered ortho-metric

//=== shear/warp by specified skew_period ====
//============================================
void BS7_SkewH( int ASkew, bool AFill, tBS_Val *AData, int AW, int AH ) { if (!AW) { AW = ixW; AH = ixH; }; // -- hrz skew
  int qShear = BS7_SkewShear(ASkew,AH); if (!qShear) return;
  ASkew = abs(ASkew); int qInc = qShear<0?1:-1, qRowsRem = AH; do { // loop by rows
    int qRowsCpy = ASkew; if (qRowsCpy>qRowsRem) qRowsCpy = qRowsRem; int qCpy = qRowsCpy*AW;
    if (qShear>0) { memmove( AData+qShear, AData, (qCpy-qShear)*cBS_Sz ); // positive shear
      if (AFill) { tBS_Val *qDat = AData+qShear; for (int qr=0;qr<qRowsCpy;qr++) { BS_VFil( qShear, qDat-qShear, qDat[0] ); qDat += AW; }; };
    } else if (qShear<0) { memmove( AData, AData-qShear, (qCpy+qShear)*cBS_Sz ); // negative shear
      if (AFill) { tBS_Val *qDat = AData+AW+qShear; for (int qr=0;qr<qRowsCpy;qr++) { BS_VFil( -qShear, qDat, qDat[-1] ); qDat += AW; }; }; };
    qShear += qInc; AData += qCpy; qRowsRem -= qRowsCpy; } while (qRowsRem); } // of loop by rows

void BS7_SkewV( int ASkew, bool AFill, tBS_Val *AData, int AW, int AH ) { if (!AW) { AW = ixW; AH = ixH; }; // -- vrt skew
  int qShear = BS7_SkewShear(-ASkew,AW); if (!qShear) return;
  ASkew = abs(ASkew); int qInc = qShear<0?1:-1, qColsRem = AW; do { // loop by cols
    int qColsCpy = ASkew; if (qColsCpy>qColsRem) qColsCpy = qColsRem;
    if (qShear>0) { tBS_Val *qDat = AData+(AH-qShear-1)*AW; // positive shear
      do { BS_Cpy( qColsCpy, qDat, qDat+qShear*AW ); qDat -= AW; } while (qDat>=AData);
      if (AFill) for (int qr=1;qr<qShear;qr++) BS_Cpy( qColsCpy, AData, AData+qr*AW );
    } else if (qShear<0) { tBS_Val *qDat = AData-qShear*AW, *qDatE = AData+(AH-1)*AW; // negative shear
      do { BS_Cpy( qColsCpy, qDat, qDat+qShear*AW ); qDat += AW; } while (qDat<=qDatE);
      if (AFill) for (int qr=-1;qr>qShear;qr--) BS_Cpy( qColsCpy, qDatE, qDatE+qr*AW ); };
    qShear += qInc; AData += qColsCpy; qColsRem -= qColsCpy; } while (qColsRem); } // of loop by cols

tBS_Val *BS7_SkewAlloc( int ASkewH, int ASkewV, const tBS_Val *AData, int *AWH, int ATag ) { const int cW = ixW, cH = ixH; // -- extend & skew image
  int      qShearH = BS7_SkewShear(abs(ASkewH),cH), qShearV = BS7_SkewShear(abs(ASkewV),cW); //qShearH = qShearV = __max(qShearH,qShearV)*2;
  int      qW = cW; if (qShearH) qW = BS_Align( qW+qShearH*2+1 ); AWH[0] = qW; 
  int      qH = cH; if (qShearV) qH = BS_Align( qH+qShearV*2+1 ); AWH[1] = qH; 
  tBS_Val *qData = new tBS_Val[qW*qH]; BS_VFil( qW*qH, qData, ATag&0xFFFF ); if (AData) {
    const int cGapH0 = (qW-cW)/2, cGapV0 = (qH-cH)/2; BS_RowsCpy( cW, AData, qW, qData+cGapH0+cGapV0*qW, cW, cH ); 
    if (ATag & cBS7_SkewAllocMarg) { // complement rectangular margins
      if (qShearH) { tBS_Val *qP0 = qData+cGapV0*qW; for (int qi0=0;qi0<cH;qi0++) { BS_VFil( cGapH0, qP0, qP0[cGapH0] ); qP0 += qW; };
        const int cGapH1 = qW-cW-cGapH0; tBS_Val *qP1 = qData+cGapV0*qW+qW-cGapH1; for (int qi1=0;qi1<cH;qi1++) { BS_VFil( cGapH1, qP1, qP1[-1] ); qP1 += qW; }; }; // if ShearH
      if (qShearV) { tBS_Val *qP0 = qData+cGapV0*qW; for (int qi0=1;qi0<=cGapV0;qi0++) BS_Cpy( qW, qP0, qP0-qi0*qW );
        const int cGapV1 = qH-cH-cGapV0; tBS_Val *qP1 = qData+(qH-cGapV1-1)*qW; for (int qi1=1;qi1<=cGapV1;qi1++) BS_Cpy( qW, qP1, qP1+qi1*qW ); }; // if ShearV
    };
    if (ATag & cBS7_SkewAllocBeve) { // complement blurred bevels
      // typically use when ShearH~=qShearV; specific addition (~>3pix) depends on implementation of _Bevels
      const int cGap = __max(cGapH0,cGapV0)+4; BS_Bevels( qW, qH, qData, cGap, cGap+1 ); };
    const bool cFill = (ATag&cBS7_SkewAllocFill)!=0; // fill in _Skew?
    BS7_SkewH( ASkewH, cFill, qData, AWH[0], AWH[1] ); BS7_SkewV( ASkewV, cFill, qData, AWH[0], AWH[1] );
  } else BS_ZFil( qW*qH, qData ); return qData; }

void BS7_SkewRelease( int ASkewH, int ASkewV, tBS_Val *AData, const int *AWH, tBS_Val *ADataSkew ) { if (AData) { // -- release skewed image
  BS7_SkewV( -ASkewV, false, ADataSkew, AWH[0], AWH[1] ); BS7_SkewH( -ASkewH, false, ADataSkew, AWH[0], AWH[1] ); // copy bk DataSkew to Data
  BS_RowsCpy( AWH[0], ADataSkew+(AWH[0]-ixW)/2+(AWH[1]-ixH)/2*AWH[0], ixW, AData, ixW, ixH ); }; delete ADataSkew; }

void BS7_SkewFillMgn( bool AforSkewV, tBS_Val *AData, int AV, int AW, int AH ) { if (!AW) { AW = ixW; AH = ixH; }; // -- hrz skew
  if (AforSkewV) { BS_VFil( AW, AData, AV ); BS_VFil( AW, AData+(AH-1)*AW, AV ); }
            else { BS_VScat( AH, AData, AW, AV ); BS_VScat( AH, AData+AW-1, AW, AV ); }; }
