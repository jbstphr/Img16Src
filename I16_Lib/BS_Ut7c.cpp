// I16_Ut7c.cpp : ribs-detection's services
// by S.B.

#include "stdafx.h"
#include "math.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut7.h"
#include "BS_Ut7c.h"

// Standard entry code for exporting procedures in this module
// Define 'N', 'DbgSfx' in ~stack_frame; redefine 'DbgPfx<>nil'.
#define mcExportProcEntry const int cN = AW*AH; char *lDbgSfx = NULL, qDbgPfx[88];\
  if (ADbgPfx) { strcpy( qDbgPfx, ADbgPfx ); ADbgPfx = qDbgPfx; lDbgSfx = strend(qDbgPfx); }

//======= Polyfit and associated funcs =======
//============================================

// polyfit
struct tLagrange1{ tLL rXi, rYi, rX2i, rXYi; int rN; // LaGrange polynomial of 1st order
  inline void Add( int x, int y ) { rXi += x; rYi += y; rX2i += x*x; rXYi += x*y; rN++; }
  inline float A1() { tLL qDen = rN*rX2i-rXi*rXi; return qDen ? (float)(rN*rXYi-rXi*rYi)/qDen : 0; }
  inline float A0( float A1 ) { return (rYi-A1*rXi)/rN; } };
int BS7c_Polyfit1( int AN, const tBS_Val *AX, const tBS_Val *AY, float *AA ) { // -- returns A0,A1 so that "AYi ~= AXi*A1+A0"
  struct tLagrange1 qLaG; memset( &qLaG, 0, sizeof(qLaG) ); for (int qi=0;qi<AN;qi++) qLaG.Add( AX[qi], AY[qi] );
  AA[1] = qLaG.A1(); AA[0] = qLaG.A0( AA[1] ); return (int)(qLaG.rXi/AN); } // return average

int BS7c_AscendRowPf1( int AN, tBS_Val *AX, tBS_Val *AY, float *AA, int AGap ) { // -- exclude XYpoints, which decrease A1
  float qA[2]; int qXC = BS7c_Polyfit1( AN, AX, AY, qA ), qN = 0;
  for (int qj=0;qj<AN;qj++) { 
    int qX = AX[qj], qY = AY[qj], qDY = qY-(int)(qX*qA[1]+qA[0]), qDX = qX-qXC;
    if (qDX<AGap) { if (qDY> AGap) continue; } else
    if (qDX>AGap) { if (qDY<-AGap) continue; } else 
                  { if (abs(qDY)>AGap) continue; };
    AX[qN] = qX; AY[qN] = qY; qN++; }; // for n
  if (AA) memcpy( AA, qA, sizeof(qA) ); return qN; }

// -- Polyfit associates
//static const tBS_Val cAscendAbsent = (tBS_Val)0x8000;
static const tBS_Val cAscendAbsent = (tBS_Val)0;
static const int cAscendArc    = mcASz(cBS7_AscendDirs),
                 cAscendMag    = 100,
                 cAscendGainLg = 9;
void BS7c_AscArcRowLimbo( int AN, int ACX, tBS_Val *ADat, float A0, float A1 ) { // -- fill row by limbo_metric (by A0,A1; center X at 'CX')
  const int cA0 = (int)(A0+0.49999f), cA1 = (int)(A1*(1<<10)); for (int qi=0;qi<AN;qi++) ADat[qi] = ((qi-ACX)*cA1>>10)+cA0; }

static float inline mcAscArcRowInc( int AV ) {
  #ifdef BS7_ASCENDDIRLESS
    const float cArcIncs[cAscendArc] = {-1,-0.33f, 0, 0.33f, 1, 3 };
  #else
    const float cArcIncs[cAscendArc] = {-3,-1,-0.33f, 0, 0.33f, 1, 3 };
  #endif
  int qPie = AV/cAscendMag; if (qPie>=cAscendArc-1) return cArcIncs[cAscendArc-1];
  return mcIpl( AV-qPie*cAscendMag, cAscendMag, cArcIncs[qPie], cArcIncs[qPie+1] ); }
int BS7c_AscArcRowSkew( int AN, int ACX, tBS_Val *ADat, float *ATmp ) { // -- row's skew by Y
  float *qTmp = ATmp?ATmp:new float[AN];
  float qSkew = 0; for (int qi0=0;qi0<AN;qi0++) { qSkew += mcAscArcRowInc( ADat[qi0] ); qTmp[qi0] = qSkew; };
  qSkew = qTmp[ACX]; for (int qi1=0;qi1<AN;qi1++) ADat[qi1] = (int)(qTmp[qi1]-qSkew);
  if (qTmp!=ATmp) delete qTmp; return __max(ADat[0],ADat[AN-1]); }

// prepare inp of AscendFrmPf
static int iCPB_ArcToMetAscend( tBS_Val *ADat ) {
  #ifdef BS7_ASCENDDIRLESS
    const int cIdx[]= {2,3,1,4,0,5};
  #else
    const int cIdx[] = {3,2,4,1,5,0,6};
  #endif
  int qMi = cIdx[0], qMn = ADat[qMi], qMx = qMn;
  for (int qi=1;qi<cAscendArc;qi++) { int idx = cIdx[qi];
    int v = ADat[idx]; if (qMn>v) { qMn = v; qMi = idx; } else if (qMx<v) qMx = v; };
  return qMx==qMn ? cAscendAbsent : qMi*cAscendMag; }
void BS7c_ArcToMet( const tBS_Val *AMsk, tBS_Val **AArc, tBS_Val *AMet, int AW, int AH ) { // -- prepare inp of AscendFrmPf
  const int cN = AW*AH; BS_ZFil( cN, AMet );
  BS_CPB( iCPB_ArcToMetAscend, cCPBNDD|7, cN, AArc, AMet );
  if (AMsk) for (int qi=0;qi<cN;qi++) if (AMsk[qi]<9) AMet[qi] = (tBS_Val)cAscendAbsent; }

// ascending polyfit of frame rows
static bool ixAscArcFrmPfDbg; static int icpl_AscArcFrmPf( int AL, tBS_Val *ADat ) { // polyfit a row
  const int cA = BS_Align( AL );
  tBS_Val  *qBufX = ADat+cA, *qBufY = qBufX+cA;
  int       qN = 0, qJ0 = 0, qJ1;
  for (int qi=0;qi<AL;qi++) { int qV = ADat[qi]; if (qV!=cAscendAbsent) {
    qBufX[qN] = qi-(AL>>1); qBufY[qN] = qV; qN++;
    if (!qJ0) qJ0 = qi; qJ1 = qi; }; };
  ADat[0] = ADat[1] = -1; // initialize output
  if (qN<5 || qJ1-qJ0<60) return 2; // don't analyze rows with too short mask's segment
  float     cGainScale = (float)(cAscendArc*cAscendMag)/qJ1,
            cA1mn = 0.5f*cGainScale, cA1mx = 1.0f*cGainScale;
  // calc polyfit' coeffs for a row
  int    cNbk = qN;
  float  qAbk[2]; qAbk[1] = 0;
  for (int qiter=0;qiter<5;qiter++) {
    float *qApt = NULL, qA[2]; int qNN = BS7c_AscendRowPf1( qN, qBufX, qBufY, qA, 8 );
    const bool cA1ok = (qA[1]>cA1mn) && (qA[1]<cA1mx);
    if (cA1ok && (qN-qNN<cNbk/9+9)) qApt = qA;
    // else if ((qNN<cNbk/3) && (qAbk[1]>1e-3f)) qApt = qAbk;
    if (qApt) { int qV;
      qV = (int)(qApt[0]+0.4999f); ADat[0] = mcRng( qV, -0x7FFF, 0x7FFF );
      qV = (int)(qApt[1]*(1<<cAscendGainLg)); ADat[1] = mcRng( qV, -0x7FFF, 0x7FFF );
      break; };
    if (cA1ok) memcpy( qAbk, qA, sizeof(qA) ); qN = qNN;
  }; // for iter
  if (!ixAscArcFrmPfDbg) return 2;
  int qV; // debug_output when AscendFrmPfDbg=true 
  qV = (int)(cA1mn*(1<<cAscendGainLg)); ADat[2] = mcRng( qV, -0x7FFF, 0x7FFF );
  qV = (int)(cA1mx*(1<<cAscendGainLg)); ADat[3] = mcRng( qV, -0x7FFF, 0x7FFF );
  BS_ZFil( AL-4, ADat+4 ); for (int qi=0;qi<qN;qi++) { int qx = qBufX[qi]+(AL>>1); if (qx>1) ADat[qx] = qBufY[qi]; };
  return AL; }
static float iGatherPfAve( int AW, int AH, const tBS_Val *ACol ) {
  int qZ = 0, qN = 0; for (int qi=0;qi<AH;qi++) { int qV = ACol[qi*AW]; if (qV>0) { qZ += qV; qN++; }; };
  return qN<(AH>>3) ? 0 : (float)qZ/qN; }
void BS7c_AscArcFrmPf( tBS_Val *ADat, tBS7_AscArcFrmPf *AAFP,
  const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- ascending polyfit of frame rows
  if (0 && ADbgPfx) { // inp limbo_metrics's frame "as-is" 
    strcpy( lDbgSfx, "A.raw" ); BS_Write( ADbgPfx, cN, ADat ); };
  //BS_VScat( AH, ADat+0,AW, -2 ); BS_VScat( AH, ADat+1,AW, -2 );
  const int cHp0 = (int)(AH*0.2f), cHp1 = AH-(int)(AH*0.15f), cHp = cHp1-cHp0;
  BS_ZFil( cHp0*AW, ADat ); BS_ZFil( (AH-cHp1)*AW, ADat+cHp1 );
  ixAscArcFrmPfDbg = ADbgPfx!=NULL;
  BS_CPLCRows( icpl_AscArcFrmPf, AW, cHp, ADat+cHp0*AW );
  if (1 && ADbgPfx) { // rows' polyfit in columns ~0..1
    if (!ixAscArcFrmPfDbg) { // if cpl_AscendFrmPf didn't write already it's debug-info into row's pixels
      for (int qr=cHp0;qr<cHp1;qr++) { tBS_Val *qRow = ADat+qr*AW;
        BS7c_AscArcRowLimbo( AW-2, AW>>1, qRow+2, qRow[0], (float)qRow[1]/(1<<cAscendGainLg) ); };
    }; strcpy( lDbgSfx, "B.raw" ); BS_Write( ADbgPfx, cN, ADat ); };
  // return frame polyfit's stats
  tBS_Val    *qDat = ADat+cHp0*AW;
  AAFP->rA[0] = iGatherPfAve( AW, cHp, qDat+0 ),
  AAFP->rA[1] = iGatherPfAve( AW, cHp, qDat+1 )/(1<<cAscendGainLg); }

// polyfit_skew of a frame by A0,A1
static void ppr_PfSkewFrm( tBS_SpawnThreadRec *AP ) { // skew section of rows
  int cW = mcVi(AP->rP[0]), cH = cW>>16; cW &= 0x7FFF; int cN = cW*cH;
  tBS_Val *cSrc = (tBS_Val*)AP->rP[1], *qDst = (tBS_Val*)AP->rP[2]; int *cSkewOfs = (int*)AP->rP[3];
  int cLst = AP->rLst; for (int qr=AP->rFst;qr<cLst;qr++) { const int cO = qr*cW;
    for (int qc=0;qc<cW;qc++) { int qOdst = cO+qc, qOsrc = qOdst+cSkewOfs[qc]; if (qOsrc>=0 && qOsrc<cN) qDst[qOdst] = cSrc[qOsrc]; };
  }; }
static int iSkewFillMgnVal( const tBS_Val *ADat, int AW, int AH ) { const int cN = AW*AH; int qFill = ADat[0];
  const int cFillCandidates[] = {AW-1,cN-AW,cN-1}; for (int qi=0;qi<mcASz(cFillCandidates);qi++) {
  int qV = ADat[cFillCandidates[qi]]; if (qFill<qV) qFill = qV; }; return qFill; }
void BS7c_PfSkewFrm( bool AHrzToArc, tBS_Val *ADat, float A0, float A1, int AFill,
  const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- polyfit_skew of a frame by A0,A1
  if (0 && ADbgPfx) { strcpy( lDbgSfx, "A.raw" ); BS_Write( ADbgPfx, cN, ADat ); }; // inp frame
  tBS_Val *qDat = new tBS_Val[cN]; BS_VFil( cN, qDat,
    (AFill<-32768)||(AFill>32767)?iSkewFillMgnVal(ADat,AW,AH):AFill ); // fill with ~ribcage_frame pixel
  tBS_Val *qRowLimb = new tBS_Val[AW]; BS7c_AscArcRowLimbo( AW, AW>>1, qRowLimb, A0, A1 ); // fill row by limbo_metric
  tBS_Val *qRowSkew = new tBS_Val[AW]; BS_Cpy( AW, qRowLimb, qRowSkew ); BS7c_AscArcRowSkew( AW, AW>>1, qRowSkew ); // row's skew by Y
  const int cWmul = AHrzToArc?-AW:AW; // skew sections of rows
  int     *qSkewOfs = new int[AW]; for (int qi=0;qi<AW;qi++) qSkewOfs[qi] = qRowSkew[qi]*cWmul;
  BS_SpawnThreads( 0, 0, AH, ppr_PfSkewFrm, mcIv((AH<<16)|AW), ADat, qDat, qSkewOfs );
  BS_Cpy( cN, qDat, ADat );
  if (1 && ADbgPfx) { BS_Cpy( AW, qRowLimb, qDat+0 ); // out frame
    BS_Cpy( AW, qRowSkew, qDat+AW ); strcpy( lDbgSfx, "B.raw" ); BS_Write( ADbgPfx, cN, qDat ); };
  delete qDat; delete qRowLimb; delete qRowSkew; delete qSkewOfs; }

//=== Apps of Descend periphery to center ====b
//============================================
// plain fast_descend
static int icpl_FastDescend( int AL, tBS_Val *ADat ) { mcWSMAd16( AL, ADat ); return AL; }
void BS7c_FastDescend( const tBS_Val *AOrgConvex, int ARadBlr, tBS_Val **ARes, 
  int AW, int AH ) { // -- plain fast_descend
  const int cN = AW*AH; BS_DataDup( 8, cN, AOrgConvex, ARes );
  BS_CPLA( icpl_FastDescend, 8, cBS_AllDirsT, AW, AH, ARes );
  BS_CPB( mccpb_Mn<8>, cCPBNDD|8, cN, ARes ); // #0:fast_descend
  BS_MovAveX( ARadBlr, AW, AH, ARes[0] ); } // #0:fast_descend_blr

// excess of slow descend
static int icpl_DescendXcs( int AL, tBS_Val *ADat ) { tBS_Val *qDat = ADat+BS_Align( AL );
  BS_Cpy( AL, ADat, qDat ); mcWSMAd16( AL, qDat ); BS_SubA( AL, qDat, ADat ); return AL; }
void BS7c_DescendXcs( int ANDD, const int *ADirs, const tBS_Val *AOrgConvex, tBS_Val **ARes, tBS_Val **ATs,
  int AW, int AH ) { // -- excess of slow descend 
  const int cN = AW*AH; BS_DataDup( ANDD, cN, AOrgConvex, ATs );
  BS_CPLA( icpl_DescendXcs, ANDD, ADirs, AW, AH, ATs ); // Ts:plain descends
  tBS_CPBProc qProc; switch (ANDD) {
    case 3: qProc = mccpb_Mn<3>; break; case 4: qProc = mccpb_Mn<4>; break; case 5: qProc = mccpb_Mn<5>; break; 
    case 6: qProc = mccpb_Mn<6>; break; case 7: qProc = mccpb_Mn<7>; break; case 8: qProc = mccpb_Mn<8>; break;
    default: ErrMsg( "BS7_DescendXcs NDD" ); return; };
  BS_CPB( qProc, ANDD|cCPBNDD, cN, ATs, ARes[0] ); }

// excess of skewed slow_descend
void BS7c_PfSkewDescendXcs( const tBS_Val *AOrgConvex, float A0, float A1, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- excess of skewed slow_descend
  const int cFill = 99999; // ~autodetect (here we deal with ~initial_convex)
  BS_Cpy( cN, AOrgConvex, ARes[0] );
  if (ADbgPfx) strcpy( lDbgSfx, "Askw" ); BS7c_PfSkewFrm( false, ARes[0], A0, A1, cFill, ADbgPfx, AW, AH );
  //  ADbgPfx = NULL; // !!!!
  #ifdef BS7_ASCENDDIRLESS
    const int cDirs[] = {8,1,2};
  #else
    const int cDirs[] = {7,8,1,2,3};
  #endif
  BS7c_DescendXcs( mcASz(cDirs), cDirs, ARes[0], ARes, ATs, AW, AH ); // out#0:by 3..5 dirs, #1:by 1..3 dirs
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "Bxcs" ); BS_DatasWrt( ADbgPfx, cN, ARes, 0, 2 ); };
  if (1 && ADbgPfx) strcpy( lDbgSfx, "Cskw0" ); BS7c_PfSkewFrm( true, ARes[0], A0, A1, cFill, ADbgPfx, AW, AH ); } //


//=== Treatment of maxmag & associates ====
//=========================================
// estimate horz_blur(main_ribcage_lap)
static int icpl_PfSkewMaxMag( int AL, tBS_Val *ADat ) { tBS_Val *qDat = ADat+BS_Align( AL );
  BS_MovAve1D( AL, 20, ADat, qDat ); BS_Cpy( AL, qDat, ADat ); return AL; }; 
void BS7c_PfSkewMaxMag( const tBS_Val *ALap, float A0, float A1, bool ATheFinal, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- estimate horz_blur(main_ribcage_lap)
  const int cFill = 0; // as here we deal with a laplacian 
  tBS_Val *qRes = ARes[0]; BS_Cpy( cN, ALap, qRes );
  if (ADbgPfx) strcpy( lDbgSfx, "Askw" ); BS7c_PfSkewFrm( false, qRes, A0, A1, cFill, ADbgPfx, AW, AH );
  tBS_Val *qBu0 = ATs[0]; if (ATheFinal) BS_Cpy( cN, qRes, qBu0 );
  BS_CPLCRows( icpl_PfSkewMaxMag, AW, AH, qRes );
  ADbgPfx = NULL; // !!!!
  if (1 && ATheFinal) {
    tBS_Val *qBuf = ATs[1], *cDs[2] = {qRes,qBuf};
    const int cShearY[] = {-7,7}, cShearN = mcASz(cShearY);
    for (int qskw=0;qskw<cShearN;qskw++) { BS_Cpy( cN, qBu0, qBuf ); //
      BS7_SkewV( cShearY[qskw], true, qBuf, AW, AH ); 
      if (0 && ADbgPfx) { sprintf( lDbgSfx, "Bskw%da.raw", qskw ); BS_Write( ADbgPfx, cN, qBuf ); };
      BS_CPLCRows( icpl_PfSkewMaxMag, AW, AH, qBuf );
      BS7_SkewV( -cShearY[qskw], true, qBuf, AW, AH );
      if (0 && ADbgPfx) { sprintf( lDbgSfx, "Bskw%db.raw", qskw ); BS_Write( ADbgPfx, cN, qBuf ); };
      BS_CPB( mccpb_MxA<2>, 2|cCPBNDD, cN, cDs );
    }; // for skw
  }; // if TheFinal
  if (ADbgPfx) strcpy( lDbgSfx, "Cskw" ); BS7c_PfSkewFrm( true, qRes, A0, A1, cFill, ADbgPfx, AW, AH ); }

struct tMMRWr { // - MaxMag's refinement wrapper
  // -- intrinsics of MMRWr 
  int  rW, rH, rCYA[2]; // average Ycounts (|negative|,|positive|)
  char rDbgPfx[80], *rDbgSfx;
  void iDbgW( char APfx, int AIdx, const char *ANam, const tBS_Val *ADat ) { if (!rDbgSfx) return; // debug file-write
    if (AIdx<0) sprintf( rDbgSfx, "%c%s.raw", APfx, ANam ); else sprintf( rDbgSfx, "%c%d%s.raw", APfx, AIdx, ANam );
    BS_Write( rDbgPfx, rW*rH, ADat ); }
  void iMxA2( tBS_Val **ADs2 ) { BS_CPB( mccpb_MxA<2>, 2|cCPBNDD, rW*rH, ADs2 ); } // max of 2 arrays

  // Ycounts
  static int icpl_CntY( int AL, tBS_Val *ADat ) { mcIvlC( AL, ADat ); return AL; }
  void iCntY( tBS_Val *ADat ) { BS_CPLCCols( icpl_CntY, rW, rH, ADat ); }
  void iCntYA( const tBS_Val *ACntY ) { const int cN = rW*rH; // calculate "CYA"
    int qA0 = 0, qA1 = 0, qC0 = 0, qC1 = 0; for (int qi=0;qi<cN;qi++) { int qV = ACntY[qi]; if (!qV) continue;
    if (qV<-11 && qV>-66) { qA0 -= qV; qC0++; } else if (qV>9 && qV<55) { qA1 += qV; qC1++; }; }; // for i
    rCYA[0] = qC0?(qA0+(qC0>>1))/qC0:0; rCYA[1] = qC1?(qA1+(qC1>>1))/qC1:0; } // of iCntYA

  // Ymetric (calc from Ycounts using "CYA"; a metric of fixed narrow range)
  static void itppr_MetY( tBS_SpawnThreadRec *AP ) { tMMRWr *cSelf = (tMMRWr*)AP->rP[0];
    tBS_Val *qDat = (tBS_Val*)AP->rP[1]; const int cLg = 16, cLst = AP->rLst,
      cA0 = cSelf->rCYA[0], cSc0 = cA0?(1<<cLg)/cA0:0, cA1 = cSelf->rCYA[1], cSc1 = cA1?(1<<cLg)/cA1:0;
    for (int qi=AP->rFst;qi<cLst;qi++) { int qV = qDat[qi]; if (!qV) continue;
      if (qV<0) { qV = cA0-abs(cA0+qV); qDat[qi] = qV<=0?0:-((qV*cSc0)>>(cLg-4)); continue; };
      qV = cA1-abs(cA1-qV); qDat[qi] = qV<=0?0:(qV*cSc1)>>(cLg-4); }; } // for i
  void iMetY( tBS_Val *ADat ) { // inp:Ycounts, out:Ymetric (0..15)
    BS_SpawnThreads( 0, 0, rW*rH, itppr_MetY, this, ADat ); }

  // Ymoment (calc from Ymetric)
  static int icpl_MomY( int AL, tBS_Val *ADat ) {
    tBS_Val *qDat = ADat+BS_Align( AL+1 ); BS_Cpy( AL, ADat, qDat+1 ); BS_Abs( AL, qDat );
    qDat[0] = qDat[1]; qDat[AL+1] = qDat[AL]; // simplify loop's condition
    int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; 
      if (qV>0) { if (qL>=0) { qL++; continue; }; goto lab1; };
      if (qV<0) { if (qL<=0) { qL--; continue; }; goto lab1; };
      if (!qL) continue; lab1:
      int qjL = qi-abs(qL), qML = (qDat[qjL]+qDat[qjL+1]),
          qjR = qi-1, qMR = (qDat[qjR+1]+qDat[qjR+2]);
      for (;;) {
        ADat[qjL++] = qML*qML; if (qjL>qjR) break; if (qML>0) qML--;
        ADat[qjR--] = qMR*qMR; if (qjL>qjR) break; if (qMR>0) qMR--; };
      qL = mcRng(qV,-1,1); }; return AL; } // for i
    void iMomY( tBS_Val *ADat ) { BS_CPLCCols( icpl_MomY, rW, rH, ADat ); }

  // Xmetric (calc from Ymetric)
  static int icpl_MetX( int AL, tBS_Val *ADat ) { // called along ~horz direction
    int qJ0 = 0;    while ((!ADat[qJ0]) && (qJ0<AL)) qJ0++;
    int qJ1 = AL-1; while ((!ADat[qJ1]) && (qJ1>qJ0)) qJ1--; ADat[++qJ1] = 0; // simplify loop's condition
    const int cJc = (qJ0+qJ1)>>1, cJm = (qJ1-qJ0)/5;
    int qL = 0, qZ = 0; for (int qi=qJ0;qi<=qJ1;qi++) { int qV = ADat[qi]; 
      if (cJm-abs(cJc-qi)<0) qV = mcRng( qV, -1, 1 );
      if (qV>0) { if (qZ<0) goto lab1; qL++; qZ += qV; continue; };
      if (qV<0) { if (qZ>0) goto lab1; qL++; qZ += qV; continue; };
      if (!qL) continue; lab1: BS_VFil( qL, ADat+qi-qL, qZ ); qL = qV?1:0; qZ = qV;
    }; return AL; } // for i, of icpl_IvlZ
  void iMetX( tBS_Val *ADat, tBS_Val **ATs ) { const int cN = rW*rH; // inp:Ymetric, out:Xmetric ~"-16000..16000"
    BS_Cpy( cN, ADat, ATs[0] ); BS_CPLCRows( icpl_MetX, rW, rH, ATs[0] );
    const int cShearY[] = {-7,7}; for (int qskw=0;qskw<mcASz(cShearY);qskw++) {
      BS_Cpy( cN, ADat, ATs[1] );
      BS7_SkewV( cShearY[qskw], true, ATs[1], rW, rH ); 
      BS_CPLCRows( icpl_MetX, rW, rH, ATs[1] );
      BS7_SkewV(-cShearY[qskw], true, ATs[1], rW, rH ); 
      iMxA2( ATs ); }; // for skw
    BS_Cpy( cN, ATs[0], ADat ); } // of iMetX

  // a version of laplacian (split by ~"signed_count")
  static int icpl_IvlCnZ( int AL, tBS_Val *ADat ) { mcIvlCnZ( AL, ADat ); return AL; }
  void iLapVersionByCnt( const tBS_Val *ALap, int ACntMn, int AReplace, int ANDD, const int *ADirs, tBS_Val **ADs ) { //
    const int cN = rW*rH; tBS_Val *qT = ADs[0]; switch (AReplace) {
      case -1: for (int qi=0;qi<cN;qi++) qT[qi] = ALap[qi]<0?0:1; break; // count pixels >=0
      case  1: for (int qi=0;qi<cN;qi++) qT[qi] = ALap[qi]>0?0:1; break; // count pixels <=0
    };
    BS_DataDup( ANDD-1, cN, qT, ADs+1 );
    BS_CPLA( icpl_IvlCnZ, ANDD, ADirs, rW, rH, ADs );
    tBS_CPBProc cProcs[] = {mccpb_Mx<2>,mccpb_Mx<3>}; if (ANDD>mcASz(cProcs)) {
      ErrMsg( "SplitLapVersion(NDD>3)" ); return; }; // no much dirs in transient optimisation
    if (ANDD>0) BS_CPB( cProcs[ANDD-1], cCPBNDD|ANDD, cN, ADs );
    switch (AReplace) {
      case -1: for (int qi=0;qi<cN;qi++) qT[qi] = qT[qi]>ACntMn ? ALap[qi]:__min(ALap[qi],-1); break;
      case  1: for (int qi=0;qi<cN;qi++) qT[qi] = qT[qi]>ACntMn ? ALap[qi]:__max(ALap[qi], 1); break;
    }; } // of iLapVersionBySign

  // a version of laplacian (split by ~"signed_value_threshold")
  void iLapVersionByThr( const tBS_Val *ALap, int AThrMn, int AReplace, tBS_Val **ADs ) { // 
    const int cN = rW*rH; tBS_Val *qT = ADs[0]; BS_Cpy( cN, ALap, qT ); switch (AReplace) {
      case -1: for (int qi=0;qi<cN;qi++) { int qV = qT[qi]; if (qV>0 && qV<AThrMn)  qT[qi] =-1; }; break;
      case  1: for (int qi=0;qi<cN;qi++) { int qV = qT[qi]; if (qV>-AThrMn && qV<0) qT[qi] = 1; }; break;
      default: for (int qi=0;qi<cN;qi++) { int qV = qT[qi]; if (abs(qV)<AThrMn) qT[qi] = -qV; };
    }; } // of iLapVersionByThr

  // max anisotropic laplacian
  void iMxLap( tBS_Val *ALap, int ARad, tBS_Val **ATs ) { const int cN = rW*rH; // 
    const int cShearY[] = {-5,-8,-11,0,11,8,5}; for (int qskw=0;qskw<mcASz(cShearY);qskw++) {
      BS_Cpy( cN, ALap, ATs[1] ); BS7_SkewV( cShearY[qskw], true, ATs[1], rW, rH ); 
      BS_MovAveX( ARad, rW, rH, ATs[1] ); BS7_SkewV( -cShearY[qskw], true, ATs[1], rW, rH ); 
      if (qskw) iMxA2( ATs ); else BS_Cpy( cN, ATs[1], ATs[0] ); }; // for skw
    BS_Cpy( cN, ATs[0], ALap ); }

  void iConstraintLap( const tBS_Val *AOrgLap, tBS_Val *ALap ) { const int cN = rW*rH; // 
    for (int qi=0;qi<cN;qi++) { int qV1 = ALap[qi]; if (!qV1) continue;
      int qV0 = AOrgLap[qi]; if ((qV0^qV1)&(1<<31)) { ALap[qi] = 0; continue; };
      if (abs(qV0)<abs(qV1)) ALap[qi] = qV0; }; }

  public: // -- publics of MMRWr
  void Init( int AW, int AH, const char *ADbgPfx ) {
    rW = AW; rH = AH; if (ADbgPfx) { strcpy( rDbgPfx, ADbgPfx ); rDbgSfx = strend(rDbgPfx); }; }

  void RibsRefine( const tBS_Val *ARibs, const tBS_Val *AMgnCon, const tBS_Val *ACrd, tBS_Val **ADs ) { // "HB*.raw"
    //  inp:  |\/|  _       inp:        _       out:        _     out:  |\/|
    // Ribs: _|  |_| |_  MgnCon: _|\/|_| |_  MgnCon: ______| |_  Ribs: _|  |_____
    // Assumes that inp_MgnCon contains mostly nodules and occacsional ribs' residuals.
    // Remove ribs' residuals in prcs_MgnCon by anisotropic spatial filters (while keep
    // nodules). Then "out_Ribs:=inp_Ribs-prcs_MgnCon". Out_Ribs should contain
    // ~untouched ribs but remove nodules.
    const int cN = rW*rH;
    tBS_Val **pB = ADs+2;
    BS_Cpy( cN, ACrd, pB[0] ); BS_Neg( cN, pB[0] ); BS_VAdd( cN, pB[0], 1<<12 );
    BS_Rng( cN, pB[0], -99999, 3333 ); BS_VMul( cN, pB[0], (float)(1<<12)/3333 );
    if (1) iDbgW( 'B', -1, "acrd", pB[0] ); // #0:crd
    // anisotropic spatial filter
    tBS_Val **pD = pB+2;
    const int cShearY[] = {0,-5,5}; for (int qskw=0;qskw<mcASz(cShearY);qskw++) { // 
      const int t = __min(qskw,1); BS_Cpy( cN, AMgnCon, pD[t] );
      BS7_SkewV( cShearY[qskw], false, pD[t], rW, rH ); 
      BS_MovAveX( 17, rW, rH, pD[t] ); // blurX(MgnCon)
      BS7_SkewV( -cShearY[qskw], false, pD[t], rW, rH ); 
      if (t) BS_CPB( mccpb_Mx<2>, 2|cCPBNDD, cN, pD ); }; // for skw; #0:maxlap
    BS_Mul( cN, pB[0], pD[0], (int)((1<<12)*0.3f) );
                                       if (1) iDbgW( 'D', -1, "a", pD[0] ); // maxlap *= 1.5..4.0
    // keep only lapacianY(mediumR,maxlap ) 
    BS_Cpy( cN, pD[0], pD[1] ); BS_MovAveY( 7, rW, rH, pD[1] );
    BS_SubA( cN, pD[1], pD[0] );       if (0) iDbgW( 'D', -1, "b", pD[0] );
    iConstraintLap( AMgnCon, pD[0] );  if (0) iDbgW( 'D', -1, "c", pD[0] ); // restrict maxlap against MgnCon
    BS_SubS( cN, AMgnCon, pD[0] );     if (1) iDbgW( 'D', -1, "d", pD[0] ); // D#0:out_MgnCon
    mcShiftLft( cN, pD[0], 1 ); BS2_MovAve2D( 1, pD[0], ADs[1], rW, rH );  // #1:blur(out_MgnCon*2)
    iConstraintLap( AMgnCon, ADs[1] ); if (0) iDbgW( 'D', -1, "e", ADs[1] ); // #1:restrict(out_MgnCon)
    // out_Ribs:=inp_Ribs-out_MgnCon
    BS_Cpy( cN, ARibs, ADs[0] ); BS_SubA( cN, ADs[1], ADs[0] ); // #0:out_Ribs
                                       if (0) iDbgW( 'F', -1, "rib", ADs[0] ); 
    // empirically decrease negatives in out_Ribs
    tBS_Val **pF = pD+1;
    BS_Cpy( cN, ADs[0], pF[0] ); BS_Rng( cN, ADs[0], 0 ); BS_Rng( cN, pF[0], -999999, 0 );
    BS_Cpy( cN, pB[0], pF[1] ); BS_VAdd( cN, pF[1], 1<<10 );
    BS_Mul( cN, pF[1], pF[0], (1<<12)+(1<<10) );
    BS_Add( cN, pF[0], ADs[0] ); iConstraintLap( ARibs, ADs[0] );
    if (1) { iDbgW( 'V', -1, "aribs", ADs[0] ); iDbgW( 'V', -1, "bmgncon", ADs[1] ); };
  } // of RibsRefine

  void MaxMagRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag, tBS_Val **ADs ) { // "HC*.raw"
    // inp:  |\/|    |\/|   out:  |\/|    |\/|     out: |\  _    _  /| 
    //  MM: _|  |_/\_|  |_   MM: _|  |____|  |_  Split: | \/ \__/ \/ |
    // 1) MaxMag has 'suggested' Ycounts for hulls,convexes (i.e. for pos,neg laps). 
    // 2) YMetric: "|Ycount-Ycount_suggested|=0":"metric=255", "...=big":"metric=0".
    // 3) XMetric: adaptive sum of YMetric along horz dir
    char qS[256]; const int cN = rW*rH;
    // Ymetric
    tBS_Val **pYMet = ADs; BS_Cpy( cN, AMaxMag, pYMet[0] );
    iCntY( pYMet[0] ); if (1) iDbgW( 'B', 0, "aycnt", pYMet[0] ); // Ycounts
    iCntYA( pYMet[0] ); // calculate "CYA"
    rCYA[0] = (int)(rCYA[0]*1.0f); // empirical attenuation coefficients
    rCYA[1] = (int)(rCYA[1]*1.0f); // for "suggested Y_count_averages"
    if (rDbgSfx) { sprintf( qS, cTPS" inf RibCage MMRef suggests CYA[neg=%d,pos=%d]",
      rCYA[0], rCYA[1] ); WriteLn( qS ); };
    iMetY( pYMet[0] ); if (0) iDbgW( 'B', 0, "bymet", pYMet[0] ); // Ymetric
    for (int qvbc1=0;qvbc1<2;qvbc1++) { const int cDirsY[]= {5};
      iLapVersionByCnt( AMaxMag, 5, qvbc1?1:-1, mcASz(cDirsY), cDirsY, pYMet+1 ); // a version of laplacian
      iCntY( pYMet[1] ); if (0) iDbgW( 'B', qvbc1+1, "aycnt", pYMet[1] ); // Ycounts
      iMetY( pYMet[1] ); if (0) iDbgW( 'B', qvbc1+1, "bymet", pYMet[1] ); // Ymetric
      iMxA2( pYMet ); }; // for vbc1
    if (1) iDbgW( 'B', -1, "Uymet", pYMet[0] );
    // Xmetric
    tBS_Val **pXMet = pYMet+1; BS_Cpy( cN, pYMet[0], pXMet[0] );
    iMetX( pXMet[0], pXMet+2 ); if (0) iDbgW( 'C', 0, "xmet", pXMet[0] ); // Xmetric
    for (int qvbc2=0;qvbc2<2;qvbc2++) { const int cDirsX[]= {1}; 
      iLapVersionByCnt( pYMet[0], 9, qvbc2?1:-1, mcASz(cDirsX), cDirsX, pXMet+1 ); // a version of laplacian
      iMetX( pXMet[1], pXMet+2 ); if (0) iDbgW( 'C', qvbc2+1, "xmet", pXMet[0] ); // Xmetric
      iMxA2( pXMet ); }; // for vbc2
    if (1) iDbgW( 'C', -1, "Uxmet", pXMet[0] );
    // re-check Xmetric for variations of "VersionByThr"
    tBS_Val **pUMet = pXMet+1; BS_Cpy( cN, pXMet[0], pUMet[0] );
    iCntY( pUMet[0] ); if (0) iDbgW( 'D', 0, "aycnt", pUMet[0] ); // Ycounts
    iMetY( pUMet[0] ); if (0) iDbgW( 'D', 0, "bymet", pUMet[0] ); // Ymetric
    iMetX( pUMet[0], pUMet+2 ); if (0) iDbgW( 'D', 0, "cxmet", pUMet[0] ); // Xmetric
    for (int qvbt=0;qvbt<3;qvbt++) {
      iLapVersionByThr( pXMet[0], 567, qvbt-1, pUMet+1 ); // a version of laplacian 888 !!!!
      iCntY( pUMet[1] ); if (0) iDbgW( 'D', qvbt+1, "aycnt", pUMet[1] ); // Ycounts
      iMetY( pUMet[1] ); if (0) iDbgW( 'D', qvbt+1, "bymet", pUMet[1] ); // Ymetric
      iMetX( pUMet[1], pUMet+2 ); if (0) iDbgW( 'D', qvbt+1, "cxmet", pUMet[1] );  // Xmetric
      iMxA2( pUMet ); }; // for vbt
    if (1) iDbgW( 'D', -1, "Uxmet", pUMet[0] );
    //
    #if 0
      iDbgW( 'E', -1, "Qa", ARibs );
      tBS_Val **pT = pUMet+1; BS_Cpy( cN, pUMet[0], pT[0] );
      tBS_Val *qM = pT[0]; for (int qi=0;qi<cN;qi++) {
        int qV0 = ARibs[qi]; if (!qV0) { qM[qi] = 0; continue; };
        int qV1 = qM[qi]; if (!qV1) continue;
        bool c0 = qV0&(1<<31), c1 = qV1&(1<<31);
        qM[qi] = c0==c1 ? 0 : (c1?1:-1);
      };
      iDbgW( 'E', -1, "Qb", pT[0] );
    #endif
    // separate and recombine masks for MaxMag_neg,_pos
    tBS_Val **pMMrf = pUMet+1; BS_DataDup( 2, cN, pUMet[0], pMMrf+0 );
    BS_Rng( cN, pMMrf[0], -99999, 0 ); BS_Neg( cN, pMMrf[0] ); BS_Rng( cN, pMMrf[1], 0 ); 
    BS_VAdd( rW*(rH/5), pMMrf[0], 999 ); BS_VAdd( rW*(rH/5), pMMrf[1], 999 );
    const int cRadBlrY = 9;
    BS_MovAveY( cRadBlrY, rW, rH, pMMrf[0] ); BS_MovAveY( cRadBlrY, rW, rH, pMMrf[1] );
    if (0) for (int qi=0;qi<2;qi++) iDbgW( 'H', qi, "msk", pMMrf[qi] );
    const int cBlackA = 0, cWhiteA = 678;
    tBS2_Cut cCutA = { {cBlackA,cWhiteA,512}, {0,0} };
    BS2_Cut( &cCutA, pMMrf[0], cN ); BS2_Cut( &cCutA, pMMrf[1], cN );
    BS_DataDup( 2, cN, AMaxMag, pMMrf+2 );
    BS_Rng( cN, pMMrf[2], -99999, 0 ); BS_Neg( cN, pMMrf[2] ); BS_Rng( cN, pMMrf[3], 0 ); 
    BS_Mul( cN, pMMrf[0], pMMrf[2], -9 ); BS_Neg( cN, pMMrf[2] ); BS_Mul( cN, pMMrf[1], pMMrf[3], -9 );
    if (1) for (int qi=0;qi<4;qi++) iDbgW( 'I', qi, "msk", pMMrf[qi] );
    BS_Add( cN, pMMrf[3], pMMrf[2] ); // #2:maxmag_rf
    BS_Rng( cN, pMMrf[2], -333, 99999 ); if (1) iDbgW( 'I', -1, "mm", pMMrf[2] ); // limit deep concaves
    // Ysplitter (masks to "determine split_by_filter", "mul split_by_filter")
    tBS_Val **pSplt = pMMrf+3; BS_Cpy( cN, pUMet[0], pSplt[0] );
    iCntY( pSplt[0] ); if (0) iDbgW( 'O', -1, "aspt", pSplt[0] );
    iMetY( pSplt[0] ); if (0) iDbgW( 'O', -1, "bspt", pSplt[0] );
    iMomY( pSplt[0] ); if (0) iDbgW( 'O', -1, "cspt", pSplt[0] );
    // const int cBlackB = 267, cWhiteB = 543; // #0:mm_sbf_det (rough but generally relevant shape)
    const int cBlackB = 111, cWhiteB = 444; // #0:mm_sbf_det (rough but generally relevant shape)
    tBS2_Cut cCutB = { {cBlackB,cWhiteB,512}, {0,0} }; BS2_Cut( &cCutB, pSplt[0], cN ); 
    // leave
    BS_Cpy( cN, pMMrf[2], ADs[0] ); if (1) iDbgW( 'V', -1, "amm", ADs[0] ); // #0:maxmag_rfn
    BS_Cpy( cN, pSplt[0], ADs[1] ); if (1) iDbgW( 'V', -1, "bspt", ADs[1] ); // #1:maxmag_sms
  } // of MaxMagRefine

  void CombineRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag, tBS_Val **ADs ) { // "HD*.raw"
    //  inp:  |\/|  _        inp:  /^^\         out:  |\/|
    // Ribs: _|  |_| |_   MaxMag: _|  |____  MaxMag: _|  |_____
    char qS[256]; const int cN = rW*rH;
    // remove overmagnitude from maxmag;normally expect "|Shrp|<|Ambi|"
    tBS_Val **pOvr = ADs; // -- to remove over_magnitude; 
    BS_Cpy( cN, AMaxMag, pOvr[0] ); iConstraintLap( ARibs, pOvr[0] ); if (0) iDbgW( 'B', -1, "amm", pOvr[0] ); // #0:constrained MaxMag
    BS_SubS( cN, AMaxMag, pOvr[0] ); if (1) iDbgW( 'B', -1, "bom", pOvr[0] ); // #0:overmagnitude
    iMxLap( pOvr[0], 17, pOvr+1 );   if (1) iDbgW( 'B', -1, "cmx", pOvr[0] ); // #0:blur_overmagnitude
    BS_SubS( cN, AMaxMag, pOvr[0] ); if (1) iDbgW( 'B', -1, "dmm", pOvr[0] ); // #0:maxmag_om 'maxmag-blur_overmagnitude'
                                  // if (1) iDbgW( 'B', -1, "drib", ARibs ); // 'HBVrib.raw' above
    BS_Cpy( cN, ARibs, pOvr[1] ); //#1:difference 'Ribs-maxmag_om'
    BS_SubA( cN, pOvr[0], pOvr[1] ); if (1) iDbgW( 'B', -1, "edif", pOvr[1] ); 
    // determine top thin band (matching ~|vertical_lap|)
    tBS_Val **pBnd = pOvr+2;
    const int c1span = 2, c2span = c1span<<1, cO2 = rW*c2span, cN2 = cN-cO2;
    BS_Cpy( cN, pOvr[0], pBnd[0] ); BS_SubA( cN2, pOvr[0]+cO2, pBnd[0] );
    BS_Abs( cN2, pBnd[0] ); BS_ZFil( cO2, pBnd[0]+cN2 ); 
    tBS2_Cut cCut = { {111,444,1<<12}, {0,0} }; BS2_Cut( &cCut, pBnd[0], cN ); // #0:maxmag_bndY0
    // loop by thin bands
    tBS_Val **pLap = pBnd+2;
    for (int qb=0;qb<3;qb++) { // Bnd#1: actual band
      const int cO = c1span*rW*qb; BS_Cpy( cN-cO, pBnd[0], pBnd[1]+cO );
      BS_ZFil( cO, pBnd[1] );               if (0) { sprintf( qS, "%da", qb ); iDbgW( 'C', -1, qS, pBnd[1] ); }; 
      const int ct = __min(1,qb); BS_Cpy( cN, pOvr[1], pLap[ct] ); // Lap#ct
      BS_Mul( cN, pBnd[1], pLap[ct], -12 ); if (0) { sprintf( qS, "%db", qb ); iDbgW( 'C', -1, qS, pLap[ct] ); };
      iMxLap( pLap[ct], 12, pLap+2 );       if (1) { sprintf( qS, "%dc", qb ); iDbgW( 'C', -1, qS, pLap[ct] ); };
      if (ct) iMxA2( pLap+0 ); }; // #0:maxlap; for b
                                     if (0) iDbgW( 'D', -1, "alap", pLap[0] ); // for b
    BS_Add( cN, pOvr[0], pLap[0] );  if (1) iDbgW( 'D', -1, "bmm", pLap[0] ); // maxmag_rfn
    iMxLap( pLap[0], 3, pLap+1 ); iConstraintLap( ARibs, pLap[0] ); 
    BS_Cpy( cN, pLap[0], ADs[0] );   if (1) iDbgW( 'V', -1, "mm", ADs[0] ); // #0:maxmag_rfn
  } // of CombineRefine
}; // end of MMRWr 

void BS7c_RibsRefine( const tBS_Val *ARibs, const tBS_Val *AMgnCon, const tBS_Val *ACrd, 
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ) { // -- Ribs -= MgnCon
  tMMRWr qMMR; memset( &qMMR, 0, sizeof(qMMR) ); qMMR.Init( AW, AH, ADbgPfx );
  qMMR.RibsRefine( ARibs, AMgnCon, ACrd, ARes ); } // invoke MMRWr.RibsRefine(...)

void BS7c_MaxMagRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag,
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ) { // -- empirically refine MaxMag
  tMMRWr qMMR; memset( &qMMR, 0, sizeof(qMMR) ); qMMR.Init( AW, AH, ADbgPfx );
  qMMR.MaxMagRefine( ARibs, AMaxMag, ARes ); } // invoke MMRWr.MaxMagRefine(...) 

void BS7c_CombineRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag,
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ) { // -- combine refines
  tMMRWr qMMR; memset( &qMMR, 0, sizeof(qMMR) ); qMMR.Init( AW, AH, ADbgPfx );
  qMMR.CombineRefine( ARibs, AMaxMag, ARes ); } // invoke MMRWr.CombineRefine(...) 

//=== Clavicle's estimation ====
//==============================
// plain directional descends
static int icpl_DirDescend( int AL, tBS_Val *ADat ) { mcWSMAd16( AL, ADat ); return AL; }
void BS7c_DirDescends( const tBS_Val *AConvex, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- plain directional descends
  BS_Cpy( cN, AConvex, ARes[0] ); 
  // non-skewed "plain directional descends"
  BS_DataDup( 8, cN, AConvex, ARes );
  BS_CPLA( icpl_DirDescend, 8, cBS_AllDirsT, AW, AH, ARes );
  // optional skewed "plain directional descends"
  if (1 && ATs) {
    const int cShear[][2] = {{-10,0},{-5,0},{0,-10},{0,-5}};
    for (int qskw=0;qskw<mcASz(cShear);qskw++) { BS_Cpy( cN, AConvex, ATs[0] ); // skew-buffer in "temps"
      BS7_SkewH( cShear[qskw][0], true, ATs[0], AW, AH ); BS7_SkewV( cShear[qskw][1], true, ATs[0], AW, AH );
      // ~"multiply single resulting frame in caller on the mask instead"
      if (1 && ADbgPfx) { sprintf( lDbgSfx, "B%da.raw", qskw ); BS_Write( ADbgPfx, cN, ATs[0] ); };
      BS_DataDup( 7, cN, ATs[0], ATs+1 ); BS_CPLA( icpl_DirDescend, 8, cBS_AllDirsT, AW, AH, ATs );
      for (int qi=0;qi<8;qi++) { // for each directional layer
        BS7_SkewV( -cShear[qskw][1], true, ATs[qi], AW, AH ); BS7_SkewH( -cShear[qskw][0], true, ATs[qi], AW, AH ); 
        tBS_Val *qDs[2] = {ARes[qi],ATs[qi]}; BS_CPB( mccpb_Mx<2>, 2|cCPBNDD, cN, qDs );
      }; }; // for i, for skw
  }; // if
  // BS_MovAve( 1, AW, AH, ARes, cBS_AllDirs, 8 ); // small orthogonal blur
  BS_MovAveX( 1, AW, AH, ARes[2] ); BS_MovAveX( 4, AW, AH, ARes[3] ); BS_MovAveX( 3, AW, AH, ARes[4] );
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "V" ); BS_DatasWrt( ADbgPfx, cN, ARes, 0, 8 ); }; }

// clavicle's estimation
static int __fastcall icpb_Clavicle( tBS_Val *ADat ) {
  tBS_Val qD[4]; for (int qi=0;qi<4;qi++) qD[qi] = abs(ADat[qi]-ADat[qi+4]);
  mcSortBAs<tBS_Val,8,2>( ADat ); ADat[0] = (ADat[0]+ADat[1])>>1;
  #if 1
    mcSortBDe<tBS_Val,4,2>( qD ); ADat[1] = mcSum<tBS_Val,2>(qD)>>1;
  #else
    ADat[1] = mcSum<tBS_Val,4>(qD)>>2;
  #endif
  return 2; }
void BS7c_Clavicle( tBS_Val **ADes, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ) { // -- clavicle's estimation
  const int cN = AW*AH; char *lDbgSfx = NULL, qDbgPfx[88]; // local debug prefix
  if (ADbgPfx) { strcpy( qDbgPfx, ADbgPfx ); ADbgPfx = qDbgPfx; lDbgSfx = strend(qDbgPfx); };
  BS_CPB( icpb_Clavicle, 8|(2<<8), cN, ADes, ARes );
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "V" ); BS_DatasWrt( ADbgPfx, cN, ARes, 0, 2  ); }; }

//======= narrow spatial filters ======
//=====================================
// splitter_mask by a spatial filter for maxmag, removed_ribs
// It is a mask of ~"mostly horz edges" (with other adhoc features; as this code requires).
struct tSbFWr { // - SplitByFilter's wrapper
  // -- intrinsics of SbFWr
  int  rW, rH, rIC, *rBuf[3]; // dimnesion, columns' buffers
  const tBS_Val *cShape, *cInpMsk; tBS_Val **rTs; // initial shape, input known <perhaps aliased> mask 0..512
  char rDbgPfx[80], *rDbgSfx; void iDbgW( const char *ANam, const tBS_Val *ADat ) { // debug file-write
    if (!rDbgSfx) return; sprintf( rDbgSfx, "%s.raw", ANam ); BS_Write( rDbgPfx, rW*rH, ADat ); }
  static void itppr_BandSttGet( tBS_SpawnThreadRec *AP ) {
    tSbFWr *cSelf = (tSbFWr*)AP->rP[0]; const tBS_Val *cDat = cSelf->rTs[0], *cInpMsk = cSelf->cInpMsk;
    int *qBufC0 = cSelf->rBuf[0], *qBufZ0 = cSelf->rBuf[1];
    const int cW = cSelf->rW, cLstRow = AP->rLst;
    for (int qr=AP->rFst;qr<cLstRow;qr++) { int qZ0 = 0, qC0 = 0, qO = qr*cW;
      for (int qc=0;qc<cW;qc++) {
        int qV = cDat[qO+qc]; if (!qV) continue;
        int qM = cInpMsk[qO+qc]; if (qM<256) { qZ0 += qV; qC0++; }; }; // for c
      qBufC0[qr] = qC0; qBufZ0[qr] = qZ0;
    }; } // for r
  static void itppr_BandSttAve( tBS_SpawnThreadRec *AP ) { tSbFWr *cSelf = (tSbFWr*)AP->rP[0];
    const int *cBufC0 = cSelf->rBuf[0], *cBufZ0 = cSelf->rBuf[1]; int *qBufA0 = cSelf->rBuf[2];
    const int cH = cSelf->rH, cLstRow = AP->rLst; for (int qr=AP->rFst;qr<cLstRow;qr++) { 
      const int cRad = cH/33+1;
      int qh0 = qr-cRad; if (qh0<0) qh0 = 0; int qh1 = qr+cRad; if (qh1>=cH) qh1 = cH-1;
      int qZ = 0, qC = 0; for (int qh=qh0;qh<=qh1;qh++) { qC += cBufC0[qh]; qZ += cBufZ0[qh]; };
      qBufA0[qr] = qC>0?qZ/qC:0;
    }; } // for r
  static int iAveThr( int AN, const tBS_Val *ADat, int AThr ) { int qZ = 0, qC = 0; for (int qi=0;qi<AN;qi++) {
    int qV = ADat[qi]; if (qV>AThr) { qZ += qV; qC++; }; }; return qZ /= qC+1; }

  public: // -- publics of SbFWr
  void Init( int AW, int AH, int AIC, const tBS_Val *AShape, const tBS_Val *AInpMsk,
    tBS_Val **ATs, const char *ADbgPfx ) { // -- initialize
    rW = AW; rH = AH; rIC = AIC; cShape = AShape, cInpMsk = AInpMsk; rTs = ATs;
    const int cBuf = BS_Align( AH ); rBuf[0]= new int[cBuf*mcASz(rBuf)];
    for (int qi=1;qi<mcASz(rBuf);qi++) rBuf[qi] = rBuf[qi-1]+cBuf;
    if (ADbgPfx) { strcpy( rDbgPfx, ADbgPfx ); rDbgSfx = strend(rDbgPfx); }; }
  void Release() { delete rBuf[0]; }; // -- release

  void SplitMskByFilter( float AKb, float AKw, tBS_Val **ARes ) { // -- determmine splitter_mask
    // step 1) calc Metric (i.e.) |Derivation(Shape)|.
    // step 2) calc local stat<s> for Metric (respectively thresholding input edges' mask).
    // step 3) equalize whole frame by "step 2". This code subtracs average Metric in flat|non_edge
    //         areas (where Mask is small).
    // step 4) invoke BS2_Cut and obtain splitter_mask.
    // Input edges' mask could be synthetic aliased 'mm_sbfdet' for MMAmbiToShrp, RRAmbiToShrp(InvCas=0)
    // or downsampled 'ribs_splt' from RRAmbiToShrp(InvCas=0) to RRAmbiToShrp(InvCas=1). Anyway it should
    // be suitable to step2.
    const int cW = rW, cH = rH, cN = cW*cH;
    switch (rIC) {
      case 0:  BS_Dr1( rW, rH, cShape, NULL, rTs[0], false ); // calc DerByY
               BS_MovAveY( cH/333+1, cW, cH, rTs[0] ); break;
      case 1:  BS_Cpy( cN, cShape, rTs[0] ); BS_DrA3x3( cW, cH, rTs[0] ); // calc DerByX,DerByY
               BS_MovAveY( cH/333+1, cW, cH, rTs[0] ); break;
      case 2:  BS_Cpy( cN, cShape, rTs[0] ); BS_Abs( cN, rTs[0] ); break; // calc Abs
      default: BS_Cpy( cN, cShape, rTs[0] ); }; // use inp Shape as-is
    if (0) iDbgW( "A", rTs[0] );
    BS_SpawnThreads( 0, 0, cH, itppr_BandSttGet, this );
    BS_SpawnThreads( 0, 0, cH, itppr_BandSttAve, this );
    tBS_Val *qDer = rTs[0];
    for (int qr=0;qr<cH;qr++) BS_VAdd( cW, qDer+cW*qr, -rBuf[2][qr] );
    if (0 && rDbgSfx) { tBS_Val *qT = rTs[1]; BS_Cpy( cN, rTs[0], qT ); // debug dump
      for (int qr=0;qr<cH;qr++) { tBS_Val *qrow = qT+cW*qr;
        for (int qi=0;qi<mcASz(rBuf);qi++) qrow[qi] = mcRg(rBuf[qi][qr]); };
       iDbgW( "B", qT ); };
    int qA0 = (int)(iAveThr(cN,qDer,0)*AKb),
        qA1 = (int)(iAveThr(cN,qDer,qA0)*AKw);
    tBS2_Cut cCut = { qA0, qA1, 512, {0,0} };
    BS_Cpy( cN, qDer, ARes[0] ); BS2_Cut( &cCut, ARes[0], cN ); // Res#0:msk
    if (1) iDbgW( "C", ARes[0] ); }
}; // end of SbFWr
static void iSplitMskByFilter( int AIC, const tBS_Val *AShape, const tBS_Val *AInpMsk, float AKb, float AKw,
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ) { // ambient_lap's splliter_msk (rather horz)
  // inp: ~ribs_ambient <to be removed>, initial mask 0..512 <aliased> to determine the
  // "split by filter; out#0:splitter_mask 0..512
  tSbFWr qSbF; memset( &qSbF, 0, sizeof(qSbF) ); qSbF.Init( AW, AH, AIC, AShape, AInpMsk, ATs, ADbgPfx );
  qSbF.SplitMskByFilter( AKb, AKw, ARes ); qSbF.Release(); }

// approach of Lap within Split_mask as ~"mx_dir_ave-ortho_dir_ave"; for BS_AllDirs!
template<int cD0,int cD1> static int __fastcall mccpb_ApproachQuad( tBS_Val *ADat ) { // quadrature(lap*split)
  int qMxA = ADat[cD0], qMxI = cD0; // for BS_AllDirs; <signed> quadrature of max dir_lap
  for (int qi=cD0+1;qi<cD1;qi++) { int qV = ADat[qi]; if (abs(qMxA)<abs(qV)) { qMxA = qV; qMxI = qi; }; };
  int qMxI0 = (qMxI-1)&7, qMxA0 = ADat[qMxI0], qMxI1 = (qMxI+1)&7, qMxA1 = ADat[qMxI1];
  if (abs(qMxA0-qMxA)>abs(qMxA1-qMxA)) { qMxA0 = qMxA1; qMxI0 = qMxI1; };
  return (qMxA-ADat[(qMxI+4)&7]+qMxA0-ADat[(qMxI0+4)&7])>>1; }
static void iApproach( int AIC, const tBS_Val *ALap, const tBS_Val *ASplit,
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // func(lap*split,InvocationCase)
  if (0 && ADbgPfx) { tBS_Val *qTs[] = {(tBS_Val*)ALap,(tBS_Val*)ASplit};
    strcpy( lDbgSfx, "A" ); BS_DatasWrt( ADbgPfx, cN, qTs, 0, mcASz(qTs) ); };
  BS_Cpy( cN, ALap, ARes[0] ); BS_Mul( cN, ASplit, ARes[0], 444 ); // #0:lap*split/nrm ??? 333
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "B.raw" ); BS_Write( ADbgPfx, cN, ARes[0] ); };
  if (AIC<0) return; // !!!!
  // calc <signed> quadrature of max directional laplacian (similar to MaxLap itself)
  static const struct { int rRad; tBS_CPBProc rProc; } cIP[3] = { // InvocationCases' parameters
    { 4, mccpb_ApproachQuad<2,7> },   // #0
    { 9, mccpb_ApproachQuad<0,8> },   // #1
    { 3, mccpb_ApproachQuad<0,8> } }; // #2
  BS_DataDup( 8, cN, ARes[0], ATs ); BS_MovAve( cIP[AIC].rRad, AW, AH, ATs, cBS_AllDirs );
  BS_CPB( cIP[AIC].rProc, cCPBNDD|8, cN, ATs, ARes[0] ); // #0:quadrature(lap*split)
  if (1 && ADbgPfx) { // BS_Cpy( cN, ALap, ARes[1] ); BS_SubA( cN, ARes[0], ARes[1] );
    strcpy( lDbgSfx, "C" ); BS_DatasWrt( ADbgPfx, cN, ARes, 0, 1 ); }; }

// approach ambient shape to exact sharp one
void BS7c_MMAmbiToShrp( const tBS_Val *AShrp, const tBS_Val *AAmbi, const tBS_Val *ASMSdet,
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- approach ambient shape to exact sharp one; "I.raw"
  // Shrp has big derX, derY (i.e. "Shrp~=LapDescend"). Res := Ambi+Msk*(Shrp-Ambi).
  // Msk's implementation assumes then Ambi has small derX, substantial derY (i.e.
  // "Ambi~=MaxMagnitudeOfRibsBasic").
  if (1 && ADbgPfx) { tBS_Val const *qTs[] = {AShrp,AAmbi,ASMSdet};
    strcpy( lDbgSfx, "Ainp" ); BS_DatasWrt( ADbgPfx, cN, (tBS_Val**)qTs, 0, mcASz(qTs) ); };
  BS_Cpy( cN, AAmbi, ATs[0] ); // TS#0:Ambi
  BS_Cpy( cN, AShrp, ATs[1] ); BS_SubA( cN, ATs[0], ATs[1] ); // TS#1:diff
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "B" ); BS_DatasWrt( ADbgPfx, cN, ATs, 0, 2 ); };
  BS_Cpy( cN, ASMSdet, ATs[2] ); // TS#2:splitter mask
  for (int qStep=0;qStep<3;qStep++) {
    tBS_Val **pTs = ATs+3;
    if (ADbgPfx) sprintf( lDbgSfx, "C%dAsms", qStep ); iSplitMskByFilter( 0, // Ambient_lap's horz_splliter_msk
      ATs[0], ATs[2], 0.5f, 0.7f, pTs, pTs+1, 1?ADbgPfx:NULL, AW, AH ); BS_Cpy( cN, pTs[0], ATs[2] );
    if (ADbgPfx) sprintf( lDbgSfx, "C%dBapr", qStep ); iApproach( 2, // func(diff*split,InvocationCase)
      ATs[1], ATs[2], pTs, pTs+1, 1?ADbgPfx:NULL, AW, AH ); // #1:approach
    BS_Add( cN, pTs[0], ATs[0] ); BS_SubA( cN, pTs[0], ATs[1] ); // TS#0:Ambi TS#1:Diff
    if (1 && ADbgPfx) { sprintf( lDbgSfx, "C%dV", qStep ); BS_DatasWrt( ADbgPfx, cN, ATs, 1, 2 ); };
  }; // for step
  tBS_Val **pNoi = ATs+2; // -- denoise
  BS_Noise( 6, AW, AH, ATs[0], pNoi[0], pNoi+1, 6 ); // #0:noise
  BS_Cpy( cN, ATs[0], pNoi[1] ); BS_Add( cN, pNoi[0], pNoi[1] ); // #0:MaxMag's approach
  if (0 && ADbgPfx) { strcpy( lDbgSfx, "Gden" ); BS_DatasWrt( ADbgPfx, cN, pNoi, 0, 2 ); };
  BS_Cpy( cN, pNoi[1], ARes[0] ); // Res#0:MaxMag's approach
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "V.raw" ); BS_Write( ADbgPfx, cN, ARes[0] ); }; } 

// approach ambient shape to <approx known> sharp one
void BS7c_RRAmbiToShrp( int AIC, const tBS_Val *AOrg, const tBS_Val *ARRib, const tBS_Val *ASplitAux, 
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ) { mcExportProcEntry; // -- approach ambient shape to <approx known> sharp one
  static const struct { int // InvocationCases' parameters
    rRadN,rTLgN,rKb100,rKw100,rRad2;} cIP[3] = {
    {   4,    6,    44,   100,    3},   // #0
    {  12,    7,     0,     0,    6} }; // #1
  if (AIC<0 || AIC>mcASz(cIP)) { ErrMsg( "BS7_RRAmbiToShrp(IC)" ); return; }
  // calc ribs_denoise, r<emoved>ribs_denoise
  tBS_Val **pRib = ATs;
  BS_Cpy( cN, AOrg, pRib[0] ); BS_SubA( cN, ARRib, pRib[0] ); // #0:ribs_org
  if (1 && ADbgPfx) { tBS_Val *qDs[] = {(tBS_Val*)AOrg,(tBS_Val*)ARRib,(tBS_Val*)ASplitAux,pRib[0]};
    strcpy( lDbgSfx, "Ainp" ); BS_DatasWrt( ADbgPfx, cN, qDs, 0, mcASz(qDs) ); };
  BS_Noise( cIP[AIC].rRadN, AW, AH, pRib[0], pRib[1], pRib+2, cIP[AIC].rTLgN ); BS_Add( cN, pRib[1], pRib[0] ); // #0:ribs_denoise
  BS_Cpy( cN, AOrg, pRib[1] ); BS_SubA( cN, pRib[0], pRib[1] ); // #1:rribs_denoise:=Org-ribs_deniose
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "Bden" ); BS_DatasWrt( ADbgPfx, cN, pRib, 0, 2 ); };
  // calc splitter mask from ribs_denoise (unless one is nt supplied in
  tBS_Val **pMsk = pRib+2;
  if (cIP[AIC].rKb100<1) // AIC==1
    BS_Cpy( cN, ASplitAux, pMsk[0] ); // SplitAux is the mask "as is"
  else { if (ADbgPfx) strcpy( lDbgSfx, "Csmf" ); iSplitMskByFilter( 1, // SplitAux is <an appoximate> mm_sbf_det
    pRib[0], ASplitAux, cIP[AIC].rKb100*0.01f, cIP[AIC].rKw100*0.01f, pMsk, pMsk+1, 1?ADbgPfx:NULL, AW, AH ); };
  // calc ~small laplacian of rribs_denoise
  tBS_Val **pLap = pMsk+1;
  BS2_MovAve2D( cIP[AIC].rRad2, pRib[1], pLap[0], AW, AH ); BS_SubS( cN, pRib[1], pLap[0] ); // #0:lap
  tBS_Val **pApr = pLap+1;
  // if ((AIC==1) BS_ZFil( cN, pApr[0] ); else ... ???
  if (ADbgPfx) strcpy( lDbgSfx, "E" ); iApproach( AIC, // func(lap*split,InvocationCase)
    pLap[0], pMsk[0], pApr, pApr+1, ADbgPfx, AW, AH );
  BS_Med( 1, AW, AH, pApr[0], pApr[1] ); BS_MovAveX( 1, AW, AH, pApr[1] ); BS_Cpy( cN, pApr[1], pApr[0] ); 
  // res := rribs_denoise-approach
  BS_Cpy( cN, pRib[1], ARes[0] ); BS_SubA( cN, pApr[0], ARes[0] );
  if (AIC==0)
    BS_Cpy( cN, pMsk[0], ARes[1] );
  if (1 && ADbgPfx) { strcpy( lDbgSfx, "V" ); BS_DatasWrt( ADbgPfx, cN, ARes, 0, 1  ); }; }

//======= helpers for ProcessLftRgtLung ======
//============================================
static int icpl_CrdRgtTopToLftBtmX( int AL, tBS_Val *ADat ) {
  const int cSpan = 4, cL = AL-cSpan;
  BS_SubA( cL, ADat+cSpan, ADat ); BS_Abs( cL, ADat );
  int qZZ = mcSum( cL, ADat ); if (qZZ<11) { const int cScLg = 20, cSc = (1<<cScLg)/cL;
    for (int qi=0;qi<cL;qi++) ADat[qi] = ((cL-qi)*cSc)>>(cScLg-12); return AL; };
  int qZ = 0; for (int qi=cL-1;qi>=0;qi--) { qZ += ADat[qi]; ADat[qi] = (qZ<<12)/qZZ; };
  BS_ZFil( cSpan, ADat+cL ); return AL; }
static int icpl_CrdRgtTopToLftBtmY( int AL, tBS_Val *ADat ) {
  AL = (int)(AL*0.4f); const int cScLg = 15, cSc = (1<<cScLg)/AL;
  for (int qi=0;qi<AL;qi++) { ADat[qi] = (int)(ADat[qi]*qi*cSc)>>cScLg; }; return AL; }
void BS7c_CrdRgtTopToLftBtm( tBS_Val *ADat, int AW, int AH ) { // -- empirical ascend RgtToLft within a lung_area
  const int cN = AW*AH;
  BS_CPLCRows( icpl_CrdRgtTopToLftBtmX, AW, AH, ADat, true );
  BS_CPLCCols( icpl_CrdRgtTopToLftBtmY, AW, AH, ADat );
  BS_ZFil( AW*(int)(AH*0.17f), ADat );
  BS_VMul( cN, ADat, 1.4f );
  BS_MovAveY( 33, AW, AH, ADat ); BS_Rng( cN, ADat, -99999, (1<<12)-1); }

// normalize convex by histogram
void BS7c_CrdOfFastDescend( float ACutDrk, float ACutBrt, bool AInv, tBS_Val *ACrd, int AW, int AH ) { // -- normalize convex by histogram
  const int cN = AW*AH; BS_Neg( cN, ACrd );
  // refer ~iHistEqu() in BS_Ut4.cpp
  const int cGam = 1<<14, cCrdMx = (1<<12)-1;
  int *qHist = (int*)BS_GamCLAHEHists( AW, AH, ACrd, 0, cGam ),
       qNP0 = cN-qHist[0];
  if (!qNP0) { ErrMsg( "BS7_CrdOfFastDescend. Only zeros in histogram" ); return; };
  int  qMn = mcHistNP<int>( (int)(qNP0*ACutDrk), qHist+1 )+1,
       qMx = cGam+mcHistNP<int>( -(int)(qNP0*ACutBrt), qHist+cGam-1 ),
       qNP1 = mcSum( qMx-qMn, qHist+qMn ), qZ = 0;
  for (int qi0=qMn;qi0<qMx;qi0++) { qZ += qHist[qi0]; qHist[qi0] = cCrdMx*qZ/qNP1; };
  memset( qHist, 0, qMn*sizeof(int) ); for (int qi1=qMx;qi1<cGam;qi1++) qHist[qi1] = cCrdMx;
  if (AInv) for (int qi=0;qi<cGam;qi++) { int qV = cCrdMx-qHist[qi]; qHist[qi] = __max(qV,0); };
  BS_GamCLAHEApply( AW, AH, ACrd, 0, cGam, qHist );
  delete qHist; }

// simple marginal connectivities
static void iiMgnConSimple( tBS_Val *ALap, int AI0, int AI1 ) { int qInc = AI0>AI1?-1:1, qZ = 0;
  for (int qi=AI0+qInc;qi!=AI1;qi+=qInc) { int l = ALap[qi]; qZ += l; ALap[qi] = l-qZ/(abs(AI0-qi)+1); }; } 
static const int cMgnConNone = 0x8000;
static int itcpl_MgnConSimple( int AL, tBS_Val *ADat ) { const int cLMn = 10;
  int qI0, qI1; tBS_Val *qDat;
  if (AL<cLMn) goto labZFil; ADat[AL] = (tBS_Val)cMgnConNone;
  qI0 = 0; while (!ADat[qI0]) qI0++; if (qI0>=AL-cLMn) goto labZFil;
  qI1 = AL-1; while (!ADat[qI1]) qI1--; if (qI1-qI0<cLMn) goto labZFil;
  qDat = ADat+BS_Align( AL ); BS_Cpy( AL, ADat, qDat );
  iiMgnConSimple( ADat, qI0, qI1+1 ); iiMgnConSimple( qDat, qI1, qI0-1 );
  for (int qi=qI0;qi<=qI1;qi++) { int la1 = qDat[qi]; if (abs(ADat[qi])>abs(la1)) ADat[qi] = la1; };
  BS_VFil( qI0, ADat, cMgnConNone ); BS_VFil( AL-qI1-1, ADat+qI1, cMgnConNone ); return AL;
  labZFil: BS_VFil( AL, ADat, cMgnConNone ); return AL; }
void BS7c_MgnConSimple( const tBS_Val *AOrg, bool ASkews, tBS_Val **ARes,
  tBS_Val **ATs, int AW, int AH ) { int cN = AW*AH; // -- marginal connectivities
  const int cShearXY[][2] = {{-12,-12},{12,12}}, cskw = ASkews ? 1+mcASz(cShearXY) : 1;
  for (int qskw=0;qskw<cskw;qskw++) { const int t = __min(qskw,1); BS_Cpy( cN, AOrg, ATs[t] );
    if (t) { BS7_SkewH( cShearXY[qskw-1][0], false, ATs[t], AW, AH ); BS7_SkewV( cShearXY[qskw-1][1], false, ATs[t], AW, AH ); };
    BS_DataDup( 7, cN, ATs[t], ATs+t+1 );
    BS_CPLA( itcpl_MgnConSimple, 8, cBS_AllDirsT, AW, AH, ATs+t );
    BS_CPB( mccpb_MnA<8>, 8|cCPBNDD, cN, ATs+t );
    if (t) { BS7_SkewV(-cShearXY[qskw-1][1], false, ATs[t], AW, AH ); BS7_SkewH(-cShearXY[qskw-1][0], false, ATs[t], AW, AH );
              BS_CPB( mccpb_MnA<2>, 2|cCPBNDD, cN, ATs ); }; }; // for skw
  BS_Cpy( cN, ATs[0], ARes[0] ); tBS_Val *qRes = ARes[0];
  for (int qi=0;qi<cN;qi++) if (qRes[qi]==(tBS_Val)cMgnConNone) qRes[qi] = 0; }

//================ in development ============
//============================================
