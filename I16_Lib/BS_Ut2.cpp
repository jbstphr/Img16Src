// I16_Ut2.cpp : CPU mid-level routines for certain dimension
// the routines here use fixed frames BS2_Dim=>BS2_W,BS2_H
// by S.B.
#include "stdafx.h"
#include "math.h"
#include "BS_Ut0.h"
#include "BS_Ut2.h"
#include "BS_Ut4.h"
#include "BS_Ut7.h"

//================= commons ==================
//============================================
static int xW, xH; // standard frame-size for local operations
void BS2_Dim( int AW, int AH ) { xW = AW; xH = AH; } // -- assigns standard frame's dimension
int BS2_W() { return xW; } int BS2_H() { return xH; }

void BS2_MimicLap( const tBS_Val *ADat, tBS_Val *ALap, int AN ) { if (!AN) AN = xW*xH; // -- mimic |Dat| and sign(Lap)
  for (int qi=0;qi<AN;qi++) { int qL = ALap[qi], qV = ADat[qi]; if (qL<0) qV = -qV; else if (!qL) qV = 0; ALap[qi] = qV; }; }
void BS2_MimicDat( const tBS_Val *ALap, tBS_Val *ADat, int AN ) { if (!AN) AN = xW*xH; // -- mimic sign(Lap) and |Dat|
  for (int qi=0;qi<AN;qi++) { int qL = ALap[qi]; if (qL>0) continue; ADat[qi] = !qL ? 0 : -ADat[qi]; }; }
void BS2_VMimicNZ( tBS_Val *ADat, int AV, int AN ) { if (!AN) AN = xW*xH; // set pix<>0 to V
  for (int qi=0;qi<AN;qi++) if (ADat[qi]) ADat[qi] = AV; }

void BS2_MovAve2D( int ARad, const tBS_Val *ASrc, tBS_Val *ADst, int AW, int AH ) { // -- moving average
  if (!AW) { AW = xW; AH = xH; }; const int cN = AW*AH; switch (ARad) {
    case 1: BS_Ave3x3( AW, AH, ASrc, ADst ); break; case 2: BS_Ave5x5( AW, AH, ASrc, ADst ); break;
    default: BS_Cpy( cN, ASrc, ADst ); if (ARad) BS_MovAve2D( ARad, AW, AH, ADst ); }; }

void BS2_MovAve2Dp( int ARad, tBS_Val *ADat, tBS_Val **ATs1, int Apyr, int AW, int AH ) { // -- 2D average within downsample
  if (ARad<=0) return; if (!AW) { AW = xW; AH = xH; };
  if (Apyr<1) { BS2_MovAve2D( ARad, ADat, ATs1[0], AW, AH ); BS_Cpy( AW*AH, ATs1[0], ADat ); return; };
  const int cW2 = AW>>Apyr, cH2 = AH>>Apyr; BS_ScaleDn( Apyr, AW, AH, ADat, ATs1[0] ); BS_MovAve2D( ARad, cW2, cH2, ADat );
  BS_ScaleUp( Apyr, cW2, cH2, ADat, ATs1[0], true ); } // _DimUpC reduces lines' aliases

void BS2_Hess( int ARad, tBS_Val **AHess, int AW, int AH ) { // -- inp#0:Lap, out#0:Hess #1:Lap_blr
  if (!AW) { AW = xW; AH = xH; }; const int cN = AW*AH; 
  BS2_MovAve2D( ARad, AHess[0], AHess[1], AW, AH ); BS_Abs( cN, AHess[0] ); BS2_MovAve2D( ARad, AHess[0], AHess[2], AW, AH ); 
  BS_Cpy( cN, AHess[1], AHess[0] ); BS_Abs( cN, AHess[0] ); BS_SubS( cN, AHess[2], AHess[0] ); }

bool BS2_Sum( float K, tBS_Val *A, float K1, const tBS_Val *A1, // -- summation without overflow
  float K2, const tBS_Val *A2, float K3, const tBS_Val *A3, float K4, const tBS_Val *A4 ) {
  if ((fabs(K-1)<1e-3f) && (fabs(K1)<1e-3f) && (fabs(K2)<1e-3f) && (fabs(K3)<1e-3f) && (fabs(K4)<1e-3f)) return false;
  const int N = xW*xH; int *I = new int[N]; memset( I, 0, N*sizeof(int) ); BS_AddMul( N, A, I, K );
  if (A1 && fabs(K1)>1e-3f) BS_AddMul( N, A1, I, K1 ); if (A2 && fabs(K2)>1e-3f) BS_AddMul( N, A2, I, K2 );
  if (A3 && fabs(K3)>1e-3f) BS_AddMul( N, A3, I, K3 ); if (A4 && fabs(K4)>1e-3f) BS_AddMul( N, A4, I, K4 );
  for (int i=0;i<N;i++) A[i] = mcRg( I[i] ); delete I; return true; }

bool BS2_MskSum( tBS_Val *Dst, tBS_Val *T, int Nrm, float K0, tBS_Val *A0, // -- sum normalized  masks
  float K1, const tBS_Val *A1, float K2, const tBS_Val *A2, float K3, const tBS_Val *A3 ) {
  const int N = xW*xH; float KZ = 0; BS_ZFil( N, Dst );
  if (A0 && fabs(K0)>1e-3f) { KZ += K0; BS_Cpy( N, A0, T ); BS_VMul( N, T, K0 ); BS_Add( N, T, Dst ); }
  if (A1 && fabs(K1)>1e-3f) { KZ += K1; BS_Cpy( N, A1, T ); BS_VMul( N, T, K1 ); BS_Add( N, T, Dst ); }
  if (A2 && fabs(K2)>1e-3f) { KZ += K2; BS_Cpy( N, A2, T ); BS_VMul( N, T, K2 ); BS_Add( N, T, Dst ); }
  if (A3 && fabs(K3)>1e-3f) { KZ += K3; BS_Cpy( N, A3, T ); BS_VMul( N, T, K3 ); BS_Add( N, T, Dst ); }
  if (KZ<1e-3f) { BS_VFil( N, Dst, Nrm ); return false; }; int qAdd = (int)((1-KZ)*Nrm);
  if (qAdd) BS_VAdd( N, Dst, qAdd ); return true; }

bool BS2_AccMim( tBS_Val *Acc, tBS_Val *T, float K0, const tBS_Val *E0, // -- add mimiclap to accumulator
  const tBS_Val *L0, float K1, const tBS_Val *E1, const tBS_Val *L1 ) { 
  const int N = xW*xH; const bool b0 = K0>1e-3f, b1 = K1>1e-3f; if (!(b0||b1)) return false;
  if (b0) { BS_Cpy( N, L0, T ); BS2_MimicLap( E0, T ); BS_VMul( N, T, K0 ); BS_Add( N, T, Acc ); };
  if (b1) { BS_Cpy( N, L1, T ); BS2_MimicLap( E1, T ); BS_VMul( N, T, K1 ); BS_Add( N, T, Acc ); }; return true; }

static void ippr_RngLapByConductor( tBS_SpawnThreadRec *AP ) {
  tBS_Val *qCnd = (tBS_Val*)AP->rP[0], *qLap = (tBS_Val*)AP->rP[1];
  const int cLst = AP->rLst; for (int qi=AP->rFst;qi<cLst;qi++) {
  int qC = qCnd[qi], qL = qLap[qi]; if (qC<abs(qL)) qLap[qi] = mcsign( qC, qL ); }; }
void BS2_RngLapByConductor( const tBS_Val *ACnd, tBS_Val *ALap ) { // -- constraint lap by cardinal conductor
  BS_SpawnThreads( 0, 0, xW*xH, ippr_RngLapByConductor, ACnd, ALap ); }

static void ippr_RidgeByConductors( tBS_SpawnThreadRec *AP ) {
  tBS_Val *qRdg = (tBS_Val*)AP->rP[0], *cCDev = (tBS_Val*)AP->rP[1],
          *cCDer0 = (tBS_Val*)AP->rP[2], *cCDer1 = (tBS_Val*)AP->rP[3];
  const int cLst = AP->rLst; for (int qi=AP->rFst;qi<cLst;qi++) {
    int cDev = cCDev[qi], cDer = __max(cCDer0[qi],cCDer1[qi]), qV = cDev-cDer;
    qRdg[qi] = qV<=0 ? 0: mcsign(qV,qRdg[qi]); }; }
void BS2_RidgeByConductors( tBS_Val *ARdg, const tBS_Val *ACDev, // -- mimic a ridge by lap & conductors
  const tBS_Val *ACDer0, const tBS_Val *ACDer1 ) { // (the conductors shall be normalized)
  BS_SpawnThreads( 0, 0, xW*xH, ippr_RidgeByConductors, ARdg, ACDev, ACDer0, ACDer1 ); }

// MxLap with extended <somehow> Limbo
static int ixMxLapXLimboK, ixMxLapXLimboTolGap;
template<int cNDDLg> static int __fastcall mccpb_MxLapXLimbo( tBS_Val *ADat ) { const int cNDD = 1<<cNDDLg;
  int qM = ADat[2], qA = qM; for (int qi=1;qi<cNDD;qi++) { int v = ADat[2+qi]; if (abs(qM)<abs(v)) qM = v; qA += v; };
  qA >>= cNDDLg; const int cK = ixMxLapXLimboK, cGap = ixMxLapXLimboTolGap;
  // calc virtual low_est, high_est relative to the basis (which is MxLapStd)
  // always "Ave<MxLapStd<MxLapExact" or "MxLapExact<MxLapStd<Ave"
  int qMstd = ADat[1]; qM = qMstd+((qM-qMstd)*cK>>10); qA = qMstd+((qA-qMstd)*cK>>10); int qLap = ADat[0];
  // low_est:=min(MxLapStd,Ave), high_est:=max(MxLapStd,Ave)
  return qM<qA ? mcRng(qLap,qM-cGap,qA+cGap) : mcRng(qLap,qA-cGap,qM+cGap); }
template<int cNDDLg> static int __fastcall mccpb_AutOscMsk( tBS_Val *ADat ) {
  const int cNDD2 = 1<<(cNDDLg-1); int qZ = 0; for (int qi=0;qi<cNDD2;qi++) {
    int q0 = ADat[qi], q1 = ADat[qi+cNDD2]; qZ += abs(q0)+abs(q1)+abs(q0-q1); }; return qZ>>(cNDDLg-1); }
void BS2_MxLapXLimbo( int ARad, float AK, int ATolGap, int AZ0, float AKAO,
  tBS_Val *ALap, tBS_Val **ATs, int AW, int AH ) { // -- MxLapStd with ~extended_limbo
  // called in BS_Ut5.cpp,_UnSharpCPU{"RadMM,KMM,TGMM,ZMM,KAO"}; BS_Ut6.cpp,_DePoisson{"1,KU,0,0,0.0f"}
  if (ARad<=0) return; if (!AW) { AW = xW; AH = xH; }; const int cN = AW*AH;
  BS_Cpy( cN, ALap, ATs[0] ); int cNDD;
  // disable BS2A_Shear here (as use byproduct directional averages in Ts[])
  if (ATolGap<0) { // Lap+=K*MxLapStd
    cNDD = BS_MxLapStd( ARad, AW, AH, ATs[0], ATs+1, 1.0f, AZ0, 0, 0.0f );
    if (AK>1e-3f){ BS_VMul( cN, ATs[0], AK ); BS_Add( cN, ATs[0], ALap ); };
  } else { // Lap:=Rng(Lap,low_est-TolGap,high_est+TolGap), _est:=K*func(MxLapStd,MxLapExact,Ave)
    cNDD = BS_MxLapStd( ARad, AW, AH, ATs[0], ATs+1, 1.0f, AZ0, 0, 0.0f );
    tBS_Val *qDs[10] = {ALap}; memcpy( qDs+1, ATs, (cNDD+1)*sizeof(void*) );
    ixMxLapXLimboK = (int)(AK*1024); ixMxLapXLimboTolGap = ATolGap;
    BS_CPB( cNDD==4?mccpb_MxLapXLimbo<2>:mccpb_MxLapXLimbo<3>, (cNDD+2)|cCPBNDD, cN, qDs );
  };
  // KA0<>0:AutoOscillator (<0:copies dbgtrc's frames to Ts^)
  bool qAODbgInf = false; if (AKAO<-1e-3f) { AKAO = -AKAO; qAODbgInf = true; };
  if (AKAO>1e-3f) { BS_CPB( cNDD==4?mccpb_AutOscMsk<2>:mccpb_AutOscMsk<3>, cNDD|cCPBNDD, cN, ATs ); BS_VMul( cN, ATs[0], AKAO );
    if (qAODbgInf) { BS_Cpy( cN, ATs[0], ATs[1] ); BS_Cpy( cN, ALap, ATs[2] ); };
    const int cTg0 = 128; const tBS2_Cut cCut = {{0,512,512-cTg0},{0,0}};
    BS2_Cut( &cCut, ATs[0] ); BS_VAdd( cN, ATs[0], cTg0 ); BS_Mul( cN, ATs[0], ALap, -9 );
  }; } // of AutoOscillator

/*
// -- |Counts| within continuous pos|neg intervals
static int icpl_CntsMx( int AL, tBS_Val *ADat ) { mcIvlC( AL, ADat ); return AL; }
static int __fastcall icpb_CntsMx( tBS_Val *ADat ) { if (!ADat[0]) return 0; // Sum of <2> max Counts
  BS_Abs( 8, ADat ); mcSortBDe<tBS_Val,8,2>( ADat ); return ADat[0]+ADat[1]; }
void BS2_CntsMxInternal( int AW, int AH, tBS_Val **ADat8 ) { const int cN = AW*AH; // -- inp|out in Dat#0
  mcShiftRgt( cN, ADat8[0], 2 ); BS_DataDup( 7, cN, ADat8[0], ADat8+1 ); 
  BS_CPLA( icpl_CntsMx, 8, cBS_AllDirsT, AW, AH, ADat8 ); BS_CPB( icpb_CntsMx, 8|cCPBNDD, cN, ADat8 ); }
static void iCntsMxInterim( int AW, int AH, tBS_Val **ADat8 ) { BS2_CntsMxInternal( AW, AH, ADat8 ); }
  //BS_Ave5x5( AW, AH, ADat8[0], ADat8[1] ); BS_Cpy( AW*AH, ADat8[1], ADat8[0] ); }
void BS2_CntsMx( const tBS_Val *ALap, tBS_Val **ATs, int Apyr, int AW, int AH ) { // -- inp:Lap, out:Sum of<~2 mx |counts|
  if (!AW) { AW = xW; AH = xH; }; BS_Cpy( AW*AH, ALap, ATs[0] ); if (Apyr<1) { iCntsMxInterim( AW, AH, ATs ); return; };
  BS_ScaleDn( Apyr, AW, AH, ATs[0], ATs[1], true ); const int cW2 = AW>>Apyr, cH2 = AH>>Apyr; iCntsMxInterim( cW2, cH2, ATs );
  mcShiftLft( cW2*cH2, ATs[0], Apyr ); BS_ScaleUp( Apyr, cW2, cH2, ATs[0], ATs[1], false ); }
*/

//========= Gamut's transformations ==========
//============================================
// cut & rescale
int *BS2_initCutLat( const tBS2_Cut *AP, int *AMM ) { // -- for CBS_XLat::Init(...)
  int qIM = AP->n[cBS2_rcuBlack], qIGa1 = AP->n[cBS2_rcuWhite], qOGa1 = AP->n[cBS2_rcuGaO]-1;
  const bool cNeg = qIM>qIGa1; if (cNeg) mcSwap( int, qIM, qIGa1 ); qIGa1 -= qIM; AMM[0] = qIM; AMM[1] = qIGa1+1;
  if (qIGa1<0) { ErrMsg( "BS2_Cut wrong params" ); return NULL; }; int *qLat = new int[qIGa1+1];
  if (cNeg) { qLat[0] = qOGa1; qLat[qIGa1] = 0; } else { qLat[0] = 0; qLat[qIGa1] = qOGa1; }; mcIplIvl<int>( qLat, qIGa1 ); return qLat; }
int BS2_initCutTD( float ATD ) { return (int)(ATD<0.99999f?-cBS_XLatTD*fabs(ATD):ATD); } // -- float => TD0,TD1
void BS2_Cut( const tBS2_Cut *AP, tBS_Val *ADat, int AN ) { // -- cut & rescale 
  int qMG[2], *qLat = BS2_initCutLat( AP, qMG ); if (!qLat) return; CBS_XLat qXLat;
  qXLat.Init( qLat, qMG[0], qMG[1], BS2_initCutTD(AP->f[cBS2_rcfTD0]), BS2_initCutTD(AP->f[cBS2_rcfTD1]) );
  qXLat.Apply( AN?AN:xW*xH, ADat ); delete qLat; }

// histogram truncation
void iiHistMM( int AGam, int ANP, const unsigned *AHist, float ACutDrk, float ACutBrt, int *AMM ) { // minmax_range by histogram
  int qCutDrk = ACutDrk<1?(int)(ANP*ACutDrk):(int)ACutDrk, qMn = mcHistNP<unsigned>( __max(qCutDrk,1), AHist );
  int qCutBrt = ACutBrt<1?(int)(ANP*ACutBrt):(int)ACutBrt, qMx = AGam+mcHistNP<unsigned>( -__max(qCutBrt,1), AHist+AGam-1 );
  AMM[0] = qMn; AMM[1] = qMx; } 
int xBS2_rhCatchMM[2]; // minmax_range by histogram
void BS2_HistCut( const tBS2_HistCut *AP, tBS_Val *ADat, const tBS_Val *ARef, tBS_Val **ATs2, int AW, int AH ) { // -- cut by histogram
  int qCatchMM[2] = {0,0}; const bool cDefaultThread = !AW; if (cDefaultThread) { AW = xW; AH = xH; };
  const int cMnI = AP->n[cBS2_rhcBlack], cMxI = AP->n[cBS2_rhcWhite], cGaI1 = abs(cMxI-cMnI),
            cMgn = AP->n[cBS2_rhcMgn], cGaO1 = AP->n[cBS2_rhcGaO]-1, cN = AW*AH; // params
  //if (1) BS_Write( "$yaa.raw", cN, ADat );
  if (cMnI<cMxI) { if (cMnI) BS_VAdd( cN, ADat, -cMnI ); } else { // negative input
    if (cMxI) BS_VAdd( cN, ADat, -cMxI ); BS_Neg( cN, ADat ); BS_VAdd( cN, ADat, cGaI1 ); BS_Rng( cN, ADat, 0, cGaI1 ); };
  //if (1) BS_Write( "$yab.raw", cN, ADat );
  const int cw = AW-cMgn*2, ch = AH-cMgn*2; unsigned *qHist = new unsigned[cGaI1+1];
  BS_GamHist( cw, ch, ADat+cMgn*AW+cMgn, AW, qHist, cGaI1+1, ARef ); // get histogram
  iiHistMM( cGaI1+1, cw*ch, qHist, AP->f[cBS2_rhfCut0], AP->f[cBS2_rhfCut1], qCatchMM ); // minmax_range by histogram
  const tBS2_Cut cC = {qCatchMM[0],qCatchMM[1],cGaO1+1,AP->f[cBS2_rhfTD0],AP->f[cBS2_rhfTD1]};
  BS2_Cut( &cC, ADat, cN ); delete qHist; // apply range
  if (cDefaultThread) memcpy( xBS2_rhCatchMM, qCatchMM, sizeof(qCatchMM) ); }

// gamma-correction
void BS2_GammaCrv( float AGamma, int AGam, unsigned *ACrv, int AGamStep ) { // -- gamma-curve for specified gamut
  if (AGamStep<=0) { AGamStep = AGam/32; AGamStep = mcRng(AGamStep,8,256); };
  const int cGa1 = AGam-1; int qi0 = 0; ACrv[0] = 0;
  for (int qi=AGamStep;;qi+=AGamStep) { if (qi>=cGa1) qi = cGa1;
    ACrv[qi] = (int)(pow((float)qi/cGa1,AGamma)*cGa1+0.49999f ); mcIplIvl( ACrv+qi0, qi-qi0 );
    if (qi>=cGa1) break; qi0 = qi; }; } // for i

void BS2_Gamma( float AGamma, int AGam, tBS_Val *ADat, int AGamStep, int AN ) { // -- apply gamma-curve
  if (AGamma<1e-3f || fabs(1-AGamma)<1e-3f) return;
  int *qCrv = new int[AGam]; BS2_GammaCrv( AGamma, AGam, (unsigned*)qCrv, AGamStep );
  CBS_XLat qXLat; qXLat.Init( qCrv, 0, AGam, -2, -2 );
  qXLat.Apply( AN?AN:xW*xH, ADat ); delete qCrv; }

//=========== Filters (conscious) ============
//============================================
static void iBS2_MovAve8( const tBS_Val *ALapOrg, int ARad, int ARadO, tBS_Val **ADs ) { // moving average by 8 dirs
  const int cW = xW, cH = xH, cN = cW*cH; 
  BS_DataDup( 8, cN, ALapOrg, ADs ); BS_MovAve( ARad, cW, cH, ADs ); // standard average by 8 dirs
  if (ARadO<=0) return; tBS_Val **qOh = ADs+8; // subtract orthogonal averages
  for (int qi1=0;qi1<8;qi1++) BS_Cpy( cN, ADs[qi1], qOh[(qi1+4)&7] ); BS_MovAve( ARadO, cW, cH, qOh );
  for (int qi2=0;qi2<8;qi2++) BS_SubA( cN, qOh[(qi2+4)&7], ADs[qi2] ); }

static void iBS2_QuadLaps( tBS_Val **ADs, tBS_Val **AQuadLaps, int ARadQ, int ARadQO ) { // quad-laps
  const int cW = xW, cH = xH, cN = cW*cH;
  for (int qe=0;qe<4;qe++) { 
    BS_Cpy( cN, ADs[qe], AQuadLaps[qe] ); BS_SubA( cN, ADs[4+qe], AQuadLaps[qe] ); };
  if (ARadQ<=0) return;
  for (int qe1=0;qe1<4;qe1++) BS_Cpy( cN, AQuadLaps[qe1], AQuadLaps[4+qe1] ); 
  BS_MovAve( ARadQO, cW, cH, AQuadLaps );
  for (int qe2=0;qe2<4;qe2++) mcSwap( tBS_Val*, AQuadLaps[qe2], AQuadLaps[4+qe2] ); 
  BS_MovAve( ARadQ, cW, cH, AQuadLaps ); }

template <int cR> static int icpl_RadonENmLap( int AL, tBS_Val *ADat ) { // edge-naive-mimic sum lap
  if (AL<cR*2+1) { BS_ZFil( AL, ADat); return AL; };
  tBS_Val *qDat = ADat+BS_Align( AL ); for (int qi=cR;qi<AL-cR;qi++)
  qDat[qi] = mcsign( abs(ADat[qi-cR]-ADat[qi+cR]), ADat[qi] ); BS_Cpy( AL-cR*2, qDat+cR, ADat+cR );
  BS_VFil( cR, ADat, ADat[cR] ); BS_VFil( cR, ADat+AL-cR, ADat[AL-cR-1] ); return AL; }
static void iBS2_ENmLaps( tBS_Val **ALap8, tBS_Val *AENm, int ARad ) { // edge-naive-mimic-laps
  const int cW = xW, cH = xH, cN = cW*cH; tBS_Val **qTs = ALap8+8;
  for (int qi=0;qi<8;qi++) BS_Cpy( cN, ALap8[qi], qTs[qi] );
  tBS_CPLProc qProc; switch (ARad) {
    case 2: qProc = icpl_RadonENmLap<2>; break; case 3: qProc = icpl_RadonENmLap<3>; break;
    case 4: qProc = icpl_RadonENmLap<4>; break; case 5: qProc = icpl_RadonENmLap<5>; break;
    default: qProc = icpl_RadonENmLap<1>; };
  BS_CPLA( qProc, 8, cBS_AllDirsT, cW, cH, qTs );
  BS_CPB( mccpb_Sum<8>, cCPBNDD|8, cN, qTs, AENm ); }
  //BS_Write( "q1.raw", cN, ALap8[0] ); BS_Write( "q2.raw", cN, qTs[0] ); BS_Write( "q3.raw", cN, AENm );

#define mcRadonLd( I0, I1 ) { qV0 = ADat[I0]; qV1 = ADat[I1]; }
#define mcRadonMM( V ) { if (qMn>V) qMn = V; else if (qMx<V) qMx = V; }
#define mcRadonShann() { int d = abs(qV0-qV1); if (d>255) d = 255; qD += d; qZ += (qV0+qV1)*d; }
static int __fastcall cpb_RadonShannon( tBS_Val *ADat ) {
  int qV0, qV1, qMn, qMx, qD = 0, qZ = 0;
  mcRadonLd( 0, 4 ); qMn = qMx = qV0;  mcRadonMM( qV1 ); mcRadonShann();
  mcRadonLd( 1, 5 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonShann();
  mcRadonLd( 2, 6 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonShann();
  mcRadonLd( 3, 7 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonShann();
  ADat[0] = qD; ADat[1] = qZ/((qD<<1)+1);
  const int cO = ADat[8]; if (abs(qMx)<abs(qMn)) qMx = qMn;
  ADat[2] = abs(qMx)<abs(cO) ? (qMx*3+cO)>>2 : (qMx+cO)>>1; return 3; }
void BS2_RadonMM( const tBS2_RadonMM *AP, const tBS_Val *ALapOrg, 
  tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp, tBS_Val *AENm,
  tBS_Val **ATs, tBS_Val **AQuadLaps ) { // -- core stats of rough Radon 
  const int cN = xW*xH;
  iBS2_MovAve8( ALapOrg, AP->n[cBS2_rrcRad], AP->n[cBS2_rrcRadO], ATs );
  if (AQuadLaps) // <quite optional> for consequent operations
    iBS2_QuadLaps( ATs, AQuadLaps, AP->n[cBS2_rrcRadQ], AP->n[cBS2_rrcRadQO] );
  if (AENm) // <optional> addition to ~Srp
    iBS2_ENmLaps( ATs, AENm, AP->n[cBS2_rrcRadENm] );
  tBS_Val *qTs[9]; memcpy( qTs, ATs, sizeof(void*)*8 ); qTs[8] = (tBS_Val*)ALapOrg;
  tBS_Val *qOut[3] = {ADev,ABlr,ASrp};
  BS_CPB( cpb_RadonShannon, (3<<8)|9, cN, qTs, qOut ); }

//============= Filters (ad-hoc) =============
//============================================

//=============== Metrics ====================
//============================================
// metrics' line_shaders NxN
// If one utilizes "the moments" then <for relevant metrics> one might add the following code-line:
//   if (qLa&1) ADat[qi-qLa--] = 0; // topology-quirk to coincide with counts' moments 
static int icpl_IvlC( int AL, tBS_Val *ADat ) { mcIvlC( AL, ADat ); return AL; } // -- counts within continuous pos|neg ivls
static int icpl_IvlCnZ( int AL, tBS_Val *ADat ) { mcIvlCnZ( AL, ADat ); return AL; } // -- counts within nonzero ivls
static int icpl_IvlA( int AL, tBS_Val *ADat ) { // -- average within nonzero ivls
  int qL = 0, qA = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qN, qV = ADat[qi]; // simplify loop's condition
    if (qV>0) { if (qL<0) { qN = 1; goto lab1; }; qL++; qA += qV; continue; };
    if (qV<0) { if (qL>0) { qN =-1; goto lab1; }; qL--; qA += qV; continue; };
    if (!qL) continue; qN = 0; lab1: qL = abs(qL); BS_VFil( qL, ADat+qi-qL, qA/qL ); qL = qN; qA = qV;
  }; return AL; } // for i
static int icpl_IvlM( int AL, tBS_Val *ADat ) { // -- max magnitude within nonzero ivls
  int qL = 0, qA = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qN, qV = ADat[qi]; // simplify loop's condition
    if (qV>0) { if (qL<0) { qN = 1; goto lab1; }; qL++; if (qA<qV) qA = qV; continue; };
    if (qV<0) { if (qL>0) { qN =-1; goto lab1; }; qL--; if (qA>qV) qA = qV; continue; };
    if (!qL) continue; qN = 0; lab1: qL = abs(qL); BS_VFil( qL, ADat+qi-qL, qA ); qL = qN; qA = qV;
  }; return AL; } // for i
static int icpl_IvlIB( int AL, tBS_Val *ADat ) { // -- interval border's
  int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; // simplify loop's condition
    if (qV>0) { if (qL<0) { qV = 1; goto lab1; }; qL++; continue; };
    if (qV<0) { if (qL>0) { qV =-1; goto lab1; }; qL--; continue; };
    if (!qL) continue; lab1: qL = abs(qL); tBS_Val *qDat = ADat+qi-qL;
    int qV0 = qDat[0], qV1 = qDat[qL-1]; BS_VFil( qL, qDat, abs(qV0)>abs(qV1)?qV0:qV1 ); qL = qV;
  }; return AL; } // for i
static void ippr_CNormByPitch( tBS_SpawnThreadRec *AP ) {
  for (int qi=AP->rFst;qi<AP->rLst;qi++) BS_VMul( mcVi(AP->rP[0]), ((tBS_Val**)AP->rP[1])[qi], ((float*)AP->rP[2])[qi] ); };

void BS2_Met( tBS_Val **ADatas, int ATag, tBS_CPLProc AProcMet, const int *ADirs, int ANDD, int AW, int AH ) { // calculate a metric along lines
  if (!AW) { AW = xW; AH = xH; }; switch (ATag&255) {
    case cBS_MetC: AProcMet = icpl_IvlC; break;
    case cBS_MetCnZ: AProcMet = icpl_IvlCnZ; break;
    case cBS_MetA: AProcMet = icpl_IvlA; break;
    case cBS_MetM: AProcMet = icpl_IvlM; break;
    case cBS_MetIB: AProcMet = icpl_IvlIB; break;
  }; // switch Tag&255
  for (int qj1=0;qj1<ANDD;qj1++) { int qd = ADirs[qj1]; if ((qd!=1) && (qd!=5)) goto lab1; };
  for (int qj2=0;qj2<ANDD;qj2++) if (ADirs[qj2]==5) // faster cartesian accessor
    BS_CPLCCols( AProcMet, AW, AH, ADatas[qj2] ); else BS_CPLCRows( AProcMet, AW, AH, ADatas[qj2] );
  return; lab1:
  BS_CPLA( AProcMet, ANDD, ADirs, AW, AH, ADatas );
  if (ATag&cBS_MetByPitch) { // normalize against directions' pixel pitch
    float qMuls[8]; for (int qi=0;qi<ANDD;qi++) qMuls[qi] = cBS_Pitches[ADirs[qi]-1];
    BS_SpawnThreads( cBS_SpawnPool|2, 0, ANDD, ippr_CNormByPitch, mcIv(AW*AH), ADatas, qMuls );
  }; } // if Tag and Pitch

// 1x1 wrapper of metrics' line_shaders _Met(_MetXYZ)
template <int cN> static int __fastcall mccpb_MeWMn8( tBS_Val *ADat ) { mcSortBAs<tBS_Val,8,cN>( ADat ); return mcAve<tBS_Val,cN>( ADat ) ; }
static int __fastcall mccpb_MeWA8( tBS_Val *ADat ) { return mcAve<tBS_Val,8>( ADat ); }
void BS2_Me8( tBS_Val *ADat, int ATag, tBS_CPLProc AProcMet, tBS_CPBProc AProcMe8, int AW, int AH ) {
  if (!AW) { AW = xW; AH = xH; }; const int cN = AW*AH; // a metric's wrapper
  tBS_Val *qLimbo[8]; BS_DatasAlloc( 7, cN, qLimbo+1 ); qLimbo[0] = ADat;
  BS_DataDup( 7, cN, ADat, qLimbo+1 );  // calculate a metric along lines
  BS2_Met( qLimbo, ATag&0xFFFF, AProcMet, cBS_AllDirsT, 8, AW, AH );
  // for (int qj=0;qj<8;qj++) { char qS[88]; sprintf( qS, "qC8_%d.raw",qj ); BS_Write( qS, cN, qLimbo[qj] ); };
  #if 0 // !!!!
  switch (ATag&255) { // ad-hoc process after line_shader
    case cBS_MetC: for (int qj=0;qj<8;qj++) BS_Abs( cN, qLimbo[qj] ); break; // if counts were signed
  }; // swich Tag&255
  #endif
  switch (ATag&0xF0000) { // pixel_shader after line_shader
    case cBS_Me8Mn1: AProcMe8 = mccpb_Mn<8>; break;
    case cBS_Me8Mn2: AProcMe8 = mccpb_MeWMn8<2>; break; 
    case cBS_Me8Mn3: AProcMe8 = mccpb_MeWMn8<3>; break; 
    case cBS_Me8Mx1: AProcMe8 = mccpb_Mx<8>; break;
    case cBS_Me8MnA: AProcMe8 = mccpb_MnA<8>; break;
    case cBS_Me8MxA: AProcMe8 = mccpb_MxA<8>; break;
    case cBS_Me8A:   AProcMe8 = mccpb_MeWA8; break;
  }; // switch Tag
  BS_CPB( AProcMe8 , 8|cCPBNDD, cN, qLimbo );
  BS_DatasRelease( 7, qLimbo+1 ); }

//============ Experiments ===================
//============================================

