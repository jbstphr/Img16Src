// I16_Ut3.cpp : CUDA mid-level routines (certain dimension either as CUDA basic & CPU mid-level)
// the routines here use fixed frames BS3_Dim=>BS3_W,BS3_H
// by S.B.
#include "stdafx.h"
#include "math.h"
#include "cuda_runtime.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut3.h"
#include "BS_Ut4.h"
#include "BS_Ut7.h"

//================= commons ==================
//============================================
static int ixW, ixH, ixDbgDmp; // standard frame-size for local operations
void BS3_Dim( int AW, int AH, int ADbgDmp ) { ixW = AW; ixH = AH; ixDbgDmp = ADbgDmp; BC_Dim( AW, AH ); } // -- assigns standard frame's dimension
int BS3_W() { return ixW; } int BS3_H() { return ixH; } int BS3_DbgDmp() { return ixDbgDmp; }

void BS3_MimicLap( const tBS_Val *ADat, tBS_Val *ALap ) { extern void KMimicLap(const tBS_Val*,tBS_Val*); KMimicLap( ADat, ALap ); } // -- mimic |Dat| and sign(Lap)
void BS3_MimicDat( const tBS_Val *ALap, tBS_Val *ADat ) { extern void KMimicDat(const tBS_Val*,tBS_Val*); KMimicDat( ALap, ADat ); } // -- mimic sign(Lap) and |Dat|
void BS3_VMimicNZ( tBS_Val *ADat, int AV ) { extern void KVMimicNZ(tBS_Val*,int); KVMimicNZ( ADat, AV); } // -- set pix<>0 to V

void BS3_MovAve2D( int ARad, const tBS_Val *ASrc, tBS_Val *ADst, tBS_Val *AT ) {
  const int cW = ixW, cH = ixH, cN = cW*cH; switch (ARad) {
    case 1: BC_Ave3x3( ASrc, ADst ); break; case 2: BC_Ave5x5( ASrc, ADst ); break;
    default: BC_Ave2D( ASrc, ADst, ARad, AT ); }; }

void BS3_MovAve2Dp( int ARad, tBS_Val *ADat, tBS_Val **ATs2, int Apyr ) { // -- 2D average within downsample
  if (ARad<=0) return; if (Apyr<1) { BS3_MovAve2D( ARad, ADat, ATs2[0], ATs2[1] ); BC_Cpy( ATs2[0], ADat ); return; };
  const int cW = ixW, cH = ixH; BC_ScaleDn( Apyr, cW, cH, ADat, ATs2[0] ); BC_Ave2D( ADat, ATs2[0], ARad, ATs2[1] );
  BC_Cpy( ATs2[0], ADat ); BC_ScaleUp( Apyr, cW>>Apyr, cH>>Apyr, ADat, ATs2[0], true ); } // _DimUpC reduces lines' aliases

void BS3_Hess( int ARad, tBS_Val **AHess ) { tBS_Val *qT = AHess[3]; // -- inp#0:Lap, out#0:Hess #1:Lap_blr
  BS3_MovAve2D( ARad, AHess[0], AHess[1], qT ); BC_Abs( AHess[0] ); BS3_MovAve2D( ARad, AHess[0], AHess[2], qT ); 
  BC_Cpy( AHess[1], AHess[0] ); BC_Abs( AHess[0] ); BC_SubS( AHess[2], AHess[0] ); }

bool BS3_Sum( float K, tBS_Val *A, float K1, const tBS_Val *A1,  // -- summation without overflow
  float K2, const tBS_Val *A2, float K3, const tBS_Val *A3, float K4, const tBS_Val *A4 ) {
  if ((fabs(K-1)<1e-3f) && (fabs(K1)<1e-3f) && (fabs(K2)<1e-3f) && (fabs(K3)<1e-3f) && (fabs(K4)<1e-3f)) return false;
  const int N = ixW*ixH; int *I = (int*)BC_AllocD( N*sizeof(int) ); cudaMemset( I, 0, N*sizeof(int) ); BC_AddMul( A, I, K );
  if (A1 && fabs(K1)>1e-3f) BC_AddMul( A1, I, K1 ); if (A2 && fabs(K2)>1e-3f) BC_AddMul( A2, I, K2 );
  if (A3 && fabs(K3)>1e-3f) BC_AddMul( A3, I, K3 ); if (A4 && fabs(K4)>1e-3f) BC_AddMul( A4, I, K4 );
  BC_RngStd( I, A ); BC_ReleaseD( I ); return true; }

bool BS3_MskSum( tBS_Val *Dst, tBS_Val *T, int Nrm, float K0, tBS_Val *A0, // -- sum normalized  masks
  float K1, const tBS_Val *A1, float K2, const tBS_Val *A2, float K3, const tBS_Val *A3 ) {
  const int N = ixW*ixH; float KZ = 0; BC_ZFil( N, Dst );
  if (A0 && fabs(K0)>1e-3f) { KZ += K0; BC_Cpy( A0, T ); BC_VMul( T, K0 ); BC_Add( T, Dst ); }
  if (A1 && fabs(K1)>1e-3f) { KZ += K1; BC_Cpy( A1, T ); BC_VMul( T, K1 ); BC_Add( T, Dst ); }
  if (A2 && fabs(K2)>1e-3f) { KZ += K2; BC_Cpy( A2, T ); BC_VMul( T, K2 ); BC_Add( T, Dst ); }
  if (A3 && fabs(K3)>1e-3f) { KZ += K3; BC_Cpy( A3, T ); BC_VMul( T, K3 ); BC_Add( T, Dst ); }
  if (KZ<1e-3f) { BC_VFil( N, Dst, Nrm ); return false; }; int qAdd = (int)((1-KZ)*Nrm+0.4999f);
  if (qAdd) BC_VAdd( Dst, qAdd ); return true; }

bool BS3_AccMim( tBS_Val *Acc, tBS_Val *T, float K0, const tBS_Val *E0, // -- add mimiclap to accumulator
  const tBS_Val *L0, float K1, const tBS_Val *E1, const tBS_Val *L1 ) { 
  const bool b0 = K0>1e-3f, b1 = K1>1e-3f; if (!(b0||b1)) return false;
  if (b0) { BC_Cpy( L0, T ); BS3_MimicLap( E0, T ); BC_VMul( T, K0 ); BC_Add( T, Acc ); };
  if (b1) { BC_Cpy( L1, T ); BS3_MimicLap( E1, T ); BC_VMul( T, K1 ); BC_Add( T, Acc ); }; return true; }

void BS3_RngLapByConductor( const tBS_Val *ACnd, tBS_Val *ALap ) { // -- constraint lap by cardinal conductor
  extern void KRngLapByConductor(const tBS_Val*,tBS_Val*); KRngLapByConductor( ACnd, ALap ); }

void BS3_RidgeByConductors( tBS_Val *ARdg, const tBS_Val *ACDev, // -- mimic a ridge by lap & conductors
  const tBS_Val *ACDer0, const tBS_Val *ACDer1 ) { // (the conductors shall be normalized)
  extern void KRidgeByConductors(tBS_Val*,const tBS_Val*,const tBS_Val*,const tBS_Val*);
  KRidgeByConductors( ARdg, ACDev, ACDer0, ACDer1 ); }

// MxLap with extended Limbo
#if 0 // reserved code for skew
  static void iMxLap16SkewCUDA( bool ADe, tBS_Val *APtrD, tBS_Val *ATmpH ) { const int cN = ixW*ixH; 
    extern void iMxLap16Skew( bool, tBS_Val*,int,int); extern void KMxLapRngZ(int,tBS_Val*);
    BC_Move( -cN*cBS_Sz, ATmpH, APtrD ); iMxLap16Skew( ADe, ATmpH, ixW, ixH ); BC_Move( cN*cBS_Sz, ATmpH, APtrD ); }
#endif
void BS3_MxLapXLimbo( int ARad, float AK, int ATolGap, int AZ0, float AKAO,
  tBS_Val *ALap, tBS_Val **ATs, tBS_Val *ATmpH ) { // -- MxLapStd with ~extended_limbo
  // called in BS_Ut5.cpp,_UnSharpCUDA{"RadMM,KMM,TGMM,ZMM,KAO"}; BS_Ut6.cpp,_DePoisson{"1,KU,0,0,0.0f"}
  extern void KMxLapXLimbo(tBS_Val*,float,int,tBS_Val**,int); extern void KMxLapAutOscMsk(tBS_Val**,int); 
  if (ARad<=0) return; BC_Cpy( ALap, ATs[0] ); int cNDD;
  // disable BS2A_Shear here (as use byproduct directional averages in Ts[])
  if (ATolGap<0) { // Lap+=K*MxLapStd
    cNDD = BC_MxLapStd( ARad, ATs[0], ATs+1, 1.0f, AZ0, 0, 0.0f );
    if (AK>1e-3f) { BC_VMul( ATs[0], AK ); BC_Add( ATs[0], ALap ); };
  } else { // Lap:=Rng(Lap,low_est-TolGap,high_est+TolGap), _est:=K*func(MxLapStd,MxLapExact,Ave)
    cNDD = BC_MxLapStd( ARad, ATs[0], ATs+1, 1.0f, AZ0, 0, 0.0f ); // BC_Write( "Qa.raw", 0, ALap ); 
    KMxLapXLimbo( ALap, AK, ATolGap, ATs, cNDD );
  };  
  // KA0<>0:AutoOscillator (<0:copies dbgtrc's frames to Ts^)
  bool qAODbgInf = false; if (AKAO<-1e-3f) { AKAO = -AKAO; qAODbgInf = true; }; if (AKAO>1e-3f) {
    KMxLapAutOscMsk( ATs, cNDD ); BC_VMul( ATs[0], AKAO );
    if (qAODbgInf) { BC_Cpy( ATs[0], ATs[1] ); BC_Cpy( ALap, ATs[2] ); };
    const int cTg0 = 128; const tBS2_Cut cCut = {{0,512,512-cTg0},{0,0}};
    BS3_Cut( &cCut, ATs[0] ); BC_VAdd( ATs[0], cTg0 ); BC_Mul( ATs[0], ALap, -9 );
  }; } // of AutoOscillator

  /*
// -- |Counts| within continuous pos|neg intervals
static void iCntsMxInterim( int AW, int AH, tBS_Val **ADat, tBS_Val **ATsH8 ) {
  BC_Move(-AW*AH*cBS_Sz, ATsH8[0], ADat[0] ); BS2_CntsMxInternal( AW, AH, ATsH8 );
  // BC_Move( AW*AH*cBS_Sz, ATsH8[0], ADat[1] ); BC_Ave5x5( ADat[1], ADat[0] ); }
  BC_Move( AW*AH*cBS_Sz, ATsH8[0], ADat[0] ); }
void BS3_CntsMx( const tBS_Val *ALap, tBS_Val **ATs, int Apyr, tBS_Val **ATsH8 ) { // inp:Lap, out:Sum of<~2 mx |counts|
  const int cW = ixW, cH = ixH; BC_Cpy( ALap, ATs[0] );
  if (Apyr<1) { iCntsMxInterim( cW, cH, ATs, ATsH8 ); return; };
  BC_ScaleDn( Apyr, cW, cH, ATs[0], ATs[1], true );
  iCntsMxInterim( cW>>Apyr, cH>>Apyr, ATs, ATsH8 );
  BC_VMul( ATs[0], (float)(1<<Apyr) ); 
  BC_ScaleUp( Apyr, cW>>Apyr, cH>>Apyr, ATs[0], ATs[1], true ); }
*/

//========= Gamut's transformations ==========
//============================================
// cut & rescale
void BS3_Cut( const tBS2_Cut *AP, tBS_Val *ADat ) { // -- cut & rescale CUDA
  int qMG[2], *qLat = BS2_initCutLat( AP, qMG ); if (!qLat) return; const int cSlotTag = 0;
  BC_XLAT_Init( cSlotTag, qLat, qMG[0], qMG[1], BS2_initCutTD(AP->f[cBS2_rcfTD0]), BS2_initCutTD(AP->f[cBS2_rcfTD1]) ); delete qLat;
  BC_XLAT_Apply( cSlotTag, ADat ); BC_XLAT_Release( cSlotTag ); }

// histogram truncation
extern void iiHistMM(int,int,const unsigned*,float,float,int*); // minmax_range by histogram
void BS3_HistCut( const tBS2_HistCut *AP, tBS_Val *ADat, const tBS_Val *ARef, tBS_Val **ATs2 ) { // -- cut by histogram
  xBS2_rhCatchMM[0] = xBS2_rhCatchMM[1] = 0; const int cW = ixW, cH = ixH,
    cMnI = AP->n[cBS2_rhcBlack], cMxI = AP->n[cBS2_rhcWhite], cGaI1 = abs(cMxI-cMnI),
    cMgn = AP->n[cBS2_rhcMgn], cGaO1 = AP->n[cBS2_rhcGaO]-1; // params
  //if (1) BC_Write( "$yaa.raw", cW*cH, ADat );
  if (cMnI<cMxI) { if (cMnI) BC_VAdd( ADat, -cMnI ); } else { // negative input
    if (cMxI) BC_VAdd( ADat, -cMxI ); BC_Neg( ADat ); BC_VAdd( ADat, cGaI1 ); BC_Rng( ADat, 0, cGaI1 ); };
  //if (1) BC_Write( "$yab.raw", cW*cH, ADat );
  const int cw = cW-cMgn*2, ch = cH-cMgn*2; unsigned *qHist = new unsigned[cGaI1+1]; if (cMgn) BC_Dim( cw, ch );
  BC_GamHist( ADat+cMgn*cW+cMgn, cW, qHist, cGaI1+1, ARef ); if (cMgn) BC_Dim( cW, cH ); // get histogram
  iiHistMM( cGaI1+1, cw*ch, qHist, AP->f[cBS2_rhfCut0], AP->f[cBS2_rhfCut1], xBS2_rhCatchMM ); // minmax_range by histogram
  const tBS2_Cut cC = {xBS2_rhCatchMM[0],xBS2_rhCatchMM[1],cGaO1+1,AP->f[cBS2_rhfTD0],AP->f[cBS2_rhfTD1]};
  BS3_Cut( &cC, ADat ); delete qHist; } // apply range

// gamma-correction
void BS3_Gamma( float AGamma, int AGam, tBS_Val *ADat, int AGamStep ) { // -- apply gamma-curve
  int *qCrv = new int[AGam]; BS2_GammaCrv( AGamma, AGam, (unsigned*)qCrv );
  BC_XLAT_Init( 0, qCrv, 0, AGam, -2, -2 ); delete qCrv;
  BC_XLAT_Apply( 0, ADat ); BC_XLAT_Release( 0 ); }

//=========== Filters (conscious) ============
//============================================
static void iBS3_MovAve8( const tBS_Val *ALapOrg, int ARad, int ARadO, tBS_Val **ADs ) { // moving average by 8 dirs
  BC_Ave8( ALapOrg, ADs, ARad ); // standard average by 8 dirs
  if (ARadO<=0) return; tBS_Val *qOh = ADs[9]; // subtract orthogonal averages
  BC_Av5( ADs[0], qOh, ARadO ); BC_SubA( qOh, ADs[0] ); BC_Av6( ADs[1], qOh, ARadO ); BC_SubA( qOh, ADs[1] );
  BC_Av7( ADs[2], qOh, ARadO ); BC_SubA( qOh, ADs[2] ); BC_Av8( ADs[3], qOh, ARadO ); BC_SubA( qOh, ADs[3] );
  BC_Av1( ADs[4], qOh, ARadO ); BC_SubA( qOh, ADs[4] ); BC_Av2( ADs[5], qOh, ARadO ); BC_SubA( qOh, ADs[5] );
  BC_Av3( ADs[6], qOh, ARadO ); BC_SubA( qOh, ADs[6] ); BC_Av4( ADs[7], qOh, ARadO ); BC_SubA( qOh, ADs[7] );}

#define mcQuadLapPlain( I, Proc1, Proc2 ) { BC_Cpy( ADs[I], qT ); BC_SubA( ADs[I+4], qT );\
  BC_Move( -cN*cBS_Sz, AQuadLapsH[I], qT ); }
#define mcQuadLapBlur2( I, Proc1, Proc2 ) { BC_Cpy( ADs[I], qT1 ); BC_SubA( ADs[I+4], qT1 ); BC_Cpy( qT1, qT2 );\
  Proc2( qT1, qT, ARadQO ); Proc1( qT, qT1, ARadQ ); BC_Move( -cN*cBS_Sz, AQuadLapsH[I], qT1 );\
  Proc1( qT2, qT, ARadQO ); Proc2( qT, qT2, ARadQ ); BC_Move( -cN*cBS_Sz, AQuadLapsH[4+I], qT2 ); }
static void iBS3_QuadLaps( tBS_Val **ADs, tBS_Val **AQuadLapsH, int ARadQ, int ARadQO ) { // quad-laps
  const int cN = ixW*ixH;
  if (ARadQ>0) { tBS_Val *qT = ADs[8], *qT1 = ADs[9], *qT2 = ADs[10];
    mcQuadLapBlur2( 0, BC_Av1, BC_Av5 ); mcQuadLapBlur2( 1, BC_Av2, BC_Av6 );
    mcQuadLapBlur2( 2, BC_Av3, BC_Av7 ); mcQuadLapBlur2( 3, BC_Av4, BC_Av8 );
  } else { tBS_Val *qT = ADs[8];
    mcQuadLapPlain( 0, BC_Av1, BC_Av5 ); mcQuadLapPlain( 1, BC_Av2, BC_Av6 );
    mcQuadLapPlain( 2, BC_Av3, BC_Av7 ); mcQuadLapPlain( 3, BC_Av4, BC_Av8 );
  }; }

static void iBS3_ENmLaps( tBS_Val **ALap8, tBS_Val *AENm, int ARad ) { // edge-naive-mimic-laps
  tBS_Val *qT = ALap8[8]; // const int cN = ixW*ixH; BC_Write( "q1.raw", cN, ALap8[0] ); BC_Write( "q2.raw", cN, qT );
  BC_ADr1( ALap8[0], qT, ARad ); BS3_MimicDat( ALap8[0], qT ); BC_Cpy( qT, AENm );
  BC_ADr2( ALap8[1], qT, ARad ); BS3_MimicDat( ALap8[1], qT ); BC_Add( qT, AENm );
  BC_ADr3( ALap8[2], qT, ARad ); BS3_MimicDat( ALap8[2], qT ); BC_Add( qT, AENm );
  BC_ADr4( ALap8[3], qT, ARad ); BS3_MimicDat( ALap8[3], qT ); BC_Add( qT, AENm );
  BC_ADr5( ALap8[4], qT, ARad ); BS3_MimicDat( ALap8[4], qT ); BC_Add( qT, AENm );
  BC_ADr6( ALap8[5], qT, ARad ); BS3_MimicDat( ALap8[5], qT ); BC_Add( qT, AENm );
  BC_ADr7( ALap8[6], qT, ARad ); BS3_MimicDat( ALap8[6], qT ); BC_Add( qT, AENm );
  BC_ADr8( ALap8[7], qT, ARad ); BS3_MimicDat( ALap8[7], qT ); BC_Add( qT, AENm ); }
  //  BC_Write( "q3.raw", cN, AENm );

extern void KRadonShannon(tBS_Val**,const tBS_Val*,tBS_Val*,tBS_Val*,tBS_Val*);
void BS3_RadonMM( const tBS2_RadonMM *AP, const tBS_Val *AOrg, 
  tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp, tBS_Val *AENm,
  tBS_Val **ATs, tBS_Val **AQuadLapsH ) { // -- core stats of rough Radon 
  int qRad = AP->n[cBS2_rrcRad]; if (qRad>16) {
    ErrMsg( "WRN BS3_RadonMM(R>16) exceeds macro in kernel.cu" ); qRad = 16; };
  iBS3_MovAve8( AOrg, qRad, AP->n[cBS2_rrcRadO], ATs );
  if (AQuadLapsH) // <quite optional> for consequent operations
    iBS3_QuadLaps( ATs, AQuadLapsH, AP->n[cBS2_rrcRadQ], AP->n[cBS2_rrcRadQO] );
  if (AENm) // <optional> addition to ~Srp
    iBS3_ENmLaps( ATs, AENm, AP->n[cBS2_rrcRadENm] );
  KRadonShannon( ATs, AOrg, ADev, ABlr, ASrp ); }

//============= Filters (ad-hoc) =============
//============================================

