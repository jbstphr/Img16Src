// I16_Ut3.h : CUDA mid-level routines (certain dimension either as CUDA basic & CPU mid-level)
// the routines here use fixed frames BS3_Dim=>BS3_W,BS3_H

#ifndef I16_UT3_H_INCLUDED
#define I16_UT3_H_INCLUDED

#include "BS_Ut2.h"
//================= commons ==================
//============================================
void BS3_Dim( int AW, int AH, int ADbgDmp=0 ); // assigns|returns standard frame's dimensions
int BS3_W(); int BS3_H(); int BS3_DbgDmp();

void BS3_MimicLap( const tBS_Val *ADat, tBS_Val *ALap ); // mimic |Dat| and sign(Lap)
void BS3_MimicDat( const tBS_Val *ALap, tBS_Val *ADat ); // mimic sign(Lap) and |Dat|
void BS3_VMimicNZ( tBS_Val *ADat, int AV ); // set pix<>0 to V
void BS3_MovAve2D( int ARad, const tBS_Val *ASrc, tBS_Val *ADst, tBS_Val *AT ); // moving average
void BS3_MovAve2Dp( int ARad, tBS_Val *ADat, tBS_Val **ATs2, int Apyr ); // 2D average within downsample
void BS3_Hess( int ARad, tBS_Val **AHess ); // inp#0:Lap, out#0:Hessian #1:Lap_blr
bool BS3_Sum( float K, tBS_Val *A, float K1, const tBS_Val *A1, // summation without overflow
  float K2=0, const tBS_Val *A2=NULL, float K3=0, const tBS_Val *A3=NULL, float K4=0, const tBS_Val *A4=NULL );
bool BS3_MskSum( tBS_Val *Dst, tBS_Val *T, int Nrm, float K0, tBS_Val *A0, // sum normalized  masks
  float K1=0, const tBS_Val *A1=NULL, float K2=0, const tBS_Val *A2=NULL, float K3=0, const tBS_Val *A3=NULL );
bool BS3_AccMim( tBS_Val *Acc, tBS_Val *T, float K0, const tBS_Val *E0, // add mimiclap to accumulator
  const tBS_Val *L0, float K1, const tBS_Val *E1, const tBS_Val *L1 );
void BS3_RngLapByConductor( const tBS_Val *ACnd, tBS_Val *ALap ); // constraint lap by cardinal conductor
void BS3_RidgeByConductors( tBS_Val *ARdg, const tBS_Val *ACDev, // mimic a ridge by lap & conductors
  const tBS_Val *ACDer0, const tBS_Val *ACDer1 ); // (the conductors shall be normalized)
  // inp Rdg - laplacian; Cnd - positive conductor ~|Lap| o r|Deviation|, CDer - negative conductors (|Der3x3|)

void BS3_MxLapXLimbo( int ARad, float AK, int ATolGap, int AZ0, float AKAO,
  tBS_Val *ALap, tBS_Val **ATs, tBS_Val *ATmpH ); // MxLapStd with ~extended_limbo
  // "TolGap<0" - increment_mode, ">=0" - filter_mode
  // refer topic ~"_MxLapXLimbo" in _ReadMe_GenC12.txt

//========= Gamut's transformations ==========
//============================================
void BS3_Cut( const tBS2_Cut *AP, tBS_Val *ADat ); // cut & rescale
void BS3_HistCut( const tBS2_HistCut *AP, tBS_Val *ADat, const tBS_Val *ARef, tBS_Val **ATs2 ); // cut by histogram 

// gamma-correction
void BS3_Gamma( float AGamma, int AGam, tBS_Val *ADat, int AGamStep=0 ); // apply gamma-curve

//=========== Filters (conscious) ============
//============================================
void BS3_RadonMM( const tBS2_RadonMM *AP, const tBS_Val *ALapOrg, 
  tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp, tBS_Val *AENm,
  tBS_Val **ATs, tBS_Val **AQuadLapsH ); // core stats of rough Radon 
  // outs: Dev:Radon' deviation, Blr:blurlap, Srp:sharplap.
  // out optional: sum(edge-naive-mimic-laps) as addition to Srp.
  // out optional: quadlaps. The concept of quad-laps is "locally robust but globally uncertain sign",
  // if RadQ=0 then output 4 quad-laps (other 4 would be identical);
  //           else output 8 quad-laps (with anisotropic blur).

#endif // I16_UT3_H_INCLUDED
