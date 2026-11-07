// I16_Ut2a.cpp : CPU mid-level routines addition

// by S.B.
#include "stdafx.h"
#include "math.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2.h"
#include "BS_Ut2a.h"

//=========== Skew/Shear ============
//===================================

//// CPU ////
// suffix '0' means 'inp'; suffix '1' means 'out'; Stride is the same for inp & out

static void iLnShear( int AFill, float AShear, int AS, int AL0, const tBS_Val *A0, int AL1, tBS_Val *A1 ) {
  float cShearFloor = (float)floor(AShear); AShear -= cShearFloor;
  int   cPixOfs = (AL0>>1)-(AL1>>1), cShearInt = (int)cShearFloor, c1 = (AL0-1)*AS;
  for (int qi=0;qi<AL1;qi++) { int qBuf4[4], q0 = qi+cPixOfs+cShearInt; // coordinate of dst[i] in src[0..]
    if (mcHerLikBuf4( AFill, AS, AL0, A0, q0, qBuf4 ))
      qBuf4[0] = (int)mcHerLikStdSpline( AFill, AShear, qBuf4 );
    *A1 = qBuf4[0]; A1 += AS; }; } // for i

typedef struct {float rShear; int rFill,rW0,rH0,rW1,rH1;} tShear; // W0=0:vertical, H0=0:horizontal
static void ippr_Shear( tBS_SpawnThreadRec *AP ) { const tShear *cSh = (const tShear*)AP->rP[0]; // spawn 'LnShear'
  tBS_Val *cA0 = (tBS_Val*)AP->rP[1], *qA1 = (tBS_Val*)AP->rP[2];
  float qShear = cSh->rShear; int cLst = AP->rLst, qH0 = cSh->rH0, qFill = cSh->rFill; 
  if (qH0) { int qW = cSh->rW1, qH1 = cSh->rH1; // H0<>0 means vertical
    for (int qi=AP->rFst;qi<cLst;qi++) iLnShear( qFill, qShear*(qi-(qW>>1)), qW, qH0, cA0+qi, qH1, qA1+qi );
  } else { int qW0 = cSh->rW0, qH = cSh->rH1, qW1 = cSh->rW1; // H0=0 means horizontal
    for (int qi=AP->rFst;qi<cLst;qi++) iLnShear( qFill, qShear*(qi-(qH>>1)), 1, qW0, cA0+qi*qW0, qW1, qA1+qi*qW1 );
  }; }
static void iShearVrt( int AFill, float AShear, int AW, int AH0, const tBS_Val *A0, int AH1, tBS_Val *A1 ) { // vertical
  const tShear cSh = { AShear,AFill,0,AH0,AW,AH1}; BS_SpawnThreads( 0, 0, AW, ippr_Shear, &cSh, A0, A1 ); }  // W0:=0
static void iShearHrz( int AFill, float AShear, int AH, int AW0, const tBS_Val *A0, int AW1, tBS_Val *A1 ) { // horizontal
  const tShear cSh = {-AShear,AFill,AW0,0,AW1,AH}; BS_SpawnThreads( 0, 0, AH, ippr_Shear, &cSh, A0, A1 ); }  // H0:=0

// standard procedure's enter to assign dimensions
#define mcShearEnter { if (fabs(AShear)<1e-3f) return;\
  if (!AW0) { AW0 = BS2_W(); AH0 = BS2_H(); }; if (!AW1) { AW1 = BS2_W(); AH1 = BS2_H(); }; }

void BS2A_ShearCPU( int AFill, float AShear, tBS_Val *ADat, tBS_Val *ATmp, // -- shear by CPU
  int AW0, int AH0, int AW1, int AH1 ) { mcShearEnter;
  if (AShear>0) { iShearHrz( AFill, AShear, AH0, AW0, ADat, AW1, ATmp ); iShearVrt( AFill, AShear, AW1, AH0, ATmp, AH1, ADat ); }
           else { iShearVrt( AFill, AShear, AW0, AH0, ADat, AH1, ATmp ); iShearHrz( AFill, AShear, AH1, AW0, ATmp, AW1, ADat ); }; }

//// CUDA ////
extern void KShear( int,float,int,int,const tBS_Val*,tBS_Val*); // Fill, Shear, W0, H0, A0, A1
static void iKShearVrt( int AFill, float AShear, int AW, int AH0, const tBS_Val *A0, int AH1, tBS_Val *A1 ) {
  BC_Dim( AW, AH1 ); KShear( AFill, AShear, 0, AH0, A0, A1 ); }
static void iKShearHrz( int AFill, float AShear, int AH, int AW0, const tBS_Val *A0, int AW1, tBS_Val *A1 ) { 
  BC_Dim( AW1, AH ); KShear( AFill, AShear, AW0, 0, A0, A1 ); }

void BS2A_ShearCUDA( int AFill, float AShear, tBS_Val *ADat, tBS_Val *ATmp, // -- shear by CUDA
  int AW0, int AH0, int AW1, int AH1 ) { mcShearEnter;
  if (AShear>0) { iKShearHrz( AFill, AShear, AH0, AW0, ADat, AW1, ATmp ); iKShearVrt( AFill, AShear, AW1, AH0, ATmp, AH1, ADat ); }
           else { iKShearVrt( AFill, AShear, AW0, AH0, ADat, AH1, ATmp ); iKShearHrz( AFill, AShear, AH1, AW0, ATmp, AW1, ADat ); }; }

//// Wrapper ////
void BS2A_Shear( int AFill, float AShear, tBS_Val *ADatH, tBS_Val *ATmpH, tBS_Val **ATmpD, // -- shear's wrapper
  int AW0, int AH0, int AW1, int AH1 ) { mcShearEnter;
  if (ATmpD) {
    BC_Move( AW0*AH0*cBS_Sz, ADatH, ATmpD[0] );
    BS2A_ShearCUDA( AFill, AShear, ATmpD[0], ATmpD[1], AW0, AH0, AW1, AH1 );
    BC_Move( -AW0*AH0*cBS_Sz, ADatH, ATmpD[0] );
  } else 
    BS2A_ShearCPU( AFill, AShear, ADatH, ATmpH, AW0, AH0, AW1, AH1 ); }
