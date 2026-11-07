// I16_Ut2a.h : CPU mid-level routines addition

#ifndef I16_UT2A_H_INCLUDED
#define I16_UT2A_H_INCLUDED

//=========== Skew/Shear ============
//===================================
// 'Shear>0':rotate 'CCW', '<0':rotate 'CW'

void BS2A_ShearCPU( int AFill, float AShear, tBS_Val *ADat, tBS_Val *ATmp,
  int AW0=0, int AH0=0, int AW1=0, int AW2=0 ); // shear by CPU
void BS2A_ShearCUDA( int AFill, float AShear, tBS_Val *ADat, tBS_Val *ATmp,
  int AW0=0, int AH0=0, int AW1=0, int AW2=0 ); // shear by CUDA
void BS2A_Shear( int AFill, float AShear, tBS_Val *ADat, tBS_Val *ATmpH, tBS_Val **ATmpD,
  int AW0=0, int AH0=0, int AW1=0, int AW2=0 ); // shear's wrapper

#endif // I16_UT2_H_INCLUDED
