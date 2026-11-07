// I16_Ut5c.h : quad-laps' treatment (CPU-only)
// The term "quadrature" here refers to "pi/4 from the radio" (not one of
// "K.F.W.Gauss's quadratures"). Initially there are 4 quad-laps from 8 basic
// dir-averages. Quad-laps give "locally robust but globally uncertain sign"
// (i.e. local zero_isoline of some 1 or 2 quad-lap<s> would be naturally
// "better" than cartesian lalacian's one. But that is not not so easy, whether
// such shape should be a convex or a dimple locally). The concept is:
//   stage a) Anisotropic blur of quad-laps (along relevant directions; rather
// in callers' codes to utilize GPU vs CPU when possible). Anisotropic blur 
// produces quad-laps.
//   stage b) Count a commonplace laplacian along quad-lap's intervals of const
// sign. This estimates the sign of quad-lap within specific intervals. As this
// code mostly "counts along lines", it hardly gets substantial surplus from GPU.
// Schematic work-flow:
//  - bufs[cQuadLaps], bufp := _QuadLapsAlloc(bufs)
//  - BS?_RadonMM(...,bufp)
//  - _QuadLaps(...,bufp)
//  - LaplacianAccumulator += _QuadLapsWait()

#ifndef I16_UT5C_H_INCLUDED
#define I16_UT5C_H_INCLUDED


#include "BS_GIni.h" // A wrapper around INI-files' (or registry's) API
int BS5c_QuadLapsEnabled(); // >0:is enabled
void BS5c_QuadLapsLoad( GIniFile *AF ); // load estimation parameters (called by I16_Init)

const int cBS5c_QuadLaps = 8; // caller reserves #pointers for quad-laps treatment
tBS_Val **BS5c_QuadLapsAlloc( int AW, int AH, tBS_Val **ABufs, // pointer to allocated buffers
  int ARO, int ARQ ); // implementation-specific params
  // Allocates up to cQuadLaps buffers (to use among BS2_RadonMM,BS3_RadonMM,BS5c_QuadLaps).
  // Caller fills these buffers with quad-laps and then invokes procedure "_QuadLaps".
  // RO - orthogonal|small radius in RadonMM ('RadO','j_RadMMO').
  // RQ - primary|big radius of quad-laps in RadonMM ('RadQ').
  
void BS5c_QuadLaps( int AW, int AH, tBS_Val **ALaps, bool ACUDA, // forks quad-laps thread
  const char *ADbgPfx, const tBS_Val *ALapOrg, const tBS_Val *AMsk ); 
  // inp Lap:LapOrg, Laps:four quad-laps (~triangle-wgt-blurred)
  // out Lap:the max of estimated correlations between the Lap & Laps#i

tBS_Val *BS5c_QuadLapsWait(); // wait for forked quad-laps thread

#endif // I16_UT5C_H_INCLUDED
