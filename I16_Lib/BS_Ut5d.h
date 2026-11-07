// I16_Ut5d.h : extended mask-recount (CPU-only; CUDA=true means inp|out in GPU-memory)
// This module was imposed by sample AI-denoise. The term "recount" here refers to
// spatial counting of "empirically flat-areas" (which are more ~isotropic against
// "object-rich-areas"). The elementary process-flow:
//   stage a) trinarization by thresholds: "-1 in flat areas", "1 in edges",
// "0 otherwize".
//   stage b) TraceShader counts intervals of const sign by 8 directions;
// PixelShader determines the "metric of flatness" (e.g. ~max of signed counts).
//   stage c) Normalize & blur the "metric of flatness".

#ifndef I16_UT5D_H_INCLUDED
#define I16_UT5D_H_INCLUDED

#include "BS_GIni.h" // A wrapper around INI-files' (or registry's) API
int BS5d_MskEnabled(); // >0:is enabled
void BS5d_MskLoad( GIniFile *AF ); // load estimation parameters (called by I16_Init)

void BS5d_Msk( int AW, int AH, bool ACUDA, const char *ADbgPfx, // forks mask-recount thread
  const tBS_Val *AImg, const tBS_Val *AmskThk, tBS_Val *AmskGDi ); // these parameters may change by implementation
  // All implementations update mskGDi (other inp frames might be used or not).
  // Various implementations might receive and recount some or another masks.

void BS5d_MskWait(); // wait for forked mask-recount thread

#endif // I16_UT5D_H_INCLUDED
