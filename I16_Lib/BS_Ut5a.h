// I16_Ut5a.h : air-estimation in CPU (within downsampled frame)
// no CUDA to avoid complication against frames's dimension in main prcess
#ifndef I16_UT5A_H_INCLUDED
#define I16_UT5A_H_INCLUDED

#include "BS_GIni.h" // A wrapper around INI-files' (or registry's) API
void BS5a_AirEstimationLoad( GIniFile *AF ); // load estimation parameters (called by I16_Init)
int BS5a_AirEstimationEnabled(); // ">0: enabled some or another mode"

void BS5a_AirEstimation( int AW, int AH, const tBS_Val *AFrm, const char *ADbgPfx, int AAdj ); // forks air-detection thread

typedef struct { // results of air-detection thread
  int      rAiA[2]; // AiA#0:abs air, #1:abs collimator|implant
  float    rAir[2]; // Air#0:relative -"-
  tBS_Val *rAiF[1]; // AiF#0:BkgPtr, where 'pixel=0' denotes certain background, '=255' denotes certain object
  int      rDbgInfFlg; // pass to UnSharpExts.rAirInfFlg
} tAirEstimationWait;
void BS5a_AirEstimationWait( tAirEstimationWait *AAEW ); // wait for forked air-detection thread
  // _AirEstimation('Frm') shall be unshanged until this point; it rets Ai...

#endif // I16_UT5A_H_INCLUDED
