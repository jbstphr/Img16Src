// I16_Ut7a.h : collimator detection (CPU-only)
#ifndef I16_UT7A_H_INCLUDED
#define I16_UT7A_H_INCLUDED

#include "BS_GIni.h" // A wrapper around INI-files' (or registry's) API

// The "Msk" here assumes <slightly blurred> shape within defaul masks' range ~"0..511".
// The workflow:
//   - PreProcess replaces collimator in Frm by ~flat area and exports;
//   - Unsharp processes Frm;
//   - PostProcess fixes aliases on collimator's border.

void BS7a_CollimatorLoad( GIniFile *AF ); // load estimation parameters (called by I16_Init)
  // load initialization profile [CO]
int BS7a_CollimatorEnabled(); // >0:is enabled
  // "[CO] nEnabled = ..."

bool BS7a_CollimatorPreProcess( int AW, int AH, tBS_Val *AFrm, // collimator-estimation
  tBS_Val *AMsk, tBS_Val **ATs, const char *ADbgPfx ); 
  // inp:Frm
  // out:Frm with collimator replaced by ~flat median
  //     Msk "0:coll", "511:nocoll"

void BS7a_CollimatorPostProcess( int AW, int AH, tBS_Val *AFrm, // collimator-fixation
  tBS_Val *AMsk, tBS_Val **ATs, const char *ADbgPfx ); 

#endif // I16_UT7A_H_INCLUDED
