// I16_Ut4.h : the routines with (~mostly) CPU|CUDA as parameter

#ifndef I16_UT4_H_INCLUDED
#define I16_UT4_H_INCLUDED
// Note1: the basic algorithms here use fixed frames; CPU:BS2_Dim=>BS2_W,BS2_H, CUDA:BS3_Dim=>BS3_W,BS3_H
// Note2: calls BS2_CLAHE_Tra... invoke internally BS2_HistCut and then define input's xBS2_rhCatchMM
const int cBS4_XLATGamLg = 9; // '=8':faster & more pixellation, '>8':slower & less pixellation

typedef struct { // beyond basic parameters of 'BS_GamCLAHE' & 'BC_GamCLAHE'
  float rSLCurvExtMn,rSLCurvExtMx; // for _Transform (contrast limitations)
} tBS4_CLAHE_Exts;
 
int BS4_CLAHE_Transform( int AB, int AGam, int ABlr, const tBS_Val *AFrmOrg,
  tBS_Val **ADs, bool ACUDA, const char *ADbgPfx, tBS4_CLAHE_Exts *AExts, const tBS_Val *ARef=NULL );
  // plain CLAHE (might produce aliases like ~median_effects) 
  // inp 'Ref<>nil':weights, '=nil':weigth=1 for each pixel
  //     'bhMap<>nil': histogram excludes pixel=0
  // out #0:CLAHE #1:~8bit_FrmOrg

void BS4_FrmNrm_CPU( int AGam, tBS_Val *AFrm, tBS_Val **ATs2, int AW=0, int AH=0 );
void BS4_FrmNrm_CUDA( int, tBS_Val*, tBS_Val** ); // normalized frame for misc CLAHE-procedures
  // following procedures need plain normalization of input FrmOrg (which shall have ~8-bit gamut)

int BS4_CLAHE_TraAveRef( int AB, int ABlr, const tBS_Val *AFrmNrm, // ref_averages (aliase-prone in areas of steady gradients)
  tBS_Val **ADs, bool ACUDA, tBS_Val *ARef0=(tBS_Val*)-1, tBS_Val *ARef1=NULL, tBS_Val *ARef2=NULL, tBS_Val *ARef3=NULL );
  // Ref='(tBS_Val*)-1' returns commonplace histogram <i.e. the counts>
  // inp FrmNrm shall be 8-bit; out: Ds#i ref_average of Ref#i

int BS4_CLAHE_TraRefAir( int AB, int ABlr, const tBS_Val *AFrmNrm, // 'TraAveRef' for voluntary dimension
  const tBS_Val *ARef, tBS_Val **ADs, int AW=0, int AH=0 );
  // CPU-only; used by BS5a_AirEstimation (which downscales frames)

#endif // I16_UT4_H_INCLUDED
