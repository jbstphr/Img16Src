// I16_Ut7b.h : ribs-detection in CPU (within downsampled frame)
// no CUDA to avoid complication against frames's dimension in main prcess
#ifndef I16_UT7B_H_INCLUDED
#define I16_UT7B_H_INCLUDED

void BS7b_RibsDetectionLoad( GIniFile *AF ); // load detection parameters (called by I16_Init)

int BS7b_RibsDetectionEnabled(); // ">0: enabled some or another mode"
  // load initialization profile [RD]

void BS7b_RibsDetection( int AW, int AH, tBS_Val *AFrm,
  tBS_Val **ATs, const char *ADbgPfx ); // ribs-detection preprocess
  // inp:Frm
  // out:Frm with collimator replaced by ~flat median
  //     Msk "0:coll", "511:nocoll"


#endif // I16_UT7B_H_INCLUDED
