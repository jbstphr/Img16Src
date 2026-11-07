// I16_Ut7.h : ad-hoc service CPU-transforms (for collimator, ribcage etc)
#ifndef I16_UT7_H_INCLUDED
#define I16_UT7_H_INCLUDED

//================= commons ==================
//============================================
void BS7_Dim( int AW, int AH, int ADbgDmp=0 ); // assigns|returns standard frame's dimensions

//======= alter dirlayer's projection ========
//============================================
typedef bool(*tBS7_ProjectionProc)(int,int,tBS_Val*); // DirIdx,PrjLen,PrjPixels; ret:"1 write back", "0:don't"
  //  DirIdx is zero-based index in passed Dirs (it is not 'one-based direction')
void BS7_Projection( tBS_CPLProc ACPLProcInp, tBS7_ProjectionProc AProjectionProc,
  int ANDD, const int *ADirs, tBS_Val **ADatas, const int *AWH=NULL ); // -- wrapper for '_ProjectionProc'
  // a] (if<>nil) invoke CPLProcInp (counts single 16-bit metric of whole line; always rets `1`)
  // b] extracts dirlayer's prj from pixel#0 of each line; invokes 'ProjectionProc'; 
  // c] if ProjectionProc rets 'true' then writes prj back 

//===== linear shear/warp by skew_period =====
//============================================
int inline BS7_SkewShear( int ASkew, int ADim ) { return ASkew?ADim/(2*ASkew):0; }; // 'mx_shear-1' for skew_period,dim
void BS7_SkewH( int ASkew, bool AFill, tBS_Val *AData, int AW=0, int AH=0 ); // hrz linear shear/warp
void BS7_SkewV( int ASkew, bool AFill, tBS_Val *AData, int AW=0, int AH=0 ); // vrt linear shear/warp
  // positive Skew (either V,H) rotates CCW; negative Skew rotates CW.
  // note: omit the Fill when next call is ~_Shrink and shear's gap is big
const int cBS7_SkewAllocFill = 1<<16, // fill in _Skew? (with lowest 16 bits)
          cBS7_SkewAllocMarg = 2<<16, // complement rectangular margins
          cBS7_SkewAllocBeve = 4<<16; // complement blurred bevels
tBS_Val *BS7_SkewAlloc( int ASkewH, int ASkewV, const tBS_Val *AData, int *AWH, int ATag=cBS7_SkewAllocFill ); // extend & skew image
  // note: it calls '_VFil( Tag&FFFFh)'; here WH[2] is only out (altered xW,xH from DSeg);
void BS7_SkewRelease( int ASkewH, int ASkewV, tBS_Val *AData, const int *AWH, tBS_Val *ADataSkew ); // release skewed image
float inline BS7_SkewTangent( int ASkew ) { return ASkew?1.0f/ASkew:0; } // change of y against 1-pix movement lft>rgt 
  // inp Skew could be either SkewH,SkewV; full tangent is '_SkewTangent(SkewH)+_SkewTangent(SkewV)'
  // note: as our Y axis goes top>btm this value is rather 'minus algebraic tangent'
void BS7_SkewFillMgn( bool AforSkewV, tBS_Val *AData, int AV, int AW=0, int AH=0 ); // forSkewV=false:fill lft&rgt cols, =true:top&btm rows
  // note: call when skew "the same frame" (without "_SkewAlloc") before "_Skew(Fill=true)"

//============ WaterShed-metrics =============
//============================================
// WaterShed-metric for a hull-shape '/^\' (also Posterior ribs projected)
inline void mcWSMP16( int AL, tBS_Val *ADat, int &Am0, int &Am1 ) { int qm0 = 0, qm1 = 0, qM = ADat[0]; for (int qi=1;qi<AL;qi++) { // Posterior WS_center
  int qV = ADat[qi]; if (qM<qV) { qM = qV; qm0 = qm1 = qi; } else if (qM==qV) qm1 = qi; }; BS_VFil( qm1-qm0+1, ADat+qm0, qM ); Am0 = qm0; Am1 = qm1; }
inline void mcWSMPd16( int AL, tBS_Val *ADat ) { int qm0, qm1; mcWSMP16( AL, ADat, qm0, qm1 ); // WS_descend (1st descend CeToPe)
  int qM0 = ADat[qm0]; while (--qm0>=0) { int qV = ADat[qm0]; if (qM0>qV) qM0 = qV; else ADat[qm0] = qM0; };
  int qM1 = ADat[qm1]; while (++qm1<AL) { int qV = ADat[qm1]; if (qM1>qV) qM1 = qV; else ADat[qm1] = qM1; }; }
inline void mcWSMPa16( int AL, tBS_Val *ADat ) { int qm0, qm1; mcWSMP16( AL, ADat, qm0, qm1 ); // WS_ascend (1st ascend PeToCe)
  int qM0 = ADat[0],    qi0 = 0;    while (++qi0<qm0) { int qV = ADat[qi0]; if (qM0<qV) qM0 = qV; else ADat[qi0] = qM0; };
  int qM1 = ADat[AL-1], qi1 = AL-1; while (--qi1>qm1) { int qV = ADat[qi1]; if (qM1<qV) qM1 = qV; else ADat[qi1] = qM1; }; }

// WaterShed-metric for a convex-shape '\_/' (also Anterior ribs projected)
inline void mcWSMA16( int AL, tBS_Val *ADat, int &Am0, int &Am1 ) { int qm0 = 0, qm1 = 0, qM = ADat[0]; for (int qi=1;qi<AL;qi++) { // Anterior WS_center
  int qV = ADat[qi]; if (qM>qV) { qM = qV; qm0 = qm1 = qi; } else if (qM==qV) qm1 = qi; }; BS_VFil( qm1-qm0+1, ADat+qm0, qM ); Am0 = qm0; Am1 = qm1; }
inline void mcWSMAd16( int AL, tBS_Val *ADat ) { int qm0, qm1; mcWSMA16( AL, ADat, qm0, qm1 ); // WSC_descend (1st descend PeToCe)
  int qM0 = ADat[0],    qi0 = 0;    while (++qi0<qm0) { int qV = ADat[qi0]; if (qM0>qV) qM0 = qV; else ADat[qi0] = qM0; };
  int qM1 = ADat[AL-1], qi1 = AL-1; while (--qi1>qm1) { int qV = ADat[qi1]; if (qM1>qV) qM1 = qV; else ADat[qi1] = qM1; }; }
//inline void mcWSMAa16( int AL, tBS_Val *ADat ) { int qm0, qm1; mcWSMA16( AL, ADat, qm0, qm1 ); // WSC_ascend (1st ascend CeToPe)
//  int qM0 = ADat[qm0]; while (--qm0>=0) { int qV = ADat[qm0]; if (qM0<qV) qM0 = qV; else ADat[qm0] = qM0; };
//  int qM1 = ADat[qm1]; while (++qm1<AL) { int qV = ADat[qm1]; if (qM1<qV) qM1 = qV; else ADat[qm1] = qM1; }; }

#endif // I16_UT7_H_INCLUDED
