// I16_Ut2.h : CPU mid-level routines with standard dimension
// the routines here typically work with default frames BS2_Dim=>BS2_W,BS2_H
// optional params ~'N' or ~'W,H' pass non-default dimensions in some routines

#ifndef I16_UT2_H_INCLUDED
#define I16_UT2_H_INCLUDED

//================= commons ==================
//============================================
void BS2_Dim( int AW, int AH ); // assigns|returns standard frame's dimensions
int BS2_W(); int BS2_H();
void BS2_MimicLap( const tBS_Val *ADat, tBS_Val *ALap, int AN=0 ); // mimic |Dat| and sign(Lap)
void BS2_MimicDat( const tBS_Val *ALap, tBS_Val *ADat, int AN=0 ); // mimic sign(Lap) and |Dat|
void BS2_VMimicNZ( tBS_Val *ADat, int AV, int AN=0 ); // set pix<>0 to V
void BS2_MovAve2D( int ARad, const tBS_Val *ASrc, tBS_Val *ADst, int AW=0, int AH=0 ); // moving average
void BS2_MovAve2Dp( int ARad, tBS_Val *ADat, tBS_Val **ATs1, int Apyr, int AW=0, int AH=0 ); // 2D average within downsample
  // note Ts1 is "^^" for compatibility with BS3_MovAve2Dp
void BS2_Hess( int ARad, tBS_Val **AHess, int AW=0, int AH=0 ); // inp#0:Lap, out#0:Hessian #1:Lap_blr
bool BS2_Sum( float K, tBS_Val *A, float K1, const tBS_Val *A1, // summation without overflow
  float K2=0, const tBS_Val *A2=NULL, float K3=0, const tBS_Val *A3=NULL, float K4=0, const tBS_Val *A4=NULL );
bool BS2_MskSum( tBS_Val *Dst, tBS_Val *T, int Nrm, float K0, tBS_Val *A0, // sum normalized  masks
  float K1=0, const tBS_Val *A1=NULL, float K2=0, const tBS_Val *A2=NULL, float K3=0, const tBS_Val *A3=NULL );
bool BS2_AccMim( tBS_Val *Acc, tBS_Val *T, float K0, const tBS_Val *E0, // add mimiclap to accumulator
  const tBS_Val *L0, float K1, const tBS_Val *E1, const tBS_Val *L1 );
void BS2_RngLapByConductor( const tBS_Val *ACnd, tBS_Val *ALap ); // constraint a lap by cardinal conductor
void BS2_RidgeByConductors( tBS_Val *ARdg, const tBS_Val *ACDev, // mimic a ridge by lap & conductors
  const tBS_Val *ACDer0, const tBS_Val *ACDer1 ); // (the conductors shall be normalized)
  // inp Rdg - laplacian; CDev - positive conductor ~|Lap| or |Deviation|, CDer - negative conductors (|Der3x3|)

void BS2_MxLapXLimbo( int ARad, float AK, int ATolGap, int AZ0, float AKAO,
  tBS_Val *ALap, tBS_Val **ATs, int AW=0, int AH=0 ); // MxLapStd with ~extended_limbo
  // "TolGap<0" - increment_mode, ">=0" - filter_mode
  // refer topic ~"_MxLapXLimbo" in _ReadMe_GenC12.txt
  // wrn: set "Shear=0" if use byproduct directional averages in TS[]
  
//========= Gamut's transformations ==========
//============================================
typedef struct{int n[3];float f[2];}tBS2_Cut; const int // cut & rescale
  cBS2_rcuBlack = 0, cBS2_rcuWhite = 1, cBS2_rcuGaO = 2, // inp_black,inp_white, out is 0..GaO-1
  cBS2_rcfTD0 = 0, // TDx<1: CBS_XLat::Init(-BS_XLatTD*|TDx|), i.e. const tangent within 0..1 (0.1 works same as -0.1),
  cBS2_rcfTD1 = 1; //   >=1: CBS_XLat::Init(TDx), i.e. the tangent is derivative at N ending points
void BS2_Cut( const tBS2_Cut *AP, tBS_Val *ADat, int AN=0 );
int *BS2_initCutLat( const tBS2_Cut *AP, int *AMM ); // for CBS_XLat::Init(...)
int BS2_initCutTD( float ATD ); // "float (which encodes meaning in sign,frac)" => "TD0,TD1"

typedef struct{int n[4];float f[4];}tBS2_HistCut; const int // cut by histogram 
  cBS2_rhcBlack = 0, cBS2_rhcWhite = 1, // inp_black,inp_white
  cBS2_rhcGaO   = 2, cBS2_rhcMgn   = 3, // out is 0..GaO-1, frame's margins
  cBS2_rhfCut0  = 0, // Cutx<1: proportion to #frame's pixels,
  cBS2_rhfCut1  = 1, //    >=1: explict int #pixels;
  cBS2_rhfTD0   = 2, cBS2_rhfTD1   = 3; // to 'cfTDx'
void BS2_HistCut( const tBS2_HistCut *AP, tBS_Val *ADat, const tBS_Val *ARef, tBS_Val **ATs2, int AW=0, int AH=0 );
extern int xBS2_rhCatchMM[]; // mn,mx

// gamma-correction
void BS2_GammaCrv( float AGamma, int AGam, unsigned *ACrv, int AGamStep=0 ); // gamma-curve for specified gamut
void BS2_Gamma( float AGamma, int AGam, tBS_Val *ADat, int AGamStep=0, int AN=0 ); // apply gamma-curve

//=========== Filters (conscious) ============
//============================================
typedef struct{int n[5];}tBS2_RadonMM; const int // MnMx denoise
  cBS2_rrcRad   = 0, cBS2_rrcRadO  = 1, // basic radius, basic ortho-radius
  cBS2_rrcRadQ  = 2, cBS2_rrcRadQO = 3, // QuadLaps' radiuses
  cBS2_rrcRadENm= 4;                    // EdgeNaiveMimic_derivation's radius
void BS2_RadonMM( const tBS2_RadonMM *AP, const tBS_Val *ALapOrg, 
  tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp, tBS_Val *AENm,
  tBS_Val **ATs, tBS_Val **AQuadLaps ); // core stats of rough Radon 
  // outs: Dev:Radon' deviation, Blr:blurlap, Srp:sharplap.
  // out optional: sum(edge-naive-mimic-laps) as addition to Srp.
  // out optional: quadlaps. The concept of quad-laps is "locally robust but globally uncertain sign",
  // if RadQ=0 then output 4 quad-laps (other 4 would be identical);
  //           else output 8 quad-laps (with anisotropic blur).

//============= Filters (ad-hoc) =============
//============================================

//=============== Metrics ====================
//============================================
const int // metrics' line_shaders NxN
          cBS_MetC = 1, cBS_MetCnZ = 2, // counts within continuous ivls pos|neg, nonzeros
          cBS_MetA = 3, // average
          cBS_MetM = 4, // max magnitude
          cBS_MetIB = 5, // interval border's
          cBS_MetByPitch = 1<<15; // normalize by pitch
void BS2_Met( tBS_Val **ADatas, int ATag, // calculate a metric along lines
  tBS_CPLProc AProcMet=NULL, const int *ADirs=cBS_AllDirsT, int ANDD=8, int AW=0, int AH=0 );
  // NDD<-1,-2 applies faster cartesian array-accesor<s>; Dirs[] shall be #1,#5 in this case

const int // 1x1 wrapper of line_shaders _Met(_Met...)
          cBS_Me8Mn1 = 1<<16, cBS_Me8Mn2 = 2<<16, cBS_Me8Mn3 = 3<<16, // return average sort[0..]
          cBS_Me8Mx1 = 4<<16, // return max
          cBS_Me8MnA = 5<<16, // mn magnitude
          cBS_Me8MxA = 6<<16, // mx magnitude
          cBS_Me8A   = 7<<16; // average
void BS2_Me8( tBS_Val *ADat, int ATag,
  tBS_CPLProc AProcMet=NULL, tBS_CPBProc AProcMe8=NULL, int AW=0, int AH=0 ); // a metric's wrapper

//============ Experiments ===================
//============================================

#endif // I16_UT2_H_INCLUDED
