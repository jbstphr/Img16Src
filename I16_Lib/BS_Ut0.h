// I16_Ut0.h : CPU general-purpose utilities
#ifndef I16_UT0_H_INCLUDED
#define I16_UT0_H_INCLUDED

#include "BS_Cmn.h"

template<typename tN>inline void mcVFil( int AL, tN *ADat, tN AV ) { for (int qi=0;qi<AL;qi++) ADat[qi] = AV; }
template<typename tN,int AL>inline void mcVFil( tN *ADat, tN AV ) { for (int qi=0;qi<AL;qi++) ADat[qi] = AV; }
template<typename tN>inline int mcSum( int AL, const tN *ADat ) { int qA = 0; for (int qi=0;qi<AL;qi++) qA += ADat[qi]; return qA; }
template<typename tN,int AL>inline int mcSum( const tN *ADat ) { int qA = 0; for (int qi=0;qi<AL;qi++) qA += ADat[qi]; return qA; }
template<typename tN>inline int mcAve( int AL, const tN *ADat ) { return mcSum<tN>(AL,ADat)/AL; } 
template<typename tN,int AL>inline int mcAve( const tN *ADat ) { return mcSum<tN,AL>(ADat)/AL; }
template<typename tN>inline tLL mcSumLL( int AL, const tN *ADat ) { tLL qA = 0; for (int qi=0;qi<AL;qi++) qA += ADat[qi]; return qA; }
template<typename tN>inline float mcAveLL( int AL, const tN *ADat ) { return (float)mcSumLL<tN>(AL,ADat)/AL; }
#define mcSortAs(V,I) {int qv=V[I];if(qv<V[I-1]){int j=0; while(qv>V[j++]); memmove(V+j,V+j-1,(I-j+1)*sizeof(V[0])); V[j-1]=qv;};} // insertion sort within '0..I'
#define mcSortAsA(V,I) {int qv=V[I];if(abs(qv)<abs(V[I-1])){int j=0; while(abs(qv)>abs(V[j++])); memmove(V+j,V+j-1,(I-j+1)*sizeof(V[0])); V[j-1]=qv;};} // insertion sort within '0..I'
  // to sort Dat[0..N-1]: invoke "for (int i=1;i<N;i++) mcSortAs(Dat,i)"
template<typename tN, int cN, int cN2> inline void mcSortBAs( tN *B ) { for (int qi=0;qi<cN2;qi++) { // bubble sort (find min, put the min in sorted order)
  tN qM = B[qi]; int qm = qi; for (int qj=qi+1;qj<cN;qj++) if (qM>B[qj]) { qM = B[qj]; qm = qj; }; tN qT = B[qi]; B[qi] = qM; B[qm] = qT; }; } // for i 
template<typename tN, int cN, int cN2> inline void mcSortBDe( tN *B ) { for (int qi=0;qi<cN2;qi++) { // bubble sort (find max, put the max in sorted order)
  tN qM = B[qi]; int qm = qi; for (int qj=qi+1;qj<cN;qj++) if (qM<B[qj]) { qM = B[qj]; qm = qj; }; tN qT = B[qi]; B[qi] = qM; B[qm] = qT; }; } // for i 
  // 'cN2=cN-1':sorts whole buffer (including last element #N-1), '=median_position+1':sorts for median

// prefix "BS_" is acronym of "Bone-Suppression-project" in ~2018 (where the array-accessors below were inroduced)

//== Ticks' mark ==
//=================
class CBS_Ticks{ int rTicks0,rTicks; public: char rPfx[8]; CBS_Ticks(); void Reset(); void Mark( const char *AMsg ); };
  // marks ticks im msec; links against 'WriteLn' (declared in ~stdafx.h)
  // 'Mark(Msg=nil)':WriteLn ticks since 'Reset'; 'Mark(Msg<>nil)':WriteLn ticks since recent Mark
  
//== 2D-iterators-service (array-accessors) ==
//============================================
const int  cBS_LnDirPass[8] = {1,2,2,2,1,2,2,2}; // #passes (index for zero-bazed 'Dir-1')
typedef struct { // iterates through straight line in 2D-array
  int rW, rH, // array's dimensions
      rMaxLen, // aligned 'diagonal(WxH)+1' (a buffer reserves at least one traling-gap-pixel)
      // iterator's internal fields; might be not initialized when 'parameter LL<0'
      rDir, // the 'angular direction'
      rPass, // index of directional-pass
      rStX, rStY, // coordinates of starting_point
      rStXInc, rStYInc, // increments of coordinates of starting_points
      rXInc, rYInc, // typical increment between neighbor straight-line-pixels 
      rPeriod; // period of straight-line-segment
} tBS_LIRec;
  // the 'angular direction' is one-based tag '1..8' ('0' is not used); negative value (not
  // implemented here) could denote reverse order of straight-line-indexes
typedef struct { tBS_LIRec rLIR; int rDir,rLC; void *rDat[4]; int rInt[2]; } tBS_SpawnLIRec;
  // it consists of 'BS_LIRec' and of few of additional ad-hoc variables 
  //  'Dir': to be used by instead of 'LIR.Dir'
  //  'LC': lines' counter; it is not 'systematic' in current implenetation (for simplicity)
  //  'Dat','Int': ad-hoc variable-containers (used within 'SpawnLIProc' only)
typedef void(*tBS_SpawnLIProc)(int,int,tBS_SpawnLIRec*,int*,tBS_Val**); // process a line
  // 1st parameter 'LL':
  //   -1/-2: enter/leave specific thread, which processes few directions
  //   -3/-4: enter/leave process of particular 'angular direction'
  //   >0: line's length
  // 2nd parameter 'Diri': current index in the 'Data'
  // 3rd parameter 'SPR' : the spawn-record
  // 4th parameter 'BufI': line-pixels' indexes (offsets); if LL>0 then Buf[LL]=0
  // 5th parameter 'Data': pointer=containers of func-parameters (passed by caller to 'SpawnLIProc')
void BS_SpawnLI( int ANDD, const int *ADirs, tBS_SpawnLIProc AProc, int AW, int AH, tBS_Val **ADatas );
  // a wrapper of 'SpawnIterator' for 2D-array 'WxH'; bits in 'Dirs':
  //   0..3: the 'angular direction'
  //   4..7: explicit 'pass_index+1'; normally it is zero and the iterator automatically goes through
  //         one or more passes (as required by specific 'angular direction')

int BS_GetProjectionBufI( int ADir, int AW, int AH, int *ABufI );
  // reindex of pix#0 of each line (project a dirlayer so that ortho-dirlayer's iterator goes lft-to-rgt)
  // BufI shall reserve 'W+H' entries; ret actual #entries for specific 'Dir,W,H'

//====== areas' consolidation pipeline =======
//============================================
// suffixes A|B|C follow consequent fits of consolidation-pipeline for <~dirty-shaped,binary> areas:
// fit#A:limbo's ray_shader, #B:texture_shader produces sum_dirlayer, #C:ray_shader of a sum_dirlayer by rows,cols
typedef int(*tBS_CPLProc)(int,tBS_Val*); // 'LnLength', 'LnPixels'; rets #pixs to keep back
void BS_CPLA( tBS_CPLProc AProc, int ANDD, const int *ADirs, int AW, int AH, tBS_Val **ADatas ); // pipeline fit#A (along dirs)
  // note: NDD<3 and Dirs=1,5 applies _CPLCRows and/or _CPLCCols 
void BS_CPLCRows( tBS_CPLProc AProc, int AW, int AH, tBS_Val *AData, bool ANonBuf=false ); // pipeline fit#C (along rows)
void BS_CPLCCols( tBS_CPLProc AProc, int AW, int AH, tBS_Val *AData ); // pipeline fit#C (along cols)

const int cCPBNDD = 1<<5; typedef int(__fastcall *tBS_CPBProc)(tBS_Val*); // switch interpretation of result of _CPBProc 
template<int cN> int __fastcall mccpb_MxA(tBS_Val *ADat) { int qM = ADat[0]; for (int qi=1;qi<cN;qi++) { int qV = ADat[qi]; if (abs(qM)<abs(qV)) qM = qV; }; return qM; }
template<int cN> int __fastcall mccpb_MnA(tBS_Val *ADat) { int qM = ADat[0]; if (!qM) return 0; for (int qi=1;qi<cN;qi++) { int qV = ADat[qi]; if (!qV) return 0; if (abs(qM)>abs(qV)) qM = qV; }; return qM; }
template<int cN> int __fastcall mccpb_Mx(tBS_Val *ADat) { int qM = ADat[0]; for (int qi=1;qi<cN;qi++) { int qV = ADat[qi]; if (qM<qV) qM = qV; }; return qM; }
template<int cN> int __fastcall mccpb_Mn(tBS_Val *ADat) { int qM = ADat[0]; for (int qi=1;qi<cN;qi++) { int qV = ADat[qi]; if (qM>qV) qM = qV; }; return qM; }
template<int cN> int __fastcall mccpb_Sum(tBS_Val *ADat) { return mcSum<tBS_Val,cN>(ADat); }
template<int cN> int __fastcall mccpb_Ave(tBS_Val *ADat) { return mcAve<tBS_Val,cN>(ADat); }
template<int cN,int cSh> int __fastcall mccpb_IplPW(tBS_Val *ADat) { int qV = ADat[0], qi = qV>>cSh; if (qi>cN-2) return ADat[cN-1]; // piece-wise interpolation (bounded)
  if (qi<1) return !qi?((ADat[1]*qV)>>cSh):0; int qV0 = ADat[qi]; return qV0+((ADat[qi+1]-qV0)*(qV-(qi<<cSh))>>cSh); } // N>1 (equ to '#cal_points+1')
//template<int cN> int __fastcall mccpb_SortAs(tBS_Val *ADat) { for (int qi=1;qi<cN;qi++) mcSortAs(ADat,qi); return cN; }

void BS_CPB( tBS_CPBProc AProc, int ANDD, int AN, tBS_Val **ADatas, void *AOut=NULL ); // pipeline fit#B; starting Data's entry<ies> are dflt out.
  // NDD '#0..CPBNDD-1:#of inp entries'; 'shr 8:NOut (#entries in Out[?]) (it is typical manipulation with a lambda-expression);
  // CPLBNDD:_CPBProc rets single value of Pix[0], otherwise _CPBProc rets #Pixs to keep back & outs Pix[].
  // if 'Out<>nil' if 'NOut=0' then 'Out replaces Datas[0]' else 'Out[?] replace Datas[?]'.

typedef int(__fastcall *tBS_CPKProc)(const tBS_Val*); // single result for a 2D-filter's kernel; the same type in "mcBS_CPKThread.tBuf"
void BS_CPK( int AR, tBS_CPKProc AProc, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); // spawn kernel-procedure DxD, where D=R*2+1

typedef int(__fastcall *tBS_CPKProcInt)(const int*); // single result for a 2D-filter's kernel
void BS_CPKInt( int AR, tBS_CPKProcInt AProc, int AW, int AH, const int *ASrc, tBS_Val *ADst ); //

//================ 1D arrays =================
//============================================
template<typename tN>void inline BS_ZFil( int AL, tN *ADst ) { memset( ADst, 0, sizeof(tN)*AL ); } // Dat:=0
template<typename tN>void inline BS_VFil( int AL, tN *ADst, int AV ) { while (--AL>=0) *ADst++ = AV; } // Dat:=V
template<typename tN>void inline BS_VScat( int AL, tN *ADst, int AStride, int AV ) { while (--AL>=0) { *ADst = AV; ADst += AStride; }; } // Dat:=V
template<typename tN>void inline BS_Abs( int AL, tN *ADat ) { for (int qi=0;qi<AL;qi++) ADat[qi] = abs(ADat[qi]); } // Dat:=|Dat|
template<typename tN>void inline BS_Neg( int AL, tN *ADat ) { for (int qi=0;qi<AL;qi++) ADat[qi] = -ADat[qi]; } // Dat:=-Dat
template<typename tN>void inline BS_VAdd( int AL, tN *ADat, int AV ) { if (AV) for (int qi=0;qi<AL;qi++) ADat[qi] += AV; } // Dat+=V
void BS_VMul( int AL, tBS_Val *ADat, float AV ); // Dat *= V

void inline BS_Add( int AL, const tBS_Val *Aadd, tBS_Val *Aacc ) { for (int qi=0;qi<AL;qi++) Aacc[qi] += Aadd[qi]; } // Acc+=Add; don't check gamut
void inline BS_AddRg( int AL, const tBS_Val *AAdd, tBS_Val *AAcc ) { for (int qi=0;qi<AL;qi++) { int qV = AAdd[qi]+AAcc[qi]; AAcc[qi] = mcRg(qV); }; } // Acc+=Add
void BS_AddMul( int AL, const tBS_Val *AAdd, int *AAcc, float AV ); // Acc += V*Add
void inline BS_SubA( int AL, const tBS_Val *ASub, tBS_Val *AAcc ) { for (int qi=0;qi<AL;qi++) AAcc[qi] -= ASub[qi]; } // Acc:=Sub-Acc
void inline BS_SubS( int AL, const tBS_Val *ASub, tBS_Val *AAcc ) { for (int qi=0;qi<AL;qi++) AAcc[qi] = ASub[qi]-AAcc[qi]; } // Acc:=Acc-Sub
void inline BS_SubSelf( int AL, tBS_Val *ADat, int AGap ) { for (int qi=AL-1;qi>=AGap;qi--) ADat[qi] -= ADat[qi-AGap]; } // Dat[i] -= Dat[i-Gap]
bool BS_Rng( int AL, tBS_Val *ADat, int AMn, int AMx=0x7FFF ); // Dat:=range(Dat); Mn<=-7FFFh, Mx>=7FFFh one-side check
bool BS_RngA( int AL, tBS_Val *ALap, int AMx, int AMxP=0x7FFF ); // -Mx<V<Mx" or "-Mx<V<MxP"

void BS_Mul( int AL, const tBS_Val *AMul, tBS_Val *ADat, int ANrm ); // Mul:=Org*Mul/Nrm; Nrm<0 denotes "shr |Nrm|"
//void BS_Div( int AL, const tBS_Val *ADen, tBS_Val *ANom, int ANrm, int AMn=-0x7FFF, int AMx=0x7FFF ); // Nom:=|Nom|*Nrm/(|Den|+1); Nrm<0 denotes "shl |Nrm|"

void inline BS_Gath( int AL, const tBS_Val *ASrc, tBS_Val *ADst, const int *ABufI ) { for (int qi=0;qi<AL;qi++) ADst[qi] = ASrc[ABufI[qi]]; }; // gathering reindex
void inline BS_Scat( int AL, const tBS_Val *ASrc, tBS_Val *ADst, const int *ABufI ) { for (int qi=0;qi<AL;qi++) ADst[ABufI[qi]] = ASrc[qi]; }; // scattering reindex
template<typename tN> void inline mcFlip( int AL, tN *ADat ) { tN *qDat= ADat+AL-1; do { tN qT = *ADat; *ADat = *qDat; *qDat = qT; } while (++ADat<--qDat); } // flip 1D-array
template <typename tN> void inline mcFlipHrz( int AW, int AH, tN *ADat, int AStride ) { for (int qi=0;qi<AH;qi++) { mcFlip<tN>( AW, ADat ); ADat += AStride; }; }

template<typename tN> void inline mcGath( int AL, const tN *ASrc, tN *ADst, int AStride ) { while (--AL>=0) { *ADst++ = *ASrc; ASrc += AStride; }; }; // gathering copy
template<typename tN> void inline mcScat( int AL, const tN *ASrc, tN *ADst, int AStride ) { while (--AL>=0) { *ADst = *ASrc++; ADst += AStride; }; }; // scattering copy
template<typename tN> int inline mcHistNP( int ANP, const tN *AHis ) { const tN *qHis = AHis; if (ANP<0) { ANP = -ANP;
  while (ANP>=0) ANP -= *qHis--; if (ANP<0) qHis++; } else { while (ANP>=0) ANP -= *qHis++; if (ANP<=0) qHis--; }; return (int)(qHis-AHis); }
template<typename tN> inline void mcShiftLft( int AL, tN *ADat, int AShift ) { for (int qi=0;qi<AL;qi++) ADat[qi] <<= AShift; }
template<typename tN> inline void mcShiftRgt( int AL, tN *ADat, int AShift ) { for (int qi=0;qi<AL;qi++) ADat[qi] >>= AShift; }

inline int mcsqrti( int AI ) {union {float rf;int ri;} qv; qv.rf = (float)abs(AI); qv.ri -= 1<<23; qv.ri >>= 1; qv.ri += 1<<29; return AI<0?-(int)qv.rf:(int)qv.rf;} // int sqrt
inline int mcsqrtu( unsigned AU ) {union {float rf;int ri;} qv; qv.rf = (float)AU; qv.ri -= 1<<23; qv.ri >>= 1; qv.ri += 1<<29; return (int)qv.rf;} // int sqrt
inline int mcsign( unsigned AU, int ASign ) { if (!ASign) return 0; if (ASign<0) return -(int)AU; return (int)AU; }
template<int U> inline int mcsign( int ASign ) { if (!ASign) return 0; if (ASign<0) return -U; return U; }

//======== Gray-scale gamut's XLAT ==========
//===========================================
void BS_GamHist( int AW, int AH, const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam, const tBS_Val *ARef=NULL ); // Hist[Dat] '+=Ref' or '+=1'
const int cBS_XLatTDLg = 8, cBS_XLatTD = 1<<cBS_XLatTDLg; // norm_one for terminal_derivation (of derivation at the left-,rgt-ends)
class CBS_XLat{ public: const int *cX; int rM,rG,rV[4]; // S-like curve's TD>0:~XLat[0]-XLat[TD], <=0:~_XLatTD/norm_one
  void Init( const int *AXLat, int AMn, int AGam, int ATD0=0, int ATD1=0 ); void Apply( int AN, tBS_Val *ADat ); };

void BS_GamXLATFit( int AGam, int *AXLAT ); // scale so that XLAT[Gam-1]:=Gam-1
unsigned *BS_GamCLAHEHists( int AW, int AH, const tBS_Val *ADat, int AB, int AGam, const tBS_Val *ARef=NULL ); // get Hists[Ny][Nx]
  // B<=0 histograms whole frame with margin B
void BS_GamCLAHEApply( int AW, int AH, tBS_Val *ADat, int AB, int AGam, int *AHists ); // apply Hists[Ny][Nx]
  // B<=0 applies to the whole frame

//============= moving averages =============
//===========================================
#define mcMovAveCore(T,L,R,Inp,Out) { int qSum = Inp[0]*(R+2); /*non-weighed 1D-moving average with denominator 'R*2+1'*/\
  for (int qi1=1;qi1<R;qi1++) qSum += Inp[qi1]; const T cDen = (R<<1)+1, *cInp = Inp;\
  for (int qi2=R+2;  --qi2;) {qSum-=Inp[0];    qSum+=cInp[R]; cInp++;*Out++=qSum/cDen;};\
  for (int qi3=L-R*2;--qi3;) {qSum-=cInp[-R-1];qSum+=cInp[R]; cInp++;*Out++=qSum/cDen;};\
  for (int qi4=R+1;  --qi4;) {qSum-=cInp[-R-1];qSum+=Inp[L-1];cInp++;*Out++=qSum/cDen;}; }
template <typename tN, int cR> void mcMovAve( int AL, const tN *AInp, tN *AOut ) { mcMovAveCore( tN, AL, cR, AInp, AOut ); } // const rad, no centinels /use BS_movave/
template <typename tN> void inline mcMovAve1( int AL, const tN *AInp, tN *AOut ) { // radius=1
  AOut[0] = (AInp[0]*2+AInp[1])/3; AL--; for (int qi=1;qi<AL;qi++) AOut[qi] = (AInp[qi-1]+AInp[qi]+AInp[qi+1])/3; AOut[AL] = (AInp[AL-1]+AInp[AL])/3; }
template <typename tN> void mcMovAve( int AL, int AR, const tN *AInp, tN *AOut ) { if (AR<1) { memcpy( AOut, AInp, AL*sizeof(tN) ); return; }; // var radius, centinels
  if (AR>=AL/2) { int qA = mcAve<tN>(AL,AInp); for (int qi1=0;qi1<AL;qi1++) AOut[qi1] = qA; return; };
  if (AR==1) { mcMovAve1( AL, AInp, AOut ); return ; }; 
  mcMovAveCore( tN, AL, AR, AInp, AOut ); }

const float cBS_Pitches[8]= {1,cPitch22,cPitch45,cPitch22,1,cPitch22,cPitch45,cPitch22};
inline int BS_SDist( int ADir, int ADist ) { int qDist = (int)(abs(ADist)/cBS_Pitches[ADir-1]+0.49999f);
  if (!qDist) qDist = 1; return ADist<0?-qDist:qDist; } // adhoc macro to interpret "_MovAve( R<0, ... )"
const int   cBS_AllDirs[8]  = {5,6,7,8,1,2,3,4}, // std limb of "conditionally radial" dirs
            cBS_AllDirsT[8] = {1,2,3,4,5,6,7,8}; // std limb of "conditionally transversal" dirs
void BS_MovAve1D( int AL, int AR, const tBS_Val *AInp, tBS_Val *AOut ); // ~mcMovAve<BS_Val,...> with sentinels
void BS_MovAve( int AR, int AW, int AH, tBS_Val **ADatas, const int *ADirs=cBS_AllDirsT, int ANDD=8 ); // threads of 2D-movave 
void BS_MovAveX( int AR, int AW, int AH, tBS_Val *AData ); void BS_MovAveY(int,int,int,tBS_Val*); // horz,vert average
void BS_MovAve2D( int AR, int AW, int AH, tBS_Val *AData, unsigned AFlg=7 ); // engineering 2D-movave (assumed to be isotropic)

//======= application spatial filters =======
//===========================================
void BS_Ave3x3( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); void BS_Ave5x5(int,int,const tBS_Val*,tBS_Val*); 
void BS_Med( int AR, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); // the medians below might be not exact ones
void BS_Dr1( int AW, int AH, const tBS_Val *ADat, tBS_Val *ADrX, tBS_Val *ADrY, bool Asgn ); // XYders (R=1); sgn=false:aboslute,=true:signed
void BS_DrA3x3( int AW, int AH, tBS_Val *ADat ); void BS_DrA3x3(int,int,tBS_Val*,tBS_Val**); // sum(abs(XYders))

void BS_MedB( int AR, int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); // 

//======= Radon-Shannon-denoise =======
//=====================================

void BS_Noise( int AR, int AW, int AH, // basic RadonShannon_noise
  const tBS_Val *AOrg, tBS_Val *ANoise, tBS_Val **ATs8, int AThrLg );
  // inp    Org: ~laplacian (around zero)
  //      ThrLg: regul_add*100 or threshold_log2
  // out  Noise: estimated noise 
  //        Ts8: directional averages (for adhoc usage)
void BS_MxLapRngZ( int AZ, int AN, tBS_Val *ALap ); // "if |Lap[i]|<Z then Lap[i]>>=2"
int BS_MxLapStd( // reduce shot_noise after RadonShannon_noise; returns #NDD of dir averages in TS^
  int AR, int AW, int AH, tBS_Val *ALap, tBS_Val **ATs, float AK, int AZ0=0, int AZ1=-1, float AShear=-0.033f );
  // call-sequence: "_MxLapRngZ(Z0); MaxLap; _MxLapRngZ(Z1)";
  // always: "Ave<MxLapStd<MxLapExact" or "MxLapExact<MxLapStd<Ave";
  // default Shear<>0 <might> slightly reduces aliases (it should be Ok for the "laplacians").
  
inline int BS_MxLapZ1dflt( int AZ0, int AZ1 ) { return AZ1<0?AZ0/3:AZ1; }

template<typename tN> void mcMedian1( int AL, tN *ADat ) { 
  tN qV0, qV1 = ADat[0], qV2 = ADat[1]; for (int qi=2;qi<AL;qi++) { qV0 = qV1; qV1 = qV2; qV2 = ADat[qi]; // advance buffer
    tN qV, qM; if (qV0>qV1) { qV = qV0; qM = qV1; } else { qM = qV0; qV = qV1; };
    if (qV>qV2) qV = __max( qM, qV2 ); ADat[qi-1] = qV; }; // for i
  ADat[0] = ADat[1]; ADat[AL-1] = ADat[AL-2]; }

//=========== down-/up-sample 2x2 ============
//============================================
void BS_DimDn( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst, bool AM=false ); // dimension down 2W*2H => W*H
void BS_DimUp( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); // dimension up W*H => 2W*2H
void BS_DimUpC( int AW, int AH, const tBS_Val *ASrc, tBS_Val *ADst ); // correction after 'dimension up' (i.e. 'unbinning2x2')
// pyramid's dn- & up-scale; M:true - max amplitude (default average); C:false - no _DimUpC (default _DimUpC)
tBS_Val *BS_ScaleDnA( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM=false ); void BS_ScaleDn(int,int,int,tBS_Val*,tBS_Val*,bool=false);
tBS_Val *BS_ScaleUpA( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC=true ); void BS_ScaleUp(int,int,int,tBS_Val*,tBS_Val*,bool=true);
 
//============== for margins ================
//===========================================
void BS_Bevels( int AW, int AH, tBS_Val *ADat, int AM0x, int AM1x, int AM0y=-1, int AM1y=-1 ); // replace 'M0' periph lines by inner; blur 'M1' periph lines; 'M0'<='M1' 
inline void BS_VBevelLftRgt( int AM, int AW, int AH, tBS_Val *ADat, int AV ) { // fill M lft,rgt margin-lines by V
  BS_VFil( AM, ADat, AV ); BS_VFil( AM, ADat+AW*AH-AM, AV ); for (int qi=1;qi<AH-1;qi++) BS_VFil( AM*2, ADat+qi*AW-AM, AV ); }
inline void BS_VBevelTopBtm( int AM, int AW, int AH, tBS_Val *ADat, int AV ) { // fill M top,btm margin-lines by V
  const int cL0 = AW*AM; BS_VFil( cL0, ADat, AV ); BS_VFil( cL0, ADat+AW*AH-cL0, AV ); }
template<int cM> inline void BS_VBevel( int AW, int AH, tBS_Val *ADat, int AV ) { // fill M rect margin-lines by V
  int cL0 = AW*cM+cM; BS_VFil( cL0, ADat, AV ); BS_VFil( cL0, ADat+AW*AH-cL0, AV ); for (int qi=cM+1;qi<AH-cM-1;qi++) BS_VFil( cM*2, ADat+qi*AW-cM, AV ); }
void BS_BevCorners( int AW, int AH, tBS_Val *AFrame, int ABevCorners, int ABevel, int AV );
void BS_DrawRect( int AW, int AH, tBS_Val *AData, const int ARect[4], int AMark ); // draw rectangle [lft,top,rgt,btm]

//====== Spatial intervals' counting =======
//==========================================
inline void mcIvlC( int AL, tBS_Val *ADat ) { // signed counts
  int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; // simplify loop's condition
  if (qV>0) { if (qL<0) { qV = qL; qL = 1; goto lab1; }; qL++; continue; }; if (qV<0) { if (qL>0) { qV = qL; qL =-1; goto lab1; }; qL--; continue; };
  if (!qL) continue; qV = qL; qL = 0; lab1: int qLa = abs(qV); // if (qLa&1) ADat[qi-qLa--] = 0; // optional topology's quirck to coincide with the moments
  BS_VFil( qLa, ADat+qi-qLa, qV ); }; }
inline void mcIvlCnZ( int AL, tBS_Val *ADat ) { // non-zero counts
  int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; // simplify loop's condition 
  if (qV) { qL++; continue; }; if (!qL) continue; BS_VFil( qL, ADat+qi-qL, qL ); qL = 0; }; }
inline void mcIvlCnMom( int AL, tBS_Val *ADat ) { // moments of non-zero intervals
  int qL = 0; ADat[AL] = 0; for (int qi=0;qi<=AL;qi++) { int qV = ADat[qi]; // simplify loop's condition 
  if (qV) { qL++; continue; }; if (qL<6) { BS_ZFil( qL, ADat+qi-qL ); qL = 0; continue; };
  tBS_Val *qDatL = ADat+qi-qL, *qDatR = ADat+qi-1; if (qL&1) *qDatL++ = 0; qL >>= 1;
  while (qL) { *qDatL++ = -qL; *qDatR-- = qL; qL--; }; }; }

template<typename tN> void mcIplIvl( int AL, tN *ADat, int AV0, int AV1 ) { // interpolate interval between V0+..V1-
  if (AL<2) { if (AL==1) ADat[0] = (AV0+AV1)/2; return; }; AV1 -= AV0; if (AV1) { AV1 = (AV1<<15)/(AL+1); 
  for (int qi=0;qi<AL;qi++) ADat[qi] = AV0+(((qi+1)*AV1)>>15); } else for (int qi=0;qi<AL;qi++) ADat[qi] = AV0; }
template<typename tN> void mcIplIvl( tN *ADat, int AI1 ) { mcIplIvl<tN>( AI1-1, ADat+1, ADat[0], ADat[AI1] ); } // interpolate Dat[1..I1-1] as Dat[0]..Dat[I1]
template<typename tN,tN AIvl> void mcIplIvls( int AL, tN *ADat, int AMxIvl=0x7FFFFFFF ) { // interpolate intervals of 'Ivl' shorter than 'MxIvl'
  ADat[AL] = (tN)~AIvl; int qC = 0; // simplify loop's condition (Dat shall locate in buffer 'L+1')
  for (int qi=0;qi<=AL;qi++) { tN qV = ADat[qi]; if (qV==AIvl) qC++; else if (qC) { if (qC==AL) return;
    int qi0 = qi-qC; tN qV0; if (!qi0) { qV0 = qV; } else { qV0 = ADat[qi0-1]; if (qi==AL) qV = qV0; };
    if (qC<AMxIvl) mcIplIvl( qC, ADat+qi0, qV0, qV ); qC = 0; }; }; } // if C<MxIvl | if C>0 | for i
template<typename tN,tN AIvl> void mcIplIvls0( int AL, tN *ADat, int AMxIvl=0x7FFFFFFF ) { // interpolate intervals of 'Ivl'
  int qC = 0; // it works in_place; Dat locates in buffer 'L'
  for (int qi=0;qi<=AL;qi++) { tN qV = qi<AL?ADat[qi]:(tN)~AIvl; if (qV==AIvl) qC++; else if (qC) { if (qC==AL) return;
    int qi0 = qi-qC; tN qV0; if (!qi0) { qV0 = qV; } else { qV0 = ADat[qi0-1]; if (qi==AL) qV = qV0; };
    if (qC<AMxIvl) mcIplIvl( qC, ADat+qi0, qV0, qV ); qC = 0; }; }; } // if C<MxIvl | if C>0 | for i

#endif // I16_UT0_H_INCLUDED