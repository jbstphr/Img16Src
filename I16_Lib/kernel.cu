// kernel.cu : impplements CUDA kernels
// by S.B.
//#include "memory.h" #include "cuda_runtime.h"
// NoteA: minimize includes here because of slow parsing of "*.cu" during "AddCudaCompilePropsDeps".
// This quirk likely relates to extremely slow check of files' existence & attributes (and gets yet
// worse under specific 3rd-party security-programs in Windows).
// NoteB: if Nsight's gutter-break-points are notoriously wrong then rely (somewhat) on "line nums"
// in the tool-windows "Debug|Windows|Breakpoints".
#include "stdafv.h"
#include "cuda_runtime.h"

///typedef short tBS_Val; // typedef unsigned short tBS_Card; // working data-type 
//__device__ inline int mcRng( int V, int Mn, int Mx ) { return V<Mn?Mn:(V>Mx?Mx:V); }
//__device__ inline int mcRg( int V ) { return mcRng( V, -0x7FFF, 0x7FFF ); } // limit of working gamut
__device__ inline int mcsign( unsigned AU, int ASign ) { if (!ASign) return 0; if (ASign<0) return -(int)AU; return (int)AU; }

//// services ////
// #define mcCache( AProc ) cudaFuncSetCacheConfig( AProc, cudaFuncCachePreferL1 )
// enter a warp with 2D indexes (here is no dedicated macros for rare 1D-operations)
#define mcEnter0() const int cx = blockDim.x*blockIdx.x+threadIdx.x; if (cx>=AW) return; const int cy = blockDim.y*blockIdx.y+threadIdx.y; if (cy>=AH) return; const int co = cy*AW+cx
#define mcEnter2(R) mcEnter0(); AInp += co; AOut += co; if ((cx<R)||(cx>=AW-R)||(cy<R)||(cy>=AH-R)) { AOut[0] = AInp[0]; return; }
//#define mcIpl(AX,AIX,AY0,AY1) ((AY0)+((AY1)-(AY0))*(AX)/(AIX)) // linear interpolation
//#define mcSwap( tN, A, B ) { tN T=(tN)(A); A=(tN)(B); B=T; }
template <typename tN> __device__ inline int mcAve( int AL, const tN *ADat ) { int qA = 0; for (int qi=0;qi<AL;qi++) qA += ADat[qi]; return qA/AL; } // simple average
template <typename tN, int AL> __device__ inline int mcSum( const tN *ADat ) { int qA = 0; for (int qi=0;qi<AL;qi++) qA += ADat[qi]; return qA; } // sum
template <typename tN, int AL> __device__ inline int mcAve( const tN *ADat ) { return mcSum<tN,AL>(ADat)/AL; } // simple average
extern int xw, xh, xBC_PMxBlkThreads; extern dim3 xG, xB; // managed in BS_Ut01.cpp

__global__ void kRowsCpy( int AW, int AH, int AWInp, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter0(); AOut[co] = AInp[cx+cy*AWInp]; }
void KRowsCpy( int AWInp, const tBS_Val *AInp, tBS_Val *AOut ) { kRowsCpy<<<xG,xB>>>( xw, xh, AWInp, AInp, AOut ); }
__global__ void kRowsCpy( int AW, int AH, const tBS_Val *AInp, int AWOut, tBS_Val *AOut ) { mcEnter0(); AOut[cx+cy*AWOut] = AInp[co]; }
void KRowsCpy( const tBS_Val *AInp, int AWOut, tBS_Val *AOut ) { kRowsCpy<<<xG,xB>>>( xw, xh, AInp, AWOut, AOut ); }

//=========== down-/up-sample 2x2 ============
//============================================
__global__ void kDimDn( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter0(); const int c2 = cy*AW*4+cx*2; AOut[co] = (AInp[c2]+AInp[c2+1]+AInp[c2+AW*2]+AInp[c2+AW*2+1])>>2; }
__global__ void kDimDnM( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter0(); const int c2 = cy*AW*4+cx*2; int qM = AInp[c2], qV = AInp[c2+1];
  if (abs(qM)<abs(qV)) qM = qV; qV = AInp[c2+AW*2]; if (abs(qM)<abs(qV)) qM = qV; qV = AInp[c2+AW*2+1]; if (abs(qM)<abs(qV)) qM = qV; AOut[co] = qM; }
void KDimDn( const tBS_Val *AInp, tBS_Val *AOut, bool AM ) { if (AM) kDimDnM<<<xG,xB>>>( xw, xh, AInp, AOut ); else kDimDn<<<xG,xB>>>( xw, xh, AInp, AOut ); } // 2W*2H => W*H
__global__ void kDimUp( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter0(); const int c2 = cy*AW*4+cx*2; AOut[c2] = AOut[c2+1] = AOut[c2+AW*2] = AOut[c2+AW*2+1] = AInp[co]; }
void KDimUp( const tBS_Val *AInp, tBS_Val *AOut ) { kDimUp<<<xG,xB>>>( xw, xh, AInp, AOut ); } // W*H => 2W*2H
__global__ void kDimUpC( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter2(1); int qo = ((cx&1)<<1)-1, qV = AInp[qo]-AInp[-qo]; qo = (((cy&1)<<1)-1)*AW; qV += AInp[qo]-AInp[-qo]; AOut[0] = (qV>>2)+AInp[0]; }
void KDimUpC( const tBS_Val *AInp, tBS_Val *AOut ) { kDimUpC<<<xG,xB>>>( xw, xh, AInp, AOut ); } // W*H => W*H

//================ 1D arrays =================
//============================================
__global__ void kAbs( int AW, int AH, tBS_Val *ADat ) { mcEnter0(); ADat[co] = abs(ADat[co]); } void KAbs( tBS_Val *ADat ) { kAbs<<<xG,xB>>>( xw, xh, ADat ); }
__global__ void kNeg( int AW, int AH, tBS_Val *ADat ) { mcEnter0(); ADat[co] = -ADat[co]; } void KNeg( tBS_Val *ADat ) { kNeg<<<xG,xB>>>( xw, xh, ADat ); }
__global__ void kVAdd( int AW, int AH, tBS_Val *ADat, int AV ) { mcEnter0(); ADat[co] += AV; } void KVAdd( tBS_Val *ADat, int AV ) { if (AV) kVAdd<<<xG,xB>>>( xw, xh, ADat, AV ); }
__global__ void kVMul( int AW, int AH, tBS_Val *ADat, float AK ) { mcEnter0(); int qv = (int)(ADat[co]*AK+0.4999f); ADat[co] = mcRg( qv ); } __global__ void kVmul( int AW, int AH, tBS_Val *ADat, float AK ) { mcEnter0(); ADat[co] = (int)(ADat[co]*AK+0.4999f); } 
void KVMul( tBS_Val *ADat, float AV ) { if (fabs(AV)>1.005f) kVMul<<<xG,xB>>>( xw, xh, ADat, AV ); else kVmul<<<xG,xB>>>( xw, xh, ADat, AV ); }
__global__ void kRng( int AW, int AH, tBS_Val *ADat, int AMn, int AMx ) { mcEnter0(); int qv = ADat[co]; if (qv<AMn) ADat[co] = AMn; else if (qv>AMx) ADat[co] = AMx; }
void KRng( tBS_Val* ADat, int AMn, int AMx ) { if (AMn<=-0x7FFF && AMx>=0x7FFF) return; kRng<<<xG,xB>>>( xw, xh, ADat, AMn, AMx ); } 
__global__ void kRngA( int AW, int AH, tBS_Val *ALap, int AMx ) { mcEnter0(); int qV = ALap[co], qA = abs(qV)-AMx; if (qA<=0) return; if (qA>AMx) qA = AMx; qA = AMx+(qA>>2); ALap[co] = qV<0?-qA:qA; };
__global__ void kRngA2( int AW, int AH, tBS_Val *ALap, int AMx, int AMxP ) { mcEnter0(); int qV = ALap[co]; if (-qV>AMx) { qV = -qV-AMx;
  if (qV>AMx) qV = AMx; ALap[co] = -AMx-(qV>>2); } else if (qV>AMxP) { qV -= AMxP; if (qV>AMxP) qV = AMxP; ALap[co] = AMxP+(qV>>2); }; }
void KRngA( tBS_Val* ALap, int AMx, int AMxP ) { if (AMxP<1 || AMxP>=0x7FFF) kRngA<<<xG,xB>>>( xw, xh, ALap, AMx ); else kRngA2<<<xG,xB>>>( xw, xh, ALap, AMx, AMxP ); }
__global__ void kRngStd( int AW, int AH, const int *ASrc, tBS_Val *ADst ) { mcEnter0(); ADst[co] = mcRg( ASrc[co] ); }
void KRngStd( const int *ASrc, tBS_Val* ADst ) { kRngStd<<<xG,xB>>>( xw, xh, ASrc, ADst ); }
__global__ void kAddRg( int AW, int AH, const tBS_Val *AAdd, tBS_Val *AAcc ) { mcEnter0(); int qv = AAdd[co]+AAcc[co]; AAcc[co] = mcRg( qv ); } void KAddRg( const tBS_Val *AAdd, tBS_Val *AAcc ) { kAddRg<<<xG,xB>>>( xw, xh, AAdd, AAcc ); }
__global__ void kAdd( int AW, int AH, const tBS_Val *Aadd, tBS_Val *Aacc ) { mcEnter0(); Aacc[co] += Aadd[co]; } void KAdd( const tBS_Val *Aadd, tBS_Val *Aacc ) { kAdd<<<xG,xB>>>( xw, xh, Aadd, Aacc ); }
__global__ void kAdd1( int AW, int AH, const tBS_Val *AAdd, int *AAcc ) { mcEnter0(); AAcc[co] += AAdd[co]; } __global__ void kSub1( int AW, int AH, const tBS_Val *AAdd, int *AAcc ) { mcEnter0(); AAcc[co] -= AAdd[co]; }
__global__ void kAddMul( int AW, int AH, const tBS_Val *AAdd, int *AAcc, float AV ) { mcEnter0(); AAcc[co] += (int)(AAdd[co]*AV+0.4999f); }
void KAddMul( const tBS_Val *AAdd, int *AAcc, float AV ) { if (fabs(AV)<1e-3f) return; if (fabs(AV-1)<1e-3f) { kAdd1<<<xG,xB>>>( xw, xh, AAdd, AAcc ); return; };
  if (fabs(AV+1)<1e-3f) { kSub1<<<xG,xB>>>( xw, xh, AAdd, AAcc ); return; }; kAddMul<<<xG,xB>>>( xw, xh, AAdd, AAcc, AV ); }
__global__ void kSubA( int AW, int AH, const tBS_Val *ASub, tBS_Val *AAcc ) { mcEnter0(); AAcc[co] -= ASub[co]; } void KSubA( const tBS_Val *ASub, tBS_Val *AAcc ) { kSubA<<<xG,xB>>>( xw, xh, ASub, AAcc ); } // Acc := Acc-Sub
__global__ void kSubS( int AW, int AH, const tBS_Val *ASub, tBS_Val *AAcc ) { mcEnter0(); AAcc[co] = ASub[co]-AAcc[co]; } void KSubS( const tBS_Val *ASub, tBS_Val *AAcc ) { kSubS<<<xG,xB>>>( xw, xh, ASub, AAcc ); } // Acc := Sub-Acc
__global__ void kMul( int AW, int AH, const tBS_Val *AMul, tBS_Val *ADat, int ANrm ) { mcEnter0(); int qV = AMul[co]*ADat[co]; ADat[co] = ANrm>0 ? qV/ANrm : qV>>(-ANrm); } void KMul( const tBS_Val *AMul, tBS_Val *ADat, int ANrm ) { kMul<<<xG,xB>>>( xw, xh, AMul, ADat, ANrm ); };
__global__ void kMxA( int AW, int AH, const tBS_Val *ACmp, tBS_Val *AAcc ) { mcEnter0(); int qV = ACmp[co]; if (abs(AAcc[co])<abs(qV)) AAcc[co] = qV; } void KMxA( const tBS_Val *ACmp, tBS_Val *AAcc ) { kMxA<<<xG,xB>>>( xw, xh, ACmp, AAcc ); };

//======== Gray-scale gamut's XLAT ==========
//===========================================
typedef struct {int *rX,rM,rG,rV[4];} tXLATe; const int cBS_XLatTD = 256; // norm_one for terminal_derivation (shall coincide with Ut0.h)
static tXLATe ixXLAT[4] = {{NULL,0,0,0,0,0,0},{NULL,0,0,0,0,0,0},{NULL,0,0,0,0,0,0},{NULL,0,0,0,0,0,0}}; // ~'array of CBS_XLat from _UT0.cpp'
void KXLAT_Release( int ASlot ) { tXLATe *qX = ixXLAT+ASlot; if (!qX->rX) return; cudaFree( qX->rX); qX->rX = NULL; } // -- release XLAT[Slot]

static void iTvk( int AIdxDir, const int *AX, int AGam, int ATD, int *AV ) { // determine terminal_derivation
  const int *cX = AIdxDir>0?AX+0:AX+AGam-1; AV[0] = cX[0]; // 'V','K' in XLATApply
  if (ATD>0) AV[1] = (cX[0]-cX[__min(ATD,AGam>>1)*AIdxDir])*cBS_XLatTD/ATD; else {
  if (AX[0]>AX[AGam-1]) AIdxDir *= -1; AV[1] = __max(-cBS_XLatTD,ATD)*AIdxDir; }; } // if descend then inverse
void KXLAT_Init( int ASlot, const int *AXLat, int AMn, int AGam, int ATD0, int ATD1 ) { KXLAT_Release( ASlot ); // -- allocate XLAT[Slot][0..Gam-1]
  tXLATe *qX = ixXLAT+ASlot; cudaMalloc( &qX->rX, AGam*sizeof(int) ); cudaMemcpy( qX->rX, AXLat, AGam*sizeof(int), cudaMemcpyHostToDevice );
  qX->rM = AMn; qX->rG = AGam; iTvk( 1, AXLat, AGam, ATD0, qX->rV+0 ); iTvk( -1, AXLat, AGam, ATD1, qX->rV+2 ); }

__global__ void kXLAT_Apply( int AW, int AH, tBS_Val *ADat, const int *AXLat, int cMn, int cGam1, int cV0, int cK0, int cV1, int cK1 ) { mcEnter0(); ADat += co; int qv = *ADat-cMn;
  if (qv<0) { qv = cV0-qv*cK0/cBS_XLatTD; *ADat = qv<=-0x7FFF?-0x7FFF:qv; return; };
  if (qv>cGam1) { qv = cV1+(qv-cGam1)*cK1/cBS_XLatTD; *ADat = qv>=0x7FFF?0x7FFF:qv; return; };
  *ADat = AXLat[qv]; }
void KXLAT_Apply( int ASlot, tBS_Val *ADat ) {
  tXLATe *qX = ixXLAT+ASlot; if (!qX->rX) return; kXLAT_Apply<<<xG,xB>>>( xw, xh, ADat, qX->rX, qX->rM, qX->rG-1, qX->rV[0], qX->rV[1], qX->rV[2], qX->rV[3] );}

__global__ void kGamHist( int AW, int AH, const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam1, const tBS_Val *ARef ) {
  const int cx = blockDim.x*blockIdx.x+threadIdx.x; if (cx>=AW) return; const int cy = blockDim.y*blockIdx.y+threadIdx.y; if (cy>=AH) return;
  int qo = cy*AStride+cx, qv = ADat[qo]; qv = mcRng( qv, 0, AGam1 ); atomicAdd( AHist+qv, ARef?ARef[qo]:1 ); }
void KGamHist( const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam, const tBS_Val *ARef ) {
  unsigned *qHist; cudaMalloc( &qHist, AGam*sizeof(unsigned) ); cudaMemset( qHist, 0, AGam*sizeof(unsigned) );
  kGamHist<<<xG,xB>>>( xw, xh, ADat, AStride, qHist, AGam-1, ARef );
  cudaMemcpy( AHist, qHist, AGam*sizeof(unsigned), cudaMemcpyDeviceToHost ); cudaFree( qHist ); }

static const int cBmn = 5; // min block's dim (as power of two)
#define mcCLAHEbl(B,Dim,N,D0) { N = Dim>>B; D0 = (Dim-(N<<B))/2; } // block's macro
typedef struct { int rNx,rNy,rX0,rY0,rSz; } tCLAHE; // hists[Ny][Nx]; skip X0 from lft, Y0 from top
static bool iCLAHE( int AW, int AH, int AB, int AGam, tCLAHE *AP ) { if (AB<cBmn || AB>cBmn+5) return false;
  mcCLAHEbl( AB, AW, AP->rNx, AP->rX0 ); mcCLAHEbl( AB, AH, AP->rNy, AP->rY0 ); AP->rSz = AGam*AP->rNx*AP->rNy; return AP->rNx>1 && AP->rNy>1; }
template<int cB> __global__ void kGamCLAHEHists( int cNx, int cNy, int cStride, const tBS_Val *ADat, unsigned cGam, unsigned *AHists, const tBS_Val *ARef ) {
  const int cy = blockDim.y*blockIdx.y+threadIdx.y, cNy0 = cy>>cB; if (cNy0>=cNy) return;
  const int cx = blockDim.x*blockIdx.x+threadIdx.x, cNx0 = cx>>cB; if (cNx0>=cNx) return;
  const int co = cy*cStride+cx; unsigned qV = ADat[co]; if (qV>=cGam) qV = cGam-1; atomicAdd( AHists+(cNy0*cNx+cNx0)*cGam+qV, ARef?ARef[co]:1 ); }
unsigned *KGamCLAHEHists( const tBS_Val *ADat, int AB, int AGam, const tBS_Val *ARef ) {
  tCLAHE qP; if (!iCLAHE( xw, xh, AB, AGam, &qP )) return NULL; // -- get Hists[Ny][Nx]
  unsigned *qHists; cudaMalloc( &qHists, qP.rSz*sizeof(unsigned) ); cudaMemset( qHists, 0, qP.rSz*sizeof(unsigned) );
  const int co = qP.rY0*xw+qP.rX0; ADat += co; if (ARef) ARef +=co; switch (AB) {
    case cBmn+0: kGamCLAHEHists<cBmn+0><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break;
    case cBmn+1: kGamCLAHEHists<cBmn+1><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break;
    case cBmn+2: kGamCLAHEHists<cBmn+2><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break;
    case cBmn+3: kGamCLAHEHists<cBmn+3><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break;
    case cBmn+4: kGamCLAHEHists<cBmn+4><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break;
    case cBmn+5: kGamCLAHEHists<cBmn+5><<<xG,xB>>>( qP.rNx, qP.rNy, xw, ADat, AGam, qHists, ARef ); break; };
  unsigned *qHistsH = new unsigned[qP.rSz];
  cudaMemcpy( qHistsH, qHists, qP.rSz*sizeof(unsigned), cudaMemcpyDeviceToHost ); 
  cudaFree( qHists ); return qHistsH; }

template<int cB> __device__ inline void mcCLAHEio( int N, int D, int P, int *IO ) { const int M = (1<<cB)+(1<<(cB-1))+D; P -= M; // index,offset-macro
  if (P<=0) { IO[0] = 0; IO[1] = ((P+M)<<cB)/M; return; }; P += 1<<cB; D = P>>cB;
  if (D>=N-2) { IO[0] = N-2; IO[1] = ((P-((N-2)<<cB))<<cB)/M; return; }; IO[0] = D; IO[1] = P-(D<<cB); }
template<int cB> __global__ void kGamCLAHEApply( int AW, int AH, int cNx, int cNy, tBS_Val *ADat, int cGam, const int *cHists, int cX0, int cY0 ) { mcEnter0(); 
  int  qIOy[2]; mcCLAHEio<cB>( cNy, cY0, cy, qIOy ); qIOy[0] *= cNx; int qIOx[2]; mcCLAHEio<cB>( cNx, cX0, cx, qIOx );
  int  qI = qIOy[0]+qIOx[0]; unsigned qV = ADat[co]; if (qV>=cGam) qV = cGam-1; const int *cHistV = cHists+qV;
  int  qV00 = cHistV[qI*cGam], qV01 = cHistV[(qI+1)*cGam], qV10 = cHistV[(qI+cNx)*cGam], qV11 = cHistV[(qI+cNx+1)*cGam];
  const int cc = (1<<(cB-1))-1;
  qV00 += ((qV01-qV00)*qIOx[1]+cc)>>cB; qV10 += ((qV11-qV10)*qIOx[1]+cc)>>cB; ADat[co] = qV00+(((qV10-qV00)*qIOy[1]+cc)>>cB); }
bool KGamCLAHEApply( tBS_Val *ADat, int AB, int AGam, int *AHistsH ) { // -- apply Hists[Ny][Nx]
  tCLAHE qP; if (!AHistsH || !iCLAHE( xw, xh, AB, AGam, &qP )) return false; int *qHists; cudaMalloc( &qHists, qP.rSz*sizeof(int) );
  cudaMemcpy( qHists, AHistsH, qP.rSz*sizeof(int), cudaMemcpyHostToDevice ); switch (AB) {
    case cBmn+0: kGamCLAHEApply<cBmn+0><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break;
    case cBmn+1: kGamCLAHEApply<cBmn+1><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break;
    case cBmn+2: kGamCLAHEApply<cBmn+2><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break;
    case cBmn+3: kGamCLAHEApply<cBmn+3><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break;
    case cBmn+4: kGamCLAHEApply<cBmn+4><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break;
    case cBmn+5: kGamCLAHEApply<cBmn+5><<<xG,xB>>>( xw, xh, qP.rNx, qP.rNy, ADat, AGam, qHists, qP.rX0, qP.rY0 ); break; };
  cudaFree( qHists ); return true; }

//=========== Skew/Shear ============
//===================================
__global__ void kShear( int AW, int AH, int AFill, float AShear, int AW0, int AH0, const tBS_Val *A0, tBS_Val *A1 ) {
  mcEnter0(); int cShearInt, qBuf4[4]; bool  qBuf4Ok; // Shear>0:"rotate CCW", <0:"CW"
  if (AH0) { const tBS_Val *cCol0 = A0+cx; // vertical; W0=W1
    AShear *= cx-(AW>>1); cShearInt = (int)floor( AShear ); // keep 'sign(Shear)'
    int q0 = cy+(AH0>>1)-(AH>>1)+cShearInt; 
    qBuf4Ok = mcHerLikBuf4( AFill, AW, AH, cCol0, q0, qBuf4 );
  } else { const tBS_Val *cRow0 = A0+cy*AW0; // horizontal; H0=H1
    AShear *= (AH>>1)-cy; cShearInt = (int)floor( AShear ); // invert 'sign(Shear)'
    int q0 = cx+(AW0>>1)-(AW>>1)+cShearInt;
    qBuf4Ok = mcHerLikBuf4( AFill,  1, AW, cRow0, q0, qBuf4 );
  };
  if (qBuf4Ok) qBuf4[0] = (int)mcHerLikStdSpline( AFill, AShear-cShearInt, qBuf4 );
  A1[co] = mcRg( qBuf4[0] ); }

void KShear( int AFill, float AShear, int AW0, int AH0, const tBS_Val *A0, tBS_Val *A1 ) {
  kShear<<<xG,xB>>>( xw, xh, AFill, AShear, AW0, AH0, A0, A1 ); } 
 
//========== Directional average  ===========
//===========================================
static struct { cudaArray *rA; cudaTextureObject_t rT; } ixTue = { NULL }; // 2D-buffer and texture for boundary checks in kernels
bool KTueInit( const tBS_Val *AInp ) { if (!ixTue.rA) { 
  cudaChannelFormatDesc qDesc = cudaCreateChannelDesc<tBS_Val>(); if (cudaMallocArray( &ixTue.rA, &qDesc, xw, xh )!=cudaSuccess) return false;
  cudaResourceDesc qTR; memset( &qTR, 0, sizeof(qTR) ); qTR.resType = cudaResourceTypeArray; qTR.res.array.array = ixTue.rA;
  cudaTextureDesc qTD; memset( &qTD, 0, sizeof(qTD) ); if (cudaCreateTextureObject( &ixTue.rT, &qTR, &qTD, NULL )!=cudaSuccess) return false; }; // if not Tue.A
  return cudaMemcpy2DToArray( ixTue.rA, 0, 0, AInp, xw*sizeof(tBS_Val), xw*sizeof(tBS_Val), xh, cudaMemcpyDeviceToDevice )==cudaSuccess; }
void KTueRelease() { if (!ixTue.rA) return; cudaFreeArray( ixTue.rA ); ixTue.rA = NULL; cudaDestroyTextureObject( ixTue.rT ); } // DestroyTextureObject(xTue)

// directional averages' macros; 'mckAv'~-'CPU'; 'mckAu'~='more sharp'
#define mcImgPch(X,Y) tex2D<tBS_Val>( AInp, cx+(X)+0.5f, cy+(Y)+0.5f ) // single texture's pixel at 'cx+X','cy+Y'
#define cAvLg  5
#define cAviC ((1<<(cAvLg-3))+0) // corrective Y-decrement to be like moving averages CPU
#define cAvi1 (1<<(cAvLg-1))
#define cAvi2 (1<<cAvLg)
template<int cR,int Dx,int Dy,int Dys> __global__ void mckAv( int AW, int AH, cudaTextureObject_t AInp, tBS_Val *AOut ) { mcEnter0();
  int qSum = mcImgPch(0,0); for (int qi=1;qi<=cR;qi++) { int qx = (qi*Dx+cAvi1-1)>>cAvLg, qy = Dys*((qi*Dy+cAvi1-1)>>cAvLg);
  qSum += mcImgPch(qx,qy)+mcImgPch(-qx,-qy); }; AOut[co] = qSum/((cR<<1)+1); }
#if 1 // ~less sharp than CPU
template<int Dx,int Dy,int Dys> inline void mcKAv( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { if (!KTueInit( AInp )) return; switch (AR) { // increment 0 or 1 by X and Y
  case  1:mckAv< 1,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  2:mckAv< 2,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  3:mckAv< 3,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  4:mckAv< 4,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  5:mckAv< 5,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  6:mckAv< 6,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  7:mckAv< 7,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  8:mckAv< 8,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  9:mckAv< 9,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 10:mckAv<10,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 11:mckAv<11,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 12:mckAv<12,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 13:mckAv<13,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 14:mckAv<14,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 15:mckAv<15,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; default:mckAv<16,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; }; }
#else // ~more sharp than CPU
template<int cR,int Dx,int Dy,int Dys> __global__ void mckAu( int AW, int AH, cudaTextureObject_t AInp, tBS_Val *AOut ) { mcEnter0();
  int qSum = mcImgPch(0,0)*3; for (int qi=1;qi<=cR;qi++) { int qx = (qi*Dx+cAvi1-1)>>cAvLg, qy = Dys*((qi*Dy+cAvi1-1)>>cAvLg),
  qV = mcImgPch(qx,qy)+mcImgPch(-qx,-qy); if (qi<cR) qV <<= 1; qSum += qV; }; AOut[co] = qSum/(((cR-1)<<2)+5); }
template<int Dx,int Dy,int Dys> inline void mcKAv( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { if (!KTueInit( AInp )) return; switch (AR) { // increment 0 or 1 by X and Y
  case  1:mckAv< 1,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  2:mckAv< 2,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  3:mckAv< 3,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  4:mckAv< 4,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  5:mckAu< 5,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  6:mckAu< 6,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  7:mckAu< 7,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case  8:mckAu< 8,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case  9:mckAu< 9,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 10:mckAu<10,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 11:mckAu<11,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 12:mckAu<12,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 13:mckAu<13,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 14:mckAu<14,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 15:mckAu<15,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; default:mckAu<16,Dx,Dy,Dys><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; }; }
#endif
// particular directional averages
void BC_Av1( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi2,    0, 1>(AInp,AOut,AR); }; // '-' #1
void BC_Av3( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi2,cAvi2, 1>(AInp,AOut,AR); }; // '\' #3
void BC_Av5( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<0,    cAvi2, 1>(AInp,AOut,AR); }; // '|' #5
void BC_Av7( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi2,cAvi2,-1>(AInp,AOut,AR); }; // '/' #7
void BC_Av2( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi2,cAvi1-cAviC, 1>(AInp,AOut,AR); }; // #1..#3
void BC_Av4( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi1-cAviC,cAvi2, 1>(AInp,AOut,AR); }; // #3..#5
void BC_Av6( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi1-cAviC,cAvi2,-1>(AInp,AOut,AR); }; // #5..#7
void BC_Av8( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcKAv<cAvi2,cAvi1-cAviC,-1>(AInp,AOut,AR); }; // #7..#1

//========== Directional |Der|  ===========
//=========================================
template<int cR,int Dx,int Dy> __global__ void mckADrR( int AW, int AH, cudaTextureObject_t AInp, tBS_Val *AOut ) { mcEnter0(); int qSum = mcImgPch(0,0);
  int qx = (cR*Dx)>>1, qy = (cR*Dy)>>1; AOut[co] = abs( mcImgPch(qx,qy)-mcImgPch(-qx,-qy) ); }
template<int Dx,int Dy> inline void mcADr( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { if (!KTueInit( AInp )) return; switch (AR) { // increment 0 or 1 by X and Y
  case 2:mckADrR<2,Dx,Dy><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 3:mckADrR<3,Dx,Dy><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  case 4:mckADrR<4,Dx,Dy><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break; case 5:mckADrR<5,Dx,Dy><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); break;
  default:mckADrR<1,Dx,Dy><<<xG,xB>>>( xw, xh, ixTue.rT, AOut ); }; }
void BC_ADr1( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<2, 0>(AInp,AOut,AR); }; // '-' #1
void BC_ADr3( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<2, 2>(AInp,AOut,AR); }; // '\' #3
void BC_ADr5( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<0, 2>(AInp,AOut,AR); }; // '|' #5
void BC_ADr7( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<2,-2>(AInp,AOut,AR); }; // '/' #7
void BC_ADr2( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<2, 1>(AInp,AOut,AR); }; // #1..#3
void BC_ADr4( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<1, 2>(AInp,AOut,AR); }; // #3..#5
void BC_ADr6( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<1,-2>(AInp,AOut,AR); }; // #5..#7
void BC_ADr8( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { mcADr<2,-1>(AInp,AOut,AR); }; // #7..#1

//======= application spatial filters =======
//===========================================
template<int cR> __device__ inline void mcLoadSquareBuf( int AW, const tBS_Val *ASrc, tBS_Val *ADst ) { const int cD = cR*2+1; // load square buffer
  for (int qi=-cR;qi<=cR;qi++) for (int qj=-cR;qj<=cR;qj++) ADst[(qi+cR)*cD+qj+cR] = ASrc[qi*AW+qj]; }
template<int cR> __global__ void kAveNxN( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter2( cR ); const int cN = cR*2+1;
  tBS_Val qB[cN*cN]; mcLoadSquareBuf<cR>( AW, AInp, qB ); AOut[0] = mcAve<tBS_Val,cN*cN>( qB ); }
void KAve3x3( const tBS_Val *AInp, tBS_Val *AOut ) { kAveNxN<1><<<xG,xB>>>( xw, xh, AInp, AOut ); } 
void KAve5x5( const tBS_Val *AInp, tBS_Val *AOut ) { kAveNxN<2><<<xG,xB>>>( xw, xh, AInp, AOut ); }

template<typename tN, int cN, int cN2> __device__ inline void mcSortBAs( tN *B ) { for (int qi=0;qi<cN2;qi++) { // bubble sort (find min, put the min in sorted order)
  tN qM = B[qi]; int qm = qi; for (int qj=qi+1;qj<cN;qj++) if (qM>B[qj]) { qM = B[qj]; qm = qj; }; tN qT = B[qi]; B[qi] = qM; B[qm] = qT; }; } // for i 
template<int cR> __global__ void kMedNxN( int AW, int AH, const tBS_Val *AInp, tBS_Val *AOut ) { mcEnter2( cR ); const int cN = (cR<<1)+1;
  tBS_Val qB[cN*cN]; mcLoadSquareBuf<cR>( AW, AInp, qB ); mcSortBAs<tBS_Val,cN*cN,(cN*cN>>1)+1>( qB ); AOut[0] = qB[cN*cN>>1]; }
bool KMed( const tBS_Val *AInp, tBS_Val *AOut, int AR ) { switch (AR) {
  case 1: kMedNxN<1><<<xG,xB>>>( xw, xh, AInp, AOut ); break; case 2: kMedNxN<2><<<xG,xB>>>( xw, xh, AInp, AOut ); break;
  case 3: kMedNxN<3><<<xG,xB>>>( xw, xh, AInp, AOut ); break; case 4: kMedNxN<4><<<xG,xB>>>( xw, xh, AInp, AOut ); break;
  case 5: kMedNxN<5><<<xG,xB>>>( xw, xh, AInp, AOut ); break; case 6: kMedNxN<6><<<xG,xB>>>( xw, xh, AInp, AOut ); break;
  case 7: kMedNxN<7><<<xG,xB>>>( xw, xh, AInp, AOut ); break; case 8: kMedNxN<8><<<xG,xB>>>( xw, xh, AInp, AOut ); break;
  case 9: kMedNxN<9><<<xG,xB>>>( xw, xh, AInp, AOut ); break; default: return false; }; return true; }

//======= Radon-Shannon-denoise =======
//=====================================
// imccpb_Noise in BS_Ut0.cpp for 'BC_Noise' in BS_Ut1.cpp
__global__ void kNoise( int AW, int AH, int cThrLg, int cReg, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3,
  const tBS_Val *A4, const tBS_Val *A5, const tBS_Val *A6, const tBS_Val *A7, tBS_Val *ANoise ) { // invoke ~"mcRadonShannNoise in BS_Ut0.h"; 
   mcEnter0(); tBS_Val qDat8[8] = { A0[co], A1[co], A2[co], A3[co], A4[co], A5[co], A6[co], A7[co] };
  // inp:eight values "Blr(Org,rad,dir]-Org"; out:"basic RadonShannon_noise.
  // "RadonDeviation~=0":max_noise (i.e subtract thick lap), ">=1 shl ThrLg":min_noise (i.e. subtract thin lap)
  int qD = 0; for (int qi2=0;qi2<4;qi2++) qD += abs(qDat8[qi2]-qDat8[qi2+4]);
  if (!qD) { ANoise[co] = qDat8[0]; return; }; int qN = 8-(qD>>cThrLg);
  if (qN==8) { ANoise[co] = mcSum<tBS_Val,8>( qDat8 )/(8+cReg); return; };
  if (qN<1) qN = 1; int qZ = 0; for (int qi4=0;qi4<qN;qi4++) { int qMnJ = 0, qMn = qDat8[qi4];
    for (int qj=qi4+1;qj<8;qj++) { int qV = qDat8[qj]; if (abs(qMn)>abs(qV)) { qMn = qV; qMnJ = qj; }; }; // for j
    if (qMnJ) qDat8[qMnJ] = qDat8[qi4]; qZ += qMn; }; ANoise[co] = qZ/(qN+cReg); } // for i4
void KNoise( int AThrLg, int AReg, tBS_Val **A8, tBS_Val *ANoise ) {
  kNoise<<<xG,xB>>>( xw, xh, AThrLg, AReg, A8[0], A8[1], A8[2], A8[3], A8[4], A8[5], A8[6], A8[7], ANoise ); }

// "if |Lap[i]|<Z then Lap[i]>>=2"
__global__ void kMxLapRngZ( int AW, int AH, tBS_Val *ALap, int AZ ) { mcEnter0(); 
  int qV = ALap[co], qA = abs(qV); if (qA<AZ) ALap[co] = (qV<0)?-(qA>>2):qA>>2; }
void KMxLapRngZ( int AZ, tBS_Val *ALap ) { if (AZ<=0) return; kMxLapRngZ<<<xG,xB>>>( xw, xh, ALap, AZ ); }

// reduce shot_noise after RadonShannon_noise
#define mcMLMx( M, A ) if (abs(M)<abs(A[co])) M = A[co]
__global__ void kMxLapStd4( int AW, int AH, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3, // mxA of 4 dir_layers
  tBS_Val *ALap ) { mcEnter0(); int qM = A0[co]; mcMLMx( qM, A1 ); mcMLMx( qM, A2 ); mcMLMx( qM, A3 ); ALap[co] = qM; } 
__global__ void kMxLapStd8( int AW, int AH, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3, // mxA of 8 dir_layers
  const tBS_Val *A4, const tBS_Val *A5, const tBS_Val *A6, const tBS_Val *A7, tBS_Val *ALap ) { mcEnter0();
  int qM0 = A0[co]; mcMLMx( qM0, A2 ); mcMLMx( qM0, A4 ); mcMLMx( qM0, A6 );
  int qM1 = A1[co]; mcMLMx( qM1, A3 ); mcMLMx( qM1, A5 ); mcMLMx( qM1, A7 );
  ALap[co] = (qM0+qM1)>>1; }
void KMxLapStd( int ANDD, tBS_Val **A, tBS_Val *ALap ) { switch (ANDD) {
  case 4:  kMxLapStd4<<<xG,xB>>>( xw, xh, A[0], A[1], A[2], A[3], ALap ); break;
  default: kMxLapStd8<<<xG,xB>>>( xw, xh, A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], ALap ); }; }

//============== Reductions =================
//===========================================
template <typename tInp,typename tAcc> __global__ void mcReduceNaiveDevice( const tInp *AInp, tAcc *AOut, int AN ) {
  extern tAcc __shared__ xAcc[]; // perform first level of reduction,
  // reading from global memory, writing to shared memory
  int tid = threadIdx.x, bid = blockIdx.x*(blockDim.x<<1)+tid;
  tAcc qAcc = (bid<AN) ? AInp[bid]:0; if (bid+blockDim.x<AN) qAcc += AInp[bid+blockDim.x];
  xAcc[tid] = qAcc; __syncthreads(); // do reduction in shared mem
  for (int qi=blockDim.x>>1;qi>0;qi>>=1) { if (tid<qi) { qAcc += xAcc[tid+qi]; xAcc[tid] = qAcc; }; __syncthreads(); };
  if (tid==0) AOut[blockIdx.x] = qAcc; } // write result for this block to global mem
template <typename tInp,typename tAcc,typename tAcZ> cudaError_t mcReduceNaiveHost(
  const tInp *AInp, tAcc *ATmp, tAcc *ATmpH, int AN, int ANt, double *ARes ) {
  cudaError_t qErr; const int cNt2 = ANt*2;
  if (AN<cNt2) { tInp *qTmpH = (tInp*)ATmpH;
    qErr = cudaMemcpy( qTmpH, AInp, AN*sizeof(tInp), cudaMemcpyDeviceToHost );
    if (qErr!=cudaSuccess) return qErr;
    tAcZ qAcZ = 0; for (int qi=0;qi<AN;qi++) qAcZ += qTmpH[qi];
    *ARes = (double)qAcZ/AN; return cudaSuccess; };
  const int cNb = (AN+cNt2-1)/cNt2;
  dim3 qGrids(cNb,1,1), qBlocks(ANt,1,1); // reduce 'at 1st level'
  int  qSharMem = ANt*sizeof(tAcc); if (ANt<=32) qSharMem <<= 1;
  mcReduceNaiveDevice<tInp,tAcc><<<qGrids,qBlocks,qSharMem>>>( AInp, ATmp, AN );
  qErr = cudaGetLastError(); if (qErr!=cudaSuccess) return qErr; // Kernel execution failed
  qErr = cudaMemcpy( ATmpH, ATmp, cNb*sizeof(tAcc), cudaMemcpyDeviceToHost );
  if (qErr!=cudaSuccess) return qErr;
  tAcZ qAcZ = 0; for (int qi=0;qi<cNb;qi++) qAcZ += ATmpH[qi];
  *ARes = (double)qAcZ/AN; return cudaSuccess; }
double KMean( const tBS_Val *AInp, tBS_Val *ATmp, tBS_Val *ATmpH ) { double qAve = 0; // (because of ambiguous comparison) avoid NAN here 
  return (mcReduceNaiveHost<tBS_Val,long,tLL>( AInp, (long*)ATmp, (long*)ATmpH, xw*xh, xBC_PMxBlkThreads, &qAve ))!=cudaSuccess?-1e37f:qAve; }

//======== kernels beside _Ut01.cpp =========
//===========================================

// BS_Ut2.cpp cpb_RadonShannon for 'BS3_RadonMM'
#define mcRadonLd( U0, U1 ) { qV0 = U0[co]; qV1 = U1[co]; }
#define mcRadonMM( V ) { if (qMn>V) qMn = V; else if (qMx<V) qMx = V; }
#define mcRadonSn() { int d = abs(qV0-qV1); if (d>255) d = 255; qD += d; qZ += (qV0+qV1)*d; }
__global__ void kRadonShannon( int AW, int AH, tBS_Val *A0, tBS_Val *A1, tBS_Val *A2, tBS_Val *A3, tBS_Val *A4,
  tBS_Val *A5, tBS_Val *A6, tBS_Val *A7, const tBS_Val *AOrg, tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp ) {
  mcEnter0(); int qV0, qV1, qMn, qMx, qD = 0, qZ = 0;
  mcRadonLd( A0, A4 ); qMn = qMx = qV0;  mcRadonMM( qV1 ); mcRadonSn();
  mcRadonLd( A1, A5 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonSn();
  mcRadonLd( A2, A6 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonSn();
  mcRadonLd( A3, A7 ); mcRadonMM( qV0 ); mcRadonMM( qV1 ); mcRadonSn();
  ADev[co] = qD; ABlr[co] = qZ/((qD<<1)+1);
  const int cO = AOrg[co]; if (abs(qMx)<abs(qMn)) qMx = qMn;
  ASrp[co] = abs(qMx)<abs(cO) ? (qMx*3+cO)>>2 : (qMx+cO)>>1; }
void KRadonShannon( tBS_Val **A8, const tBS_Val *AOrg, tBS_Val *ADev, tBS_Val *ABlr, tBS_Val *ASrp ) {
  kRadonShannon<<<xG,xB>>>( xw, xh, A8[0], A8[1], A8[2], A8[3], A8[4], A8[5], A8[6], A8[7], AOrg, ADev, ABlr, ASrp ); }

__device__ inline int mcMimic( int ADat, int ALap ) { int qL = ALap, qV = abs(ADat); if (qL<0) qV = -qV; else if (!qL) qV = 0; return qV; }
__global__ void kMimicLap( int AW, int AH, const tBS_Val *ADat, tBS_Val *ALap ) { mcEnter0(); ALap[co] = mcMimic( ADat[co], ALap[co] ); }
__global__ void kMimicDat( int AW, int AH, const tBS_Val *ALap, tBS_Val *ADat ) { mcEnter0(); ADat[co] = mcMimic( ADat[co], ALap[co] ); }
__global__ void kVMimicNZ( int AW, int AH, tBS_Val *ADat, int AV ) { mcEnter0(); if (ADat[co]) ADat[co] = AV; }
void KMimicLap( const tBS_Val *ADat, tBS_Val *ALap ) { kMimicLap<<<xG,xB>>>( xw, xh, ADat, ALap ); }
void KMimicDat( const tBS_Val *ALap, tBS_Val *ADat ) { kMimicDat<<<xG,xB>>>( xw, xh, ALap, ADat ); }
void KVMimicNZ( tBS_Val *ADat, int AV ) { kVMimicNZ<<<xG,xB>>>( xw, xh, ADat, AV ); }

__global__ void kRngLapByConductor( int AW, int AH, const tBS_Val *ACnd, tBS_Val *ALap ) {
  mcEnter0(); int qC = ACnd[co], qL = ALap[co]; if (qC<abs(qL)) ALap[co] = mcsign( qC, qL ); }
void KRngLapByConductor( const tBS_Val *ACnd, tBS_Val *ALap ) { kRngLapByConductor<<<xG,xB>>>( xw, xh, ACnd, ALap ); }

__global__ void kRidgeByConductors( int AW, int AH, tBS_Val *ARdg, const tBS_Val *ACDev, const tBS_Val *ACDer0, const tBS_Val *ACDer1 ) {
  mcEnter0();  int cDev = ACDev[co], cDer = __max(ACDer0[co],ACDer1[co]), qV = cDev-cDer; ARdg[co] = qV<=0 ? 0 : mcsign(qV,ARdg[co]); }
void KRidgeByConductors( tBS_Val *ARdg, const tBS_Val *ACDev, const tBS_Val *ACDer0, const tBS_Val *ACDer1 ) {
  kRidgeByConductors<<<xG,xB>>>( xw, xh, ARdg, ACDev, ACDer0, ACDer1 ); }

// MxLapStd with ~extended_limbo; BS_Ut3.cpp,_MxLapXLimbo
#define mcMxLapXLimbo1( A ) { int v = A[co]; if (abs(qM)<abs(v)) qM = v; qA += v; } // ~common maxlap
__device__ inline int imcMxLapXLimbo2( int L, int S, int M, int A, int K, int TG ) {
  M = S+((M-S)*K>>10); A = S+(int)((A-S)*K>>10); return M<A ? mcRng(L,M-TG,A+TG) : mcRng(L,A-TG,M+TG); }
__global__ void kMxLapXLimbo4( int AW, int AH, tBS_Val *ALap, int AK, int ATolGap, tBS_Val *AStd,
  const tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3 ) { mcEnter0();
  int qM = A0[co], qA = qM; mcMxLapXLimbo1( A1 ); mcMxLapXLimbo1( A2 ); mcMxLapXLimbo1( A3 ); qA >>= 1;
  ALap[co] = imcMxLapXLimbo2( ALap[co], AStd[co], qM, qA, AK, ATolGap ); }
__global__ void kMxLapXLimbo8( int AW, int AH, tBS_Val *ALap, int AK, int ATolGap, tBS_Val *AStd, const tBS_Val *A0, const tBS_Val *A1,
  const tBS_Val *A2, const tBS_Val *A3, const tBS_Val *A4, const tBS_Val *A5, const tBS_Val *A6, const tBS_Val *A7 ) { mcEnter0(); 
  int qM = A0[co], qA = qM; mcMxLapXLimbo1( A1 ); mcMxLapXLimbo1( A2 ); mcMxLapXLimbo1( A3 );
  mcMxLapXLimbo1( A4 ); mcMxLapXLimbo1( A5 ); mcMxLapXLimbo1( A6 ); ; mcMxLapXLimbo1( A7 ); qA >>= 2;
  ALap[co] = imcMxLapXLimbo2( ALap[co], AStd[co], qM, qA, AK, ATolGap ); }
void KMxLapXLimbo( tBS_Val *ALap, float AK, int ATolGap, tBS_Val **ATs, int ANDD ) {
  // calc virtual low_est, high_est relative to the basis (which is MxLapStd)
  // always "Ave<MxLapStd<MxLapExact" or "MxLapExact<MxLapStd<Ave"
  if (ANDD==4) kMxLapXLimbo4<<<xG,xB>>>( xw, xh, ALap, (int)(AK*1024), ATolGap, ATs[0], ATs[1], ATs[2], ATs[3], ATs[4] ); else
    kMxLapXLimbo8<<<xG,xB>>>( xw, xh, ALap, (int)(AK*1024), ATolGap, ATs[0], ATs[1], ATs[2], ATs[3], ATs[4], ATs[5], ATs[6], ATs[7], ATs[8] ); }

// AutoOscillator here reflects "mccpb_AutOscMsk{... return qZ>>(cNDDLg-1)}" in BS_Ut2.cpp
#define mcMxLapXAutOsc( A0, A1 ) { int q0 = A0[co], q1 = A1[co]; qZ += abs(q0)+abs(q1)+abs(q0-q1); } // AutoOscillator's mask
__global__ void kMxLapXAutOscMsk4( int AW, int AH, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3 ) {
  mcEnter0(); int qZ = 0; mcMxLapXAutOsc( A0, A2 ); mcMxLapXAutOsc( A1, A3 ); A0[co] = qZ>>1; } 
__global__ void kMxLapXAutOscMsk8( int AW, int AH, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3,
  const tBS_Val *A4, const tBS_Val *A5, const tBS_Val *A6, const tBS_Val *A7 ) { mcEnter0(); int qZ = 0;
  mcMxLapXAutOsc( A0, A4 ); mcMxLapXAutOsc( A1, A5 ); mcMxLapXAutOsc( A2, A6 ); mcMxLapXAutOsc( A3, A7 ); A0[co] = qZ>>2; }
void KMxLapAutOscMsk( tBS_Val **ATs, int ANDD ) {
  if (ANDD==4) kMxLapXAutOscMsk4<<<xG,xB>>>( xw, xh, ATs[0], ATs[1], ATs[2], ATs[3] ); else
  kMxLapXAutOscMsk8<<<xG,xB>>>( xw, xh, ATs[0], ATs[1], ATs[2], ATs[3], ATs[4], ATs[5], ATs[6], ATs[7] ); }

// general despeckle (refers both median and average); BS_Ut6.cpp,_DespeckleCUDA
__global__ void kDespeckle( int AW, int AH, int AThr, tBS_Val *ADat, const tBS_Val *AMed, const tBS_Val *AAve ) { mcEnter0();
  int qV = ADat[co], qV1 = AMed[co], qV2 = AAve[co]; qV1 -= qV; qV2 -= qV;
  int qD = abs(qV1)<abs(qV2)?qV1:qV2, qA = abs(qD)-AThr; if (qA<=0) return;
  if (qA>AThr) qA = AThr; qA >>= 2; if (qD>0) qA = -qA; ADat[co] = qV+qD+qA; }
void KDespeckle( int AThr, tBS_Val *ADat, const tBS_Val *AMed, const tBS_Val *AAve ) { kDespeckle<<<xG,xB>>>( xw, xh, AThr, ADat, AMed, AAve ); }

// BS_Ut6.cpp mccpb_RadonDespeckle; BS_Ut6.cpp,_DePoisson{iRadonDespeckle_CUDA}
#define mcRadonM1( A ) { int qv = A[co]; if (qMn>qv) qMn = qv; else if (qMx<qv) qMx = qv; qAv += qv; }
__global__ void kRadonDespeckle( int AW, int AH, float ALimK, tBS_Val *ADat, tBS_Val *AAve,
  tBS_Val *A0, tBS_Val *A1, tBS_Val *A2, tBS_Val *A3, tBS_Val *A4, tBS_Val *A5, tBS_Val *A6, tBS_Val *A7 ) {
  mcEnter0(); int qMn = A0[co], qMx = qMn, qAv = qMn, qNDD;
  mcRadonM1( A1 ); if (!A2) { qNDD = 2; goto labCalc; };
  mcRadonM1( A2 ); mcRadonM1( A3 ); if (!A4) { qNDD = 4; goto labCalc; };
  mcRadonM1( A4 ); mcRadonM1( A5 ); if (!A6) { qNDD = 6; goto labCalc; };
  mcRadonM1( A6 ); mcRadonM1( A7 ); qNDD = 8;
  labCalc: qAv /= qNDD; int qV0 = ADat[co];
  qMn = qAv+(int)(ALimK*(qMn-qAv)); if (qV0<=qMn) { qV0 = qMn; goto labLeave; };
  qMx = qAv+(int)(ALimK*(qMx-qAv)); if (qV0>qMx) qV0 = qMx; labLeave:
  ADat[co] = qV0-qAv; AAve[co] = qAv; }
void KRadonDespeckle( float ALimK, tBS_Val *ADat, tBS_Val *AAve, tBS_Val **AA ) {
  kRadonDespeckle<<<xG,xB>>>( xw, xh, ALimK, ADat, AAve, AA[0], AA[1], AA[2], AA[3], AA[4], AA[5], AA[6], AA[7] ); }

// BS_Ut6.cpp icpl_ScrapeGrid
// Empirically emphasizes symmetric laplacians ~"--+--", "-+++-";
// doesn't emphasize non-symmetric laplacians ~"+++--", "--+++".
__device__ inline int mcScrapeGridIdx( int AR, int ADX, int ADY, int AW ) { return AR*(ADX+ADY*AW); }
__device__ inline int mcScrapeGridAdd( int AR, int A0, int A1 ) { A0 += A1; A0 -= A0>>(3-AR); return A0; }
__device__ int iScrapeGrid1( int ADX, int ADY, int AW, const tBS_Val *ADat, int *AMu ) {
  int qv = ADat[0], qu = 0; if (!qv) goto labLeave; for (int qr=1;qr<3;qr++) {
    int q0 = ADat[mcScrapeGridIdx(-qr,ADX,ADY,AW)], q1 = ADat[mcScrapeGridIdx(qr,ADX,ADY,AW)];
    if (qv<0) { if (q0<=0 || q1<=0) continue; q0 = mcScrapeGridAdd(qr,q0,q1); qu += __max(qv,-q0); }
         else { if (q0>=0 || q1>=0) continue; q0 = mcScrapeGridAdd(qr,q0,q1); qu += __min(qv,-q0); }; }; // for
  labLeave: if (abs(*AMu)>qu) *AMu = qu; return qu; } 
__device__ int iScrapeGrid4( int AW, const tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2, const tBS_Val *A3 ) {
  int M = 0x7FFFFFFF, Z = iScrapeGrid1( 1, 0, AW, A0, &M )+iScrapeGrid1( 1, 1, AW, A1, &M )+
    iScrapeGrid1( 0, 1, AW, A2, &M )+iScrapeGrid1( 1, -1, AW, A3, &M ); return (Z-M)>>1; }
__global__ void kScrapeGrid( int AW, int AH, tBS_Val *A0, const tBS_Val *A1, const tBS_Val *A2,
  const tBS_Val *A3, const tBS_Val *A4, const tBS_Val *A5, const tBS_Val *A6, const tBS_Val *A7 ) {
  mcEnter0(); if ((cx<2)||(cx>=AW-2)||(cy<2)||(cy>=AH-2)) { A0[co] = 0; return; }; 
  // summary hypo_grid estimates as ~"{scrape(R=1)*3+scrape(R=2)*2}/5"
  A0[co] = (iScrapeGrid4(AW,A0+co,A1+co,A2+co,A3+co)*3+iScrapeGrid4(AW,A4+co,A5+co,A6+co,A7+co)*2)/5; }
void KScrapeGrid( tBS_Val **AA ) { // estimate "summary hypothetical grid"
  // inp #0..3:lap(R=1), #4..7:lap(R=2); out #0:sum_hypo_grid
  kScrapeGrid<<<xG,xB>>>( xw, xh, AA[0], AA[1], AA[2], AA[3], AA[4], AA[5], AA[6], AA[7] ); }
