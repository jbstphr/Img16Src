// I16_Ut1.cpp : CUDA interconnect utilities
// the routines here mostly operate with fixed frames BC_Dim()
// by S.B.
#include "stdafx.h"
#include "math.h"
#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include "BS_Cmn.h"
#include "BS_Ut0.h"
#include "BS_Ut1.h"
#include "BS_Ut2a.h" // Skew/Shear

//================= commons ==================
//============================================
bool xBC_Err = false;
static const char *iErr( cudaError_t AStat ) {
  if (AStat==cudaSuccess) { xBC_Err = false; return NULL; }; xBC_Err = true; return cudaGetErrorString( AStat ); }
static void ikLeave( cudaError_t AStat, const char *AProcName ) {
  const char *cErr = iErr( AStat ); if (cErr) { char qS[128]; sprintf( qS, "%s,%s", AProcName, cErr ); ErrMsg( qS ); }; }
static void ikLeave( const char *AProcName ) { ikLeave( cudaGetLastError(), AProcName ); }

int xBC_PCnt = 0, xBC_PMxBlk = 0, xBC_PMxBlkThreads = 0; // MxBlk is 'rational number' for reference
void BC_DeviceFind() {
  xBC_PCnt = xBC_PMxBlk = xBC_PMxBlkThreads = 0; const int cDevice = 0; // Set the device to be used
  ikLeave( cudaSetDevice(cDevice), "BC_DeviceFind" ); if (xBC_Err) return;
  cudaDeviceProp qP = { 0 }; ikLeave( cudaGetDeviceProperties(&qP,cDevice),"BC_DeviceFind|Prop" ); if (xBC_Err) return;
  //if (!qP.cooperativeLaunch) { ErrMsg( "BC_DeviceFind|Prop does not report Cooperative Kernel Launch" ); xBC_Err = true; return; };
  xBC_PCnt = qP.multiProcessorCount; xBC_PMxBlkThreads = qP.maxThreadsPerBlock; xBC_PMxBlk = xBC_PCnt*(qP.maxThreadsPerMultiProcessor/xBC_PMxBlkThreads);
  if (1) fprintf( stderr, cTPS" BC_DeviceFind Prc=%d*%d Blk=%d*%d\n", xBC_PCnt, qP.maxThreadsPerMultiProcessor, xBC_PMxBlk, xBC_PMxBlkThreads ); }
void BC_DeviceReset() { ikLeave( cudaDeviceReset(), "BC_DeviceReset" ); }
void BC_Sync() { ikLeave( cudaDeviceSynchronize(), "BC_Sync" ); }

//// Host|GPU transfers ////
void *BC_AllocD( int ASize ) { if (0) fprintf( stderr, cTPS" BC_AllocD(%d)\n", ASize ); void *qPtr; ikLeave( cudaMalloc(&qPtr,ASize), "BC_AllocD" ); return xBC_Err?NULL:qPtr; } // -- alloc pointers
void BC_ReleaseD( void *APtr ) { if (APtr) ikLeave( cudaFree(APtr), "BC_ReleaseD" ); } // -- dealloc pointers
void BC_Move( int ASize, void *APtrH, void *APtrD ) { ikLeave( ASize>0?cudaMemcpy(APtrD,APtrH,ASize,cudaMemcpyHostToDevice):cudaMemcpy(APtrH,APtrD,-ASize,cudaMemcpyDeviceToHost),"BC_Move" ); }

//// GPU Allocations ////
void BC_DatasAllocD( int ANDD, int AN, tBS_Val **AData ) { for (int qi=0;qi<ANDD;qi++) { AData[qi] = (tBS_Val*)BC_AllocD( AN*cBS_Sz ); if (xBC_Err) { BC_DatasReleaseD( qi, AData ); break; }; }; }
void BC_DatasReleaseD( int ANDD, tBS_Val **AData ) { for (int qi=0;qi<ANDD;qi++) { BC_ReleaseD( AData[qi] ); if (xBC_Err) break; }; }
void BC_VFil( int AL, tBS_Val *ADst, int AV ) { tBS_Val *qTmp = new tBS_Val[AL]; BS_VFil( AL, qTmp, AV ); BC_Move( AL*cBS_Sz, qTmp, ADst ); delete qTmp; }
void BC_ZFil( int AL, tBS_Val *ADst ) { cudaMemset( ADst, 0, AL*cBS_Sz ); }

//==== implicit parameters for kernel.cu =====
//============================================
extern void KTueRelease(); // DestroyTextureObject(xTue)
int  xw=0, xh=0; // standard frame-size for local operations
dim3 xG, xB; // frames' dimensions

void BC_Write( const char *AFName, int AN, const tBS_Val *ADatD ) { if (!AFName || !*AFName) return; if (!AN) AN = xw*xh;
  tBS_Val *qTmp = new tBS_Val[AN]; xBC_Err = false; BC_Move( -AN*cBS_Sz, qTmp, (void*)ADatD ); BS_Write( AFName, AN, qTmp ); delete qTmp; }

static int inline mcDivUp( int ANom, int ADen ) { int qRes = ANom/ADen; return (ANom-qRes*ADen)?qRes+1:qRes; }
void BC_Dim( int AW, int AH ) { if (xw==AW && xh==AH) return; xw = AW; xh = AH; KTueRelease(); if (AW*AH<1) return; 
  const int cBlkW = 16, cBlkH = 16; // for 2D-kernels
  xG.x = mcDivUp(xw,cBlkW); xG.y = mcDivUp(xh,cBlkH); xG.z = 1; xB.x = cBlkW; xB.y = cBlkH; xB.z = 1;
  if (0) fprintf( stderr, cTPS" BC_Dim %dx%d Grd|Blk 2D=%dx%d|%dx%d\n", xw, xh, xG.x, xG.y, xB.x, xB.y ); }

void BC_Cpy( const tBS_Val *ASrc, tBS_Val *ADst ) { if (ASrc!=ADst) ikLeave( cudaMemcpy(ADst,ASrc,xw*xh*cBS_Sz,cudaMemcpyDeviceToDevice), "BC_Cpy" ); } // -- copy frame
extern void KRowsCpy(int,const tBS_Val*,tBS_Val*); // -- copy WInp,H to W,H
void BC_RowsCpy( int AWInp, const tBS_Val *AInp, tBS_Val *AOut ) { KRowsCpy( AWInp, AInp, AOut ); ikLeave( "BC_<fromSrc>RowsCpy()" ); }
extern void KRowsCpy(const tBS_Val*,int,tBS_Val*); // -- copy W,H to WOut,H
void BC_RowsCpy( const tBS_Val *AInp, int AWOut, tBS_Val *AOut ) { KRowsCpy( AInp, AWOut, AOut ); ikLeave( "BC_RowsCpy<toDst>()" ); };

void BC_DatasWrt( const char *AFPrefix, tBS_Val **ADataD, int AI1, int AI2 ) { tBS_Val *qT[25]; memset( qT, 0, sizeof(qT) ); // -- dump frames
  int qFlags = 0xFFFFFFFF, qI1 = AI1, qI2 = AI2; if (!qI2) { qFlags = qI1; qI1 = 0; qI2 = mcASz(qT); }; const int cN = xw*xh; 
  for (int qj1=qI1;qj1<qI2;qj1++) if (ADataD[qj1] && (1<<qj1)&qFlags) { qT[qj1] = new tBS_Val[cN]; BC_Move( -cN*cBS_Sz, qT[qj1], ADataD[qj1] ); };
  BS_DatasWrt( AFPrefix, cN, qT, AI1, AI2 ); for (int qj9=qI1;qj9<qI2;qj9++) if (qT[qj9]) delete qT[qj9]; }
void BC_DatasWrt( const char *AFPrefix, tBS_Val *A0, tBS_Val *A1, tBS_Val *A2, tBS_Val *A3, tBS_Val *A4, tBS_Val *A5 ) {
  int qFlg = 0; tBS_Val *qDs[6] = {A0,A1,A2,A3,A4,A5}; for (int qi=0;qi<mcASz(qDs);qi++) if (qDs[qi]) qFlg |= 1<<qi; BC_DatasWrt( AFPrefix, qDs, qFlg, 0 ); }

//=========== down-/up-sample 2x2 ============
//============================================
extern void KDimDn(const tBS_Val*,tBS_Val*,bool); // -- dimension down 2W*2H => W*H
void BC_DimDn( const tBS_Val *AInp, tBS_Val *AOut, bool AM ) { KDimDn( AInp, AOut, AM ); ikLeave( "BC_DimDn" ); }
extern void KDimUp(const tBS_Val*,tBS_Val*); // -- dimension up W*H => 2W*2H
void BC_DimUp( const tBS_Val *AInp, tBS_Val *AOut ) { KDimUp( AInp, AOut ); ikLeave( "BC_DimUp" ); }
extern void KDimUpC(const tBS_Val*,tBS_Val*); // -- correction after 'dimension up'
void BC_DimUpC( const tBS_Val *AInp, tBS_Val *AOut ) { KDimUpC( AInp, AOut ); ikLeave( "BC_DimUpC" ); }

tBS_Val *BC_ScaleDnA( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM ) { for (int qst=1;qst<=Apyr;qst++) { // -- pyramid's dn-scale
  Aw0 >>= 1; Ah0 >>= 1; BC_Dim( Aw0, Ah0 ); BC_DimDn( ADat, ATmp, AM ); mcSwap( tBS_Val*, ADat, ATmp ); }; return ADat; } 
void BC_ScaleDn( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM ) { BC_Cpy( BC_ScaleDnA(Apyr,Aw0,Ah0,ADat,ATmp,AM), ADat ); }
tBS_Val *BC_ScaleUpA( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC ) { for (int qst=Apyr;qst>0;qst--) { // -- pyramid's up-scale 
  BC_Dim( Awp, Ahp ); BC_DimUp( ADat, ATmp ); Awp <<= 1; Ahp <<= 1;
  if (AC) { BC_Dim( Awp, Ahp ); BC_DimUpC( ATmp, ADat ); } else mcSwap( tBS_Val*, ADat, ATmp ); }; 
  if (!AC) BC_Dim( Awp, Ahp ); return ADat; }
void BC_ScaleUp( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC ) { BC_Cpy( BC_ScaleUpA(Apyr,Awp,Ahp,ADat,ATmp,AC), ADat ); }

//================ 1D arrays =================
//============================================
extern void KAbs(tBS_Val*); void BC_Abs( tBS_Val *ADat ) { KAbs( ADat ); ikLeave( "BC_Abs" ); }
extern void KNeg(tBS_Val*); void BC_Neg( tBS_Val *ADat ) { KNeg( ADat ); ikLeave( "BC_Neg" ); }
extern void KVAdd(tBS_Val*,int); void BC_VAdd( tBS_Val *ADat, int AV ) { if (AV) KVAdd( ADat, AV ); ikLeave( "BC_VAdd" ); }
extern void KVMul(tBS_Val*,float); void BC_VMul( tBS_Val *ADat, float AV ) { // -- Dat *= V 
  if (fabs(AV-1)<1E-4f) return; if (fabs(AV+1)<1E-4f) { BC_Neg( ADat ); return; };
  if (fabs(AV)<1E-4f) { BC_ZFil( xw*xh, ADat ); return; }; KVMul( ADat, AV ); ikLeave( "BC_VMul" ); }
extern void KRng(tBS_Val*,int,int); bool BC_Rng( tBS_Val *ADat, int AMn, int AMx ) {
  if ((AMn<=-0x7FFF && AMx>=0x7FFF) || (AMn>=AMx)) return false; KRng( ADat, AMn, AMx ); ikLeave( "BC_Rng" ); return true; };

extern void KRngA(tBS_Val*,int,int); bool BC_RngA( tBS_Val *ALap, int AMx, int AMxP ) { if (AMx<=0) return false; KRngA( ALap, AMx, AMxP ); ikLeave( "BC_RngA" ); return true; };
extern void KRngStd(const int*,tBS_Val*); void BC_RngStd( const int *ASrc, tBS_Val *ADst ) { KRngStd( ASrc, ADst ); ikLeave( "BC_RngStd" ); };
extern void KMul(const tBS_Val*,tBS_Val*,int); void BC_Mul( const tBS_Val *AMul, tBS_Val *ADat, int ANrm ) { KMul( AMul, ADat, ANrm ); ikLeave( "BC_Mul" ); }
extern void KMxA(const tBS_Val*,tBS_Val*); void BC_MxA( const tBS_Val *ACmp, tBS_Val *AAcc ) { KMxA(ACmp,AAcc); ikLeave("BC_MxA"); }

extern void KAdd(const tBS_Val*,tBS_Val*); void BC_Add( const tBS_Val *AAdd, tBS_Val *AAcc ) { KAdd( AAdd, AAcc ); ikLeave( "BC_Add" ); }
extern void KAddRg(const tBS_Val*,tBS_Val*); void BC_AddRg( const tBS_Val *Aadd, tBS_Val *Aacc ) { KAddRg( Aadd, Aacc ); ikLeave( "BC_AddRg" ); }
extern void KAddMul(const tBS_Val*,int*,float); void BC_AddMul( const tBS_Val*AAdd, int *AAcc, float AV ) { KAddMul( AAdd, AAcc, AV ); ikLeave( "BC_AddMul" ); }
extern void KSubA(const tBS_Val*,tBS_Val*); void BC_SubA( const tBS_Val *ASub, tBS_Val *AAcc ) { KSubA( ASub, AAcc ); ikLeave( "BC_SubA" ); }
extern void KSubS(const tBS_Val*,tBS_Val*); void BC_SubS( const tBS_Val *ASub, tBS_Val *AAcc ) { KSubS( ASub, AAcc ); ikLeave( "BC_SubS" ); }

//======== Gray-scale gamut's XLAT ==========
//===========================================
extern void KXLAT_Init(int,const int*,int,int,int,int); // -- allocate XLAT[Slot][0..Gam-1]
void BC_XLAT_Init( int ASlot, const int *AXLat, int AMn, int AGam, int ATD0, int ATD1 ) {
  KXLAT_Init( ASlot, AXLat, AMn, AGam, ATD0, ATD1 ); ikLeave( "BC_XLAT_Init" ); }
extern void KXLAT_Release(int); // -- release XLAT[Slot]
void BC_XLAT_Release( int ASlot ) { KXLAT_Release( ASlot); ikLeave( "BC_XLAT_Release" ); }
extern void KXLAT_Apply(int,tBS_Val*); // -- Dat := XLAT[Slot][Dat]
void BC_XLAT_Apply( int ASlot, tBS_Val *ADat ) { KXLAT_Apply( ASlot, ADat ); ikLeave( "BC_XLAT_Apply" ); }
extern void KGamHist(const tBS_Val*,int,unsigned*,int,const tBS_Val*); // -- Hist[Dat] '+=Ref' or '+=1'
void BC_GamHist( const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam, const tBS_Val *ARef ) {
  KGamHist( ADat, AStride, AHist, AGam, ARef ); ikLeave( "BC_GamHist" ); }
extern unsigned *KGamCLAHEHists(const tBS_Val*,int,int,const tBS_Val*); // -- get Hists[Ny][Nx]
unsigned *BC_GamCLAHEHists( const tBS_Val *ADat, int AB, int AGam, const tBS_Val *ARef) { 
  if (AB<=0) { unsigned *qHist = new unsigned[AGam]; AB *= -1; const int cw = xw, ch = xh; int qo = AB+AB*cw; // whole frame
    if (AB) BC_Dim( cw-AB*2, ch-AB*2 ); BC_GamHist( ADat+qo, cw, qHist, AGam, ARef?ARef+qo:NULL );
    if (AB) BC_Dim( cw, ch ); return xBC_Err?NULL:qHist; };
  unsigned *qHists = KGamCLAHEHists( ADat, AB, AGam, ARef ); if (!qHists) { ErrMsg( "KGamCLAHEHists params" ); return NULL; };
  ikLeave( "BC_GamCLAHEHists" ); return qHists; }
extern bool KGamCLAHEApply(tBS_Val*,int,int,int*); // -- equ to CBS_XLat.Apply
void BC_GamCLAHEApply( tBS_Val *ADat, int AB, int AGam, int *AHists ) {
  if (AB<=0) { BC_XLAT_Init( 0, AHists, 0, AGam, 0, 0 ); BC_XLAT_Apply( 0, ADat ); BC_XLAT_Release( 0 ); return; };
  if (!KGamCLAHEApply( ADat, AB, AGam, AHists )) { ErrMsg( "KGamCLAHEApply params" ); return; }; ikLeave( "BC_GamCLAHEHists" ); }

//========== Directional average  ===========
//===========================================
void BC_Ave2D( const tBS_Val *AOrg, tBS_Val *AAve, int AR, tBS_Val *ATmp ) { if (!AR) { BC_Cpy( AOrg, AAve ); return; }; // -- directional averages (adhoc-recursively)
  if (AR>16) { ErrMsg( "WRN BC_Ave2D(R>16) likely exceeds macro in kernel.cu" ); AR = 16; };
  BC_Av1( AOrg, ATmp, AR ); BC_Av5( ATmp, AAve, AR );     if (AR<3) return; int qR45 = BS_SDist( 3, AR );
  BC_Av3( AAve, ATmp, qR45 ); BC_Av7( ATmp, AAve, qR45 ); if (AR<7) return; int qR22 = BS_SDist( 2, AR );
  BC_Av2( AOrg, ATmp, qR22 ); BC_Av4( ATmp, AAve, qR22 ); BC_Av6( AAve, ATmp, qR22 ); BC_Av8( ATmp, AAve, qR22 ); }

void BC_Ave4( const tBS_Val *AOrg, tBS_Val **AAve4, int AR ) { // -- Ave4 := Ave1|3|5|7( Org )
  if (!AR) { for (int qi=0;qi<4;qi++) BC_Cpy( AOrg, AAve4[qi] ); return; }; int qR45 = BS_SDist( 3, AR );
  BC_Av1( AOrg, AAve4[0], AR ); if (xBC_Err) return; BC_Av3( AOrg, AAve4[1], qR45 ); if (xBC_Err) return;
  BC_Av5( AOrg, AAve4[2], AR ); if (xBC_Err) return; BC_Av7( AOrg, AAve4[3], qR45 ); if (xBC_Err) return; }

void BC_Ave8( const tBS_Val *AOrg, tBS_Val **AAve8, int AR ) { // -- Ave8 := Ave1|2|3|4|5|6|7|8( Org )
  if (AR<1) { for (int qi=0;qi<8;qi++) { BC_Cpy( AOrg, AAve8[qi] ); if (xBC_Err) return; }; return; };
  tBS_Val *qT[4] = {AAve8[0],AAve8[2],AAve8[4],AAve8[6]}; BC_Ave4( AOrg, qT, AR ); if (xBC_Err) return;
  if (AR<2) for (int qi=1;qi<8;qi+=2) {
    BC_Cpy( AAve8[qi-1], AAve8[qi] ); BC_Add( AAve8[(qi+1)&7], AAve8[qi] ); BC_VMul( AAve8[qi], 0.5f ); if (xBC_Err) return;
  } else { int qR22 = BS_SDist( 2, AR );
    BC_Av2( AOrg, AAve8[1], qR22 ); if (xBC_Err) return; BC_Av4( AOrg, AAve8[3], qR22 ); if (xBC_Err) return;
    BC_Av6( AOrg, AAve8[5], qR22 ); if (xBC_Err) return; BC_Av8( AOrg, AAve8[7], qR22 ); if (xBC_Err) return;
  }; }

//========== Directional |Der|  ===========
//=========================================
extern void KDerFst1(const tBS_Val*,tBS_Val*); extern void KDerFst5(const tBS_Val*,tBS_Val*); // |Der| by Dir#1,#5
void BC_DrA3x3( tBS_Val *ADat, tBS_Val **ATs ) { // -- absolute derivation (Rad=1)
  BC_ADr1( ADat, ATs[0], 1 ); if (xBC_Err) return;
  BC_ADr5( ADat, ATs[1], 1 ); if (xBC_Err) return;
  BC_Cpy( ATs[0], ADat ); BC_Add( ATs[1], ADat ); }

//======= application spatial filters =======
//===========================================
extern void KAve3x3(const tBS_Val*,tBS_Val*); void BC_Ave3x3( const tBS_Val *AInp, tBS_Val *AOut ) { KAve3x3( AInp, AOut ); ikLeave( "BC_Ave3x3" ); }
extern void KAve5x5(const tBS_Val*,tBS_Val*); void BC_Ave5x5( const tBS_Val *AInp, tBS_Val *AOut ) { KAve5x5( AInp, AOut ); ikLeave( "BC_Ave5x5" ); }
extern bool KMed(const tBS_Val*,tBS_Val*,int); void BC_Med( const tBS_Val *AInp, tBS_Val *AOut, int AR ) {
  if (KMed( AInp, AOut, AR )) ikLeave( "BC_Med" ); else ErrMsg( "KMed params" ); }

//======= Radon-Shannon-denoise =======
//=====================================
extern void KNoise(int,int,tBS_Val**,tBS_Val*); void BC_Noise( // -- basic RadonShannon_noise
  const tBS_Val *AOrg, tBS_Val *ANoise, tBS_Val **ATs8, int AR, int AThrLg ) {
  BC_Ave8( AOrg, ATs8, AR ); if (xBC_Err) return; for (int qi0=0;qi0<8;qi0++) BC_SubA( AOrg, ATs8[qi0] ); 
  // this might balance StdDev's differene (not the morphology's!) against CPU 
  // for (int qi1=0;qi1<8;qi1++) BC_VMul( ATs8[qi1], 0.97f );
  KNoise( AThrLg%100, AThrLg/100, ATs8, ANoise ); ikLeave( "BC_Noise" ); }

int BC_MxLapStd( // -- reduce shot_noise after RadonShannon_noise
  int AR, tBS_Val *ALap, tBS_Val **ATs, float AK, int AZ0, int AZ1, float AShear ) {
  extern void KMxLapRngZ(int,tBS_Val*); extern void KMxLapStd(int,tBS_Val**,tBS_Val*);
  // call-sequence: "_MxLapRngZ(Z0); MaxLap; _MxLapRngZ(Z1)"
  if (AR<1) return 0; KMxLapRngZ( AZ0, ALap ); 
  const float cEps = 1e-3f; const bool cJoin = (AK>cEps)&&(AK<1-cEps);
  if (cJoin) BC_Cpy( ALap, ATs[9] ); int cNDD; switch (AR) {
    case  1: { BC_Ave4( ALap, ATs, AR ); if (xBC_Err) return 0;
      KMxLapStd( 4, ATs, ALap ); KMxLapRngZ( BS_MxLapZ1dflt(AZ0,AZ1), ALap );
      if (AK<cEps) { BC_Ave3x3( ALap, ATs[8] ); BC_Add( ATs[8], ALap ); BC_VMul( ALap, 0.5f ); };
      cNDD = 4; }; break;
    default: {
      BS2A_ShearCUDA( 0, AShear, ALap, ATs[0] );
      BC_Ave8( ALap, ATs, AR ); if (xBC_Err) return 0;
      KMxLapStd( 8, ATs, ALap ); KMxLapRngZ( BS_MxLapZ1dflt(AZ0,AZ1), ALap );
      BS2A_ShearCUDA( 0, -AShear, ALap, ATs[0] );
      if (AK<cEps) { BC_Ave3x3( ALap, ATs[8] ); BC_Cpy( ATs[8], ALap ); };
      cNDD = 8; }; }; // switch R
  if (AK<cEps) return cNDD; BC_VMul( ALap, AK );
  if (cJoin) { BC_VMul( ATs[9], 1-AK ); BC_Add( ATs[9], ALap ); };
  return cNDD; } // returns #NDD of dir averages in TS^

//=============== Reductions ================
//===========================================
extern double KMean(const tBS_Val*,tBS_Val*,tBS_Val*); double BC_Mean( const tBS_Val *AInp, tBS_Val *ATmp, tBS_Val *ATmpH ) {
  double qRes = KMean( AInp, ATmp, ATmpH ); if (qRes>-1e36) return qRes; ikLeave( "BC_Mean"); return 0; }
