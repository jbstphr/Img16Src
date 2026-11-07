// I16_Ut1.h : CUDA interconnect utilities
// I16_Ut1.h : CUDA interconnect utilities
// the routines here mostly operate with fixed frames BC_Dim()

#ifndef I16_UT1_H_INCLUDED
#define I16_UT1_H_INCLUDED

//================= commons ==================
//============================================
extern bool xBC_Err; void BC_DeviceFind();
void BC_DeviceReset();
void BC_Sync();

//// Host|GPU transfers ////
void *BC_AllocD( int ASize );
void BC_ReleaseD( void *APtr );
void BC_Move( int ASize, void *APtrH, void *APtrD ); // size in bytes >0:HostToDev, <0:DevToHost
void BC_Write( const char *AFName, int AN, const tBS_Val *ADat ); // dump frame

//// GPU Allocations ////
void BC_DatasAllocD( int ANDD, int AN, tBS_Val **ADataD ); // alloc pointers
void BC_DatasReleaseD( int ANDD, tBS_Val **ADataD ); // dealloc pointers
void BC_ZFil( int AL, tBS_Val *ADst );
void BC_VFil( int AL, tBS_Val *ADst, int AV );

//==== implicit parameters for kernel.cu =====
//============================================
void BC_Dim( int AW, int AH ); // Frame's dimension in DSeg: sets xw,xh & assigns xGrd,xBlk
// -- below operates with standard dimension & perhaps other implicit parameters for kernel.cu -- //
// -- procedures with 'Stride' need '_Dim( Winner, Hinner )' -- //

void BC_Cpy( const tBS_Val *ASrc, tBS_Val *ADst ); // copy frame
void BC_RowsCpy( int AWInp, const tBS_Val *AInp, tBS_Val *AOut ); // copy WInp,H to W,H
void BC_RowsCpy( const tBS_Val *AInp, int AWOut, tBS_Val *AOut ); // copy W,H to WOut,H

void BC_DatasWrt( const char *AFPrefix, tBS_Val **ADataD, int AI1=0, int AI2=8 ); // dump frames
void BC_DatasWrt( const char *AFPrefix, tBS_Val *A0, tBS_Val *A1=NULL, tBS_Val *A2=NULL, tBS_Val *A3=NULL, tBS_Val *A4=NULL, tBS_Val *A5=NULL );

//=========== down-/up-sample 2x2 ============
//============================================
void BC_DimDn( const tBS_Val *AInp, tBS_Val *AOut, bool AMA=false ); // dimension down 2W*2H => W*H
void BC_DimUp( const tBS_Val *AInp, tBS_Val *AOut ); // dimension up W*H => 2W*2H
void BC_DimUpC( const tBS_Val *AInp, tBS_Val *AOut ); // correction after 'dimension up'
// pyramid's dn- & up-scale; M:true - max amplitude (default average); C:false - no _DimUpC (default _DimUpC)
tBS_Val *BC_ScaleDnA( int Apyr, int Aw0, int Ah0, tBS_Val *ADat, tBS_Val *ATmp, bool AM=false ); void BC_ScaleDn(int,int,int,tBS_Val*,tBS_Val*,bool=false);
tBS_Val* BC_ScaleUpA( int Apyr, int Awp, int Ahp, tBS_Val *ADat, tBS_Val *ATmp, bool AC=true ); void BC_ScaleUp(int,int,int,tBS_Val*,tBS_Val*,bool=true);

//================ 1D arrays =================
//============================================
void BC_Abs( tBS_Val *ADat ); // Dat:=|Dat|
void BC_Neg( tBS_Val *ADat ); // Dat:=-Dat
void BC_VAdd( tBS_Val *ADat, int AV ); // Dat+=V
void BC_VMul( tBS_Val *ADat, float AV ); // Dat*=V
bool BC_Rng( tBS_Val *ADat, int AMn, int AMx=0x7FFF ); // Dat:=range(Dat)
bool BC_RngA( tBS_Val *ADat, int AMx, int AMxP=0x7FFF ); // Mx<V<Mx" or "-Mx<V<MxP"
void BC_RngStd( const int *ASrc, tBS_Val *ADst ); // Dst:=standard_range(Src)
void BC_Mul( const tBS_Val *AMul, tBS_Val *ADat, int ANrm );
void BC_MxA( const tBS_Val *ACmp, tBS_Val *AAcc );

void BC_Add( const tBS_Val *AAdd, tBS_Val *AAcc ); // Acc += Add
void BC_AddRg(const tBS_Val*,tBS_Val*); // Acc += Add
void BC_AddMul( const tBS_Val *AAdd, int *AAcc, float AV ); // Acc += V*Add
void BC_SubA( const tBS_Val *ASub, tBS_Val *AAcc ); // Acc:=Sub-Acc
void BC_SubS( const tBS_Val *ASub, tBS_Val *AAcc ); // Acc:=Acc-Sub

//======== Gray-scale gamut's XLAT ==========
//===========================================
void BC_XLAT_Init( int ASlot, const int *AXLat, int AMn, int AGam, int ATD0=0, int ATD1=0 ); // allocate XLAT[Slot][0..Gam-1]
void BC_XLAT_Release( int ASlot ); // release XLAT[Slot]
void BC_XLAT_Apply( int ASlot, tBS_Val *ADat ); // Dat := XLAT[Slot][Dat]
void BC_GamHist( const tBS_Val *ADat, int AStride, unsigned *AHist, int AGam, const tBS_Val *ARef=NULL ); // Hist[Dat] '+=Ref' or '+=1'
unsigned *BC_GamCLAHEHists( const tBS_Val *ADat, int AB, int AGam, const tBS_Val *ARef ); // get Hists[Ny][Nx]
void BC_GamCLAHEApply( tBS_Val *ADat, int AB, int AGam, int *AHists ); // equ to CBS_XLat.Apply
  
//========== Directional average  ===========
//===========================================
// 'conditionally transversal' to the standard limbo in Ut0.cpp
// const tBS_Val *AInp, tBS_Val *AOut, int AR
void BC_Av1( const tBS_Val*,tBS_Val*,int); void BC_Av3( const tBS_Val*,tBS_Val*,int); // '-' #1 '\' #3
void BC_Av5( const tBS_Val*,tBS_Val*,int); void BC_Av7( const tBS_Val*,tBS_Val*,int); // '|' #5 '/' #7
void BC_Av2( const tBS_Val*,tBS_Val*,int); void BC_Av4( const tBS_Val*,tBS_Val*,int); // #1..#3 #3..#5
void BC_Av6( const tBS_Val*,tBS_Val*,int); void BC_Av8( const tBS_Val*,tBS_Val*,int); // #5..#7 #7..#1
void BC_Ave2D( const tBS_Val *AOrg, tBS_Val *AAve, int AR, tBS_Val *ATmp ); // directional averages (adhoc-recursively)
void BC_Ave4( const tBS_Val *AOrg, tBS_Val **AAve4, int AR ); // Ave4 := Ave1|3|5|7( Org )
void BC_Ave8( const tBS_Val *AOrg, tBS_Val **AAve8, int AR ); // Ave8 := Ave1|2|3|4|5|6|7|8( Org )
  // "R<-2" returns skewed averages by |R|

//========== Directional |Der|  ===========
//=========================================
// const tBS_Val *AInp, tBS_Val *AOut, int AR
void BC_ADr1( const tBS_Val*,tBS_Val*,int); void BC_ADr3( const tBS_Val*,tBS_Val*,int); // '-' #1 '\' #3
void BC_ADr5( const tBS_Val*,tBS_Val*,int); void BC_ADr7( const tBS_Val*,tBS_Val*,int); // '|' #5 '/' #7
void BC_ADr2( const tBS_Val*,tBS_Val*,int); void BC_ADr4( const tBS_Val*,tBS_Val*,int); // #1..#3 #3..#5
void BC_ADr6( const tBS_Val*,tBS_Val*,int); void BC_ADr8( const tBS_Val*,tBS_Val*,int); // #5..#7 #7..#1
void BC_DrA3x3( tBS_Val *ADat, tBS_Val **ATs ); // absolute derivation (Rad=1)

//======= application spatial filters =======
//===========================================
void BC_Ave3x3( const tBS_Val *AInp, tBS_Val *AOut ); void BC_Ave5x5( const tBS_Val*, tBS_Val* ); 
void BC_Med( const tBS_Val *AInp, tBS_Val *AOut, int AR ); // the medians might be not exact ones for large radiuses

//======= Radon-Shannon-denoise =======
//=====================================
void BC_Noise( // basic RadonShannon_noise
  const tBS_Val *AOrg, tBS_Val *ANoise, tBS_Val **ATs8, int AR, int AThrLg ); 
  // inp    Org: ~laplacian (around zero)
  //      ThrLg: regul_add*100 or threshold_log2
  // out  Noise: estimated noise 
  //        Ts8: directional averages (for adhoc usage)

int BC_MxLapStd( // reduce shot_noise after RadonShannon_noise; returns #NDD of dir averages in TS^
  int AR, tBS_Val *ALap, tBS_Val **ATs, float AK, int AZ0=0, int AZ1=-1, float AShear=-0.033f ); 
  // call-sequence: "_MxLapRngZ(Z0); MaxLap; _MxLapRngZ(Z1)";
  // always: "Ave<MxLapStd<MxLapExact" or "MxLapExact<MxLapStd<Ave";
  // default Shear<>0 <might> slightly reduces aliases (it should be Ok for the "laplacians").

//=============== Reductions ================
//===========================================
double BC_Mean( const tBS_Val *AInp, tBS_Val *ATmp, tBS_Val *ATmpH );

#endif // I16_UT1_H_INCLUDED
