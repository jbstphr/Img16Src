// I16_Ut6.h : CPU|CUDA ad-hoc utilities (with explicit dimension)
// This header is for ~"despeckle", when CUDA "resembles CPU", but doesn't match
// exactly. Sample situation is:
//  a) Suppress defecive lines before unsharp by strong anisotropic CPU-despeckle
// within <logically optional> ~virtual horz-bands .
//  b) Within unsharp, call commonplace 2D isotropic despeckle (either CPU or
// CUDA) as ~"cartesian 2D-median-filter with threshold<s>".
//  c) We would prefer C-code of both procedures in the same file (as CPU-code
//for "a" and "b" is quite similar (though params of _xyzCPU(...) differ
//against _xyzCUDA(...). 

#ifndef I16_UT6_H_INCLUDED
#define I16_UT6_H_INCLUDED

//// despeckle (for certain bands) ////
// CPU despeckle (~same code for general_denoise & horz_defects_correction).
// This approach compares a pixel against both 'median' and 'mean'.
void BS6_DespeckleScalars( int AThr, int AStrong ); // despeckle scalar params
  // Strong=false:for ~general_denoise
  //       =true :for ~horz_defects_correction
  // general approach accounts median and isotropic_average with ~similar importance
  // strong approach accounts less median but more anisotropic average (by vert_neighbors,ignores horz_neighbors)
  
void BS6_Despeckle( int AW, int AH, int AH0, int AH1, // despeckle within a band
  tBS_Val *ADat, int ARadMed, tBS_Val *AMed, int ARadAve, tBS_Val *AAve ); 
  // calculate median & average within band [0..H[
  // despeckle within band [H0..H1[

// CUDA general despeckle (as 'Strong=false' for CPU ).
void BS6_DespeckleCUDA( int AW, int AH, tBS_Val *ADat, // -- general despeckle
  int ARadMed, tBS_Val *AMed, int ARadAve, tBS_Val *AAve );

//// statistical suppress of thin & long cartesian laplacians ////
// This approach subtracts thin & long laplacians (optionally divided on a mask, obtained from 'cut
// |Dr3x3A|'). TLLs incapsulate either "objects' edges" and "slow gradients", which are ~proportional
// to the module of derivation <with ~similar coefficients in a practice>.
// Typical TLL:+++++   After       :ooooo   After    :o++oo   where   :29920   If TLL divides on an
// in presence:-----   subtract    :ooooo   subtract :oo--o   |Der3x3|:01881   <empirically normalized>
// of a <tilt>:+++++   TLL by lucky:ooooo   real TTL :o+ooo   could   :07000   mask by |Der3x3|, then we 
// antiscatter:-----   RadX, RadY, :ooooo   of tilted:o-ooo   look ~as:06000   could avoid excessive
// grid       :+++++   Attenuation :ooooo   grid     :oo+oo   this    :00500   removal of object's edges.
typedef struct { int n[3]; float f[1]; } tThLoLap; const int
  cBS6_tllMode = 0, // '=1':set optional normalized mask, '=2':'execute'
  cBS6_tllCGN0 = 1, // gamut's cut0 for optional normalized mask
  cBS6_tllCGN1 = 2, //         -"-1 
  cBS6_tllRadX = 1, // 
  cBS6_tllRadY = 2, // 
  cBS6_tlfAtt  = 0; // attenuation
  // Mode=1 
  // Mode=2 "plain suppression"; RadX>RadY:supress blemish rows, RadY>RadYXsupress blemish cols
void BS6_ThLoLapReset();
void BS6_ThLoLapSuppress( const tThLoLap *ATLL, // suppress thin & long laplacians
  int AW, int AH, tBS_Val *ADat, tBS_Val **ATsH, tBS_Val **ATsD, const char *ADbgPfx );
  // changes Dat; records the change in Ts#0

//// Combined Poisson's denoise ////
void BS6_DePoisson( int AW, int AH, const tUnSharp *AUS, tUnSharpExts *AUSExts,
  tBS_Val *AFrm, tBS_Val **ATsH, tBS_Val **ATsD, const char *ADbgPfx ); // Poisson_noise removal

//// Radon's despeckle //// 
//void BS6_RadonDespeckle( int AW, int AH, tBS_Val *ADat, tBS_Val *AAve, tBS_Val **ATs,
//  int AR, int AFlg, float ALimK = 1.25f ); // Radon's despeckle CPU
//void BS6_RadonDespeckle_CUDA( int AW, int AH, tBS_Val *ADat, tBS_Val *AAve, tBS_Val **ATs,
//  int AR, int AFlg, float ALimK = 1.25f ); // Radon's despeckle CUDA


//// Unsharp's upscale //// 
// emprical blur/sharp is easier than ~'bicubic spline' for correctnes of edges
//int BS6_ScaleUpBeforeUnsharp( int APyr, int AW, int AH, tBS_Val *ADat, tBS_Val **ADs ); // upscale before unsharp
  // upscale might keep #NBuUps few bufs WxH at Ds+0 (caller shall increment Ds)
//int BS6_ScaleDnAfterUnsharp( int APyr, int AW, int AH, tBS_Val *ADat, tBS_Val **ADs ); // downscale after unsharp
  // downscale uses upscale's bufs WxH at Ds-#NBuUps; it returns #bufs at Ds+0 to print (if debug trace is enabled)

//void BS6_QTest( int AW, int AH, tBS_Val *A0, tBS_Val *A1=NULL, tBS_Val *A2=NULL );
#endif // I16_UT6_H_INCLUDED