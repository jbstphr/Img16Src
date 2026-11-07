// I16_Ut7c.h : ribs-detection's services
#ifndef I16_UT7C_H_INCLUDED
#define I16_UT7C_H_INCLUDED

//======= Polyfit and associated funcs =======
//============================================
int BS7c_Polyfit1( int AN, const tBS_Val *AX, const tBS_Val *AY, float *AA ); // returns A0,A1 so that "AYi ~= AXi*A1+A0"
int BS7c_AscendRowPf1( int AN, tBS_Val *AX, tBS_Val *AY, float *AA, int AGap=8 ); // exclude XYpoints, which decrease A1
// Polyfit associates
#define BS7_ASCENDDIRLESS
#ifdef BS7_ASCENDDIRLESS
  const int cBS7_AscendDirs[] = {7,8,1,2,3,4};
#else
  const int cBS7_AscendDirs[] = {6,7,8,1,2,3,4};
#endif
void BS7c_AscArcRowLimbo( int AN, int ACX, tBS_Val *ADat, float A0, float A1 ); // fill row by limbo_metric (by A0,A1; center X at 'CX')
int BS7c_AscArcRowSkew( int AN, int ACX, tBS_Val *ADat, float *ATmp=NULL ); // row's skew by Y
void BS7c_ArcToMet( const tBS_Val *AMsk, tBS_Val **AArc, tBS_Val *AMet, int AW, int AH ); // prepare inp of AscendFrmPf

typedef struct { float rA[2]; } tBS7_AscArcFrmPf;
void BS7c_AscArcFrmPf( tBS_Val *ADat, tBS7_AscArcFrmPf *AAFP,
  const char *ADbgPfx, int AW, int AH ); // ascending polyfit of frame rows
  // inp: raw limbo_metrics (which should be ascending); i.e. indexes ~of arc's pies;
  // out: the fileds of R^; contents of Dat^ become ~meaningless
void BS7c_PfSkewFrm( bool AHrzToArc, tBS_Val *ADat, float A0, float A1, int AFill,
  const char *ADbgPfx, int AW, int AH ); // polyfit_skew of a frame by A0,A1
  // Convex=false:"arc to ~hrz_line or hrz_line to convex", =true:"convex to ~hrz_line or hrz_line to arc"
  // A0,A1 specify against the "_BS9_AscendDirs"

//=== Apps of Descend periphery to center ====
//============================================
void BS7c_FastDescend( const tBS_Val *AOrgConvex, int ARadBlr, tBS_Val **ARes,
  int AW, int AH ); // plain fast_descend
  // unlike practical slow_descend it applies 8 directions; it is more sharp than "clavicle's fast_descend"
void BS7c_DescendXcs( int ANDD, const int *ADirs, const tBS_Val *AOrgConvex, tBS_Val **ARes, tBS_Val **ATs,
  int AW, int AH ); // excess of slow_descend (i.e. ~"Org-slow_descend")
  // This procedure doesn't utilize linear shear/warp. NDD=7 is valid rather on ~top|clavicle and on
  // the right (where ribs might connect outer ribcage either ~horizontally or ~vertically). NDD=5 or
  // =3 is valid rather where ~vertical connectivity with outer ribcage denotes ~vessels:
  //  - on the left (the ribs connect outer ribcage ~horizontally there);
  //  - near ~bottom|heart (vessels connect ~vertically to bottom's framing shadow).
  // Out#0:xcs_accurate[0..NDD-1], #1:xcs_rough[1..NDD-2]; Ts:NDD descends.
void BS7c_PfSkewDescendXcs( const tBS_Val *AOrgConvex, float A0, float A1, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ); // excess of skewed slow_descend
  // Invoke _PfSkewFrm and then _DescendXcs(NDD~=3). 
  // This procedure doesn't utilize linear shear/warp. 

//=== Treatment of maxmag & associates ====
//=========================================
void BS7c_PfSkewMaxMag( const tBS_Val *ALap, float A0, float A1, bool ATheFinal, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ); // estimate horz_blur(main_ribcage_lap)
  // TheFinal=true utilizes linear shear/warp internally (alike extended_mask_recount in I16_Ut5d.h).
  // inp: main_ribcage_lap (aliased). This means posterior ribs with up-arc "/^\".
  // out: main_ribcage_lap (blurred along pf-curvature).
  // Clavicle_remains: Org-main_ribcage_lap (clavicle and anterior ribs with down-arc "\_/").

void BS7c_RibsRefine( const tBS_Val *ARibs, const tBS_Val *AMgnCon, const tBS_Val *ACrd,
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ); // "Ribs -= MgnCon"
  // Res#0:ribs_rfn 1:mgncon_rfn (remove margin connections from primay ribs)
  // invoke MMRWr.RibsRefine(...); 

void BS7c_MaxMagRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag,
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ); // empirically refine MaxMag
  // Res#0:maxmag_rfn #1:maxmag_sms
  // invoke MMRWr.MaxMagRefine(...); 

void BS7c_CombineRefine( const tBS_Val *ARibs, const tBS_Val *AMaxMag,
  tBS_Val **ARes, const char *ADbgPfx, int AW, int AH ); // combine refines
  // Res#0:maxmag_rfn
  // invoke MMRWr.CombineRefine(...); 

//=== Clavicle's estimation ====
//==============================
void BS7c_DirDescends( const tBS_Val *AConvex, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ); // plain directional descends.
  // Ts<>nil utilizes linear shear/warp internally.
  // inp: convex of Clavicle_remains.
  // out: 8 directional_descends of Clavicle_remains.
void BS7c_Clavicle( tBS_Val **ADes, tBS_Val **ARes, tBS_Val **ATs,
  const char *ADbgPfx, int AW, int AH ); // clavicle's estimation
  // inp: 8 directional_descends of Clavicle_remains.
  // out #0:fast_descend, #1:cross_directional_excess.
  // The cross_directional_excess:
  // - nearly estimates ribs on the ~top|clavicle;
  // - estimates ~ribs on the ~right;
  // - nearly estimates vessels on the ~left.

//======= narrow spatial filters ======
//=====================================
void BS7c_MMAmbiToShrp( const tBS_Val *AShrp, const tBS_Val *AAmbi, const tBS_Val *ASMSdet,
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ); // approach ambient shape to exact sharp one
  // inp:exact "sharp shape", "ambient shape", auxilliary split mask; out#0:ambi_shrp
  // Unlike RRAmbiToShrp the target sharp shape here is "known exactly". Hence this apporach is ~simple.

void BS7c_RRAmbiToShrp( int AIC, const tBS_Val *AOrg, const tBS_Val *ARRib, const tBS_Val *ASplitAux,
  tBS_Val **ARes, tBS_Val **ATs, const char *ADbgPfx, int AW, int AH ); // approach ambient shape to <approx known> sharp one
  // inp:Org, removed_ribs, aux_splt (migth be nil)
  // out#0:removed_ribs_shrp, (when SplitAux=nil) #1:ribs_splt
  // Unlike MMAmbiToShrp the target sharp shape here isn't "known exactly". Hence this approach needs
  // complex directional filters.
  // "InvocationCase=0:LRLung_RRAmbiToShrp", "=1,2:ProcessLftRgtFROI"

//======= helpers for ProcessLftRgtLung ======
//============================================

void BS7c_CrdRgtTopToLftBtm( tBS_Val *ADat,
  int AW, int AH ); // empirical ascend RgtToLft within a lung_area

void BS7c_CrdOfFastDescend( float ACutDrk, float ACutBrt, bool AInv, tBS_Val *ACrd,
  int AW, int AH ); // normalize convex by histogram

void BS7c_MgnConSimple( const tBS_Val *AOrg, bool ASkews, tBS_Val **ARes,
  tBS_Val **ATs, int AW, int AH ); // simple marginal connectivities
  // "Skews=1" gives more accuracy & periheral aliases, which caller should remove
  // later by ~"Mul(lung_area_msk,Res[0])".

//================ in development ============
//============================================

/*
//void BS7_CrdA( bool AVrt, tBS_Val *ADat, int AW=0, int AH=0 ); // coordinates of cardinal metric

//void BS7_NearShape( const tBS_Val *AOrg, const tBS_Val *ADat, int ARad, int ANDD, const int *ADirs,
//  tBS_Val *ARes, tBS_Val **ATs, int AW=0, int AH=0 ); // -- nearest to Org dir_ave(Dat) 

void BS7_RadonRibCageLaps( const tBS_Val *AOrg, tBS_Val **AAv8, tBS_Val **ALaps, int AW=0, int AH=0 ); // ribcage laps
void BS7_MgnCon( const tBS_Val *AOrg, tBS_Val **ALap8, tBS_Val **ADev8, int AW=0, int AH=0 ); // --
void BS9_MatchAreas( const tBS_Val *ASrc, tBS_Val **ARess, const int *AWH=NULL );

void BS9_Trinaries( int ANDD, tBS_Val **ALaps, const int *AWH=NULL );
void BS9_Trinaries( int ARad, tBS_Val *ADat, const int *AWH=NULL );
*/
#endif // I16_UT7C_H_INCLUDED
