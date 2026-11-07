// I16_Ut5.h : CPU|CUDA unsharpen's core
// frequent change needs 'synchronization|equivalence' between CPU & CUDA

#ifndef I16_UT5_H_INCLUDED
#define I16_UT5_H_INCLUDED
// unsharp
typedef struct {
  bool  rEnabled;
  int   rDeG[1]; // des_Thr
  float fDeG[4]; // des_Gamma des_Kgr des_KnP des_KSnP
  int   rDeR[1]; // des_RadR
  float fDeR[3]; // des_KR des_KU des_Kmf
  float fDeD[4]; // des_KdSS "A+a* A+l* L+a* L+l*"
  int   rDeN[4]; // des_RadN des_LgN des_RadN2 des_ZN2
  float fDeN[2]; // des_KN des_KN2
  int   rDeS[3]; // des_RadSnP des_r0SnP des_r1SnP
  int   rMgn[3]; // mgn_Msk mgn_Bev mgn_Cor
  int   rA[6]; // a_RadHF a_RadMM a_TGMM a_ZMM a_Lim a_LimP
  float fA[5]; // a_KHF a_KMM a_KAO a_KCPU a_DnSpDi
  float fB[3]; // b_GammaSpa b_GammaCLA b_GammaAtt
  int   rC[3]; // c_RadBas c_RadBig c_RadMM
  float fC[3]; // c_RefBig c_DelBig c_KMM
  int   rD[2]; // d_Blk d_Blr
  float fD[3]; // d_BodyOne d_ExtMn d_ExtMx
  int   rE[4]; // e_Blk e_Blr e_Lim e_Zero
  float fF[4]; // f_Att f_EmLthn f_EmLedg f_MxLbias
  float fFK[3][3]; // f_K
  int   rFC[4][5]; // f_C
  int   rH[2]; // h_Lim h_LimP
  int   rJ[4]; // j_RadMM j_RadMMO j_RadENm j_Lim
  float fJ[2]; // j_KMM j_KENm
  int   rL[3]; // l_RadHess l_RadMM l_Lim
  float fL[3]; // l_KEdg l_KRdg l_KMM
  float fN[3]; // n_Z
  float fSA[4]; // s_Za
  float fSB[3]; // s_Zb
  float fU[3]; // u_BlOuDG u_Air0Mn u_Air1Mn
} tUnSharp; // control record
  // des_Thr   ~11; despeckle's Threshold
  // des_Gamma  ~1; balance_gamma before despeckle
  // des_KG   ~1.0; Org += _KG*SrapeGrid(Org)
  // des_KnP  ~0.6; Org += _KnP*ScrapenP(Org) <apply at ~denoise_begin>
  // des_KSnP ~0.2; <apply at ~denoise_end> 
  // des_RadR  ~0; rad of Radon_Radial_Despeckle 
  // des_KR  ~1.2; "RadonDespeckle(...,KR)"
  // des_KU  ~3.0; _MxLapXLimbo(-KU,...)
  // des_Kmf ~1.0; Denoise_Res:=Ave+Kmf*Mid_Freq
  // des_KdSS ~0.5 0.0 3.0 -1.0; scrap_std dispatch "A+a* A+l* L+a* L+l*"
  // des_RadN  ~4; rad of Radon_Shannon_denoise
  // des_RadN2 ~9; rad2 of morpho_residuals
  // des_RadSnP ~203; large blur mask to determine Scrape_NonProcessed
  // des_r0SnP,r1SnP ~111,-3456; frame+=KSnP*SnP*maskSnP(frame)
  // des_LgN   ~5; log2 of shannon_noise's threshold
  // des_ZN2   ~0; MxLapRngZ
  // des_KN  ~0.5; "-=KN*shannon_noise"
  // des_KN2 ~0.5; shannon_noise-=KN2*morpho_residuals
  // mgn_Msk   ~0; >0(~1) forces org_zeros to be also final_zeros
  //               <0: do nothing
  // mgn_Bev ~0; size of corner_bevel by X,Y
  // mgn_Cor ~0; location of corner_bevels
  // a_RadHF ~1; HiFblr:=blr(RadHF,Inp), HiFrem:=Inp-HiFblr
  // a_RadMM ~0:do nothing; _MxLapXLimbo(Rad,TolerangeGap=TGMM);
  // a_TGMM  -1:no_tolerance_gap
  // a_ZMM   ~0; _MxLapStd(Z0=ZMM,Z1=-1)
  // a_Lim =0 don't truncate HiFrem; ~999
  // a_KHF ~0.5; operate with 'HiFblr+HiFrem*KHF' (~>1.5 would cause ~non-robust CLAHE)
  // a_KMM ~0.5; Lap+=KMM*MaxMag(Lap,RadMM)
  // a_KAO ~6.0; Auto_Oscillator; HiF ~= func(HiF,|MaxLap|)
  // a_KCPU ~0.0; multiply on |K|; K<0:when>_in_CPU; K>0:when>_in_CUDA,
  // a_DnSpDi ~0.75; spatial dilation if downsample activates
  // b_GammaSpa ~1.0; Gamma for spatial metrics
  // b_GammaCLA ~3.0; Gamma for CLAHE
  // b_GammaAtt ~1.0; Gamma for Attenutaion mask
  // c_RadBas,c_RadBig - radiuses of LapBase,LapBig
  // c_RadMM ~3; MaxMag(LapBas,RadMM); 0:do nothing
  // c_RefBig - reference average of |LapBig|; 0-don't use
  //   it might change <by fixed formulas> ~e_Blr, etc
  // c_DelBig ~0.3; 'delta=1' for c_RefBig
  // c_KMM  ~0.5; LapBas+=KMM*MaxMag(LapBas,RadMM)
  // d_Blk =5..9
  // d_Blr =0 no blur; ~22
  // d_BodyOne =0 don't use; ~0.1 absolute min of a tissue (~max background)
  // d_ExtMn,Mx ~0.5; SLikeCurve gamut's extensions
  // e_Blk =5..9
  // e_Blr =0 no blur; ~7
  // e_Lim  ~123 (against e_Blk=6); if "count>e_Lim" then "mskAtt=0"
  // e_Zero ~256 (against MskAve=512) decrement before multiplications
  // f_Att     ~5.0; MskGDi += MskAtt*Att
  // f_EmLthn  ~0.10; lapAcc+=EmLthn*Mimic(LapBas,mskThn)
  // f_EmLedg  ~0.05; lapAcc+=EmLedg*Mimic(LapBig,mskEdg)
  // f_MxLbias ~1.15; _MxLapStd(~mskGDi or mskEdg,_MxLbias)
  // f_K - denormalized_masks; |Der3x3|*,|Der3x3_CLAHE|*,|LapBig|*
  //    a ~GSh; global shape of object (calcs MskThk,MskEdg)
  //    b ~GDi; global dilation of high freqs
  //    c ~GThn; derivations mimic laplacians <edges>
  // f_C - gamut_cuts of normalized_masks; rg0 rg1 tg0 tg1 rad
  //    a ~mskThn; ~'Thin'; merge LapAda,LapMM to LapAdaEff
  //    b ~mskThk; ~'Thick'; merge LapBas,LapBig to LapAda
  //    c ~mskGDi; ~'general dilation'; dilate HiFrem, LapAdaMM
  //    d ~mskEdg; ~'edges for MimicLap()'
  // h_Lim =0 don't truncate LapAda; ~999
  // j_RadMM ~8; RadonMinMax's radus 
  // j_RadMMO ~0; OrthoLap's radius; >0 reduces auras (for ~NDT)
  // j_RadENm ~1; EdgeNaiveMimic's radius
  // j_Lim =0 don't truncate LapAdaEff; ~99
  // j_KMM ~1.25; ponder LapMM*KMM against LapAda by mskThn
  // j_KENm ~0.5; lapSrp += lapENm*j_ENm (EdgeNaive's mimic Laplacian)
  // l_RadHess ~3; Hessian's radius
  // l_RadMM ~3; MaxMag(LapAcc,RadMM); 0:do nothing
  // l_Lim =0 don't truncate LapAcc; ~999
  // l_KEdg ~0.5; += edge*l_KEdg
  // l_KRdg ~0.1; += ridge*l_KRdg
  // l_KMM  ~0.5; LapAcc+=KMM*MaxMag(LapAcc,RadMM)
  // n_Z ~0 0.4 0.2; mskGDiRes := (1-n_Z0-n_Z1-n_Z2) +
  //       n_Z0*mskThk +
  //       n_Z1*mskGDi +
  //       n_Z2*mskEdg.
  // s_Za ~0.5 0.3 9.9 0.0; BodySum :=
  //       (1-s_Za0-s_Za1)*HiFblr + -- Ab2.raw
  //       s_Za0*BlrAda +           -- Hb1.raw
  //       s_Za1*GammaCLA +         -- B1.raw
  //       s_Za2*HiFrem +           -- Ab0.raw
  //       s_Za3*LapMid.            -- Ha0.raw
  // s_Zb ~17.0 22.0, 0.0; LapSum := 1.0*LapAcc +
  //       s_Zb0*CLAInc +
  //       s_Zb1*LapAda +
  //       s_Zb2*LapAdaEff.
  // u_BlOuDG ~0.2l; BlackOutDnGray
  // u_Air0Mn u_Air1Mn ~3e-3 1e-5; DetectedAir0Mn DetectedAir1Mn

typedef struct { int rBCSi; bool rIsFloat; float rKmn,rKmx; int rAddrOffs; } tUnSharpAEn;
  // adjustment's entity (adjustment's constraints for field '@UnSharp+AddrOffs')
  // BCSi: '0':change this field by tag B, '1':by tag C, '2':by tag S
  // Kmn,Kmx: multiplers of default value (at 'tag=8') for min (at 'tag=1'), max (at 'tag=15')
  // AddrOffs: field's offset from begin of the record (in bytes)
typedef struct { int rN; tUnSharpAEn rE[32]; } tUnSharpAdjs;
  // collection of adjustment's entities; created by ~"App_Init():UnSharpLoad()";
  // used by "_UnSharpRecAdjust"

typedef struct {
  // extended (transient, ad-hoc) data from UnSharpRecAdjust to _UnSharp_
  int   rAdj,          // keeps original TBCS
        rModeId;       // adhoc copy of global 'xModeId'
  // _UnSharp_'s reports about run-time facts
  int   rAirEn,        // BS5a_AirEstimationEnabled
        rQuadEn,       // BS5c_QuadLapEnabled
        rMskEn;        // BS5d_MskEnabled
  // run-time stats, changeable params inside _UnSharp_
  //int   rSttFlg;       // work-flow 'state flags'
  //                     // 1:c_RadBas>c_RadBig
  //                     // 2:auto-change the roles of LapBas,LapBig
  //                     // 4:f_Cb[0]>f_Cb[1]
  float rSttAveBig,    // average of original <denormalized> |LapBig|
        rSttNrmDelBig; // normalized delta of |LapBig| (it is ~rSttAveBig/c_RefBig)
  int   rSttCLABlr,    // actual d_Blr (changed by NrmDelBig)
        rStackFrmBuf;  // #used frames
  bool  rDnSample;     // spatial downsample 2x2 (for operations after "stage A")
  int   rAirInfFlg;    // shortest info from AirEstimation
} tUnSharpExts;
  // call-sequence in ~VaPrcFlu.cpp_jUnSharp():
  //   ~BS5::_UnSharpRecAdjust( out:tUnSharp, out:tUnSharpExts )
  //   high-level adhoc alteratios of tUnSharp,tUnSharpExts
  //   preprocess (like GrayScale cut)
  //   ~BS5::_UnSharp( inp:tUnSharp, inp:tUnSharpExts ):
  //     iUnSharpRecReConsider( out:tUnSharp, out:tUnSharpExts )
  //     stageA, DnSample, StageB..., UpSample

void BS5_UnSharpRecAdjust( int AW, int AH, const tUnSharp *AUSOrg, const tUnSharpAdjs *AUSAdjs,
   int AModeId, int AAdj, tUnSharp *AUS, tUnSharpExts *AUSExts ); // adjusts control-record, fills transient parameters

void BS5_UnSharp_CPU( tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH, tBS_Val *AFrm, tBS_Val **ATs,
  const char *ADbgPfx, char *lDbgSfx ); // UnSharp CPU

void BS5_UnSharp_CUDA( tUnSharp *AUS, tUnSharpExts *AUSExts, int AW, int AH, tBS_Val *AFrm, tBS_Val **ATsH, tBS_Val **ATsD,
  const char *ADbgPfx, char *lDbgSfx ); // UnSharp CUDA
#endif // I16_UT5_H_INCLUDED
