========================================================================
            16-bit Image-Process (for Windows, Ubuntu)
                     The library Img16Prc.so|.dll
========================================================================

////////////////////////////
// External documentation //
////////////////////////////
1. The introduction
1.1 The concept
1.2. App-specific dictionary
  The dictionary for computing
  The dictionary for presentation enhancement
2. Presentation enhancement
2.1. Input & output
2.2. The stages
2.3. CLAHE at stages D, E
2.4. Generation of normalized masks at stage F
2.5. LapAdaptive at stage H
2.6. Rough Radon's MinMax and complementary statistics at stage J
2.7. The "summa summarum"
2.8. Ad-hoc shift of isolines=0



/////////////////////////////////////////////
// Appendix A: Debug & Environment's notes //
/////////////////////////////////////////////

// Topic A1: 
///////////////////////////////////////////////////

// Topic A2: IDE settings (MS Visual Studio) 
///////////////////////////////////////////////////

        General|configuration type "dll", extension "dll"
        General|character set "Multi-Byte"
typical General|output directory "..\Tst\"
typical General|target name "Img16Prc" 
typical Linker|output file "$(OutDir)\$(TargetName).dll"
typical Debugging|working command "..\Tst\Debug64\Img16Prc.exe", ".\Release64\Img16Prc.exe"
typical Debugging|working directory "..\Tst" (should be absolute paㅐth actually)

Typical debug-flow:
  Build outs:  ~"..\Tst\Img16Prc.dll,Img16Prc.lib,...".
  Execute:     ~"..\Tst\Debug64\Img16Prc.exe 1500 1500 x1500x1500\0037.raw 0888".
  Work-dir:    (absolute path of) ~"..\Tst\" .
  Inp files:   ~"x1500x1500\????.raw".
  Out filess:  ~"x1500x1500\????_nn_0888.raw".
  Debug dumps: ~"$Tst_00_0888\*.raw".



///////////////////////////////
// Appendix B: the work-flow //
///////////////////////////////

// Topic B1: the semantics of principal transient frames
////////////////////////////////////////////////////////

  0          1           2            The notions:
A HiFrem     HiF         HiFblr     - high freqs preliminaries
                                    - optional dnsample
B GammaSpa   GammaCLA    GammaAtt   - gamma-correction;
                                      |Der3x3| in #3
C BlrBas     BlrBig      |LapBig|   - basic spatial metrics;
                                      LapBas in #3, LapBig in #4
D CLAHE      mskBodyOne  CLAInc     - basic_CLAHE
E mskAtt                            - attenuation mask by segmentation_CLAHE
F mskThn     mskThk      mskGDi     - normalized masks;
                                      mskEdg in #3
H LapMid     BlrAda      LapAda     - LapAdaptive
J LapMM      LapAdaEff              - Radon & LapAdaEff
L mskEdg_blr                        - the edges; LapAcc+=Edg
M                                   - the ridges; LapAcc+=Rdg  
Q BkgPtr                            - <wait for> air/background estimation;
                                      it outputs early anyway
N                                   - global dilation;
                                      LapAdaEff*=mskGDiRes, HiF:=mskGDiRes
S BodySum    LapSum      ResSum     - summa summarum
U Res                               - histogram cut out
V ResBkGM                           - zero_background mask
  HiFrem     mskGDiRes   ResFrm     - optional upsample (stack of stage "A")
X                        ResFrm     - ResFrm+=HiFrem*mskGDiRes

"usage of n_Z": mskGDiRes :=
           (1-n_Z0-n_Z1-n_Z2)+
           n_Z0*mskThk+     Fd1.raw
           n_Z1*mskGDi+     Fd2.raw
           n_Z2*mskEdg.     Fd3.raw 

"usage of s_Za": BodySum :=
 (1-s_Za0-s_Za1)*HiFblr +   Ab2.raw
           s_Za0*BlrAda +   Hb1.raw
           s_Za1*GammaCLA + B1.raw
           s_Za2*HiFrem +   Ab0.raw
           s_Za3*LapMid.    Ha0.raw  // if s_Za3<0 then "LapMid_blr H1.raw"

"usage of s_Zb": LapSum :=
         s_Zb0*CLAInc +     Dc2.raw
         s_Zb1*LapAda +     Hb2.raw
         s_Zb2*LapAdaEff +  Nb_LapAdaEff_J1.raw
         1.0*LapAcc.        S_lapacc.raw


// Topic B2: Trace files
////////////////////////

A0inp         - I16_Invoke if (xShiftRgt) mcShiftRgt
CbCut         - jGrayCut 
CcGamma       - <rarely used alone>
  jUnSharp (..., DbgPfx="e" )
e0Inp0        - the input as-is; likely it is CbCut
            = after BS6_Despeckle or BS6_DespeckleCUDA
e1Desp0       - <trans> median's residual
e1Desp1       - <trans> average's residual
e1Desp2       - <trans> result's residual
            = stage A; high freqs preliminaries
eA0           - HiF <optional out>
eA1           - HiFrem
eA2           - HiFblr
            = stage B; gamma-correction
eB0           - GammaSpa
eB1           - GammaCLA
eB2           - GammaAtt
eB3           - |Der3x3|_nrm (of GammaCLA; 'cGSDer3x3')
            = stage C; basic sptial metrics
eC0           - BlrBas        
eC1           - BlrBig
eC2           - |LapBig|_nrm ('cGSLapBig')
eC3           - LapBas
eC4           - LapBig
            = stage D: CLAHE
eDc0          - plain_CLAHE slightly blurred
eDc1          - BodyOne (or just 8..9bit GammaCLA)
eDc2          - CLAInc (i.e. plain_CLAHE-GammaCLA)
eDc3          - |Der3x3CLA|_nrm (of plain_CLAHE,'cGSDer3x3CLA')
            = stage E: attenuate |LapBig| to LabBigEff
eEa_GSAtt0    - <opt> CLAHE_FrmNrm(GammaAtt)
eEa_GSAtt2    - GSAtt basic 
eEb_GSAdd0    - GSAtt0+GSAdd2
eEb_GSAdd2    - <trans> k0*|Der3x3|+k1*|Der3x3CLA|+k2*|LapBig|
eEc_mskAtt    - norm_mask(GSAdd0)
eEd_mskAttB1  - <opt> Ed_mskAtt*BodyOne
            = stage F: normalized masks
eFa0          - |Der3x3| of B1
eFa1          - |Der3x3CLA| of Eb0
eFb0          - GSh (global shape to make MskThk,MskEdg)
eFb1          - <opt> GDi (shape of global dilation, could be set later as ~GSh)
eFb2          - <opt> GThn (shape of thin_mask, could be derived later from GSh)
eFc0          - MskEdg:=GSh-MskThk
eFc1          - MskThk:=blur(GSh)
eFd0          - mskThn:=cut(GThn)
eFd1          - mskThk:=cut(GSh_blr)
eFd2          - mskGDi:=cut(GDi)
eFd3          - mskEdg:=cut(MskEdg)
eFd4          - mskEmL
eFe4          - mskEmL*Edg
            = stage H: LapAdaptive
eHa0          - LapMid:=Ca0-Ca1
eHa1          - LapMid_dil:=Fc1*Ha0
eHb1          - BlrAda:=Ca1+Ha1
eHb2          - LapAda:=B0-Hb1
eHc2          - <opt> Rng(Hb2+eFd3*h_Iso)
eHd2          - <opt> LapAda *= transient_mask
eHd3          - <opt> transient_mask = k0*mskThk+k1*mskGDi+k2*qmskEdg
            = stage J: Radon & LapAdaEff
eJa0          - Blr*KMM (~LapMM*KMM) 
eJa1          - Srp (~MinLap) 
eJa2          - Dev  (Radon deviation)
eJa3          - ENm (~EdgeNaive mimic lap)
eJb1          - LapAdaEff (~Blr+(Srp-Blr)*mskThn )
eJc1          - LapAdaEff += ENm_conducted*KENm
eJc3          - ENm_conducted*KENm (it is limit_ENm_by_Dev)
eJd           - dilate local peaks in LapAdaEff (small effect)
            = stage K: global dilation;
ekaAHiFdi1    - HiFrem_dil:=A1*mskGDi
ekbJgdi1      - LapAdaEff_dil:=Jc1*mskGDi
            = stage N: edge
eNa0          - MimicDat(LapEdg,mskEdg), where LapEdg:=RadHess<0?LapAdaEff:LapMM
eNb0          - ~Hessian
eNb1          - Na0 blurred
eNc0          - Edg:=MimicDat(Nb1,Nb0)
            = stage S: summa summarum
eS_lapacc     - mskEmL+Edg+QuadLaps
eS0           - BodySum:=(a0-a1)*A2+a0*Hb1+a1*B1+a2*HiFrem_dil
eS1           - LapSum:=b0*Jb1+b1*Ed2+Nb0
eS2           - ResSum:=S0+S1
            = stage U: histogram cut out
eU0           - histogram cut out, compensate most of negatives
Zres          - final




//////////////////////////////////////////////
// Appendix C: algorithms' excerption notes //
//////////////////////////////////////////////

// Topic C1: _MxLapXLimbo
/////////////////////////

= invocation ~:
 ~ in BS_Ut5.cpp,_UnSharpCPU{"RadMM,KMM,TGMM,ZMM,KAO"};
 ~ in BS_Ut6.cpp,_DePoisson{"1,KU,0,0,0.0f"}

= declaration ~: _MxLapXLimbo ( // MxLapStd with ~extended_limbo
  int Rad, float K, int TolGap, int Z0, float KAO, *Lap, **Ts )

= step_1: modify Lap; calc RadonDeviation:
  - "TolGap<0" - increment_mode "Lap+=K*MxLapStd";
  -      ">=0" - filter_mode "Rng(Lap,low_est-TolGap,high_est+TolGap)".

= the filter_mode utilizes the facts that: 
  - "MxLapStd" is somewhere between "Ave" & "MxLapExact";
  - Ave & MxLapExact are ~byproducts of MxLapStd;
  - low_est := min(Ave,MxLapExact), low_est := MxLapStd+(low_est-MxLapStd)*K;
  - high_est := max(Ave,MxLapExact), high_est := MxLapStd+(high_est-MxLapStd)*K;.
 
= step_2: optional auto-oscillator:
  - utilizes the fact that RadonDeviation is ~byproduct of MxLapStd;
  - "KAO=0" - skip this step;
  -   "<>0" - "Lap*=AutoOscillatorMsk", where
    AutoOscillatorMsk:="_Cut(RadonDeviation*|KAO|,0..512)=>128..512",

= Debug-trace frames in Ts (when "KAO<0"):
  - #0: AutoOscillatorMsk;
  - #1: RadonDeviation*|KAO|;
  - #2: the Lap after step_1 (before "*=AutoOscillatorMsk").
  Note: calling code sets negative "KAO" when ~"debug_name_prefix<>nil";
        "KAO in control_record ~tUnSharp" is positive.

= The notes:
  - The extended_limbo could mean 16 directions. But in current implementation this
  is just reuse of 8 directional averages (obtained as byproduct of MxLapStd).

  - It uses byproduct directional averages in TS and hence invokes _MxLapStd(Shear=0).



//// Appendix D: Collection of ad-hoc command-lines
///////////////////////////////////////////////////

256 256 x0256MRI/3o.raw 90 0888 24
256 256 x0256MRI/5o.raw 90 0888 24
1024 1024 x1024ifr/0011.raw 100 0888 24
324 256 x0324x0256ifr/01210.raw 102 0888 24
1280 960 x1280x0960bar/001.raw 200 0888 24
1472 2172 x1472x2172ndt/2050.raw 70 0888 24
1472 2172 x1472x2172ndt/3003.raw 71 0888 24
