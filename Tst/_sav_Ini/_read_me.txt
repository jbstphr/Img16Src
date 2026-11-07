  -- Typical folders -
Work ~"d:\MeGePrcs" 
Dbg  ~"..\" 

  -- Process-flow --
reduce dimension by PyrMxXY
shift grayscale right
[Gray Cut]
[COllimator]
[Thin Long Laplacians] (anti-scatter)
[Ribs Detection]
[UnSharp]
  [Air Estimation]
  [Quadrature Laplacians]
  [Mask re-Estimation]
shift grayscale left
restore dimension from PyrMxXY


  -- Typical parameters --
Work    | Dbg    | Key name     |
--------+--------+--------------+---------
[Settings]
 0        -1,1     nPrintLabel   Trace few text-logs in left-up corner by very dark|bright pixels.
 0        -1,1     nCUDA         0=auto.
      2..3         nPyrMxXY      Max downsample by X, by Y. It modifies process' dimensions to be
                                 a divident of respective power of two.
      0..1         nShrShu       Shift gray-scale values right before process (and restore back);
                               "=0" when (signed) inp gamut <="1 shl 14",
                               "=1" when (signed) inp gamut is "1 shl 15",
                               "=2" when unsigned inp gamut is "1 shl 16".
 0        1        nIniSetDflt   Overwrite certain sections by hard-coded defaults from "VaPrcFlu.inl".

[Gray Cut]
     0,1           nInvert       Set as suitable.
                   nGamInpOut    Inp,out gamuts. Current hard-code needs (~always) "16384 16384".
                   fCutDrkBrt    Histogram's cut; set as suitable ~"0.05 0.01".
     0             fGamma        [UnSharpen] does the most of gammas right before or right after
                                 "despeckle" as if "fGamma=0" (i.e. as if "=1" tehnically). However
                                 [Ribs Detection] acts before [US] and needs realistic gamma to be
                                 set here. Then hard-coded logic in [US] autonatically reshuffles
                                 its gammas when this "fGamma<>0".

                                 
[TLL]
     0             a,b,c,d       There are fits "a..d=bEn nRadX nRadY fAtt" for anti-scatter grid. 
                                 "bEn=0" disables respective fit.
