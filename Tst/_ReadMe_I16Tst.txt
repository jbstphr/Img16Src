========================================================================
            16-bit Image-Process (for Windows, Ubuntu)
               The directory for debug-work-flow 
========================================================================

= The files' iterator:
  - Debug64/Img16Prc<.exe>.
  - Release64/Img16Prc<.exe>.


= Images' sub-directories xWWWWxHHHHsfx/*:
  - The DataSeries typically have signed 16-bit pixels.
  - Transient debug-trace-images './$a??_TBCS/*" also are "signed".
  - Certain DataSeries might be unsigned 16-bit pixels (e.g. x256MRI/*).


= Headerless 16-bit DataSeries' viewer:
  - Img16DatSer.exe (should work smooth under WinE in Ubuntu).
  - Img16DatSer.ini (transient initialization profile; could be
                     deleted at any moment ~without_a_harm).
  - Dedicated documentation "Img16DatSer.pdf".
  
  The notes:
    [Grayscale combo-box] typically has to be "signed".
    Unsigned pixels (e.g. "x256MRI/*") need [grayscale combo-box] "unsigned".


= Processing library:
  Img16Prc.so|.dll


= Processing library's initialization profiles Img16PrcNN.ini:
  - _sav_Ini/Img16PrcNN.ini - the templates (serve here and also for the
    application_evaluation directory ~"e:/Img16Prc/").
  - ./Img16PrcNN.ini - working copies for standard debug-work-flow
    and ad-hoc editing.

  The notes:
  A) Copy templates to this directory as/if/when required during debug-
     work-flow.

  B) If "nCUDA=0" then Img16Prc.so|.dll uses CUDA when "BC_DeviceFind():=true".
     Typically set "nCUDA=" to "-1" or "1". Then Img16Prc.so|.dll explicitly
     uses CPU or CUDA.
     _sav_Ini/Img16PrcNN.ini        ./Img16PrcNN.ini 
     [Settings]                     [Settings]
     nCUDA=0                        nCUDA=-1

  C) If "nIniSetDflt=0" then Img16Prc.so|.dll sets only non-existing keys
     at once (and further uses keys from ./Img16PrcNN.ini).
     Typically change "nIniSetDflt" to "1". Then Img16Prc.so|.dll always
     ovewrites unsharp-parameters by hard-coded initialization-values from
     "Img16Prc.inl tUnSharp cUSNN = {..." as ~:
     _sav_Ini/Img16PrcNN.ini         ./Img16PrcNN.ini
     [Settings]                      [Settings]
     nIniSetDflt=0                   nIniSetDflt=1
     [US]                            [US]
     nEnabled=1                      nEnabled=1
     <empty>                         nDeG1=...; ...
                                     fDeG4=... ... ... ...; ...
     

= Service batch files:
  - _clean.sh|.bat - delete processing results in xWWWWxHHHHsfx/*.
  - _cpyto_Img16Prc.sh|.bat - copy templates to application_evaluation
    directory  (~"e:/Img16Prc/").
