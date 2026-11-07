# Img16Src
            16-bit Image-Process (for Windows, Ubuntu)
                   The root source directory 

This repository utilizes "Occam axe" to avoid external dependencies and possible associations with any external (cpy-right/lefts). Here is just sole "2D-Imaging with Redundant Coordinates" (mainly utilizing "rough Radon Transform"). I cleanly state that this is completely personal work in "free personal time" for personal cognitive interest; without any external funding or specific external assistance from another person or organization. Possible future patent-application of the author might (try to) utilize "rough Radon Transform" for blind denoise. This code utilizes only "linear combinations of directional moving averages".

The explicit references are:
  - https://zenodo.org/records/18210101 "...2D-imaging in redundant coordinates",
  - K&R C RTL,
  - hermite-spline from Wiki,
  - moving average from Wiki.

The theoretical background of this repository might be associated with "geometrical measurement" by H.Federer, L.Simon and their colleagues from mainly in ~1970th; ~<1985. However this code utilizes exactly "adhoc non-non-regular" metrics;  hence such association could only illustrate certain possibilities but can't prove anything about formal correctness of this code/algorithms.

This directory contains ~all the code for "application-functionalities"  (beside K&R C RTL). There are

  = The library Img16Prc.so|.dll:
  
    "I16_Lib/*" - library's source
    "I16_Cmn/*" - common files to include in calling CPP-project
    
  = Process_execution_example & process_files_iterator Img16Prc<.exe>:
  
    "I16_Prc/*"
  
  = The directory for coding-debug-work-flow:
  
    "Tst/*"

  = The documentation:
  
  - I16_Prc/_ReadMe_I16Prc.pdf - systematic
  - I16_Prc/_ReadMe_I16Prc.txt - for reference during workflow
  - I16_Lib/_ReadMe_I16Lib.txt - for reference during workflow
   
   
Auxilliary associated directories/repositories:
 
  = 16-bit Headerless Images DataSeries Viewer MFC:
  
    "Img16DatSer/*"
    - It is suitable for Windows and Linux under WinE.
    - It deals with clipboard, drag&drop, WindowsExplorer&Nautilus.

  = The directory for application-debug-work-flow:
  
     like "e:\Img16Prc\*" or "/mnt/e/Img16Prc".
    - "*.bat/*.sh" for process group of input images by Img16Prc<.exe>
    - "Img16Prcs.ljm" to process single opened image im ImageJ's Java-AWT
