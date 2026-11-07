========================================================================
            16-bit Image-Process (for Windows, Ubuntu)
          Exec example & files' iterator Img16Prc<.exe>.
========================================================================

This file-iterator executes ~various I16Prc.so|.dll.
Refer _ReadMe_I16Prc.rtf.

  // Quick notes:
  ///////////////
Note1: this code ~frequently invokes "flush(stdout)"; a I16Prc.so|.dll shall do 
       the same (as the stdout here might be routinely redirected|rebuffered
       by a caller, wich executes this utility as files' iterator).
Note2: code always deals with ASCII, x64.


  // A collection of command-lines:
  /////////////////////////////////
768 768 c:/I16/CAL0768x0768/x??????.bin 3 0 24
768 768 c:/I16/CAL0768x0768/A??_?????.bin 4 0 24