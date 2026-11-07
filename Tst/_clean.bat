echo off
rem - delete processing results
del *.raw
del/S/Q *_??.raw
del/S/Q *_????.raw
del *.exp
del *.lib
del *.pdb
for /f %i in ('dir /a:d /b $*') do rd /s /q %i")