@echo off
rem -- clean transient files of MS-Build
del/S/Q "*.ncb"
del/S/Q "*.aps"
del/S/Q/AH "*.suo"
rmdir/S/Q "Debug"
rmdir/S/Q "Release"
rmdir/S/Q "X64"
rmdir/S/Q ".vs"

