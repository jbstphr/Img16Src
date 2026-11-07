@echo off
rem -- It creates reference sub-directories "./l_*", which serve
rem -- only to open editor-windows in VSCodeExplorer (avoiding buggy
rem -- multi-root workspaces. ~Don't use *.bat if you don't use VSCode.
echo -- create symlink sub-dir for VSCode-explorer
@echo on
rmdir l_VPF_Cmn
mklink /D l_I16_Cmn "..\I16_Cmn"


