# -- It creates reference sub-directories "./l_*", which serve
# -- only to open editor-windows in VSCodeExplorer (avoiding buggy
# -- multi-root workspaces. ~Don't use *.sh if you don't use VSCode.
eacho "-- create symlink sub-dir for VSCode-explorer"
rm l_I16_Cmn
ln -s "../I16_Cmn" l_I16_Cmn
