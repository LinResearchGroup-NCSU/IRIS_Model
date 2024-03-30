#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Feb 15 22:24:24 2024
# File Name: cmd.extract_chains.sh
# Description: 
#########################################################################
#!/bin/bash

PDB_ID=$1

ExtractFromPDBFiles.pl --overwrite --mode Chains --chains A,B,R --CombineChains yes $PDB_ID.pdb
mv ${PDB_ID}ExtractedChains.pdb ${PDB_ID}_tmp.pdb
