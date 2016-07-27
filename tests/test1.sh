#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test1
#_______________________________________________________________________________

../src/seqorrmap --pdbIn 1LKV.pdb --corr test1.dat --firstResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --corrOut test1_corr.dat --tclOut test1.vmd > >(tee test1.out) 2> >(tee test1.err >&2)

