#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test3
#_______________________________________________________________________________

../src/seqorrmap --pdbIn 1LKV.pdb --corr test3.dat --firstResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --rLoCutoff 0.4 --rUpCutoff 40.0 --sLoCutoff 2 --sUpCutoff 50 --corrOut test3_corr.dat --tclOut test3.vmd > >(tee test3.out) 2> >(tee test3.err >&2)

