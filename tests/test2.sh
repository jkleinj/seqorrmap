#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test2
#_______________________________________________________________________________

../src/seqorrmap --pdbIn 1LKV.pdb --corr test2.dat --firstResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --rLoCutoff 0.4 --rUpCutoff 40.0 --corrOut test2_corr.dat --tclOut test2.vmd > >(tee test2.out) 2> >(tee test2.err >&2)

