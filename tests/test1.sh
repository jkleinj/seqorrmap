#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test1
#_______________________________________________________________________________

../src/seqorrmap --pdbIn 1LKV.pdb --corr test1.dat --firstStrResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --corrOut test1_corr.dat --tclOut test1.vmd > >(tee test1.o) 2> >(tee test1.e >&2)

