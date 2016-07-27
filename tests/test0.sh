#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test0
#_______________________________________________________________________________

../src/seqorrmap --pdbIn 1LKV.pdb --corr test0.dat --firstResidue 217 --corrOut test0_corr.dat --tclOut test0.vmd > >(tee test0.out) 2> >(tee test0.err >&2)

