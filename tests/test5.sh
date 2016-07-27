#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test5
#_______________________________________________________________________________

valgrind ../src/seqorrmap --pdbIn 1LKV.pdb --corr test5.dat --firstResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --rLoCutoff 0.4 --rUpCutoff 40.0 --sLoCutoff 2 --sUpCutoff 50 --regionStart 230 --regionEnd 320 --gappedSeq test5.fasta --corrOut test5_corr.dat --tclOut test5.vmd > >(tee test5.out) 2> >(tee test5.err >&2)

