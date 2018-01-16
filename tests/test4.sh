#! /bin/bash
#_______________________________________________________________________________
# seqorrmap : test4
#_______________________________________________________________________________

valgrind --leak-check=full  --show-reachable=yes ../src/seqorrmap --pdbIn 1LKV.pdb --corr test4.dat --firstStrResidue 217 --cLoCutoff 0.4 --cUpCutoff 1.4 --rLoCutoff 0.4 --rUpCutoff 40.0 --sLoCutoff 2 --sUpCutoff 50 --select1 '217-224 228' --corrOut test4_corr.dat --tclOut test4.vmd > >(tee test4.o) 2> >(tee test4.e >&2)

