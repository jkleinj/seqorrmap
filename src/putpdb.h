/*===============================================================================
 putpdbs.h : Write alignments of PDB structures
 Copyright (C) 2006 Jens Kleinjung
 see README file for more information 
================================================================================*/

#ifndef PUTPDBS_H
#define PUTPDBS_H

#include "pdb_structure.h"

void print_pdb(FILE *pdboutfile, Str *pdb, float cutoff_radius, int atom);

#endif
