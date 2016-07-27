/*==============================================================================
putpdbs.c : Routines for printing PDB structures and Ramchandran plots
Copyright (C) 2004 Jens Kleinjung
see README file for more information 
==============================================================================*/

#include <stdlib.h>
#include <stdio.h>

#include "putpdb.h"

/*----------------------------------------------------------------------------*/
/* calculate rmsd of atom pair */
__inline__ static float calc_rmsd(Atom *atom0, Atom *atom1)
{
    return sqrt(pow(atom0->pos.x - atom1->pos.x, 2)
              + pow(atom0->pos.y - atom1->pos.y, 2)
              + pow(atom0->pos.z - atom1->pos.z, 2));
}

/*----------------------------------------------------------------------------*/
/* write superimposed PDB files */
void print_pdb(FILE *pdboutfile, Str *str, float cutoff_radius, int atom)
{
	unsigned int i;

	for (i = 0; i < str->nAtom; ++i) {
		if (calc_rmsd(&str->atom[atom], &str->atom[i]) <= cutoff_radius) {
			fprintf(pdboutfile, "%30s", str->atom[i].description);
			fprintf(pdboutfile, "%8.3f%8.3f%8.3f\n",
				str->atom[i].pos.x, str->atom[i].pos.y, str->atom[i].pos.z);
		}
	}
	fprintf(pdboutfile, "TER\n");
}

