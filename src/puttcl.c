/*==============================================================================
puttcl.c : write Tcl scripts
Copyright (C) 2013 Jens Kleinjung
see README file for more information 
==============================================================================*/

#include "puttcl.h"

void print_tcl(Arg *arg, Str *pdbca, Corrlist *corrlist)
{
	unsigned int i, m;
	int *p_idx1 = 0; /* pointer to atom index */
	int *p_idx2 = 0;

	/* print header */
    fprintf(arg->tclOutFile,	"# seqorrmap TCL script for input files:\n#\t%s\n#\t%s\n#\n",
				arg->corrInFileName, arg->pdbInFileName);

	/* load structure */
    fprintf(arg->tclOutFile,	"# load molecule\n"
								"mol new %s\n"
								"mol rename top %s\n"
								"# basic selections and initial view\n"
								"set mymol [atomselect top \"all\"]\n"
								"set mysel [atomselect top \"alpha\"]\n",
				arg->pdbInFileName, arg->pdbInFileName);


	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
			/*for all mappings */
		for (m = 0; m < corrlist->corr[i].nCombimap; ++ m) {
			if (corrlist->corr[i].combimap[m].selected) {
				p_idx1= &(corrlist->corr[i].combimap[m].idx_pdbca[0]);
				p_idx2= &(corrlist->corr[i].combimap[m].idx_pdbca[1]);

				fprintf(arg->tclOutFile,	"# show correlation %d\n"
											"set corrsel0 [atomselect top \"resid %d and name %s and chain %s\"]\n"
											"set corrsel1 [atomselect top \"resid %d and name %s and chain %s\"]\n"
											"set pos0 [lindex [$corrsel0 get {x y z}] 0]\n"
											"set pos1 [lindex [$corrsel1 get {x y z}] 0]\n",
					i,
					pdbca->atom[*p_idx1].residueNumber,
					pdbca->atom[*p_idx1].atomName,
					pdbca->atom[*p_idx1].chainIdentifier,
					pdbca->atom[*p_idx2].residueNumber,
					pdbca->atom[*p_idx2].atomName,
					pdbca->atom[*p_idx2].chainIdentifier);

				if (corrlist->corr[i].corrval < 0.0)
					fprintf(arg->tclOutFile,	"draw color blue\n");
				if ((corrlist->corr[i].corrval >= 0.0) && (corrlist->corr[i].corrval < 0.2))
					fprintf(arg->tclOutFile,	"draw color gray\n");
				if ((corrlist->corr[i].corrval >= 0.2) && (corrlist->corr[i].corrval < 0.4))
					fprintf(arg->tclOutFile,	"draw color purple\n");
				if ((corrlist->corr[i].corrval >= 0.4) && (corrlist->corr[i].corrval < 0.6))
					fprintf(arg->tclOutFile,	"draw color red\n");
				if ((corrlist->corr[i].corrval >= 0.6) && (corrlist->corr[i].corrval < 0.8))
					fprintf(arg->tclOutFile,	"draw color orange\n");
				if ((corrlist->corr[i].corrval >= 0.8) && (corrlist->corr[i].corrval < 1.0))
					fprintf(arg->tclOutFile,	"draw color yellow\n");
				if (corrlist->corr[i].corrval >= 1.0)
					fprintf(arg->tclOutFile,	"draw color white\n");

				fprintf(arg->tclOutFile,	"draw line $pos0 $pos1 width 2\n");
			}
		}
	}
}
