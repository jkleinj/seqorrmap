/*==============================================================================
mapcorr.c : map correlations
Copyright (C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "mapcorr.h"

/*____________________________________________________________________________*/
/** adjust residue numbering of correlations
	according to firstResidue and cumulative gaps */
void adjust_corresnums(Arg *arg, Corrlist *corrlist, Seq *seq)
{
	unsigned int i, j;

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* for the residue pair */
		for (j = 0; j < 2; ++ j) {
			if (arg->gapped) {
				Error("Gapped mapping is not yet fully implemented!\n");
				/* flag up gapped target sequence positions */
				if (seq->dgap[corrlist->corr[i].corresnum[j]]) {
					corrlist->corr[i].adjcorresnum[j] = INT_MIN;
				/* or assign adjusted residue number */
				} else {
				/* if sequence gapped, adjust for number of gaps and firstResidue shift */
					corrlist->corr[i].adjcorresnum[j] = corrlist->corr[i].corresnum[j] - seq->cgap[i] + arg->firstStrResidue;
				}
				/* transform according to firstStrResidue and firstSeqResidue */
			} else {
				corrlist->corr[i].adjcorresnum[j] = corrlist->corr[i].corresnum[j] + arg->firstStrResidue - arg->firstSeqResidue;
			}
		}
	/*
#ifdef DEBUG
		fprintf(stderr, "%s:%d: %d : %d -> %d | %d -> %d\n",
			__FILE__, __LINE__,
			i,
			corrlist->corr[i].corresnum[0],
			corrlist->corr[i].adjcorresnum[0],
			corrlist->corr[i].corresnum[1],
			corrlist->corr[i].adjcorresnum[1]);
#endif
	*/
	}
}

/*____________________________________________________________________________*/
/** initialise mapping of correlations */
void init_map(Corrlist *corrlist, int nChain)
{
	unsigned int i, m, j;

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* allocate memory for multiple mappings of one correlation onto PDB chains */
		corrlist->corr[i].map = safe_malloc(nChain * sizeof(Map));

		/* for all potential mappings */
		for (m = 0; m < nChain; ++ m) {
			/* for the correlated residue pair */
			for (j = 0; j < 2; ++ j) {
				corrlist->corr[i].map[m].idx_pdbca[j] = INT_MIN;
				corrlist->corr[i].map[m].idx_pdbaa[j] = INT_MIN;
				corrlist->corr[i].map[m].rmsdca = FLT_MAX;
			}
		}
		corrlist->corr[i].nMap = 0;
	}
}

/*____________________________________________________________________________*/
/** map correlations */
void map_corrs(Arg *arg, Str *pdbca, Corrlist *corrlist, Seq *seq, int allocated)
{
	unsigned int i, a, j;
	int *p_nMap; /* pointer to nMap */
	unsigned int nMatch[2];

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		fprintf(arg->logoOutFile, "\tcorrelation %d:\n", i);

		p_nMap = &(corrlist->corr[i].nMap);
		nMatch[0] = 0;
		nMatch[1] = 0;
		/* for all atoms */
		for (a = 0; a < pdbca->nAtom; ++ a) {
			/* for both correlation residues */
			for (j = 0; j < 2; ++ j) {
				if (pdbca->atom[a].residueNumber == corrlist->corr[i].adjcorresnum[j] &&
					pdbca->sequence.res[a] == corrlist->corr[i].adjcorresnam[j][0]) { 
					corrlist->corr[i].map[*p_nMap].idx_pdbca[j] = a; /* assign atom index */
					++ nMatch[j]; /* count individual residue match */
				
					fprintf(arg->logoOutFile, "\t\tresidue%d %d, chain %c, index %d\n",
									j+1, pdbca->atom[a].residueNumber,
									pdbca->atom[a].chainIdentifier[0], a);
				}
			}
			/* if both residues have been matched at least once */
			if ((nMatch[0] > 0) && (nMatch[1] > 0)) {
				++ (*p_nMap);
				nMatch[0] = 0;
				nMatch[1] = 0;
			}
		}
	}
}

/*____________________________________________________________________________*/
/** initialise mapping of combined correlations */
void init_combimap(Corrlist *corrlist, int combiChain)
{
	unsigned int i, m, j;

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* allocate memory for multiple mappings of one correlation onto PDB chains */
		corrlist->corr[i].combimap = safe_malloc(combiChain * sizeof(Map));

		/* for all potential mappings */
		for (m = 0; m < combiChain; ++ m) {
			/* for the correlated residue pair */
			for (j = 0; j < 2; ++ j) {
				corrlist->corr[i].combimap[m].idx_pdbca[j] = INT_MIN;
				corrlist->corr[i].combimap[m].idx_pdbaa[j] = INT_MIN;
				corrlist->corr[i].combimap[m].rmsdca = FLT_MAX;
			}
		}
		corrlist->corr[i].nCombimap = 0;
	}
}

/*____________________________________________________________________________*/
/** map correlations */
void combine_corrs(Arg *arg, Str *pdbca, Corrlist *corrlist, Seq *seq)
{
	unsigned int i, m, n;
	int *p_idx0; /* pointer to atom index */
	int *p_idx1;
	int nCombimap = 0;

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* for all mappings */
		for (m = 0, nCombimap = 0; m < corrlist->corr[i].nMap; ++ m) {
			p_idx0 = &(corrlist->corr[i].map[m].idx_pdbca[0]);
			/* and their combinations */
			for (n = 0; n < corrlist->corr[i].nMap; ++ n) {
				p_idx1 = &(corrlist->corr[i].map[n].idx_pdbca[1]);
				/* assign combined mappings */
				corrlist->corr[i].combimap[nCombimap].idx_pdbca[0] = *p_idx0;			
				corrlist->corr[i].combimap[nCombimap].idx_pdbca[1] = *p_idx1;			
/*
#ifdef DEBUG
				fprintf(stderr, "===> %d %d %d: %d %d: %d %d\n", 
						i, m, n,
						corrlist->corr[i].nMap, nCombimap,
						*p_idx0, *p_idx1);
#endif
*/
				++ nCombimap;
			}
		}
		corrlist->corr[i].nCombimap = nCombimap;
	}
}

/*____________________________________________________________________________*/
/** initialise selection of correlation list */
void init_corrlist_selection (Corrlist *corrlist)
{
	unsigned int i, m;

	corrlist->ssel = 0;
	corrlist->csel = 0;
	corrlist->rsel = 0;
	corrlist->lsel = 0;
	corrlist->selected = 0;

	for (i = 0; i < corrlist->nCorr; ++ i) {
		for (m = 0; m < corrlist->corr[i].nCombimap; ++ m) {
			corrlist->corr[i].combimap[m].ssel = 0;
			corrlist->corr[i].combimap[m].csel = 0;
			corrlist->corr[i].combimap[m].rsel = 0;
			corrlist->corr[i].combimap[m].lsel = 0;
			corrlist->corr[i].combimap[m].selected = 0;
		}
	}
}

/*____________________________________________________________________________*/
/* check whether int value is within cutoff */
__inline__ static int within_cutoff_int(int value, int loCutoff, int upCutoff)
{
	if ((loCutoff <= value) && (value <= upCutoff))
		return 1;
	else
		return 0;
}

/*____________________________________________________________________________*/
/* check whether int value is within cutoff */
__inline__ static int within_cutoff_float(float value, float loCutoff, float upCutoff)
{
	if ((loCutoff <= value) && (value <= upCutoff))
		return 1;
	else
		return 0;
}

/*____________________________________________________________________________*/
/* check whether selected and system chain identifiers match */
__inline__ static int within_system_chain(int systemChain, char *selChain)
{
    if (systemChain == (int)selChain[0])
        return 1;
    else
        return 0;
}

/*____________________________________________________________________________*/
/* check whether selected and PDB chain identifiers match */
__inline__ static int within_pdb_chain(char *molChain, char *selChain)
{
    if (molChain[0] == selChain[0])
        return 1;
    else
        return 0;
}


/*____________________________________________________________________________*/
__inline__ static int within_list(Arg *arg, Str *pdbca, Corrlist *corrlist, Sel *sel, 
					unsigned int i, unsigned int m)
{
	unsigned int j, k;
	int nSel[2];
	int *p_idx = 0; /* pointer to atom index */
	char *p_chain = 0; /* pointer to chain */

	nSel[0] = 0;
	nSel[1] = 0;

	for (j = 0; j < 2; ++ j) {
		/* selection list is empty: select all residues */
		if (! sel[j].nSelection) {
			++ nSel[j];
		/* selection list is not empty: select specified residues */
		} else {
			p_idx = &(corrlist->corr[i].combimap[m].idx_pdbca[j]);

			/* for all selection groups */
			for (k = 0; k < sel[j].nSelection; ++ k) {
				/* use program-internal chain definition if specified */
				if (arg->systemChain) {
					p_chain = &(pdbca->atom[*p_idx].systemChainID);
				/* else use PDB chain definition */
				} else {
					p_chain = &(pdbca->atom[*p_idx].chainIdentifier[0]);
				}

				/* if this correlated residue pair is within the range of selected residues
					and the selected chain identifier matches that of the structure as well,
					select this pair (further down) */
					/* 'j' refers to the selection list and selects automatically
						the first residue/chain for the first selection list and
						the second residue/chain for the second selection list */
/*
#ifdef DEBUG
				fprintf(stderr, "====> %d:%d:%d:%d: is %d | sel %d %d\tis %c | sel %c\n",
								i, m, j, k,
								corrlist->corr[i].adjcorresnum[j],
								sel[j].residueSelection[k][0],
								sel[j].residueSelection[k][1],
								*p_chain,
								sel[j].chainSelection[k][0]);
# endif
*/
				within_system_chain(*p_chain,
									&(sel[j].chainSelection[k][0]));
				/* verify residue */
				if (within_cutoff_int(corrlist->corr[i].adjcorresnum[j], 
										sel[j].residueSelection[k][0], 
										sel[j].residueSelection[k][1]) &&
					within_system_chain(*p_chain,
										&(sel[j].chainSelection[k][0]))) {
					++ nSel[j];
					break;
				}
			}
		}
	}

	/* if both correlation residues have been selected,
		either explicitely in one of the two selection lists
		or implicitly by the absence of either selection list,
		flag this correlation up as being selected */
	if ((nSel[0] > 0) && (nSel[1] > 0)) {
		fprintf(arg->logoOutFile, "\tselected correlation %d, mapping %d\n", i, m);
		++ corrlist->corr[i].combimap[m].lsel;
/*
#ifdef DEBUG
		fprintf(stderr, "===> +\n");
#endif
*/
		return 1;
	} else {
		return 0;
	}
}

/*____________________________________________________________________________*/
/** select correlations */
void corrlist_selection(Arg *arg, Str *pdbaa, Str *pdbca, Corrlist *corrlist, Seq *seq, Sel *sel)
{
	unsigned int i, m;
	int *p_idx0; /* pointer to atom index */
	int *p_idx1;

	if (! arg->silent) fprintf(stdout, "\tapplying cutoff filters\n");

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* for all mappings */
		for (m = 0; m < corrlist->corr[i].nCombimap; ++ m) {
			p_idx0 = &(corrlist->corr[i].combimap[m].idx_pdbca[0]);
			p_idx1 = &(corrlist->corr[i].combimap[m].idx_pdbca[1]);

			/*____________________________________________________________________________*/
			/* sequence (distance) cutoff filter */
/*
#ifdef DEBUG
			fprintf(stderr, "=> %d %d %d\n",
					(corrlist->corr[i].corresnum[1] - corrlist->corr[i].corresnum[0]),
					arg->sLoCutoff, arg->sUpCutoff);
#endif
*/
			if (within_cutoff_int((corrlist->corr[i].corresnum[1] - corrlist->corr[i].corresnum[0]), 
									arg->sLoCutoff, arg->sUpCutoff)) {
				++ (corrlist->corr[i].combimap[m].ssel);
				++ (corrlist->ssel);
			}
			/*____________________________________________________________________________*/
			/* correlation (strength) cutoff filter */
/*
#ifdef DEBUG
			fprintf(stderr, "==> %f %f %f\n",
						corrlist->corr[i].corrval,
						arg->cLoCutoff, arg->cUpCutoff);
#endif
*/
			if (within_cutoff_float(corrlist->corr[i].corrval,
									arg->cLoCutoff, arg->cUpCutoff)) {
				++ corrlist->corr[i].combimap[m].csel;
				++ corrlist->csel;
			}
			/*____________________________________________________________________________*/
			/* CA atom distance */
/*
#ifdef DEBUG
			fprintf(stderr, "==> %f %f %f\n",
						corrlist->corr[i].map[m].rmsdca,
						arg->rLoCutoff, arg->rUpCutoff);
#endif
*/
			corrlist->corr[i].combimap[m].rmsdca = v_rmsd(&(pdbca->atom[*p_idx0].pos), &(pdbca->atom[*p_idx1].pos));
			if (within_cutoff_float(corrlist->corr[i].combimap[m].rmsdca, arg->rLoCutoff, arg->rUpCutoff)) {
					++ corrlist->corr[i].combimap[m].rsel;
					++ corrlist->rsel;
			}
			/*____________________________________________________________________________*/
			/* selection list cutoff filter */
			if (within_list(arg, pdbca, corrlist, sel, i, m)) {
				++ corrlist->corr[i].combimap[m].lsel;
				++ corrlist->lsel;
			}

			/*____________________________________________________________________________*/
			/* combine all cutoff criteria */
			if (corrlist->corr[i].combimap[m].ssel && corrlist->corr[i].combimap[m].csel && \
			    corrlist->corr[i].combimap[m].rsel && corrlist->corr[i].combimap[m].lsel) {
				++ corrlist->corr[i].combimap[m].selected;
				++ corrlist->selected;
			}
		}
	}

	/*____________________________________________________________________________*/
	/* report selection statistics */
	if (! arg->silent) fprintf(stdout,	"\t\tsequence distance: %d\n"
										"\t\tcorrelation strength: %d\n"
										"\t\tRMSD: %d\n"
										"\t\tselection: %d\n",
								corrlist->ssel,
								corrlist->csel,
								corrlist->rsel,
								corrlist->lsel);
	if (! arg->silent) fprintf(stdout,	"\tcombined: selected %d mappings from %d correlation(s)\n",
								corrlist->selected, corrlist->nCorr);
}

/*____________________________________________________________________________*/
/* print correlation */
void print_corr(Arg *arg, Str *pdbca, Corrlist *corrlist)
{
	unsigned int i, m, n;

	int *p_idx1 = 0; /* pointer to atom index */
	int *p_idx2 = 0;

	arg->corrOutFile = safe_open(arg->corrOutFileName, "w");
	fprintf(arg->corrOutFile, "# seqResNr1 strResNr1 strResNe1 strChaNe1 | seqResNr2 strResNr2 strResNe2 strChaNe2 | seqDist rmsdCA/A corrVal\n");

	/* for all correlations */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* for all mappings */
		for (m = 0; m < corrlist->corr[i].nCombimap; ++ m) {
			if (corrlist->corr[i].combimap[m].selected) {
				/* and their combinations */
				for (n = 0; n < corrlist->corr[i].nCombimap; ++ n) {
					if (corrlist->corr[i].combimap[n].selected) {
						p_idx1= &(corrlist->corr[i].combimap[m].idx_pdbca[0]);
						p_idx2= &(corrlist->corr[i].combimap[n].idx_pdbca[1]);

						fprintf(arg->corrOutFile, "%d\t%d\t%s\t%s\t%s\t|\t%d\t%d\t%s\t%s\t%s\t|\t%d\t%f\t%f\n",
							corrlist->corr[i].corresnum[0],
							pdbca->atom[*p_idx1].residueNumber,
							pdbca->atom[*p_idx1].atomName,
							pdbca->atom[*p_idx1].residueName,
							pdbca->atom[*p_idx1].chainIdentifier,
							corrlist->corr[i].corresnum[1],
							pdbca->atom[*p_idx2].residueNumber,
							pdbca->atom[*p_idx2].atomName,
							pdbca->atom[*p_idx2].residueName,
							pdbca->atom[*p_idx2].chainIdentifier,
							(corrlist->corr[i].corresnum[1] - corrlist->corr[i].corresnum[0]),
							corrlist->corr[i].combimap[m].rmsdca,
							corrlist->corr[i].corrval);
					}
				}
			}
		}
	}

	fclose(arg->corrOutFile);
}

