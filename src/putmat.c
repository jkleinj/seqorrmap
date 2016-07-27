/*==============================================================================
putmat.c : Write corelation matrices 
Copyright (C) 2013 Jens Kleinjung
see README file for more information 
==============================================================================*/

#include "putmat.h"

/*____________________________________________________________________________*/
/* the correlation matrix is scaled to [-1 1] to fit with the Mutual Information
	values of the GSA Tools */
void print_corrmat(Arg *arg, Corrlist *corrlist, double **corrmat)
{
	unsigned int i;
	double corrval_sc = 0.; /* scaled [-1 1] correlation value */

	/* for all correlations */
	/* (all matrix values have already been zeroed,
		therefore it suffices to set the scaled correlation values) */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* scale correlation value */
		if (corrlist->corr[i].corrval >= 0) {
			corrval_sc = corrlist->corr[i].corrval / abs(corrlist->corrvalmax);
		} else {
			corrval_sc = corrlist->corr[i].corrval / abs(corrlist->corrvalmin);
		}
		/* assign correlation value to symmetric matrix */
		assert(corrlist->corr[i].corresnum[0] < (corrlist->corresmax + 1));
		assert(corrlist->corr[i].corresnum[1] < (corrlist->corresmax + 1));
		corrmat[corrlist->corr[i].corresnum[0]][corrlist->corr[i].corresnum[1]] = corrval_sc;
		corrmat[corrlist->corr[i].corresnum[1]][corrlist->corr[i].corresnum[0]] = corrval_sc;
	}

	print_mat2D_double(arg->corrmatOutFileName, corrmat, corrlist->corresmax, corrlist->corresmax);
}

/*____________________________________________________________________________*/
/* the complementary correlation matrix is intended to represent correlations
	as distances, meaning stronger correlations correspond to shorter distances;
	it is scaled to the range 1 - [0 1 */
void print_corrcompmat(Arg *arg, Corrlist *corrlist, double **corrcompmat)
{
	unsigned int i;
	double corrval_sc = 0.; /* scaled 1 - [0 1] correlation value */

	/* for all correlations */
	/* (all matrix values have already been zeroed,
		therefore it suffices to set the scaled and complemented correlation values) */
	for (i = 0; i < corrlist->nCorr; ++ i) {
		/* scale correlation value and assign complement value */
		corrval_sc = 1. - ((corrlist->corr[i].corrval - corrlist->corrvalmin) / abs(corrlist->corrvalmax));

		/* assign correlation value to symmetric matrix */
		corrcompmat[corrlist->corr[i].corresnum[0]][corrlist->corr[i].corresnum[1]] = corrval_sc;
		corrcompmat[corrlist->corr[i].corresnum[1]][corrlist->corr[i].corresnum[0]] = corrval_sc;
	}

	print_mat2D_double(arg->corrcompmatOutFileName, corrcompmat, corrlist->corresmax, corrlist->corresmax);
}
