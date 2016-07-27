/*===============================================================================
getcorr.h : read correlations 
Copyright (C) 2013 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETCORR_H
#define GETCORR_H

#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arg.h"
#include "error.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* structures */
/* correlation list */
/* information about correlations that have been mapped to PDB i residues
 * multiple mappings are allowed, for example in multi-chain structures */
typedef struct
{
	int idx_pdbca[2]; /* atom indices of matched pdbca residues */ 
	float rmsdca; /* distance between CA atoms */
	int idx_pdbaa[2]; /* atom indices of matched pdbaa residues */ 
	int ssel; /* within sequence distance cutoffs */
	int csel; /* within correlation strength cutoffs */
	int rsel; /* within radial distance cutoffs */
	int lsel; /* within listed selections (command line) */
	int selected; /* selected correlation */
} Map;

typedef struct
{
	int corresnum[2]; /* residue numbering correlation pair list */
	int gapcorresnum[2]; /* residue numbering including gaps */
	int adjcorresnum[2]; /* adjusted numbering (considering structure and sequence) */
	char (*adjcorresnam)[2];
	Map *map; /* atom index of matched pdbaa residue */ 
	Map *combimap; /* atom index of matched pdbaa residue */ 
	float corrval; /* correlation strength */
	int nMap; /* times this correlation has been mapped */
	int nCombimap;
} Corr;

/* array of (input) correlations */
typedef struct
{
	Corr *corr; /* correlations in list */
	int nCorr; /* number of correlations in list */
	int corresmax; /* highest residue number in correlation list (lowest is 0) */
	double corrvalmin; /* lowest correlation value in correlation list */
	double corrvalmax; /* highest correlation value in correlation list */
	int ssel; /* within sequence distance cutoffs */
	int csel; /* within correlation strength cutoffs */
	int rsel; /* within radial distance cutoffs */
	int lsel; /* within listed selections (command line) */
	int selected; /* selected correlation */
} Corrlist;

/*____________________________________________________________________________*/
/* prototypes */
void read_corresnums(Arg *arg, Corrlist *corrlist);

#endif
