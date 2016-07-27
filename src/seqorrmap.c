/*=============================================================================
seqorrmap: map seqorr contacts 
Copyright (C) 2012-2013 Jens Kleinjung
==============================================================================*/

#include "seqorrmap.h"

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	unsigned int i = 0;
	unsigned int j = 0;

	/* data structures */
	Arg arg; /* command line arguments */
	Argpdb argpdb; /* data structure for PDB command line arguments */
	Sel sel[2]; /* residue selections (two selection lists) */
	Str pdbca; /* PDB structure CA only*/
	Str pdbaa; /* PDB structure all atoms*/
	Corrlist corrlist; /* correlation list */
	Seq seq; /* gapped target sequence */
	double **corrmat = 0; /* correlation matrix */
	double **corrcompmat = 0; /* correlation complement matrix */
	int combiChain = 0; /* number of chain combinations (upper limit for mappings) */

    /*____________________________________________________________________________*/
    /** parse command line arguments */
	parse_args(argc, &(argv[0]), &arg, &argpdb);
	if (! arg.silent) fprintf(stdout, "\nResidue selections\n");
	/* selection lists */
	parse_selects(&arg, arg.select[i], &(sel[i]));

    /*____________________________________________________________________________*/
    /** read input structure */
	if (! arg.silent) fprintf(stdout, "\nInput structure\n");
	/* pdbca is used for all calculations referring to residue pair correlations */
	argpdb.coarse = 1;
	read_structure(&arg, &argpdb, &pdbca); /* CA only */
	/* pdbaa is used for shortest atom distance calculations between residue pairs */
	argpdb.coarse = 0;
	read_structure(&arg, &argpdb, &pdbaa); /* all atom */

    /*____________________________________________________________________________*/
    /** read correlation data */
	if (! arg.silent) fprintf(stdout, "\nInput residue pair correlations\n");
   	read_corresnums(&arg, &corrlist);

    /*____________________________________________________________________________*/
    /** if seqorr correlation list is for entire alignment, read gapped target sequence */
	if (arg.gapped) {
		fprintf(stderr, "Gapped sequence not yet implemented!\n");
		exit(1);
	}

    /*____________________________________________________________________________*/
    /** adjust numbering of correlation list to gaps and firstResidue */
	if (! arg.silent) fprintf(stdout, "\nAdjusting residue numbering\n");
	adjust_corresnums(&arg, &corrlist, &seq);

    /*____________________________________________________________________________*/
    /** map correlations */
	arg.logoOutFile = safe_open(arg.logoOutFileName, "w");
	arg.logeOutFile = safe_open(arg.logeOutFileName, "w");

	if (! arg.silent) fprintf(stdout, "\nMapping correlations onto structure\n");
	fprintf(arg.logoOutFile, "\nMapping correlations onto structure\n");
	init_map(&corrlist, pdbca.nChain);
	map_corrs(&arg, &pdbca, &corrlist, &seq, pdbca.nChain);

	if (! arg.silent) fprintf(stdout, "\nCombining mapped correlations\n");
	combiChain = pdbca.nChain * pdbca.nChain;
	init_combimap(&corrlist, combiChain);
	combine_corrs(&arg, &pdbca, &corrlist, &seq);

    /*____________________________________________________________________________*/
    /** select and print correlations */
	if (! arg.silent) fprintf(stdout, "\nSelecting and printing mapped correlations\n");
	fprintf(arg.logoOutFile, "\nSelecting and printing mapped correlations\n");
	init_corrlist_selection(&corrlist);
	corrlist_selection(&arg, &pdbaa, &pdbca, &corrlist, &seq, &(sel[0]));
	print_corr(&arg, &pdbca, &corrlist);

	fclose(arg.logoOutFile);
	fclose(arg.logeOutFile);

    /*____________________________________________________________________________*/
	/* print TCL script */
	if (! arg.silent) fprintf(stdout, "\nCorrelation view is printed to '%s'.\n", arg.tclOutFileName);
	arg.tclOutFile = safe_open(arg.tclOutFileName, "w");
	print_tcl(&arg, &pdbca, &corrlist);
	fclose(arg.tclOutFile);

    /*____________________________________________________________________________*/
	/* print correlation matrix */
	if (! arg.silent) fprintf(stdout, "\nCorrelation matrix is printed to '%s'\n", arg.corrmatOutFileName);
	corrmat = alloc_mat2D_double(corrmat, (corrlist.corresmax + 1), (corrlist.corresmax + 1));
	init_mat2D_double(corrmat, (corrlist.corresmax + 1), (corrlist.corresmax + 1), 0.);
	print_corrmat(&arg, &corrlist, corrmat);

    /*____________________________________________________________________________*/
	/* print correlation complement matrix */
	if (! arg.silent) fprintf(stdout, "\nCorrelation complement matrix is printed to '%s'\n", arg.corrmatOutFileName);
	corrcompmat = alloc_mat2D_double(corrcompmat, (corrlist.corresmax + 1), (corrlist.corresmax + 1));
	init_mat2D_double(corrcompmat, (corrlist.corresmax + 1), (corrlist.corresmax + 1), 0.);
	print_corrcompmat(&arg, &corrlist, corrcompmat);

    /*____________________________________________________________________________*/
	/* free memory */
	/* PDB */
	free(pdbca.atom);
	free(pdbca.atomMap);
	free(pdbca.sequence.res);
	free(pdbca.sequence.name);

	free(pdbaa.atom);
	free(pdbaa.atomMap);
	free(pdbaa.sequence.res);
	free(pdbaa.sequence.name);

	/* free arguments */
	free(arg.select);

	/* corrlist */
	for (i = 0; i < corrlist.nCorr; ++ i) {
		free(corrlist.corr[i].map);
		free(corrlist.corr[i].combimap);
		free(corrlist.corr[i].adjcorresnam[j]);
	}
	free(corrlist.corr);

	/* selection */
	for (i = 0; i < 2; ++ i) {
		free(sel[i].chainSelection);
		free(sel[i].residueSelection);
	}

	/* correlation matrix */
	free_mat2D_double(corrmat, (corrlist.corresmax + 1));
	free_mat2D_double(corrcompmat, (corrlist.corresmax + 1));

    /*____________________________________________________________________________*/
    fprintf(stdout, "\nClean termination\n\n");
	return 0;
}

