/*==============================================================================
getcorr.c :  read correlations
Copyright (C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "getcorr.h"

/*____________________________________________________________________________*/
/** compile pattern correlation format */
static void compile_pattern_corrform(regex_t *regex)
{
    char matchPattern[256];

    /* pattern for correlation format */
    strcpy (matchPattern, "[[:digit:]]+[\t][[:digit:]]+[\t][[:print:]]{1}[\t][[:digit:]]+[[\t][[:digit:]]+[\t][[:print:]]{1}[\t][[:print:]]+");

    assert(regcomp(regex, matchPattern, REG_EXTENDED) == 0); 
}

/*____________________________________________________________________________*/
/** match pattern */
static int match_pattern(regex_t *regex, char *searchPattern)
{
    return regexec(regex, searchPattern, 0, NULL, 0); 
}

/*___________________________________________________________________________*/
void read_corresnums(Arg *arg, Corrlist *corrlist)
{
    unsigned int allocated = 64; 
	char line[80];
	regex_t corrform; /* regular expression of correlation format */

    /*___________________________________________________________________________*/
    /* compile correlation format (match pattern) */
    compile_pattern_corrform(&corrform);

    /*___________________________________________________________________________*/
	/* initialise / allocate */
    corrlist->nCorr = 0; /* number of residue pair correlations */
    corrlist->corr = safe_malloc(allocated * sizeof(Corr));
	corrlist->corresmax = 0;
	corrlist->corrvalmin = DBL_MAX;
	corrlist->corrvalmax = DBL_MIN;

    /*___________________________________________________________________________*/
    /* read input data file */
    arg->corrInFile = safe_open(arg->corrInFileName, "r");

	while(fgets(line, 80, arg->corrInFile) != 0) {
/*
#ifdef DEBUG
		fprintf(stderr, "%s:%d: %d input: %s", __FILE__, __LINE__, corrlist->nCorr, line);
#endif
*/
		/* allocate space for residue names */
		corrlist->corr[corrlist->nCorr].adjcorresnam = safe_malloc(2 * sizeof(char [2]));

		/* assign correlations */
		if (match_pattern(&corrform, line) == 0) {
			assert((sscanf(&(line[0]), "%d%d%s%d%d%s%f\n",
				&(corrlist->corr[corrlist->nCorr].corresnum[0]),
				&(corrlist->corr[corrlist->nCorr].gapcorresnum[0]),
				&(corrlist->corr[corrlist->nCorr].adjcorresnam[0][0]),
				&(corrlist->corr[corrlist->nCorr].corresnum[1]),
				&(corrlist->corr[corrlist->nCorr].gapcorresnum[1]),
				&(corrlist->corr[corrlist->nCorr].adjcorresnam[1][0]),
				&(corrlist->corr[corrlist->nCorr].corrval)) == 7) && \
				"input data format has to be: [int] [int] [char] [int] [int] [char] [float]");
			corrlist->corr[corrlist->nCorr].adjcorresnam[0][1] = '\0';
			corrlist->corr[corrlist->nCorr].adjcorresnam[1][1] = '\0';
/*
#ifdef DEBUG
			fprintf(stderr, "%s:%d: %d assigned: %d\t%d\t%s\t%d\t%d\t%s\t%f\n",
				__FILE__, __LINE__,
				corrlist->nCorr,
				corrlist->corr[corrlist->nCorr].corresnum[0],
				corrlist->corr[corrlist->nCorr].gapcorresnum[0],
				corrlist->corr[corrlist->nCorr].adjcorresnam[0],
				corrlist->corr[corrlist->nCorr].corresnum[1],
				corrlist->corr[corrlist->nCorr].gapcorresnum[1],
				corrlist->corr[corrlist->nCorr].adjcorresnam[1],
				corrlist->corr[corrlist->nCorr].corrval);
#endif
*/

			/* determine maximal residue number in correlation list (lowest is 0) */
			corrlist->corresmax = corrlist->corresmax >= corrlist->corr[corrlist->nCorr].corresnum[0] ? corrlist->corresmax : corrlist->corr[corrlist->nCorr].corresnum[0];
			corrlist->corresmax = corrlist->corresmax >= corrlist->corr[corrlist->nCorr].corresnum[1] ? corrlist->corresmax : corrlist->corr[corrlist->nCorr].corresnum[1];

			/* determine minimal correlation value in correlation list */
			corrlist->corrvalmin = corrlist->corrvalmin <= corrlist->corr[corrlist->nCorr].corrval ? corrlist->corrvalmin : corrlist->corr[corrlist->nCorr].corrval;

			/* determine maximal correlation value in correlation list */
			corrlist->corrvalmax = corrlist->corrvalmax >= corrlist->corr[corrlist->nCorr].corrval ? corrlist->corrvalmax : corrlist->corr[corrlist->nCorr].corrval;

			/* count correlations */
			++ corrlist->nCorr;
		} else {
			fprintf(stderr, "Skipping incorrectly formatted input line: %s", line);
		}

		/*___________________________________________________________________________*/
		/* allocate more memory if needed */
        if (corrlist->nCorr == allocated) {
            allocated += 64;
            corrlist->corr = safe_realloc(corrlist->corr, (allocated * sizeof(Corr)));
		}
	}

	fclose(arg->corrInFile);

    /*___________________________________________________________________________*/
	assert((corrlist->nCorr > 0) && "need at least one correlation; check format of your input data");

    /*___________________________________________________________________________*/
    if (! arg->silent) fprintf(stdout, "\tCORR file: %s\n"
										"\tCORR file content:\n"
										"\t\tnCorr = %d\n",
						arg->corrInFileName, corrlist->nCorr);	

	regfree(&corrform);
}

