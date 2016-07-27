/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <limits.h>
#include <safe.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "argpdb.h"
#include "error.h"

/*____________________________________________________________________________*/
/* structures */
/* parameters for command line arguments */
typedef struct  
{
	FILE *pdbInFile;
	char *pdbInFileName;
	FILE *corrInFile;
	char *corrInFileName;
	FILE *gappedseqInFile;
	char *gappedseqInFileName;
	float rLoCutoff;
	float rUpCutoff;
	float cLoCutoff;
	float cUpCutoff;
	int sLoCutoff;
	int sUpCutoff;
	int firstSeqResidue;
	int firstStrResidue;
	int silent;
	int strict_aa;
	int strict_nr;
	int multiModel;
	FILE *pdbOutFile;
	char *pdbOutFileName;
	FILE *corrOutFile;
	char *corrOutFileName;
	FILE *tclOutFile;
	char *tclOutFileName;
	FILE *logoOutFile;
	char *logoOutFileName;
	FILE *logeOutFile;
	char *logeOutFileName;
	FILE *corrmatOutFile;
	char *corrmatOutFileName;
	FILE *corrcompmatOutFile;
	char *corrcompmatOutFileName;
	char **select;
	int gapped;
	int systemChain;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb);

#endif
