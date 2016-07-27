/*===============================================================================
seqorrmap.h : map seqorr correlations onto structure 
Copyright (C) 2012 Jens Kleinjung
================================================================================*/

#ifndef SEQORRMAP_H
#define SEQORRMAP_H

#include <assert.h>
#include <ctype.h>
/*#include <gsl/gsl_math.h>*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "arg.h"
#include "error.h"
#include "getcorr.h"
#include "getpdb.h"
#include "getseq.h"
#include "mapcorr.h"
#include "matrix.h"
#include "putmat.h"
#include "puttcl.h"
#include "safe.h"
#include "select.h"

/*____________________________________________________________________________*/
/* rmsd of aligned atom pairs */
typedef struct
{
	float rmsd; /* rmsd value of atom pair */
	int pos0; /* flags for match type */
	int pos1;
} Atompair;

#endif
