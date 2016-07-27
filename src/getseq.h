/*===============================================================================
getseq.h : read alignment of FASTA sequences
Copyright (C) 2004 Jens Kleinjung and John Romein
Read the COPYING file for license information.
================================================================================*/

#ifndef GETSEQ_H
#define GETSEQ_H

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "arg.h"
#include "safe.h"
#include "seq.h"

/*____________________________________________________________________________*/
void read_sequence(Arg *arg, Seq *seq);

#endif
