/*==============================================================================
mapcorr.h : map correlations
Copyright (C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef MAPCORR_H
#define MAPCORR_H

#include <assert.h>
#include <safe.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "arg.h"
#include "getcorr.h"
#include "getpdb.h"
#include "select.h"
#include "seq.h"
#include "vector.h"

/*___________________________________________________________________________*/
/* prototype */
void adjust_corresnums(Arg *arg, Corrlist *corrlist, Seq *seq);
void init_map(Corrlist *corrlist, int nChain);
void init_combimap(Corrlist *corrlist, int combiChain);
void map_corrs(Arg *arg, Str *str, Corrlist *corrlist, Seq *seq, int allocated);
void combine_corrs(Arg *arg, Str *str, Corrlist *corrlist, Seq *seq);
void init_corrlist_selection(Corrlist *corrlist);
void corrlist_selection(Arg *arg, Str *pdbaa, Str *pdbca, Corrlist *corrlist, Seq *seq, Sel *sel);
void print_corr(Arg *arg, Str *pdbca, Corrlist *corrlist);

#endif

