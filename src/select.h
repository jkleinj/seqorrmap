/*==============================================================================
select.h : selection list 
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef SELECT_H
#define SELECT_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "arg.h"
#include "safe.h"

/*____________________________________________________________________________*/
/* structures */
/* parameters for command line arguments */
typedef struct  
{
	char (*chainSelection)[1]; /* selected chain identifier */
	int (*residueSelection)[2]; /* selected residue number */
	int nSelection;
} Sel;

/*____________________________________________________________________________*/
/* prototypes */
void parse_selects(Arg *arg, char *select, Sel *sel);

#endif
