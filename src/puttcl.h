/*===============================================================================
puttcl.h : Write TCL scripts
Copyright (C) 2013 Jens Kleinjung
see README file for more information 
================================================================================*/

#ifndef PUTTCL_H
#define PUTTCL_H

#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "getcorr.h"
#include "pdb_structure.h"

void print_tcl(Arg *arg, Str *pdbca, Corrlist *corrlist);

#endif
