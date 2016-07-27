/*===============================================================================
putmat.h : Write correlation matrices
Copyright (C) 2013 Jens Kleinjung
see README file for more information 
================================================================================*/

#ifndef PUTMAT_H
#define PUTMAT_H

#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "getcorr.h"
#include "matrix.h"

void print_corrmat(Arg *arg, Corrlist *corrlist, double **corrmat);
void print_corrcompmat(Arg *arg, Corrlist *corrlist, double **corrcompmat);

#endif
