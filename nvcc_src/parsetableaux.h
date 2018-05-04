#ifndef PARSETABLEAUX_H
#define PARSETABLEAUX_H
/*****************************************************************************
 * 
 * File:    parsetableaux.c
 * Author:  Alex Stivala
 * Created: January 2010
 *
 * Declaratinos and definitions for tableaux/dist matrix parsing functions
 *
 * $Id: parsetableaux.h 3328 2010-02-12 06:39:40Z alexs $
 *
 *****************************************************************************/

#include "saparams.h"


#define MAX_LINE_LEN 2048


/* Index into 2d m x n array stored in contiguous memory */
#define INDEX2D(i,j,m,n) ( ((i)*(n) + (j)) )

/* parse a 2-char encoded tableau */
int parse_tableau(FILE *fp, size_t dim, int n, char *tab);

/* parse distance matrix */
int parse_distmatrix(FILE *fp, size_t dim, int n, float *dmat, int discard);

/* read whole database into memory */
int read_database(FILE *fp, char **tableaux, float **distmatrices, 
                  char **large_tableaux, float **large_distmatrices,
                  int **orders, char **names,
                  int **large_orders, char **large_names,
                  int *num_large_read);
#endif /* PARSETABLEAUX_H */
