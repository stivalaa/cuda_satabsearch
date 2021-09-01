/*****************************************************************************
 * 
 * File:    parsetableaux.c
 * Author:  Alex Stivala
 * Created: January 2010
 *
 *
 * Functions to parse tableaux and distance matrices.
 * Equivalent to the FORTRAN subroutines RDTABD and RDISTM.
 * Note that unlike FORTRAN, the 2-character encoding of tableaux is not
 * convenient in C, we will still use 2-dimensional arrays (actual 2D arrays
 * not C style arrays of pointers) but encode each 2-character tableaux
 * code in 1 char according to:
 *
 * top 4 bits:     P    0
 *                 R    1
 *                 O    2
 *                 L    3
 *
 * bottom 4 bits:  E    0 
 *                 D    1 
 *                 S    2
 *                 T    3
 *
 * so OT for example is 0x22, LE is 0x30, etc.
 *                 
 *
 * The SSE types are encoded in one byte as:
 *
 * 00 e     beta strand
 * 01 xa    alpha helix
 * 02 xi    pi helix
 * 03 xg    3_10 helix
 *
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parsetableaux.h"

/*
 * encode the given 2-char SSE type code as per the scheme in header
 * comment
 *
 * Parameters:
 *    ssecode - 2 char SSE code
 * Return value:
 *   single byte encoding of the SSE code
 */
static char encode_ssetype(char ssetype[2])
{
  char ssecode;
  
  if (ssetype[0] == 'e')
    ssecode = 0;
  else
    switch (ssetype[1])
    {
      case 'a':
        ssecode = 1;
        break;
      case 'i':
        ssecode = 2;
        break;
      case 'g':
        ssecode = 3;
        break;
      default:
        fprintf(stderr, "Bad helix type %c\n", ssetype[1]);
        exit(1);
        break;
    }
  return ssecode;
}


/*
 * encode the given 2-char tableaux code as per the scheme in header
 * comment
 *
 * Parameters:
 *     tabcode - 2 char tableaux code
 * Return value
 *     single byte encoding of the tableaux code
 */
static char encode_tabcode(char tabcode[2])
{
  char high,low,codeval;

  switch (tabcode[0])
  {
    case 'P':
      high = 0;
      break;
    case 'R':
      high = 1;
      break;
    case 'O':
      high = 2;
      break;
    case 'L':
      high = 3;
      break;
    case '?':
      high = 4; /* sometimes get ?? when some error computing angles */
      break;
    default:
      fprintf(stderr, "invalid tableaux code %c\n", tabcode[0]);
      exit(1);
      break;
  }
  
  switch (tabcode[1]) 
  {
    case 'E':
      low = 0;
      break;
    case 'D':
      low = 1;
      break;
    case 'S':
      low = 2;
      break;
    case 'T':
      low = 3;
      break;
    case '?':
      low = 4; /* sometimes get ?? when some error computing angles */
      break;
    default:
      fprintf(stderr, "invalid tableaux code %c\n", tabcode[1]);
      exit(1);
      break;
  }

  codeval = (high << 4) | low;
  return codeval;
}

/*
 * Read tableau database entry in discrete tableaux format.
 *
 * The format of the 'database' is a text file with an entry for each
 * structure.
 * The first line of an entry is the identifier and
 * order of tableau (i.e. dimension of square array), then
 * each subsequent row is a row of the tableau, lower triangle
 * only (since it is symmetric).
 * The diagonal entries are meaningless (self-angle) in tableaux,
 * and are included instead to specify the SSE type, with
 * the following codes:
 *
 * e     beta strand
 * xa    alpha helix
 * xi    pi helix
 * xg    3_10 helix
 *
 * Width of identifier is 8 chars, blank padded on right,
 * width of order is 4 digits, blank padded on left.
 * There is a single space between identifier and order.
 * Each entry in tableau is two characters, with a space betwen
 * each on a line, and one line
 * per row of matrix.
 *
 * E.g.:
 *
 * /local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii
 *  T F
 * D1UBIA_    6
 * e  
 * OT e  
 * LE RT xa 
 * RT LE RT e  
 * LE RD LE OT e  
 * PE RT LE OT PE e  
 *
 * The tableau is followed by a SSE distance matrix, not read
 * by this subroutine, rather handled by a call to  parse_distmatrix()
 *
 * Parameters:
 *      fp - filehandle to read from
 *      dim  - max dimension  allowed (tab must be allocated this size)
 *      n    - order of tableau to read
 *      name - name for error messages only
 *      tab  - (OUT) full 2d tableaux  (allocated by caller as dim*dim chars)
 *
 * Return value:
 *      Order of tableau read, or
 *      -n (n>1) tableau order n is too large to read
 */
int parse_tableau(FILE *fp, size_t dim, int n, const char *name, char *tab)
{
  char buf[MAX_LINE_LEN];
  char tabcode[2];
  char encval;
  int i,j,order;
  int toolarge = 0;
  
  order = n;
  if (n > (int)dim)
  {
    fprintf(stderr, "Tableau %s order %d is too large (max is %d)\n", 
            name, n, (int)dim);
    order = -n;
    toolarge = 1;
  }

  for (i = 0; i < n; i++)
  {
    fgets(buf, MAX_LINE_LEN, fp);
    if (!toolarge)
      for (j = 0; j <= i; j++)
      {
        tabcode[0] = buf[j*3];
        tabcode[1] = buf[j*3+1];
        if (i == j)
          encval = encode_ssetype(tabcode);
        else
          encval = encode_tabcode(tabcode);
        tab[INDEX2D(i,j,dim,dim)] = encval;
        tab[INDEX2D(j,i,dim,dim)] = encval;
      }
  }
  return order;
}


/*
 * Read numeric SSE distance matrix.
 *
 * Each row is a row of the distance matrix, lower triangle
 * only (since it is symmetric).
 * The diagonal entries are meaningless (self-distance)
 * and are included instead to specify the SSE type, with
 * the following codes:
 *
 * 0.000 beta strand
 * 1.000 alpha helix
 * 2.000 pi helix
 * 3.000 3_10 helix
 *
 * Each entry in matrix is in Angstroms format
 * F6.3 with a space between each on a line, and one line
 * per row of matrix.
 *
 * E.g.:
 *
 *   0.000 
 *   4.501  0.000 
 *  11.662 10.386  1.000 
 *  16.932 17.644  9.779  3.000 
 *  10.588 13.738 11.815 10.527  0.000 
 *  15.025 18.692 17.143 15.341  6.466  0.000 
 *  15.298 17.276 16.276 20.075 13.264 11.610  3.000 
 *   7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000 
 *
 * 
 * This subroutine is used to read the distance matrix from the 
 * database, where the identifier and order and tableau have already
 * been read by parse_tableaux().
 *
 * Parameters:
 *      fp - filehandle to read from
 *      dim  - max dimension  allowed (tab must be allocated this size)
 *      n  - order of distance matrix to read
 *      dmat  - (OUT) full 2d tableaux array allocated by caller as
 *                     dim*dim floats
 *      discard - if true, read data and don't store (used when too large)
 *
 * Return value:
 *      Order of tableau read, or negative on error
 *
 */
int parse_distmatrix(FILE *fp, size_t dim, int n, float *dmat, int discard)
{
  char buf[MAX_LINE_LEN];
  int i,j;
  float d;
  
  for (i = 0; i < n; i++)
  {
    fgets(buf, MAX_LINE_LEN, fp);
    if (!discard)
      for (j = 0; j <= i; j++)
      {
        d = strtof(&buf[j*7], NULL);
        dmat[INDEX2D(i,j,dim,dim)] = d;
        dmat[INDEX2D(j,i,dim,dim)] = d;
      }
  }
  return n;
}

/*
 * Read the entire tableaux+distmatrix database from file into memory,
 * with ones with order <= MAXDIM_GPU and ones with order > MAXDIM_GPU
 * in separte memory allocations.
 *
 * Parameters:
 *    fp - open FILE* of database file
 *    tableaux - (out) all tableaux, allocated here
 *    distmatrices - (out) all distmatrices, allocated here
 *    large_tableaux (out) tableux > MAXDIM_GPU, allocated here
 *    large_distmatrices (out) distmatrices > MAXDIM_GPU, allocated here
 *    orders - (out) order of <=MAXDIM_GPU tableaux&distmatrix, allocated here
 *    names - (out) name of each <=MAXDIM_GPU structure, allocated here
 *    large_orders - (out) orders of >MAXDIM_GPU tableaux+distmatrices
 *    large_names - (out) names or >MAXDIM_GPU structures
 *    num_large_read - (out) number of large structures read
 *
 * Return value: total number of tableaux+distmatrices read, -ve on error.
 *               this is large_num_read + num_read
 *
 */
int read_database(FILE *fp, char **tableaux, float **distmatrices, 
                  char **large_tableaux, float **large_distmatrices,
                  int **orders, char **names,
                  int **large_orders, char **large_names,
                  int *num_large_read)
{
  const int INITIAL_NUM = 1000; /* allocate this many to begin with */
  char *curtab, *large_curtab;
  float *curdmat, *large_curdmat;
  int *curord, *large_curord;
  char *curname, *large_curname;
  int num_read = 0, large_num_read = 0;
  char loctab[MAXDIM_GPU*MAXDIM_GPU];
  float locdmat[MAXDIM_GPU*MAXDIM_GPU];
  char large_loctab[MAXDIM*MAXDIM];
  float large_locdmat[MAXDIM*MAXDIM];
  int num_skipped = 0;
  int order,read_order;
  char name[LABELSIZE+1];
  
  *num_large_read = 0;

  if (!(*tableaux = (char *)malloc(INITIAL_NUM*MAXDIM_GPU*MAXDIM_GPU*sizeof(char)))) 
  {
    fprintf(stderr, "malloc tableaux failed\n");
    exit(1);
  }
  if (!(*distmatrices = (float *)malloc(INITIAL_NUM*MAXDIM_GPU*MAXDIM_GPU*sizeof(float))))
  {
    fprintf(stderr, "malloc distmatrices failed\n");
    exit(1);
  }
  if (!(*orders = (int *)malloc(INITIAL_NUM*sizeof(int))))
  {
    fprintf(stderr, "malloc orders failed\n");
    exit(1);
  }
  if (!(*names = (char *)malloc(INITIAL_NUM*(LABELSIZE+1)*sizeof(char))))
  {
    fprintf(stderr, "malloc names failed\n");
    exit(1);
  }
  curname = *names;
  curord = *orders;
  curtab = *tableaux;
  curdmat = *distmatrices;

  if (!(*large_tableaux = (char *)malloc(INITIAL_NUM*MAXDIM*MAXDIM*sizeof(char)))) 
  {
    fprintf(stderr, "malloc large_tableaux failed\n");
    exit(1);
  }
  if (!(*large_distmatrices = (float *)malloc(INITIAL_NUM*MAXDIM*MAXDIM*sizeof(float))))
  {
    fprintf(stderr, "malloc large_istmatrices failed\n");
    exit(1);
  }
  if (!(*large_orders = (int *)malloc(INITIAL_NUM*sizeof(int))))
  {
    fprintf(stderr, "malloc orders failed\n");
    exit(1);
  }
  if (!(*large_names = (char *)malloc(INITIAL_NUM*(LABELSIZE+1)*sizeof(char))))
  {
    fprintf(stderr, "malloc names failed\n");
    exit(1);
  }
  large_curname = *large_names;
  large_curord = *large_orders;
  large_curtab = *large_tableaux;
  large_curdmat = *large_distmatrices;

  while (!feof(fp))
  {
    if (fscanf(fp, "%8s %d\n", name, &order) != 2)
      break; /* normal for reading past last entry in db, detects EOF */
    
    if (order <= MAXDIM_GPU && num_read >= INITIAL_NUM)
    {
      /* past the initial allocation, realloc for this one */
      if (!(*tableaux = (char *)realloc(*tableaux, (num_read+1) * MAXDIM_GPU*MAXDIM_GPU*sizeof(char)))) 
      {
        fprintf(stderr, "realloc tableaux failed\n");
        exit(1);
      }
      curtab = *tableaux + num_read * MAXDIM_GPU*MAXDIM_GPU;
      if (!(*distmatrices = (float *)realloc(*distmatrices, (num_read+1) * MAXDIM_GPU*MAXDIM_GPU*sizeof(float))))
      {
        fprintf(stderr, "realloc distmatrices failed\n");
        exit(1);
      }
      curdmat = *distmatrices + num_read * MAXDIM_GPU*MAXDIM_GPU;
      if (!(*orders = (int *)realloc(*orders, (num_read+1)*sizeof(int))))
      {
        fprintf(stderr, "realloc orders failed\n");
        exit(1);
      }
      curord = *orders + num_read;
      if (!(*names = (char *)realloc(*names, (num_read+1)*(LABELSIZE+1)*sizeof(char))))
      {
        fprintf(stderr, "realloc names failed\n");
        exit(1);
      }
      curname = *names + num_read * (LABELSIZE+1);
    }
    else if (order > MAXDIM_GPU && large_num_read >= INITIAL_NUM)
    {
      /* past the initial allocation, realloc for this one */
      if (!(*large_tableaux = (char *)realloc(*large_tableaux, (large_num_read+1) * MAXDIM*MAXDIM*sizeof(char)))) 
      {
        fprintf(stderr, "realloc large_tableaux failed\n");
        exit(1);
      }
      large_curtab = *large_tableaux + large_num_read * MAXDIM*MAXDIM;
      if (!(*large_distmatrices = (float *)realloc(*large_distmatrices, (large_num_read+1) * MAXDIM*MAXDIM*sizeof(float))))
      {
        fprintf(stderr, "realloc large_distmatrices failed\n");
        exit(1);
      }
      large_curdmat = *large_distmatrices + large_num_read * MAXDIM*MAXDIM;
      if (!(*large_orders = (int *)realloc(*large_orders, (large_num_read+1)*sizeof(int))))
      {
        fprintf(stderr, "realloc large_orders failed\n");
        exit(1);
      }
      large_curord = *large_orders + large_num_read;
      if (!(*large_names = (char *)realloc(*large_names, (large_num_read+1)*(LABELSIZE+1)*sizeof(char))))
      {
        fprintf(stderr, "realloc large_names failed\n");
        exit(1);
      }
      large_curname = *large_names + large_num_read * (LABELSIZE+1);
    }

    
    if (order <= MAXDIM_GPU)
      read_order = parse_tableau(fp, MAXDIM_GPU, order, name, loctab);
    else
      read_order = parse_tableau(fp, MAXDIM, order, name, large_loctab);

    if (read_order < -1)  /* tableau too large*/
    {
      parse_distmatrix(fp, order <= MAXDIM_GPU ? MAXDIM_GPU : MAXDIM,
                       -(read_order), NULL, 1); /* read+discard the distmatrix */
      num_skipped++;
      continue;
    }

/* stop strncpy trunaction warning, this usage is actually what we want */
#pragma GCC diagnostic ignored "-Wstringop-truncation"

    if (order <= MAXDIM_GPU)
    {
      strncpy(curname, name, LABELSIZE);
      *curord = order;
      memcpy(curtab, &loctab, sizeof(loctab));
      parse_distmatrix(fp, MAXDIM_GPU, read_order, locdmat, 0);
      memcpy(curdmat, &locdmat, sizeof(locdmat));
      curname += (LABELSIZE+1);
      curord++;
      curtab += MAXDIM_GPU*MAXDIM_GPU;
      curdmat += MAXDIM_GPU*MAXDIM_GPU;
      num_read++;
    }
    else
    {
      strncpy(large_curname, name, LABELSIZE);
      *large_curord = order;
      memcpy(large_curtab, &large_loctab, sizeof(large_loctab));
      parse_distmatrix(fp, MAXDIM, read_order, large_locdmat, 0);
      memcpy(large_curdmat, &large_locdmat, sizeof(large_locdmat));
      large_curname += (LABELSIZE+1);
      large_curord++;
      large_curtab += MAXDIM*MAXDIM;
      large_curdmat += MAXDIM*MAXDIM;
      large_num_read++;
    }
  }

  if (num_skipped > 0) 
  {
    fprintf(stderr, "WARNING: skipped %d tableaux of order > %d\n",
            num_skipped, MAXDIM);
  }

  *num_large_read = large_num_read;
  return num_read + large_num_read;
}

/*
 * Read the entire list of query structures (tableaux+distmatrix)
 * from file into memory,
 *
 * Parameters:
 *    fp - open FILE* of database file
 *    tableaux - (out) all tableaux, allocated here
 *    distmatrices - (out) all distmatrices, allocated here
 *    orders - (out) order of tableaux&distmatrix, allocated here
 *    names - (out) name of each structure, allocated here
 *
 * Return value: total number of tableaux+distmatrices read, -ve on error.
 *
 */
int read_queries(FILE *fp, char **tableaux, float **distmatrices, 
                 int **orders, char **names)
{
  const int INITIAL_NUM = 1000; /* allocate this many to begin with */
  char *curtab;
  float *curdmat;
  int *curord;
  char *curname;
  int num_read = 0;
  char loctab[MAXDIM*MAXDIM];
  float locdmat[MAXDIM*MAXDIM];
  int num_skipped = 0;
  int order,read_order;
  char name[LABELSIZE+1];
  

  if (!(*tableaux = (char *)malloc(INITIAL_NUM*MAXDIM*MAXDIM*sizeof(char)))) 
  {
    fprintf(stderr, "malloc tableaux failed\n");
    exit(1);
  }
  if (!(*distmatrices = (float *)malloc(INITIAL_NUM*MAXDIM*MAXDIM*sizeof(float))))
  {
    fprintf(stderr, "malloc distmatrices failed\n");
    exit(1);
  }
  if (!(*orders = (int *)malloc(INITIAL_NUM*sizeof(int))))
  {
    fprintf(stderr, "malloc orders failed\n");
    exit(1);
  }
  if (!(*names = (char *)malloc(INITIAL_NUM*(LABELSIZE+1)*sizeof(char))))
  {
    fprintf(stderr, "malloc names failed\n");
    exit(1);
  }
  curname = *names;
  curord = *orders;
  curtab = *tableaux;
  curdmat = *distmatrices;


  while (!feof(fp))
  {
    if (fscanf(fp, "%8s %d\n", name, &order) != 2)
      break; /* normal for reading past last entry in db, detects EOF */
    
    if (order <= MAXDIM && num_read >= INITIAL_NUM)
    {
      /* past the initial allocation, realloc for this one */
      if (!(*tableaux = (char *)realloc(*tableaux, (num_read+1) * MAXDIM*MAXDIM*sizeof(char)))) 
      {
        fprintf(stderr, "realloc tableaux failed\n");
        exit(1);
      }
      curtab = *tableaux + num_read * MAXDIM*MAXDIM;
      if (!(*distmatrices = (float *)realloc(*distmatrices, (num_read+1) * MAXDIM*MAXDIM*sizeof(float))))
      {
        fprintf(stderr, "realloc distmatrices failed\n");
        exit(1);
      }
      curdmat = *distmatrices + num_read * MAXDIM*MAXDIM;
      if (!(*orders = (int *)realloc(*orders, (num_read+1)*sizeof(int))))
      {
        fprintf(stderr, "realloc orders failed\n");
        exit(1);
      }
      curord = *orders + num_read;
      if (!(*names = (char *)realloc(*names, (num_read+1)*(LABELSIZE+1)*sizeof(char))))
      {
        fprintf(stderr, "realloc names failed\n");
        exit(1);
      }
      curname = *names + num_read * (LABELSIZE+1);
    }

    
    read_order = parse_tableau(fp, MAXDIM, order, name, loctab);

    if (read_order < -1)  /* tableau too large*/
    {
      parse_distmatrix(fp,  MAXDIM,
                       -(read_order), NULL, 1); /* read+discard the distmatrix */
      num_skipped++;
      continue;
    }

/* stop strncpy trunaction warning, this usage is actually what we want */
#pragma GCC diagnostic ignored "-Wstringop-truncation"

    strncpy(curname, name, LABELSIZE);
    *curord = order;
    memcpy(curtab, &loctab, sizeof(loctab));
    parse_distmatrix(fp, MAXDIM, read_order, locdmat, 0);
    memcpy(curdmat, &locdmat, sizeof(locdmat));
    curname += (LABELSIZE+1);
    curord++;
    curtab += MAXDIM*MAXDIM;
    curdmat += MAXDIM*MAXDIM;
    num_read++;
  }

  if (num_skipped > 0) 
  {
    fprintf(stderr, "WARNING: skipped %d query tableaux of order > %d\n",
            num_skipped, MAXDIM);
  }

  return num_read;
}
