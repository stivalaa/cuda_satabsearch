/*****************************************************************************
 * 
 * File:    cudaSaTabsearch.cu
 * Author:  Alex Stivala
 * Created: January 2010
 *
 * $Id: cudaSaTabsearch.cu 3600 2010-05-03 07:17:03Z alexs $
 *
 * CUDA host code for simulated annealing tableau matching (discrete).
 * This is a CUDA implemenation of the FORTRAN subroutine TSAMTD.
 * Since the GPU has limited memory (and specifically, very limited
 * per block shared memory), we split the database into 'small' and
 * 'large' structures. The small ones can run on the GPU in shared memory,
 * the large ones cannot so we either have to not use shared memory
 * (OK, but a bit slower) or run them on the host.
 * When runnign on the host, we can simultaneously run the GPU and
 * host in separate threads. For multiple GPU cards, CUDA also requires
 * that there is a separate host thread for each GPU, so this program
 * is multithreaded: each thread is either for a separate GPU or for
 * running the same kernel (but compiled for host) on the host CPU.
 *
 * Usage: cudaSaTabsearch [-c] [-q dbfile] [-r restarts] < inputfile
 *
 * -c : run on host CPU not GPU card
 *
 * -q : query list mode: instead of reading query data on stdin
 *      just as in the original Fortran version tlocsd, a list
 *      of query sids to be read from the database is read on stdin (one per
 *      line),
 *      and db filenaame is specified on command
 *      line. In this mode options are assumed as LORDER=T, LTYPE=T,
 *      LSOLN=N. The output is still to stdout, but each query following
 *      immediately from the previous (can parse using the  header comment
 *      niformation lines as separators.
 *
 * -r restarts: number of restarts (iterations of cooling schedule).
 *              Should be a multiple of blocksize. Defaults to 128.
 *
 * The 'database' to search is an ASCII file of  tableaux
 * (Omega matrices) in format described in rdtabd.f.
 *
 * The results are printed to stdout as 
 *
 * name score
 *
 *
 * Both the name of the database file to read, and the actual
 * query tableau are read from stdin. 
 * The first line is the name
 * of the database file.
 * The second line is for options. There are currently 3 logical
 * options, for SSE type constraint (only allow SSEs of same type ot
 * match) and ordering constraint (disallow out of sequence order 
 * matches). The third is to output not just the scores but also solution
 * vector values.
 * They are single character logical values (T or F).
 * First is type, second is order, third is solution output,
 * separated by one space.
 *
 * The subsequent lines are a single tableau in the same format as
 * each tableau entry in the database i.e.:
 *
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
 * Following the tableau is the distance matrix.
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
 * 
 * E.g.:
 * 
 * /local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii
 *  T T F
 * D1UBIA_    8
 * e  
 * OT e  
 * LE RT xa 
 * PD OS RD xg 
 * RT LE RT LS e  
 * LE RD LE LS OT e  
 * RT LS LS RD PE OS xg 
 * PE RT LE RD OT PE RT e  
 *  0.000 
 *  4.501  0.000 
 *  1.662 10.386  1.000 
 * 16.932 17.644  9.779  3.000 
 * 10.588 13.738 11.815 10.527  0.000 
 * 15.025 18.692 17.143 15.341  6.466  0.000 
 * 15.298 17.276 16.276 20.075 13.264 11.610  3.000 
 *  7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000 
 *
 * 
 *
 *****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include <multithreading.h>
#include <cutil_inline.h>
#include "MersenneTwister.h"
#include "parsetableaux.h"
#include "cudaSaTabsearch_kernel.h"


/*****************************************************************************
 *
 * Type definitions
 *
 *****************************************************************************/

/* dbIndex_t is for the query list mode, an array of these gives for each
   query the index in the appropriate ('small' or 'large' according to the
   large flag) tableaux and distmatrix db arrays */
typedef struct dbIndex_s 
{
    bool large;  /* true if query is 'large' (>MAXDIM_GPU) structure */
    int  index;  /* index in tableaux and distmatrix db list, or 'large'
                    tableaux and distmatrix db list if large is true */
} dbIndex_t;

/* searchParams_t is a struct for parameter to tableau search functions
   dcelared as CUT_THREADROUTINE to be callable as threads */
typedef struct searchParams_s
{
    int ltype; int lorder; int lsoln; /* type,order,soln flags */
    int maxstart;           /* number of restarts */
    int maxdim;             /*dimension of tableaux, distmatrices here */
    int num_queries;        /* number of queries; 0 if not query list mode */
    int single_query_qid; /* if >=0, do only the one at this index */
    dbIndex_t *query_dbindex_list; /* if num_queries>0, the query db index */
    char qtab[MAXDIM*MAXDIM];     /* if num_queries==0, the query tableau */
    float qdmat[MAXDIM*MAXDIM];   /*                    the query distmatrix*/
    char qid[LABELSIZE+1];        /*                    the query identifier*/
    int qn;                       /*                    the query order */
    char *qssetypes;              /*                    the query SSE types*/

    int dbsize;             /* number of entries in the db */
    char *tableaux;         /* the tableaux database */
    float *distmatrices;    /* the distance matrices database */
    int   *orders;          /* orders of entries in db */
    char  *names;           /* names of entries in db */
    
} searchParams_t;


/*****************************************************************************
 *
 * Globals
 *
 *****************************************************************************/

static char dbfile[MAX_LINE_LEN];   /* database file name */
static bool use_gpu = true;   /* use the GPU */
static bool use_shared_memory = true; /* use GPU shared mem for db structs */
static char *tableaux, *large_tableaux; /* small and large tableaux */
static float *distmatrices, *large_distmatrices; /* same for dist.matrices*/
static int *orders, *large_orders; /* and for orders */
static char *names, *large_names;  /* and names */
static bool querydbmode = false;   /* use list of query ids in db */
static char *queryid_list = NULL;  /* this is the list of query ids */
static dbIndex_t *query_dbindex_list = NULL; /* and their indices in db */

static int maxstart = DEFAULT_MAXSTART; /* number of restarts */


/*****************************************************************************
 *
 * This part contains code adapted from CUDA SDK 2.3 
 *
 *****************************************************************************/
/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation and 
 * any modifications thereto.  Any use, reproduction, disclosure, or distribution 
 * of this software and related documentation without an express license 
 * agreement from NVIDIA Corporation is strictly prohibited.
 * 
 */




///////////////////////////////////////////////////////////////////////////////
// Common host and device function 
///////////////////////////////////////////////////////////////////////////////
//ceil(a / b)
extern "C" int iDivUp(int a, int b){
    return ((a % b) != 0) ? (a / b + 1) : (a / b);
}

//floor(a / b)
extern "C" int iDivDown(int a, int b){
    return a / b;
}

//Align a to nearest higher multiple of b
extern "C" int iAlignUp(int a, int b){
    return ((a % b) != 0) ?  (a - a % b + b) : a;
}

//Align a to nearest lower multiple of b
extern "C" int iAlignDown(int a, int b){
    return a - a % b;
}



///////////////////////////////////////////////////////////////////////////////
// Reference MT front-end
///////////////////////////////////////////////////////////////////////////////
extern "C" void initMTRef(const char *fname);
extern "C" void RandomHost(void);


///////////////////////////////////////////////////////////////////////////////
// Data configuration
///////////////////////////////////////////////////////////////////////////////
//const int    PATH_N = 24000000;
//const int N_PER_RNG = iAlignUp(iDivUp(PATH_N, MT_RNG_COUNT), 2);
//const int    RAND_N = MT_RNG_COUNT * N_PER_RNG;

const unsigned int SEED = 777;

int InitRandomNumberGenerator(void) {



    fprintf(stderr,"Loading CPU and GPU twisters configurations...\n");
/*
        const char *raw_path = cutFindFilePath("MersenneTwister.raw", argv[0]);
        const char *dat_path = cutFindFilePath("MersenneTwister.dat", argv[0]);
*/
        const char *raw_path = "data/MersenneTwister.raw";
        const char *dat_path = "data/MersenneTwister.dat"; 

        if (use_gpu) {
            loadMTGPU(dat_path);
            seedMTGPU(SEED);
        }

        initMTRef(raw_path);

    return 0;
}



/*****************************************************************************
 *
 * End of CUDA SDK 2.3 Mersenne Twister code
 *
 *****************************************************************************/

/*
 * tabsearch_host_thread - run the tableau search kernel on host CPU
 *
 * Started as a thread by cutStartThread in main
 *
 * Parameters:
 *   params - paramter block for thread. See comments on searchParams_t defn.
 *
 * Return value: None.
 *
 */
static CUT_THREADPROC tabsearch_host_thread(searchParams_t *params)
{
  /* extern declartions of host version of gpu constant memory */
  extern int c_qn_host;    // query structure size
  extern char c_qtab_host[MAXDIM*MAXDIM];  // query tableau
  extern float c_qdmat_host[MAXDIM*MAXDIM];  // query distance matrix
  extern char c_qssetypes_host[MAXDIM]; // main diagonal of c_qn


  unsigned int hTimer;
  double runtime;
  int *ssemaps;
  int i,j;
  char qid[LABELSIZE+1];
  int *scores;

  int query_count = (params->num_queries == 0 || params->single_query_qid >= 0
                     ? 1 : params->num_queries);

  cudaExtent tableaux_extent = {params->maxdim, params->maxdim,
                                params->dbsize};
  cudaPitchedPtr tableaux_pp = {params->tableaux, params->maxdim,
                                params->maxdim, params->dbsize};
  cudaExtent distmatrices_extent = {params->maxdim*sizeof(float), 
                                    params->maxdim,
                                    params->maxdim};
  cudaPitchedPtr distmatrices_pp = {params->distmatrices, 
                                    params->maxdim*sizeof(float),
                                    params->maxdim,
                                    params->maxdim};


  /* allocate space for output */
  if (!(scores = (int *)malloc(params->dbsize*sizeof(int))))
  {
    fprintf(stderr, "malloc scores failed\n");
    return;
  }
  if (!(ssemaps = (int *)malloc(params->dbsize*MAXDIM*sizeof(int))))
  {
    fprintf(stderr, "malloc ssemaps failed\n");
    return;
  }

  for (int qi = 0; qi < query_count; qi++)
  {
    if (params->query_dbindex_list)
    {
      dbIndex_t *dbindex_entry =  params->single_query_qid >= 0 ? 
        &params->query_dbindex_list[params->single_query_qid] :
        &params->query_dbindex_list[qi];
      int qdbi = dbindex_entry->index;

      if (dbindex_entry->large) /* query in 'large' struct db */
      {
        strncpy(qid, large_names+qdbi*(LABELSIZE+1), LABELSIZE);
        c_qn_host = large_orders[qdbi];
        memcpy(c_qtab_host, large_tableaux+qdbi*MAXDIM*MAXDIM,
               MAXDIM*MAXDIM*sizeof(char));
        memcpy(c_qdmat_host, large_distmatrices+qdbi*MAXDIM*MAXDIM,
               MAXDIM*MAXDIM*sizeof(float));
        /* NB the query in constant memory is MAXDIM not MAXDIM_GPU 
           since constant memory larger than shared memory. */
        // set the qssetypes vector as main diagonal of the query tableau
        for (i = 0; i < c_qn_host; i++)
          c_qssetypes_host[i] = (large_tableaux+qdbi*MAXDIM*MAXDIM)[INDEX2D(i,i,MAXDIM,MAXDIM)];
      }
      else /* query in 'small' struct db */
      {
        strncpy(qid, names+qdbi*(LABELSIZE+1), LABELSIZE);
        c_qn_host = orders[qdbi];
        
        /* NB the query in constant memory is MAXDIM not MAXDIM_GPU 
           since constant memory larger than shared memory.
           This means we need to reformat the matrices into the larger 
             size if they are in the smaller class */
        for (i = 0; i < orders[qdbi]; i++)
        {
          for (j = i + 1; j < orders[qdbi]; j++)
          {
            char tabcode = (tableaux+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,j,MAXDIM_GPU,MAXDIM_GPU)];
            c_qtab_host[INDEX2D(i,j,MAXDIM,MAXDIM)] = tabcode;
            c_qtab_host[INDEX2D(j,i,MAXDIM,MAXDIM)] = tabcode;
            float dist = (distmatrices+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,j,MAXDIM_GPU,MAXDIM_GPU)];
            c_qdmat_host[INDEX2D(i,j,MAXDIM,MAXDIM)] = dist;
            c_qdmat_host[INDEX2D(j,i,MAXDIM,MAXDIM)] = dist;
          }
        }
        // set the qssetypes vector as main diagonal of the query tableau
        for (i = 0; i < c_qn_host; i++)
          c_qssetypes_host[i] = (tableaux+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,i,MAXDIM_GPU,MAXDIM_GPU)];
      }
    }
    else
    {
      strncpy(qid, params->qid, LABELSIZE);
      c_qn_host = params->qn;
      memcpy(c_qtab_host, params->qtab, sizeof(c_qtab_host));
      memcpy(c_qdmat_host, params->qdmat, sizeof(c_qdmat_host));
      memcpy(c_qssetypes_host, params->qssetypes, sizeof(c_qssetypes_host));
    }
    
    printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
           params->ltype ? 'T' : 'F' , 
           params->lorder ? 'T' : 'F' , 
           params->lsoln ? 'T' : 'F');
    printf("# QUERY ID = %-8s\n", qid);
    printf("# DBFILE = %-80s\n", dbfile);
      
    fprintf(stderr, "Executing simulated annealing tableaux match kernel on host for query %s...\n", qid);
    cutilCheckError( cutCreateTimer(&hTimer) );
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );
    sa_tabsearch_host(params->dbsize,
                      params->lorder, 
                      params->lsoln,
                      params->maxstart,
                      tableaux_pp, tableaux_extent,
                      params->orders,
                      distmatrices_pp, distmatrices_extent,
                      scores,
                      ssemaps);
    cutilCheckError( cutStopTimer(hTimer) );
    runtime = cutGetTimerValue(hTimer);
    fprintf(stderr,  "host execution time %f ms\n", runtime);
    fprintf(stderr,  "%f million iterations/sec\n", (params->dbsize * (params->maxstart * MAXITER) / (runtime/1000)) / 1.0e6);
    
    for (i = 0; i < params->dbsize; i++)
    {
      printf("%-8s  %d\n", params->names+i*(LABELSIZE+1), scores[i]);
      if (params->lsoln)
        for (int k = 0; k < c_qn_host; k++)
          if (ssemaps[i*MAXDIM + k] >= 0)
            printf("%3d %3d\n", k+1, ssemaps[i*MAXDIM + k]+1);
    }
  }
  free(scores);
  if (params->lsoln)
    free(ssemaps);
}



/*
 * copyQueryToConstantMemory() - copy the query data to device constant memory
 *
 *
 * Parameters:
 *   qi - the query index of the query to copy. 
 *        Otherwise (query_dbinex_list is NULL), these used:
 *   qn -query order
 *   qtab - query tableau  (in/out: may be set here)
 *   qdmat - query distance matrix  (in/out: may be set here)
 *   qssetypes - query SSE types vector (in/out: may be set here)
 *   qid - query id (in/out: may be set here)
 *   c_qn_symbol - name of the c_qn constant ("q_qn" or "c_qn_noshared")
 *   c_qtab_symbol - name fo the c_qtab constant
 *   c_qdmat_synmbol - name of the c_qdmat constant
 *   c_qssetypes_symbol - name of the c_qssetypes constant
 *   
 *
 * Uses the global variables query_dbindex_list, tableaux, etc.
 *
 * Return value: None.
 *
 */
static void copyQueryToConstantMemory(int qi, 
                                      int qn, char *qtab, float *qdmat,
                                      char *qssetypes, char *qid,
                                      const char *c_qn_symbol,
                                      const char *c_qtab_symbol,
                                      const char *c_qdmat_symbol,
                                      const char *c_qssetypes_symbol)
{
  unsigned int hTimer;
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );
  if (query_dbindex_list)
  {
    int qdbi = query_dbindex_list[qi].index;
    if (query_dbindex_list[qi].large)
    {
      strncpy(qid, large_names+qdbi*(LABELSIZE+1), LABELSIZE);
      // set the qssetypes vector as main diagonal of the query tableau
      for (int i = 0; i < large_orders[qdbi]; i++)
        qssetypes[i] = (large_tableaux+qdbi*MAXDIM*MAXDIM)[INDEX2D(i,i,MAXDIM,MAXDIM)];
      /* copy query structure to constant memory on device */
      cutilSafeCall( cudaMemcpyToSymbol(c_qn_symbol, &large_orders[qdbi], sizeof(int)) );
      /* NB the query in constant memory is MAXDIM not MAXDIM_GPU 
         since constant memory larger than shared memory. */
      cutilSafeCall( cudaMemcpyToSymbol(c_qtab_symbol, large_tableaux+qdbi*MAXDIM*MAXDIM, MAXDIM*MAXDIM*sizeof(char)) );
      cutilSafeCall( cudaMemcpyToSymbol(c_qdmat_symbol, large_distmatrices+qdbi*MAXDIM*MAXDIM, MAXDIM*MAXDIM*sizeof(float)) );
      cutilSafeCall( cudaMemcpyToSymbol(c_qssetypes_symbol, qssetypes, MAXDIM*sizeof(char)) );
    }
    else /* query is in the 'small' structure dbase */
    {
      strncpy(qid, names+qdbi*(LABELSIZE+1), LABELSIZE);
      // set the qssetypes vector as main diagonal of the query tableau
      for (int i = 0; i < orders[qdbi]; i++)
        qssetypes[i] = (tableaux+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,i,MAXDIM_GPU,MAXDIM_GPU)];
      /* copy query structure to constant memory on device */
      cutilSafeCall( cudaMemcpyToSymbol(c_qn_symbol, &orders[qdbi], sizeof(int)) );
      /* NB the query in constant memory is MAXDIM not MAXDIM_GPU 
         since constant memory larger than shared memory.
         This means we need to reformat the matrices into the larger 
         size if they are in the smaller class */
      for (int i = 0; i < orders[qdbi]; i++)
      {
        for (int j = i + 1; j < orders[qdbi]; j++)
        {
          char tabcode = (tableaux+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,j,MAXDIM_GPU,MAXDIM_GPU)];
          qtab[INDEX2D(i,j,MAXDIM,MAXDIM)] = tabcode;
          qtab[INDEX2D(j,i,MAXDIM,MAXDIM)] = tabcode;
          float dist = (distmatrices+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,j,MAXDIM_GPU,MAXDIM_GPU)];
          qdmat[INDEX2D(i,j,MAXDIM,MAXDIM)] = dist;
          qdmat[INDEX2D(j,i,MAXDIM,MAXDIM)] = dist;
        }
      }
      cutilSafeCall( cudaMemcpyToSymbol(c_qtab_symbol, qtab, MAXDIM*MAXDIM*sizeof(char)) );
      cutilSafeCall( cudaMemcpyToSymbol(c_qdmat_symbol, qdmat, MAXDIM*MAXDIM*sizeof(float)) );
      
      cutilSafeCall( cudaMemcpyToSymbol(c_qssetypes_symbol, qssetypes, MAXDIM*sizeof(char)) );
    }
  }
  else // single query mode - copy to constant memory
  {
    cutilSafeCall( cudaMemcpyToSymbol(c_qn_symbol, &qn, sizeof(qn)) );
    cutilSafeCall( cudaMemcpyToSymbol(c_qtab_symbol, qtab, MAXDIM*MAXDIM*sizeof(char)) );
    cutilSafeCall( cudaMemcpyToSymbol(c_qdmat_symbol, qdmat, MAXDIM*MAXDIM*sizeof(float)) );
    cutilSafeCall( cudaMemcpyToSymbol(c_qssetypes_symbol, qssetypes, MAXDIM*sizeof(char)) );
  }
  cutilCheckError( cutStopTimer(hTimer) );
  float qtime = cutGetTimerValue(hTimer);
  fprintf(stderr, "Copying query to constant memory (%s) took %f ms\n", 
          c_qtab_symbol,
          qtime);
}


static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s [-c] [-q dbfile]\n", progname);
  fprintf(stderr, "  -c : run on host CPU not GPU card\n");
  fprintf(stderr, "  -q dbfile : database is read from dbfile, list of query\n"
          "              ids is read from stdin\n");
  fprintf(stderr, "   -r restarts : number of restarts. Default %d\n",
          DEFAULT_MAXSTART);
  exit(1);
}


int main(int argc, char *argv[])
{
  CUTThread threadID[MAX_THREADS];
  int num_threads = 0;
  int exit_status = 0;
  char buf[MAX_LINE_LEN];
  char qtab[MAXDIM*MAXDIM];
  float qdmat[MAXDIM*MAXDIM];
  int qn;
  char qid[LABELSIZE+1];
  int ltype=0,lorder=0,lsoln=0;
  char cltype,clorder,clsoln;
  FILE *dbfp;
  unsigned int hTimer;
  int total_dbsize, large_dbsize, gpu_dbsize;
  double dbtime,runtime;
  cudaPitchedPtr d_tableaux;
  cudaPitchedPtr d_distmatrices;
  int *d_orders;
  int *scores = NULL;
  int *ssemaps = NULL;
  int *d_scores;
  int *d_ssemaps;
  cudaError_t cuda_errcode;
  int i,j;
  char qssetypes[MAXDIM];
  int c;
  char *queryptr = NULL;
  int num_queries = 0;
  int large_query_count = 0;

  while ((c = getopt(argc, argv, "cq:r:")) != -1)
  {
    switch (c)
    {
      case 'c':
        use_gpu = false;
        break;

      case 'q':
        querydbmode = true;
        strncpy(dbfile, optarg, sizeof(dbfile)-1);
        break;

      case 'r':
        maxstart = atoi(optarg);
        break;
        
      default:
        usage(argv[0]);
        break;
    }
  }

  if (querydbmode)
  {
    cltype = 'T'; ltype = 1;
    clorder = 'T'; lorder = 1;
    clsoln = 'F'; lsoln = 0;
    if (!(queryid_list = (char *)malloc(LABELSIZE+1)))
    {
      fprintf(stderr, "malloc queryid_list failed\n");
      exit(1);
    }
    queryptr = queryid_list;
    while (!feof(stdin))
    {
      if (num_queries > 0)
      {
        if ((!(queryid_list = (char *)realloc(queryid_list, (num_queries+1)*(LABELSIZE+1)))))
        {
          fprintf(stderr, "realloc queryid_list failed\n");
          exit(1);
        }
      }
      if (!fgets(buf, MAX_LINE_LEN, stdin))
        break;
      strncpy(queryptr, buf, LABELSIZE);
      queryptr[LABELSIZE-1] = '\0';
      if (queryptr[strlen(queryptr)-1] == '\n')
        queryptr[strlen(queryptr)-1] = '\0';
      queryptr += (LABELSIZE+1);
      num_queries++;
    }
  }
  else
  {
    if (fscanf(stdin, "%s\n", dbfile) != 1)
    {
      fprintf(stderr, "ERROR reading dbfilename from stdin\n");
      exit(1);
    }
    if (fscanf(stdin, "%c %c %c\n", &cltype, &clorder, &clsoln) != 3)
    {
      fprintf(stderr, "ERROR reading options from stdin\n");
      exit(1);
    }
    if (cltype == 'T')
      ltype = 1;
    if (clorder == 'T')
      lorder = 1;
    if (clsoln == 'T')
      lsoln = 1;
    
    if (fscanf(stdin, "%8s %d\n", qid, &qn) != 2)
    {
      fprintf(stderr, "ERROR parsing query tableau header from stdin\n");
      exit(1);
    }
    if (parse_tableau(stdin, MAXDIM, qn, qtab) < 0)
    {
      fprintf(stderr, "ERROR parsing query tableau from stdin\n");
      exit(1);
    }
    if (parse_distmatrix(stdin, MAXDIM, qn, qdmat, 0) < 0)
    {
      fprintf(stderr, "ERROR parsing query distance matrix from stdin\n");
      exit(1);
    }
  }

  if (!ltype)
  {
    fprintf(stderr, "WARNING: LTYPE is always set to T\n");
    ltype = 1; cltype = 'T';
  }
  
  if (!(dbfp = fopen(dbfile, "r")))
  {
    fprintf(stderr, "ERROR opening db file %s\n", dbfile);
    exit(1);
  }

  fprintf(stderr, "Loading database...\n");
  cutilCheckError( cutCreateTimer(&hTimer) );
  cutilCheckError( cutResetTimer(hTimer) );
  cutilCheckError( cutStartTimer(hTimer) );
  total_dbsize = read_database(dbfp, &tableaux, &distmatrices, 
                               &large_tableaux, &large_distmatrices,
                               &orders, &names,
                               &large_orders, &large_names,
                               &large_dbsize);
  if (total_dbsize < 0)
  {
    fprintf(stderr, "ERROR loading database\n");
    exit(1);
  }
  gpu_dbsize = total_dbsize - large_dbsize;
  cutilCheckError( cutStopTimer(hTimer) );
  dbtime = cutGetTimerValue(hTimer);
  fprintf(stderr, "Loaded %d db entries (%d order > %d) in %f ms\n", 
          total_dbsize, large_dbsize, MAXDIM_GPU, dbtime);
          

  if (querydbmode)
  {
    /* Convert the list of query sids to list of indices in db for later
       rapid lookup.
       TODO: we should build a hash table rather than this highly 
       inefficient linear search for each query id, but it's only
       done once and db not that big... 
    */
    fprintf(stderr, "Building query index list...\n");
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );
    if (!(query_dbindex_list = (dbIndex_t *)malloc(num_queries*sizeof(dbIndex_t))))
    {
      fprintf(stderr, "malloc query_dbindex_list failed\n");
      exit(1);
    }
    for (i = 0; i < num_queries; i++)
    {
/*      fprintf(stderr, "zzz %s\n", queryid_list+i*(LABELSIZE+1)); */
      bool found = false;
      for (j = 0; j < gpu_dbsize; j++) /* search 'small' structure dbase */
      {
        if (!strcasecmp(queryid_list+i*(LABELSIZE+1),names+j*(LABELSIZE+1)))
        {
          query_dbindex_list[i].large = false;
          query_dbindex_list[i].index = j;
          found = true;
          break;
        }
      }
      if (!found)
      {
        for (j = 0; j < large_dbsize; j++) /* search 'large' structure dbase*/
        {
          if (!strcasecmp(queryid_list + i*(LABELSIZE+1),
                           large_names + j*(LABELSIZE+1)))
          {
            query_dbindex_list[i].large = true;
            query_dbindex_list[i].index = j;
            large_query_count++;
            found = true;
            break;
          }
        }
      }
      if (!found)
      {
        fprintf(stderr, "ERROR: query %s not found\n", queryid_list+i*(LABELSIZE+1));
        exit(1);
      }
    }
    cutilCheckError( cutStopTimer(hTimer) );
    fprintf(stderr, "Built query index (%d queries (%d large)) in %f ms\n",
            num_queries, large_query_count, cutGetTimerValue(hTimer));
  }
  else
  {
    num_queries = 0; 
    query_dbindex_list = NULL;
    // set the qssetypes vector as main diagonal of the query tableau
    for (i = 0; i < qn; i++)
      qssetypes[i] = qtab[INDEX2D(i,i,MAXDIM,MAXDIM)];
  }
    
  /* TODO allow multiple GPUs (need one thread for each) */

  if (use_gpu)
  {
/*
    int devnum = cutGetMaxGflopsDeviceId();
    fprintf(stderr, "using max gflops device %d: ", devnum);
*/
    /* If there is a compute capability 2 device ("Fermi"
       architecture) (or higher) then use that, and do NOT use shared
       memory as it is faster to just rely on the new "NVIDIA Parallel
       DataCache (TM)" -- just use global memory for all (small and large)
       structures
    */

    int devnum, deviceCount, gflops,max_gflops=0, sel_devnum;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
      fprintf(stderr, "There is no device supporting CUDA.\n");
      exit(1);
    }
    fprintf(stderr, "found %d CUDA devices\n", deviceCount);
    for (devnum = 0; devnum < deviceCount; devnum++)
    {  
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, devnum);
      if (deviceProp.major >= 2)
      {
        fprintf(stderr,
          "found Fermi architecture (compute capability %d.%d) device %d: %s\n",
                deviceProp.major, deviceProp.minor, devnum, deviceProp.name);
        sel_devnum = devnum;
        use_shared_memory = false;
        break;
      }
      else
      {
        gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
        fprintf(stderr, "device %d: %s\n", devnum,
                deviceProp.name);
        if (gflops > max_gflops)
        {
          max_gflops = gflops;
          sel_devnum = devnum;
          use_shared_memory = true;
        }
      }
    }
    
    fprintf(stderr, "using device %d: ", sel_devnum);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, sel_devnum);
    fprintf(stderr, "%s\n", deviceProp.name);
    cudaSetDevice( sel_devnum );
  }

  InitRandomNumberGenerator();

  fprintf(stderr, "maxstart = %d\n", maxstart);

  if (use_gpu)
  {
    /* setup execution configuration parameters */
    /* TODO optimize for different architectures (automatically) */
    dim3 dimGrid(128);          // blocks
    dim3 dimBlock(128);         // threads per block


    fprintf(stderr, "Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);

    if (dimBlock.x * dimGrid.x > MT_RNG_COUNT)
    {
      fprintf(stderr, "ERROR: can only have a max of %d threads in configuration\n", MT_RNG_COUNT);
      exit(1);
    }

    fprintf(stderr, "using shared memory for small db structs: %s\n",
            use_shared_memory ? "YES" : "NO");

    /* first do the 'small' db structures on GPU (with shared memory) */

    fprintf(stderr, "Copying database to device...\n");
    cutilCheckError( cutResetTimer(hTimer) );
    cutilCheckError( cutStartTimer(hTimer) );

    cudaExtent tableaux_extent = make_cudaExtent(MAXDIM_GPU, MAXDIM_GPU,
                                                 gpu_dbsize);
    cutilSafeCall( cudaMalloc3D(&d_tableaux, tableaux_extent) );
    fprintf(stderr, "d_tableaux.pitch == %u xsize == %u ysize == %u\n", d_tableaux.pitch, d_tableaux.xsize, d_tableaux.ysize);

    cudaExtent distmatrices_extent = make_cudaExtent(MAXDIM_GPU*sizeof(float), MAXDIM_GPU, gpu_dbsize);
    cutilSafeCall( cudaMalloc3D(&d_distmatrices, distmatrices_extent) );
    fprintf(stderr, "d_distmatrices.pitch == %u xsize == %u ysize == %u\n", d_distmatrices.pitch, d_distmatrices.xsize, d_distmatrices.ysize);

    cutilSafeCall( cudaMalloc((void **)&d_orders, gpu_dbsize*sizeof(int)) );


    cudaMemcpy3DParms copyParams = { 0 };
    // srcPtr is tricky: need to give pitch of row, #elements in row,
    // then height, omitting 3rd dimension (doesn't seem to be documented)
    // (I found this info on 28/1/2010 at 
    // http://sites.google.com/site/cudaiap2009/cookbook-1).
    // Note pitch of row on host is just MAXDIM_GPU, we don't need padding here
    copyParams.srcPtr = make_cudaPitchedPtr((void*)tableaux, MAXDIM_GPU, MAXDIM_GPU, MAXDIM_GPU);

    fprintf(stderr, "srcPtr.pitch == %u\n", copyParams.srcPtr.pitch);

    copyParams.dstPtr = d_tableaux;
    copyParams.extent = tableaux_extent;
    copyParams.kind = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&copyParams) );

    
    cudaMemcpy3DParms copyParams2 = { 0 };
    copyParams2.srcPtr = make_cudaPitchedPtr((void*)distmatrices,
                                             MAXDIM_GPU*sizeof(float),
                                             MAXDIM_GPU, MAXDIM_GPU);
    fprintf(stderr, "distmatrices srcPtr.pitch == %u\n", copyParams2.srcPtr.pitch);
    copyParams2.dstPtr = d_distmatrices;
    copyParams2.extent = distmatrices_extent;
    copyParams2.kind = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&copyParams2) );

    cutilSafeCall( cudaMemcpy(d_orders, orders, gpu_dbsize*sizeof(int),
                              cudaMemcpyHostToDevice) );

    cutilCheckError( cutStopTimer(hTimer) );
    dbtime = cutGetTimerValue(hTimer);
    fprintf(stderr, "Copied %d entries to GPU in %f ms\n", gpu_dbsize, dbtime);

    /* allocate space for output */
    cutilSafeCall( cudaMalloc((void **)&d_scores, gpu_dbsize*sizeof(int)));
    if (!(scores = (int *)malloc(gpu_dbsize*sizeof(int))))
    {
      fprintf(stderr, "malloc scores failed\n");
      goto bye;
    }
    if (lsoln)
    {
      cutilSafeCall( cudaMalloc((void **)&d_ssemaps, gpu_dbsize*MAXDIM*sizeof(int)));
      if (!(ssemaps = (int *)malloc(gpu_dbsize*MAXDIM*sizeof(int))))
      {
        fprintf(stderr, "malloc ssemaps failed\n");
        goto bye;
      }
    }
      
    int query_count = (num_queries == 0 ? 1 : num_queries);
    for (int qi = 0; qi < query_count; qi++)
    {
      if (use_shared_memory)
        copyQueryToConstantMemory(qi, qn, qtab, qdmat, qssetypes, qid,
                                  "c_qn", "c_qtab", "c_qdmat", "c_qssetypes");
      else
        copyQueryToConstantMemory(qi, qn, qtab, qdmat, qssetypes, qid,
                                  "c_qn_noshared_small",
                                  "c_qtab_noshared_small", 
                                  "c_qdmat_noshared_small", 
                                  "c_qssetypes_noshared_small");
        

      printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
             cltype, clorder, clsoln);
      printf("# QUERY ID = %-8s\n", qid);
      printf("# DBFILE = %-80s\n", dbfile);


      /* launch thread to do large db structs on host */
      searchParams_t host_params;
      host_params.ltype = ltype;
      host_params.lorder = lorder;
      host_params.lsoln = lsoln;
      host_params.maxstart = maxstart;
      host_params.num_queries = num_queries;
      host_params.query_dbindex_list = query_dbindex_list;
      host_params.single_query_qid = qi; 
      memcpy(host_params.qtab, qtab, sizeof(qtab));
      memcpy(host_params.qdmat, qdmat, sizeof(qdmat));
      memcpy(host_params.qid, qid, sizeof(qid));
      host_params.qn = qn;
      host_params.qssetypes = qssetypes;
      host_params.maxdim = MAXDIM;
      host_params.dbsize = large_dbsize;
      host_params.tableaux = large_tableaux;
      host_params.distmatrices = large_distmatrices;
      host_params.orders = large_orders;
      host_params.names = large_names;

//XXX      threadID[num_threads++] = cutStartThread((CUT_THREADROUTINE)tabsearch_host_thread, &host_params);


      fprintf(stderr, "Executing simulated annealing tableaux match kernel (%sshared memory) on GPU for qid %s...\n", use_shared_memory ? " " : "no ", qid);
      cutilSafeCall( cudaThreadSynchronize() );
      cutilCheckError( cutResetTimer(hTimer) );
      cutilCheckError( cutStartTimer(hTimer) );
      if (use_shared_memory)
        sa_tabsearch_gpu<<<dimGrid,dimBlock>>>(gpu_dbsize,
                                               lorder, 
                                               lsoln,
                                               maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                               d_ssemaps);
      else
        sa_tabsearch_gpu_noshared_small<<<dimGrid,dimBlock>>>(gpu_dbsize,
                                               lorder, 
                                               lsoln,
                                               maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                               d_ssemaps);

      cuda_errcode = cudaGetLastError();
      if (cuda_errcode != cudaSuccess)
      {
        fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cuda_errcode));
        exit_status = 1;
        goto bye;
      }

      cutilSafeCall( cudaThreadSynchronize() );
      cutilCheckError( cutStopTimer(hTimer) );
      runtime = cutGetTimerValue(hTimer);
      fprintf(stderr,  "GPU execution time %f ms\n", runtime);
      fprintf(stderr,  "%f million iterations/sec\n", ((float)gpu_dbsize * ((float)maxstart * (float)MAXITER) / (runtime/1000)) / 1.0e6);
      
      /* Get results from device */
      cutilSafeCall( cudaMemcpy(scores, d_scores, gpu_dbsize*sizeof(int), 
                                cudaMemcpyDeviceToHost) );
      if (lsoln)
        cutilSafeCall( cudaMemcpy(ssemaps, d_ssemaps, 
                                  gpu_dbsize*MAXDIM*sizeof(int),
                                  cudaMemcpyDeviceToHost) );

      /* Wait for host thread */
//XXX      cutWaitForThreads(threadID, num_threads);
//XXX      --num_threads;



      /* TODO we could reduce wasted time waiting by running all host 
         (large db) queries in the one thread instead of matching up
         with GPU query in this loop (actuall, more like the other way
         around usually, the GPU ends up idle while host is still runnign since
         the latter is so much slower even though it has very few
         db entries unlike GPU) */

      for (i = 0; i < gpu_dbsize; i++)
      {
        printf("%-8s  %d\n", names+i*(LABELSIZE+1), scores[i]);
        if (lsoln)
          for (int k = 0; k < qn; k++)
            if (ssemaps[i*MAXDIM + k] >= 0)
              printf("%3d %3d\n", k+1, ssemaps[i*MAXDIM + k]+1);
      }
    }

    cutilSafeCall( cudaFree(d_tableaux.ptr) );
    cutilSafeCall( cudaFree(d_distmatrices.ptr) );
    cutilSafeCall( cudaFree(d_orders) );
    cutilSafeCall( cudaFree(d_scores) );
    free(scores); scores = NULL;
    if (lsoln)
    {
      cutilSafeCall( cudaFree(d_ssemaps) );
      free(ssemaps); ssemaps = NULL;
    }

    /* now do the 'large' db structures on GPU (not using shared memory) */
    if (large_dbsize > 0)
    {
      fprintf(stderr, "Copying large structure database to device...\n");
      cutilCheckError( cutResetTimer(hTimer) );
      cutilCheckError( cutStartTimer(hTimer) );

      tableaux_extent = make_cudaExtent(MAXDIM, MAXDIM, large_dbsize);
      cutilSafeCall( cudaMalloc3D(&d_tableaux, tableaux_extent) );
      fprintf(stderr, "d_tableaux.pitch == %u xsize == %u ysize == %u\n", d_tableaux.pitch, d_tableaux.xsize, d_tableaux.ysize);

      distmatrices_extent = make_cudaExtent(MAXDIM*sizeof(float), MAXDIM, large_dbsize);
      cutilSafeCall( cudaMalloc3D(&d_distmatrices, distmatrices_extent) );
      fprintf(stderr, "d_distmatrices.pitch == %u xsize == %u ysize == %u\n", d_distmatrices.pitch, d_distmatrices.xsize, d_distmatrices.ysize);

      cutilSafeCall( cudaMalloc((void **)&d_orders, large_dbsize*sizeof(int)) );


      cudaMemcpy3DParms copyParamsl = { 0 };
      // srcPtr is tricky: need to give pitch of row, #elements in row,
      // then height, omitting 3rd dimension (doesn't seem to be documented)
      // (I found this info on 28/1/2010 at 
      // http://sites.google.com/site/cudaiap2009/cookbook-1).
      // Note pitch of row on host is just MAXDIM_GPU, we don't need padding here
      copyParamsl.srcPtr = make_cudaPitchedPtr((void*)large_tableaux, MAXDIM, MAXDIM, MAXDIM);

      fprintf(stderr, "srcPtr.pitch == %u\n", copyParamsl.srcPtr.pitch);

      copyParamsl.dstPtr = d_tableaux;
      copyParamsl.extent = tableaux_extent;
      copyParamsl.kind = cudaMemcpyHostToDevice;
      cutilSafeCall( cudaMemcpy3D(&copyParamsl) );


      cudaMemcpy3DParms copyParams2l = { 0 };
      copyParams2l.srcPtr = make_cudaPitchedPtr((void*)large_distmatrices,
                                                MAXDIM*sizeof(float),
                                                MAXDIM, MAXDIM);
      fprintf(stderr, "distmatrices srcPtr.pitch == %u\n", copyParams2l.srcPtr.pitch);
      copyParams2l.dstPtr = d_distmatrices;
      copyParams2l.extent = distmatrices_extent;
      copyParams2l.kind = cudaMemcpyHostToDevice;
      cutilSafeCall( cudaMemcpy3D(&copyParams2l) );

      cutilSafeCall( cudaMemcpy(d_orders, large_orders, large_dbsize*sizeof(int),
                                cudaMemcpyHostToDevice) );

      cutilCheckError( cutStopTimer(hTimer) );
      dbtime = cutGetTimerValue(hTimer);
      fprintf(stderr, "Copied %d large entries to GPU in %f ms\n", large_dbsize, dbtime);

      /* allocate space for output */
      cutilSafeCall( cudaMalloc((void **)&d_scores, large_dbsize*sizeof(int)));
      if (!(scores = (int *)malloc(large_dbsize*sizeof(int))))
      {
        fprintf(stderr, "malloc scores failed\n");
        goto bye;
      }
      if (lsoln)
      {
        cutilSafeCall( cudaMalloc((void **)&d_ssemaps, large_dbsize*MAXDIM*sizeof(int)));
        if (!(ssemaps = (int *)malloc(large_dbsize*MAXDIM*sizeof(int))))
        {
          fprintf(stderr, "malloc ssemaps failed\n");
          goto bye;
        }
      }

      for (int qi = 0; qi < query_count; qi++)
      {
        copyQueryToConstantMemory(qi, qn, qtab, qdmat, qssetypes, qid,
                                  "c_qn_noshared", "c_qtab_noshared", 
                                  "c_qdmat_noshared", "c_qssetypes_noshared");

        printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
               cltype, clorder, clsoln);
        printf("# QUERY ID = %-8s\n", qid);
        printf("# DBFILE = %-80s\n", dbfile);



        fprintf(stderr, "Executing simulated annealing tableaux match kernel (no shared memory) on GPU for qid %s...\n",qid);
        cutilSafeCall( cudaThreadSynchronize() );
        cutilCheckError( cutResetTimer(hTimer) );
        cutilCheckError( cutStartTimer(hTimer) );
        sa_tabsearch_gpu_noshared<<<dimGrid,dimBlock>>>(large_dbsize,
                                                        lorder, 
                                                        lsoln,
                                                        maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                               d_ssemaps);
        cuda_errcode = cudaGetLastError();
        if (cuda_errcode != cudaSuccess)
        {
          fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cuda_errcode));
          exit_status = 1;
          goto bye;
        }

        cutilSafeCall( cudaThreadSynchronize() );
        cutilCheckError( cutStopTimer(hTimer) );
        runtime = cutGetTimerValue(hTimer);
        fprintf(stderr,  "GPU (no shared memory) execution time %f ms\n", runtime);
        fprintf(stderr,  "%f million iterations/sec\n", ((float)large_dbsize * ((float)maxstart * (float)MAXITER) / (runtime/1000)) / 1.0e6);

        /* Get results from device */
        cutilSafeCall( cudaMemcpy(scores, d_scores, large_dbsize*sizeof(int), 
                                  cudaMemcpyDeviceToHost) );
        if (lsoln)
          cutilSafeCall( cudaMemcpy(ssemaps, d_ssemaps, 
                                    large_dbsize * MAXDIM * sizeof(int),
                                    cudaMemcpyDeviceToHost) );

        for (i = 0; i < large_dbsize; i++)
        {
          printf("%-8s  %d\n", large_names+i*(LABELSIZE+1), scores[i]);
          if (lsoln)
            for (int k = 0; k < qn; k++)
              if (ssemaps[i*MAXDIM + k] >= 0)
                printf("%3d %3d\n", k+1, ssemaps[i*MAXDIM + k]+1);
        }
      }
    }
  }
  else
  {
    /* running on host CPU */

    searchParams_t host_params;
    host_params.ltype = ltype;
    host_params.lorder = lorder;
    host_params.lsoln = lsoln;
    host_params.maxstart = maxstart;
    host_params.num_queries = num_queries;
    host_params.single_query_qid = -1;
    host_params.query_dbindex_list = query_dbindex_list;

    memcpy(host_params.qtab, qtab, sizeof(qtab));
    memcpy(host_params.qdmat, qdmat, sizeof(qdmat));
    memcpy(host_params.qid, qid, sizeof(qid));
    host_params.qn = qn;
    host_params.qssetypes = qssetypes;
    host_params.maxdim = MAXDIM_GPU;
    host_params.dbsize = gpu_dbsize;

    /* first do small structure db */
    host_params.tableaux = tableaux;
    host_params.distmatrices = distmatrices;
    host_params.orders = orders;
    host_params.names = names;
    
    tabsearch_host_thread(&host_params);

    /* then large structure db */
    if (large_dbsize > 0)
    {
      host_params.maxdim = MAXDIM;
      host_params.dbsize = large_dbsize;
      host_params.tableaux = large_tableaux;
      host_params.distmatrices = large_distmatrices;
      host_params.orders = large_orders;
      host_params.names = large_names;
      tabsearch_host_thread(&host_params);
    }
  }

bye:
  /* cleanup and exit */
  free(tableaux);
  free(distmatrices);
  free(orders);
  free(names);
  free(scores);
  free(large_tableaux);
  free(large_distmatrices);
  free(large_names);
  free(large_orders);
  if (lsoln)
    free(ssemaps);
  cutilCheckError( cutDeleteTimer( hTimer) );
  if (use_gpu)
  {
    if (large_dbsize > 0)
    {
      cutilSafeCall( cudaFree(d_tableaux.ptr) );
      cutilSafeCall( cudaFree(d_distmatrices.ptr) );
      cutilSafeCall( cudaFree(d_orders) );
      cutilSafeCall( cudaFree(d_scores) );
      if (lsoln)
        cutilSafeCall( cudaFree(d_ssemaps) );
    }
    cudaThreadExit();
  }
  exit(exit_status);
}
