/*****************************************************************************
 * 
 * File:    cudaSaTabsearch.cu
 * Author:  Alex Stivala
 * Created: January 2010
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
 * name rawscore norm2score z-score p-value
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
 * The subsequent lines are multiple query structures (tableau and
 * distance matrix for each), separated by a blank line. I.e. the
 * same format as the database.
 *
 * The tableau is in the same format as
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
 * D1AE6H1   13
 * e  
 * PD e  
 * OT OS e  
 * LS LS RD xg 
 * LE LE RT RT e  
 * RT RT LE LE OT e  
 * RT RT LE LE OT PE e  
 * LE LE RT RT PE OT PE e  
 * RT OT PE LS OT LE PE OT e  
 * PE PE OT LS PE RT OT PE OT e  
 * OT OT PE RD OS LE PE OT PE OT xg 
 * OT RT PE RD OT PE PE OT PE OT PD e  
 * PD PE OT LS LE RT RT LE OT PE OT OT e  
 *  0.000 
 * 19.130  0.000 
 *  8.850 13.371  0.000 
 * 14.608 29.221 15.945  3.000 
 * 12.469 19.135 11.231 16.008  0.000 
 * 18.479 21.128 16.959 21.982  6.730  0.000 
 * 16.153 22.704 13.140 13.210  6.909 10.709  0.000 
 * 20.850 24.610 16.558 16.527 10.946 12.552  4.935  0.000 
 * 15.604 18.394  8.791 14.366 11.402 16.188  8.316  9.609  0.000 
 * 13.949 13.565  5.751 17.301 10.771 15.314 10.725 12.661  4.876  0.000 
 * 24.234 12.620 19.140 31.786 17.166 14.790 20.733 21.202 20.224 16.665  3.000 
 *  9.731 17.355  9.936 16.942  3.841  8.797 10.327 14.566 13.021 11.226 17.023  0.000 
 * 16.856  5.985 12.706 27.454 15.011 16.146 19.829 22.156 17.744 13.079  9.541 12.996  0.000 
 *
 *****************************************************************************/

#define CUDASATABSEARCH_MAIN 1

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include <multithreading.h>
#include <helper_cuda.h>
#include <helper_timer.h>
#include <curand_kernel.h>
#include "parsetableaux.h"
#include "cudaSaTabsearch_kernel.h"
#include "cudaGetDeviceConstantAddresses.h"
#include "gumbelstats.h"



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

/* queryData_t is a struct containing the data for a single query structure */
typedef struct queryData_s {
    char qtab[MAXDIM*MAXDIM];     /* the query tableau */
    float qdmat[MAXDIM*MAXDIM];   /* the query distmatrix*/
    char qid[LABELSIZE+1];        /* the query identifier*/
    int qn;                       /* the query order */
    char *qssetypes;              /* the query SSE types*/
} queryData_t;

/* searchParams_t is a struct for parameter to tableau search functions
   dcelared as CUT_THREADROUTINE to be callable as threads */
typedef struct searchParams_s
{
    int ltype; int lorder; int lsoln; /* type,order,soln flags */
    unsigned long maxstart;           /* number of restarts */
    unsigned long maxdim;             /*dimension of tableaux, distmatrices here */
    int num_queries;        /* number of queries */
    int single_query_qid; /* if >=0, do only the one at this index */
    dbIndex_t *query_dbindex_list; /* the query db index (query list mode) */
                                   /* OR query data (notq query list mode) : */
    char  *q_tableaux;                     /* query tableaux */
    float *q_distmatrices;                 /* query distance matrices */
    int   *q_orders;                       /* sizes of query tableaux */
    char  *q_names;                        /* names of queries */

    unsigned long dbsize;   /* number of entries in the db */
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

static char *query_tableaux;       /* query tableaux */
static float *query_distmatrices;; /* query dist.matrices*/
static int *query_orders;          /* query orders */
static char *query_names;          /* query names */



/*
 * init_rng()
 *
 * Initialize CURAND pseudrandom number generator
 * See CUDA Toolkit 4.1 CURAND Guide (p.21)
 *
 * Parameters:
 *    state - CURAND state for random number generation
 *
 */
__global__ void init_rng(curandState *state)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  /* give each therad same seed, different sequence number, no offset */
  curand_init(1234, tid, 0, &state[tid]);
}




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


  StopWatchInterface *hTimer = NULL;
  double runtime;
  int *ssemaps;
  int i,j;
  char qid[LABELSIZE+1];
  int *scores;
  double norm2score,zscore,pvalue;
  char qssetypes[MAXDIM];

  int query_count = (params->query_dbindex_list && params->single_query_qid >= 0
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
      strncpy(qid, params->q_names+qi*(LABELSIZE+1), LABELSIZE);
      c_qn_host = params->q_orders[qi];
      memcpy(c_qtab_host, params->q_tableaux+qi*MAXDIM*MAXDIM, sizeof(c_qtab_host));
      memcpy(c_qdmat_host, params->q_distmatrices+qi*MAXDIM*MAXDIM, sizeof(c_qdmat_host));
      // set the qssetypes vector as main diagonal of the query tableau
      for (i = 0; i < params->q_orders[qi]; i++)
        qssetypes[i] = (params->q_tableaux+qi*MAXDIM*MAXDIM)[INDEX2D(i,i,MAXDIM,MAXDIM)];
      memcpy(c_qssetypes_host, qssetypes, sizeof(c_qssetypes_host));
    }
    
    printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
           params->ltype ? 'T' : 'F' , 
           params->lorder ? 'T' : 'F' , 
           params->lsoln ? 'T' : 'F');
    printf("# QUERY ID = %-8s\n", qid);
    printf("# DBFILE = %-80s\n", dbfile);
      
    fprintf(stderr, "Executing simulated annealing tableaux match kernel on host for query %s...\n", qid);
    sdkCreateTimer(&hTimer) ;
    sdkResetTimer(&hTimer) ;
    sdkStartTimer(&hTimer) ;
    int state = 0; /*unused*/
    sa_tabsearch_host(params->dbsize,
                      params->lorder, 
                      params->lsoln,
                      params->maxstart,
                      tableaux_pp, tableaux_extent,
                      params->orders,
                      distmatrices_pp, distmatrices_extent,
                      scores,
                      ssemaps,
                      &state);
    sdkStopTimer(&hTimer);
    runtime = sdkGetTimerValue(&hTimer);
    fprintf(stderr,  "host execution time %f ms\n", runtime);
    fprintf(stderr,  "%f million iterations/sec\n", (params->dbsize * (params->maxstart * MAXITER) / (runtime/1000)) / 1.0e6);
    
    for (i = 0; i < params->dbsize; i++)
    {
/*      printf("%-8s  %d\n", params->names+i*(LABELSIZE+1), scores[i]); */
      norm2score = norm2(scores[i], c_qn_host, params->orders[i]);
      zscore = z_gumbel(norm2score, gumbel_a, gumbel_b);
      pvalue = pv_gumbel(zscore);
      printf("%-8s %d %g %g %g\n", params->names+i*(LABELSIZE+1),
             scores[i], norm2score, zscore, pvalue);
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
 *   c_qn_addr - address of c_qn device constant (q_qn or c_qn_noshared)
 *   c_qtab_addr - address of c_qtab device constant
 *   c_qdmat_addr - address of c_qdmat device constant
 *   c_qssetypes_addr - address c_qssetypes device constant
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
                                      int *c_qn_addr,
                                      char *c_qtab_addr,
                                      float *c_qdmat_addr,
                                      char *c_qssetypes_addr)
{
  StopWatchInterface *hTimer = NULL;
  sdkCreateTimer(&hTimer) ;
  sdkResetTimer(&hTimer) ;
  sdkStartTimer(&hTimer) ;
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
      checkCudaErrors( cudaMemcpy(c_qn_addr, &large_orders[qdbi], sizeof(int), cudaMemcpyHostToDevice) );
      /* NB the query in constant memory is MAXDIM not MAXDIM_GPU 
         since constant memory larger than shared memory. */
      checkCudaErrors( cudaMemcpy(c_qtab_addr, large_tableaux+qdbi*MAXDIM*MAXDIM, MAXDIM*MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
      checkCudaErrors( cudaMemcpy(c_qdmat_addr, large_distmatrices+qdbi*MAXDIM*MAXDIM, MAXDIM*MAXDIM*sizeof(float), cudaMemcpyHostToDevice) );
      checkCudaErrors( cudaMemcpy(c_qssetypes_addr, qssetypes, MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
    }
    else /* query is in the 'small' structure dbase */
    {
      strncpy(qid, names+qdbi*(LABELSIZE+1), LABELSIZE);
      // set the qssetypes vector as main diagonal of the query tableau
      for (int i = 0; i < orders[qdbi]; i++)
        qssetypes[i] = (tableaux+qdbi*MAXDIM_GPU*MAXDIM_GPU)[INDEX2D(i,i,MAXDIM_GPU,MAXDIM_GPU)];
      /* copy query structure to constant memory on device */
      checkCudaErrors( cudaMemcpy(c_qn_addr, &orders[qdbi], sizeof(int), cudaMemcpyHostToDevice) );
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
      checkCudaErrors( cudaMemcpy(c_qtab_addr, qtab, MAXDIM*MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
      checkCudaErrors( cudaMemcpy(c_qdmat_addr, qdmat, MAXDIM*MAXDIM*sizeof(float), cudaMemcpyHostToDevice) );
      
      checkCudaErrors( cudaMemcpy(c_qssetypes_addr, qssetypes, MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
    }
  }
  else // single query mode - copy to constant memory
  {
    fprintf(stderr, "XXX c_qn_addr = %p , qn = %d\n", c_qn_addr,  qn);
    checkCudaErrors( cudaMemcpy(c_qn_addr, &qn, sizeof(qn), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(c_qtab_addr, qtab, MAXDIM*MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(c_qdmat_addr, qdmat, MAXDIM*MAXDIM*sizeof(float), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy(c_qssetypes_addr, qssetypes, MAXDIM*sizeof(char), cudaMemcpyHostToDevice) );
  }

  sdkStopTimer(&hTimer) ;
  float qtime = sdkGetTimerValue(&hTimer);
  fprintf(stderr, "Copying query to constant memory took %f ms\n", 
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
//  CUTThread threadID[MAX_THREADS];
//  int num_threads = 0;
  int exit_status = 0;
  char buf[MAX_LINE_LEN];
  char qtab[MAXDIM*MAXDIM];
  float qdmat[MAXDIM*MAXDIM];
  int qn;
  char qid[LABELSIZE+1];
  int ltype=0,lorder=0,lsoln=0;
  char cltype,clorder,clsoln;
  FILE *dbfp;
  StopWatchInterface *hTimer = NULL;
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
  double norm2score, zscore, pvalue;

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
      //fprintf(stderr, "XXX num_queries = %d\n", num_queries);
      if (num_queries > 0)
      {
        if ((!(queryid_list = (char *)realloc(queryid_list, (num_queries+1)*(LABELSIZE+1)))))
        {
          fprintf(stderr, "realloc queryid_list failed\n");
          exit(1);
        }
        queryptr = queryid_list + num_queries*(LABELSIZE+1);
      }
      if (!fgets(buf, MAX_LINE_LEN, stdin))
        break;
      strncpy(queryptr, buf, LABELSIZE);
      queryptr[LABELSIZE-1] = '\0';
      if (queryptr[strlen(queryptr)-1] == '\n')
        queryptr[strlen(queryptr)-1] = '\0';
      //fprintf(stderr, "XXX queryptr = '%s'\n", queryptr);
      queryptr += (LABELSIZE+1);
      num_queries++;
    }
  }
  else /* not querydbmode: read query tableaux+distmatrices on stdin */
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
    
    num_queries = read_queries(stdin, &query_tableaux, &query_distmatrices,
                               &query_orders, &query_names);
    if (num_queries < 0) {
      fprintf(stderr, "ERROR loading query structures from stdin\n");
      exit(1);
    } else if (num_queries == 0) {
      fprintf(stderr, "ERROR: no query structures found on stdin\n");
      exit(1);
    }
    fprintf(stderr, "Read %d query structures\n", num_queries);
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
  sdkCreateTimer(&hTimer) ;
  sdkResetTimer(&hTimer) ;
  sdkStartTimer(&hTimer) ;
  total_dbsize = read_database(dbfp, &tableaux, &distmatrices, 
                               &large_tableaux, &large_distmatrices,
                               &orders, &names,
                               &large_orders, &large_names,
                               &large_dbsize);
  fclose(dbfp);
  if (total_dbsize < 0)
  {
    fprintf(stderr, "ERROR loading database\n");
    exit(1);
  }
  gpu_dbsize = total_dbsize - large_dbsize;
  sdkStopTimer(&hTimer) ;
  dbtime = sdkGetTimerValue(&hTimer);
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
    sdkResetTimer(&hTimer) ;
    sdkStartTimer(&hTimer) ;
    if (!(query_dbindex_list = (dbIndex_t *)malloc(num_queries*sizeof(dbIndex_t))))
    {
      fprintf(stderr, "malloc query_dbindex_list failed\n");
      exit(1);
    }
    for (i = 0; i < num_queries; i++)
    {
      //fprintf(stderr, "zzz %s\n", queryid_list+i*(LABELSIZE+1));  //XXX
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
    sdkStopTimer(&hTimer);
    fprintf(stderr, "Built query index (%d queries (%d large)) in %f ms\n",
            num_queries, large_query_count, sdkGetTimerValue(&hTimer));
  }
  else
  {
    query_dbindex_list = NULL;
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
    checkCudaErrors( cudaGetDeviceCount(&deviceCount) );
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
          "found modern architecture (compute capability %d.%d) device %d: %s\n",
                deviceProp.major, deviceProp.minor, devnum, deviceProp.name);
        sel_devnum = devnum;
        use_shared_memory = true; 
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
    fprintf(stderr, "totalGlobalMem              = %g GB\n"
                    "sharedMemPerBlock           = %g KB\n"
                    "warpSize                    = %d\n"
                    "maxThreadsPerBlock          = %d\n"
                    "clockRate                   = %g MHz\n"
                    "totalConstMem               = %g KB\n"
                    "multiProcessorCount         = %d\n"
                    "maxThreadsPerMultiProcessor = %d\n"
                    "sharedMeMPerMultiprocessor  = %d KB\n"
                    "maxBlocksPerMultiProcessor  = %d\n",
            (double)deviceProp.totalGlobalMem / (1024*1024*1024),
            (double)deviceProp.sharedMemPerBlock / 1024,
            deviceProp.warpSize,
            deviceProp.maxThreadsPerBlock,
            (double)deviceProp.clockRate / 1000,
            (double)deviceProp.totalConstMem / 1024,
            deviceProp.multiProcessorCount,
            deviceProp.maxThreadsPerMultiProcessor,
            (double)deviceProp.sharedMemPerMultiprocessor / 1024,
            deviceProp.maxBlocksPerMultiProcessor);

    cudaSetDevice( sel_devnum );
  }


  fprintf(stderr, "maxstart = %d\n", maxstart);

  srand48(1234);

  if (use_gpu)
  {
    /* setup execution configuration parameters */
    /* TODO optimize for different architectures (automatically) */
    const int blocks = 128;
    const int NUM_THREADS = 128;
    dim3 dimGrid(blocks);          // blocks
    dim3 dimBlock(NUM_THREADS);         // threads per block


    fprintf(stderr, "Execution configuration: Grid = (%d,%d,%d) Block = (%d,%d,%d)\n", dimGrid.x,dimGrid.y,dimGrid.z, dimBlock.x,dimBlock.y,dimBlock.z);



    fprintf(stderr, "using shared memory for small db structs: %s\n",
            use_shared_memory ? "YES" : "NO");

    /* first do the 'small' db structures on GPU (with shared memory) */

    fprintf(stderr, "Copying database to device...\n");
    sdkResetTimer(&hTimer) ;
    sdkStartTimer(&hTimer) ;

    curandState *devStates;
  /* allocate space on device for random number generator state */
    int rc;
  if ((rc = cudaMalloc((void **)&devStates, 
                       blocks*NUM_THREADS*sizeof(curandState))) != cudaSuccess)
  {
    fprintf(stderr, "cudaMalloc devStates failed %d\n", rc);
    exit(1);
  }
  
  /* initialize device random number generator */
  sdkStartTimer(&hTimer) ;
  init_rng<<<dimGrid, dimBlock>>>(devStates);
  if ((rc = cudaGetLastError()) != cudaSuccess)
  {
    fprintf(stderr, "init_rng kernel error %d\n", rc);
  }
  cudaDeviceSynchronize();
  if ((rc = cudaGetLastError()) != cudaSuccess)
  {
    fprintf(stderr, "init_rng sync error %d\n", rc);
  }
  sdkStopTimer(&hTimer) ;
  fprintf(stderr, "Initialized device RNG with %d states (%d KB) in %f ms\n",
          blocks*NUM_THREADS, 
          blocks*NUM_THREADS*sizeof(curandState)/1024,
          sdkGetTimerValue(&hTimer));

    cudaExtent tableaux_extent = make_cudaExtent(MAXDIM_GPU, MAXDIM_GPU,
                                                 gpu_dbsize);
    checkCudaErrors( cudaMalloc3D(&d_tableaux, tableaux_extent) );
    fprintf(stderr, "d_tableaux.pitch == %u xsize == %u ysize == %u\n", d_tableaux.pitch, d_tableaux.xsize, d_tableaux.ysize);

    cudaExtent distmatrices_extent = make_cudaExtent(MAXDIM_GPU*sizeof(float), MAXDIM_GPU, gpu_dbsize);
    checkCudaErrors( cudaMalloc3D(&d_distmatrices, distmatrices_extent) );
    fprintf(stderr, "d_distmatrices.pitch == %u xsize == %u ysize == %u\n", d_distmatrices.pitch, d_distmatrices.xsize, d_distmatrices.ysize);

    checkCudaErrors( cudaMalloc((void **)&d_orders, gpu_dbsize*sizeof(int)) );


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
    checkCudaErrors( cudaMemcpy3D(&copyParams) );

    
    cudaMemcpy3DParms copyParams2 = { 0 };
    copyParams2.srcPtr = make_cudaPitchedPtr((void*)distmatrices,
                                             MAXDIM_GPU*sizeof(float),
                                             MAXDIM_GPU, MAXDIM_GPU);
    fprintf(stderr, "distmatrices srcPtr.pitch == %u\n", copyParams2.srcPtr.pitch);
    copyParams2.dstPtr = d_distmatrices;
    copyParams2.extent = distmatrices_extent;
    copyParams2.kind = cudaMemcpyHostToDevice;
    checkCudaErrors( cudaMemcpy3D(&copyParams2) );

    checkCudaErrors( cudaMemcpy(d_orders, orders, gpu_dbsize*sizeof(int),
                              cudaMemcpyHostToDevice) );

    sdkStopTimer(&hTimer) ;
    dbtime = sdkGetTimerValue(&hTimer);
    fprintf(stderr, "Copied %d entries to GPU in %f ms\n", gpu_dbsize, dbtime);

    /* allocate space for output */
    checkCudaErrors( cudaMalloc((void **)&d_scores, gpu_dbsize*sizeof(int)));
    if (!(scores = (int *)malloc(gpu_dbsize*sizeof(int))))
    {
      fprintf(stderr, "malloc scores failed\n");
      goto bye;
    }
    if (lsoln)
    {
      checkCudaErrors( cudaMalloc((void **)&d_ssemaps, gpu_dbsize*MAXDIM*sizeof(int)));
      if (!(ssemaps = (int *)malloc(gpu_dbsize*MAXDIM*sizeof(int))))
      {
        fprintf(stderr, "malloc ssemaps failed\n");
        goto bye;
      }
    }

    const_addr_t const_addr;
    for (int qi = 0; qi < num_queries; qi++)
    {
      if (!querydbmode) {
      strncpy(qid,  query_names+qi*(LABELSIZE+1), LABELSIZE);
      qn = query_orders[qi];
      // set the qssetypes vector as main diagonal of the query tableau
      for (i = 0; i < query_orders[qi]; i++)
        qssetypes[i] = (query_tableaux+qi*MAXDIM*MAXDIM)[INDEX2D(i,i,MAXDIM,MAXDIM)];

      } else {
        qn = orders[qi];
      }
      if (use_shared_memory) {
        get_device_constant_addresses(&const_addr);
        copyQueryToConstantMemory(qi, qn,
                                  querydbmode ? qtab : query_tableaux+qi*MAXDIM*MAXDIM,
                                  querydbmode ? qdmat : query_distmatrices+qi*MAXDIM*MAXDIM,
                                  qssetypes,
                                  qid,
                                  const_addr.c_qn_addr, const_addr.c_qtab_addr,
                                  const_addr.c_qdmat_addr,
                                  const_addr.c_qssetypes_addr);

//      checkCudaErrors( cudaMemcpy(const_addr.c_qn_addr, &qn, sizeof(qn), cudaMemcpyHostToDevice) ); fprintf(stderr,"qn=%d\n",qn); //XXX

      }
      else {
        get_device_constant_addresses_noshared_small(&const_addr);    
        copyQueryToConstantMemory(qi, qn,
                                  querydbmode ? qtab : query_tableaux+qi*MAXDIM*MAXDIM,
                                  querydbmode ? qdmat : query_distmatrices+qi*MAXDIM*MAXDIM,
                                  qssetypes,
                                  qid,
                                  const_addr.c_qn_noshared_small_addr,
                                  const_addr.c_qtab_noshared_small_addr, 
                                  const_addr.c_qdmat_noshared_small_addr, 
                                  const_addr.c_qssetypes_noshared_small_addr);
      }
        

      printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
             cltype, clorder, clsoln);
      printf("# QUERY ID = %-8s\n", qid);
      printf("# DBFILE = %-80s\n", dbfile);


      fprintf(stderr, "Executing simulated annealing tableaux match kernel (%sshared memory) on GPU for qid %s...\n", use_shared_memory ? " " : "no ", qid);
      checkCudaErrors( cudaDeviceSynchronize() );

      sdkResetTimer(&hTimer) ;
      sdkStartTimer(&hTimer) ;
      if (use_shared_memory)
      {
        int xxx_qn=-1; checkCudaErrors( cudaMemcpy(&xxx_qn, const_addr.c_qn_addr, sizeof(qn), cudaMemcpyDeviceToHost) ); fprintf(stderr,"xxx_qn=%d\n",xxx_qn); //XXX

        sa_tabsearch_gpu<<<dimGrid,dimBlock>>>(gpu_dbsize,
                                               lorder, 
                                               lsoln,
                                               maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                               d_ssemaps,
                                               devStates);
      }
      else 
      {
        sa_tabsearch_gpu_noshared_small<<<dimGrid,dimBlock>>>(gpu_dbsize,
                                               lorder, 
                                               lsoln,
                                               maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                                              d_ssemaps,
                                                              devStates);
      }

      cuda_errcode = cudaGetLastError();
      if (cuda_errcode != cudaSuccess)
      {
        fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cuda_errcode));
        exit_status = 1;
        goto bye;
      }

      checkCudaErrors( cudaDeviceSynchronize() );
      sdkStopTimer(&hTimer) ;
      runtime = sdkGetTimerValue(&hTimer);
      fprintf(stderr,  "GPU execution time %f ms\n", runtime);
      fprintf(stderr,  "%f million iterations/sec\n", ((float)gpu_dbsize * ((float)maxstart * (float)MAXITER) / (runtime/1000)) / 1.0e6);
      
      /* Get results from device */
      checkCudaErrors( cudaMemcpy(scores, d_scores, gpu_dbsize*sizeof(int), 
                                cudaMemcpyDeviceToHost) );
      if (lsoln)
        checkCudaErrors( cudaMemcpy(ssemaps, d_ssemaps, 
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
/*        printf("%-8s  %d\n", names+i*(LABELSIZE+1), scores[i]); */
        norm2score = norm2(scores[i], qn, orders[i]);
        zscore = z_gumbel(norm2score, gumbel_a, gumbel_b);
        pvalue = pv_gumbel(zscore);
        printf("%-8s %d %g %g %g\n", 
               names+i*(LABELSIZE+1), scores[i], norm2score, zscore, pvalue);
        if (lsoln)
          for (int k = 0; k < qn; k++)
            if (ssemaps[i*MAXDIM + k] >= 0)
              printf("%3d %3d\n", k+1, ssemaps[i*MAXDIM + k]+1);
      }
    }

    checkCudaErrors( cudaFree(d_tableaux.ptr) );
    checkCudaErrors( cudaFree(d_distmatrices.ptr) );
    checkCudaErrors( cudaFree(d_orders) );
    checkCudaErrors( cudaFree(d_scores) );
    free(scores); scores = NULL;
    if (lsoln)
    {
      checkCudaErrors( cudaFree(d_ssemaps) );
      free(ssemaps); ssemaps = NULL;
    }

    /* now do the 'large' db structures on GPU (not using shared memory) */
    if (large_dbsize > 0)
    {
      fprintf(stderr, "Copying large structure database to device...\n");
      sdkResetTimer(&hTimer) ;
      sdkStartTimer(&hTimer) ;

      tableaux_extent = make_cudaExtent(MAXDIM, MAXDIM, large_dbsize);
      checkCudaErrors( cudaMalloc3D(&d_tableaux, tableaux_extent) );
      fprintf(stderr, "d_tableaux.pitch == %u xsize == %u ysize == %u\n", d_tableaux.pitch, d_tableaux.xsize, d_tableaux.ysize);

      distmatrices_extent = make_cudaExtent(MAXDIM*sizeof(float), MAXDIM, large_dbsize);
      checkCudaErrors( cudaMalloc3D(&d_distmatrices, distmatrices_extent) );
      fprintf(stderr, "d_distmatrices.pitch == %u xsize == %u ysize == %u\n", d_distmatrices.pitch, d_distmatrices.xsize, d_distmatrices.ysize);

      checkCudaErrors( cudaMalloc((void **)&d_orders, large_dbsize*sizeof(int)) );


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
      checkCudaErrors( cudaMemcpy3D(&copyParamsl) );


      cudaMemcpy3DParms copyParams2l = { 0 };
      copyParams2l.srcPtr = make_cudaPitchedPtr((void*)large_distmatrices,
                                                MAXDIM*sizeof(float),
                                                MAXDIM, MAXDIM);
      fprintf(stderr, "distmatrices srcPtr.pitch == %u\n", copyParams2l.srcPtr.pitch);
      copyParams2l.dstPtr = d_distmatrices;
      copyParams2l.extent = distmatrices_extent;
      copyParams2l.kind = cudaMemcpyHostToDevice;
      checkCudaErrors( cudaMemcpy3D(&copyParams2l) );

      checkCudaErrors( cudaMemcpy(d_orders, large_orders, large_dbsize*sizeof(int),
                                cudaMemcpyHostToDevice) );

      sdkStopTimer(&hTimer) ;
      dbtime = sdkGetTimerValue(&hTimer);
      fprintf(stderr, "Copied %d large entries to GPU in %f ms\n", large_dbsize, dbtime);

      /* allocate space for output */
      checkCudaErrors( cudaMalloc((void **)&d_scores, large_dbsize*sizeof(int)));
      if (!(scores = (int *)malloc(large_dbsize*sizeof(int))))
      {
        fprintf(stderr, "malloc scores failed\n");
        goto bye;
      }
      if (lsoln)
      {
        checkCudaErrors( cudaMalloc((void **)&d_ssemaps, large_dbsize*MAXDIM*sizeof(int)));
        if (!(ssemaps = (int *)malloc(large_dbsize*MAXDIM*sizeof(int))))
        {
          fprintf(stderr, "malloc ssemaps failed\n");
          goto bye;
        }
      }

      for (int qi = 0; qi < num_queries; qi++)
      {
        get_device_constant_addresses_noshared(&const_addr);    
        copyQueryToConstantMemory(qi, qn,
                                  querydbmode ? qtab : query_tableaux+qi*MAXDIM*MAXDIM,
                                  querydbmode ? qdmat : query_distmatrices+qi*MAXDIM*MAXDIM,
                                  qssetypes,
                                  qid,
                                  const_addr.c_qn_noshared_addr,
                                  const_addr.c_qtab_noshared_addr, 
                                  const_addr.c_qdmat_noshared_addr,
                                  const_addr.c_qssetypes_noshared_addr);


        printf("# cudaSaTabsearch LTYPE = %c LORDER = %c LSOLN = %c\n",
               cltype, clorder, clsoln);
        printf("# QUERY ID = %-8s\n", qid);
        printf("# DBFILE = %-80s\n", dbfile);



        fprintf(stderr, "Executing simulated annealing tableaux match kernel (no shared memory) on GPU for qid %s...\n",qid);
        checkCudaErrors( cudaDeviceSynchronize() );
        sdkResetTimer(&hTimer) ;
        sdkStartTimer(&hTimer) ;

        int xxx_qn_noshared=-1; checkCudaErrors( cudaMemcpy(&xxx_qn_noshared, const_addr.c_qn_noshared_addr, sizeof(qn), cudaMemcpyDeviceToHost) ); fprintf(stderr,"xxx_qn_noshared=%d\n",xxx_qn_noshared); //XXX
        sa_tabsearch_gpu_noshared<<<dimGrid,dimBlock>>>(large_dbsize,
                                                        lorder, 
                                                        lsoln,
                                                        maxstart,
                                               d_tableaux, tableaux_extent,
                                               d_orders,
                                               d_distmatrices, distmatrices_extent,
                                               d_scores,
                                                        d_ssemaps,
                                                        devStates);
        cuda_errcode = cudaGetLastError();
        if (cuda_errcode != cudaSuccess)
        {
          fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cuda_errcode));
          exit_status = 1;
          goto bye;
        }

        checkCudaErrors( cudaDeviceSynchronize() );
        sdkStopTimer(&hTimer) ;
        runtime = sdkGetTimerValue(&hTimer);
        fprintf(stderr,  "GPU (no shared memory) execution time %f ms\n", runtime);
        fprintf(stderr,  "%f million iterations/sec\n", ((float)large_dbsize * ((float)maxstart * (float)MAXITER) / (runtime/1000)) / 1.0e6);

        /* Get results from device */
        checkCudaErrors( cudaMemcpy(scores, d_scores, large_dbsize*sizeof(int), 
                                    cudaMemcpyDeviceToHost) );
        if (lsoln)
          checkCudaErrors( cudaMemcpy(ssemaps, d_ssemaps, 
                                    large_dbsize * MAXDIM * sizeof(int),
                                    cudaMemcpyDeviceToHost) );

        for (i = 0; i < large_dbsize; i++)
        {
/*          printf("%-8s  %d\n", large_names+i*(LABELSIZE+1), scores[i]); */
          norm2score = norm2(scores[i], qn, large_orders[i]);
          zscore = z_gumbel(norm2score, gumbel_a, gumbel_b);
          pvalue = pv_gumbel(zscore);
          printf("%-8s %d %g %g  %g\n", 
                 large_names+i*(LABELSIZE+1), scores[i], norm2score,
                 zscore, pvalue);
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
    host_params.q_tableaux = query_tableaux;
    host_params.q_distmatrices = query_distmatrices;
    host_params.q_orders = query_orders;
    host_params.q_names = query_names;
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
  sdkDeleteTimer( &hTimer);
  if (use_gpu)
  {
    if (large_dbsize > 0)
    {
      checkCudaErrors( cudaFree(d_tableaux.ptr) );
      checkCudaErrors( cudaFree(d_distmatrices.ptr) );
      checkCudaErrors( cudaFree(d_orders) );
      checkCudaErrors( cudaFree(d_scores) );
      if (lsoln)
        checkCudaErrors( cudaFree(d_ssemaps) );
    }
    cudaDeviceReset(); /* replaces deprecated cudaThreadExit() */
  }
  exit(exit_status);
}
