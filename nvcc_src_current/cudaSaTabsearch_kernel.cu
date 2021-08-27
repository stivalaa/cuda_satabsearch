/*****************************************************************************
 * 
 * File:    cudaSaTabsearch_kernel.cu
 * Author:  Alex Stivala
 * Created: January 2010
 *
 * $Id: cudaSaTabsearch_kernel.cu 4755 2013-11-20 03:46:36Z astivala $
 *
 * CUDA kernel for simulated annealing tableau matching (discrete).
 * This is a CUDA implemenation of the FORTRAN subroutine TSAMTD.
 * 
 * if CUDA preprocessor symbol is defined, this is the CUDA kernel version.
 * __DEVICE_EMULATION__ may also be defined for this case (nvcc -deviceemu)
 * in which case device emulation mode is being used
 *
 * Otherwise (CUDA symbol not defined), this builds a host (single threaded)
 * version.
 *
 * if DEBUG (in which case either __DEVICE_EMULATION__ must be defined,
 * or CUDA must not be defined),  is defined then verbose stderr output
 * is generated, and various assertions and checks are compiled in.
 *
 * If CUDA5_DEBUG is defined, then debugging via printf() from device
 * is used.
 *
 * If USE_SHARED_MEMORY is defined, then each block copies the tableau
 * and distance matrix it is operating on frmot he db in global memory
 * into the block shared memory and uses it there, to take advantage
 * of faster (but very small) shared memory. Not using this allows
 * larger structures to be used.  Note that, even when this is not
 * defined, the ssetypes and maxscores vectors are kept in shared
 * memory (this does not limit the maximum size of db structures).
 *
 *****************************************************************************/


#undef CUDA5_DEBUG
#undef TESTING

#if defined(CUDA)
//#include <math_functions.h>
#include <helper_cuda.h>
#include <curand_kernel.h>
#endif
#include "saparams.h"
#include "cudaGetDeviceConstantAddresses.h"

#if defined(__DEVICE_EMULATION__) || !defined(CUDA) || defined(CUDA5_DEBUG)
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#endif

#if !defined(CUDA)
#include <math.h>
#include <driver_types.h> /* for struct cudaPitchedPtr, cudaExtent */
#undef __constant__
#define __constant__
#undef __shared__
#define __shared__ static
#undef __global__
#define __global__ 
#undef __device__
#define __device__
#define curandState int
#endif

#define EPS 1.1e-7  /* epsilon for making sure rand is < 1.0 */

#if defined(CUDA) 
#if defined(USE_SHARED_MEMORY)
#define MAXDIM_KERNEL MAXDIM_GPU  /*for shared memory,restrict to MAXDIM_GPU*/
#elif defined(SMALL_MAXDIM)
#define MAXDIM_KERNEL MAXDIM_GPU
#else
#define MAXDIM_KERNEL MAXDIM
#endif
#else
#define MAXDIM_KERNEL MAXDIM      /* otherwise, use largest maxdim */
#endif



/*****************************************************************************
 * 
 * __constant__ memory
 *
 * The query tableau and distance matrix is loaded into constant memory.
 * These are MAXDIM not MAXDIM_KERNEL is constant memory is larger than
 * shared memory so not so restircted (at the moment, on e.g. GTX285,
 * constant memory is 64K but shared memory is only 16K per block).
 *
 *****************************************************************************/

/*tricky: we compile the kernel without CUDA defined for host version
 *annd compile three versions of kernel:
 * (1) with shared memory, maxdim restricted to small so fits in shared
 * (2) without shared memory, large maxdim
 * (3) without shared memory, small maxdim (for Fermi architecture, don't
 *     use shared memory at all, cacheing faster, but still want to run
 *     with small maxdim to not waste so much memory and be consistent
 */

#if !defined(CUDA)
int c_qn_host;    // query structure size
char c_qtab_host[MAXDIM*MAXDIM];  // query tableau
float c_qdmat_host[MAXDIM*MAXDIM];  // query distance matrix
char c_qssetypes_host[MAXDIM]; // main diagonal of c_qn
#endif

/* FIXME because of hack having 3 different versions of everything,
   there is too much __constant__ memory here, but if compiled without
   debug acdtually there is no warning/error and it just gets wrong results
   at runtim for the noshared kernel if __constant__ is used, so 
   changed to __device__ for CUDA 5 just to get it to work.
   It doesn't seem to make any different to speed to use __device__ instead
   of __constant__ though so possibly can just leave it as is.
*/
/* __constant__ */ __device__ int c_qn;    // query structure size
/* __constant__ */ __device__ char c_qtab[MAXDIM*MAXDIM];  // query tableau
/* __constant__ */ __device__ float c_qdmat[MAXDIM*MAXDIM];  // query distance matrix
/* __constant__ */ __device__ char c_qssetypes[MAXDIM]; // main diagonal of c_qn

/* __constant__ */ __device__ int c_qn_noshared;    // query structure size
/* __constant__ */ __device__ char c_qtab_noshared[MAXDIM*MAXDIM];  // query tableau
/* __constant__ */ __device__ float c_qdmat_noshared[MAXDIM*MAXDIM];  // query distance matrix
/* __constant__ */ __device__ char c_qssetypes_noshared[MAXDIM]; // main diagonal of c_qn


/* __constant__ */ __device__ int c_qn_noshared_small;    // query structure size
/* __constant__ */ __device__ char c_qtab_noshared_small[MAXDIM*MAXDIM];  // query tableau
/* __constant__ */ __device__ float c_qdmat_noshared_small[MAXDIM*MAXDIM];  // query distance matrix
/* __constant__ */ __device__ char c_qssetypes_noshared_small[MAXDIM]; // main diagonal of c_qn


#if defined(CUDA)
#if defined(USE_SHARED_MEMORY)
void get_device_constant_addresses
#elif defined(SMALL_MAXDIM)
void get_device_constant_addresses_noshared_small
#else
void get_device_constant_addresses_noshared
#endif
(const_addr_t *const_addr)
{
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qn_addr, c_qn) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qtab_addr, c_qtab) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qdmat_addr, c_qdmat) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qssetypes_addr, c_qssetypes) );
  
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qn_noshared_addr, c_qn_noshared) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qtab_noshared_addr, c_qtab_noshared) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qdmat_noshared_addr, c_qdmat_noshared) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qssetypes_noshared_addr, c_qssetypes_noshared) );

  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qn_noshared_small_addr, c_qn_noshared_small) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qtab_noshared_small_addr, c_qtab_noshared_small) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qdmat_noshared_small_addr, c_qdmat_noshared_small) );
  checkCudaErrors( cudaGetSymbolAddress((void **)&const_addr->c_qssetypes_noshared_small_addr, c_qssetypes_noshared_small) );

}
#endif


#if !defined(CUDA)
/* tricky - we redefined these symbols to the host versions are different */
#define c_qn c_qn_host
#define c_qtab c_qtab_host
#define c_qdmat c_qdmat_host
#define c_qssetypes c_qssetypes_host
#define tscord tscord_cpu
#define deltasd deltasd_cpu
#define thinit thinit_cpu
#define randtypeind randtypeind_cpu
#define icopy icopy_cpu
#define tmscord tmscord_cpu
#define debug_dump_tableau debug_dump_tableu_cpu
#define debug_dump_distmatrix debug_dump_distmatrix_cpu
#else
/* nasty hack -- can't seem to be able to get constant memory in a seaparte
   file, so forced to have 3 different versions of the constant for the
   3 kernel versions. */
#if !defined(USE_SHARED_MEMORY)
#if defined(SMALL_MAXDIM)
#define c_qn c_qn_noshared_small
#define c_qtab c_qtab_noshared_small
#define c_qdmat c_qdmat_noshared_small
#define c_qssetypes c_qssetypes_noshared_small
#define tscord tscord_noshared_small
#define deltasd deltasd_noshared_small
#define thinit thinit_noshared_small
#define randtypeind randtypeind_noshared_small
#define icopy icopy_noshared_small
#define tmscord tmscord_noshared_small
#define debug_dump_tableau debug_dump_tableu_noshared_small
#define debug_dump_distmatrix debug_dump_distmatrix_noshared_small
#else
#define c_qn c_qn_noshared
#define c_qtab c_qtab_noshared
#define c_qdmat c_qdmat_noshared
#define c_qssetypes c_qssetypes_noshared
#define tscord tscord_noshared
#define deltasd deltasd_noshared
#define thinit thinit_noshared
#define randtypeind randtypeind_noshared
#define icopy icopy_noshared
#define tmscord tmscord_noshared
#define debug_dump_tableau debug_dump_tableu_noshared
#define debug_dump_distmatrix debug_dump_distmatrix_noshared
#endif
#endif
#endif

/*****************************************************************************
 * 
 * __device__ functions: callable on GPU only, inlined
 *
 *****************************************************************************/


/* Index into 2d m x n array stored in contiguous memory */
//#define INDEX2D(i,j,m,n) ( ((i)*(n) + (j)) )

/* Get char* to (i,j) element of into 2d array A strored in condiguous
 * memory with pitch (CUDA version of Fortran stride or leading
 * dimension). NB we don't index like INDEX2D but return address as char*
 * which must then be cast to appropriate type, since the address/pitch
 * computations in CUDA are always done in units of bytes, so don't want
 * C address computation using size of actual type */
//#define GET2D(A,i,j,pitch,type) ( *((type *)((char *)(A) + (i)*(pitch) + (j)*sizeof(type) )) )
#define GET2D(A,i,j,pitch,type)  ( ((type *)((A) + ((i)*(pitch))))[(j)] )


/*
 * debug_dump_tableau() - in debug/evmulation build dump tableau
 * 
 * Parameters:
 * tab  (input) 2d char array
 *        tableau one structure. symmetric
 * tab_pitch (input) size_t
 *       pitch of tab
 * n (input) int
 *     dimension of tableau
 *
 * Return value: 
 *   None.
 */
__device__ void debug_dump_tableau(char *tab, size_t tab_pitch, int n)
{
  int i,j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j <= i; j++)
    {
      printf("%02X ", GET2D(tab, i, j, tab_pitch, char));
    }   
    printf( "\n");
  }
}
/*
 * debug_dump_distmatrix() - in debug/evmulation build dump distmatrix
 * 
 * Parameters:
 * dmat  (input) 2d float array
 *        SSE distance matrix for one structure. symmetric
 * dmat_pitch (input) size_t
 *       pitch of dmat1
 * n (input) int
 *     dimenstion of distance matrix
 *
 * Return value: 
 *   None.
 */
__device__ void debug_dump_distmatrix(char *dmat, size_t dmat_pitch, int n)
{
  int i,j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j <= i; j++)
    {
      printf( "%5.3f ", GET2D(dmat, i, j, dmat_pitch, float));
    }   
    printf( "\n");
  }
}



/*
 *
 * tscord - Tableau (discrete) matching score function
 *
 *    Return the tableau matching score between two tableau entries
 *    x and y.
 *    The score is 2 if the tableau entries are equal, 1 if they are
 *    equal in only one position, else -2.
 *
 * Parameters:
 *     x, y - the two two-char tableau codes encoded as 4 bits per char
 *            as per parsetableaux.c
 *
 * Return value:
 *     tableau matching score for x and y
 */
__device__ int tscord(char x, char y)
{
  char xhigh,xlow,yhigh,ylow;
  int score;
  
  xhigh = (x & 0xF0);
  xlow =  (x & 0x0F);
  yhigh = (y & 0xF0);
  ylow =  (y & 0x0F);

  score = ( xhigh == yhigh ? (xlow == ylow ? 2 : 1) :
            (xlow == ylow ? 1 : -2) );
/*
  if (xhigh == yhigh)
  {
    if (xlow == ylow)
      score = 2;
    else
      score = 1;
  }
  else if (xlow == ylow)
    score = 1;
  else
    score = -2;
*/
  return score;
}

/*
 * Compute the score for a given SSE matching between two structures
 * given their tableaux (discrete version), and distnace matrices.
 *
 * The score computed is
 *
 * \sum{i=1,j=1}^{N_A} \sum{j=1,k=1}^{N_B} \zeta(T_{ik},T_{kl}) x_{ik}x{jl}
 *
 * in the QIP formulation where x_{ij} is the binary indicator variable
 * indication SSE i in A matched with SSE j in B. 
 *
 * But actually here we are representing the matching with the ssemap
 * vector so can much more efficiently compute this in only
 * O(N_A^2) with 2 nested loops over the ssemap vector rather than 
 * requring O(N_A^2 N_B^2) with 4 nested loops in the naive implentation
 * of the score computation using indicator variables (required only
 * for using a general purpose QP solver, can do it more efficiently here).
 * 
 * Furthermore, we can actually halve the computation since the tableaux
 * matrices are symmetric by only iterating from k = i .. N_A 
 * inside the outer loop i = 1 .. N_A.
 * 
 * Parameters:
 * tab1  (input) encoded as two 4-bit char code
 *        Tableau for one structure. Symmetric.
 *
 * tab1_pitch (input) size_t
 *         pitch of tab1
 *
 * n1     (input) INTEGER
 *        Dimension of tab1 array.
 *
 * tab2  (input) encoded as two 4-bit char code
 *        Tableau for second structure. Symmetric.
 *
 * tab2_pitch (input) size_t
 *        pitch of tab2
 *
 * n2     (input) INTEGER
 *        Dimension of tab2 matrix.
 *
 * dmat1  (input) 2d float array
 *        SSE distance matrix for one structure. symmetric
 *
 * dmat1_pitch (input) size_t
 *       pitch of dmat1
 *
 * dmat2  (input) 2d float array
 *        SSE distance matrix for second structure. symmetric
 *
 * dmat2_pitch (input) size_t
 *         pitch of dmat2
 *
 * ssemap (input) int vector, dimension(n1)
 *        SSE map vector of dimension n1. Each ssemap(i) is the SSE index
 *        in tab2 that SSE i in tab1 is matched with.
 *
 *
 * Return value:
 *   The tableau matching score for given mapping by ssemap.
 *
 */
__device__ int tmscord(char *tab1, size_t tab1_pitch, int n1,
                       char *tab2, size_t tab2_pitch, int n2,
                       char *dmat1, size_t dmat1_pitch, 
                       char *dmat2, size_t dmat2_pitch,
                       int ssemap[])
{
  int i,j,k,l;
  int score;
  
  score = 0;
  for (i = 0; i < n1; i++)
  {
    for (k = i + 1; k  < n1; k++)
    {
      j = ssemap[i];
      l = ssemap[k];
      /* only add to score when both are mapped to something, and */
      /* diagonal entries are SSE type not angle so don't use them either */
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) || defined(CUDA5_DEBUG)
      assert(j == -1 || i != k && j != l);
#endif
      if (j >= 0 && l >= 0)
      {
        /*
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))        
        fprintf(stderr, "%d %d %d %d %X %X\n", i,j,k,l,
                GET2D(tab1,i,k,tab1_pitch,char), 
                GET2D(tab2,j,l,tab2_pitch,char));
        fprintf(stderr, "%d %d %d %d %f %f (%f)\n",i,j,k,l,
                GET2D(dmat1,i,k,dmat1_pitch,float) ,
                GET2D(dmat2,j,l,dmat2_pitch,float) ,
                fabsf(GET2D(dmat1,i,k,dmat1_pitch,float) - GET2D(dmat2,j,l,dmat2_pitch,float)) );
#endif
        */
        /* don't add score when difference between SSE distances
           exceeds threshold */
        if (fabsf(GET2D(dmat1,i,k,dmat1_pitch,float) - GET2D(dmat2,j,l,dmat2_pitch,float)) <= MXSSED)
        {
          score += tscord(GET2D(tab1,i,k,tab1_pitch,char), GET2D(tab2,j,l,tab2_pitch,char));
        }
      }
    }
  }
  return score;
}


/*
 * deltasd - 
 *
 * Compute the difference in score from due to removing a particular
 * matching of two SSEs and replacing it with a new one.
 * We can do this in O(N_A) time rather than the O(N_A^2) required for
 * computing the score from scratch as in tmscord.
 *
 * 
 * Parameters:
 * tab1  (input) encoded as two 4-bit char code
 *        Tableau for one structure. Symmetric.
 *
 * tab1_pitch (input) size_t
 *         pitch of tab1
 *
 * n1     (input) INTEGER
 *        Dimension of tab1 array.
 *
 * tab2  (input) encoded as two 4-bit char code
 *        Tableau for second structure. Symmetric.
 *
 * tab2_pitch (input) size_t
 *        pitch of tab2
 *
 * n2     (input) INTEGER
 *        Dimension of tab2 matrix.
 *
 * dmat1  (input) 2d float array
 *        SSE distance matrix for one structure. symmetric
 *
 * dmat1_pitch (input) size_t
 *       pitch of dmat1
 *
 * dmat2  (input) 2d float array
 *        SSE distance matrix for second structure. symmetric
 *
 * dmat2_pitch (input) size_t
 *         pitch of dmat2
 *
 * ssemap (input) int vector, dimension(n1)
 *        SSE map vector of dimension n1. Each ssemap(i) is the SSE index
 *        in tab2 that SSE i in tab1 is matched with.
 *
 * sse_i (input) int 
 *        SSE in tab1 that is being replaced with a new matchig
 *
 * old_j (input) int
 *        SSE in tab2 of old matching
 *
 * new_j (input) int
 *        SSE in tab2 of new matching (matched to new_i)
 *
 *
 * Return value:
 *   The difference to add to the current score due to replacing
 *    the old_i <-> old_k matching with the new_i <-> new_k matching.
 *
 */
__device__ int deltasd(char *tab1, size_t tab1_pitch, int n1,
                       char *tab2, size_t tab2_pitch, int n2,
                       char *dmat1, size_t dmat1_pitch,
                       char *dmat2, size_t dmat2_pitch,
                       int ssemap[],
                       int sse_i,
                       int old_j, int new_j)
{
  int k,l;
  int delta = 0;
  float dmat1_i_k;

#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))
  fprintf(stderr,"aaa %d %d %d \n", sse_i, old_j, new_j);
#endif

  for (k = 0; k < n1; k++)
  {
    l = ssemap[k];
    if (l >= 0)
    {
      dmat1_i_k = GET2D(dmat1,sse_i,k,dmat1_pitch,float);
      if (old_j >= 0 && l != old_j && k != sse_i && fabsf(dmat1_i_k - GET2D(dmat2,old_j,l,dmat2_pitch,float)) <= MXSSED)
        delta -= tscord(GET2D(tab1,sse_i,k,tab1_pitch,char), GET2D(tab2,old_j,l,tab2_pitch,char));
      
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))
      fprintf(stderr,"yyy %d %d %d %d %d\n", sse_i, old_j, new_j,k,l);
#endif
      if (new_j >= 0 && l != new_j && k != sse_i && fabsf(dmat1_i_k - GET2D(dmat2,new_j,l,dmat2_pitch,float)) <= MXSSED)
        delta += tscord(GET2D(tab1,sse_i,k,tab1_pitch,char), GET2D(tab2,new_j,l,tab2_pitch,char));
    }
  }
  return delta;
}


 /*
  * Build the initial mapping of the two structurs for heruristic
  * tableaux matching algoriths.
  *
  * we make an initial matching where we just go along
  * the sequence set match of same SSEs e.g. if 1st in query is helix,
  * match that to first helix in db struture, and so on.
  * (Unless LTYPE flag not set, then we don't care about SSE types and 
  * just go along sequence of SSEs).
  * Then compute the score.
  *
  * Parameters:
  *
  *
  * ssetypes1 (input) char vector length n1
  *        vector of SSE types in structure 1
  *
  * n1     (input) INTEGER
  *        Dimension of tab1 matrix
  *
  * ssetypes2 (input) char vector length n2
  *        vector of SSE types in structure 2
  *
  * n2     (input) INTEGER
  *        Dimension of tab2 matrix.
  *
  * lorder (input) LOGICAL
  *        if true, penalize matches between SSEs not maintaining sequence
  *        order between the tableaux i.e. if i < k and j >= l for i,k
  *        indices in tab1 and j,l indices in tab2.
  *
  *        
  * ssemap (output) INTEGER vector, dimension (n1)
  *        solution SSE map vector of dimension n1. 
  *        Each ssemap(i) is the SSE index
  *        in tab2 that SSE i in tab1 is matched with.
  *
  * revmap (output) INTEGER vector, dimension(n2)
  *     reverse ssemap: revmap(j) for j index in tab2 is the index i
  *     in tab1 that matches that sse i.e. if ssemap(i) = j then
  *     revmap(j) = ssemap(i) and vice versa, for quick lookup of what
  *     is matched so we can easily check that one-to-one mapping maintained
  *
  * state (input/output) State for RNG
  *
  * Return value:
  *        on exit, status of the computation
  *        =  0 : successful exit
  *        =  1 : cannot setup intial ssemap with both lorder and ltype
  */
__device__ int thinit(char ssetypes1[], int n1,
                      char ssetypes2[], int n2,
                      int lorder,
                      int ssemap[], int revmap[],
                      curandState *state)
{
  int i,j;
  int info = 0;
  float randnum;

#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))
  for(int k = 0; k < n1; k++)
    fprintf(stderr,"%02X ", ssetypes1[k]);
  fprintf(stderr, "\n");
  for(int k = 0; k < n2; k++)
    fprintf(stderr,"%02X ", ssetypes2[k]);
  fprintf(stderr, "\n");
#endif

  /* initialize ssemap to all -1 meaning no match for each sse */
  for (i = 0; i < n1; i++)
    ssemap[i] = -1;
  for (j = 0; j < n2; j++)
    revmap[j] = -1;

  /* initial SSE map set by matching along sequence, only matching SSEs
   * of the saem type if LFTYPE flag is set.
   */
  
  j = 0;
  for (i = 0; i < n1; i++)
  {
#if defined(CUDA)
    randnum = curand_uniform(state);
#else
    randnum = drand48();
#endif
#if defined(CUDA) && defined( CUDA5_DEBUG )
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
    printf("tid = %d thinit randnum = %f\n", tid ,randnum);
#endif
    if (randnum < INIT_MATCHPROB)
    {
      while (j < n2 && ssetypes1[i] != ssetypes2[j])
        j++;
      if (j >= n2) 
      {
        /* not all SSEs in tab1 are mapped, but that's OK */
        info = 0;
        return info;
      }
      else
      {
        ssemap[i] = j;
        revmap[j] = i;
        j++;
      }
    }
  }
  return info;
}

/*
 * find the index of first SSE of same type in tableaux that is not
 * already mapped or -1 if not found
 *
 * Parameters:
 *
 * ssetypesvec (input) char vector length n
 *        vector of SSE types in structure 
 *
 * n     (input) INTEGER
 *        Dimension of tableaux, legnth of ssetypesvec
 *
 * startind (input) INTEGER
 *        SSE index to start at in tab
 *
 * ssetype (input) CHARACTER*2
 *        SSE type as two charcter string 'xa' etc.
 *
 * smap    (input) INTEGER vector, dimension(n1)
 *         each smap(i) is index in other tableau it is already mapped
 *         to, or 0 for not mapped.
 *
 * endind (input) INTEGER
 *         last SSE index to consider in tab 
 *
 * state (input/output) State for RNG
 */
__device__ int randtypeind(char ssetypesvec[], int n,
                           int startind, char ssetype, int smap[], int endind,
                           curandState *state)
{
  int i,indi,rti;
  int indlist[MAXDIM];
  float randnum;
  unsigned int randidx;

  i = startind;
  indi = 0;
  rti = -1;
  for (i = startind; i < endind; i++)
  {

#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) || defined(CUDA5_DEBUG)
  assert(i >= 0);
  assert(i < n);
#endif

    if (ssetypesvec[i] == ssetype && smap[i] < 0)
      indlist[indi++] = i;
  }

  if (indi == 1)
    rti = indlist[0];
  else if (indi > 1) 
  {
#if defined(CUDA)
    randnum = curand_uniform(state);
#else
    randnum = drand48();
#endif
    randidx = (unsigned int)((randnum - EPS) * indi);
    rti = indlist[randidx];
  }
  return rti;
}


/*
 * integer vector copy y <- x
 *
 */
__device__ void icopy(int n, int x[], int y[])
{
  int i;
  for (i = 0; i < n; i++)
    y[i] = x[i];
}

/*****************************************************************************
 * 
 * __global__ functions: GPU kernels, callable from host
 *
 *****************************************************************************/

/*
 * CUDA GPU kernel for tableau matching using simulated annealing.
 *
 * We make an initial matching where we just go along
 * the sequence set match of same SSEs e.g. if 1st in query is helix,
 * match that to first helix in db struture, and so on.
 * (Unless LTYPE flag not set, then we don't care about SSE types and 
 * just go along sequence of SSEs).
 * Then compute the score.
 *
 * Then we use simulated annealing to improve the score. At each 
 * iteration a random SSE is chosen to be remapped to a random
 * other SSE (obeying constraints that are set) or
 * mapped to no SSE in the other structure.
 *
 *
 * "Embarrasingly parallel" version: just do all the loops in here,
 * each thread does a different database structure.
 * The query tableau and distance matrix is placed in constant memory
 * for faster access (constant memory is cached but very limited size:
 * we certainly can't put the whole db of structures there for instance).
 *
 * Parameters:
 *
 * dbsize (input) INTEGER
 *        number of strucures in database
 *
 * lorder (input) LOGICAL
 *        if true, penalize matches between SSEs not maintaining sequence
 *        order between the tableaux i.e. if i < k and j >= l for i,k
 *        indices in tab1 and j,l indices in tab2.
 *
 * lsoln  (input) LOGICAL
 *         if true, return the SSE mapping for the best solution found.
 *
 * maxtart (input) INTEGER
 *         number of restarts (iteratinos of cooling schedule).
 *         Should be a multiple of blocksize.
 *
 * d_qdmat  (input) float array, dimension (n1,n1)
 *        SSE distance matrix for query structure. symmetric
 *
 * d_qdmat_pitch (input) size_t
 *         pitch of d_qdmat
 *
 * d_tableaux (input) pointer to char arrays, CUDA Pitched Pointer
 *        Pointer to database of tableaux
 *
 * tableaux_extent (input) cudaExtent
 *        Extent structure for d_tableaux
 *
 * d_ordrers ( input) pointer to ints
 *        Pointer to database of orders (order of each db tableau)
 *
 * d_distmatcies (input) pitched pointer to float arrays
 *        Pointer to database of distance matrices
 *
 * distmatrices_extent (input) cudaExtent
 *       Extent structure for d_distmatrices
 *
 * outscore  (output) INTEGER vector, dimension (dbsize)
 *        scores of matching query with each db structure
 *
 * outssemap (output) INTEGER array, dimension (dbszie, n1)
 *        solution SSE map vector of dimension n1 for each db structure
 *        Each ssemap(d,i) is the SSE index
 *        in dbentry d that SSE i in query is matched with.
 *
 * state - (in/out) state for RNG
 */
#if defined(CUDA)
#if defined(USE_SHARED_MEMORY)
__global__ void sa_tabsearch_gpu
#elif defined(SMALL_MAXDIM)
__global__ void sa_tabsearch_gpu_noshared_small
#else
__global__ void sa_tabsearch_gpu_noshared
#endif
#else
void sa_tabsearch_host
#endif
                                (int dbsize,
                                 int lorder,
                                 int lsoln,
                                 int maxstart,
                                 cudaPitchedPtr d_tableaux,
                                 cudaExtent tableaux_extent,
                                 int *d_orders,
                                 cudaPitchedPtr d_distmatrices,
                                 cudaExtent distmatrices_extent,
                                 int *outscore,
                                 int *outssemap,
                                 curandState *state)
{
 
  /*
   * 
   * __shared__ memory
   *
   * Each block of threads copies one database tableau and distance matrix
   * from the global memory into shared memory. Each thread in the block
   * runs the simulated annealing schedule (with different RNG) on the quey
   * and this shared tableu+distmatrix, so the 'restarts' are pallelized
   * within the block. 
   *
   * Note the shared memory is very restriced in size (16K) so we can
   * only fit limited size structures in it.
   *
   */
#if defined(USE_SHARED_MEMORY)
  __shared__ char s_tab[MAXDIM_KERNEL*MAXDIM_KERNEL];
  __shared__ float s_dmat[MAXDIM_KERNEL*MAXDIM_KERNEL];
#endif
  __shared__ char s_ssetypes[MAXDIM_KERNEL]; // TODO maybe shouldn't use this in shared
  __shared__ int s_maxscores[128]; // FIXME should be max threads in block
  __shared__ int s_maxscore_threadid;

  /*
   * automatic (register and local) memory
   */    

//  const int THREAD_N = blockDim.x * gridDim.x;  // total number of threads
#if defined(CUDA)
  const int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
  const int blockid = blockIdx.x;            // block id
  const int gridDimx = gridDim.x;            // number of blocks in grid
  const int blockDimx  = blockDim.x;         // number of threads in block
  const int threadIdxx = threadIdx.x;        // thread id in the block
  curandState localState = state[tid];/* cache state in fast local memory */
#else
  const int tid = 0;
  const int blockid = 0;
  const int gridDimx = 1;
  const int blockDimx = 1;
  const int threadIdxx = 0;
  curandState localState = 0;/*unused*/
#endif

#ifdef CUDA5_DEBUG
  if(tid==0)
    printf("c_qn = %d\nc_qn_noshared = %d\nc_qn_noshared_small = %d\n",
           c_qn, c_qn_noshared, c_qn_noshared_small);
#endif

  int revmap[MAXDIM_KERNEL];  /* reverse ssemap: revmap(j) for j index in
                         tab2 is the index i in tab1 that matches that
                         sse i.e. if ssemap(i) = j then revmap(j) =
                         ssemap(i) and vice versa, for quick lookup of
                         what is matched so we can easily check that
                         one-to-one mapping maintained revmap has
                         dimension (n2) */
    
  int bestmap[MAXDIM]; /* best ssemap feound. this has dimenion (n1) */
  int ssemap[MAXDIM];
  int maxscore,score,newscore;
  int iter;
  float temp;
  float randnum;
  int startj,endj,k,oldj,newj;
  int ssei;
  char *tab1;
  char *tab2;
  char *dmat1; // we use char* not float* to do pitched pointer arithmetic
  char *dmat2;
  int n1,n2;
  int restart;
  int dbi;
  //int iState;
  size_t tab1_pitch, tab2_pitch,dmat1_pitch,dmat2_pitch;
  int i,j;
  int blockmaxscore;
  int delta;


#if defined(__DEVICE_EMULATION__)
  fprintf(stderr, "running in device emulation mode\n");
  fprintf(stderr, "sizeof(int) == %d\n", sizeof(int));
#if defined(USE_SHARED_MEMORY)
  fprintf(stderr, "using shared memory\n");
#else
  fprintf(stderr, "NOT using shared memory\n");
#endif
  fprintf(stderr, "MAXDIM_KERNEL = %d\n", MAXDIM_KERNEL);
#endif

#if !defined(CUDA)
  fprintf(stderr, "running on host\n");
#endif
  

  n1 = c_qn;
  tab1 = c_qtab;
  tab1_pitch = MAXDIM; /* NB MAXDIM not MAXDIM_KERNEL, see comments on c_qtab */
  dmat1 = (char*)c_qdmat;
  dmat1_pitch = MAXDIM * sizeof(float);


  // each of the gridDim.x blocks does as many as needed to do whole database
  for (dbi = blockid; dbi < dbsize; dbi += gridDimx) 
  {

    n2 = d_orders[dbi];

    // get the tableau aray for db entry index dbi using pitched pointer
    char *d_tableauxPtr = (char *)d_tableaux.ptr;
    size_t tableauxPitch = d_tableaux.pitch;
    size_t tableauxSlicePitch = tableauxPitch * tableaux_extent.height;
    char *tableauxSlice = d_tableauxPtr + dbi * tableauxSlicePitch;
    tab2 = tableauxSlice;
    tab2_pitch = tableauxPitch;

    // and similarly for distmatrices (2d float arrays)
    char *d_distmatricesPtr = (char *)d_distmatrices.ptr;
    size_t distmatricesPitch = d_distmatrices.pitch;
    size_t distmatricesSlicePitch = distmatricesPitch * distmatrices_extent.height;
    char *distmatricesSlice = d_distmatricesPtr + dbi * distmatricesSlicePitch;
    dmat2 = distmatricesSlice;
    dmat2_pitch = distmatricesPitch;

    
    // set the s_ssetypes vector as main diagonal of this db instance tableau
    // in parallel (each thread in block does one element)
    for (j = threadIdxx; j < n2; j += blockDimx)
      s_ssetypes[j] = GET2D(tab2,j,j,tab2_pitch,char); // use global not shared so no sync required

#if defined(USE_SHARED_MEMORY)
    //
    // parallel copy (each thread in block does as many elements as needed)
    // of the db entry for this block into the shared memory for the block
    // we'll have each thread do one row of the copy (may leave threads idle
    // since likely to have more threads in block than rows in tableau).
    //
    for (i = threadIdxx; i < n2; i += blockDimx)
      for (j = 0; j < n2; j++)
      {
        *(s_tab + i*MAXDIM_KERNEL + j) = GET2D(tab2,i,j,tab2_pitch,char);
        *(s_dmat + i*MAXDIM_KERNEL + j) = GET2D(dmat2,i,j,dmat2_pitch,float);
      }
    tab2_pitch = MAXDIM_KERNEL;  /* pitch is now leading dimension in shared */
    dmat2_pitch = MAXDIM_KERNEL*sizeof(float);
#else
    /* not using shared memory, just point the s_* variables to the 
       global memory */
    char *s_tab = tab2;
    char *s_dmat = dmat2;
#endif

#if defined(CUDA)
    // sync point so all threads have loaded into shared memory
    __syncthreads();
#endif

#if defined(CUDA5_DEBUG)
    if (tid == 0) {
      printf("maxstart = %d lsoln = %d\n", maxstart, lsoln);
      printf( "tab1 (n = %d):\n", n1);
    debug_dump_tableau(tab1, tab1_pitch, n1);
      printf( "dmat1:\n");
    debug_dump_distmatrix(dmat1, dmat1_pitch, n1);
    printf("c_q_ssetypes:");
    for (int l = 0; l < n1; l++)
      printf("%02X ", c_qssetypes[l]);
    printf("\n");
    printf( "tab2 (n = %d):\n", n2);
    debug_dump_tableau(s_tab, tab2_pitch, n2);
      printf( "s_dmat:\n");
    debug_dump_distmatrix((char *)s_dmat, dmat2_pitch, n2);
    printf("s_ssetypes:");
    for (int l = 0; l < n2; l++)
      printf("%02X ", s_ssetypes[l]);
    printf("\n");

    }
#endif

    maxscore = -99999;


    // each of the blockDim.x threads in the block does as many iterations
    // as need to get to maxstart restarts
    for (restart = 0; restart < maxstart; restart += blockDimx)
    {
      /* setup initial mapping */
      thinit(c_qssetypes, n1, s_ssetypes, n2, lorder, ssemap, revmap,
             &localState);

      score =  tmscord(tab1, tab1_pitch, n1, s_tab, tab2_pitch, n2, 
                       dmat1, dmat1_pitch, 
                       (char *)s_dmat, dmat2_pitch,
                       ssemap);
      if (score > maxscore)
      {
        maxscore = score;
        icopy(n1, ssemap, bestmap);
      }

      temp = TEMP0;

      for (iter = 0; iter < MAXITER; iter++)
      {
        /* generate neighbour state by picking random SSE in tab1 and
           moving its mapping to a radnom SSE in tab2, maintaining 
           constraints */
#if defined(CUDA)
        randnum = curand_uniform(&localState);
#else
        randnum = drand48();
#endif
        ssei = ((randnum - EPS) * n1);

#if defined (DEBUG) && !defined(CUDA)
//        fprintf(stderr, "xxx %f %d\n", randnum, ssei);
#endif
#if defined(CUDA) && defined(CUDA5_DEBUG)
        if (tid==0)
//          printf( "xxx %f %d\n", randnum, ssei);
#endif


        if (lorder)
        {
          startj = ssemap[ssei];
          k = ssei;
          while (startj < 0 && k >= 0)
          {
            startj = ssemap[k];
            k--;
          }
          if (startj < 0)
            startj = n2;
          if (ssei == n1-1)
            endj = n2;
          else if (ssemap[ssei+1] < 0)
          {
            endj = -1;
            k = 1;
            while (endj == -1 && ssei + k < n1)
            {
              endj = ssemap[ssei + k];
              k++;
            }
          }
          else
            endj = ssemap[ssei+1];
        }
        else
        {
          startj = 0;
          endj = n2;
        }
        newj = randtypeind(s_ssetypes, n2, startj, 
                           c_qssetypes[ssei],
                           revmap, endj, &localState);
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) || defined(CUDA5_DEBUG)
#ifdef CUDA
        if(tid==0)
#endif
        {
          printf("%d %d %d %d %d %d %d\n", tid , restart, iter, ssei, startj, endj, newj);
          printf( "%d ssemap: ", tid);
          for (int q = 0; q < n1; q++)
            printf( "%d ", ssemap[q]);
          printf( "\n");
        }
#endif
        oldj = ssemap[ssei];
        delta = deltasd(tab1, tab1_pitch, n1, s_tab, tab2_pitch, n2,
                        dmat1, dmat1_pitch,
                        (char*)s_dmat, dmat2_pitch,
                        ssemap, ssei, oldj, newj);
//#undef TESTING
#ifdef TESTING 
#if defined(__DEVICE_EMULATION__) || !defined(CUDA) || defined(CUDA5_DEBUG)
        int revnewmap[MAXDIM_KERNEL],ssenewmap[MAXDIM];
        icopy(n1, ssemap, ssenewmap);
        icopy(n2, revmap, revnewmap);
          if (newj > -1)
          {
            ssenewmap[ssei] = newj;
            if (oldj > -1)
              revnewmap[oldj] = -1;
            revnewmap[newj] = ssei;
          }
          else
          {
            /* the SSE was removed from the matching */
            if (oldj > -1)
            {
              revnewmap[ssenewmap[ssei]] = -1;
              revnewmap[oldj] = -1;
            }
            ssenewmap[ssei] = -1;
          }
        int fullscore = tmscord(tab1, tab1_pitch, n1, s_tab, tab2_pitch, n2,
                                dmat1, dmat1_pitch,
                                (char*)s_dmat, dmat2_pitch,
                                ssenewmap);
//        fprintf(stderr, "zzz %d %d %d\n", delta, score+delta,fullscore);
        assert(score + delta == fullscore);
#endif
#endif

        newscore = score + delta;
        if (newscore > maxscore)
        {
          maxscore = newscore;
          if (lsoln)
          {
            icopy(n1, ssemap, bestmap);
            if (newj > -1)
              bestmap[ssei] = newj;
            else
              bestmap[ssei] = -1;
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA)) 
            fprintf(stderr, "NNN %d %d %d\n", ssei, oldj, newj);
            fprintf(stderr, "%d bestmap: ", tid);
            for (int q = 0; q < n1; q++)
              fprintf(stderr, "%d ", bestmap[q]);
            fprintf(stderr, "\n");
#endif
          }
        }

#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))
        fprintf(stderr, "%d %d %d %f %d %d %f\n",tid, restart, iter, temp, score, newscore,
                expf((float)delta / temp));
#endif
#if defined(CUDA)
        randnum = curand_uniform(&localState);
#else
        randnum = drand48();
#endif
        if (expf((float)delta / temp) > randnum)
        {
          /* accept the move, update ssemap and revmap accordingly */
          score = newscore;
          if (newj > -1)
          {
            ssemap[ssei] = newj;
            if (oldj > -1)
              revmap[oldj] = -1;
            revmap[newj] = ssei;
          }
          else
          {
            /* the SSE was removed from the matching */
            if (oldj > -1)
            {
              revmap[ssemap[ssei]] = -1;
              revmap[oldj] = -1;
            }
            ssemap[ssei] = -1;
          }
        }

        temp *= ALPHA;

      }
    }
    
    s_maxscores[threadIdxx] = maxscore;

#if defined(CUDA)
    // synchronization point: now we need to find max score over each thread in
    // block for that block's db structure.
    __syncthreads();
#endif

    // reduction (MAX) over threads in block to get max score for
    // TODO make this a proper reduction operation instead of a 
    // loop in a single thread
    if (threadIdxx == 0)
    {
      s_maxscore_threadid = 0;
      blockmaxscore = s_maxscores[0];
      for (i = 1; i < blockDimx; i++)
      {
        if (s_maxscores[i] > blockmaxscore)
        {
          blockmaxscore = s_maxscores[i];
          s_maxscore_threadid = i;
        }
      }
#if defined(__DEVICE_EMULATION__) || (defined(DEBUG) && !defined(CUDA))
      fprintf(stderr, "%d says maxscore is %d for %d\n", tid, blockmaxscore,dbi);
#endif
      outscore[dbi] = blockmaxscore;
    }

    if (lsoln)
    {
#if defined(CUDA)
    // synchronization point: need to wait for threadid 0 to have found max
      __syncthreads();
#endif
      // Now that we have the best score and the thread that found it,
      // THAT thread only will put its bestssemap as the output SSE map
      if (threadIdxx == s_maxscore_threadid)
        icopy(n1, bestmap, outssemap + dbi * MAXDIM);
    }
  }
  state[tid] = localState; /* copy back new state from local cache */
}


