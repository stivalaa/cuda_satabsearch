#ifndef SAPARAMS_H
#define SAPARAMS_H
/*****************************************************************************
 * 
 * File:    saparams.h
 * Author:  Alex Stivala
 * Created: January 2010
 *
 * Constants for CUDA simulated annealing tableau search
 *
 *
 *****************************************************************************/

/* maximum size of tableaux / distance matrices that can be read */
#define MAXDIM 111

/* maximum size of tableaux / distance matrices handled by GPU (shared mem) */
#define MAXDIM_GPU 96  /* TODO conditional compilation for this */
/* Fermi (compute capability 2.0) has 48 KB shared but only 16 KB for 
   earlier architectures. For larger shared memory (Fermi arch.) need
   -arch sm_20 compile option and cuda/3.0 not cuda/2.3 */
/*#define MAXDIM_GPU 96*/ /* 48 KB shared memory on Fermi*/

/* max length of structure labels (d1ubia_ etc.) */
#define LABELSIZE 8

/* SSE distance difference threshold */
#define MXSSED (float)4.0

/* number of iterations in cooling schedule */
#define MAXITER 100

/* Initial temperature */
#define TEMP0 (float)10.0

/* factor to multiply temperature by at each step */
#define ALPHA (float)0.95

/* default number of restarts if not specified */
#define DEFAULT_MAXSTART 128

/* probability that we should try to find an initial matching for each SSE */
#define INIT_MATCHPROB 0.5

/* max number of host threads (one for each GPU, one for host) we can have */
#define MAX_THREADS 8

#endif /* SAPARAMS_H */

