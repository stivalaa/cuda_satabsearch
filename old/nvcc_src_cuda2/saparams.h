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
 * $Id: saparams.h 3545 2010-03-31 00:40:50Z alexs $
 *
 *****************************************************************************/

/* maximum size of tableaux / distance matrices that can be read */
#define MAXDIM 111

/* maximum size of tableaux / distance matrices handled by GPU (shared mem) */
#define MAXDIM_GPU 32

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

