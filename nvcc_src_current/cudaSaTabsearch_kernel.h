#ifndef CUDASATABSEARCH_KERNEL_H
#define CUDASATABSEARCH_KERNEL_H
/*****************************************************************************
 * 
 * File:    cudaSaTabsearch_kernel.h
 * Author:  Alex Stivala
 * Created: January 2010
 *
 * $Id: cudaSaTabsearch_kernel.h 4751 2013-11-20 01:43:29Z astivala $
 *
 *****************************************************************************/

#include "saparams.h"

__global__ void sa_tabsearch_gpu(int dbsize,
                                 int lorder,
                                 int lsoln,
                                 int maxstart,
                                 cudaPitchedPtr d_tableaux,
                                 cudaExtent tableaux_extent,
                                 int *d_orders,
                                 cudaPitchedPtr d_distmatrices,
                                 cudaExtent distmatrices_extent,
                                 int *outscore,
                                 int ssemap[],
                                 curandState *state);

__global__ void sa_tabsearch_gpu_noshared(int dbsize,
                                          int lorder,
                                          int lsoln,
                                          int maxstart,
                                          cudaPitchedPtr d_tableaux,
                                          cudaExtent tableaux_extent,
                                          int *d_orders,
                                          cudaPitchedPtr d_distmatrices,
                                          cudaExtent distmatrices_extent,
                                          int *outscore,
                                          int ssemap[],
                                          curandState *state);

__global__ void sa_tabsearch_gpu_noshared_small(int dbsize,
                                          int lorder,
                                          int lsoln,
                                          int maxstart,
                                          cudaPitchedPtr d_tableaux,
                                          cudaExtent tableaux_extent,
                                          int *d_orders,
                                          cudaPitchedPtr d_distmatrices,
                                          cudaExtent distmatrices_extent,
                                          int *outscore,
                                                int ssemap[],
                                                curandState *state);

void sa_tabsearch_host          (int dbsize,
                                 int lorder,
                                 int lsoln,
                                 int maxstart,
                                 cudaPitchedPtr d_tableaux,
                                 cudaExtent tableaux_extent,
                                 int *d_orders,
                                 cudaPitchedPtr d_distmatrices,
                                 cudaExtent distmatrices_extent,
                                 int *outscore,
                                 int ssemap[],
                                 int *state);


extern int MAXSTART; /* number of restarts */

#endif /* CUDASATABSEARCH_KERNEL_H */
