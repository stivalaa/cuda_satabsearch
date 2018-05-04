
/*****************************************************************************
 *
 * This file contains code adapted from CUDA SDK 2.3 
 *
 * The Mersenne Twister RNG kernel from the CUDA SDK, modified so that
 * instead of standalone generating numbers and storing them,
 * it is called to get next number like the standard C library rand() etc,.
 * from the device only (not a __global__ kernel any more).
 *
 * $Id: MersenneTwister_kernel.cu 3350 2010-02-18 00:32:08Z alexs $
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

#include <stdio.h>
#include <cutil_inline.h>
#include "MersenneTwister.h"



__device__ static mt_struct_stripped ds_MT[MT_RNG_COUNT];
/* static mt_struct_stripped h_MT[MT_RNG_COUNT]; */


///////////////////////////////////////////////////////////////////////////////
//
// Initialize RNG state vecotr
//
///////////////////////////////////////////////////////////////////////////////
__device__ void InitRandomGPU(int *iState0, unsigned int mt[MT_NN])
{
  const int      tid = blockDim.x * blockIdx.x + threadIdx.x;
  int iState;

  int iRng = tid;

  //Load bit-vector Mersenne Twister parameters
  mt_struct_stripped config = ds_MT[iRng];

  //Initialize current state
  mt[0] = config.seed;
  for(iState = 1; iState < MT_NN; iState++)
    mt[iState] = (1812433253U * (mt[iState - 1] ^ (mt[iState - 1] >> 30)) + iState) & MT_WMASK;
  *iState0 = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Get next random number in sequence, for this thread.
//
///////////////////////////////////////////////////////////////////////////////
__device__ float RandomGPU(
  int *iState,
  unsigned int mt[MT_NN]
){
    const int      tid = blockDim.x * blockIdx.x + threadIdx.x;
//    const int THREAD_N = blockDim.x * gridDim.x;

    int iState1, iStateM;
    unsigned int mti, mti1, mtiM, x;

    int iRng = tid;
    //Load bit-vector Mersenne Twister parameters
    mt_struct_stripped config = ds_MT[iRng];
    mti1 = mt[0];

    //iState1 = (iState +     1) % MT_NN
    //iStateM = (iState + MT_MM) % MT_NN
    iState1 = *iState + 1;
    iStateM = *iState + MT_MM;
    if(iState1 >= MT_NN) iState1 -= MT_NN;
    if(iStateM >= MT_NN) iStateM -= MT_NN;
    mti  = mti1;
    mti1 = mt[iState1];
    mtiM = mt[iStateM];
    
    x    = (mti & MT_UMASK) | (mti1 & MT_LMASK);
    x    =  mtiM ^ (x >> 1) ^ ((x & 1) ? config.matrix_a : 0);
    mt[*iState] = x;
    *iState = iState1;
    
    //Tempering transformation
    x ^= (x >> MT_SHIFT0);
    x ^= (x << MT_SHIFTB) & config.mask_b;
    x ^= (x << MT_SHIFTC) & config.mask_c;
    x ^= (x >> MT_SHIFT1);
    
    //Convert to (0, 1] float and return
    float newrand = ((float)x + 1.0f) / 4294967296.0f;
#ifdef __DEVICE_EMULATION__
    fprintf(stderr, "%d RandomGPU %f\n", tid, newrand);
#endif
    return newrand;
}


#ifdef UNUSED

////////////////////////////////////////////////////////////////////////////////
// Transform each of MT_RNG_COUNT lanes of NPerRng uniformly distributed 
// random samples, produced by RandomGPU(), to normally distributed lanes
// using Cartesian form of Box-Muller transformation.
// NPerRng must be even.
////////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265358979f
__device__ void BoxMuller(float& u1, float& u2){
    float   r = sqrtf(-2.0f * logf(u1));
    float phi = 2 * PI * u2;
    u1 = r * __cosf(phi);
    u2 = r * __sinf(phi);
}

__global__ void BoxMullerGPU(float *d_Random, int NPerRng){
    const int      tid = blockDim.x * blockIdx.x + threadIdx.x;
    const int THREAD_N = blockDim.x * gridDim.x;

    for(int iRng = tid; iRng < MT_RNG_COUNT; iRng += THREAD_N)
        for(int iOut = 0; iOut < NPerRng; iOut += 2)
            BoxMuller(
                d_Random[iRng + (iOut + 0) * MT_RNG_COUNT],
                d_Random[iRng + (iOut + 1) * MT_RNG_COUNT]
            );
}
#endif
