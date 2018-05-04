
/*****************************************************************************
 *
 * This file contains code adapted from CUDA SDK 2.3 
 *
 * The Mersenne Twister RNG kernel from the CUDA SDK, modified so that
 * instead of standalone generating numbers and storing them,
 * it is called to get next number like the standard C library rand() etc,.
 * from the device only (not a __global__ kernel any more).
 *
 * $Id: MersenneTwister_host.cu 3350 2010-02-18 00:32:08Z alexs $
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


/* __device__ static mt_struct_stripped ds_MT[MT_RNG_COUNT]; */
static mt_struct_stripped h_MT[MT_RNG_COUNT];



//Load twister configurations
void loadMTGPU(const char *fname){
    FILE *fd = fopen(fname, "rb");
    if(!fd){
        printf("initMTGPU(): failed to open %s\n", fname);
        printf("TEST FAILED\n");
        exit(0);
    }
    if( !fread(h_MT, sizeof(h_MT), 1, fd) ){
        printf("initMTGPU(): failed to load %s\n", fname);
        printf("TEST FAILED\n");
        exit(0);
    }
    fclose(fd);
}

//Initialize/seed twister for current GPU context
void seedMTGPU(unsigned int seed){
    int i;
    //Need to be thread-safe
    mt_struct_stripped *MT = (mt_struct_stripped *)malloc(MT_RNG_COUNT * sizeof(mt_struct_stripped));

    for(i = 0; i < MT_RNG_COUNT; i++){
        MT[i]      = h_MT[i];
        MT[i].seed = seed;
    }
    CUDA_SAFE_CALL( cudaMemcpyToSymbol("ds_MT", MT, sizeof(h_MT)) );

    free(MT);
}
