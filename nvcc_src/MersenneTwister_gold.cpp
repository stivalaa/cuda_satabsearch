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

/* modified by ADS: added RandomHost()
   $Id: MersenneTwister_gold.cpp 3310 2010-02-05 08:22:12Z alexs $
*/

#ifdef _WIN32
#   pragma warning( disable : 4996 ) // disable deprecated warning 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MersenneTwister.h"
#include "dci.h"




static mt_struct MT[MT_RNG_COUNT];
static uint32_t state[MT_NN];



extern "C" void initMTRef(const char *fname){
    
    FILE *fd = fopen(fname, "rb");
    if(!fd){
        printf("initMTRef(): failed to open %s\n", fname);
        printf("TEST FAILED\n");
        exit(0);
    }

    for (int i = 0; i < MT_RNG_COUNT; i++){
        //Inline structure size for compatibility,
        //since pointer types are 8-byte on 64-bit systems (unused *state variable)
        if( !fread(MT + i, 16 /* sizeof(mt_struct) */ * sizeof(int), 1, fd) ){
            printf("initMTRef(): failed to load %s\n", fname);
            printf("TEST FAILED\n");
            exit(0);
        }
    }

    fclose(fd);

    const int iRng = 0;
    unsigned int seed = 777;
    MT[iRng].state = state;
    sgenrand_mt(seed, &MT[iRng]);
}


extern "C" float RandomHost()
{
  const int iRng = 0;



  return ((float)genrand_mt(&MT[iRng]) + 1.0f) / 4294967296.0f;
}

extern "C" void RandomRef(
    float *h_Random,
    int NPerRng,
    unsigned int seed
){
    int iRng, iOut;

    for(iRng = 0; iRng < MT_RNG_COUNT; iRng++){
        MT[iRng].state = state;
        sgenrand_mt(seed, &MT[iRng]);

        for(iOut = 0; iOut < NPerRng; iOut++)
           h_Random[iRng * NPerRng + iOut] = ((float)genrand_mt(&MT[iRng]) + 1.0f) / 4294967296.0f;
    }
}


#define PI 3.14159265358979f
void BoxMuller(float& u1, float& u2){
    float   r = sqrtf(-2.0f * logf(u1));
    float phi = 2 * PI * u2;
    u1 = r * cosf(phi);
    u2 = r * sinf(phi);
}

extern "C" void BoxMullerRef(float *h_Random, int NPerRng){
    int i;

    for(i = 0; i < MT_RNG_COUNT * NPerRng; i += 2)
        BoxMuller(h_Random[i + 0], h_Random[i + 1]);
}
