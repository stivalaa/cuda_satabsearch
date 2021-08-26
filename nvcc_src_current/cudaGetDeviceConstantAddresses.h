#ifndef CUDAGETDEVICEONSTANTADDRESSES_H
#define CUDAGETDEVICEONSTANTADDRESSES_H
/*****************************************************************************
 * 
 * File:    cudaDeviceGetConstantAdresses.h
 * Author:  Alex Stivala
 * Created: November 2013
 *
 * Definition of struct and prototype for function to get addresses
 * of all device constant variables since in CUDA 5.0 can no longer
 * using string symbol names in cudaMemcpyToSymbol().
 *
 * $Id: cudaGetDeviceConstantAddresses.h 4751 2013-11-20 01:43:29Z astivala $
 *
 *****************************************************************************/


/* structure to contain addresses of device constants */
typedef struct const_addr_s {
  int    *c_qn_addr;
  char   *c_qtab_addr;
  float  *c_qdmat_addr;
  char   *c_qssetypes_addr;

  int    *c_qn_noshared_addr;
  char   *c_qtab_noshared_addr;
  float  *c_qdmat_noshared_addr;
  char   *c_qssetypes_noshared_addr;
  
  int    *c_qn_noshared_small_addr;
  char   *c_qtab_noshared_small_addr;
  float  *c_qdmat_noshared_small_addr;
  char   *c_qssetypes_noshared_small_addr;
} const_addr_t;

/* return struct with addresses of device constants */
void get_device_constant_addresses(const_addr_t *const_addr);
void get_device_constant_addresses_noshared(const_addr_t *const_addr);
void get_device_constant_addresses_noshared_small(const_addr_t *const_addr);



#endif /* CUDAGETDEVICEONSTANTADDRESSES_H */
