###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala
# Created   : January 2010
#
# $Id: Makefile 284 2009-02-23 04:06:18Z alexs $
#
# Makefile for building the the CUDA simulated annealing tableau matching code
#
#
# The default target "all" will compile and link the program.
#
# The "clean" target will remove object and executable files.
#
# set MODE=DEVICEEMU to build with device emulation and debugging on. 
# default is to build with optimizations on and no debug or profile, for GPU
# execution.
# set MODE=GPUDEBUG to build with debugging on GPU with cuda-gdb enabled.
# Otherwise,
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile, for GPU
# execution.
#
###############################################################################


# Requires the CUDA SDK 
CUDA_INCLUDE = /usr/local/cuda/2.3/cuda/include
CUDA_SDK_ROOT = /usr/local/cuda/2.3/cuda/sdk
CUDA_SDK_INCLUDE = $(CUDA_SDK_ROOT)/C/common/inc
CUDA_SDK_UTIL_LIB = $(HOME)/cuda/sdk/common/lib  # libcutil must be built by user


CUDA_KERNELS  = cudaSaTabsearch_kernel.cu
CUDA_HOSTCODE = cudaSaTabsearch.cu MersenneTwister_host.cu
CSRCS         = parsetableaux.c genmtrand.c 
CPPSRCS       = MersenneTwister_gold.cpp
HDRS          = MersenneTwister.h parsetableaux.h cudaSaTabsearch_kernel.h MersenneTwister_kernel.cu
CUDA_OBJS     = cudaSaTabsearch_kernel.o cudaSaTabsearch_kernel_noshared.o \
                cudaSaTabsearch_kernel_noshared_small.o
CUDA_HOST_OBJS= $(CUDA_HOSTCODE:.cu=.o)
CPP_OBJS      = $(CPPSRCS:.cpp=.o)
COBJS         = $(CSRCS:.c=.o)
OBJS          = $(CUDA_OBJS) $(COBJS) $(CUDA_HOST_OBJS) $(CPP_OBJS) cudaSaTabsearch_host.o
SRCS          = $(CSRCS) $(CPPSRCS) $(CUDA_KERNELS) $(CUDA_HOSTCODE)


CC           = g++ 
NVCC         = nvcc

CUDA_CPPFLAGS   = -I$(CUDA_SDK_INCLUDE)
HOST_CPPFLAGS   = 
CDEBUG     = -g -DDEBUG
OPTFLAGS   = -O3 
NVCCOPTFLAGS = --ptxas-options=-v --use_fast_math
CFLAGS     = $(OPTFLAGS) -Wall
#              the following warnings are not implied by -Wall
CFLAGS     += -Wextra -Wfloat-equal  \
              -Wundef -Wshadow \
              -Wpointer-arith -Wcast-qual -Wcast-align\
              -Wwrite-strings \
              -Wmissing-declarations -Wunreachable-code
ifeq ($(MODE), DEVICEEMU)
  NVCFLAGS   = $(CDEBUG) -deviceemu
else
  ifeq ($(MODE),GPUDEBUG)
    CFLAGS += -g -G
    NVCFLAGS = -g -G
  else
    NVCFLAGS   = $(OPTFLAGS) $(NVCCOPTFLAGS)
  endif
  ifeq ($(MODE),DEBUG)
    CFLAGS = $(CDEBUG)
    NVCFLAGS = $(CDEBUG) -G
  endif
endif


# NVIDIA CUDA compute capability - we do not require double precision
GPU_ARCH = sm_10

LD         = nvcc
LDFLAGS    = -arch $(GPU_ARCH)
ifeq ($(MODE),DEVICEEMU)
  LDFLAGS += -deviceemu
else
  ifeq ($(MODE),GPUDEBUG)
    LDFLAGS += -g -G
  endif
  ifeq ($(MODE),DEBUG)
    LDFLAGS += -g -G
   endif
endif

LDLIBPATH  = -L$(CUDA_SDK_UTIL_LIB)
LDLIBS     = -lcutil -lm

MAKEDEPEND = gcc -MM  $(CPPFLAGS)
DEPENDFILE = .depend

# Program to build TAGS file for EMACS
MAKETAGS   = etags

all: cudaSaTabsearch

cudaSaTabsearch: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)


# tricky: we compile the kernel without CUDA defined for host version
cudaSaTabsearch_host.o: cudaSaTabsearch_kernel.cu
	$(CC) -x c++ -c $(CFLAGS) $(HOST_CPPFLAGS) -I$(CUDA_INCLUDE) -o $@  $<

cudaSaTabsearch_constant_host.o: cudaSaTabsearch_constant.cu
	$(CC) -x c++ -c $(CFLAGS) $(HOST_CPPFLAGS) -I$(CUDA_INCLUDE) -o $@  $<

# annd compile three versions of kernel:
#  (1) with shared memory, maxdim restricted to small so fits in shared
#  (2) without shared memory, large maxdim
#  (3) without shared memory, small maxdim (for Fermi architecture, don't
#      use shared memory at all, cacheing faster, but still want to run
#      with small maxdim to not waste so much memory and be consistent
#      with others).
cudaSaTabsearch_kernel.o : cudaSaTabsearch_kernel.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -DCUDA -DUSE_SHARED_MEMORY -arch $(GPU_ARCH) -o $@ $<

cudaSaTabsearch_kernel_noshared.o : cudaSaTabsearch_kernel.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -DCUDA -arch $(GPU_ARCH) -o $@ $<

cudaSaTabsearch_kernel_noshared_small.o : cudaSaTabsearch_kernel.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -DCUDA -DSMALL_MAXDIM -arch $(GPU_ARCH) -o $@ $<


depend: $(SRCS) $(HDRS)
	$(MAKEDEPEND) $(SRCS) $(HDRS) > $(DEPENDFILE)

# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS) cudaSaTabsearch


%.o : %.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -DCUDA -arch $(GPU_ARCH) -o $@ $<

%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CUDA_CPPFLAGS) -o $@ $<

include $(DEPENDFILE)


