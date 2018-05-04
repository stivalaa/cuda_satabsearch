#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=1
#PBS -N cpu_cudaSaTabsearch
#PBS -l walltime=0:10:0

uname -a >&2

module load cuda/5.0.35
cd $PBS_O_WORKDIR
set CONV_RSH = ssh



/usr/bin/time ./cudaSaTabsearch -c -r4096  < d2phlb1.input

