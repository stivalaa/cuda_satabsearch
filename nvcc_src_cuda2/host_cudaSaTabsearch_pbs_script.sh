#!/bin/bash

# pbs launching script:
 	
#PBS -l nodes=1
#PBS -N cpu_cudaSaTabsearch
#PBS -l walltime=10:0:0

uname -a >&2

module load cuda/2.3
cd $PBS_O_WORKDIR
set CONV_RSH = ssh


/usr/bin/time ./cudaSaTabsearch -c -r4096  < $HOME/phd/qptabsearch/data/d1ubia_.input

