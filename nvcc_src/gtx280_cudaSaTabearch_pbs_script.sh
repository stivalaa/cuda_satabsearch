#!/bin/bash

# pbs launching script:
 	
# run on GTX 280 (nv1)
#PBS -l nodes=nv1
#PBS -q gpu
#PBS -N gpucudaSaTabsearch_gtx280
#PBS -l walltime=1:0:0

uname -a >&2

module load cuda
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ae6h1.input
#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ubia_.input
#/usr/bin/time ./cudaSaTabsearch < d2phlb1.input
#/usr/bin/time ./cudaSaTabsearch < d1ubia_.input
/usr/bin/time ./cudaSaTabsearch -r4096  < $HOME/phd/qptabsearch/data/d1ubia_.input

