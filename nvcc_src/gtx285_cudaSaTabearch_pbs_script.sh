#!/bin/bash

# pbs launching script:
 	
# run on GTX 285 (nv3)
#PBS -l nodes=nv3
#PBS -q gpu
#PBS -N gpucudaSaTabsearch_gtx285
#PBS -l walltime=1:0:0

uname -a >&2

module load cuda/2.3
cd $PBS_O_WORKDIR
set CONV_RSH = ssh

#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ae6h1.input
#/usr/bin/time ./cudaSaTabsearch < 2qp2-1.input
#/usr/bin/time ./cudaSaTabsearch < d2phlb1.input
#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ubia_.input
/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d2phlb1.input
#/usr/bin/time ./cudaSaTabsearch -r4096  < $HOME/phd/qptabsearch/data/d1ubia_.input

