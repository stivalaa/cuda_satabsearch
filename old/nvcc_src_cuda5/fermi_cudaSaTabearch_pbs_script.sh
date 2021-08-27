#!/bin/bash

# pbs launching script:
 	
# run on Fermi architecture (Tesla M2070) 

#PBS -q gpu
#PBS -N gpucudaSaTabsearch_fermi
#PBS -l walltime=0:0:20

uname -a >&2

module load cuda/5.0.35
cd $PBS_O_WORKDIR
set CONV_RSH = ssh


#/usr/bin/time ./cudaSaTabsearch -r4096  < d1ubia_.input

#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ae6h1.input
#/usr/bin/time ./cudaSaTabsearch -r128  < d1ubia_.input
#/usr/bin/time ./cudaSaTabsearch -r4096 < 2qp2-1.input
/usr/bin/time ./cudaSaTabsearch -r4096 < d2phlb1.input
#/usr/bin/time ./cudaSaTabsearch < d1ubia_.input

