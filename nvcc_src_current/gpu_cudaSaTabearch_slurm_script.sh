#!/bin/bash
#SBATCH --job-name="gpu_cudaSaTabsearch"
#SBATCH --time=0-00:00:20
#SBATCH --partition=gpu
#SBATCH --output=gpu_cudaSaTabsearch-%j.out
#SBATCH --error=gpu_cudaSaTabsearch-%j.err

echo -n "started at: " >&2; date >&2

uname -a >&2

module load cudatoolkit

#/usr/bin/time ./cudaSaTabsearch -r4096  < d1ubia_.input

#/usr/bin/time ./cudaSaTabsearch < $HOME/phd/qptabsearch/data/d1ae6h1.input
#/usr/bin/time ./cudaSaTabsearch -r128  < d1ubia_.input
#/usr/bin/time ./cudaSaTabsearch -r4096 < 2qp2-1.input
/usr/bin/time ./cudaSaTabsearch -r4096 < d2phlb1.input
#/usr/bin/time ./cudaSaTabsearch < d1ubia_.input

echo -n "ended at: " >&2; date >&2
