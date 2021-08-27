#!/bin/bash
#SBATCH --job-name="host_cudaSaTabsearch"
#SBATCH --time=0-00:10:00
#SBATCH --partition=slim
#SBATCH --output=host_cudaSaTabsearch-%j.out
#SBATCH --error=host_cudaSaTabsearch-%j.err

echo -n "started at: " >&2; date >&2
uname -a >&2

module load cudatoolkit

/usr/bin/time ./cudaSaTabsearch -c -r4096  < d2phlb1.input

echo -n "ended at: " >&2; date >&2
