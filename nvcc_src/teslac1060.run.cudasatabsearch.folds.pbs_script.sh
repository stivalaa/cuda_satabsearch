#!/bin/bash
#
# run.tlocsd.folds.pbs_script.sh
#
# Usage: qsub run.tlocsd.folds.pbs_script.sh
#
# run tlocsd
# on the manually defind fold queries
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: teslac1060.run.cudasatabsearch.folds.pbs_script.sh 3342 2010-02-15 04:11:04Z alexs $
#

# run on Tesla C1060 (nv4)
#PBS -l nodes=nv4
#PBS -q gpu
#PBS -N folds_teslac1060_gpucudaSaTabsearch
#PBS -l walltime=1:0:0

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

TIME=/usr/bin/time
TLOCSD=./cudaSaTabsearch
INPUT_DIR=../../data

QUERYID_LIST="d1ubia_ betagrasp d1ae6h1 d1tima_ d1bhne_ d1h6rb_ d1tttb1 d2phlb1 d1f6dc_"

tlocsdopts=""


echo "# Run as: " $0 $* 
echo "# at: " `date`
echo "# on: " `uname -a`

for query in ${QUERYID_LIST}
do
    queryfile=${INPUT_DIR}/${query}.input
    $TIME ${TLOCSD} ${tlocsdopts}  < ${queryfile} > ${query}.out 2> ${query}.err
done

