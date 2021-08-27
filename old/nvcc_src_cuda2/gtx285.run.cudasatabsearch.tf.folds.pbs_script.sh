#!/bin/bash
#
# run.tlocsd.folds.pbs_script.sh
#
# Usage: qsub gtx285.run.cudasatabsearch.tf.folds.pbs_script.sh
#
# run cudaSaTabsearch
# on the manually defind fold queries with the ORDER constraint disabled
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: gtx285.run.cudasatabsearch.tf.folds.pbs_script.sh 3342 2010-02-15 04:11:04Z alexs $
#

# run on GTX285 (nv3)
#PBS -l nodes=nv3
#PBS -q gpu
#PBS -N folds_gtx285_gpucudaSaTabsearch
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
    queryfile=${INPUT_DIR}/${query}.tf.input
    $TIME ${TLOCSD} ${tlocsdopts}  < ${queryfile} > ${query}.tf.out 2> ${query}.tf.err
done

