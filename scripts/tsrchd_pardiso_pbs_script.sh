#!/bin/bash
#
# File:    tsrchd_pardiso_pbs_script.sh
# Author:  Alex Stivala
# Created: July 2009
#
# PBS script for submitting QP tableau search jobs on tango.vpac.org
# requires PATH and PYTHONPATH already set up in environment 
#
# $Id: tsrchd_pardiso_pbs_script.sh 2941 2009-11-15 03:27:36Z astivala $
 	
#PBS -N d1ae6h1_tsrchd_pardiso
#PBS -l walltime=23:0:0
#PBS -l nodes=1
#PBS -v MKL_NUM_THREADS=1



QUERY=d1ae6h1
INPUT_DIR=${HOME}/phd/qptabsearch/data
OUTDIR=.
TSRCHD=${HOME}/phd/qptabsearch/src/tsrchd_pardiso
TIME=/usr/bin/time

cd $PBS_O_WORKDIR
set CONV_RSH = ssh


$TIME $TSRCHD -t < ${INPUT_DIR}/${QUERY}.input > ${OUTDIR}/${QUERY}.out 2> ${OUTDIR}/${QUERY}.err


