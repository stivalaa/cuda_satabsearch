#!/bin/bash
#
# File:    tsrchn_pardiso_pbs_script.sh
# Author:  Alex Stivala
# Created: July 2009
#
# PBS script for submitting QP tableau search (numeric) jobs on tango.vpac.org
# requires PATH and PYTHONPATH already set up in environment 
#
# $Id: tsrchn_pardiso_pbs_script.sh 2906 2009-11-06 00:05:34Z astivala $
 	
#PBS -N tsrchn_pardiso
#PBS -l walltime=23:0:0
#PBS -l nodes=1
#PBS -v MKL_NUM_THREADS=1



QUERY=d1f6dc_
INPUT_DIR=${HOME}/phd/qptabsearch/data
OUTDIR=.
TSRCHN=${HOME}/phd/qptabsearch/src/tsrchn_pardiso
TIME=/usr/bin/time

cd $PBS_O_WORKDIR
set CONV_RSH = ssh


$TIME $TSRCHN -t < ${INPUT_DIR}/${QUERY}.omega.input > ${OUTDIR}/${QUERY}.tsrchn.out 2> ${OUTDIR}/${QUERY}.tsrchn.err


