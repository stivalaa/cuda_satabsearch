#!/bin/bash
#
# File:    tsrchd_umfpack_pbs_script.sh
# Author:  Alex Stivala
# Created: July 2009
#
# PBS script for submitting QP tableau search jobs on tango.vpac.org
# requires PATH and PYTHONPATH already set up in environment 
#
# Run from qptabsearch/ directory, i.e. :
#
#       qsub scripts/tsrchd_umfpack_pbs_script.sh
#
# $Id: tsrchd_umfpack_pbs_script.sh 2726 2009-08-05 04:56:11Z astivala $
 	
#PBS -N tsrchd_umfpack
#PBS -l walltime=16:0:0
#PBS -l nodes=1
#PBS -v MKL_NUM_THREADS=1


QUERY=d1f6dc_
OUTDIR=tango_results_umfpack

TIME=/usr/bin/time
cd $PBS_O_WORKDIR
set CONV_RSH = ssh


$TIME src/tsrchd_sparse < data/${QUERY}.input > ${OUTDIR}/${QUERY}.out 2> ${OUTDIR}/${QUERY}.err


