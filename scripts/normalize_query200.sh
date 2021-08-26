#!/bin/sh
#
# File:    normalize_query200.sh
# Author:  Alex Stivala
# Created: October 2008
#
# normalize_query200.sh - build all normalized scores files for query200 set
#
# Usage: normalize_query200.sh 
#
# Run from the query200 results directory: processes files under there
# makeing norm1/ norm2/ norm3/ subdirectories,
# using also tableaux and contact map input files under data/
# via scripts normalize_tabmatch.py 
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_query200.sh 2092 2009-03-09 23:19:14Z astivala $

TABMATCH_QUERY200_DIR=.
TABLEAUXDB=/home/alexs/tableauxdb/ASTRAL-1.75/tableauxdistmatrixdb.ascii

NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing QP tableau search query200 results with norm ${norm}..."
    outdir=${TABMATCH_QUERY200_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_QUERY200_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
#        echo $qid
        normalize_tabmatch.py ${norm} ${qid} ${TABLEAUXDB} < ${infile} > ${outfile}
    done
done

