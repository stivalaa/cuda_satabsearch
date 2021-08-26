#!/bin/sh
#
# File:    normalize_multiquery200.sh
# Author:  Alex Stivala
# Created: February 2010
#
# normalize_multiquery200.sh - build  normalized scores for query200 set
#
# Usage: normalize_multiquery200.sh 
#
# Run from the query200 results directory: processes files under there
# makeing norm1/ norm2/ norm3/ subdirectories,
# using also tableaux and contact map input files under data/
# via scripts normalize_tabmatch.py 
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_multiquery200.sh 3345 2010-02-16 03:42:05Z alexs $

TABMATCH_QUERY200_DIR=.
TABLEAUXDB=/home/alexs/tableauxdb/ASTRAL-1.75/tableauxdistmatrixdb.ascii

NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing  query200 results with norm ${norm}..."
    outdir=${TABMATCH_QUERY200_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_QUERY200_DIR}/*.out
    do
        outfile=${outdir}/`basename $infile .out`.out
        normalize_tabmatch.py ${norm} -m ${TABLEAUXDB} < ${infile} > ${outfile}
    done
done

