#!/bin/sh
#
# File:    normalize_astral40multiquery.sh
# Author:  Alex Stivala
# Created: February 2010
#
# normalize_astral40multiquery.sh - build  normalized scores for ASTRAL40
#
# Usage: normalize_astral40multiquery.sh 
#
# Run from the ASTRAL 40% all-all results directory: processes files under there
# making norm1/ norm2/ norm3/ subdirectories,
# using also tableaux and contact map input files 
# via scripts normalize_tabmatch.py 
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_astral40multiquery.sh 3928 2010-07-15 01:22:09Z alexs $

TABMATCH_ASTRAL40_DIR=.
TABLEAUXDB=/home/alexs/tableauxdb/ASTRAL-sel-gs-bib-40-1.75/tableauxdistmatrixdb.sorted.ascii

NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing  ASTRAL40 results with norm ${norm}..."
    outdir=${TABMATCH_ASTRAL40_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_ASTRAL40_DIR}/*.out
    do
        outfile=${outdir}/`basename $infile .out`.out
        normalize_tabmatch.py ${norm} -m ${TABLEAUXDB} < ${infile} > ${outfile}
    done
done

