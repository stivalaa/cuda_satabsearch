#!/bin/sh
#
# File:    normalize_fischer_multiquery.sh
# Author:  Alex Stivala
# Created: October 2008
#
# normalize_fischer_multiquery.sh - build all normalized scores files for Fischer db
#
# Usage: normalize_fischer_multiquery.sh 
#
# Run from the fischer/ subdirectory: processes files under there
# makeing norm1/ norm2/ norm3/ subdirectories,
# using also tableaux and contact map input files under data/
# via script normalize_tabmatch.py 
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_fischer_multiquery.sh 3909 2010-07-11 04:45:25Z alexs $

TABMATCH_FISCHER_DIR=.

TABLEAUX_FISCHER_DIR=${HOME}/phd/qptabsearch/data/fischer_db


NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing tableau search Fischer results with norm ${norm}..."
    outdir=${TABMATCH_FISCHER_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    infile=fischer.out
    outfile=${outdir}/fischer.out
    normalize_tabmatch.py ${norm} -m ${TABLEAUX_FISCHER_DIR} < ${infile} > ${outfile}
done

