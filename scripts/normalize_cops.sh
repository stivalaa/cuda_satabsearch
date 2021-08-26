#!/bin/sh
#
# File:    normalize_cops_multiquery.sh
# Author:  Alex Stivala
# Created: October 2008
#
# normalize_cops.sh - build all normalized scores files for COPS db
#
# Usage: normalize_cops.sh 
#
# Run from the cops/ subdirectory: processes files under there
# makeing norm1/ norm2/ norm3/ subdirectories,
# using also tableaux files for COPS benchmark data set
# via script normalize_tabmatch.py 
#
# The scripts/ directory (containing this scrpipt and the abovementioned ones)
# must be in the PATH
#
# $Id: normalize_cops.sh 3635 2010-05-12 06:48:14Z alexs $

TABMATCH_COPS_DIR=.

TABLEAUXDB=/home/alexs/tableauxdb/COPS/COPS.tableaux_db_and_queries.ascii


NORMTYPES="1 2 3"

for norm in ${NORMTYPES}
do
    echo "normalizing tableau search COPS results with norm ${norm}..."
    outdir=${TABMATCH_COPS_DIR}/norm${norm}
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir}
    fi
    for infile in ${TABMATCH_COPS_DIR}/*.out
    do
        qid=`basename ${infile} .out`
        outfile=${outdir}/${qid}.out
        normalize_tabmatch.py ${norm} ${qid} ${TABLEAUXDB} < ${infile} > ${outfile}
    done
done

