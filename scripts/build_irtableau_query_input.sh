#!/bin/sh
#
# File:    build_irtableau_query_input.sh
# Author:  Alex Stivala
# Created: March 2010
#
# build_irtableau_query_input.sh - 
#                          build irtableau input files for a list of query ids
#
# Usage: build_irtableau_query_input.sh outdir <query.list
#
#   outdir is name of diretory which is created, and each tableau vector
#   in ASCII format for use with irtableau is
#   created as a separate file in that directory
#

# location of tableaux vector db file for all ASTRAL SCOP pdbstyle files
ALL_VECTOR_DB=/home/alexs/tableauxdb/ASTRAL-1.75/tableaux.vectordb.ascii

# tableaux vector db file specified in input to search
TABLEAUX_DB=/home/alexs/tableauxdb/jASTRAL-sel-gs-bib-95-1.75/tableaux.vectordb.ascii

if [ $# -ne 1 ]; then
    echo "Usage: $0  outdir < querylist" 2>&1
    exit 1
fi
outdir=$1

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

while read scopsid
do
    echo "${TABLEAUX_DB}" > ${outdir}/${scopsid}.input
    fgrep ${scopsid} $ALL_VECTOR_DB >> ${outdir}/${scopsid}.input
done

