#!/bin/sh
#
# File:    build_query_input.sh
# Author:  Alex Stivala
# Created: November 2008
#
# build_query_input.sh - build tsrchd input files for a list of query ids
#
# Usage: build_query_input.sh outdir <query.list
#
#   outdir is name of diretory which is created, and each tableau 
#   in ASCII format for use with tsrchd_sparse etc. is
#   created as a separate file in that directory, in format for input
#   for use with qptabmatch_allpairs.py for example
#
#   To stdout is written the ASCII format db of all the tableaux+dist matrices
#   (just all the ones written to outdir concatenated together with
#   blank line between each). 
#
#

# location of ASTRAL SCOP pdbstyle files
#ASTRAL_DIR=/usr/local/ASTRAL/pdbstyle-sel-gs-bib-95-1.75
ASTRAL_DIR=/usr/local/ASTRAL/pdbstyle-1.75

# tableaux+distmatrix db file
#TABLEAUX_DB=/home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75/omegadistmatrixdb.ascii
#TABLEAUX_DB=/home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75/tableauxdistmatrixdb.sorted.ascii
TABLEAUX_DB=/home/alexs/tableauxdb/minlen4-ASTRAL-sel-gs-bib-95-1.75/tableauxdistmatrixdb.ascii

if [ $# -ne 1 ]; then
    echo "Usage: $0  outdir < querylist" 2>&1
    exit 1
fi
outdir=$1

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

# pytableaucreate.py options
#tabopts="-35 -f -t dssp -p none"
tabopts="-m4 -35 -f -t dssp -p none"
#tabopts="-n -35 -f -t dssp -p none"

while read scopsid
do
    div=`echo $scopsid | cut -c3-4`
    pdbfile=${ASTRAL_DIR}/${div}/${scopsid}.ent
    echo "${TABLEAUX_DB}" > ${outdir}/${scopsid}.input
    echo "T T F" >> ${outdir}/${scopsid}.input  # options: type,order,output
    pytableaucreate.py ${tabopts} ${pdbfile} >> ${outdir}/${scopsid}.input
    # append distance matrix, removing identifier on first line
    pytableaucreate.py -d ${tabopts} ${pdbfile} | awk 'NR > 1' >> ${outdir}/${scopsid}.input 
done

