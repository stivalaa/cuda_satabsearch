#!/bin/sh
###############################################################################
#
# qptabmatch_.nh3d.sh - run the QP tableau matching on Nh3D data set
#
# File:    qptabmatch_nh3d.sh
# Author:  Alex Stivala
# Created: September 2008
#
# Run QP tableau matching on 
# the Nh3D data set (Thiruv et al 2005 BMC Struct. Biol. 5:12)
# with the 73 queries defined in Pelta et al 2008 BMC Bioinformatics 9:161
#
# Usage:
#     qptabmatch_nh3d.sh indir dbfile outdir
#
#     indir is directory containing the tableaux+dismatrices built with
#     build_nh3d_db.sh
#
#     dbfile is the nh3d tableau+distmatrices database file built with
#     build_nh3d_db.sh
#
#     outdir is diretory to place corresponding output from tsrchd_sparse
#     created if it does not exist
#     WARNNG: .out files in outdir overwritten if they exist
#
# Environment variables:
#
#   PATH must contain the location of tsrchd_sparse
#
# $Id: qptabmatch_n3hd.sh 1887 2008-09-13 00:52:50Z astivala $
# 
###############################################################################


# List of query CATH identifiers, from the Additional File 1 spreadsheet
# for Pelta et al 2008
QUERY_LIST="1.10.1040 1.10.1320 1.10.533 1.10.645 1.20.1280 1.20.210 1.20.5 1.20.840 2.10.25 2.10.260 2.10.270 2.10.90 2.170.16 2.170.230 2.170.290 2.170.40 2.30.110 2.30.18 2.30.230 2.30.29 2.30.40 2.40.155 2.40.160 2.40.180 2.40.340 2.40.50 2.60.130 2.60.260 2.60.420 2.60.90 2.70.100 2.70.180 2.70.220 2.70.98 3.10.105 3.10.170 3.10.270 3.10.330 3.10.400 3.20.120 3.20.140 3.20.19 3.20.70 3.20.90 3.30.1530 3.30.1690 3.30.240 3.30.559 3.30.560 3.30.60 3.30.990 3.40.1210 3.40.1380 3.40.225 3.40.720 3.60.100 3.60.120 3.60.20 3.60.40 3.60.90 3.90.1280 3.90.1300 3.90.1350 3.90.1580 3.90.510 3.90.850 4.10.1080 4.10.1090 4.10.220 4.10.260 4.10.480 4.10.540 4.10.790"

if [ $# -ne 3 ]; then
    echo "Usage: $0 indir dbfile outdir" 2>&1
    exit 1
fi

indir=$1
dbfile=$2
outdir=$3

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

for qid in ${QUERY_LIST}
do
    tmpfile=/tmp/nh3d.${qid}.$$
    echo ${dbfile} >$tmpfile
    echo "T T F" >> $tmpfile
    tfile=${indir}/`echo $qid | tr -d .`.tableaudistmatrix
    cat ${tfile} >> $tmpfile
    tsrchd_sparse < $tmpfile > ${outdir}/${qid}.out
    rm $tmpfile
done

