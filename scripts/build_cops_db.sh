#!/bin/sh
#
# File:    build_cops_db.sh
# Author:  Alex Stivala
# Created: May 2010
#
#
# build_cops_db.sh - build tableaux database for COPS benchmark data set
#
# Usage: build_cops_db.sh querydir dbfile
#
#         querydir is directory to put query input tableaux into
#         dbfile   is basename of tableaux database to create, will create
#                  dbfile.tableaux.pickle, dbfile.distmatrix.pickle and
#                  dbfile.tableauxdistmatrixdb.ascii
#
# Builds tableaux for queries and database for the COPS benchmark data set
# (Frank et al. 1999 "COPS Benchmark: interactive analysis of database 
# search methods" Bioinformatics 26(4):574-575) available from
# http://benchmark.services.came.sbg.ac.at/
#
# Requires the buildtableauxdb.py and pytableaucreate.py and convdb2.py
#  scripts in PATH.
#
# WARNING: dbfile and files in querydir are overwritten if they exist.
#
# $Id: build_cops_db.sh 3632 2010-05-12 02:07:26Z alexs $

COPS_ROOT=${HOME}/cops-benchmark-2009-6-full
COPS_PDB_QUERIES=${COPS_ROOT}/queries/pdb
COPS_PDB_DB=${COPS_ROOT}/database/pdb

if [ $# -ne 2 ]; then
    echo "Usage: $0 querydir dbfile" >&2
    exit 1
fi

QUERYDIR=$1
DBFILE=$2

OPTIONS="-p none -35 -t dssp"

if [ ! -d $QUERYDIR ]; then
    mkdir $QUERYDIR
fi

tableaux_pickle=${DBFILE}.tableaux.pickle
distmatrix_pickle=${DBFILE}.distmatrix.pickle
tableauxdb=${DBFILE}.tableauxdb.ascii

for query in ${COPS_PDB_QUERIES}/*.pdb
do
  qid=`basename $query .pdb`
  qfile=${QUERYDIR}/${qid}.input
  echo $tableauxdb > $qfile
  echo "T T F" >> $qfile # options: type, order, output
  pytableaucreate.py -f -b $OPTIONS $query >> $qfile
done


buildtableauxdb.py $OPTIONS $COPS_PDB_DB $tableaux_pickle
buildtableauxdb.py -d $OPTIONS $COPS_PDB_DB $distmatrix_pickle

convdb2.py $tableaux_pickle $distmatrix_pickle > $tableauxdb

