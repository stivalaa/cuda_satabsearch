#!/bin/bash
###############################################################################
#
# genpermutedqueries.sh - generate permuted folding pattern queries
#
# File:    genpermutedqueries.sh
# Author:  Alex Stivala
# Created: March 2009
#
#
# Generate 5 random permutations of the tableaux of each of the 8 
# structures we use for testing. The idea is then to run tsrchd_sparse
# with the ordering constraint disabled, and check that the actual
# structures are still matched (even though the ordering of SSEs in
# the query structure has been disturbed randomly).
#
#
# Usage:
#     genpermutedqueries.sh outdir
#
#     outdir is the directory to write the query files to. It is created
#     if it does not exist. WARNING: files overwritten if they do exist.
#     The files are named d1ubia_.1.input d1ubia_.1.permutation etc.
#     where the .permutation file is the output of pytableaucreate, which
#     contains the permutation = 8,5,1,7,6,3,2,4 (for example) line
#     showing the permutation used (needed to reconstruct alignments later
#     if desired with ssepermutationremap.py) and .input is the input
#     to tsrchd_sparse containging the database name and options (and
#     not the permutation).
#
# Uses the Python scripts pytableaucreate.p to create tableaux for input
# to the FORTRAN tsrchd_sparse program
#
# Environment variables:
#
#   PATH must contain the location of the Python scripts, ie where this
#   script itself is and the ptgraph/ directory with pytableaucreate.py etc.,
#   and the location of tsrchd_sparse.
#   The dssp program must also be in the PATH.
#
#   PYTHONPATH must contain the directory containing the ptsecstruct.py 
#   and other Python modules used by the Python scripts.
#
# $Id: genpermutedqueries.sh 2120 2009-03-21 01:04:39Z astivala $
# 
###############################################################################


# Root of ASTRAL divided PDB style hierarchy
ASTRAL_ROOT=/local/charikar/ASTRAL/pdbstyle-1.73

# Tableau+distmatrix database file
TABLEAUXDB=/local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii

# number of permutations of each structure to make
NUM_PERMUTATIONS=5

# list of the structures we use as queries
STRUCTS="d1ubia_ d1tttb1 d1ae6h1 d1bhne_ d1h6rb_ d2phlb1 d1tima_ d1f6dc_"


if [ $# -ne 1 ]; then
    echo "Usage: $0 outdir" >&2
    exit 1
fi

outdir=$1
if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

for struct in $STRUCTS
do
    pnum=1
    while [ $pnum -le $NUM_PERMUTATIONS ]
    do
        div=`echo $struct | cut -c3-4`
        pdbfile=${ASTRAL_ROOT}/${div}/${struct}.ent
        pfile=${outdir}/${struct}.${pnum}.permutation
        qfile=${outdir}/${struct}.${pnum}.input
        pytableaucreate.py -u -b -35 -f -t dssp -p none $pdbfile > $pfile
        echo $TABLEAUXDB > $qfile
        echo "T F F" >> $qfile     # options: type,order,output
        awk 'NR > 1' <  $pfile >> $qfile
        pnum=`expr $pnum + 1`
    done
done

