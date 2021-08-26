#!/bin/sh
#
# File:    build_nh3d_db.sh
# Author:  Alex Stivala
# Created: September 2008
#
# build_nh3d_db.sh - build tableaux database for Nh3D data set
#
# Usage: build_nh3d_db.sh outdir 
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
# builds database of tableaux, using pytableaycreate.py,
# for the Nh3D data set (Thiruv et al 2005 BMC Struct. Biol. 5:12)
#

# location of Nh3D data set, PDB format files
NH3D_PDB_DIR=/local/charikar/Nh3D/v3.0


if [ $# -ne 1 ]; then
    echo "Usage: $0  outdir" 2>&1
    exit 1
fi
outdir=$1

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

# pytableaucreate.py options
tabopts="-35 -f -t dssp -p none"

first=1
for pdbfile in ${NH3D_PDB_DIR}/*.pdb
do
    cathid=`basename ${pdbfile} .pdb`
    # dodgy: remove periods so that CATH id fits in 8 chars... 
    # hopefully will get no duplicates...
    cathid=`echo ${cathid} | tr -d .`
    if [ $first -eq 0 ]; then
        echo
    else
        first=0
    fi
    pytableaucreate.py ${tabopts} -i ${cathid} ${pdbfile} | tee ${outdir}/${cathid}.tableaudistmatrix
    # append distance matrix, removing identifier on first line
    pytableaucreate.py -d ${tabopts} -i ${cathid} ${pdbfile} | awk 'NR > 1' |tee -a ${outdir}/${cathid}.tableaudistmatrix 
done

