#!/bin/sh
#
# File:    build_nh3d_cm.sh
# Author:  Alex Stivala
# Created: September 2008
#
# build_nh3d_cm.sh - build contact maps for Nh3D data set
#
# Usage: build_nh3d_cm.sh outdir 
#
#   outdir is name of diretory which is created, and each contact map
#   in ASCII format for use with MSVNS4MaxCMO (Pelta et al 2008) 
#   or other program using this format of contact matrix is
#   created as a separate file in that directory, in format for input
#   for use with msvns4maxcmo_allall.py for example
#
# builds contact maps, using pconpy.py,
# for the Nh3D data set (Thiruv et al 2005 BMC Struct Biol 5:12)
#
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

# pconpy.py options
threshold=7.0
pconpyopts="--cmaplist --threshold=${threshold} --seq_separation=2"

for pdbfile in ${NH3D_PDB_DIR}/*.pdb
do
    cathid=`basename ${pdbfile} .pdb`
    pconpy.py ${pconpyopts} --pdb=${pdbfile} --output=${outdir}/${cathid}.cm_a${threshold}
done

