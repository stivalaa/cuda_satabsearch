#!/bin/sh
#
# File:    build_tops_files_pdb.sh
# Author:  Alex Stivala
# Created: April 2010
#
# build_tops_files_pdb.sh - build TOPS files from hierarchy of PDB files
#
# Usage: build_tops_files_pdb.sh pdb_root outdir 
#
#   outdir is name of diretory which is created, and a .tops file
#   created as a separate file in that directory, for each .ent file
#   in the PDB hierarchy.
#
#   pdb_root is the root of the PDB hierarchy.
#
#
# $Id: build_tops_files_pdb.sh 3575 2010-04-22 00:35:26Z alexs $

# location of TOPS directory, contains tops.def etc.
# Note all the .dssp and .tops files are temporarily created here,
# (tops.def has these specifications)
TOPS_ROOT=$HOME/Tops

if [ $# -ne 2 ]; then
    echo "Usage: $0 pdb_root outdir" 2>&1
    exit 1
fi

pdb_root=$1
outdir=$2

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

cd $TOPS_ROOT

for ent in `find $pdb_root -name \*.ent`
do
  pdbfile=`basename $ent`
  pdbcode=`expr substr $pdbfile 4 4`
  cp $ent $pdbfile
  dssp $pdbfile > ${pdbcode}.dssp
  ${TOPS_ROOT}/bin/Tops $pdbcode
  mv ${pdbcode}.tops ${outdir}/${pdbcode}.tops
  rm ${pdbcode}.dssp
  rm $pdbfile
done
