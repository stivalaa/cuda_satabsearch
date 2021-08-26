#!/bin/sh
#
# File:    build_tops_files.sh
# Author:  Alex Stivala
# Created: March 2009
#
# build_tops_files.sh - build TOPS files from hierarchy of ASTRAL files
#
# Usage: build_tops_files.sh astral_root outdir 
#
#   outdir is name of diretory which is created, and a .tops file
#   created as a separate file in that directory, for each .ent file
#   in the ASTRAL SCOP hierarchy.
#
#   astral_root is the root of the ASTRAL SCOP pdbstyle hierarchy.
#
#
# $Id: build_tops_files.sh 3617 2010-05-05 06:13:27Z alexs $

# location of TOPS directory, contains tops.def etc.
# Note all the .dssp and .tops files are temporarily created here,
# (tops.def has these specifications)
TOPS_ROOT=$HOME/Tops

if [ $# -ne 2 ]; then
    echo "Usage: $0 astral_root outdir" 2>&1
    exit 1
fi

astral_root=$1
outdir=$2

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

cd $TOPS_ROOT

for ent in `find $astral_root -name \*.ent`
do
  # TOPS can only cope with 4 letter PDB codes, so we have to name 
  # input files that way
  sid=`basename $ent .ent`
  pdbcode=`echo $sid | cut -c2-5`
  cp $ent pdb${pdbcode}.ent
  dssp pdb${pdbcode}.ent > ${pdbcode}.dssp
  # the -C ALL options ensures all chains in one file, otherwise we get
  # multiple files for genetic domain style e.g. d1pid.1 files (multi chains
  # in 1 domain)
  ${TOPS_ROOT}/bin/Tops -C ALL $pdbcode
  mv ${pdbcode}ALL.tops ${outdir}/${sid}.tops
  rm ${pdbcode}.dssp
  rm pdb${pdbcode}.ent
done
