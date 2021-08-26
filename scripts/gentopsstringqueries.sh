#!/bin/bash
###############################################################################
#
# gentopsstringqueries.sh - generate TOPS string queries
#
# File:    gentopsstringqueries.sh
# Author:  Alex Stivala
# Created: March 2009
#
#
# Generate TOPS strings for each of the 8 
# structures we use for testing. 
#
#
# Usage:
#     gentopsstringqueries.sh outdir
#
#     outdir is the directory to write the query files to. It is created
#     if it does not exist. WARNING: files overwritten if they do exist.
#
# $Id: gentopsstringqueries.sh 2122 2009-03-23 22:39:40Z astivala $
# 
###############################################################################


# location of TOPS directory, contains tops.def etc.
# Note all the .dssp and .tops files are temporarily created here,
# (tops.def has these specifications)
TOPS_ROOT=/local/charikar/astivala/biosoftware/Tops

# location of tops_comparison directory, contains jars/translation.jar etc.
TOPS_COMPARISON_ROOT=/local/charikar/astivala/biosoftware/tops_comparison

# Root of ASTRAL divided PDB style hierarchy
ASTRAL_ROOT=/local/charikar/ASTRAL/pdbstyle-1.73

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

cd $TOPS_ROOT

for struct in $STRUCTS
do
  div=`echo $struct | cut -c3-4`
  pdbfile=${ASTRAL_ROOT}/${div}/${struct}.ent
  # TOPS can only cope with 4 letter PDB codes, so we have to name 
  # input files that way
  sid=$struct
  pdbcode=`echo $sid | cut -c2-5`
  cp ${pdbfile} pdb${pdbcode}.ent
  dssp ${pdbfile} > ${pdbcode}.dssp
  ${TOPS_ROOT}/bin/Tops $pdbcode
  topsfile=${outdir}/${sid}.tops
  mv ${pdbcode}.tops ${topsfile}
  rm ${pdbcode}.dssp
  rm pdb${pdbcode}.ent

  java -cp ${TOPS_COMPARISON_ROOT}/jars/translation.jar tops.translation.Tops2String $topsfile $sid > ${outdir}/${sid}.topsstring
done

