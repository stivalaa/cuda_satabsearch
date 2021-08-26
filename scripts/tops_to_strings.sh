#!/bin/sh
#
# File:    tops_to_strings.sh
# Author:  Alex Stivala
# Created: March 2009
#
# tops_to_strings.sh - build TOPS strings from database of TOPS files
#
# Usage: tops_to_strings.sh tops_db_dir 
#
#   top_db_dir is the name of the directory containing TOPS file
#   (as built with build_tops_files)
#
# The TOPS strings, one per line, are written to stdout
#
# $Id: tops_to_strings.sh 3689 2010-05-18 00:52:29Z alexs $

# location of tops_comparison directory, contains jars/translation.jar etc.
TOPS_COMPARISON_ROOT=$HOME/tops_comparison

if [ $# -ne 1 ]; then
    echo "Usage: $0 tops_db_dir" 2>&1
    exit 1
fi

tops_db_dir=$1


for topsfile in `find $tops_db_dir -maxdepth 1 -name \*.tops`
do
    #sid=`basename $topsfile .tops`
    sid=`basename $topsfile .pdb.tops`
    pdbcode=`echo $sid | cut -c2-5`
    # for some reason even though TOPS only handles PDB ids e.g. 1NDD not
    # SCOP sids e.g. d1ndda_, Tops2String crashes when when the PDB id
    # is used, regardless of what we use as the 'string name' on the command
    # line  so we have to replace it in the .tops file with a SCOP sid 
    # e.g. 
    # DOMAIN_NUMBER 0 1ndd 1 1 74
    # becomes
    # DOMAIN_NUMBER 0 d1ndda_ 1 1 74
    # then we can use anythnign as the 'string name' but it gets junk added
    # to the end, which can can anyway remove later (topscompreout2col.sh).
    temp_tops_file=/var/tmp/topstmp.$$.tops
    sed "s/DOMAIN_NUMBER 0 [a-z0-9.]* \(.*\) \(.*\) \(.*\)/DOMAIN_NUMBER 0 ${sid} \1 \2 \3/" < ${topsfile} > ${temp_tops_file}
    java -cp ${TOPS_COMPARISON_ROOT}/jars/translation.jar tops.translation.Tops2String $temp_tops_file $sid
    rm ${temp_tops_file}
done
