#!/bin/bash
###############################################################################
#
# mergeoutput.sh - merge two sets of tableau search output into single file
#
# File:    mergeoutput.sh
# Author:  Alex Stivala
# Created: December 2009
#
#
# For all the .out files in outdir1 and outdir2, assumed to be named
# with as queryid.out e.g. d4ubpb_.out and with output in two-column
#
#   dbid score
#
# format, create a single merged output with 3 columns
#
#   queryid dbid score1 score2
#
# where score1 and score2 are score for matching queryid with dbid
# according to output in dir1 and dir2 respectively.
# Used for large-scale comparison of differences in scores between two
# methods (in R etc.)
#
# Output is to stdout.
#
# Usage:
#     mergeoutput.sh outdir1 outdir2
#
# $Id: mergeoutput.sh 3007 2009-12-04 06:31:37Z alexs $
# 
###############################################################################

if [ $# -ne 2 ]; then
    echo Usage: $0 outdir1 outdir2  >&2
    exit 1
fi

dir1=$1
dir2=$2

for qpfile in ${dir1}/*.out
do
    queryid=`basename $qpfile .out`
    safile=${dir2}/${queryid}.out
    sortedqpfile=`mktemp`
    sortedsafile=`mktemp`
    grep -v '^#' $qpfile | sort > $sortedqpfile
    grep -v '^#' $safile | sort > $sortedsafile
    join $sortedqpfile $sortedsafile  \
            | grep -v ERROR           \
            | awk "{print \"$queryid\", \$0}" 
    rm $sortedqpfile 
    rm $sortedsafile
done


