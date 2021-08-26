#!/bin/sh
#
# getrank.sh - get rank of a hit from the output file
#
# $Id: getrank.sh 1790 2008-08-04 02:17:59Z astivala $
#
#    Usage: getrank.sh  query_sid  tabsearch.outputfile
#
#    query_sid is the SCOP id (e.g. 'd1ubia_') of the query domain
#    
#    tabsearch.outputfile is the output from tabsearchqpml / tsrchn / tsrchd
#    which is a text file where each line is identifier then whitespace
#    then score, sorted by score from most negative to least negative e.g.
#
#    d1xksa_ -35.99999999
#    d3sila_ -35.99999999
#    ....
#    d2mhua_ -0.499999999
#
#    ie this means the best hit as the top of the file, and worst at bottom.
#
# Output is rank in file as fraction of line number/total line and percentage
# e.g. 
#
# 944/15174 (6%)
#
# Uses GNU sort options.

if [ $# -ne 2 ]; then
    echo "Usage: $0 query_sid  tabsearch_outputfile" >&2
    exit 1
fi


query_sid=$1
resultsfile=$2
totallines=`grep -v '^#' ${resultsfile} | wc -l`
line=`grep -v '^#' ${resultsfile} | sort -k2,2n | grep -n ${query_sid} | cut -d: -f1`
percentage=`echo "(${line} / ${totallines}) * 100" | bc -l`
printf "%d/%d (%.0f%%)\n" $line $totallines $percentage

