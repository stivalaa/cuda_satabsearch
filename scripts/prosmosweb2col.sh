#!/bin/sh
#
# File:    prosmosweb2col
# Author:  Alex Stivala
# Created: March 2010
#
# prosmosweb2col.sh - Convert ProSMoS web server html results to sam
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: prosmosweb2col queryid prosmoswebout.html 
#
# queryid is the query id to put as QUERY ID = comment line in output
# prosmoswebout.html is the html saved from the ProSMoS results webpage
#   see http://prodata.swmed.edu/ProSMoS/index.html
# 
# Output has two columns, database id and PROSMOS score
#
# Output is to stdout.
#
# $Id: prosmosweb2col.sh 3485 2010-03-17 04:48:13Z alexs $
#

if [ $# -ne 2 ]; then
    echo "Usage: $0 queryid prosmosweboutput.html" >&2
    exit 1
fi

queryid=$1
htmlfile=$2

tmpfile1=/var/tmp/p1$$
tmpfile2=/var/tmp/p2$$

awk 'BEGIN { i = 1 } /^\* [0-9]*/{print i,$3,$NF;i=i+1}' < $htmlfile > $tmpfile1
awk 'BEGIN {i=1} / d[0-9]..... <http/ {print i,$1;i+=1}' < $htmlfile > $tmpfile2

echo "# QUERYID = " $queryid
join $tmpfile1 $tmpfile2  | awk '{print $4,$3}'

rm $tmpfile1 $tmpfile2

