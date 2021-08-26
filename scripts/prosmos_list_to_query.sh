#!/bin/sh
#
# File:    prosmos_list_to_query.sh
# Author:  Alex Stivala
# Created: March 2009
#
# prosmst_list_to_query.sh - convert ProSMoS .list file to first go at query
#
# Usage: prosmost_list_to_query.sh  listfile 
#
# Takes the .list file created by the fetchmatrix.pl script and convert
# it to format for ProSMoS query. All this does is reformatting of the
# data (and removing some) - the .query file MUST be manually edited to
# be useful. (Actually, may be simpler to edit the a copy of list file first,
# to remove unwanted SSEs,
# then run this to convert format)
# 
# Writes query file format to stdout
#
# $Id: prosmos_list_to_query.sh 2138 2009-03-26 23:54:00Z astivala $
#

if [ $# -ne 1 ]; then
    echo "Usage: $0 listfile" 2>&1
    exit 1
fi
listfile=$1
num_sses=`awk '/^[0-9]+[ ]*\*/ {print}' $listfile | wc -l`
i=1
while [ $i -le $num_sses ]; do
    printf "%d" $i
    if [ $i -lt $num_sses ]; then
        printf " "
    else
        printf '\n'
    fi
    i=`expr $i + 1`
done
awk '/^[0-9]+[ ]*[EH]/{print substr($2,1,1)}' $listfile  | tr '\n' ' '
echo
awk '/^[0-9]+[ ]*\*/{rownum=$1 ; row=""; for (i=1;i<=length($2);i++){row = row substr($2,i,1) " "} ;      numspaces = rownum * 2 - 2 ; for(i=0; i<numspaces;i++) printf(" ");  printf("%s\n",row)}' $listfile


