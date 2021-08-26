#!/bin/bash
###############################################################################
#
# get_superfamily_list.sh - Get list of superfamilies in top scoring hits
#
# File:    get_superfamily_list.sh
# Author:  Alex Stivala
# Created: March 2009
#
# For the top hits according to structural search output
# on stdin (in two column name score format), get list of unique superfamilies.
#
# Usage:
#     get_superfamily_list.sh [-n num_hits] 
#
#     -n num_hits : number of top hits to use (default all)
#
# Input is tsrchd output file on stdin.  
# Output is list of SCOP superfamily sccs identifier (e.g. d.15.1) on stdout.
#
# Environment variables:
#
#   PATH must contain the location of the Python scripts, ie where this
#   script itself is and the ptgraph/ directory with scopdominfo.py etc.
#
#   PYTHONPATH must contain the directory containing the ptsecstruct.py 
#   and other Python modules used by the Python scripts.
#
#
# Uses GNU utilities options (and bash)
#
# Relies on the output format generated by scopdominfo.py
#
# $Id: get_superfamily_list.sh 2104 2009-03-16 06:47:45Z astivala $
# 
###############################################################################

num_hits=0

while getopts 'n:' opt
do
    case $opt in
    n)
       num_hits="$OPTARG"
       ;;
    ?)
    echo "Usage: $0  [-n num_hits]" >&2
    exit 1
    ;;
    esac
done
shift $(($OPTIND - 1))


if [ $# -ne 0 ]; then
    echo "Usage: $0 [-n num_hits" >&2
    exit 1
fi

if [ $num_hits -gt 0 ]; then
    filter="head -n${num_hits}"
else
    filter="cat"
fi

grep -v '^#' | sort -k2,2n | $filter | cut -d' ' -f1 | scopdominfo.py |  \
    grep '^d' | cut -d\( -f2 | cut -d \) -f1 | sort | uniq