#!/bin/sh
#
# File:    yakusaout2col.sh
# Author:  Alex Stivala
# Created: March 2010
#
# yakusaout2col.sh - Convert yakusa default output format to same
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: yakusaout2col.sh [-q] < yakusaoutput.yakusa.out 
#
# -q : add QUERY ID = line
# 
# Output has two columns, database id and SHEBA m-score.
# The query id is put in a comment line at top of file, it is assumed
# to be the same in every line of yakusa -A output since that mode
# runs one query against a db.
#
# Output is to stdout.
#
# Uses the output format from YAKUSA, see documentation at
# http://bioserv.rpbs.jussieu.fr/Yakusa/download/README_yakusa
# for more information
# E.g:
#
# Protein rank: 1 score:     118.48 Z-score:      24.29 name: d1u6ra1 : 0000 SCOP/ASTRAL domain d1u6ra1 [11960CHAIN A
#
# $Id: yakusaout2col.sh 3397 2010-03-05 04:28:21Z alexs $
#
# Uses GNU head options (-n -1)

outputqueryid=0
if [ $# -gt 1 ]; then
   echo "usage: $0 [-q] < yakusaoutput" >&2
   exit 1
elif [ $# -eq 1 ]; then
  if [ $1 = "-q" ]; then
    outputqueryid=1
  else
    echo "usage: $0 [-q] < yakusaoutput" >&2
    exit 1
  fi
fi

awk -v outputqueryid=$outputqueryid '
     /^Protein rank:/ { score = $7;
                        if (score  == "inf") score = 99999;
                        print $9, score;
                      } 
     /^Description query :/ { if (outputqueryid == 1) 
                              printf("# QUERY ID = %s\n", $7);
                            }
     /^Query: / { printf("# %s\n", $0);} 
     /^Database: / {printf("# %s\n", $0)}'
