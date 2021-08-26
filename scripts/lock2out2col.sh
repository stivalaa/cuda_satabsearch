#!/bin/sh
#
# File:    LOCK2out2col.sh
# Author:  Alex Stivala
# Created: March 2010
#
# locks2out2col.sh - Convert LOCK2 default output format to same
#                    format as output by tsrchd_sparse etc. which can
#                    be processed with tsevalfn.py etc.
#
# Usage: LOCK2out2col.sh [-q] < LOCK2output.LOCK2.out 
#
# -q : add QUERY ID = line
# 
# Output has two columns, database id and LOCK2 score.
# The query id is put in a comment line at top of file, it is assumed
# to be the same in every line of LOCK2 -A output since that mode
# runs one query against a db.
#
# Output is to stdout.
#
# Uses the output format from LOCK2 in the FoldMiner package
# http://motif.stanford.edu/distributions/foldminer/FoldMinerDistribution.tar.gz
#
#
# $Id: lock2out2col.sh 3522 2010-03-24 05:50:47Z alexs $
#

outputqueryid=0
if [ $# -gt 1 ]; then
   echo "usage: $0 [-q] < LOCK2output" >&2
   exit 1
elif [ $# -eq 1 ]; then
  if [ $1 = "-q" ]; then
    outputqueryid=1
  else
    echo "usage: $0 [-q] < LOCK2output" >&2
    exit 1
  fi
fi

awk -v outputqueryid=$outputqueryid '
     /^\*\* Target =/ { path = $4;
                        splitlen = split(path, splitpath, "/");
                        entname = splitpath[splitlen];
                        target = substr(entname, 1, 7)
                      }
     /^final score:/ { score = $3;
                        print target, score;
                      } 
     /^\*\* Query =/ { path = $4;
                        if (outputqueryid == 1 && !donequery) {
                          splitlen = split(path, splitpath, "/");
                          entname = splitpath[splitlen];
                          queryid = substr(entname, 1, 7)
                          printf("# QUERY ID = %s\n", queryid);
                          donequery = 1;
                        }
                      }

    '
