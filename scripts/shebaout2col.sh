#!/bin/sh
#
# File:    shebaout2col.sh
# Author:  Alex Stivala
# Created: Novermber 2008
#
# shebaout2col.sh - Convert sheba -A summary column format to 2-column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: shebaout2col.sh < shebaoutput.A.out 
# 
# Output has two columns, database id and SHEBA m-score.
# The query id is put in a comment line at top of file, it is assumed
# to be the same in every line of sheba -A output since that mode
# runs one query against a db.
#
# Output is to stdout.
#
# Uses the output format from SHEBA 3.1.1, see documentation at
# http://rex.nci.nih.gov/RESEARCH/basic/lmb/mms/sheba.htm 
# for more information
#
# $Id: shebaout2col.sh 2031 2008-11-22 04:57:16Z astivala $
#
# Uses GNU head options (-n -1)


awk '/ pdb1   na       pdb2   nb   id    m   %ma    %mb /,/^$/' | awk 'NR > 1' | head -n -1 | awk 'NR==1 {printf "# QUERYID = %s\n",$1} {printf "%s    %s\n",$3,$6}'
