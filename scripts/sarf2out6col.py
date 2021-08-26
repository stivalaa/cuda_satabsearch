#!/usr/bin/env python
#
# File:    sarf2out6col.sh
# Author:  Alex Stivala
# Created: April 2010
#
# sarf2out6col.sh - Convert SARF2  output format to 6 column format
#                   queryid dbid querylen dblen Nres RMSD
#                  
#
# Usage: sarf2out6col.sh < sarf2_results_file
#
# The input file is read fomr stdin
# Output is to stdout.
#
#
# This is quite tricky since it is not really space-delimited, but (sort of)
# fixed column where fields can end up too large so no space between them,
# and yet don't always START in same column.
#
# $Id: sarf2out6col.py 3695 2010-05-18 07:32:14Z alexs $
#


import sys,os


if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))


querypdbid = None


for line in sys.stdin:
    queryid = line[:8]
    if line[11] == "*": # sometimes we get '***' for some reason
       querylen = 0
    else:
       querylen = int(line[11:14].lstrip())
    dbid = line[16:24]
    if line[27] == "*": # sometimes we get '***' for some reason
       dblen = 0
    else:
       dblen = int(line[27:30].lstrip())
    nres = int(line[32:35].lstrip())
    if line[35] == "*": # sometimes we get '***' for some reason
        sys.stderr.write('WARNING: bad RMSD for %s - %s\n' % (queryid,dbid))
        rmsd = 9999
    else:
        rmsd = float(line[35:40].lstrip())
    sys.stdout.write("%s %s %3d %3d %3d %6.2f\n"
                     % (queryid, dbid, querylen, dblen, nres, rmsd))
    
