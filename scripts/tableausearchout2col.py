#!/usr/bin/env python
#
# File:    tableausearchout2col.py
# Author:  Alex Stivala
# Created: March 2009
#
# tableausearchout2col.py - Convert Arun's TableauSearch (TableauComparer)
#                   output format to same 
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: tableausearchout2col.py < tableausearchoutput.A.out 
# 
# Output has two columns, database id and tableau comparer score
#
# Input is TableauComparer output (search.scores) on stdin, e.g.:
#/local/charikar/TableauSearchDB/d1u3ya_.ent.angles   Score-of-comparison:    -149.2
#/local/charikar/TableauSearchDB/d1geea_.ent.angles   Score-of-comparison:    -593.7
#
# Output is to stdout.
#
# $Id: tableausearchout2col.py 2088 2009-03-07 02:30:49Z astivala $
#


import os,sys

for line in sys.stdin:
    splitline = line.split()
    fname = splitline[0]
    dbid = os.path.splitext(os.path.splitext(os.path.basename(fname))[0])[0]
    score = splitline[-1]
    sys.stdout.write('%s    %s\n' % (dbid, score))
