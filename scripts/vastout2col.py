#!/usr/bin/env python
#
# File:    vastout2col.sh
# Author:  Alex Stivala
# Created: November 2008
#
# vastout2col.sh - Convert VAST .gibbs output format to 2 column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: vastout2col.sh < domain.gibbs
# 
# Output has two columns, database id and VAST Pcli score
#
# Output is to stdout.
#
# Uses the output format from  VAST (Gibrat et al 1996; Madej et al 1995),
# available from
# http://migale.jouy.inra.fr/outils/mig/vast/
#
# $Id: vastout2col.py 3603 2010-05-04 04:47:51Z alexs $
#

import os,sys
from itertools import groupby

value_header = False
dbid = None

scorelist = []  # list of (targetpdbid,Pcli) tuples

for line in sys.stdin:
    splitline = line.split()
    if len(splitline) > 1 and splitline[1] == 'Nclique=':
        dbid = splitline[0]
        value_header = False
    elif splitline[0] == 'Nres' and splitline[6] == 'Pcli':
        value_header = True
    elif value_header:
        Pcli = splitline[6]
        scorelist.append((dbid, Pcli))
        value_header = False

# for reasons I don't entirely understand 
# there are sometimes two or more entries for the same target
# with differing Pcli and other values. We will also choose the one
# with highest Pcli
single_scorelist = []
targetpdbid_group_iter = groupby(sorted(scorelist), lambda t : t[0])
for (targetpdbid, targetpdbid_iter) in targetpdbid_group_iter:
    maxPcli = max([Pcli for (pdbid,Pcli) in targetpdbid_iter])
    single_scorelist.append((targetpdbid, maxPcli))


for (targetpdbid, Pcli) in single_scorelist:
    sys.stdout.write('%s    %s\n' % (targetpdbid, Pcli))

