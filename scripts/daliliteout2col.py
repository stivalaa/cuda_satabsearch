#!/usr/bin/env python
#
# File:    daliliteout2col.sh
# Author:  Alex Stivala
# Created: March 2010
#
# daliliteout2col.sh - Convert DaliLite .dccp output format to 2 column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: daliliteout2col.sh < dccpfile 
#
# The input file is read fomr stdin
#
# Output has two columns, database id and DaliLite Z-score
#
# Output is to stdout.
#
#
# We also convert pdbid with chain format e.g. 1atnA to format
# for Fischer evaluation e.g. 1atn_a 
# and those that did not have chain specified in Fischer data set,
# do now have for DaliLite, so remove them (so we need the Fischer
# pdb list)
#
# $Id: daliliteout2col.py 3575 2010-04-22 00:35:26Z alexs $
#


import sys,os
from itertools import groupby

from fischer_tables import FISCHER_ID_FOLD_DICT

def daliid_to_fischerid(daliid):
    """
    Convert a DaliLite id with chain e.g. 1atnA
    to format for Fischer dat set e.g. 1atn_a, those without
    chains in Fischer data set
    e.g. 1cew have chain removed i.e. 1cewA becomes 1cew
    NB all inputs (Dali ids) have chain

    Parameters:
       daliid -DaliLite identifier, with chain on end (no dliemiter)

    Return value:
       PDB identifier with chain after _
    """
    if FISCHER_ID_FOLD_DICT.has_key(daliid[:4]):
        return daliid[:4]
    else:
        return daliid[:4] + '_' + daliid[4]

def usage(progname):
    sys.stderr.write("Usage: " + progname + " < dccpdata\n")
    sys.exit(1)


if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))


querypdbid = None


scorelist = []  # list of (targetpdbid,zscore) tuples
for line in sys.stdin:
    splitline = line.split()
    if len(splitline) > 0 and splitline[0] == 'DCCP':
        if len(splitline) == 10:
            targetpdbid = splitline[9]
            zscore = splitline[5]
            if querypdbid == None:
                querypdbid = splitline[8]
        else: # sometimes fields 2 and 3 get stuck together
            targetpdbid = splitline[8]
            zscore = splitline[4]
        scorelist.append((targetpdbid, zscore))

# for reasons I don't understand (and can't find any documentation on the
# dccp file format) there are often two or more entries for the same target
# with differing Z-scores and other values. We will also choose the one
# with highest Z-score
single_scorelist = []
targetpdbid_group_iter = groupby(sorted(scorelist), lambda t : t[0])
for (targetpdbid, targetpdbid_iter) in targetpdbid_group_iter:
    maxzscore = max([zscore for (pdbid,zscore) in targetpdbid_iter])
    single_scorelist.append((daliid_to_fischerid(targetpdbid), maxzscore))


querypdbid = daliid_to_fischerid(querypdbid)

sys.stdout.write('# QUERY ID = ' + querypdbid + '\n')
for (targetpdbid, zscore) in single_scorelist:
    sys.stdout.write('%s    %s\n' % (targetpdbid, zscore))

