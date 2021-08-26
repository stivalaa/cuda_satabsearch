#!/usr/bin/env python
#
# File:    sarf2fischerid.py
# Author:  Alex Stivala
# Created: April 2010
#
# sarf2fischerid.py - Convert SARF2 id never with chain to fischer id
#
# Usage: sarf2fischerid.py < sarfmultiqeryout
#
# The input file is 2 column from sarf2_6col_to_score.py
#  (multiple query) on stdin
# Output is to stdout.
#
#
# We convert pdbid never with chain format e.g. 8ilb to format
# for Fischer evaluation e.g. 8ilb_a for those that have chain
# spcified in Fischer data set, otherwise leave the chain off.
#
# $Id: sarf2fischerid.py 3599 2010-05-03 04:22:27Z alexs $
#


import sys,os
from itertools import groupby

from tsevalutils import iter_searchresult
from fischer_tables import FISCHER_ID_FOLD_DICT

# dict of { pdbid : pdbid_chain } where pdbid is each identifier
# in Fischer data set WITHOUT chainid and pdbid_chain has the chainid
# if it is specified in Fischer (oterhwise same as pdbid)

FISCHER_CHAINID_DICT = dict( [ (d[:4], d) for d in  FISCHER_ID_FOLD_DICT.keys() ] )

def sarfid_to_fischerid(sarfid):
    """
    We convert pdbid never with chain format e.g. 8ilb to format
    for Fischer evaluation e.g. 8ilb_a for those that have chain
    spcified in Fischer data set, otherwise leave the chain off
    NB no inputs have chain

    Parameters:
       sarf -sarf identifier, no chain

    Return value:
       PDB identifier with chain after _ only if chain specified in Fischer

    Uses global FISCHER_CHAIN_ID dict
    """
    return FISCHER_CHAINID_DICT[sarfid]


def usage(progname):
    sys.stderr.write("Usage: " + progname + " < sarf2colout\n")
    sys.exit(1)


if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))


query_groupby_iter = groupby(sorted(iter_searchresult(sys.stdin,multiquery=True) ), lambda t : t[0])
for (queryid, result_iter) in query_groupby_iter:
    sys.stdout.write("# QUERY ID = %s\n" % sarfid_to_fischerid(queryid))
    for (queryid, score,dbid) in result_iter:
        sys.stdout.write("%s    %f\n" % (sarfid_to_fischerid(dbid), score))

