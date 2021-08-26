#!/usr/bin/env python
#
# File:    sheba2fischerid.py
# Author:  Alex Stivala
# Created: April 2010
#
# sheba2fischerid.py - Convert SHEBA id always with chain to fischer id
#
# Usage: sheba2fischerid.py < shebamultiqeryout
#
# The input file is 2 column from shebaout2col (multiple query) on stdin
# Output is to stdout.
#
#
# We convert pdbid always with chain format e.g. 8ilb_a to format
# for Fischer evaluation e.g. 8ilb for those that do not have chain
# spcified in Fischer data set, otherwise leave the chain on.
#
# $Id: sheba2fischerid.py 3596 2010-05-03 02:38:39Z alexs $
#


import sys,os
from itertools import groupby

from tsevalutils import iter_searchresult
from fischer_tables import FISCHER_ID_FOLD_DICT

def shebaid_to_fischerid(shebaid):
    """
    We convert pdbid always with chain format e.g. 8ilb_a to format
    for Fischer evaluation e.g. 8ilb for those that do not have chain
    spcified in Fischer data set, otherwise leave the chain on.
    NB all inputs have chain

    Parameters:
       sheba -sheba identifier, with chain on end after underscore

    Return value:
       PDB identifier with chain after _ only if chain specified in Fischer
    """
    if FISCHER_ID_FOLD_DICT.has_key(shebaid[:4]):
        return shebaid[:4]
    else:
        return shebaid

def usage(progname):
    sys.stderr.write("Usage: " + progname + " < sheba2colout\n")
    sys.exit(1)


if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))


query_groupby_iter = groupby(sorted(iter_searchresult(sys.stdin,multiquery=True) ), lambda t : t[0])
for (queryid, result_iter) in query_groupby_iter:
    sys.stdout.write("# QUERY ID = %s\n" % shebaid_to_fischerid(queryid))
    for (queryid, score,dbid) in result_iter:
        sys.stdout.write("%s    %d\n" % (shebaid_to_fischerid(dbid), score))

