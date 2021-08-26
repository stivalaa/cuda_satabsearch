#!/usr/bin/env python
#
# File:    dalilitequery2col.sh
# Author:  Alex Stivala
# Created: March 2010
#
# dalilitequery2col.sh - Convert DaliLite .dccp output format to 2 column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: dalilitequery2col.sh < querysidlist
#
# The input file is read as pdbid.dccp in cwd
# where the pdbid is constructed from
# the query id SCOP sid e.g. d1u6ra2 becomes 1u6rA, for each query sid
# in the querysidlist (one sid per line) on stdin.
#
# Output has two columns, database id and DaliLite Z-score
#
# Output is to stdout.
#
# $Id: dalilitequery2col.py 3473 2010-03-16 04:06:23Z alexs $
#


import sys,os
from itertools import groupby

from Bio.SCOP import *

from tsevalutils import filter_domains_astral_nrpercent
from pathdefs import SCOP_DIR,SCOP_VERSION
from pdbid2scopsid import pdbid_to_scopsid


def usage(progname):
    sys.stderr.write("Usage: " + progname + " < querysidlist\n")
    sys.exit(1)


if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))



# read SCOP and ASTRAL data
sys.stderr.write('reading SCOP data...\n')
scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)
astral = Astral(dir_path=SCOP_DIR,version=SCOP_VERSION,scop=scop)
nrpercent = 95 # Always use 95% nr subset. TODO make this an option
all_domains = scop.getRoot().getDescendents('domain')
if nrpercent != None:
    all_domains = filter_domains_astral_nrpercent(all_domains,
                                                  scop, astral,
                                                  nrpercent)
all_scopsids_dict = dict( [(d.sid,True) for d in all_domains] )


for querysid in sys.stdin:
    querysid = querysid.rstrip()
    querypdbid = querysid[1:5] + querysid[5].upper()
    dccpfilename = querypdbid + ".dccp"

    scorelist = []  # list of (targetsid,zscore) tuples
    for line in open(dccpfilename):
        splitline = line.split()
        if len(splitline) > 0 and splitline[0] == 'DCCP':
            if len(splitline) == 10:
                targetpdbid = splitline[9]
                zscore = splitline[5]
            else: # sometimes fields 2 and 3 get stuck together
                targetpdbid = splitline[8]
                zscore = splitline[4]
            targetsid = pdbid_to_scopsid(targetpdbid, all_scopsids_dict)
            scorelist.append((targetsid, zscore))

    # for reasons I don't understand (and can't find any documentation on the
    # dccp file format) there are often two or more entries for the same target
    # with differing Z-scores and other values. We will also choose the one
    # with highest Z-score
    single_scorelist = []
    targetsid_group_iter = groupby(sorted(scorelist), lambda t : t[0])
    for (targetsid, targetsid_iter) in targetsid_group_iter:
        maxzscore = max([zscore for (sid,zscore) in targetsid_iter])
        single_scorelist.append((targetsid, maxzscore))


    sys.stdout.write('# QUERY ID = ' + querysid + '\n')
    for (targetsid, zscore) in single_scorelist:
        sys.stdout.write('%s    %s\n' % (targetsid, zscore))

