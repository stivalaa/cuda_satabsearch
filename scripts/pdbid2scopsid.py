#!/usr/bin/env python
###############################################################################
#
# pdbid2scopsid.py - Convert PDB identifier with chain id to SCOP sid
#
# File:    pdbid2scopsid.py
# Author:  Alex Stivala
# Created: March 2010
#
# $Id: pdbid2scopsid.py 3471 2010-03-15 23:13:09Z alexs $
# 
###############################################################################

"""
Usage: pdbid2scopsid.py < pdbidlist

Convert a list of PDB identifiers (one per line)
with chain identifiers as used by DaliLite
e.g. 1qlpA to SCOP sid e.g. d1qlpa_

The reverse (SCOP sid to PDB id with chain) can be done with just

 sed 's/.$/\U&/g'

but this is more complicated since we need to get the right character
on the end (_ or a domain id) and also in the latter case there
is ambiguity that we need
to resolve e.g. 1u6rA could be either d1u6ra1 or d1u6ra2.
A further complication arises with the 'genetic domain' type ASTRAL
structures e.g. d1mtp.1 (g1mtp.1 in the sequence and id SCOP files, another
complication (see the diff to Bio/SCOP/__init__.py for this)
which is 1mtpA and 1mtpB in the PDB DaliLite scheme (2 chains, both in
same domain).
"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt

 
from Bio.SCOP import *

from tsevalutils import filter_domains_astral_nrpercent

from pathdefs import SCOP_DIR,SCOP_VERSION


def pdbid_to_scopsid(pdbid, all_scopsids_dict):
    """
    Convert a PDB id with chain to SCOP sid, as per description in
    module header docstring

    Parameters:
       pdbid - pdb id with chain e.g. 1u6rA
       all_scopsids_dict - dict (scopsid, True) of all SCOP sids to check

    Return value:
       scop sid e.g. d1u6ra2
    """
    scopsid = 'd' + pdbid.lower() + '_'
    if not all_scopsids_dict.has_key(scopsid):
        sid_list = [sid for sid in all_scopsids_dict.keys() if
                    sid[1:6] == scopsid[1:6]]
        if len(sid_list) < 1:
            sid_list = [sid for sid in all_scopsids_dict.keys() if
                        sid[1:5] == scopsid[1:5]]
        if len(sid_list) > 0:
            scopsid = sorted(sid_list)[0] # always take lowest domain id FIXME
        else:
            scopsid="UNKNOWN"
    return scopsid

    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + 
                     " < pdbidlist\n")
    sys.exit(1)


def main():
    """
    main for pdbid2scopsid
    see usage in module docstring
    """

    
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

    
    for pdbid in sys.stdin:
        scopsid = pdbid_to_scopsid(pdbid.rstrip(), all_scopsids_dict)
        sys.stdout.write(scopsid + '\n')

if __name__ == "__main__":
    main()
