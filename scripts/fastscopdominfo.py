#!/usr/bin/env python
###############################################################################
#
# fastscopdominfo.py - Report information folds and classes of a list of SCOP sids
#
# File:    fastscopdominfo.py
# Author:  Alex Stivala
# Created: March 2010
#
# $Id: fastscopdominfo.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Report information on the folds and superfamilies and classes of a list
of SCOP domain identifiers (sids).
scopdominfo.py does this from Bio.SCOP, but this is quite slow to load
so this version uses a cached (pickle dictionary) table built by
build_fastscopdominfo_cache.py

See usage in docstring for main()
"""

import sys,os
import pickle


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def write_dom_info(scopsid_list, scopdominfo_dict):
   """
   Write information about the list of SCOP sids (domain identifiers)
   in the scopsid_list to fh. For each domain write the superfamily sccs,
   superfamily description and fold description.
   Delimiter is '|'


   Parameters:
     scopsid_list - list of SCOP sids (domain ids)
     scopdominfo_dict -
      dict {sid: (superfamily_sccs, superfamily_description, fold_sccs, fold_description)}
      where
     superfamily_sccs is SCOP sccs identifying the superfamily for the domain
     superamily_description is SCOP dessription of the superfamily
     fold_description is the SCOP descriptino of the fold the domain is in
   Return value:
      None.
   """
   for sid in scopsid_list:
       entry = scopdominfo_dict[sid]
       sf_sccs = entry[0]
       sf_desc = entry[1]
       fold_sccs =entry[2]
       fold_desc = entry[3]
       sys.stdout.write("%s | %s | %s | %s\n" %(sid, sf_sccs, sf_desc, fold_desc))

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + 
                     "cachefile < domainidlist\n")
    sys.exit(1)

    
def main():
    """
    main for fastscopdominfo.py

    Usage: fastscopdominfo.py cachefile  < domainidlist

    cachefile is the filename of the cache (pickled) file built by
    build_fastscopdominfo_cache.py
    
    The list of SCOP domain ids (sids) is read from stdin
    Output is written to stdout.
    """
    if len(sys.argv) != 2:
        usage(os.path.basename(sys.argv[0]))

    pickle_filename = sys.argv[1]
    scopdominfo_dict = pickle.load(open(pickle_filename))
    scopsid_list = sys.stdin.read().split('\n')[:-1]
    write_dom_info(scopsid_list, scopdominfo_dict)

            
if __name__ == "__main__":
    main()
