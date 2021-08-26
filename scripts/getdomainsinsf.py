#!/usr/bin/env python
###############################################################################
#
# getdomainsinsf.py - Get a domain in each superfamily specfied on stdin.
#
# File:    getdomainsinsf.py
# Author:  Alex Stivala
# Created: February 2010
#
# $Id: getdomainsinsf.py 3322 2010-02-11 05:46:13Z alexs $
# 
###############################################################################

"""
For each SCOP sccs superfamily string (e.g. 'd.58.1') read from stdin,
output a domain in the 95% nr ASTRAL subset.

See usage in docstring for main()

SCOP and ASTRAL data is obtained using the Bio.SCOP library (Casbon et
al 2006 'A high level interface to SCOP and ASTRAL implemented in
Python' BMC Bioinformatics 7:10) and depends on having the data
downloaded, in SCOP_DIR (defined below).

Downloaded SCOP files from

http://scop.mrc-lmb.cam.ac.uk/scop/parse/index.html

and ASTRAL files (in scopseq-1.73) from

http://astral.berkeley.edu/scopseq-1.73.html

The files downlaoded are:

/local/charikar/SCOP/:
dir.cla.scop.txt_1.73
dir.des.scop.txt_1.73
dir.hie.scop.txt_1.73

/local/charikar/SCOP/scopseq-1.73:
astral-scopdom-seqres-all-1.73.fa
astral-scopdom-seqres-sel-gs-bib-95-1.73.id

Other files there are indices built by Bio.SCOP when first used.
"""

import sys,os
import getopt
import random
 
from Bio.SCOP import *

from pathdefs import SCOP_DIR,SCOP_VERSION


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def get_domain_for_each_sf(sccs_list, scop, astral):
   """
   For each superfamily named by sccs in the sccs_list,
   return a domain sid in that
   superfamily in the 95% nr ASTRAL subset.

   Parameters:
      sf_list - list of Bio.SCOP superfamily objects
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance

   Return value:
      list of SCOP sids, one for each superfamily.
   """
   
   # Bio.SCOP actually doesn't seem to have a facility to look up by
   # sccs so we'll build a dictionary ourselves of all superfamilies
   # keyed by sccs
   all_superfamilies = scop.getRoot().getDescendents('sf')
   sccs_dict = dict([(sf.sccs, sf) for sf in all_superfamilies])

   domain_sids = []
   for sccs in sccs_list:
      sf = sccs_dict[sccs]
      domain_list = [ dom for dom in sf.getDescendents('domain')
                      if astral.isDomainInId(dom, 95) ]
#      sys.stderr.write('xxx ' + str(domain_list))
      if len(domain_list) > 0:
         domain = random.choice(domain_list)
         domain_sids.append(domain.sid)
         
   return domain_sids

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + "\n")
    sys.exit(1)

    
def main():
    """
    main for getdomainsinsf.py

    Usage: getdomainsinsf.py 


    List of SCOP superfamily ids (sccs) is read from stdin.
    The list of SCOP domain ids (sids) is printed to stdout.
    """
    global verbose
    verbose = False
    
    use_nonredundant = False

    if len(sys.argv) != 1:
       usage(os.path.basename(sys.argv[0]))
       

    # read SCOP and ASTRAL data
    scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)
    astral = Astral(dir_path=SCOP_DIR,version=SCOP_VERSION,scop=scop)
       
    sccs_list = sys.stdin.read().split('\n')[:-1]
    sid_list = get_domain_for_each_sf(sccs_list, scop, astral)
    for sid in sid_list:
       sys.stdout.write(sid + '\n')
    
            
if __name__ == "__main__":
    main()
