#!/usr/bin/env python
###############################################################################
#
# scopsuperfamilyinfo.py - Report information folds and classes of a
#                          list of SCOP sccs identifiers
#
# File:    scopsuperfamilyinfo.py
# Author:  Alex Stivala
# Created: March 2009
#
# $Id: scopsuperfamilyinfo.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Report information on the folds, superfamilies and classes of a list
of SCOP superfamily identifiers in the form of sccs (SCOP Concise
Classification Strings) identifeirs
(e.g. d.58.1)

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
 
from Bio.SCOP import *

from pathdefs import SCOP_DIR,SCOP_VERSION


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def write_scopsuperfamily_info(scopsccs_list, fh, scop):
   """
   Write information about the list of SCOP sccs ids
   in the scopsccs_list to fh. For each sccs id write the
   superfamily and fold description.

   Parameters:
      sccsid_list - list of SCOP superfamilies (sccs ids)
      fh - open (write) filehandle to write to
      scop - previously built Bio.SCOP Scop instance

   Return value:
      None.
   """
   # Bio.SCOP actually doesn't seem to have a facility to look up by
   # sccs so we'll build a dictionary ourselves of all superfamilies
   # keyed by sccs
   all_superfamilies = scop.getRoot().getDescendents('sf')
   sccs_dict = dict([(sf.sccs, sf) for sf in all_superfamilies])

   for sccs in scopsccs_list:
       sf = sccs_dict[sccs]
       fold = sf.getAscendent('fold')
       fh.write('%s\t%s\t%s\n' % (sf.sccs, sf.description, fold.description))
   

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
                     " < sccslist\n")
    sys.exit(1)

    
def main():
    """
    main for scopsuperfamilyinfo.py

    Usage: scopsuperfamilyinfo.py  < sccslist

    
    The list of SCOP sccs superfamily strings (e.g. 'd.58.1')
    is read from stdin
    Output is written to stdout.
    """
    if len(sys.argv) != 1:
        usage(os.path.basename(sys.argv[0]))


    # read SCOP data
    scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)

    sccs_list = sys.stdin.read().split('\n')[:-1]
    write_scopsuperfamily_info(sccs_list, sys.stdout, scop)

            
if __name__ == "__main__":
    main()
