#!/usr/bin/env python
###############################################################################
#
# scomdominfo.py - Report information folds and classes of a list of SCOP sids
#
# File:    scomdominfo.py
# Author:  Alex Stivala
# Created: November 2008
#
# $Id: scopdominfo.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Report information on the folds, superfamilies and classes of a list
of SCOP domain identifiers (sids).

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


def write_scopdom_info(scopsid_list, fh, scop):
   """
   Write information about the list of SCOP sids (domain identifiers)
   in the scopsid_list to fh. For each domain write the fold and class,
   then write stats about number of different folds represented 
   and the number of domains in each class.

   Parameters:
      scopsid_list - list of SCOP sids (domain ids)
      fh - open (write) filehandle to write to
      scop - previously built Bio.SCOP Scop instance

   Return value:
      None.
   """
   superfamily_count = {} # dict of {sf_sunid : count} counting domains in eac superfamily
   fold_count= {} # dict of {fold_sunid : count} counting domains in each fold
   class_count={} # dict of {class_sunid : count} counting domains in each class
   for sid in scopsid_list:
       scop_dom = scop.getDomainBySid(sid)
       scop_superfamily = scop_dom.getAscendent('superfamily')
       scop_fold = scop_dom.getAscendent('fold')
       scop_class = scop_dom.getAscendent('class')
       if superfamily_count.has_key(scop_superfamily.sunid):
          superfamily_count[scop_superfamily.sunid] += 1
       else:
          superfamily_count[scop_superfamily.sunid] = 1
       if fold_count.has_key(scop_fold.sunid):
           fold_count[scop_fold.sunid] += 1
       else:
           fold_count[scop_fold.sunid] = 1
       if class_count.has_key(scop_class.sunid):
           class_count[scop_class.sunid] += 1
       else:
           class_count[scop_class.sunid] = 1
       fh.write('%s\t(%s) %s\t%s\t%s\n' % (sid,  scop_superfamily.sccs,scop_superfamily.description, scop_fold.description, scop_class.description))
       
   num_domains = len(scopsid_list)
   num_superfamilies = len(superfamily_count)
   num_folds = len(fold_count)
   num_classes = len(class_count)
   fh.write('Totals: %d domains\t%d superfamilies\t%d folds\t%d classes\n' % 
             (num_domains, num_superfamilies, num_folds, num_classes))
   fh.write('Class distribution:\n')
   for (class_sunid, count) in class_count.iteritems():
       fh.write('\t%s:\t%d\n' % (scop.getNodeBySunid(class_sunid).description,
                               count))


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
                     " < domainidlist\n")
    sys.exit(1)

    
def main():
    """
    main for scomdominfo.py

    Usage: scomdominfo.py  < domainidlist

    
    The list of SCOP domain ids (sids) is read from stdin
    Output is written to stdout.
    """
    if len(sys.argv) != 1:
        usage(os.path.basename(sys.argv[0]))


    # read SCOP data
    scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)

    scopsid_list = sys.stdin.read().split('\n')[:-1]
    write_scopdom_info(scopsid_list, sys.stdout, scop)

            
if __name__ == "__main__":
    main()
