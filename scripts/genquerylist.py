#!/usr/bin/env python
###############################################################################
#
# genquerylist.py - Generate a list of SCOP domain ids to query database with
#
# File:    genquerylist.py
# Author:  Alex Stivala
# Created: November 2008
#
# $Id: genquerylist.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Generate a list of SCOP domain ids that are representative for use as
database queries. For a given number of queries to generate, ensure
that each class is represeneted in proportion to its total representatino
in the SCOP database.

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


def generate_query_list(num_queries, use_nonredundant, scop, astral):
   """
   Return a list of randomly chosen num_queries SCOP sids for folds that are
   representative of the distribution of folds in the SCOP classes
   all-alpha, all-beta, alpha+beta, alpha/beta, multidomain, membrane/cell
   surface and small proteins.
   Only one (randomly) chosen domain per fold is returned, i.e. each
   fold in the list is uniquely represented.
   

   Parameters:
      num_queries - number of domain ids to generate
      use_nonredundant - if True, use only from ASTRAL 95% nr subset not
                         whole SCOP.
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance

   Return value:
      list of SCOP sids representing domains.
   """
   class_sunids = [ 46456,   # all alpha
                    48724,   # all beta
                    51349,   # alpha/beta
                    53931,   # alpha+beta
                  ]
#                    56572,   # multi-domain (alpha and beta)
#                    56835,   # membrane and cell surface proteins
#                    56992 ]  # small proteins

   num_folds_per_class = []  # correspdoning to class_sunids, num folds in class
   for class_id in class_sunids:
       scop_class = scop.getNodeBySunid(class_id)
       assert(scop_class.type == 'cl') 
       fold_count = len(scop_class.getDescendents('fold'))
       if verbose:
           sys.stderr.write('%d folds in class %d (%s)\n'  %
                            (fold_count, class_id, scop_class.description))
       num_folds_per_class.append(fold_count)
   total_folds = sum(num_folds_per_class)

   if num_queries > total_folds:
       sys.stderr.write(
           'WARNING: There are only %d folds, num_queries changed to %d\n'
                         % (total_folds, total_folds))
       num_queries = total_folds

   num_folds_per_class = [ int(round((float(n)/float(total_folds)) *
                                      num_queries))
                           for n in num_folds_per_class ]
   if verbose:
       for i in xrange(len(class_sunids)):
           sys.stderr.write('%d folds will be represented for class %d (%s)\n'
           % (num_folds_per_class[i], class_sunids[i], 
              scop.getNodeBySunid(class_sunids[i]).description))

   fold_list = []
   for i in xrange(len(class_sunids)):
        scop_class = scop.getNodeBySunid(class_sunids[i])
        folds = random.sample(scop_class.getDescendents('fold'),
                              num_folds_per_class[i])
        fold_list += folds

   domain_sids = []
   for fold in fold_list:
        if use_nonredundant:
            domain_list = [ dom for dom in fold.getDescendents('domain')
                                if astral.isDomainInId(dom, 95) ]
        else:
            domain_list = fold.getDescendents('domain')
        if len(domain_list) > 0:
            domain = random.choice(domain_list)
            domain_sids.append(domain.sid)

   return domain_sids[:num_queries]

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-vn] "
                     " <number_of_queries>\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -n use ASTRAL 95% nonredundant subset not all SCOP\n')
    sys.exit(1)

    
def main():
    """
    main for genquerylist.py

    Usage: genquerylist.py [-vn]  <number_of_quries>

    
    -v turns on debug output to stderr

    -f use only the ASTRAL 95% sequence identity subset not all of SCOP
    
    The list of SCOP domain ids (sids) is printed to stdout.
    """
    global verbose
    verbose = False
    
    use_nonredundant = False

    try:
        opts,args = getopt.getopt(sys.argv[1:], "vn?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        elif opt == "-n": # use ASTRAL 95% nr subset
            use_nonredundant = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    num_queries  = int(args[0])

    # read SCOP and ASTRAL data
    if verbose:
        sys.stderr.write('Reading SCOP data...\n')
    scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)
    astral = Astral(dir_path=SCOP_DIR,version=SCOP_VERSION,scop=scop)

    sid_list = generate_query_list(num_queries, use_nonredundant, scop, astral)
    for sid in sid_list:
        sys.stdout.write(sid + '\n')
    
            
if __name__ == "__main__":
    main()
