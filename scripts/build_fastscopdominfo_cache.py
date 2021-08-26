#!/usr/bin/env python
###############################################################################
#
# build_fastscopdominfo_cache.py - build pickle file for cached SCOP info
#
# File:    build_fastscopdominfo_cache.py
# Author:  Alex Stivala
# Created: March 2010
#
# $Id: scopdominfo.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Build cache (Python pickled dictionary) of information on the folds
and superfamilies SCOP domain identifiers (sids).

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
import pickle

from Bio.SCOP import *

from pathdefs import SCOP_DIR,SCOP_VERSION

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def build_scopdominfo_dict(scop):
    """
    Build dictionary with
    information about superfamily and class of all SCOP domains
    
    Parameters:
      scop - previously built Bio.SCOP Scop instance
      
    Return value:
     dict {sid: (superfamily_sccs, superfamily_description, fold_sccs,fold_description)}
      where
     superfamily_sccs is SCOP sccs identifying the superfamily for the domain
     superamily_description is SCOP dessription of the superfamily
     fold_description is the SCOP descriptino of the fold the domain is in
    """
    scopdominfo_dict = {}
    for scop_dom in scop.getDomains():
        sid = scop_dom.sid
        scop_superfamily = scop_dom.getAscendent('superfamily')
        scop_fold = scop_dom.getAscendent('fold')
        scop_class = scop_dom.getAscendent('class')
        scopdominfo_dict[sid] = (scop_superfamily.sccs,
                                 scop_superfamily.description,
                                 scop_fold.sccs,
                                 scop_fold.description)

    return scopdominfo_dict


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " cachefile\n")
    sys.exit(1)

    
def main():
    """
    main for scomdominfo.py

    Usage: scomdominfo.py  cachefile

    cachefile is the file to create the pickled domain info dictionary as
    WARNING: overwritten if it exists
    """
    if len(sys.argv) != 2:
        usage(os.path.basename(sys.argv[0]))

    pickle_filename = sys.argv[1]
        
    sys.stderr.write("Reading SCOP Data...")
    scop = Scop(dir_path=SCOP_DIR,version=SCOP_VERSION)
    sys.stderr.write("done\n")

    sys.stderr.write("Building domain info cache...")
    scopdominfo_dict = build_scopdominfo_dict(scop)
    sys.stderr.write("done. Got %d domain descriptions\n" %
                     len(scopdominfo_dict))

    sys.stderr.write("Writing cache to file %s...\n" % pickle_filename)
    fh = open(pickle_filename, "w")
    pickle.dump(scopdominfo_dict, fh)
    fh.close()
    sys.stderr.write("done\n")
    
    
            
if __name__ == "__main__":
    main()
    
