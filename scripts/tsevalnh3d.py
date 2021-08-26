#!/usr/bin/env python
###############################################################################
#
# tsevalnh3d.py - Evaluate tableaux search against Nh3D data set
#
# File:    tsevalnh3d.py
# Author:  Alex Stivala
# Created: June 2008
#
# Evaluate QP tableau search for all against all matching
# for the Nh3D data set at CATH architecture level
# as per Pelta et all 2008 BMC Bioinformatics 9:161
#
# $Id: tsevalnh3d.py 1956 2008-10-07 00:28:16Z astivala $
# 
###############################################################################

# OBSOLETE - use rocrnh3d.py and rocauc.r now

"""
Evaluate false negatives using CATH architecture :
domains that are not included above the cuttoff but have same architecture id
according to CATH are false negatives.
Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()
"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt

from cathmap import CATHMAP
from tsevalutils import compute_auc,parse_searchresult,eval_fn

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-rv] "
                     " <query_id> <tabsearch.outputfile>\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -r higher scores are better (for MSVNS4MaxCMO)\n')
    sys.exit(1)

    
def main():
    """
    main for tsevalnh3d.py

    Usage: tsevalnh3d.py [-rv]  <query_id> <tabsearch.outputfile>

    
    -r higher scores are better (rather than default of more negative
       scores better). MSVNS4MaxCMO requires this (QP tableau search
       gives lower (more negative) scores for better matches)
    -v turns on debug output to stderr

    <query_id> is the CATH id (e.g. 1.10.1290) of the query structure
    
    <tabsearch.outputfile> is the output from the tabsearchqpml.file,
    which is a text file where each line is identifier then whitespace
    then score, sorted by score from most negative to least negative e.g.

    1101290 -35.99999999

    ie this means the best hit as the top of the file, and worst at bottom.
    May be specified as - for stdin.

    The table of positive and false negative rates is printed to stdout.
    """
    global verbose
    verbose = False
    reverseflag = False

    try:
        opts,args = getopt.getopt(sys.argv[1:], "rv?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        elif opt == "-r": # reversed: higher scores better
            reverseflag = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    query_id = args[0]
    tabsearch_file = args[1]

    # get the gold standard as the list of 'compreseed' CATH ids
    # that have same architecture as query
    query_id_split = query_id.split('.')
    query_class = query_id_split[0]
    query_arch = query_id_split[1]
    goldstd_ids = []
    for (compressed_id, cathid) in CATHMAP.iteritems():
        cathid_split = cathid.split('.')
        cathid_class = cathid_split[0]
        cathid_arch = cathid_split[1]
        if cathid_class == query_class and cathid_arch == query_arch:
            goldstd_ids.append(compressed_id)

    if verbose:
        sys.stderr.write('parsing search results...\n')
    if tabsearch_file == '-':
        tabsearch_fh = sys.stdin
    else:
        tabsearch_fh = open(tabsearch_file)
    (searchresult,commentlist) = parse_searchresult(tabsearch_fh, reverseflag)
    if tabsearch_file != '-':
        tabsearch_fh.close()

    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('# results from:\n')
    for line in commentlist:
        sys.stdout.write('#  ')
        sys.stdout.write(line) # identifying information about search run
    sys.stdout.write('\n')
    eval_fn(goldstd_ids, searchresult, reverseflag)
    
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
