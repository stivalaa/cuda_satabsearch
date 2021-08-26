#!/usr/bin/env python
###############################################################################
#
# tszscore.py - compute Z scores for QP tableau search output
#
# File:    tszscore.py
# Author:  Alex Stivala
# Created: March 2009
#
# See usage in docstring for main()
#
# Imports from tsevalutils.py, so the directory containing it (ptgraph/)
# must be in the PYTHONPATH
# Also requires numpy.
#
# $Id: tszscore.py 2166 2009-03-30 03:59:01Z astivala $
# 
###############################################################################

import sys,os
import getopt

from tsevalutils import parse_searchresult
from numpy import mean,std

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def usage(progname):
    """ print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " [-no]\n")
    sys.stderr.write('  -n negate scores\n')
    sys.stderr.write('  -o take log10 of scores\n')
    sys.exit(1)

    
#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def main():
    """
    main for tszscore.py


    Given list of scores from QP tableau search (tsrchd_sparse etc.) output,
    compute Z-scores for hits to the query, and output sorted by
    descending Z-score.
    
    Usage:
       tszscore.py   < tsrchd_output

       -n negate scores
       -o take log10 of scores
       
       Input is from stdin, tab-delimited in the format
       
           pdbid score

       Output is to stdout, in the same format as input i.e.

           pdbid Z-score
    """
    negateflag = False
    logflag = False

    try:
        opts,args = getopt.getopt(sys.argv[1:], "no?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-n": # negate scores
            negateflag = True
        elif opt == "-o": # take log10 of scores
            logflag = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))


    (searchresult,commentlist) = parse_searchresult(sys.stdin,
                                                    negateflag,
                                                    logflag)

    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('# Z-scores computed from:\n')
    for line in commentlist:
        sys.stdout.write('#  ')
        sys.stdout.write(line) # identifying information about search run
    sys.stdout.write('\n')

    scores = [s[0] for s in searchresult]
    mu = mean(scores)
    sigma = std(scores)
    # Z-score for s is (s - mu) / sigma
    zscores = [ ( (s[0] - mu) / sigma, s[1]) for s in searchresult ]

    zscores.sort(reverse=True)
    for zs in zscores:
        sys.stdout.write("%s %5.3f\n" % (zs[1], zs[0]))
    
            
if __name__ == "__main__":
    main()
