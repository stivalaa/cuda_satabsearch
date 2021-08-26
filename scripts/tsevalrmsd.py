#!/usr/bin/env python
###############################################################################
#
# tsevalrmsd.py - Evaluate tableaux search TPR/FPR against superposition
#
# File:    tsevalrmsd.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id: tsevalrmsd.py 2079 2009-03-03 07:43:11Z astivala $
# 
###############################################################################

"""
Evaluate true/false positives against superposition of
the matched substructures at each cutoff (rank/score).
Domains that are included below the cutoff and
with superposition RMSD below some (constant) threshold
at the cutoff are considered true positives, and false positives if they
are included below the cutoff score but RMSD is above the (constant) threshold.

Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt

from tsevalutils import compute_auc,parse_searchresult,eval_fn

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# RMSD threshold equal to or below which structures considered
# adequately superimposed (Angstroms)
RMSD_THRESHOLD = 3.5


#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class RMSDRecord:
    """
    Dummy class for containing data from line of superimposessemap.py
    output parsed by parse_superimposessemap() i.e.

    identifier, score, nsses nres, rmsd
    """
    pass


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def parse_superimposessemap(fh):
    """
    Parse the output of the superimposessemap.py script:
    format is one result per line, fields whitespace delimited:
    
    identifier score num_sses_matched num_aligned_points rmsd

    e.g.

    d1t10a_ -40.9999 8 16 16.93

    num_aligned_points is number of points used in the superposition,
    RMSD is the RMS deviation of those points (in Angstroms).

    Parameters:
       fh - open (read) filehandle to read superimposessemap.py output from
    Return value:
       list of RMSDRecord objects, each having the fields
          identifier - identifier of structure
          score - score of matching query to structure
          nsses - number of SSEs matched by maximlly similar subtableau finding
          nres - number of points used in superposition of query and structure
          rmsd - RMSD of superposition with the nres points
    """
    rr_list = []
    for line in fh:
        if line[0] == '#':
            continue # skip over comment lines
        rr = RMSDRecord()
        (rr.identifier, rr.score, rr.nsses, rr.nres, rr.rmsd) = line.split()
        rr.score = float(rr.score)
        rr.nsses = int(rr.nsses)
        rr.nres = int(rr.nres)
        rr.rmsd = float(rr.rmsd)
        rr_list.append(rr)
    return rr_list

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-v] "
                     " <tabsearch.outputfile>\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('superimposessemap.py output is read from stdin\n')
    sys.exit(1)

    
def main():
    """
    main for tsevalfn.py

    Usage: tsevalfn.py [-v]  <tabsearch.outputfile>
    
    -v turns on debug output to stderr

    superimposessemap.py output is read from stdin

    <tabsearch.outputfile> is the output from the tabsearchqpml.file,
    which is a text file where each line is identifier then whitespace
    then score, sorted by score from most negative to least negative e.g.

    d1xksa_ -35.99999999
    d3sila_ -35.99999999
    ....
    d2mhua_ -0.499999999

    ie this means the best hit as the top of the file, and worst at bottom.

    The superimposessemap.py output read on stdin
    format is one result per line, fields whitespace delimited:
    
    identifier score num_aligned_points rmsd

    e.g.

    d1t10a_ -40.9999 16 16.93

    num_aligned_points is number of points used in the superposition,
    RMSD is the RMS deviation of those points (in Angstroms).

    The table of positive and false negative rates is printed to stdout.
    """
    global verbose
    verbose = False
    

    try:
        opts,args = getopt.getopt(sys.argv[1:], "v?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    tabsearch_file = args[0]

    if verbose:
        sys.stderr.write('parsing supersimposessemap.py output...')

    rr_list = parse_superimposessemap(sys.stdin)
    if verbose:
        sys.stderr.write('parsed ' + str(len(rr_list)) + ' records\n')


    # get dict of identifiers where number of SSEs matched less than query size
    querysize = 8 #FIXME 8 for d1ubia_
    remove_dict = dict([(rr.identifier,True)
                        for rr in rr_list if rr.nsses < querysize])
        
    
    # select only those not less than query size and
    # with RMSD below threshold as gold standard positives
    goldstd_domains=[rr.identifier for rr in rr_list if not remove_dict.has_key(rr.identifier) and rr.rmsd <= RMSD_THRESHOLD]
    if verbose:
        sys.stderr.write('got ' + str(len(goldstd_domains)) +
                         ' structures below RMSD threshold ' +
                         str(RMSD_THRESHOLD) + '\n')
    sys.stderr.write(str([str(d) for d in goldstd_domains]))

    if verbose:
        sys.stderr.write('parsing search results...')
    tabsearch_fh = open(tabsearch_file)
    (searchresult,commentlist) = parse_searchresult(tabsearch_fh, negateflag=True)
    tabsearch_fh.close()
    if verbose:
        sys.stderr.write('parsed ' + str(len(searchresult)) + ' records\n')

    # set score to 0 for those where number of SSEs matched less than query size
    oldsearchresult = list(searchresult)
    searchresult = []
    for i in range(len(oldsearchresult)):
        # set score to 0 for identifiers matching fewer than query size SSEs
        if remove_dict.has_key(oldsearchresult[i][1]):
            searchresult.append((0.0, oldsearchresult[i][1]))
        else:
            searchresult.append(oldsearchresult[i])

    if verbose:
        sys.stderr.write('set score to 0 for ' + str(len(remove_dict)) +
                         ' structures matching fewer than ' +
                         str(querysize) + ' SSEs\n')
    searchresult.sort()
            

    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('# results from:\n')
    for line in commentlist:
        sys.stdout.write('#  ')
        sys.stdout.write(line) # identifying information about search run
    sys.stdout.write('\n')
    eval_fn(goldstd_domains, searchresult)
    
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
