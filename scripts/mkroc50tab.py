#!/usr/bin/env python
###############################################################################
#
# mkroc50tab.py - Make table of ROC50 scores from multiquery output file
#
# File:    mkroc50tab.py
# Author:  Alex Stivala
# Created: December 2008
#
# $Id: mkroc50tab.py 3009 2009-12-08 03:01:48Z alexs $
# 
###############################################################################

"""
Build a table of ROC50 scores from a multiquery output file from a db
search program (in the two-column (dbid score) format), with a comment
('#' in first column) in this format:

# QUERYID = querid

e.g.

# QUERYID = d1ubia_

to delimit results for each separate query. The blastout2col.sh etc.
scripts create this format.

NB the results (within each query) are assumed to be sorted from
best score to worst score i.e. the ordering of the output from the 
program is used, rather than the numeric values of the scores as such.

This way the SCOP database only has to be loaded once, rather than
running evalroc50.py once for each query - loading the database is the
slowest part so this is much more efficient.

See usage in docstring for main()

"""

import sys,os
import getopt
from itertools import groupby
 
from Bio.SCOP import *

from pathdefs import SCOP_DIR,SCOP_VERSION
from tsevalutils import (iter_searchresult, get_goldstd_domains,
                         tsevalutils_set_verbose,is_true_positive)


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def compute_roc50(searchresult, goldstd_domains,
                 use_nonredundant, 
                 scop, astral):
    """
    Compute the ROC50 score by means of the Mann-Whitney U statistic
    (Wilcoxon rank-sum test) evaluated up to the first 50 false positives
    and then divided by product of number of false positives (ie 50)
    and the total number of true positives possible in the data.

    The ROC50 score is the area under the ROC curve plotted until 50 true
    negatives are found. This is computed by processing scores (from best
    to worst) marking them as true postivies (same fold/superfamily as query)
    or false positives (different fold/superfamily from query) until 50 false
    positives found.  For each of the first 50 false positives, the number
    of true positives with a better simliarity score is calculatd, and this
    is divided by the number of false postives (50). This is then divided
    by the total number of possible true postivies in the database
    (i.e. the number of mumbers in the SCOP fold/superfamily).
    
    ROC50 was defined by Gribskov & Robinson 1996 'Use of receiver operating
    characteristic (ROC) analysis to evaluate sequence matching'
    Computers Chem. 20(1):25-33
    and also used in e.g. Hulsen et al. 2006 'Testing statistical significance
    scores of sequence comparison methods with structural similarity'
    BMC Bioinformatics 7:444 from where this definition was taken.
    
    Note however we do not compute it directly this way, but by the
    rank-sum computation mentioned above.
    
    Parameters:
       searchresult -list of (score, domainid) tuples parsed from pogram to eval
                     This is assumed to be sorted from best to worst score
       goldstd_domains - list of gold standard domains to eval against
       use_nonredundant - Bool if True use ASTRAL 95% nr subset
       scop - previously built Bio.SCOP Scop instance
       astral - previously build Bio.SCOP Astral instance

    Return value:
       The ROC50 score for the searchresults evaluated against the 
       goldstd_domains.
    """
    NUM_FALSEPOS = 50 # could change this, but then not ROC50

    # convert goldstd list of domainids to dictionary keyed by domainid
    # for fast lookup as we iterate through search results.
    # The dictionary is { domainid : True } (we don't have a value,
    # just need to quickly test for presence of domainid in gold std
    # postive list).
    goldstd_pos_dict = dict([(sid, True) for
                       sid in [scopdom.sid for scopdom in goldstd_domains]])

    total_tp_possible = len(goldstd_domains) - 1 # don't count query itself
#    print 'bbb total tp possible',total_tp_possible

    labelled_ranks = [] # list of (sid, score, rank, is_pos)
    i = 0
    fp_count = 0
    while fp_count < NUM_FALSEPOS:
        try:
            (score, domainid) = searchresult[i]
        except IndexError:
            sys.stderr.write('WARNING: only %d search results, only got %d false positives\n' % (len(searchresult), fp_count))
            break
#        print 'zzz',score,domainid

        istruepos = is_true_positive(domainid, goldstd_pos_dict)
        if not istruepos:
            fp_count += 1
#            if fp_count == 1:
#                print 'eee first fp',domainid,'at ',i
        labelled_ranks.append((domainid, score, i, istruepos))
        i += 1

    if fp_count == 0:
        # if no false positives found, then set roc50 to 1.0 as it means
        # the results are 'perfect' up to the number of scores reported.
        sys.stderr.write('WARNING: found zero false positives, set roc50=1.0\n')
        return 1.0

    tp_count = len([sid for (sid,score,rank,istruepos) in labelled_ranks 
                    if istruepos])

#    print 'ddd fp_count',fp_count
#    print 'fff tp_count',tp_count
    assert(tp_count == len(labelled_ranks) - fp_count)
    
#    for tup in labelled_ranks:
#        print 'xxx',tup


#     ########################################################################
#     # naive 'direct' computation 
#     total = 0
#     for i in xrange(len(labelled_ranks)):
#         count = 0
#         if not labelled_ranks[i][3]:
#             for j in xrange(i):
#                 if labelled_ranks[j][3]:
#                     count += 1
#             total += count
#     print 'nnn total = ',total
#     print 'ooo ', float(total) / ( float(fp_count) * float(total_tp_possible) )
#     #########################################################################

    fp_rank_sum = sum([rank+1 for (sid,score,rank,istruepos)
                       in labelled_ranks if not istruepos])
    mann_whitney_U = fp_rank_sum - fp_count*(fp_count+1)/2
#    print 'uuuu U =',mann_whitney_U

    tp_rank_sum = len(labelled_ranks)*(len(labelled_ranks)+1)/2 - fp_rank_sum
    mann_whitney_U2 = tp_rank_sum - tp_count*(tp_count+1)/2
#    print 'vvvv U2 = ',mann_whitney_U2

    assert(mann_whitney_U + mann_whitney_U2 == fp_count * tp_count)

    auc50 = float(mann_whitney_U) / float(fp_count*total_tp_possible)


#     #########################################################################
#     # average precision calculation
#     total = 0
#     for i in xrange(len(labelled_ranks)):
#         higher_tp_count = 1
#         if labelled_ranks[i][3]:
#             for j in xrange(i):
#                 if labelled_ranks[j][3]:
#                     higher_tp_count += 1
#             pos_count = i+1
#             prec = float(higher_tp_count) / float(pos_count)
#             print 'qqqq i prec = ',i,prec
#             total += prec
#     avgprec = total / float(len(labelled_ranks)) # paper says this; seems wrong
# #    avgprec = total / float(tp_count)
#     print 'pppp avgprec = ',avgprec
#     #########################################################################


    return auc50



#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [a|u][-vep]\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -u evaluate at superfamily rather than fold level\n')
    sys.stderr.write('  -a evaluate at family rather than fold level\n')
    sys.stderr.write('  -e version: SCOP version (default %s)\n' %  SCOP_VERSION)
    sys.stderr.write('  -p percent: ASTRAL sequence percent nonredundant subset to use\n')
                     
    sys.exit(1)

    
def main():
    """
    main for mkroc50tab.py

    Usage: mkroc50tab.py [-vep] [a|u]

    
    -v turns on debug output to stderr

    -u use domains from the same superfamily as the query rather than from
       the same fold as the query as the gold standard.

    -a use domains from same family as the query rather tahn from the
       same fold as the query as the gold standard.

    -e version: SCOP version to use (default 1.73)

    -p percent: ASTRAL sequence nonredundant subset percentage to use
                (e.g. 95). If not specified, use all domains.

    input on stdin is the two column (identifier score) output of the
    program to evaluate for that query as described in module header
    docstring.

    output on stdout is the table of ROC50 values (two column (queryid roc50)) 
    for all queries concatenated together, with header for use in R
    read.table(header=TRUE).
    """
    global verbose
    verbose = False
    
    negateflag = False
    use_prosmos_dataset = False
    use_nonredundant = False
    use_superfamily = False
    use_family = False
    bottom_score = None
    logflag = False
    nrpercent = None
    scop_version = SCOP_VERSION

    try:
        opts,args = getopt.getopt(sys.argv[1:], "auve:p:?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
            tsevalutils_set_verbose(True) # tsevalutils.py verbose flag
        elif opt == "-u": # evaluate at superfamily not fold level
            use_superfamily = True
        elif opt == "-a": # evaluate at family not fold level
            use_family = True
        elif opt == "-p": # sequence identity percentage for nr subset
            nrpercent = int(arg)
            use_nonredundant = True
        elif opt == "-e": # SCOP version number
            scop_version = arg
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))

    if use_family and use_superfamily:
        sys.stderr.write('-a (family) and -u (superfamily) are mutually exclusive\n')
        usage(sys.argv[0])

    level = 'fold'
    if use_family:
        level = 'family'
    elif use_superfamily:
        level = 'superfamily'

    # read SCOP and ASTRAL data
    if verbose:
        sys.stderr.write('Reading SCOP data (version %s)...\n' % scop_version)
    scop = Scop(dir_path=SCOP_DIR,version=scop_version)
    astral = Astral(dir_path=SCOP_DIR,version=scop_version,scop=scop)


    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('queryid  roc50\n')
    sys.stdout.write('#-------------\n')

    # get list of iterables each for same queryid.
    # iter_searchresult() is isterable of tuples (queryid, score, domainid)
    # groupby() requires that the iterable already has identical consecutive
    # queryids (first element of tuple) - iter_searchresult() should yield this
    query_group_iter = groupby(iter_searchresult(sys.stdin, multiquery=True,
                                                 skip_self_query=True),
                               lambda t : t[0])
    
    total = 0
    num = 0
    error_count = 0
    for (queryid, result_iter) in query_group_iter:
        if verbose:
            sys.stderr.write('processing query id %s\n' % queryid)

        try:
            goldstd_domains = get_goldstd_domains(queryid, 
                                                  level,
                                                  use_nonredundant,
                                                  scop, astral,
                                                  nrpercent)
        except AttributeError:
            # AttributeError: 'NoneType' object has no attribute 'getAscendent'
            # when we get a None from looking up query domain identifier
            sys.stderr.write('ERROR: query sid %s not found in SCOP, skipping\n' % queryid)
            error_count += 1
            continue
        
        if len(goldstd_domains) < 2: # At least 1 for query itself
            sys.stderr.write('ERROR: query sid %s has no domains but self in same %s in ASTRAL SCOP subset, skipping\n' % (queryid,level))
            error_count += 1
            continue


        roc50 = compute_roc50([(score,dbid) 
                               for (queryid,score,dbid) in result_iter], 
                              goldstd_domains,
                              use_nonredundant, scop, astral)
        total += roc50
        num += 1
        sys.stdout.write('%s    %f\n'  % (queryid, roc50))

    avg = total / num
    sys.stdout.write('# AVERAGE = %f\n' % (avg))
            
if __name__ == "__main__":
    main()
