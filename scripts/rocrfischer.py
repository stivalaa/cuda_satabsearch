#!/usr/bin/env python
#
#
# rocrfischer.py - Output score+label table for results on Fischer dataset
#
# File:    rocrfischer.py
# Author:  Alex Stivala
# Created: September 2008
#
# Evaluate structure search for all against all matching
# for the Fischer data set (Fischer et al 1996 Pac. Symp. Biocomput. 300-318))
# as per Pelta et all 2008 BMC Bioinformatics 9:161
# Output a table of scores and true class labels (binary: in or not in same
# class/fold) , for use with R CRAN package ROCR
#
# $Id: rocrfischer.py 3603 2010-05-04 04:47:51Z alexs $
# 
#

"""
Write scores from structure search method and actual class labels
(0/1 for same/different fold as query).
Using Fischer Table II as the gold standard at fold or class level

Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os,glob
import getopt
from itertools import groupby

from tsevalutils import parse_searchresult,iter_searchresult
from fischer_tables import FISCHER_ID_FOLD_DICT,FISCHER_FOLD_IDLIST_DICT,FISCHER_ID_CLASS_DICT,FISCHER_CLASS_IDLIST_DICT

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-cmnv] [-z score] "
                     "  <outdir>\n")
    sys.stderr.write('  -c class level not fold level evaluation\n')
    sys.stderr.write('  -m read multiquery file on stdin\n')
    sys.stderr.write('  -n negate scores (so that most -ve is best)\n')
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -z score : assign identifiers not present in the output a score of score\n')
    sys.exit(1)

    
def main():
    """
    main for rocfischer.py

    Usage: rocfischer.py [-cmnv] [-z score]   <outdir>

    
    -c evaluate at class level rather than default fold level
    -v turns on debug output to stderr
    -n negate all scores (so that most -ve is best)
    -m read multiquery file (all query results in one file) from stdin
       instead of separate .out file for each query
    -z score : any identifiers that are not present in the output but
               are in the gold standard data are given the specified score.
               This is for programs that do not assign a score to all
               domains, but only those above some threshold or the top n,
               or just cannot assign a score to some domains for some reason.
               This would generally be specified as some score lower 
               than all other scores.


    <outdir> is the directory containing output files as generated
    by tsrchd_sparse (qptabmatch_allall.py) or msvns4maxcmo_allall.py 

    The table of scores and labels is printed to stdout.
    """
    global verbose
    verbose = False
    use_class = False
    negateflag = False
    multiquery = False
    bottom_score = None
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "vz:ncm?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True # this module only
        elif opt == "-n": # negate scores
            negateflag = True
        elif opt == '-c': # class not fold level evaluation
            use_class = True
        elif opt == '-m': # read multiquery file
            multiquery = True
        elif opt == "-z": # score to give to domains that have no score
            bottom_score = float(arg)
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    outdir = args[0]

    # build dict of search results, one for each query
    # and corresponding dict of gold standard results

    searchresult_dict = {}
    goldstd_dict = {}

    if multiquery:
        # get list of iterables each for same queryid.
        # iter_searchresult() is isterable of tuples (queryid, score, domainid)
        # groupby() requires that the iterable already has identical consecutive
        # queryids (first element of tuple) - iter_searchresult() should yield this
        # and does, when there is only one instance of the QUERY ID for each query,
        # but new version of cudaSaTabsearch does all queries in small structure
        # db, then in large structure db, so two QUERY ID for each query in the
        # output, so we sort by queryid first before groupby.
        query_group_iter = groupby(sorted(
            iter_searchresult(sys.stdin, multiquery=True,
                              skip_self_query=True,
                              negateflag=negateflag) ),
                                   lambda t : t[0])

        for (query_id, result_iter) in query_group_iter:
            try:
                if use_class:
                    goldstd_ids = FISCHER_CLASS_IDLIST_DICT[FISCHER_ID_CLASS_DICT[query_id.lower()]]
                else:
                    goldstd_ids = FISCHER_FOLD_IDLIST_DICT[FISCHER_ID_FOLD_DICT[query_id.lower()]]
            except KeyError:
                if verbose:
                    sys.stderr.write('skipped ' + query_id + '\n')
                continue
            searchresult_dict[query_id.lower()] = [(score, dbid.lower()) for (domainid,score,dbid)  in result_iter]
            goldstd_dict[query_id.lower()] = goldstd_ids
    else:
        # each query is in a separate .out file
        for result_file in glob.glob(os.path.join(outdir, '*.out')):
            query_id = os.path.splitext(os.path.basename(result_file))[0].lower()
            result_fh = open(result_file)
            (searchresult,commentlist) = parse_searchresult(result_fh, negateflag)
            result_fh.close()
            try:
                if use_class:
                    goldstd_ids = FISCHER_CLASS_IDLIST_DICT[FISCHER_ID_CLASS_DICT[query_id]]
                else:
                    goldstd_ids = FISCHER_FOLD_IDLIST_DICT[FISCHER_ID_FOLD_DICT[query_id]]
            except KeyError:
                if verbose:
                    sys.stderr.write('skipped ' + query_id + '\n')
                continue
            searchresult_dict[query_id] = searchresult
            goldstd_dict[query_id] = goldstd_ids


    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('score    label\n')
    sys.stdout.write('#-------------\n')
    qcount = 0
    for query_id in searchresult_dict.iterkeys():
        sys.stdout.write('# %s\n' % query_id) # XXX
        slcount=0
        qcount += 1
        searchresult = searchresult_dict[query_id]
        goldstd_pos_dict = dict([(domainid,True) for domainid in 
                                 goldstd_dict[query_id.lower()]])
        if len(searchresult) == 0:
            sys.stderr.write("warning: no results for query %s\n" % query_id)
        for (score, domainid) in searchresult:
            # skip self-query
            if domainid.lower() == query_id.lower():
                continue
            if goldstd_pos_dict.has_key(domainid.lower()):
                label = 1
            else:
                label = 0
            sys.stdout.write('%20.8f  %d\n' % (score,label))
            slcount += 1

        if bottom_score != None:
            lowscore_domains = 0
            for domid in FISCHER_ID_FOLD_DICT.keys():
                if domid.lower() == query_id.lower():
                    continue #skip self-query
                if domid.lower() not in [d.lower() for (s,d) in searchresult]:
                    lowscore_domains += 1
                    if goldstd_pos_dict.has_key(domid.lower()):
                        label = 1
                    else:
                        label = 0
                    sys.stdout.write('%20.8f  %d\n' % (bottom_score,label))
                    slcount += 1
            if verbose and lowscore_domains > 0:
                sys.stderr.write("(queryid %s): set score to %f for %d domains\n" % (query_id, bottom_score, lowscore_domains))
        if verbose:
            sys.stderr.write('wrote %d (score,label) pairs for query %s\n' % (slcount, query_id))
                    
    if verbose:
        sys.stderr.write("processed %d queries\n" % qcount)

    # some methods (actually only VAST as far as I've found) actually
    # give NO results for some queries, which is a major hassle specially
    # for StAR which requires all methods to have all results.
    # so we'll just give the bottom score to all those matchings
    if qcount < len(list(FISCHER_ID_FOLD_DICT.iterkeys())):
        sys.stderr.write("WARNING: only %d of %d queries have results\n" % (qcount, len(list(FISCHER_ID_FOLD_DICT.iterkeys()))))
        for qid in FISCHER_ID_FOLD_DICT.iterkeys():
            if qid not in searchresult_dict.iterkeys():
                sys.stderr.write("query %s has no results\n" %qid)
                if bottom_score != None:
                    sys.stdout.write('# %s\n' % qid) # XXX
                    if use_class:
                        goldstd_ids = FISCHER_CLASS_IDLIST_DICT[FISCHER_ID_CLASS_DICT[qid.lower()]]
                    else:
                        goldstd_ids = FISCHER_FOLD_IDLIST_DICT[FISCHER_ID_FOLD_DICT[qid.lower()]]
                    goldstd_dict[qid.lower()] = goldstd_ids
                    goldstd_pos_dict = dict([(domainid,True) for domainid in 
                                             goldstd_dict[qid.lower()]])
                    lowscore_domains = 0
                    for domid in FISCHER_ID_FOLD_DICT.keys():
                            if domid.lower() == qid.lower():
                                continue #skip self-query
                            lowscore_domains += 1
                            if goldstd_pos_dict.has_key(domid.lower()):
                                label = 1
                            else:
                                label = 0
                            sys.stdout.write('%20.8f  %d\n' % (bottom_score,label))
                    if verbose and lowscore_domains > 0:
                        sys.stderr.write("(queryid %s): set score to %f for %d domains\n" % (qid, bottom_score, lowscore_domains))
                    
        
    
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()

