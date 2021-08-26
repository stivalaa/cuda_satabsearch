#!/usr/bin/env python
#
#
# rocrcops.py - Output score+label table for results on COPS dataset
#
# File:    rocrcops.py
# Author:  Alex Stivala
# Created: May 2010
#
# Evaluate structure search for COPS benchmark data set
# (Frank et al. 1999 "COPS Benchmark: interactive analysis of database 
# search methods" Bioinformatics 26(4):574-575) available from
# http://benchmark.services.came.sbg.ac.at/
#
# $Id: rocrcops.py 3635 2010-05-12 06:48:14Z alexs $
# 
#

"""
Write scores from structure search method and actual class labels
(0/1 for same/different fold as query).
Using COPS benchmark true positives as gold standard.

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


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# The true positives are stored in a text file, not embedded here

COPS_DIR = "/home/alexs/phd/qptabsearch/data/COPS/"
COPS_TP_FILE = COPS_DIR + "cops.truepositives"
COPS_QUERYLIST = COPS_DIR + "cops.querylist"
COPS_DBLIST = COPS_DIR + "cops.dblist"

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def parse_cops_tp_file(fname):
    """
    Parse the COPS true positives file (copied from COPS readme.txt file).
    Each (whitespace delmited) line starts with query id then has
    (exactly 6) true positives for that query.
    Comment lines starting with # are ignored.

    For conveniene for methods that mess with case of identifiers,
    everything is converted to lowercase here.
    
    Parameters:
       filename - name of COPS true positives text file to parse

    Return value:
       dict { queryid : tp_list } where queryid is query structure name
          and tp_list is list of the true positives structure names for
          that query.
    """
    tp_dict = {}
    for line in open(fname):
        if line[0] == '#':
            continue
        sline = line.split()
        if len(sline) < 6:
            sys.stderr.write('bad line in COPS tp file: %s\n' % line)
            continue
        tp_dict[sline[0].lower()] = [s.lower()  for s in sline[1:]]
    return tp_dict

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-mnv] [-z score] "
                     "  <outdir>\n")
    sys.stderr.write('  -m read multiquery file on stdin\n')
    sys.stderr.write('  -n negate scores (so that most -ve is best)\n')
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -z score : assign identifiers not present in the output a score of score\n')
    sys.exit(1)

    
def main():
    """
    main for rocrcops.py

    Usage: rocrcops.py [-mnv] [-z score]   <outdir>

    
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

    cops_truepositives = parse_cops_tp_file(COPS_TP_FILE)
    cops_dblist = [line.rstrip().lower() for line in open(COPS_DBLIST)]
    cops_querylist = [line.rstrip().lower() for line in open(COPS_QUERYLIST)]
    
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
                goldstd_ids = cops_truepositives[query_id.lower()]
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
                goldstd_ids = cops_truepositives[query_id]
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
            for domid in cops_dblist:
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
    if qcount < len(cops_querylist):
        sys.stderr.write("WARNING: only %d of %d queries have results\n" % (qcount, len(cops_querylist)))
        for qid in cops_querylist:
            if qid not in searchresult_dict.iterkeys():
                sys.stderr.write("query %s has no results\n" %qid)
                if bottom_score != None:
                    sys.stdout.write('# %s\n' % qid) # XXX
                    goldstd_ids = cops_truepositives[qid.lower()]
                    goldstd_dict[qid.lower()] = goldstd_ids
                    goldstd_pos_dict = dict([(domainid,True) for domainid in 
                                             goldstd_dict[qid.lower()]])
                    lowscore_domains = 0
                    for domid in cops_dblist:
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

