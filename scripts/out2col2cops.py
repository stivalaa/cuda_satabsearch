#!/usr/bin/env python
###############################################################################
#
# out2col2copy.py - Convert 2 column format (score,domainid) multiquery
#                    format to format 
#                    for COPS benchmark web server.
#
# See usage comment in docstring for main()
#
# File:    out2col2cops.py
# Author:  Alex Stivala
# Created: May 2010
#
# $Id: out2col2cops.py 3697 2010-05-18 08:05:56Z alexs $
# 
###############################################################################

"""
Converts the 2 column output from qptabsearch, SA tab search, etc.
to format for COPS benchmark web server

http://benchmark.services.came.sbg.ac.at/

The structure of the result file for the automatic analysis of the
benchmark results is very simple. For each query of the benchmark
there must be a line with the results of the predictor. The format of
these individual result lines is as follows:

<query> <hit #1> <hit #2> <hit #3> ...

Where query is the identifier of the query protein and the individual
hits are the identifiers of the respective results from the hitlist in
order of alignment quality. The different identifiers may be separated
either by one or more white space characters or one or more tabulator
signs.

Examples:

c2d4eD1 c3ek1A1 c1t90D1 c3efvA1 c2qe0D1 c1uxvA1 c2jg7H1 c1ad3B1 c1ez0A1 c2vroB1
c2oz5A_ c2q47B_ c2gx5C_ c1de5B2 c1d5rA1


Output is to stdout.


"""

import sys,os
import getopt
import random
from itertools import groupby

from tsevalutils import iter_searchresult



#-----------------------------------------------------------------------------
#
# constants
#
#-----------------------------------------------------------------------------

COPS_DIR = "/home/alexs/phd/qptabsearch/data/COPS/"
COPS_QUERYLIST = COPS_DIR + "cops.querylist"
COPS_DBLIST = COPS_DIR + "cops.dblist"

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [zn] < 2colout\n")
    sys.stderr.write('  -n negate scores\n')
    sys.exit(1)


def main():
    """
    main for out2col2cops.py

    Usage: out2col2cops.py [-n]

    -n negate all scores (so that most -ve is best)

    Input is on stdin,
    format is identifier and score

    Output for COPS web server as described in module info comment.
    """
    negateflag = False
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "nz:?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-n": # negate scores
            negateflag = True
        else:
            usage(os.path.basename(sys.argv[0]))


    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))

    cops_querylist = [line.rstrip().lower() for line in open(COPS_QUERYLIST)]
    cops_dblist = [line.rstrip() for line in open(COPS_DBLIST)]
    
    # get list of iterables each for same queryid.
    # iter_searchresult() is isterable of tuples (queryid, score, domainid)
    # groupby() requires that the iterable already has identical consecutive
    # queryids (first element of tuple) - iter_searchresult() should yield this
    # and does, when there is only one instance of the QUERY ID for each query,
    # but new version of cudaSaTabsearch does all queries in small structure
    # db, then in large structure db, so two QUERY ID for each query in the
    # output, so we sort by queryid first before groupby.
    done_queries_dict = {}
    query_group_iter = groupby(sorted(
        iter_searchresult(sys.stdin, multiquery=True,
                          negateflag=negateflag) ),
                               lambda t : t[0])

    for (queryid, result_iter) in query_group_iter:
        done_queries_dict[queryid.lower()] = True
        sys.stdout.write(queryid[:5].lower() + queryid[5:].upper() + ' ')
        sys.stdout.write(' '.join([dbid[:5].lower() + dbid[5:].upper()
                                   for (tqueryid,score,dbid) in
                                   sorted(result_iter,
                                          key  = lambda t : t[1],
                                          reverse = True)]))
        sys.stdout.write('\n')

    # output a query id and no results for any COPS query not in output
    # the web server actually insists on some results there,
    # so put a random db element (dodgy)
    for query in cops_querylist:
        if not done_queries_dict.has_key(query):
            sys.stderr.write('WARNING: no results for %s, putting a random one\n' % query)
            sys.stdout.write(query[:5].lower() + query[5:].upper() + '\t' + random.sample(cops_dblist, 1)[0] +  '\n')

if __name__ == "__main__":
    main()
