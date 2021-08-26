#!/usr/bin/env python
###############################################################################
#
# multi2colout2single.py - Convert multiquery output to multiple single query
#
# File:    mkroctabs.py
# Author:  Alex Stivala
# Created: July 2010
#
# $Id: multi2colout2single.py 3866 2010-07-05 06:00:40Z alexs $
# 
###############################################################################

"""
Convert output file from a db search program (in the two-column (dbid
score) format), with a comment ('#' in first column) in this format:

# QUERYID = querid

or

# QUERY ID = querid

e.g.

# QUERY ID = d1ubia_

to delimit results for each separate query, to similar format but
where each query is in a separate file named by the query id.

See usage in docstring for main()

"""

import sys,os
from itertools import groupby
from tsevalutils import iter_searchresult


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " + progname + "  <outdir>\n")
    sys.stderr.write(" <outdir> is directory create output files in\n")
    sys.stderr.write(" WARNING: they will be overwritten if they exist\n")
    sys.stderr.write("input is multiquery file on stdin\n")
    sys.exit(1)

    
def main():
    """

    Usage: multi2colout2single.py <outdir>

         <outdir> is directory to crate output files. WARNING: they will
                  be overwritten if they exist.

         input is multiquery input on stdin

    """
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    outdir = sys.argv[1]

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
                          skip_self_query=False,
                          negateflag=False) ),
                               lambda t : t[0])
    
    for (query_id, result_iter) in query_group_iter:
        outfile = outdir + os.path.sep + query_id.lower() + os.path.extsep + "out"
        outfile_fh = open(outfile, 'w')
        for (domainid,score,dbid) in result_iter:
            outfile_fh.write(dbid + "    " + str(score) + "\n")
        outfile_fh.close()

            
if __name__ == "__main__":
    main()



