#!/usr/bin/env python
###############################################################################
#
# soln2ssemap.py - Convert solution x vector from tsrchd/tsrchn to SSE mapping
#
# File:    solns2ssemap.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id: soln2ssemap.py 4135 2010-09-07 06:51:43Z alexs $
# 
###############################################################################

"""
Parse the output of the FORTRAN-77 tsrchn/tsrchd programs with the LSOLN
logical option set so that solution x vector is output as well as score
(objective function value at the solution), and convert to a correspondance
between SSEs in each structure (sequentially numbered from 1 along chains).

Note we need to know the SSEs (or at least the number of them) so this
script actually gets the name of the tableaux database from the
tsrchd/tsrchnd output file it parses, and reads that information from the
database.

This all is rather cumbersome and inefficient, but as we actually only
want the matching for comparatively few 'hits' (the best scoring ones)
it would be wasteful to have to do it for every match in the FORTRAN code
itself, this way the searching stays fast and we just have this slow
and cumbersome step to get the matching (and eg generate a PyMOL script
to show the correspondences with colours in PyMOL from it with ssemap2pml.py).

Requires the Numeric library:
   http://sourceforge.net/projects/numpy

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
from time import strftime,localtime

import numpy.oldnumeric as Numeric

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

EPS = 1e-1 # epsilon for closesness to 0/1

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class SearchSoln:
    """
    SearchSoln is just a dummy class for containing the search results
    with solution vectors, returned by parse_searchsoln()
    """
    pass

class QuerySoln:
    """
    QuerySoln is a dummy class for containign result from individual query,
    in SearchSoln.
    """
    pass

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def parse_searchsoln(fh, domid=None):
    """
    Parse the output of tsrchn/tsrchd with LSOLN set to .true.
    so each structure has first structure name and
    scores (raw,nor2m,zscore,pvalue), then
    x vector one element per line e.g. 
    
    # TSRCHD LTYPE = T LORDER = T LSOLN = T
    # QUERY ID = D1KI9A_ 
    # DBFILE = /local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii                     
    # Mon Aug  4 12:34:07 2008

    d1xksa_  41 1.51852 0.903567  0.161557
      0.0000
      1.0000
      0.0000
      0.0000
      0.0000
      ...
    d3sila_  21 0.823529 -1.27278 0.943443
      0.0000
      ...
      etc.

    Parameters:
        fh - open (read) filehandle to parse from
        domid - (default None) if not None, get only this identifier,
                not whole file.
        
    Return value:
        search_soln - dummy class SearchSoln containing:
             queryid - query identifier parsed from comments
             dbfile - filename of database file parsed from comments
             query_soln_list - list of dummy class QuerySoln containing:
                 domid - id of domain in db
                 score - score of queryid against domid
                 soln - list of values in solution vector
             comment_lines - list of comment lines read
    """
    search_soln = SearchSoln()
    query_soln_list = []
    query_soln = None
    search_soln.comment_lines = []
    for line in fh:
        if line[0] == '#':
            sline = line[1:].split('=')
            if sline[0].lstrip().rstrip() == 'QUERY ID':
                search_soln.queryid = sline[1].lstrip().rstrip().lower()
            elif sline[0].lstrip().rstrip() == 'DBFILE':
                search_soln.dbfile = sline[1].lstrip().rstrip()
            search_soln.comment_lines.append(line)
            continue
        if len(line.split()) > 1:  # new identifier
            if query_soln:
                query_soln_list.append(query_soln)   # finsished with prev one
                if domid:
                    query_soln = None # appended here, don't do again
                    break # if only getting one, we have finished
            splitline = line.split()
            if len(splitline) != 5:
                sys.stderr.write('bad line: ' + line + '\n')
                continue
            domainid = splitline[0]
            if (domid and domid != domainid):
                continue # only interested in domid, skip others
            score_str = splitline[1] # raw score
#            score_str = splitline[4] # now use pvalue
            if score_str.lower() == 'nan' or score_str == '********':
                # if we get any NaN values then then sort() gets completely
                # screwed up and the evaluation is all wrong, so skip them.
                sys.stderr.write('skipping NaN: ' + line + '\n')
                query_soln = None
                continue
            score = float(score_str)
            query_soln = QuerySoln()
            query_soln.domid = domainid
            query_soln.score = score
            query_soln.soln = []
        else:
            # line should be an element of the solution vector
            try:
                xval = float(line)
            except ValueError:
                sys.stderr.write('bad float value: ' + line + '\n')
                # skip whole query of a value is bad
                if query_soln:
                    sys.stderr.write('skipping ' + query_soln.domid + '\n')
                    query_soln = None
                continue
            if query_soln:
                query_soln.soln.append(xval)
    if query_soln:
        query_soln_list.append(query_soln)
        
    search_soln.query_soln_list = query_soln_list
    return search_soln


def parse_tableauxdb_sizes(fh):
    """
    Parse the dimensinos of tableaux from the ascii tableauxdb

    Parameters:
        fh - open (read) filehandle for ascii tableauxdb file (numeric or
             discrete)
    Return value:
        dictionary { domid : dim } mapping domain identifiers to tableau
        dimension (numer of sses)
    """
    dimdict = {}
    line = fh.readline()
    while line != "":
        (domid, sdim) = line.split()
        dim = int(sdim)
        dimdict[domid.lstrip().rstrip()] = dim
        i = 0
        while i < dim * 2:
            line = fh.readline()  # read through the tableau and dist matrix
            i += 1
        line = fh.readline() # read the blank line between entries
        line = fh.readline()
    return dimdict


def soln_list_to_matrix(soln_vec, n1, n2):
    """
    Convert the solution vector in the form of list of floats to matching
    matrix where (i,j) is 1.0 for SSE i matching SSE j.
    This requires dimensions of the two tableaux.

    Parameters:
       soln_vec - solution vector as list of floats, length n1*n2+n1+n2
       n1 - dimension of tableau (number of SSEs) for 1st structure
       n2 - dimension of tableau (number of SSEs) for 2nd structure

    Return value:
       Numeric.array matrix dimension n1 x n2 of matching matrix.
    """
    x = Numeric.array(soln_vec[:n1*n2])
    matchmat = Numeric.reshape(x, (n1, n2))
    return matchmat

def matrix_to_tuplelist(matchmat):
    """
    Convert the matrching matrix in the form of Numeric.array
    matrix where (i,j) is 1.0 for SSE i matching SSE j
    to a list of tuples (k,l) where k is sequential (from 1) SSE number
    in first structure and l is sequential (from 2) SSE number in second
    structure.

    Parameters:
        matchmat - Numeric.array matching matrix (n1 x n2)
    Return value:
        list of tuples as described above.
    """
    matchlist = []
    i = 0
    while i < Numeric.shape(matchmat)[0]:
        j = Numeric.argmax(matchmat[i])
        if matchmat[i,j] >= 1.0 - EPS:
            matchlist.append((i+1,j+1))
        i += 1
    return matchlist


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-d domainid] [-q querysize]\n")
    sys.stderr.write("  -d domainid: only output result for hit against domainid\n")
    sys.stderr.write("  -q querysize: specify size of query tableau\n")
    sys.exit(1)

    
def main():
    """
    main for soln2ssemap.py

    Usage: sol2ssemap.py [ -d domainid ] [ -q querysize ]

    -d domainid: only output result for query against that domain (SCOP)  id,
                 not all.

    -q querysize: specify size of query tableau, for when the query structure
                  is not itself in the database. NB if this does not match
                  actual size of query tableau, results undefined
                  (probably get ValueError on Numeric.reshape()).
                 
    Input is on stdin, the output of tsrchn/tsrchd with LSOLN set to .true.
    so each structure has first structure name and
    scores (raw, norm2, zscore, pvalue), then
    x vector one element per line e.g. 
    
    # TSRCHD LTYPE = T LORDER = T LSOLN = T
    # QUERY ID = D1KI9A_ 
    # DBFILE = /local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii                     
    # Mon Aug  4 12:34:07 2008
    d1xksa_  41 1.51852 0.903567  0.161557    
      0.0000
      1.0000
      0.0000
      0.0000
      0.0000
      ...
    d3sila_  21 0.823529 -1.27278 0.943443    
      0.0000
      ...
      etc.

    Note the header 'comment' information is actually required, to get
    the query id since we need to get the SSEs for that structure
    (as well as the structures matched against it), we read this
    from the tableaux database, the location of which is also in
    the information at the top of the file we parse.

    Output format is identifier and score (as per input), then
    for each matching a line containing
    i and j separated by a space,
    one per line (with blank line before next id) e.g.:

    # TSRCHD LTYPE = T LORDER = T LSOLN = T
    # QUERY ID = D1KI9A_ 
    # DBFILE = /local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii
    # Mon Aug  4 12:34:07 2008
    d1wiua_     -23.0000
    1 1
    3 2
    8 4
    9 5
    11 6
    14 9

    Note we copy the QUERY ID and DBFILE and other information to output for use
    in later processing.

    """
    global verbose
    verbose = False

    dbdomid = None
    qsize = None
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "d:q:")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-d": # domain id specified, only get this one
            dbdomid = arg
        elif opt == "-q": # specify query tableau szie
            qsize = int(arg)
        else:
            usage(os.path.basename(sys.argv[0]))
        
    if len(args) != 0:
        usage(os.path.basename(sys.argv[0]))

    search_soln = parse_searchsoln(sys.stdin, dbdomid)
    dim_dict = parse_tableauxdb_sizes(open(search_soln.dbfile))
    n1 = None
    try:
        n1 = dim_dict[search_soln.queryid]
    except KeyError:
        if not qsize:
            sys.stderr.write("Query structure " + search_soln.queryid +
                             " not found in database. Use -q to specify" +
                             " query tableau size.\n")
            sys.exit(1)
    if qsize and n1 and qsize != n1:
        sys.stderr.write("WARNING: -q " + str(qsize) + " does not match "
                         "size of query " + search_soln.queryid +
                         " in database (" + str(n1) + "). Using " +
                         " user-specified value " + str(qsize) + "\n")
        n1 = qsize
    elif qsize and not n1:
        n1 = qsize
    sys.stdout.write('# ' + os.path.basename(sys.argv[0]) +
                     ' processed:\n')
    sys.stdout.write('#\n')
    for cline in search_soln.comment_lines:
        sys.stdout.write(cline)
    sys.stdout.write('#\n')
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    sys.stdout.write('# on ' + timestamp + '\n')
    sys.stdout.write('#\n')    
    for query_soln in search_soln.query_soln_list:
        if ((not dbdomid) or (query_soln.domid == dbdomid)):
            n2 = dim_dict[query_soln.domid]
            matchmat = soln_list_to_matrix(query_soln.soln, n1, n2)
            matchlist = matrix_to_tuplelist(matchmat)
            sys.stdout.write('%s %g\n' % (query_soln.domid,query_soln.score))
            for (i,j) in matchlist:
                sys.stdout.write(str(i) + ' ' + str(j) + '\n')
            sys.stdout.write('\n')
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
