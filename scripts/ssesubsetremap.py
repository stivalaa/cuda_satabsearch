#!/usr/bin/env python
###############################################################################
#
# ssesubsetremap.py - Remap SSE mapping output of solns2ssemap.py for the
#                     case that one of the tableaux was a subset of the SSEs
#                     in the structure.
#
# See usage comment in docstring for main()
#
# File:    ssesubsetremap.py
# Author:  Alex Stivala
# Created: November 2008
#
# $Id: ssesubsetremap.py 3953 2010-07-21 04:34:49Z alexs $
# 
###############################################################################

"""
In the case that one of the tableaux is not a whole structure but a
substructure (-s option on pytableaucreate.py (-e option on 
qptabmatchstructs.sh)) we cannot just use the mapping produced
by soln2ssemap.py (q.v.) as it just numbers SSEs 1,2,3,... along
the tableau. If one of the tableaux actually reprsents for example
SSEs 4,5,6,7,8 (or 2,6,9,10,11,12,18, need not be sequential) then
we need to map the sequential 1-based numbering of the tableau columns
from soln2ssemap back to the actual SSE sequential nubmers in the structure
according to the mapping provided by the SSE subset list used on 
pytableaucreate.py (and provided to this script).

This all is rather cumbersome and inefficient, but as we actually only
want the matching for comparatively few 'hits' (the best scoring ones)
it would be wasteful to have to do it for every match in the FORTRAN code
itself, this way the searching stays fast and we just have this slow
and cumbersome step to get the matching (and eg generate a PyMOL script
to show the correspondences with colours in PyMOL from it with ssemap2pml.py).

"""

import sys,os
import getopt
from time import strftime,localtime

from parsessemap import parse_ssemap,SearchMap,QuerySSEMap

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
    
    sys.stderr.write("Usage: " +progname + " [-5] [-d domainid] sse_num_list\n")
    sys.exit(1)

    
def main():
    """
    main for ssesubsetremap.py

    Usage: sol2ssemap.py [-5] [-d domainid] sse_num_list

    -d domainid: only output result for query against that domain (SCOP)  id,
                 not all.

    -5: use new 5-column (name rawscore norm2score zscore pvalue) format not
        old 2-column (name rawscore) format

    sse_num_list is the comma-separted list of sequential SSE numbers
    represented by the tableau of the first ('query') structure, 
    as supplied to pytableaucreate.py -s option and 
    qptabmatscructs.py -e option
    Note that this must be the same list used to 
    to produce the mapping, i.e. that the tableaux query
    was built with, otherwise it won't mean anything.

    Input is on stdin, the output of soln2ssemap.py,
    format is identifier and score, then
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

    Output is the same format, but the SSE numbers in the first column
    have been remapped according to the sse_num_list; i.e. the numbers
    in that column of the input as used as indices into the sse_num_list
    and the resulting numbers from that list are the correspdoning output.
    """
    global verbose
    verbose = False

    dbdomid = None
    use5col = False
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "5d:")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-d": # domain id specified, only get this one
            dbdomid = arg
        elif opt == "-5": # use 5 column score format
            use5col = True
        else:
            usage(os.path.basename(sys.argv[0]))
        
    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    sse_id_list_str = args[0].split(',')
    sse_id_list = []
    sse_id_uniq_dict = {} # { id : True } just for checking all unique
    for sse_id_str in sse_id_list_str:
        if sse_id_str.isdigit():
            if sse_id_uniq_dict.has_key(int(sse_id_str)):
                sys.stderr.write("duplicate SSE sequential number "  +
                                 sse_id_str + "\n")
                usage(sys.argv[0])
            sse_id_uniq_dict[int(sse_id_str)] = True
            sse_id_list.append(int(sse_id_str))
        else:
            sys.stderr.write("not a valid SSE sequential number '" +
                             sse_id_str + "'\n")
            usage(sys.argv[0])
    sse_id_list.sort() # ensure SSEs are in order

    search_maps = parse_ssemap(sys.stdin, use5col)

    sys.stdout.write('# generated by ' + os.path.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]) +'\n')
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    sys.stdout.write('# on ' + timestamp + '\n')
    sys.stdout.write('# from:\n')
    for cline in search_maps.comment_lines:
        sys.stdout.write(cline)
    sys.stdout.write('#\n')    
    for query_ssemap in search_maps.query_ssemap_list:
        if ((not dbdomid) or (query_ssemap.domid == dbdomid)):
            sys.stdout.write('%s %g\n' % (query_ssemap.domid,query_ssemap.score))
            for (i,j) in query_ssemap.sse_map:
                iprime = sse_id_list[i-1]
                sys.stdout.write(str(iprime) + ' ' + str(j) + '\n')



            
if __name__ == "__main__":
    main()
