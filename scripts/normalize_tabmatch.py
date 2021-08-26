#!/usr/bin/env python
###############################################################################
#
# normalize_tabmatch.py - normalize tableau match output scores by protein size
#
# File:    normalize_tabmatch.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Given list of scores from QP tableau search (tsrchd_sparse etc.) output,
# normalize the scores by proteni size.
#
# Usage:
#     normalize_tabmatch.py normtype <queryid | -m> <db_file | db_directory> 
#
# normtype is the normalization function: 1,2, or 3
# queryid is the identifier for the query structure, all the scores are
# assumed to be for matching this to the identifier on each line of the input.
# If -m is specified instead of a queryid, it is a multiquery input file with
# multiple queries, each with a '# QUERY ID = queyrid' comment line to
# delimit each query: the output file is then also a multiquery file.
# Input is from stdin, tab-delimited in the format
# 
#   pdbid score
#
# Also reads .tableuxdistmatrix files from db_directory, as used by
# qptabmatch_allparirs.py, tsrchd_sparse etc., to get 
# orders of tableaux (number of SSEs)
# required for normalization.
# If db_file is specified instead of db_directory, just reads tableaux/
# distmatrix db from that file.
#
# Output is to stdout, in the same format as input i.e.
#
#  pdbid normscore
#
# where normscore is the normalized score.
#
# $Id: normalize_tabmatch.py 3345 2010-02-16 03:42:05Z alexs $
# 
###############################################################################

import sys,os,glob
from itertools import groupby
from norms import norm1,norm2,norm3
from tsevalutils import iter_searchresult

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

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
        dimdict[domid.lstrip().rstrip().upper()] = dim
        i = 0
        while i < dim * 2:
            line = fh.readline()  # read through the tableau and dist matrix
            i += 1
        line = fh.readline() # read the blank line between entries
        line = fh.readline()
    return dimdict


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <normtype> <queryid|-m> <db_directory>\n")
    sys.exit(1)


def main():
    """
    main for normalize_tabmatch.py
    """
    multiquery = False
    
    if len(sys.argv) != 4:
        usage(os.path.basename(sys.argv[0]))

    normtype = int(sys.argv[1])
    queryid = sys.argv[2]
    if queryid == "-m":
        multiquery = True
    db_directory = sys.argv[3]
    if normtype < 1 or normtype > 3:
        sys.stderr.write('normtype must be 1,2, or 3\n')
        usage(sys.argv[0])

    # get numbers of SSEs from .tableauxdistmatrix files in db_diurectory
    # build dict of {name : size}
    size_dict = {}
    if (not os.path.isdir(db_directory)):
        # all tableaux/distmatrix in one db file, not dirctory of them
        dbfile = db_directory
        size_dict = parse_tableauxdb_sizes(open(dbfile))
    else:
        input_list = glob.glob(os.path.join(db_directory, '*.tableaudistmatrix'))
        if len(input_list) == 0:
            # try .angles files instead
            input_list = glob.glob(os.path.join(db_directory, '*.angles'))
            for dbfile in input_list:
                qid = os.path.splitext(os.path.splitext(os.path.basename(dbfile))[0])[0].upper()
                firstline = open(dbfile).readline()
                if len(firstline) < 1:
                    # can happen if no SSEs
               #     sys.stderr.write('skipped ' + dbfile + '\n')
                    continue
                qsize=int(firstline)
                size_dict[qid] = qsize
        else:
            for dbfile in input_list:
                idline = open(dbfile).readline()
                if len(idline) < 2:
                    # can happen if no SSEs
    #                sys.stderr.write('skipped ' + dbfile + '\n')
                    continue
                qid = idline[:8].lstrip().rstrip().upper()
                qsize = int(idline[8:])
                size_dict[qid] = qsize
 
    # build list of (name , score)
    scorelist = []
    commentlines = []
    if multiquery:
        # get list of iterables each for same queryid.
        # iter_searchresult() is isterable of tuples (queryid, score,
        # domainid) groupby() requires that the iterable already has
        # identical consecutive queryids (first element of tuple) -
        # iter_searchresult() should yield this
        query_group_iter = groupby(iter_searchresult(sys.stdin,multiquery=True,
                                                     skip_self_query=False,
                                                     negateflag=False),
                                   lambda t : t[0])

    else:
        for line in sys.stdin:
            if line[0] == '#':
                commentlines.append(line)
                continue
            s = line.split()
            if len(s) != 2:
                sys.stderr.write('skipped line ' + line)
                continue
            score_str = s[1]
            if score_str.lower() == 'nan' or score_str == '********':
                sys.stderr.write('skipped line ' + line)
                continue
            scorelist.append((s[0], float(score_str)))
        query_group_iter = [ ( queryid, [ (queryid,score,domainid)
                                        for (domainid,score) in scorelist ] ) ]

    for (queryid, result_iter) in query_group_iter:
        querysize = size_dict[queryid.upper()]

        sys.stdout.write('# normalization type ' + str(normtype) + '\n')
        sys.stdout.write('# generated by ' + os.path.basename(sys.argv[0]) +'\n')
        sys.stdout.write('# from:\n')
        for line in commentlines:
            sys.stdout.write(line)
        sys.stdout.write('#\n')

        if multiquery:
            sys.stdout.write('# QUERY ID = ' + queryid + '\n')
            
        skipcount = 0
        for (queryid, score, pdbid) in result_iter:
            try:
                dsize = size_dict[pdbid.upper()]
            except KeyError:
                skipcount += 1

            if normtype == 1:
                normscore = norm1(score, querysize, dsize)
            elif normtype == 2:
                normscore = norm2(score, querysize, dsize)
            elif normtype == 3:
                normscore = norm3(score, querysize, dsize)
            else:
                raise ValueError('unknown norm type ' + str(normtype) + '\n')
            sys.stdout.write('%s\t%20.8f\n' % (pdbid, normscore))



if __name__ == "__main__":
    main()

