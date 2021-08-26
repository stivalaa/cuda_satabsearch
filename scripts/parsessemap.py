###############################################################################
#
# parsessemap.py - Functions to parse soln2ssemap.py output
#
# File:    parsessemap.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id: parsessemap.py 3120 2009-12-24 04:42:47Z alexs $
# 
###############################################################################

"""
Parse output of soln2ssepmap.py
"""
import sys,os



#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class SearchMap:
    """
    SearchMap is just a dummy class for containing the search results
    with solution vectors, returned by parse_ssemap()
    """
    pass

class QuerySSEMap:
    """
    QuerySSEMap is a dummy class for containign result from individual query,
    in SearchMap.
    """
    pass

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def parse_ssemap(fh):
    """
    Parse the output of soln2ssemap.py;
    identifier and score (as per input), then
    for each matching a line containing
    i and j separated by a space,
    one per line (with blank line before next id) e.g.:

    # soln2ssemap.py processed:
    #
    # TSRCHD LTYPE = T LORDER = F LSOLN = T
    # QUERY ID = D1KI9A_ 
    # DBFILE = /home/astivala/tableauxdb.ascii
    # Tue Aug  5 12:31:50 2008
    #
    # on 05Aug2008 15:36:05
    #
    d1wiua_     -23.0000
    1 1
    3 2
    8 4
    9 5
    11 6
    14 9

    Note the header 'comment' information is required, we get the QUERY ID
    from it.

    Parameters:
        fh - open (read) filehandle to parse from
        
    Return value:
        search_maps - dummy class SearchMap containing:
             queryid - query identifier parsed from comments
             query_ssemap_list - list of dummy class QuerySSEMap containing:
                 domid - id of domain in db
                 score - score of queryid against domid
                 sse_map - list of (i,j) SSE sequential index tuples
             comment_lines - list of comment lines read

    """
    search_maps = SearchMap()
    query_ssemap_list = []
    query_ssemap = None
    search_maps.comment_lines = []
    for line in fh:
        if line[0] == '#':
            sline = line[1:].split('=')
            if sline[0].lstrip().rstrip() == 'QUERY ID':
                search_maps.queryid = sline[1].lstrip().rstrip().lower()
            search_maps.comment_lines.append(line)
            continue
        elif len(line) < 2:
            continue #skip blank lines
        elif not line.split()[0].isdigit(): # new identifier
            if query_ssemap:
                query_ssemap_list.append(query_ssemap) # finsished with prev one
            splitline = line.split()
            if len(splitline) != 2:
                sys.stderr.write('bad line: ' + line + '\n')
                continue
            domainid = splitline[0]
            score_str = splitline[1]
            query_ssemap = QuerySSEMap()
            query_ssemap.domid = domainid
            query_ssemap.sse_map = []
            query_ssemap.score = float(score_str)
        else:
            # line should be two SSE indices spearated by space
            try:
                sline = line.split()
                i = int(sline[0])
                j = int(sline[1])
            except ValueError:
                sys.stderr.write('bad line: ' + line + '\n')
                # skip whole query of a value is bad
                if query_ssemap:
                    sys.stderr.write('skipping ' + query_ssemap.domid + '\n')
                    query_ssemap = None
                continue
            query_ssemap.sse_map.append((i,j))
    if query_ssemap:
        query_ssemap_list.append(query_ssemap)
        
    search_maps.query_ssemap_list = query_ssemap_list
    return search_maps

