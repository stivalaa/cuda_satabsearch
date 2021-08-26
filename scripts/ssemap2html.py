#!/usr/bin/env python
###############################################################################
#
# ssemap2html.py - Convert ssemap format (e.g. soln2ssemap.py output) to HTML
#                     in the structure.
#
# See usage comment in docstring for main()
#
# File:    ssemap2html.py
# Author:  Alex Stivala
# Created: November 2008
#
# $Id: ssemap2html.py 2038 2008-11-25 05:39:26Z astivala $
# 
###############################################################################

"""
Converts the ssemap format from soln2ssemap.py or ssesubsetreamp.py
or tlocsd to HTML formatted text, with each dbid given link to
pro-origami webserver database prebuilt cartoon for that ASTRAL SCOP
sid, with the selected SSEs indicated to be highlighted, and link to
SCOP entry, also SCOP superfamily sccs and fold description.

Output is to stdout.

The cache file is a python pickle dictionary:
     scopdominfo_dict -
      dict {sid: (superfamily_sccs, superfamily_description, fold_sccs, fold_description)}
      where
     superfamily_sccs is SCOP sccs identifying the superfamily for the domain
     superamily_description is SCOP dessription of the superfamily
     fold_sccs is SCOP sccs of the fold it is in
     fold_description is the SCOP descriptino of the fold the domain is in

"""

import sys,os
import getopt
import pickle

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
    
    sys.stderr.write("Usage: " +progname + " cachefile\n")
    sys.exit(1)


def querymap_cmp(qmap1, qmap2):
    """
    compare two QuerySSEMap objects by score for use in sorting.
    uses absolute value so works with QP tableau search (-ve) and
    heuristic tableau search (+ve) scores
    """
    if abs(qmap1.score) < abs(qmap2.score):
        return -1
    elif abs(qmap1.score) > abs(qmap2.score):
        return  1
    else:
        return 0
    
def main():
    """
    main for ssemap2html.py

    Usage: ssemap2html.py cachefile


    cachefile is the filename of the cache (pickled) file built by
    build_fastscopdominfo_cache.py

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

    Ouput is HTML format text with each dbid given link to pro-origami
    webserver database prebuilt cartoon for that ASTRAL SCOP sid, with
    the selected SSEs indicated to be highlighted, and link to SCOP
    entry, also superfamily sccs id with link to SCOP and fold
    description.

    """
    dbdomid = None
    
    if len(sys.argv) != 2:
        usage(os.path.basename(sys.argv[0]))

    pickle_filename = sys.argv[1]
    scopdominfo_dict = pickle.load(open(pickle_filename))

    search_maps = parse_ssemap(sys.stdin)

    print '<html>'
    print '<link rel="stylesheet" href="/pro-origami/style.css" />'
    print '<div id="qpresults">'
    print '<table>'

    print '<tr><th>score<th>cartoon<th>SCOP entry<th>superfamily<th>fold</tr>'

    for query_ssemap in sorted(search_maps.query_ssemap_list,cmp=querymap_cmp,reverse=True):
        if len(query_ssemap.sse_map) == 0:
            sseseqnums = None
            sseseqnums_str = "none"
        else:
            sseseqnums = [j for (i,j) in query_ssemap.sse_map]
            sseseqnums_str = ','.join([str(s) for s in sseseqnums])

        sid = query_ssemap.domid
        entry = scopdominfo_dict[sid]
        sf_sccs = entry[0]
        sf_desc = entry[1]
        fold_sccs = entry[2]
        fold_desc = entry[3]

        sys.stdout.write("<tr><td>%s <td><a href=\"/pro-origami/cgi-bin/podbget.cgi?pdbcode=%s&format=SVGINTERACTIVE&selsses=%s\">%s</a><td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid=%s\">%s</a><td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs=%s\">%s %s<td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs=%s\">%s</tr>\n"
                         % (query_ssemap.score,
                            query_ssemap.domid,
                            sseseqnums_str,
                            query_ssemap.domid, query_ssemap.domid,
                            query_ssemap.domid,
                            sf_sccs, sf_sccs, sf_desc,
                            fold_sccs, fold_desc)
                        )

    print '</table>'
    print '</div>'
    print '</html>'

            
if __name__ == "__main__":
    main()
