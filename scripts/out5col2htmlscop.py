#!/usr/bin/env python
###############################################################################
#
# out5col2htmlscop.py - Convert 5 column format
#                       (domainid,rawscore,norm2score,zscore,pvalue)
#                       to HTML with SCOP information and links to SCOP
#
# See usage comment in docstring for main()
#
# File:    out5col2htmlscop.py
# Author:  Alex Stivala
# Created: March 2010
#
# $Id: out5col2htmlscop.py 2038 2008-11-25 05:39:26Z astivala $
# 
###############################################################################

"""
Converts the 5 column output from qptabsearch, SA tab search, etc.
HTML formatted text, with each dbid given link to
pro-origami webserver database prebuilt cartoon for that ASTRAL SCOP
sid, with the selected SSEs indicated to be highlighted, and link to
SCOP entry, also SCOP superfamily sccs and fold description.

Output is to stdout.

The cache file is a python pickle dictionary:
     scopdominfo_dict -
      dict {sid: (superfamily_sccs, superfamily_description, fold_sccs fold_description)}
      where
     superfamily_sccs is SCOP sccs identifying the superfamily for the domain
     superamily_description is SCOP dessription of the superfamily
     fold_sccs is SCOP sccs of the fold it is in
     fold_description is the SCOP descriptino of the fold the domain is in

"""

import sys,os
import pickle



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


def main():
    """
    main for out5col2htmlscop.py

    Usage: out5col2htmlscop.py cachefile


    cachefile is the filename of the cache (pickled) file built by
    build_fastscopdominfo_cache.py

    Input is on stdin, the output of soln2ssemap.py,
    format is:

    domainid rawscore norm2score zscore pvalue

    Ouput is HTML format text with each dbid given link to pro-origami
    webserver database prebuilt cartoon for that ASTRAL SCOP sid, with
    the selected SSEs indicated to be highlighted, and link to SCOP
    entry, also superfamily sccs id with link to SCOP and fold
    description.

    """
    
    if len(sys.argv) != 2:
        usage(os.path.basename(sys.argv[0]))

    pickle_filename = sys.argv[1]
    scopdominfo_dict = pickle.load(open(pickle_filename))


    print '<html>'
    print '<link rel="stylesheet" href="/pro-origami/style.css" />'
    print '<div id="qpresults">'
    print '<table>'

    print '<tr><th>score<th>cartoon<th>SCOP entry<th>superfamily<th>fold</tr>'
#    print '<tr><th>p-value<th>cartoon<th>SCOP entry<th>superfamily<th>fold</tr>'

    for line in sys.stdin:
        splitline = line.split()
        sid = splitline[0]
        pvalue = splitline[4]
        score = splitline[1]
        entry = scopdominfo_dict[sid]
        sf_sccs = entry[0]
        sf_desc = entry[1]
        fold_sccs = entry[2]
        fold_desc = entry[3]

        sys.stdout.write("<tr><td>%s <td><a href=\"/pro-origami/cgi-bin/podbget.cgi?pdbcode=%s&format=SVGINTERACTIVE\">%s</a><td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid=%s\">%s</a><td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs=%s\">%s %s<td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs=%s\">%s</tr>\n"
                         % (score,
                            sid,
                            sid, sid,
                            sid,
                            sf_sccs, sf_sccs, sf_desc,
                            fold_sccs, fold_desc)
                        )

    print '</table>'
    print '</div>'
    print '</html>'

            
if __name__ == "__main__":
    main()
