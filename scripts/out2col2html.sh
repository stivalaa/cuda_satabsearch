#!/bin/sh
#
# File:    out2col2html.sh
# Author:  Alex Stivala
# Created: October 2009
#
# out2col2html.sh - Convert 2-colum (dbid score) format to HTML format
#
# Usage: out2col2html.sh < qptabsearchoutput
# 
# Output from qp tableau search has two columns, database id and matching score.
# This converts to HTML formatted text, with each dbid given link to
# pro-origami webserver database prebuilt cartoon for that ASTRAL SCOP sid,
# and to SCOP entry.
#
# Output is to stdout.
#
# $Id: out2col2html.sh 3090 2009-12-20 06:06:14Z alexs $
#

echo '<html>'
echo '<link rel="stylesheet" href="/pro-origami/style.css" />'
echo '<div id="qpresults">'
echo '<table>'
awk '{printf("<tr><td>%s <td><a href=\"/pro-origami/cgi-bin/podbget.cgi?pdbcode=%s&format=SVGINTERACTIVE\">%s</a><td><a href=\"http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid=%s\">SCOP entry for %s</a></tr>\n"),$2,$1,$1,$1,$1}'
echo '</table>'
echo '</div>'
echo '</html>'


