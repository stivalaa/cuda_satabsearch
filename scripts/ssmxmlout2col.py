#!/usr/bin/env python
#
# File:    ssmxmlout2col.sh
# Author:  Alex Stivala
# Created: November 2008
#
# ssmxmlout2col.sh - Convert SSM webserver XML output to 2 column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: ssmxmlout2col.sh < domain.xml
# 
# Output has two columns, database id and SSM Pcli score
#
# Output is to stdout.
#
# Uses the XML output format from  the SSM webserver
# http://www.ebi.ac.uk/msd-srv/ssm/
#
# Developed with SSM v2.36 output
# 
# $Id: ssmxmlout2col.py 2103 2009-03-16 05:15:19Z astivala $
#

import os,sys
from xml.dom import minidom

if (len(sys.argv) != 1):
    sys.stderr.write('Usage: ' + sys.argv[0] + '\n')
    sys.exit(1)
    
doc = minidom.parse(sys.stdin)
matches = doc.getElementsByTagName("Match")
for match in matches:
    qscore = [child for child in match.childNodes
              if child.nodeType == child.ELEMENT_NODE and
                 child.nodeName == "Q-score"][0]
    qval = qscore.firstChild.data
    target = [child for child in match.childNodes
              if child.nodeType == child.ELEMENT_NODE and
                 child.nodeName == "Target"][0]
    name =  [child for child in target.childNodes
              if child.nodeType == child.ELEMENT_NODE and
                 child.nodeName == "name"][0]
    sid = name.firstChild.data
    sys.stdout.write('%s    %s\n' % (sid, qval))    
        
        

