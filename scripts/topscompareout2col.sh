#!/bin/sh
#
# File:    topscompareout2col.sh
# Author:  Alex Stivala
# Created: March 2009
#
# topscompareout2col.sh - Convert tops_comparison output format to 2-column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: topscompareout2col.sh < topscompareoutput
#
#
# Output has two columns, database id and tops_comparison compressino score
#
# Output is to stdout.
#
# $Id: topscompareout2col.sh 2161 2009-03-29 03:27:50Z astivala $
#

awk '$2 != "probe" {print substr($2,1,7),$1}'
