#!/bin/sh
#
# File:    mkdiscretetab.sh
# Author:  Alex Stivala
# Created: August 2008
#
# mkdiscretetab.sh - make table of AUC and time for tableaux on some folds
#
# Usage: mkdiscretetab.sh
#
# Output is to stdout.
#
# Each row contains:
#
#            Fold & SCOP sid & \# SSEs & AUC & time\\
# e.g.
#            Immunoglobulin      & \texttt{d1ae6h1}  & 13 & 0.93 & 9:20  \\
#
# The mktabrow.sh script is called to make each row.
# The table is sorted by #SSEs
#
# Uses options specific to GNU sort
#
# $Id: mkdiscretetab.sh 1790 2008-08-04 02:17:59Z astivala $
#

RTABS_DISCRETE="d1ubia_.tsrchd.tt.rtab d1ae6h1.tsrchd.tt.rtab d1bhne_.tsrchd.tt.rtab d1h6rb_.tsrchd.tt.rtab d1tttb1.tsrchd.tt.rtab d1tima_.tsrchd.tt.rtab d2phlb1.tsrchd.tt.rtab d1f6dc_.tsrchd.tt.rtab"

cat <<EOF
{\begin{tabular}{llrrr} 
\toprule
Fold & SCOP sid & \# SSEs & AUC & time \\\\
\midrule
EOF

for fold in $RTABS_DISCRETE ; do
    mktabrow.sh `basename ${fold} .rtab`
done | sort -t '&' -k 3,3n

cat <<EOF
\botrule
\end{tabular}}
EOF
