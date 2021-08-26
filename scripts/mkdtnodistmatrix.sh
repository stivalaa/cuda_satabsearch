#!/bin/sh
#
# File:    mkdtnodistmatrix.sh
# Author:  Alex Stivala
# Created: August 2008
#
# mkdtnodistmatrix.sh - make table of AUC and time for tableaux on some folds
#
# Usage: mkdtnodistmatrix.sh
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
# $Id: mkdtnodistmatrix.sh 1857 2008-09-05 07:34:16Z astivala $
#

RTABS_DISCRETE="d1ubia_.nodistmatrix.tsrchd.tt.rtab d1ae6h1.nodistmatrix.tsrchd.tt.rtab d1bhne_.nodistmatrix.tsrchd.tt.rtab d1h6rb_.nodistmatrix.tsrchd.tt.rtab d1tttb1.nodistmatrix.tsrchd.tt.rtab d1tima_.nodistmatrix.tsrchd.tt.rtab d2phlb1.nodistmatrix.tsrchd.tt.rtab d1f6dc_.nodistmatrix.tsrchd.tt.rtab"

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
