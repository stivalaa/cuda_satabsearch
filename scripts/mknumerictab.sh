#!/bin/sh
#
# File:    mknumerictab.sh
# Author:  Alex Stivala
# Created: August 2008
#
# mknumerictab.sh - make table of AUC and time for Omega matrix on some folds
#
# Usage: mknumerictab.sh
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
# $Id: mknumerictab.sh 1790 2008-08-04 02:17:59Z astivala $
#

cat <<EOF
{\begin{tabular}{llrrr} 
\toprule
Fold & SCOP sid & \# SSEs & AUC & time \\\\
\midrule
EOF

for fold in *.orderpen0.typepen0.out ; do
    mktabrow.sh `basename ${fold} .out`
done | sort -t '&' -k 3,3n

cat <<EOF
\botrule
\end{tabular}}
EOF
