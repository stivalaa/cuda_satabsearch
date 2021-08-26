#!/bin/sh
#
# File:    starauctable2tex.sh
# Author:  Alex Stivala
# Created: May 2010
#
# starauctable2tex.sh -convert output of star2auctable.py to LaTeX format
#            
# Usage: starauctable2tex.sh
#
# Input on on stdin is stdout of star2auctable.py
# Output is to stdout
#
# Uses options specific to GNU sort
#
# $Id: starauctable2tex.sh 3682 2010-05-17 06:05:18Z alexs $
#

if [ $# -ne 0 ]; then
    echo "Usage: $0" 2>&1
    exit 1
fi


cat <<EOF
{\begin{tabular}{lrr}  \hline
Method(s) & $\Delta\mathrm{AUC}$ & p-value \\\\
\hline
EOF

TAB=`printf "\t"`
sort -k4,4n -t"$TAB" | awk -vFS="$TAB" '{printf("%-40s & %s & %s \\\\\n",$1,$4,$3)}'


cat <<EOF
\hline
\end{tabular}}
EOF

