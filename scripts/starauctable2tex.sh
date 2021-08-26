#!/bin/sh
#
# File:    starauctable2tex.sh
# Author:  Alex Stivala
# Created: May 2010
#
# starauctable2tex.sh -convert output of star2auctable.py to LaTeX format
#            
# Usage: starauctable2tex.sh [-n]
#
#  -n : do not include the p-value column
#
# Input on on stdin is stdout of star2auctable.py
# Output is to stdout
#
# Uses options specific to GNU sort
#
# $Id: starauctable2tex.sh 3885 2010-07-07 07:45:07Z alexs $
#

if [ $# -gt 1 ]; then
    echo "Usage: $0" 2>&1
    exit 1
fi
incpvalue=1
if [ $# -eq 1 ]; then
    if [ $1 = "-n" ]; then
        incpvalue=0
    else
        echo "Usage: $0" 2>&1
        exit 1
    fi
fi
     

if [ $incpvalue -eq 0 ]; then
    cat <<EOF
{\begin{tabular}{lr}  \hline
Method(s) & $\Delta\mathrm{AUC}$ \\\\
\hline
EOF
else
    cat <<EOF
{\begin{tabular}{lrr}  \hline
Method(s) & $\Delta\mathrm{AUC}$ & p-value \\\\
\hline
EOF
fi

TAB=`printf "\t"`
if [ $incpvalue -eq 0 ]; then
    sort -k4,4n -t"$TAB" | awk -vFS="$TAB" '{printf("%-40s & %s  \\\\\n",$1,$4)}' | sed 's/IRTableau-F77/IR Tableau/'
else
    sort -k4,4n -t"$TAB" | awk -vFS="$TAB" '{printf("%-40s & %s & %s \\\\\n",$1,$4,$3)}' | sed 's/IRTableau-F77/IR Tableau/'
fi


cat <<EOF
\hline
\end{tabular}}
EOF

