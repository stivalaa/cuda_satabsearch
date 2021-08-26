#!/bin/sh
#
# File:    mkfoldstab2.sh
# Author:  Alex Stivala
# Created: September 2008
#
# mkfoldstab2.sh - make table of AUC and time for tableaux on some folds
#                  for two variations of method (with+without distance 
#                  difference constraints)
#
# Usage: mkfoldstab2.sh
#
# Output is to stdout.
#
# Each row contains:
#
#   Fold & SCOP sid & \# SSEs & AUC (without dist) & time (without dist) & AUC (with dist) & time (with dist) \\
# e.g.
# $\beta$-grasp        & \texttt{d1ubia\_}    &   8 & 0.80 & 0 h 47 m & 0.92 & 0 h 32 m 
#
# The mktabrow.sh script is called to make each row.
# The table is sorted by #SSEs
#
# Uses options specific to GNU sort
# Note may have problems with
# multiple backslash escapes used to format LaTeX table
# - this works on room0219pos09.cs.mu.oz.au
# but not charikar.cs.mu.oz.au (both Linux)!
#
# $Id: mkfoldstab2.sh 2204 2009-04-06 00:23:10Z astivala $
#

RTABS_DISCRETE="d1ubia_.tsrchd.tt.rtab d1ae6h1.tsrchd.tt.rtab d1bhne_.tsrchd.tt.rtab d1h6rb_.tsrchd.tt.rtab d1tttb1.tsrchd.tt.rtab d1tima_.tsrchd.tt.rtab d2phlb1.tsrchd.tt.rtab d1f6dc_.tsrchd.tt.rtab"

RTABS_NODISTMATRIX_DISCRETE="d1ubia_.nodistmatrix.tsrchd.tt.rtab d1ae6h1.nodistmatrix.tsrchd.tt.rtab d1bhne_.nodistmatrix.tsrchd.tt.rtab d1h6rb_.nodistmatrix.tsrchd.tt.rtab d1tttb1.nodistmatrix.tsrchd.tt.rtab d1tima_.nodistmatrix.tsrchd.tt.rtab d2phlb1.nodistmatrix.tsrchd.tt.rtab d1f6dc_.nodistmatrix.tsrchd.tt.rtab" 

ftab_discrete_tmpfile=/tmp/ftab$$
ftab_nodistmatrix_discrete_tmpfile=/tmp/ftabnd$$

cat <<EOF
{\begin{tabular}{llrrrrr} 
\hline
 &&&\multicolumn{4}{c}{distance information} \\\\
 &&&\multicolumn{2}{c}{without} & \multicolumn{2}{c}{with} \\\\
 Fold & SCOP sid & \# SSEs & AUC & time & AUC & time \\\\

\hline
EOF

for fold in $RTABS_DISCRETE ; do
    mktabrow.sh `basename ${fold} .rtab`
done | sort -t '&' -k2,2 > ${ftab_discrete_tmpfile}

for fold in $RTABS_NODISTMATRIX_DISCRETE ; do
    mktabrow.sh `basename ${fold} .rtab`
done | sort -t '&' -k2,2 > ${ftab_nodistmatrix_discrete_tmpfile}

join -t\& -12 -22  ${ftab_nodistmatrix_discrete_tmpfile} ${ftab_discrete_tmpfile} | sed 's/\\\\//g' |awk -F\& '{print $2 "&" $1 "&" $3 "&" $4 "&" $5 "&"  $8 "&" $9 "\\\\"}' | sort -t '&' -k3,3n

cat <<EOF
\hline
\end{tabular}}
EOF

rm ${ftab_discrete_tmpfile} ${ftab_nodistmatrix_discrete_tmpfile}

