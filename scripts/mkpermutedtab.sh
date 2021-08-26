#!/bin/sh
#
# File:    mkpermutedtab.sh
# Author:  Alex Stivala
# Created: Mach 2009
#
# mkpermutedtab.sh - make table of average over permutations AUC of the folds
#
# Usage: mkpermutedtab.sh
#
# Output is to stdout.
#
# Each row contains:
#
#            Fold & SCOP sid & average AUC \\
# e.g.
#            Immunoglobulin      & \texttt{d1ae6h1}  & 0.93 \\
#
# $Id: mkpermutedtab.sh 2153 2009-03-28 05:20:23Z astivala $
#

AVERAGES="d1ubia_.average_auc.txt d1tttb1.average_auc.txt d1ae6h1.average_auc.txt d1bhne_.average_auc.txt d1h6rb_.average_auc.txt d2phlb1.average_auc.txt  d1tima_.average_auc.txt d1f6dc_.average_auc.txt"

cat <<EOF
{\begin{tabular}{llr} 
\hline
Fold & SCOP sid & Average AUC \\\\
\hline
EOF

for query in $AVERAGES ; do

    basefilename=${query}
    dotindex=`expr index ${basefilename} '.'`
    dotindex=`expr ${dotindex} - 1`
    scopsid=`expr substr ${basefilename} 1  ${dotindex}`

    case ${scopsid} in
        d1ubia_)
        fold='$\beta$-grasp'
        ;;
        d1ae6h1)
        fold='Immunoglobulin'
        ;;
        d1tima_)
        fold='Tim-barrel'
        ;;
        d1bhne_)
        fold='Plait (ferredoxin)'
        ;;
        d1h6rb_)
        fold='GFP-like'
        ;;
        d1tttb1)
        fold='Key-barrel'
        ;;
        d2phlb1)
        fold='Jelly-roll'
        ;;
        d1f6dc_)
        fold='NAD-binding fold'
        ;;
        *)
        fold=`echo ${scopsid} | sed 's/_/\\\_/g'`
    esac
    scopsid=`echo ${scopsid} | sed 's/_/\\\_/g'`
    avgauc=`awk '{print $NF}' ${query}`
    printf '%-20s & %-20s & %4.2f ' "${fold}" "\texttt{${scopsid}}" ${avgauc}
    echo '\\\'
done

cat <<EOF
\hline
\end{tabular}}
EOF


