#!/bin/sh
#
# File:    mkfischer3dtabs_all.sh
# Author:  Alex Stivala
# Created: October 2008
#
# mkfischer3dtab.sh - make all LaTeX tables of results for MSVNS and
#                     QP tableau search on Fischer and Nh3D data sets
#            
# Uses the mkfischer3dtab.sh script to make ecah table
# Muast be run at directoryu containg results/ and maxcmo_results/ dirs
#
# Usage: mkfischer3dtabs_all.sh
#
# $Id: mkfischer3dtabs_all.sh 2996 2009-11-30 05:25:05Z alexs $
#

if [ $# -ne 0 ]; then
    echo "Usage: $0 " 2>&1
    exit 1
fi

tmpfile=/var/tmp/f3tab.$$

for dataset in fischer nh3d
do
    for level in fold class arch
    do
        if [ \( $dataset = "fischer" -a $level = "arch" \) -o \( $dataset = "nh3d" -a $level = "fold" \) ]; then
            continue
        fi
        outfile=results/${dataset}.${level}.textab
        cat <<EOF >${outfile}
\begin{tabular}{llrrrr} 
\hline
 & & & standard & \multicolumn{2}{c}{95\% confidence interval} \\\\
Method & Normalization & AUC & error & lower & upper \\\\
\hline
EOF
        mkfischer3dtab.sh $dataset $level >  $tmpfile
        tabcolmax.sh $tmpfile 3 >> ${outfile}
        cat <<EOF >>${outfile}
\hline
\end{tabular}
EOF
   done
done
rm ${tmpfile}
