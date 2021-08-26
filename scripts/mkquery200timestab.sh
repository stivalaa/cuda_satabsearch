#!/bin/sh
#
# File:    mkquery200timestab.sh
# Author:  Alex Stivala
# Created: March 2009
#
# mkquery200timestab.sh - make LaTeX table of AUC values and timne for 
#             different QP tableau search implementations on query200 data set
#            
#
# Usage: mkquery200timestab.sh
#
# Uses sumtimes.sh script to get total time from all .err files
#
# Uses options specific to GNU sort
#
# $Id: mkquery200timestab.sh 3909 2010-07-11 04:45:25Z alexs $
#

if [ $# -ne 0 ]; then
    echo "Usage: $0" 2>&1
    exit 1
fi

cat <<EOF
\begin{tabular}{lrrrr} 
\hline
Method & min. SSE length & AUC & time & speedup  \\\\
\hline
EOF

baseline_s=0

for statsfile in ../tango_results_umfpack/query200/norm2/*.stats ../tango_results_pardiso/query200/norm2/*.stats ../tango_results_ma57/query200/norm2/*.stats ../tango_results_umfpack/query200-minlen4/norm2/*.stats  #../tango_results/localsearch/query200/norm2/*.stats
do
    resdir=`echo ${statsfile} | cut -d/ -f2`
    if [ ${resdir} = "tango_results_pardiso" ]; then
        solver="PARDISO"
    elif [ ${resdir} = "tango_results_ma57" ]; then
        solver="MA57"
    elif [ ${resdir} = "tango_results_umfpack" ]; then
        solver="UMFPACK"
    else
        solver="simulated annealing"
    fi
    resdir2=`echo ${statsfile} | cut -d/ -f3`
    if [ ${resdir2} = "query200-minlen4" ]; then
        minlen=4
    else
        minlen=1
    fi
    auc=`fgrep 'RROC          AUC' ${statsfile} | cut -d= -f2`
    errdir=`dirname ${statsfile}`
    errdir=`dirname ${errdir}`
#    echo XXX $errdir
    hms=`sumtimes.sh ${errdir}/*.err`
    h=`echo $hms | cut -d' ' -f1`
    m=`echo $hms | cut -d' ' -f3`
    s=`echo $hms | cut -d' ' -f5`
    total_s=`expr $h \* 3600 + $m \* 60 + $s`
    if [ $baseline_s -eq 0 ]; then
        baseline_s=$total_s
    fi
    speedup=`echo "$baseline_s / $total_s" | bc -l`
    printf '%-22s & %2d & %5.2f & %s & %8.2f \\\\\n' "${solver}" "${minlen}" ${auc} "${hms}"  ${speedup}
done | sort -t'&'  -k4,4nr

cat <<EOF
\hline
\end{tabular}
EOF
