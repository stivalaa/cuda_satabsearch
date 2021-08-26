#!/bin/sh
#
# File:    mksfcounttab.sh
# Author:  Alex Stivala
# Created: Mach 2009
#
# mksfcounttab.sh - make table of  superfamily counts
#
# Usage: mksfcounttab.sh
#
# Output is to stdout.
# Must be run from other_results directory after make has been run
# to build .sflist files in all subsidiary directories that are used.
#
# Each row contains:
#
#                 &          & Number of superfamilies \\
#            Fold & SCOP sid & ProSMoS & SSM & TOPS & QP tabsearch & SCOP & P/R & S/R & T/R & Q/R \\
# e.g.
#   Immunoglobulin & \texttt{d1ae6h1}  & 27 & 1 & 2 & 3 & & & & & \\
#
# Note the SCOP column (number explicitly mentioned in SCOP ) and following
# depend on the .scopmentioned files WHICH MUST BE GENERATED MANUALLY
# (see README), so if any of the .sflist files actally change, 
# the .scopmentioned files must be manually updated.
#
# $Id: mksfcounttab.sh 2181 2009-03-31 03:17:24Z astivala $
#

RESULTS_DIRS="ProSMoS SSM tops/folds ../results"
SF_LISTS="d1ubia_.sflist d1tttb1.sflist d1ae6h1.sflist d1bhne_.sflist d1h6rb_.sflist d2phlb1.sflist  d1tima_.sflist d1f6dc_.sflist"



# write common lines between two files to stdount
# Parameters:
#     file1 - filename of first file
#     file2 - filename of secondfile
commonlines() {
    file1=$1
    file2=$2
    tmpfile1=/var/tmp/msf$$.1
    tmpfile2=/var/tmp/msf$$.2
    sort $file1 > $tmpfile1
    sort $file2 > $tmpfile2
    comm -12 $tmpfile1 $tmpfile2
    rm ${tmpfile1} ${tmpfile2}
}


cat <<EOF
{\begin{tabular}{llrrrrrrrrr} 
\hline
EOF
printf "Fold & SCOP sid "
for resdir in ${RESULTS_DIRS} ; do
    if [ `basename ${resdir}` = "results" ]; then
        method="QP tableau search"
    else
        method=`echo ${resdir} | cut -d/ -f1`
        if [ "${method}" != "TableauSearch" -a "${method}" != "ProSMoS" ]; then
            method=`echo "${method}" | tr a-z A-Z`
        fi
    fi
    case "${method}" in
        ProSMoS)
            methodabbrev='P'
            ;;
        SSM)
            methodabbrev='S'
            ;;
        TOPS)
            methodabbrev='T'
            ;;
        "QP tableau search")
            methodabbrev='Q'
            ;;
    esac
    printf "& ${methodabbrev} "
done
printf " & R & P/R & S/R & T/R & Q/R "
cat <<EOF
\\\\
\hline
EOF

for sflist in ${SF_LISTS} ; do
    scopsid=`basename ${sflist} .sflist`
    case ${scopsid} in
        d1ubia_)
        fold='$\beta$-grasp'
        ;;
        d1ae6h1)
        fold='Immunoglobulin'
        ;;
        d1tima_)
        fold='TIM-barrel'
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
    printf '%-20s & %-20s ' "${fold}" "\texttt{${scopsid}}"
    for resdir in ${RESULTS_DIRS} ; do
        sfcount=`wc -l ${resdir}/${sflist} | cut -d' ' -f1`
        printf '& %3d ' ${sfcount}
    done
    scopmentioned=`basename ${sflist} .sflist`.scopmentioned
    scopcount=`wc -l ${scopmentioned} | awk '{print $1}'`
    printf " & %3d " ${scopcount} 
    for resdir in ${RESULTS_DIRS} ; do
        sf_in_r_count=`commonlines ${resdir}/${sflist} ${scopmentioned} | wc -l`
        printf '& %s ' "${sf_in_r_count}/${scopcount}"
    done
    echo '\\\'
done

cat <<EOF
\hline
\end{tabular}}
EOF


