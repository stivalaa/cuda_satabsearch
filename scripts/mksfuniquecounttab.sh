#!/bin/sh
#
# File:    mksfuniquecounttab.sh
# Author:  Alex Stivala
# Created: Mach 2009
#
# mksfuniquecounttab.sh - make table of superfamily counts unique to each method
#
# Usage: mksfuniquecounttab.sh
#
# Output is to stdout.
# Must be run from other_results directory after make has been run
# to build .sflist files in all subsidiary directories that are used.
#
# Each row contains:
#
#            Fold & SCOP sid & Pu  & Su  & Tu  & Qu  & R & Pu/R & Su/R & Tu/R & Qu/R \\
# e.g.
#   Immunoglobulin & \texttt{d1ae6h1}  & 27 & 1 & 2 & 3 & & & & & \\
#
# Note the SCOP column (number explicitly mentioned in SCOP ) and following
# depend on the .scopmentioned files WHICH MUST BE GENERATED MANUALLY
# (see README), so if any of the .sflist files actally change, 
# the .scopmentioned files must be manually updated.
#
# WARNING: this script creates the $scopsid.$methodabbrev.uniquesfinscop files
# in the cwd, will overwrite if they exist
#
# $Id: mksfuniquecounttab.sh 2186 2009-04-01 00:03:31Z astivala $
#

RESULTS_DIRS="ProSMoS SSM tops/folds ../results"
SF_LISTS="d1ubia_.sflist d1tttb1.sflist d1ae6h1.sflist d1bhne_.sflist d1h6rb_.sflist d2phlb1.sflist  d1tima_.sflist d1f6dc_.sflist"


# write list of lines in one file and not in all others to stoud
# Parameters:
#     ufile - filename of file to find unique lines in
#     subsequence paramters - filenames of other files to combine all
#                             lines together, results are lines
#                             in ufile and not in this combined list
uniquelines() {
    ufile=$1
    shift 1
    others=$*
    tmpufile=/var/tmp/msu$$.u
    tmpothers=/var/tmp/msu$$.others
    sort $ufile > $tmpufile
    cat $others | sort | uniq > $tmpothers
    comm -23 $tmpufile $tmpothers
    rm $tmpufile $tmpothers
}


# write common lines between two files to stdout
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

#
# write abbreviation for method name based on results directory
# Paramters:
#     results directory name
getmethodabbrev() {
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
    echo $methodabbrev
}


#############################################################################
#
# Main
#
############################################################################# 

cat <<EOF
{\begin{tabular}{llrrrrrrrrr} 
\hline
EOF
printf "Fold & SCOP sid "
for resdir in ${RESULTS_DIRS} ; do
    methodabbrev=`getmethodabbrev ${resdir}`
    printf "& ${methodabbrev}u "
done
printf " & R & Pu/R & Su/R & Tu/R & Qu/R "
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
    sid=${scopsid}
    scopsid=`echo ${scopsid} | sed 's/_/\\\_/g'`
    printf '%-20s & %-20s ' "${fold}" "\texttt{${scopsid}}"
    for resdir in ${RESULTS_DIRS} ; do
        others=""
        for oresdir in ${RESULTS_DIRS} ; do
            if [ "$oresdir" != "$resdir" ] ; then
                others="$others ${oresdir}/${sflist}"
            fi
        done
        uniquesfcount=`uniquelines ${resdir}/${sflist} ${others} | wc -l`
        printf '& %3d ' ${uniquesfcount}
    done

    scopmentioned=`basename ${sflist} .sflist`.scopmentioned

    # build the .uniquesfinscop files with sccs ids found only by
    # each method, that are in the .scopmentioned list
    for resdir in ${RESULTS_DIRS} ; do
        others=""
        for oresdir in ${RESULTS_DIRS} ; do
            if [ "$oresdir" != "$resdir" ] ; then
                others="$others ${oresdir}/${sflist}"
            fi
        done
        methodabbrev=`getmethodabbrev ${resdir}`
        tmpuniquesf=/var/tmp/${sid}.${methodabbrev}.uniquesf
        uniquelines ${resdir}/${sflist} ${others} > $tmpuniquesf
        uniquesfinscop=${sid}.${methodabbrev}.uniquesfinscop
        commonlines ${tmpuniquesf} ${scopmentioned} > $uniquesfinscop
        rm $tmpuniquesf
    done

    # count the total number of superfamilies that are uniquely
    # found for each method, and which is also mentioned in SCOP:
    # this is just the total number in all .uniquesfinscop files
    # since they are disjoint sets (an entry in one cannot be in any
    # other since they are those found ONLY by each method)
    uniquescopcount=`wc -l ${sid}.*.uniquesfinscop | grep total | awk '{print $1}'`
    printf '& %3d '  $uniquescopcount

    # now count the unique superfamilies and unique superfamilies mentioned
    # in scop in the files just generated
    # for each method and write them in the table
    for resdir in ${RESULTS_DIRS} ; do
        methodabbrev=`getmethodabbrev ${resdir}`
        uniquesfinscop=${sid}.${methodabbrev}.uniquesfinscop
        uniquesf_in_r_count=`wc -l ${uniquesfinscop} | awk '{print $1}'`
        printf '& %s ' "${uniquesf_in_r_count}/${uniquescopcount}"
    done

    echo '\\\'
done

cat <<EOF
\hline
\end{tabular}}
EOF

