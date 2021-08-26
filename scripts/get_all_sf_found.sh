#!/bin/sh
#
# File:    get_all_sf_found.sh
# Author:  Alex Stivala
# Created: Mach 2009
#
# get_all_sf_found.sh - make list of all superfamilies found by all methods
#
# Usage: get_all_sf_found.sh
#
# Output is to the scopsid.allsffound files in the cwd
# WARNING: these files are overwritten if they exist
# Must be run from other_results directory after make has been run
# to build .sflist files in all subsidiary directories that are used.
#
#
# $Id: get_all_sf_found.sh 2169 2009-03-30 05:33:43Z astivala $
#

RESULTS_DIRS="ProSMoS SSM tops/folds ../results"
SF_LISTS="d1ubia_.sflist d1tttb1.sflist d1ae6h1.sflist d1bhne_.sflist d1h6rb_.sflist d2phlb1.sflist  d1tima_.sflist d1f6dc_.sflist"



for sflist in ${SF_LISTS} ; do
    scopsid=`basename ${sflist} .sflist`
    tmpfile=/var/tmp/gsf$$.${scopsid}
    cat /dev/null > ${tmpfile}
    for resdir in ${RESULTS_DIRS} ; do
        cat ${resdir}/${sflist} >> ${tmpfile}
    done
    sort ${tmpfile} | uniq > ${scopsid}.allsffound
    rm ${tmpfile}
done

