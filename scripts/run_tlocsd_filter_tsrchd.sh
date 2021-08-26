#!/bin/bash
###############################################################################
#
# run_tlocsd_filter_tsrchd.sh - Run tlocsd heuristic (fast) then use top
#                               10% hits only as db for slower tsrchd
#
# File:    run_tlocsd_filter_tsrchd.sh
# Author:  Alex Stivala
# Created: November 2009
#
# For the top 10% (or other) hits according to tlocsd simulated annealing
# use that as database for slower QP tableau matching.
# 
# Usage:
#     run_tlocsd_filter_tsrchd.sh [-p percent] [-e sse_num_list] [-q] querystruct dbdir
#
#     -q: do not use the ordering constraint (allow nonsequential  matchings)
#
#     -e sse_num_list: list of SSE sequential numbers to select from
#                      query struct rather than whole structure
#
#     -p percent_hits : percentage of top hits to use (default 10)
#
#     querystruct is a structure in PDB format for the query
#     dbdir is the tableaux+distmatrix db directory. It must contain
#     the files:
#               distmatrixdb.pickle
#               tableauxdb.pickle
#               tableauxdistmatrixdb.ascii
#
#     where the .pickle files were created with buildtableauxdb.py
#     and the .ascii file is created from them with convdb2.py
#
# Output to stdout is the output from tsrchd
#
# The script works by first building the query tableau+distmatrix
# with pytableaucreate.py then using this as input to the fast
# tlocsd program. The top 10% (or n%) of the hits are then used
# to select only those structures from the database to build a
# new database with only those, and then tsrchd_pardiso run on that
# reduced database.
# TODO: this is just a temporary hack until it is done properly by
# having an option on tsrchd etc. to just examine the listed structures.
# Output is then the output from tsrchd run on the reduced database.
#
# Environment variables:
#
#   PATH must contain the location of the Python scripts, ie where this
#   script itself is and the ptgraph/ directory with pytableaucreate.py etc.,
#   and the location of tsrchd_sparse.
#   The dssp program must also be in the PATH.
#
#   PYTHONPATH must contain the directory containing the ptsecstruct.py 
#   and other Python modules used by the Python scripts.
#
# $Id: run_tlocsd_filter_tsrchd.sh 2110 2009-03-18 05:58:44Z astivala $
# 
###############################################################################

TIME=/usr/bin/time
TLOCSD=tlocsd
TSRCHD=tsrchd_pardiso

# write tableau and distance matrix for tsrchd to stdout
# Parameters:
#     pdbfile - filename of PDB file
#     use_hk - if 1, use the HH and KK codes
#     extra_opts - extra options to add to pytableaucreate.py
writetableau() {
    pdbfile=$1
    use_hk=$2
    extra_opts="$3"
    tabopts="-b -35 -f -t dssp -p none ${extra_opts}"
    if [ $use_hk -eq 1 ]; then
        tabopts="${tabopts} -k"
    fi
    pytableaucreate.py ${tabopts} ${pdbfile}
}



use_ordering=1
sse_num_list=''
sse_num_list_opt=''
percent_hits=10
run_mustang=0

while getopts 'qe:p:' opt
do
    case $opt in
    q) use_ordering=0
    ;;
    e) sse_num_list="$OPTARG"
       sse_num_list_opt="-e ${sse_num_list}"
    ;;
    p)
       percent_hits="$OPTARG"
       ;;
    ?)
    echo "Usage: $0 [-q] [-e sse_num_list] [-p percent_hits] query_pdb_file dbdir" >&2
    exit 1
    ;;
    esac
done
shift $(($OPTIND - 1))


if [ $# -ne 2 ]; then
    echo "Usage: $0 [-q] [-e sse_num_list] [-p percent_hits] query_pdb_file dbdir" >&2
    exit 1
fi

querystruct=$1
dbdir=$2

db_size=`grep -c '^d' ${dbdir}/tableauxdistmatrixdb.ascii`
num_hits=`echo "$db_size * $percent_hits / 100" | bc`
#echo 'xxx num_hits =' $num_hits

queryinput1_tmp=/var/tmp/rtft1$$
tableauxdb_tmp=/var/tmp/rtftdb$$
queryinput2_tmp=/var/tmp/rtft2$$

echo ${dbdir}/tableauxdistmatrixdb.ascii >${queryinput1_tmp}
if [ $use_ordering -eq 0 ]; then
    echo "T F F" >> ${queryinput1_tmp}     # options: type,order,output
else
    echo "T T F" >> ${queryinput1_tmp}     # options: type,order,output
fi
extra_opts="${sse_num_list_opt}"
writetableau ${querystruct} 0 "" >> ${queryinput1_tmp}

trap "rm ${queryinput1_tmp} ${queryinput2_tmp} ${tableauxdb_tmp}" 0

${TIME} ${TLOCSD} < ${queryinput1_tmp} | sort -k2,2nr | head -${num_hits} |
        cut -d' ' -f1 | 
        convdb2.py -l ${dbdir}/tableauxdb.pickle ${dbdir}/distmatrixdb.pickle > ${tableauxdb_tmp}

# build input for tsrchd as same as original but using new (reduced) database
awk "NR == 1 {print(\"${tableauxdb_tmp}\")} NR > 1 {print}" < ${queryinput1_tmp} > ${queryinput2_tmp}

${TIME} ${TSRCHD} < ${queryinput2_tmp} 


