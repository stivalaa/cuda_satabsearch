#!/bin/bash
#
# gtx285.run.cudasatabsearch.querylist.pbs_sript.sh
#
# Usage: qsub gtx285.run.cudasatabsearch.querylist.pbs_sript.sh
#
# run tlocsd
# on list of sids in database
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: gtx285.run.cudasatabsearch.querylist.pbs_sript.sh 3356 2010-02-19 04:30:27Z alexs $
#

# run on GTX285 (nv3)
#PBS -l nodes=nv3
#PBS -q gpu
#PBS -N gtx285_qlist_gpucudaSaTabsearch
#PBS -l walltime=10:0:0

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

TIME=/usr/bin/time
TLOCSD=./cudaSaTabsearch

tlocsdopts="-q /home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75/tableauxdistmatrixdb.sorted.ascii -r8192"


echo "# Run as: " $0 $* 
echo "# at: " `date`
echo "# on: " `uname -a`

$TIME ${TLOCSD} ${tlocsdopts}  <<EOF
d1ndda_
d1wnda_
d1rrva_
d1muka_
d2byla1
d1a32a_
d1mjta_
EOF

