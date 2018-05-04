#!/bin/bash
#
# host.run.cudasatabsearch.querylist.pbs_sript.sh
#
# Usage: qsub host.run.cudasatabsearch.querylist.pbs_sript.sh
#
# run tlocsd
# on list of sids in database
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: host_run.cudasatabsearch.querylist.pbs_sript.sh 4764 2013-11-20 22:17:21Z astivala $
#

#PBS -q serial
#PBS -N host_qlist_gpucudaSaTabsearch
#PBS -l walltime=50:0:0


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

TIME=/usr/bin/time
TLOCSD=./cudaSaTabsearch

# the (compressed) tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii file 
# can be downloaded from
# http://munk.csse.unimelb.edu.au/~astivala/satabsearch/tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii.gz
# (then uncompress with gunzip)
tlocsdopts="-c -q ${HOME}/tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii -r8192"


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

