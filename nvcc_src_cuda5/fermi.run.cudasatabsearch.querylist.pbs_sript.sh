#!/bin/bash
#
# fermi.run.cudasatabsearch.querylist.pbs_sript.sh
#
# Usage: qsub fermi.run.cudasatabsearch.querylist.pbs_sript.sh
#
# run tlocsd
# on list of sids in database
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: fermi.run.cudasatabsearch.querylist.pbs_sript.sh 4760 2013-11-20 04:55:07Z astivala $
#

# run on Fermi architecture (Tesla M2070) 
#PBS -q gpu
#PBS -N fermi_qlist_gpucudaSaTabsearch
#PBS -l walltime=1:0:0

module load cuda/5.0.35

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

TIME=/usr/bin/time
TLOCSD=./cudaSaTabsearch

# the (compressed) tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii file 
# can be downloaded from
# http://munk.csse.unimelb.edu.au/~astivala/satabsearch/tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii.gz
# (then uncompress with gunzip)
tlocsdopts="-q ${HOME}/tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii -r8192"


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

