#!/bin/bash
#
# host.run.cudasatabsearch.querylist.slurm_sript.sh
#
# Usage: sbatch host.run.cudasatabsearch.querylist.slurm_sript.sh
#
# run tlocsd
# on list of sids in database
# saving results and stderr (for $TIME)
# in cwd (so run from results/tlocsd/ dir)
# and running one at a time.
#
# $Id: host_run.cudasatabsearch.querylist.pbs_sript.sh 4764 2013-11-20 22:17:21Z astivala $
#
#SBATCH --job-name="host_qlist_cudaSaTabsearch"
#SBATCH --ntasks=1
#SBATCH --time=0-48:00:00
#SBATCH --partition=slim
#SBATCH --output=host_qlist_cudaSaTabsearch-%j.out
#SBATCH --error=host_qlist_cudaSaTabsearch-%j.err

TIME=/usr/bin/time
TLOCSD=./cudaSaTabsearch

# the (compressed) tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii file 
# can be downloaded from
# http://munk.cis.unimelb.edu.au/~astivala/satabsearch/tableauxdistmatrixdb-sel-gs-bib-95-1.75.sorted.ascii.gz
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

