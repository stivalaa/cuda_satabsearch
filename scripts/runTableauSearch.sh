#!/bin/sh
#
# runTableauSearch.sh
#
# run Arun's TableauSearch 
# on a given PDB file against database (hardcoded in exercutable)
# saving results and stderr (for $TIME)
#
# Uses other scripts to build query angles file, convert output
# to format for evaluation with tsevalfn.py etc. (also runs the latter)
#
# Usage: runTableauSearch.sh querypdbfile outdir
#
# Puts output in outdir, with names starting with the querypdbfile baesname
# WARNING: overwrites output files
#
# $Id: runTableauSearch.sh 2958 2009-11-19 04:10:07Z astivala $
#

TIME=/usr/bin/time
TABLEAUSEARCH=/home/alexs/phd/TableauxCompare/bin/TableauComparer

if [ $# -ne 2 ]; then
    echo "Usage: $0 querypdbfile outdir" >&2
    exit 1
fi
infile=$1
outdir=$2
queryid=`basename $1`

anglesfile=${outdir}/${queryid}.angles
outfile=${outdir}/${queryid}.TableauSearch.out

pytableaucreate.py -e -35 -t dssp -p none ${infile} > ${anglesfile}

# NB TableauSearch always writes to search.scores in output dir
$TIME ${TABLEAUSEARCH} ${anglesfile} ${outdir} > ${outdir}/${queryid}.err 2>&1
tableausearchout2col.py < ${outdir}/search.scores > ${outfile}

cat ${outdir}/${queryid}.err

# run tsevalfn.py if it is a SCOP identifier
if [ `expr substr ${queryid} 1 1` = 'd' ]; then
    scopsid=`basename ${queryid} .ent`
    rtabfile=${outdir}/${scopsid}.TableauSearch.rtab
    tsevalfn.py $scopsid ${outfile} > ${rtabfile}
    grep AUC ${rtabfile}
fi

