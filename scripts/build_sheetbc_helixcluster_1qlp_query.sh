#!/bin/sh
#
# build the tableau + distance matrix for tsrchd_sparse for the 
# serpin B/C sheet. 
#
# $Id: build_sheetbc_helixcluster_1qlp_query.sh 2108 2009-03-18 05:49:01Z astivala $

# tableaux+distmatrix db file
TABLEAUX_DB=/local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii

echo "${TABLEAUX_DB}" 
echo "T T F" # options: type,order,output
pytableaucreate.py -pnone -b -f -35 -i sheetbch -tdssp -s2,3,4,5,26,25,15,14,12,13,18,24 /local/charikar/pdb/pdb/ql/pdb1qlp.ent.gz 
