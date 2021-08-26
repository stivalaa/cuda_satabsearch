#!/bin/sh
#
# build the tableau + distance matrix for tsrchd_sparse for the 
# serpin B/C sheet. 
#
# $Id: build_sheetbc_1qlp_query.sh 2105 2009-03-17 01:30:33Z astivala $

# tableaux+distmatrix db file
TABLEAUX_DB=/local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii

echo "${TABLEAUX_DB}" 
echo "T T F" # options: type,order,output
pytableaucreate.py -pnone -b -f -35 -i sheetbc -tdssp -s2,26,25,15,14,12,13,18,24 /local/charikar/pdb/pdb/ql/pdb1qlp.ent.gz 
