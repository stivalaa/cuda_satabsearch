#!/bin/sh
#
# build the tableau + distance matrix for tsrchd_sparse for the 
# beta-grasp query (as the 4 largest strands and 1 alpha hexli in 
# ubiquitin structure)
#
# $Id: build_betagrasp_query.sh 2908 2009-11-06 05:33:18Z astivala $

# tableaux+distmatrix db file
TABLEAUX_DB=${HOME}/tableauxdistmatrixdb.ascii

echo "${TABLEAUX_DB}" 
echo "T T F" # options: type,order,output
pytableaucreate.py -bf -35 -tdssp -p none -i BGRASP -s2,1,8,5,3 ${HOME}/pdb/d1ubia_.ent  
