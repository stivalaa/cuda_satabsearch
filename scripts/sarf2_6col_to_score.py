#!/usr/bin/env python
#
# File:    sarf2_6col_to_score.sh
# Author:  Alex Stivala
# Created: April 2010
#
# sarf2_6col_to_score.sh - Convert SARF2  6 column format from sarf2out6col
#                          to 2 column format for tsevalfn.py etc.
#                  
#
# Usage: sarf2out6col < results.sarf2 | sarf2_6col_to_score.py 
#
# The input in 6 column (id1 id2 len1 len2 nres rmsd) is read from stdin
# Output is to stdout.
#
#
# $Id: sarf2_6col_to_score.py 3695 2010-05-18 07:32:14Z alexs $
#


import sys,os
from itertools import groupby


def Q_ssm(nres, rmsd, n1, n2):
    """
    Compute the SSM 'Q' score for an alignment of two proteins with n1
    and n2 residues, with nres residues aligned and RMSD value of rmsd
    (Angstroms).

    This score is defined in

    Kirssinel, E. & Henrick, K. 2004 'Secondary-structure matching (SSM), a new
    tool for fast protein structure alignment in three dimensions'
    Acta Crystallographica D60:2256-2268

    Parameters:
         nres - number of residues in alignment
         rmsd - root mean square deviation of aligned residues (Angstroms)
         n1 - number of residues in protein 1
         n2 - number of residues in protein 2

    Return value: Q score for alignment
        
    """
    R0 = 3.0 # Krissinel & Henrick p. 2262
    return nres**2 / ( (1 + (rmsd / R0)**2) * n1 * n2)

#
# main
#

if len(sys.argv) != 1:
    usage(os.path.basename(sys.argv[0]))

sarflist = [] # list of (id1,id2,size1,size2,nres,rmsd)
for line in sys.stdin:
    (id1, id2, size1, size2, nres, rmsd) = line.split()
    size1 = int(size1)
    size2 = int(size2)
    nres = int(nres)
    rmsd = float(rmsd)
    sarflist.append((id1, id2, size1, size2, nres, rmsd))

query_group_iter = groupby(sorted(sarflist), lambda t : t[0])

for (queryid, result_iter) in query_group_iter:
    sys.stdout.write("# QUERY ID = %s\n" % queryid)
    for (id1, id2, size1, size2, nres, rmsd) in result_iter:
        if size1 == 0 or size2 ==0: # somtimes we get junk results from SARF
            sys.stderr.write("WARNING: bad results for %s - %s\n" % (id1, id2))
            score = 0
        else:
            score = Q_ssm(nres, rmsd, size1, size2)
        sys.stdout.write("%s %f\n" % (id2, score))





    

