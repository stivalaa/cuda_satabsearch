#!/usr/bin/env python
###############################################################################
#
# norms.py - functions to normalize tableau match score by protein size
#
# File:    norms.py
# Author:  Alex Stivala
# Created: September 2008
#
# $Id: norms.py 3585 2010-04-29 03:56:03Z alexs $
# 
###############################################################################
"""
 Contains function to normalize tableua matching scores by number
 of SSEs, analogous to the normalizatino functions used in Pelta et al 2008
 BMC Bioinformatics 9:161

 Note comments on these functions refer to #SSEs for normalizing tableau
 match score, but actually they can equally well be used for normalizing
 MAX-CMO scores where size is number of contacts, as original use
 in Pelta et al 2008.

 Also note using these functions with QP tableau search (size is
 number of SSES) does not 'normalize' the results into [0,1] either,
 since max score is actually n(n-1) where n is #SSEs in smaller structure.
 TODO: should have different normalization functions for tableau search

"""
import sys,os,glob
import numpy.oldnumeric as Numeric



def norm1(score, size1, size2):
    """
    Normalization similar to norm1 in Pelta et al 2008 (from Lancia et al
    2006) for MAX-CMO.

    norm1(struct1,struct2) = tabmatch_score(struct1,struct2) / 
                               min(#sses(struct1), #sses(struct2))

    Parameters:
        score  - tableau match score for the two structures
        size1  - number of SSEs in one structure
        size2  - number of SSEs in other structure

    Return value:
        normalized score as described above.
    """
    return score / float(min(size1, size2))
#    n = min(size1,size2)
#    if n < 2:
#        n = 2
#    return score / float(n*(n-1))

    
def norm2(score, size1, size2):
    """
    Normalization similar to norm2 in Pelta et al 2008 (from Xie & Sahinidis
    2006) for MAX-CMO.

    norm1(struct1,struct2) = 2*tabmatch_score(struct1,struct2) / 
                               (#sses(struct1) + #sses(struct2))

    Parameters:
        score  - tableau match score for the two structures
        size1  - number of SSEs in one structure
        size2  - number of SSEs in other structure

    Return value:
        normalized score as described above.
    """
    return 2 * score / float(size1 + size2)


def norm3(score, size1, size2):
    """
    Normalization similar to norm3 in Pelta et al 2008 for MAX-CMO.

    norm1(struct1,struct2) = 0 if #SSE differnce > 75%
                             norm1(struct1,struc2) otherwise

    Parameters:
        score  - tableau match score for the two structures
        size1  - number of SSEs in one structure
        size2  - number of SSEs in other structure

    Return value:
        normalized score as described above.
    """
    if float(abs(size1 - size2)) / float(max(size1,size2)) > 0.75:
        return 0
    else:
        return norm1(score, size1, size2)
    

