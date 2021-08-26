#!/usr/bin/env python
###############################################################################
#
# convdbpacked2ascii.py - Convert PTTableauPacked tableaux db to vector form
#
# File:    convdbpacked2vector.py
# Author:  Alex Stivala
# Created: March2010
#
#
# Usage:
#    convdbpacked2vector.py inputdb > outputfile
#
# (output is to stdout)
#
# $Id: convdbpacked2vector.py 3528 2010-03-25 05:33:30Z alexs $
#
###############################################################################

import sys
import pickle
import numpy
from pttableau import PTTableauPacked


"""
This script converts tableaux database in pickled PTTableauPacked format
built by buildtableauxdb.py to the vector feature format used
FORTRAN-77 implementation of IR Tableau, tableau searching by
vector cosine similarity as described by:

Zhang, Bailey, Konagurthu, Ramamohanarao 2010 'A fast indexing approach
to protein structure comparison' BMC Bioinformatics 11(Suppl 1):S46D
8th Asia-Pacific Bioinformatics Conference (APBC 2010)

The format of the input is pickled PTTableauPacked objects built by
buildtableauxdb.py.

The output 'database' is an ASCII file where each line contains
the structure identifier and the 32-element integer vector representing
the tableau for that structure reduced to the 'bag of words' vector
as described in the paper cited above.


The output is to stdout.

"""


# TYPE_CODE_DICT gives an index 0..32 for each combinatino of a pair
# of SSE type codes ('ee' is two strands, 'ex' is strand/helix, etc.)
# and tableau code
TYPE_CODE_DICT = \
{
    ('xx','PE') : 0,
    ('xx','PD') : 1,
    ('xx','RD') : 2,
    ('xx','RT') : 3,
    ('xx','OT') : 4,
    ('xx','OS') : 5,
    ('xx','LS') : 6,
    ('xx','LE') : 7,

    ('xe','PE') : 8,
    ('xe','PD') : 9,
    ('xe','RD') : 10,
    ('xe','RT') : 11,
    ('xe','OT') : 12,
    ('xe','OS') : 13,
    ('xe','LS') : 14,
    ('xe','LE') : 15,

    ('ex','PE') : 16,
    ('ex','PD') : 17,
    ('ex','RD') : 18,
    ('ex','RT') : 19,
    ('ex','OT') : 20,
    ('ex','OS') : 21,
    ('ex','LS') : 22,
    ('ex','LE') : 23,

    ('ee','PE') : 24,
    ('ee','PD') : 25,
    ('ee','RD') : 26,
    ('ee','RT') : 27,
    ('ee','OT') : 28,
    ('ee','OS') : 29,
    ('ee','LS') : 30,
    ('ee','LE') : 31
}

def tableau_to_feature_vector(tableau):
    """
    Given a tableau in PTTableau format, convert to feature vector
    (numpy.array vector of integers) of dimension 32, in which each
    element is the count of the number of occurrences of one of the
    32 (4 different SSE combinatinos and 8 differet two-character
    tableau codes) possible types, as described by the Zhou et al 2010
    IR Tableau paper.

    Paramteters:
       tableau - tableau in PTTableau format

    Return value:
       numpy.array shape (32) vector of integers

    Uses global constant TYPE_CODE_DICT to get index for the type/tabcode
    combination
    """
    fvec = numpy.zeros(32, dtype='i')
    for i in xrange(len(tableau)):
        for j in xrange(i+1,len(tableau)):
            if tableau[(i,i)][0] == 'e':
                ssetypecode = 'e'
            else:
                ssetypecode = 'x'
            if tableau[(j,j)][0] == 'e':
                ssetypecode += 'e'
            else:
                ssetypecode += 'x'
            tabcode = tableau[(i,j)]
            if tabcode == '??':  # sometimes cannot get angle, ignore
                continue
            vecindx = TYPE_CODE_DICT[(ssetypecode, tabcode)]
            fvec[vecindx] += 1
    return fvec
    
    
def usage(prog):
    """
    print usage message and exit
    """
    sys.stderr.write("Usage: " + prog + " inputdb > outputfile\n")
    sys.exit(1)
    

def main():
    """
    main for convdbpacked2vector.py - load PTTableauPacked pickle and 
    output feature vectors in ASCII
    """
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    dbfile = sys.argv[1]
    sys.stderr.write('loading tableaux... ')
    db = pickle.load(open(dbfile))
    sys.stderr.write('done\n')
    for pdbid,dbtablist in db.iteritems():
        tabnum = 0
        while tabnum < len(dbtablist):
            tableau = dbtablist[tabnum]
            n = len(tableau)
            name = pdbid
            if len(dbtablist) > 1:
                name += str(tabnum)
            sys.stdout.write('%s ' % name)
            fvec = tableau_to_feature_vector(tableau)
            sys.stdout.write(" ".join(str(x) for x in fvec))
            sys.stdout.write('\n')
            tabnum += 1

            
if __name__ == "__main__":
    main()
