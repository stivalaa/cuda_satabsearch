#!/usr/bin/env python
###############################################################################
#
# slrtabs2star.py - convert .slrtab files into format for StAR
#
# File:    slrtabs2star.py
# Author:  Alex Stivala
# Created: April 2010
#
# $Id: slrtabs2star.py 3593 2010-05-02 05:02:35Z alexs $
###############################################################################
"""
 slrtabs2star.py - convert .slrtab files into format for StAR

 Usage: slrtabs2star.py [-v] posfile negfile < listfile


   posfile is the postivies filename to create
   negfile is the negatives filename to create

   -v specifies verbose output to stderr
 
 listfile read from stdin is a two-column tab-delimited file with method namd
 in first column and slrtab filename in second column. The script
 reads each of the slrtab files named and creates the postivies and 
 negatives file for StAR.

 The listfile looks like e.g.

"QP Tableau Search" /astivala/phd/qptabsearch/results/query200/norm2/query200_roc.slrtab
 
 NB tab delimiter not space

 Converts the .slrtab files, one for each method,
 created by rocrfischer.py, tsevalfn.py etc.
 into the positives.dat and negatives.dat files used by
 StAR.
 Both files are tab-delimited, with names on first line,
 and scores for the true negative and positive classes (for negatives
 and postivies file respectively) for each
 classifier (according to order of names on first line) on each
 subsequent line. (See examples in StAR documentation).

 WARNING: the positives file and negatives files are overwritten
 if they exist.

 The reference for StAR is:

 Vergara, Normabuena, Ferrada, Slater and Melo 'StAR: a simple tool for
 the statistical comparison of ROC curves' BMC Bioinformatics 2008 9:265

"""


import sys,os
import getopt

from tsevalutils import iter_slrtab


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage " + progname + " [-v] posfile negfile < listfile\n")
    sys.exit(1)

def main():
    """
    main for slrtabs2star.py
    see usage message in header comment
    """
    verbose = False
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "v?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-v": # verbose
            verbose = True
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    posfilename = args[0]
    negfilename = args[1]

    (namelist, slrtabfilelist) = zip(*[line.split('\t') for line in sys.stdin])  # trciky use of zip and * to unzip list
    slrtabfilelist = [x[:-1] for x in slrtabfilelist] # remove newlines on end

    posfile_fh = open(posfilename, "w")
    negfile_fh = open(negfilename, "w")

    posscores = [] # list of lists: each list is scores for each method in pos class
    negscores = [] # similarly for negative class scores
    for slrtabfile in slrtabfilelist:
        if verbose:
            sys.stderr.write("Reading results from file %s..." % slrtabfile)
        slrlist = list(iter_slrtab(open(slrtabfile))) # (score,label) list
        posscores.append([sl[0] for sl in slrlist if sl[1] == 1])
        negscores.append([sl[0] for sl in slrlist if sl[1] == 0])
        assert(len(posscores[-1]) + len(negscores[-1]) == len(slrlist))
        if verbose:
            sys.stderr.write(" %d entries (%d pos, %d neg)\n" % (len(slrlist),len(posscores[-1]),len(negscores[-1])))
        
    if verbose:
        sys.stderr.write("writing output to %s and %s..." %(posfilename, negfilename))
    
    posfile_fh.write('\t'.join(namelist) + '\n')
    negfile_fh.write('\t'.join(namelist) + '\n')

    numpos = len(posscores[0]) # FIXME may be different lengths
    for i in xrange(numpos):
        for j in xrange(len(namelist)):
            posfile_fh.write(str(posscores[j][i]))
            if j < len(posscores) - 1:
                posfile_fh.write('\t')
        posfile_fh.write('\n')

    numneg = len(negscores[0]) # FIXME may be different lengths
    for i in xrange(numneg):
        for j in xrange(len(namelist)):
            negfile_fh.write(str(negscores[j][i]))
            if j < len(negscores) - 1:
                negfile_fh.write('\t')
        negfile_fh.write('\n')


    posfile_fh.close()
    negfile_fh.close()
    if verbose:
        sys.stderr.write("done\n")

    
if __name__ == "__main__":
    main()
